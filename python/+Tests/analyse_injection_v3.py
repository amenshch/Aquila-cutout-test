#!/usr/bin/env python
"""analyse_injection_v3.py -- analyse the getsf extractions of the two
population-matched injection runs.

The two runs place identical models at identical positions and differ only in
whether the cloud background under each injected core was flattened before
insertion.  Comparing them therefore isolates the contribution of the cloud's
own fluctuations to the error in the recovered mass.

Two questions are answered separately.

1.  Does the injected population look like the real one?  This is the question
    the placement method was built to solve, and it is now answerable
    like-for-like, because the injected cores have been through getsf and
    carry their own PEAK^BGF, PEAK^SBF and AFWHM, measured exactly as the real
    sources' were.  Until now it could only be checked against the truth
    table, whose background column is a disk average over a different scale.

2.  How well is the mass recovered and corrected?  Each recovered source is
    matched to its injected core by position, and the reported mass and the
    two corrected masses are compared with the true Bonnor-Ebert mass from the
    truth table.

Usage:
    python analyse_injection_v3.py [--dir inj2] [--out-prefix inj_v3]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '.')
sys.path.insert(0, 'python')
import mass_correction_pipeline_6 as mc                    # noqa: E402
from getsf_columns import GetsfTable                       # noqa: E402
from scipy.stats import ks_2samp, spearmanr               # noqa: E402

DIST = 260.0
MATCH_PIX = 10.0
PIX_ARCSEC = 3.0
RUNS = [('unflattened', 'cmf2.00_n600_s2029'),
        ('flattened', 'cmf2.00_n600_s2029_flat')]


def read_truth(path):
    """id, tag, x_pix, y_pix, local_Sigma, SD_emb, T_BE, rho_BE, M_BE,
    R_BE_as, FWHM_SD, scale."""
    out = []
    for ln in open(path):
        if not ln.strip() or ln.lstrip()[0] in '#!':
            continue
        f = ln.split()
        out.append(dict(tag=f[1], x=float(f[2]), y=float(f[3]),
                        sigma=float(f[4]), sd_emb=float(f[5]),
                        mbe=float(f[8]), rbe=float(f[9]),
                        fwhm=float(f[10]) if len(f) > 11 else np.nan))
    return out


def col(t, name):
    """First occurrence of a column, for catalogs that merge several files."""
    v = np.asarray(t.col(name), float)
    if len(v) == t.nrows:
        return v
    return v.reshape(t.nrows, -1)[:, 0]


def analyse(run_dir, label, target):
    truth = read_truth([os.path.join(run_dir, f) for f in os.listdir(run_dir)
                        if f.endswith('_truth.txt')][0])
    cats = [os.path.join(run_dir, 'Aquila.s.sources.ok.cat'),
            os.path.join(run_dir, 'Aquila.s.sources.ok.add.cat')]
    table, _ = mc.correct_cloud(cats, DIST)

    t = GetsfTable(cats[0])
    xs, ys = col(t, 'XCO_P'), col(t, 'YCO_P')
    no_all = col(t, 'NO')
    afw = col(t, 'AFWHM03')          # 13.5" column-density band in this setup
    pos = {int(n): (x, y) for n, x, y in zip(no_all, xs, ys)}
    size = {int(n): a for n, a in zip(no_all, afw)}

    tx = np.array([r['x'] for r in truth])
    ty = np.array([r['y'] for r in truth])

    rows = []
    used = {}
    for i, no in enumerate(table['NO']):
        no = int(no)
        if no not in pos:
            continue
        xp, yp = pos[no]
        d = np.hypot(tx - xp, ty - yp)
        j = int(np.argmin(d))
        if d[j] >= MATCH_PIX:
            continue
        if j in used and used[j][0] <= d[j]:
            continue
        used[j] = (d[j], i, no)
    for j, (d, i, no) in used.items():
        rows.append(dict(
            mbe=truth[j]['mbe'], fwhm_model=truth[j]['fwhm'],
            m_init=table['M_init'][i], m3=table['M_corr'][i],
            m4=table['M_corr_4D'][i],
            sigma=table['sigma'][i], conc=table['conc'][i],
            afwhm=size.get(no, np.nan),
            contrast=1.0 + table['conc'][i] * 0,   # filled below
            in_dom=bool(table['IN_DOMAIN'][i]),
            peak_src=np.nan))
    # contrast from the catalog columns, as for the real sources
    g = mc.load_getsf(cats)
    okm = mc.source_mask(g)
    cno = {int(n): k for k, n in enumerate(g['NO'])}
    for r, (j, (d, i, no)) in zip(rows, used.items()):
        k = cno.get(int(no))
        if k is not None and g['Nbg'][k] > 0:
            r['contrast'] = 1.0 + g['peak'][k] / g['Nbg'][k]
            r['sigma'] = g['Nbg'][k]

    n_inj = len(truth)
    n_rec = len(rows)
    print('\n%s' % ('=' * 74))
    print('%s injection: %d cores injected, %d recovered and matched (%.0f%%)'
          % (label.capitalize(), n_inj, n_rec, 100.0 * n_rec / n_inj))
    print('%s' % ('=' * 74))

    S = np.array([r['sigma'] for r in rows])
    C = np.array([r['contrast'] for r in rows])
    A = np.array([r['afwhm'] for r in rows])
    Mr = np.array([r['m_init'] for r in rows])
    good = np.isfinite(S) & (S > 0) & np.isfinite(C) & (C > 1)

    Mt, St, Ct, At = target
    print('\nDoes the recovered injected population match the real one?')
    print('  quantity                          injected      real   '
          'KS / probability')
    for nm, vi, vr in (('background column Sigma (cm^-2)', S[good], St),
                       ('contrast PEAK^SBF/PEAK^BGF', C[good], Ct),
                       ('reported SED mass (M_sun)', Mr[good], Mt),
                       ('measured size AFWHM (arcsec)',
                        A[np.isfinite(A) & (A > 0)], At)):
        k = ks_2samp(np.log10(vi), np.log10(vr))
        print('  %-32s %9.4g %9.4g   %.3f / %.2g'
              % (nm, np.median(vi), np.median(vr), k.statistic, k.pvalue))
    print('  rank correlation of contrast on column: injected %+.3f, '
          'real %+.3f'
          % (spearmanr(np.log10(S[good]), np.log10(C[good])).statistic,
             spearmanr(np.log10(St), np.log10(Ct)).statistic))

    print('\nMass recovery and correction, ratio to the true mass M_BE')
    print('  %-26s %5s %9s %9s %10s %11s'
          % ('treatment', 'n', 'median', 'scatter', 'within 2x', 'med |log10|'))
    a = np.array([[r['mbe'], r['m_init'], r['m3'], r['m4']] for r in rows])
    dom = np.array([r['in_dom'] for r in rows])
    res = {}
    for subset, m in (('all matched', np.ones(len(a), bool)),
                      ('inside the domain', dom)):
        if m.sum() == 0:
            continue
        print('  %s (n=%d)' % (subset, int(m.sum())))
        for nm, jj in (('uncorrected', 1), ('three observables', 2),
                       ('four observables', 3)):
            r = a[m, jj] / a[m, 0]
            gg = np.isfinite(r) & (r > 0)
            lr = np.log10(r[gg])
            print('    %-24s %5d %9.3f %9.2f %9.0f%% %11.3f'
                  % (nm, gg.sum(), np.median(r[gg]), 10 ** lr.std(),
                     100 * np.mean((r[gg] > 0.5) & (r[gg] < 2)),
                     np.median(np.abs(lr))))
            res[(subset, nm)] = (np.median(r[gg]), 10 ** lr.std())
    return rows, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='inj2')
    ap.add_argument('--target', default='python/joint_target_Aquila_v3.txt')
    a = ap.parse_args()
    mc.VERBOSE = 0

    d = np.loadtxt(a.target)
    target = (d[:, 0], d[:, 1], d[:, 2], d[:, 3])
    print('real reference: %d Aquila sources' % len(d))

    out = {}
    for label, sub in RUNS:
        p = os.path.join(a.dir, sub)
        if not os.path.isdir(p):
            print('missing: %s' % p)
            continue
        out[label] = analyse(p, label, target)

    if len(out) == 2:
        print('\n%s' % ('=' * 74))
        print('Flattened minus unflattened: the effect of the cloud '
              'fluctuations')
        print('%s' % ('=' * 74))
        ru, _ = out['unflattened']
        rf, _ = out['flattened']
        print('  cores recovered: unflattened %d, flattened %d'
              % (len(ru), len(rf)))
        for subset in ('all matched', 'inside the domain'):
            for nm in ('uncorrected', 'three observables',
                       'four observables'):
                u = out['unflattened'][1].get((subset, nm))
                f = out['flattened'][1].get((subset, nm))
                if u and f:
                    print('  %-18s %-20s median %.3f -> %.3f, '
                          'scatter %.2f -> %.2f'
                          % (subset, nm, u[0], f[0], u[1], f[1]))


if __name__ == '__main__':
    main()
