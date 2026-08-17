#!/usr/bin/env python
"""analyse_quad_injection.py -- analyse a four-way injection set.

The four runs share one placement: identical models at identical positions,
differing only in whether the cloud background under each core was flattened
and whether the cores that overlap a neighbour were removed.  Every core
therefore appears in each run in which it survives the isolation filter, and
the four can be compared core by core:

    as placed              cloud fluctuations present, neighbours present
    flattened              cloud fluctuations removed, neighbours present
    isolated               cloud fluctuations present, neighbours removed
    flattened and isolated both removed

The differences separate the two contributions to the error in the recovered
mass.  Comparing the flattened pair at fixed isolation gives the contribution
of the cloud fluctuations; comparing the isolated pair at fixed flattening
gives the contribution of blending between injected cores.

Comparisons are made on the cores recovered in BOTH members of a pair, never
sample against sample.  The distinction is not pedantic: an earlier comparison
of a flattened pair made sample against sample suggested the fluctuations
biased the recovered mass high by 39%, and the effect vanished once the same
cores were compared in both, the apparent bias having been sampling.

Usage:
    python analyse_quad_injection.py [--dir inj3] [--seed s2030]
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '.')
sys.path.insert(0, 'python')
import mass_correction_pipeline_6 as mc                    # noqa: E402
from getsf_columns import GetsfTable                       # noqa: E402
from scipy.stats import ks_2samp, spearmanr, wilcoxon      # noqa: E402

DIST = 260.0
MATCH_PIX = 10.0
VARIANTS = [('as placed', ''), ('flattened', '_flat'),
            ('isolated', '_iso'), ('flattened+isolated', '_flat_iso')]


def _col(t, name):
    v = np.asarray(t.col(name), float)
    return v if len(v) == t.nrows else v.reshape(t.nrows, -1)[:, 0]


def read_truth(path):
    out = []
    for ln in open(path):
        if not ln.strip() or ln.lstrip()[0] in '#!':
            continue
        f = ln.split()
        out.append(dict(key=(f[1], int(round(float(f[2]))),
                             int(round(float(f[3])))),
                        tag=f[1], x=float(f[2]), y=float(f[3]),
                        sigma=float(f[4]), mbe=float(f[8]),
                        rbe=float(f[9]),
                        fwhm=float(f[10]) if len(f) > 11 else np.nan))
    return out


def load_run(run_dir):
    """{(tag, x, y): recovered quantities} for one run, keyed so that the same
    core can be identified across runs without relying on row order."""
    truth = read_truth(glob.glob(os.path.join(run_dir, '*_truth.txt'))[0])
    cats = [os.path.join(run_dir, 'Aquila.s.sources.ok.cat'),
            os.path.join(run_dir, 'Aquila.s.sources.ok.add.cat')]
    table, _ = mc.correct_cloud(cats, DIST)
    t = GetsfTable(cats[0])
    xs, ys, nos = _col(t, 'XCO_P'), _col(t, 'YCO_P'), _col(t, 'NO')
    afw = _col(t, 'AFWHM03')
    pos = {int(n): (x, y, a) for n, x, y, a in zip(nos, xs, ys, afw)}
    tx = np.array([r['x'] for r in truth])
    ty = np.array([r['y'] for r in truth])
    mc.load_getsf._dist = DIST
    g = mc.load_getsf(cats)
    cno = {int(n): k for k, n in enumerate(g['NO'])}
    best = {}
    for i, no in enumerate(table['NO']):
        no = int(no)
        if no not in pos:
            continue
        xp, yp, aa = pos[no]
        d = np.hypot(tx - xp, ty - yp)
        j = int(np.argmin(d))
        if d[j] >= MATCH_PIX:
            continue
        k = truth[j]['key']
        if k in best and best[k]['d'] <= d[j]:
            continue
        c = cno.get(no)
        best[k] = dict(d=d[j], mbe=truth[j]['mbe'], fwhm=truth[j]['fwhm'],
                       m_init=table['M_init'][i], m3=table['M_corr'][i],
                       m4=table['M_corr_4D'][i],
                       in_dom=bool(table['IN_DOMAIN'][i]), afwhm=aa,
                       sigma=(g['Nbg'][c] if c is not None else np.nan),
                       contrast=(1.0 + g['peak'][c] / g['Nbg'][c]
                                 if c is not None and g['Nbg'][c] > 0
                                 else np.nan))
    return truth, best


def stats_line(name, r):
    g = np.isfinite(r) & (r > 0)
    if g.sum() == 0:
        return '    %-22s no sources' % name
    lr = np.log10(r[g])
    return ('    %-22s n=%3d  median %.3f  scatter %.2f  within 2x %3.0f%%  '
            'median |log10| %.3f'
            % (name, g.sum(), np.median(r[g]), 10 ** lr.std(),
               100 * np.mean((r[g] > 0.5) & (r[g] < 2)), np.median(np.abs(lr))))


def compare(A, B, nameA, nameB, what):
    common = sorted(set(A) & set(B))
    if len(common) < 5:
        print('  %s: only %d cores in common' % (what, len(common)))
        return
    mbe = np.array([A[k]['mbe'] for k in common])
    print('\n  %s  (%d cores recovered in both)' % (what, len(common)))
    for label, D in ((nameA, A), (nameB, B)):
        print('   %s' % label)
        for nm, key in (('uncorrected', 'm_init'),
                        ('three observables', 'm3'),
                        ('four observables', 'm4')):
            print(stats_line(nm, np.array([D[k][key] for k in common]) / mbe))
    print('   paired difference in log10 of the ratio, %s minus %s'
          % (nameA, nameB))
    for nm, key in (('uncorrected', 'm_init'), ('three observables', 'm3'),
                    ('four observables', 'm4')):
        ra = np.array([A[k][key] for k in common]) / mbe
        rb = np.array([B[k][key] for k in common]) / mbe
        g = np.isfinite(ra) & np.isfinite(rb) & (ra > 0) & (rb > 0)
        dl = np.log10(ra[g]) - np.log10(rb[g])
        try:
            p = wilcoxon(dl).pvalue
        except Exception:
            p = np.nan
        sa, sb = np.log10(ra[g]).std(ddof=1), np.log10(rb[g]).std(ddof=1)
        print('    %-22s median %+.3f dex (p=%.2g); scatter %.2f versus %.2f'
              % (nm, np.median(dl), p, 10 ** sa, 10 ** sb))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='inj3')
    ap.add_argument('--stem', default='cmf2p00_n600_s2030', nargs='+',
                    help='one or more run stems; several are pooled, the key '
                         'of each core being prefixed by its stem so that '
                         'cores of different seeds are never confused')
    ap.add_argument('--real', default='python/joint_target_Aquila_v3.txt')
    a = ap.parse_args()
    mc.VERBOSE = -1

    stems = a.stem if isinstance(a.stem, list) else [a.stem]
    runs = {}
    for label, suf in VARIANTS:
        truth_all, best_all = [], {}
        for stem in stems:
            p = os.path.join(a.dir, stem + suf)
            if not os.path.isdir(p):
                print('missing: %s' % p)
                continue
            tr, bs = load_run(p)
            truth_all += tr
            # prefix the key with the stem: two seeds can place the same model
            # on the same pixel, and pooling on the bare key would merge them
            for k, v in bs.items():
                best_all[(stem,) + k] = v
        if truth_all:
            runs[label] = (truth_all, best_all)
    if len(stems) > 1:
        print('pooling %d seeds: %s\n' % (len(stems), ', '.join(stems)))

    print('%-22s %9s %11s %11s' % ('run', 'injected', 'recovered', 'fraction'))
    for label, _ in VARIANTS:
        if label not in runs:
            continue
        truth, best = runs[label]
        print('%-22s %9d %11d %10.0f%%'
              % (label, len(truth), len(best), 100.0 * len(best) / len(truth)))

    # ---- population match, on the run as placed ------------------------
    if 'as placed' in runs:
        d = np.loadtxt(a.real)
        Mr, Sr, Cr, Ar = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
        best = runs['as placed'][1]
        S = np.array([v['sigma'] for v in best.values()])
        C = np.array([v['contrast'] for v in best.values()])
        A = np.array([v['afwhm'] for v in best.values()])
        g = np.isfinite(S) & (S > 0) & np.isfinite(C) & (C > 1)
        print('\nRecovered injected cores against the real Aquila sources')
        print('  %-30s %10s %10s %14s'
              % ('quantity (median)', 'injected', 'real', 'KS / probability'))
        for nm, vi, vr in (('background column (cm^-2)', S[g], Sr),
                           ('contrast', C[g], Cr),
                           ('measured size AFWHM (arcsec)',
                            A[np.isfinite(A) & (A > 0)], Ar)):
            ks = ks_2samp(np.log10(vi), np.log10(vr))
            print('  %-30s %10.4g %10.4g %8.3f / %.2g'
                  % (nm, np.median(vi), np.median(vr), ks.statistic, ks.pvalue))
        print('  rank correlation of contrast on column: injected %+.3f, '
              'real %+.3f'
              % (spearmanr(np.log10(S[g]), np.log10(C[g])).statistic,
                 spearmanr(np.log10(Sr), np.log10(Cr)).statistic))

    # ---- the two paired comparisons ------------------------------------
    print('\n%s\nEffect of the cloud fluctuations (flattening), at fixed '
          'isolation\n%s' % ('=' * 74, '=' * 74))
    if 'as placed' in runs and 'flattened' in runs:
        compare(runs['as placed'][1], runs['flattened'][1],
                'as placed', 'flattened', 'with neighbours present')
    if 'isolated' in runs and 'flattened+isolated' in runs:
        compare(runs['isolated'][1], runs['flattened+isolated'][1],
                'isolated', 'flattened+isolated', 'with neighbours removed')

    print('\n%s\nEffect of blending between injected cores (isolation), at '
          'fixed flattening\n%s' % ('=' * 74, '=' * 74))
    if 'as placed' in runs and 'isolated' in runs:
        compare(runs['as placed'][1], runs['isolated'][1],
                'as placed', 'isolated', 'with the cloud fluctuations present')
    if 'flattened' in runs and 'flattened+isolated' in runs:
        compare(runs['flattened'][1], runs['flattened+isolated'][1],
                'flattened', 'flattened+isolated',
                'with the cloud fluctuations removed')


if __name__ == '__main__':
    main()
