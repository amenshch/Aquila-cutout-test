#!/usr/bin/env python
"""paired_injection_v3.py -- compare the flattened and unflattened injection
runs on the cores recovered in BOTH, and test whether the recovered mass
tracks the local cloud fluctuation.

The two runs place identical models at identical positions and differ only in
whether the cloud background under each injected core was flattened before
insertion.  Comparing them on the cores recovered in both makes the comparison
exactly paired: the same model, the same position, the same true mass, so the
difference is the cloud fluctuation and nothing else.

Two questions are answered.

1.  Paired effect of the fluctuations.  For every core recovered in both runs,
    the ratio of the recovered mass to the true Bonnor-Ebert mass is formed in
    each run, and the two are compared core by core.

2.  Does the recovered mass track the local fluctuation?  In the unflattened
    run the uncorrected mass exceeds the truth at the median, which would
    follow if the recovered cores sit preferentially on positive fluctuations
    of the cloud.  The truth table records, per core, the column density of
    the map averaged over a disk of twice the truncation radius,
    local_Sigma, and the column the model was built for, SD_emb.  Their ratio
    is the fluctuation the core sits on, and the test is whether the error in
    the recovered mass correlates with it.

Usage:
    python paired_injection_v3.py [--dir inj2]
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
from scipy.stats import spearmanr, wilcoxon                # noqa: E402

DIST = 260.0
MATCH_PIX = 10.0
RUNS = [('unflattened', 'cmf2.00_n600_s2029'),
        ('flattened', 'cmf2.00_n600_s2029_flat')]


def _col(t, name):
    v = np.asarray(t.col(name), float)
    return v if len(v) == t.nrows else v.reshape(t.nrows, -1)[:, 0]


def read_truth(path):
    out = []
    for ln in open(path):
        if not ln.strip() or ln.lstrip()[0] in '#!':
            continue
        f = ln.split()
        out.append(dict(x=float(f[2]), y=float(f[3]), sigma=float(f[4]),
                        sd_emb=float(f[5]), mbe=float(f[8]),
                        rbe=float(f[9]),
                        fwhm=float(f[10]) if len(f) > 11 else np.nan))
    return out


def load_run(run_dir):
    """Return {truth index: dict of recovered quantities} for one run."""
    truth = read_truth(glob.glob(os.path.join(run_dir, '*_truth.txt'))[0])
    cats = [os.path.join(run_dir, 'Aquila.s.sources.ok.cat'),
            os.path.join(run_dir, 'Aquila.s.sources.ok.add.cat')]
    table, _ = mc.correct_cloud(cats, DIST)
    t = GetsfTable(cats[0])
    xs, ys, nos = _col(t, 'XCO_P'), _col(t, 'YCO_P'), _col(t, 'NO')
    pos = {int(n): (x, y) for n, x, y in zip(nos, xs, ys)}
    tx = np.array([r['x'] for r in truth])
    ty = np.array([r['y'] for r in truth])
    best = {}
    for i, no in enumerate(table['NO']):
        no = int(no)
        if no not in pos:
            continue
        xp, yp = pos[no]
        d = np.hypot(tx - xp, ty - yp)
        j = int(np.argmin(d))
        if d[j] >= MATCH_PIX:
            continue
        if j in best and best[j]['d'] <= d[j]:
            continue
        best[j] = dict(d=d[j], mbe=truth[j]['mbe'],
                       m_init=table['M_init'][i], m3=table['M_corr'][i],
                       m4=table['M_corr_4D'][i],
                       in_dom=bool(table['IN_DOMAIN'][i]),
                       fluct=truth[j]['sigma'] / truth[j]['sd_emb'],
                       fwhm=truth[j]['fwhm'])
    return truth, best


def stats(r, name):
    g = np.isfinite(r) & (r > 0)
    lr = np.log10(r[g])
    return ('    %-24s n=%3d  median %.3f  scatter %.2f  within 2x %3.0f%%  '
            'median |log10| %.3f'
            % (name, g.sum(), np.median(r[g]), 10 ** lr.std(),
               100 * np.mean((r[g] > 0.5) & (r[g] < 2)),
               np.median(np.abs(lr))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='inj2')
    a = ap.parse_args()
    mc.VERBOSE = -1

    runs = {}
    for label, sub in RUNS:
        runs[label] = load_run(os.path.join(a.dir, sub))
    truth = runs['unflattened'][0]
    U, F = runs['unflattened'][1], runs['flattened'][1]

    common = sorted(set(U) & set(F))
    print('cores injected                     : %d' % len(truth))
    print('recovered in the unflattened run   : %d' % len(U))
    print('recovered in the flattened run     : %d' % len(F))
    print('recovered in BOTH (paired sample)  : %d' % len(common))
    print('recovered only without flattening  : %d' % len(set(U) - set(F)))
    print('recovered only with flattening     : %d' % len(set(F) - set(U)))

    if not common:
        return
    mbe = np.array([U[j]['mbe'] for j in common])
    print('\n=== Paired sample: the same %d cores in both runs ===' % len(common))
    print('  ratio of the mass obtained to the true mass M_BE')
    for label, D in (('unflattened', U), ('flattened', F)):
        print('  %s' % label)
        for nm, k in (('uncorrected', 'm_init'),
                      ('three observables', 'm3'),
                      ('four observables', 'm4')):
            r = np.array([D[j][k] for j in common]) / mbe
            print(stats(r, nm))

    print('\n  paired difference, per core, in log10 of the ratio')
    print('  %-22s %10s %10s %10s' % ('treatment', 'median', 'mean', 'p'))
    for nm, k in (('uncorrected', 'm_init'), ('three observables', 'm3'),
                  ('four observables', 'm4')):
        ru = np.array([U[j][k] for j in common]) / mbe
        rf = np.array([F[j][k] for j in common]) / mbe
        g = np.isfinite(ru) & np.isfinite(rf) & (ru > 0) & (rf > 0)
        dl = np.log10(ru[g]) - np.log10(rf[g])
        try:
            p = wilcoxon(dl).pvalue
        except Exception:
            p = np.nan
        print('  %-22s %10.3f %10.3f %10.2g'
              % (nm, np.median(dl), np.mean(dl), p))
    print('  positive means the unflattened run reports more mass than the')
    print('  flattened one for the same core; p is the Wilcoxon signed-rank')
    print('  probability that the paired differences are centred on zero.')

    # ---- does the error track the local fluctuation? -------------------
    print('\n=== Does the recovered mass track the cloud fluctuation? ===')
    print('  fluctuation = local_Sigma / SD_emb, the column of the map')
    print('  averaged over twice the truncation radius, divided by the column')
    print('  the model was built for.')
    for label, D in (('unflattened', U), ('flattened', F)):
        j = [k for k in D if np.isfinite(D[k]['fluct']) and D[k]['fluct'] > 0]
        fl = np.array([D[k]['fluct'] for k in j])
        r = np.array([D[k]['m_init'] / D[k]['mbe'] for k in j])
        g = np.isfinite(r) & (r > 0)
        rho = spearmanr(np.log10(fl[g]), np.log10(r[g]))
        print('  %-12s n=%3d  median fluctuation %.3f  '
              'rank correlation with log10(M_reported/M_BE) %+.3f (p=%.2g)'
              % (label, g.sum(), np.median(fl[g]), rho.statistic, rho.pvalue))
    allf = np.array([t['sigma'] / t['sd_emb'] for t in truth
                     if t['sd_emb'] > 0])
    recf = np.array([U[j]['fluct'] for j in U])
    print('  median fluctuation, all %d injected cores      : %.3f'
          % (len(allf), np.median(allf)))
    print('  median fluctuation, the %d recovered cores     : %.3f'
          % (len(recf), np.median(recf)))


if __name__ == '__main__':
    main()
