#!/usr/bin/env python
"""analyse_aquila_sim.py -- analyse the original simulated Aquila field, the
one whose injected population motivated this work.

The simulation placed 5622 model cores in the Aquila background and the field
was extracted with getsf, exactly as the observations were.  Two questions are
asked, and they are separate.

1.  WAS THE INJECTED POPULATION LIKE THE REAL ONE?  The claim to be tested is
    that the injected cores were too bright to be consistent with the sources
    visible in the observed images.  The test compares the contrast of the
    recovered simulated sources, PEAK^SBF / PEAK^BGF, with that of the real
    Aquila sources, both measured by getsf in the same way, so the comparison
    is like for like and needs no model-side assumption.

2.  WHAT DOES THE SELECTION FUNCTION LOOK LIKE ON A LARGE SAMPLE?  With 5622
    injected cores against the 1086 core instances of the population-matched
    injections, the completeness can be measured far more precisely, albeit
    for a brighter population.  Each injected core is matched to the nearest
    extracted source, and the recovered fraction is tabulated against the
    injected mass, size and central column density.

Note the waveband layout of this extraction: 070, 160, 165, 250, 350, 500
micron, so the column-density band is 03, not 02 as in G. Zhang's whole-field
Aquila catalogue.

Usage:
    python analyse_aquila_sim.py [--dir asim] [--real python/joint_target_Aquila_v3.txt]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '.')
sys.path.insert(0, 'python')
from getsf_columns import GetsfTable                    # noqa: E402
import mass_correction_pipeline_6 as mc                 # noqa: E402
from scipy.stats import ks_2samp, spearmanr             # noqa: E402

MATCH_PIX = 10.0
SIM_CAT = 'aquila-sim.s.sources.ok.cat'
SIM_ADD = 'aquila-sim.s.sources.ok.add.cat'
TRUTH = 'truth.table.0250.dat'


def col(t, name):
    v = np.asarray(t.col(name), float)
    return v if len(v) == t.nrows else v.reshape(t.nrows, -1)[:, 0]


def read_truth(path):
    """XPIX, YPIX, RAD_AS, AFWHM, MASS, C_COLUMN, PKMJYSR, BACKSTD, SNRATIO."""
    out = []
    for ln in open(path):
        if not ln.strip() or ln.lstrip()[0] in '#!':
            continue
        f = ln.split()
        out.append((float(f[3]), float(f[4]), float(f[7]), float(f[10]),
                    float(f[11]), float(f[15]), float(f[19]), float(f[20]),
                    float(f[21])))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='asim')
    ap.add_argument('--real', default='python/joint_target_Aquila_v3.txt')
    a = ap.parse_args()
    mc.VERBOSE = -1

    T = read_truth(os.path.join(a.dir, TRUTH))
    tx, ty = T[:, 0], T[:, 1]
    print('simulated cores injected: %d' % len(T))
    print('  injected mass (M_sun)          : 5th %.3f, median %.3f, 95th %.3f'
          % tuple(np.percentile(T[:, 4], [5, 50, 95])))
    print('  injected AFWHM (arcsec)        : 5th %.1f, median %.1f, 95th %.1f'
          % tuple(np.percentile(T[:, 3], [5, 50, 95])))
    print('  injected central column (cm^-2): 5th %.2e, median %.2e, 95th %.2e'
          % tuple(np.percentile(T[:, 5], [5, 50, 95])))

    cats = [os.path.join(a.dir, SIM_CAT), os.path.join(a.dir, SIM_ADD)]
    t = GetsfTable(cats[0])
    xs, ys, nos = col(t, 'XCO_P'), col(t, 'YCO_P'), col(t, 'NO')
    mc.load_getsf._dist = 260.0
    g = mc.load_getsf(cats)
    ok = mc.source_mask(g) & (g['Nbg'] > 0) & (g['peak'] > 0)
    print('\nextracted sources: %d, of which %d pass the quality and '
          'observable cuts' % (t.nrows, int(ok.sum())))

    # match every extracted source to the nearest injected core
    idx = {int(n): k for k, n in enumerate(nos)}
    Cs, Ss, Ms, matched = [], [], [], np.zeros(len(T), bool)
    best = {}
    for k in np.where(ok)[0]:
        no = int(g['NO'][k])
        j0 = idx.get(no)
        if j0 is None:
            continue
        d = np.hypot(tx - xs[j0], ty - ys[j0])
        j = int(np.argmin(d))
        if d[j] >= MATCH_PIX:
            continue
        if j in best and best[j][0] <= d[j]:
            continue
        best[j] = (d[j], k)
    for j, (d, k) in best.items():
        matched[j] = True
        Cs.append(1.0 + g['peak'][k] / g['Nbg'][k])
        Ss.append(g['Nbg'][k])
        Ms.append(g['mass'][k])
    Cs, Ss, Ms = np.array(Cs), np.array(Ss), np.array(Ms)
    print('injected cores recovered and matched: %d of %d (%.0f%%)'
          % (matched.sum(), len(T), 100.0 * matched.mean()))

    # ---- 1. was the population like the real one? ----------------------
    d = np.loadtxt(a.real)
    Mr, Sr, Cr, Ar = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    print('\nRecovered simulated sources against the real Aquila sources')
    print('  %-34s %11s %11s %14s'
          % ('quantity (median)', 'simulated', 'real', 'KS / probability'))
    for nm, vi, vr in (('contrast PEAK^SBF/PEAK^BGF', Cs, Cr),
                       ('background column (cm^-2)', Ss, Sr),
                       ('reported SED mass (M_sun)', Ms, Mr)):
        gg = np.isfinite(vi) & (vi > 0)
        ks = ks_2samp(np.log10(vi[gg]), np.log10(vr))
        print('  %-34s %11.4g %11.4g %8.3f / %.1e'
              % (nm, np.median(vi[gg]), np.median(vr), ks.statistic, ks.pvalue))
    gg = np.isfinite(Cs) & (Cs > 1)
    print('  fraction with contrast above 2: simulated %.0f%%, real %.0f%%'
          % (100 * np.mean(Cs[gg] > 2), 100 * np.mean(Cr > 2)))
    print('  rank correlation of contrast on column: simulated %+.3f, '
          'real %+.3f'
          % (spearmanr(np.log10(Ss[gg]), np.log10(Cs[gg])).statistic,
             spearmanr(np.log10(Sr), np.log10(Cr)).statistic))

    # ---- 2. selection function on a large sample -----------------------
    print('\nCompleteness of the extraction, against injected properties')
    for nm, v, edges, fmt in (
            ('injected mass (M_sun)', T[:, 4],
             [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 100.0], '%7.3g'),
            ('injected AFWHM (arcsec)', T[:, 3],
             [0, 20, 30, 40, 60, 100, 1e4], '%7.0f'),
            ('injected central column (cm^-2)', T[:, 5],
             [1e21, 3e21, 1e22, 3e22, 1e23, 1e25], '%7.1e'),
            ('signal-to-noise of the injected peak', T[:, 8],
             [0, 0.5, 1, 2, 5, 10, 1e4], '%7.2f')):
        print('  %s' % nm)
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (v >= lo) & (v < hi)
            if m.sum() >= 20:
                print(('    ' + fmt + ' - ' + fmt + ' %7d injected %7d '
                       'recovered %7.0f%%')
                      % (lo, hi, m.sum(), matched[m].sum(),
                         100 * matched[m].mean()))


if __name__ == '__main__':
    main()
