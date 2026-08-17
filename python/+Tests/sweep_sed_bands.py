#!/usr/bin/env python
"""sweep_sed_bands.py -- choose the band-selection criteria of the spectral
energy distribution fit by measuring them against known truth.

Two settings govern which bands are fitted:

    SNR_MIN_PEAK, SNR_MIN_TOTAL  a band is admitted only if its peak and total
                                 fluxes exceed their own uncertainties by these
                                 factors
    MIN_SED_BANDS                the fit is attempted only if at least this
                                 many bands are admitted

Neither can be settled from first principles.  A stricter threshold rejects
poorly measured bands but also discards sources; a looser one keeps sources but
admits fluxes that carry no information.  On the population-matched injections
the true mass of every recovered core is known, so the trade-off can simply be
measured: for each setting the script reports how many cores survive, how
accurately their masses are recovered, and how many are catastrophically wrong.

The same sweep is run on the real Aquila catalogue, where the truth is not
known, to show how many real sources each setting would retain.

Usage:
    python sweep_sed_bands.py [--dir inj3] [--seeds 2030 2031 2032 2033]
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

MATCH_PIX = 10.0
DIST = 260.0


def _col(t, name):
    v = np.asarray(t.col(name), float)
    return v if len(v) == t.nrows else v.reshape(t.nrows, -1)[:, 0]


def truth_of(run_dir):
    out = []
    for ln in open(glob.glob(os.path.join(run_dir, '*_truth.txt'))[0]):
        if not ln.strip() or ln.lstrip()[0] in '#!':
            continue
        f = ln.split()
        out.append((float(f[2]), float(f[3]), float(f[8])))
    return np.array(out)


def measure(run_dirs):
    """Reported mass and true mass of every matched core, at the current
    module-level thresholds.  load_getsf is re-run so the thresholds bite."""
    mt, mx = [], []
    for d in run_dirs:
        T = truth_of(d)
        cats = [os.path.join(d, 'Aquila.s.sources.ok.cat'),
                os.path.join(d, 'Aquila.s.sources.ok.add.cat')]
        t = GetsfTable(cats[0])
        xs, ys, nos = _col(t, 'XCO_P'), _col(t, 'YCO_P'), _col(t, 'NO')
        mc.load_getsf._dist = DIST
        g = mc.load_getsf(cats)
        ok = mc.source_mask(g) & np.isfinite(g['mass']) & (g['mass'] > 0)
        idx = {int(n): k for k, n in enumerate(nos)}
        best = {}
        for k in np.where(ok)[0]:
            j0 = idx.get(int(g['NO'][k]))
            if j0 is None:
                continue
            dd = np.hypot(T[:, 0] - xs[j0], T[:, 1] - ys[j0])
            j = int(np.argmin(dd))
            if dd[j] >= MATCH_PIX or (j in best and best[j][0] <= dd[j]):
                continue
            best[j] = (dd[j], g['mass'][k])
        for j, (dd, m) in best.items():
            mt.append(T[j, 2]); mx.append(m)
    return np.array(mt), np.array(mx)


def real_count(cats):
    mc.load_getsf._dist = DIST
    g = mc.load_getsf(cats)
    return int((mc.source_mask(g) & np.isfinite(g['mass']) & (g['mass'] > 0)).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='inj3')
    ap.add_argument('--seeds', nargs='+',
                    default=['2030', '2031', '2032', '2033'])
    ap.add_argument('--variant', default='_flat')
    ap.add_argument('--real',
                    default='cats/Aquila.s.sources.ok.cat='
                            'Aquila.s.sources.ok.add.cat')
    a = ap.parse_args()
    mc.VERBOSE = -1

    runs = [os.path.join(a.dir, 'cmf2p00_n600_s%s%s' % (s, a.variant))
            for s in a.seeds]
    runs = [r for r in runs if os.path.isdir(r)]
    print('injection runs used: %d (%s)' % (len(runs), a.variant or 'as placed'))

    print('\nRecovery of the true mass, over the injected cores that survive '
          'each setting.')
    print('Ratio is the reported mass divided by the true Bonnor-Ebert mass; '
          'scatter is the\nstandard deviation of its base-10 logarithm, as a '
          'multiplicative factor;\n"wrong by >5x" counts cores whose ratio '
          'lies outside 0.2 to 5.\n')
    print('  %-6s %-6s %7s %9s %9s %11s %11s %9s'
          % ('S/N', 'bands', 'cores', 'median', 'scatter', 'within 2x',
             'wrong by>5x', 'real'))
    for snr in (1.0, 1.5, 2.0, 3.0):
        for nb in (2, 3, 4):
            mc.SNR_MIN_PEAK = snr
            mc.SNR_MIN_TOTAL = snr
            mc.MIN_SED_BANDS = nb
            mt, mx = measure(runs)
            nreal = real_count([a.real]) if os.path.exists(a.real) else -1
            if len(mt) == 0:
                continue
            r = mx / mt
            g = np.isfinite(r) & (r > 0)
            lr = np.log10(r[g])
            n5 = int((np.abs(lr) > np.log10(5)).sum())
            print('  %-6.1f %-6d %7d %9.3f %9.2f %10.0f%% %6d (%3.0f%%) %9s'
                  % (snr, nb, g.sum(), np.median(r[g]), 10 ** lr.std(),
                     100 * np.mean((r[g] > 0.5) & (r[g] < 2)),
                     n5, 100.0 * n5 / g.sum(),
                     nreal if nreal >= 0 else '-'))


if __name__ == '__main__':
    main()
