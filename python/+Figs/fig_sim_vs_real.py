#!/usr/bin/env python
"""fig_sim_vs_real.py -- compare the original simulated Aquila population with
the real Aquila sources, in the two planes that decide whether an injected
population resembles the observed one.

Left panel   reported SED mass against measured size AFWHM.
Right panel  contrast PEAK^SBF / PEAK^BGF against background column density.

Both quantities in both panels are measured by getsf, on the simulated field
and on the observed field respectively, so the comparison is like for like and
requires no model-side assumption.  Only sources passing the same quality and
observable cuts are shown.

The simulated sources are the 1545 injected cores of the original simulation
that were recovered and matched to their truth-table entry; the real sources
are the 1030 Aquila sources of the reference catalogue.  For reference, the
population-matched injection built with the placement method described in the
paper is overlaid as contours where its file is available.

Usage:
    python fig_sim_vs_real.py [-o fig_sim_vs_real.pdf]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                          # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '.')
sys.path.insert(0, 'python')
from analyse_aquila_sim import read_truth, col           # noqa: E402
from getsf_columns import GetsfTable                     # noqa: E402
import mass_correction_pipeline_6 as mc                  # noqa: E402
try:
    from aa_plotstyle import clean_log
except Exception:
    def clean_log(ax):
        return


def simulated(dirname='asim'):
    """Recovered simulated sources: size, mass, background column, contrast."""
    T = read_truth(os.path.join(dirname, 'truth.table.0250.dat'))
    tx, ty = T[:, 0], T[:, 1]
    cats = [os.path.join(dirname, 'aquila-sim.s.sources.ok.cat'),
            os.path.join(dirname, 'aquila-sim.s.sources.ok.add.cat')]
    t = GetsfTable(cats[0])
    xs, ys, nos = col(t, 'XCO_P'), col(t, 'YCO_P'), col(t, 'NO')
    afw = col(t, 'AFWHM03')          # 250 um band in this extraction
    mc.load_getsf._dist = 260.0
    g = mc.load_getsf(cats)
    ok = mc.source_mask(g) & (g['Nbg'] > 0) & (g['peak'] > 0) & (g['mass'] > 0)
    idx = {int(n): k for k, n in enumerate(nos)}
    best = {}
    for k in np.where(ok)[0]:
        j0 = idx.get(int(g['NO'][k]))
        if j0 is None:
            continue
        d = np.hypot(tx - xs[j0], ty - ys[j0])
        j = int(np.argmin(d))
        if d[j] >= 10 or (j in best and best[j][0] <= d[j]):
            continue
        best[j] = (d[j], k, j0)
    A, M, S, C = [], [], [], []
    for j, (d, k, j0) in best.items():
        a = afw[j0]
        if not (np.isfinite(a) and 0 < a < 1e3):
            continue
        A.append(a); M.append(g['mass'][k]); S.append(g['Nbg'][k])
        C.append(1.0 + g['peak'][k] / g['Nbg'][k])
    return map(np.array, (A, M, S, C))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', default='fig_sim_vs_real.pdf')
    ap.add_argument('--real', default='python/joint_target_Aquila_v3.txt')
    ap.add_argument('--matched', default='dry_r95.txt',
                    help='placement table of the population-matched '
                         'injection, drawn as contours if present')
    a = ap.parse_args()
    mc.VERBOSE = -1

    As, Ms, Ss, Cs = simulated()
    d = np.loadtxt(a.real)
    Mr, Sr, Cr, Ar = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    print('simulated sources plotted: %d ; real sources: %d' % (len(As), len(Ar)))

    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3))

    ax[0].scatter(Ar, Mr, s=4, c='0.55', alpha=0.35, lw=0,
                  label='real Aquila, %d' % len(Ar))
    ax[0].scatter(As, Ms, s=5, c='crimson', alpha=0.45, lw=0,
                  label='simulated, %d' % len(As))
    ax[0].set_xscale('log'); ax[0].set_yscale('log')
    ax[0].set_xlabel(r'AFWHM  (arcsec)', fontsize=8)
    ax[0].set_ylabel(r'reported SED mass  ($M_\odot$)', fontsize=8)

    ax[1].scatter(Sr, Cr - 1, s=4, c='0.55', alpha=0.35, lw=0)
    ax[1].scatter(Ss, Cs - 1, s=5, c='crimson', alpha=0.45, lw=0)
    ax[1].set_xscale('log'); ax[1].set_yscale('log')
    ax[1].set_xlabel(r'$\Sigma_{\rm cloud}$  (cm$^{-2}$)', fontsize=8)
    ax[1].set_ylabel(r'contrast $-$ 1  =  PEAK$^{\rm SRC}$/PEAK$^{\rm BGF}$',
                     fontsize=8)

    if os.path.exists(a.matched):
        # the first column of the placement table is the node tag, so the
        # numeric columns are read by name-order after dropping it:
        # 1 x_pix, 2 y_pix, 3 SD_emb, 4 local_Sigma, 5 contrast, 6 M_BE,
        # 7 M_SED, 8 FWHM_SD, 9 R_BE_as
        m = np.loadtxt(a.matched, usecols=range(1, 10))
        m_sig, m_con, m_sed, m_fw = m[:, 3], m[:, 4], m[:, 6], m[:, 7]
        ax[1].scatter(m_sig, m_con - 1, s=4, c='steelblue', alpha=0.5, lw=0)
        ax[0].scatter(m_fw, m_sed, s=4, c='steelblue', alpha=0.5, lw=0,
                      label='population-matched, %d' % len(m))

    for x in ax:
        clean_log(x)
        x.tick_params(labelsize=7)
    ax[0].legend(fontsize=6.5, loc='upper left')
    fig.tight_layout()
    fig.savefig(a.o)
    print('wrote %s' % a.o)


if __name__ == '__main__':
    main()
