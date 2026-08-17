#!/usr/bin/env python
"""make_fig_grid_coverage3.py -- coverage of the extended production grid in
the plane of Bonnor-Ebert mass against intrinsic, unconvolved projected size.

Unlike its predecessor, this script reads the grid file written by
make_bes_grid.py directly, so no quantity is re-derived: M_BE, the intrinsic
FWHM, xi_max and the critical-mass ratio a_BE are all taken from their own
columns.  The intrinsic size is denoted H_BE, for consistency with M_BE and
T_BE, and is the FWHM the sphere would present before convolution with any
beam.

Points are the grid nodes, colored by their dimensionless truncation radius
xi_max, which measures how centrally concentrated the sphere is: xi_max below
the critical 6.451 is a sub-critical, nearly uniform sphere, above it a
concentrated, gravitationally unstable one.  Open symbols mark nodes below the
contrast floor, which are tabulated but not computed.  The grey curves join,
for each embedding column density, the nodes that are exactly critical.

Usage:
    python make_fig_grid_coverage3.py [grid file] [-o fig_grid_coverage.pdf]
"""
import argparse
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                          # noqa: E402
import matplotlib.colors as mcolors                       # noqa: E402
from matplotlib.colors import Normalize                   # noqa: E402
from matplotlib.lines import Line2D                       # noqa: E402
import aa_plotstyle                                       # noqa: E402,F401
from aa_plotstyle import clean_log, no_top_minor          # noqa: E402

PC_TO_AU = 206265.0

XI_CRIT = 6.451
DIST_PC = 260.0


def read_grid(path):
    head = open(path).readline().split()
    rows = [l.split() for l in open(path).readlines()[1:] if l.strip()]
    col = lambda n: np.array([float(r[head.index(n)]) for r in rows])
    return dict(sd=col('SD_emb'), M=col('M_BE'), H=col('FWHM_pc'),
                xi=col('xi_max'), det=col('detect') > 0, aBE=col('a_BE'),
                hole=col('hole') > 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('grid', nargs='?', default='bes_grid_final3.txt')
    ap.add_argument('-o', default='fig_grid_coverage.pdf')
    a = ap.parse_args()
    g = read_grid(a.grid)
    print('%d nodes, %d detectable, %d flagged as digging a hole'
          % (len(g['M']), int(g['det'].sum()), int(g['hole'].sum())))

    fig, ax = plt.subplots(figsize=(4.6, 4.2))

    # the top few per cent of plasma is a pale yellow that reads poorly on a
    # white background, so the colormap is truncated below it
    cmap = mcolors.LinearSegmentedColormap.from_list(
        'plasma_dark', plt.cm.plasma(np.linspace(0.0, 0.88, 256)))
    xi = g['xi']
    norm = Normalize(vmin=np.log10(xi.min()), vmax=np.log10(xi.max()))

    d = g['det']
    ax.scatter(g['H'][~d], g['M'][~d], s=10, facecolors='none',
               edgecolors='0.72', linewidths=0.4, rasterized=False, zorder=2)
    sc = ax.scatter(g['H'][d], g['M'][d], c=np.log10(xi[d]), cmap=cmap,
                    norm=norm, s=10, alpha=0.85, linewidths=0,
                    rasterized=False, zorder=3)

##    sc = ax.scatter(g['H'][d], g['M'][d], c=np.log10(xi[d]), cmap=_plasma_dark, norm=norm,
##                    s=10, alpha=0.85, linewidths=0, zorder=2)

    # The critical curves are drawn in a grey ramp, dark for the lowest cloud
    # column and light for the highest.  Grey separates them from the plasma
    # colormap at both of its ends, which red does not, and the ramp still
    # distinguishes the individual curves from one another.
    SD = np.unique(g['sd'])
    greys = plt.cm.Greys(np.linspace(0.95, 0.45, len(SD)))
    for c, s_ in zip(greys, SD):
        m = g['sd'] == s_
        H, M, A = g['H'][m], g['M'][m], g['aBE'][m]
        xs, ys = [], []
        for h in np.unique(H):
            k = H == h
            if k.sum() < 2:
                continue
            j = np.argmin(np.abs(np.log(A[k])))
            if abs(np.log10(A[k][j])) < 0.3:
                xs.append(h); ys.append(M[k][j])
        if len(xs) > 2:
            o = np.argsort(xs)
            ax.plot(np.array(xs)[o], np.array(ys)[o], color=c, lw=1.3,
                    zorder=4)
    leg = Line2D([0], [0], color=greys[len(greys) // 2], lw=1.3)
    ax.legend([leg],
              [r'critical ($\xi_{\max}=6.451$), per $\Sigma_{\rm cloud}$'],
              loc='upper left', fontsize=7.5, frameon=False)

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'$H_{\rm BE}$ (pc)')
    ax.set_ylabel(r'$M_{\rm BE}$ ($M_\odot$)')
    ax.set_xlim(2e-3, 0.5)
    ax.set_ylim(1e-3, 60)
    clean_log(ax)

    secax = ax.secondary_xaxis('top', functions=(lambda x: x * PC_TO_AU,
                                                 lambda x: x / PC_TO_AU))
    secax.set_xlabel(r'$H_{\rm BE}$ (AU)')
    # the secondary axis carries its own scale in AU, so minor ticks from the
    # primary pc axis must not be drawn on top of it
    no_top_minor(ax)

    for _a in (ax, secax):
        _a.tick_params(which='major', direction='in', length=5.0, width=0.8)
        _a.tick_params(which='minor', direction='in', length=3.5, width=0.6)
    ax.tick_params(which='both', top=False, right=True)

    cb = fig.colorbar(sc, ax=ax, pad=0.01)
    cb.set_label(r'$\log_{10}\,\xi_{\max}$')
    cb.ax.tick_params(which='both', direction='in', length=4.0)
    
    fig.savefig(a.o, bbox_inches='tight', pad_inches=0)
    print('wrote %s' % a.o)

if __name__ == '__main__':
    main()
