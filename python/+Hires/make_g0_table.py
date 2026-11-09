#!/usr/bin/env python3
"""make_g0_table.py -- the table of model temperatures with which hicorrect2 estimates the radiation field of a cloud.

For each of the nine fields of the grid, the median fitted temperature of the models is computed in bins of 0.2 dex of the fitted
surface density, from 10^21.0 to 10^23.2 cm^-2 (no model is denser), exactly as g0_requirement.py does it: the median over the cells of the table
harvest_calib3_G0=<field>.txt whose surface density lies in the bin, weighted by their numbers of pixels. The table is printed in
the syntax of hicorrect2.py (TTAB, NaN where a bin has no models) and of hicorrect2.for (ttab, -1 where a bin has no models), to be
pasted into the two codes. It needs to be recomputed only when the grid of radiation fields is recomputed.

Usage:

    python3 make_g0_table.py <directory of the harvest tables>
"""

import glob
import os
import sys

import numpy as np

STEP, XMIN, XMAX = 0.2, 21.0, 23.2


def read_harvest(path):
    """Columns x (log10 fitted surface density), y (fitted temperature, K) and the number of pixels of each cell."""
    a = []
    for line in open(path):
        c = line.split()
        if not c or c[0].startswith("!"):
            continue
        a.append((float(c[0]), float(c[1]), float(c[3])))
    a = np.array(a)
    return a[:, 0], a[:, 1], a[:, 2]


def wmedian(v, w):
    """Median of the values v with the weights w, as in g0_requirement.py."""
    k = np.argsort(v)
    v, w = v[k], w[k]
    c = np.cumsum(w)
    c = (c - 0.5 * w) / c[-1]
    return np.interp(0.5, c, v)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    lad = []
    for p in glob.glob(os.path.join(sys.argv[1], "harvest_calib3_G0=*.txt")):
        tag = os.path.basename(p)[len("harvest_calib3_G0="):-4]
        lad.append((float(tag.replace("p", ".")), read_harvest(p)))
    lad.sort(key=lambda t: t[0])
    edges = np.arange(XMIN, XMAX + 0.5 * STEP, STEP)
    med = np.full((len(edges) - 1, len(lad)), np.nan)
    for k in range(len(edges) - 1):
        for j, (g, (hx, hy, hn)) in enumerate(lad):
            m = (hx >= edges[k]) & (hx < edges[k + 1]) & (hn > 0)
            if m.any():
                med[k, j] = wmedian(hy[m], hn[m])

    gs = [g for g, _ in lad]
    print("# fields of the grid: " + ", ".join("%g" % g for g in gs))
    print("\n# Python (hicorrect2.py):")
    print("TTAB = np.array([")
    for k in range(len(edges) - 1):
        print("    [" + ", ".join("%6.3f" % v if np.isfinite(v) else "np.nan" for v in med[k])
              + "],  # %.1f-%.1f" % (edges[k], edges[k + 1]))
    print("])")
    print("\n! Fortran (hicorrect2.for), stored by bins (the first index), as gtab:")
    print("       data ttab /")
    vals = ["%.3fd0" % v if np.isfinite(v) else "-1.0d0" for v in med.T.ravel()]
    for j in range(len(gs)):
        print("! G_0 = %g" % gs[j])
        row = vals[j * med.shape[0]:(j + 1) * med.shape[0]]
        for i in range(0, len(row), 8):
            last = (j == len(gs) - 1) and (i + 8 >= len(row))
            print("     &             " + (", " if (j or i) else "  ") + ", ".join(row[i:i + 8]) + (" /" if last else ""))


if __name__ == "__main__":
    main()
