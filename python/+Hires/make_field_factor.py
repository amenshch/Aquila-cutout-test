#!/usr/bin/env python3
"""make_field_factor.py -- the field factor g(Sigma_F, G_0) of hicorrect2 (GTAB), from the grid of radiation fields.

For each of the nine radiation fields it reads the histogram written by harvest_calib3.py (harvest_calib3_G0=*.txt.hist.npz), which
holds for every pixel of every model the fitted surface density at 13.5 arcsec, x = log10 Sigma_F, and the factor that restores the
true surface density, f = Sigma_true / Sigma_F.

The pixels are binned in the true surface density, t = x + log10 f, in bins of 0.05 dex. Binning in x instead fails above about
10^22 cm^-2: most pixels of a model lie on the plateau of its uniform embedding cloud, the embedding surface densities step by a
factor of 2, and each field shifts these plateaus to a different x, so that neighboring bins of x alternate between hundreds of
thousands of pixels and a few hundred, and two fields compared in one bin of x hold unrelated populations. In bins of t the
populations of all fields coincide, and the mean of x and the mean of ln f in each bin vary smoothly. Each bin gives one point
(mean x, mean ln f) of the factor of that field as a function of the fitted surface density. At each row of the table, a straight
line is fitted to these points with Gaussian weights of width 0.1 dex in x about the row, times the square root of the number of
pixels of the bin, and its value at the row is taken. The field factor is ln g = ln f(G_0) - ln f(10) at every row. Up to
10^22.85 cm^-2 every field has enough pixels; above it g is held at its value there.

It prints the table in the form of hicorrect2.py (GTAB, rows from GX0 = 20.95 in steps of GDX = 0.1) and in Fortran DATA form.

Usage:

    python3 make_field_factor.py <directory of the histograms>
"""

import os
import sys
import numpy as np

TAGS = ["0p98", "2p12", "4p61", "10", "21p7", "47p1", "102p4", "222p2", "482p5"]
GX0, GDX, NROW = 20.95, 0.1, 24
XHOLD = 22.85          # the last row at which g is measured; rows above it take its value
WIDTH = 0.1            # dex, the width of the Gaussian window of the local linear fit
TSTEP = 0.05           # dex, the bins of the true surface density
NMIN = 200             # bins of t with fewer pixels are not used


def points(path):
    """Mean x, mean ln f and number of pixels in bins of the true surface density t = x + log10 f."""
    d = np.load(path)
    shape = (len(d["xedge"]) - 1, len(d["yedge"]) - 1, len(d["zedge"]) - 1, len(d["fedge"]) - 1)
    ix, iy, iz, jf = np.unravel_index(d["index"], shape)
    x = 0.5 * (d["xedge"][ix] + d["xedge"][ix + 1])
    lnf = np.log(np.sqrt(d["fedge"][jf] * d["fedge"][jf + 1]))
    w = d["count"].astype(float)
    t = x + lnf / np.log(10.0)
    edges = np.arange(np.floor(t.min() / TSTEP) * TSTEP, t.max() + TSTEP, TSTEP)
    k = np.clip(np.digitize(t, edges) - 1, 0, len(edges) - 2)
    n = np.bincount(k, w, len(edges) - 1)
    mx = np.bincount(k, w * x, len(edges) - 1) / np.maximum(n, 1.0)
    mf = np.bincount(k, w * lnf, len(edges) - 1) / np.maximum(n, 1.0)
    ok = n >= NMIN
    return mx[ok], mf[ok], n[ok]


def smooth(px, pf, pn, r):
    """ln f at x = r, from a local linear fit to the points with Gaussian weights about r."""
    kw = np.sqrt(pn) * np.exp(-0.5 * ((px - r) / WIDTH) ** 2)
    a = np.vstack([np.ones(len(px)), px - r]).T * np.sqrt(kw)[:, None]
    return np.linalg.lstsq(a, pf * np.sqrt(kw), rcond=None)[0][0]


def build(hdir):
    """The table of g: rows of x, columns of the fields."""
    data = [points(os.path.join(hdir, "harvest_calib3_G0=%s.txt.hist.npz" % t)) for t in TAGS]
    rows = GX0 + GDX * np.arange(NROW)
    lnf = np.array([[smooth(*p, min(r, XHOLD)) for p in data] for r in rows])
    table = np.exp(lnf - lnf[:, [TAGS.index("10")]])
    return rows, table


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    rows, table = build(sys.argv[1])
    print("GTAB = np.array([")
    for r, row in zip(rows, table):
        print("    [" + ", ".join("%.3f" % v for v in row) + "],  # %.2f" % r)
    print("])")
    print()
    print("! Fortran DATA statement for hicorrect2.for, column by column (one column per field)")
    for k, t in enumerate(TAGS):
        print("!                   G_0 = %s" % t.replace("p", "."))
        vals = ["%.3fd0" % v for v in table[:, k]]
        for j in range(0, NROW, 8):
            lead = "       data          gtab / " if (k == 0 and j == 0) else "     &             , "
            tail = " /" if (k == len(TAGS) - 1 and j + 8 >= NROW) else ""
            print(lead + ", ".join(vals[j:j + 8]) + tail)


if __name__ == "__main__":
    main()
