#!/usr/bin/env python3
"""field_factor_errors.py -- the errors of the correction on the grid of radiation fields (Table tab:fieldR of the paper).

For each of the nine radiation fields it reads the histogram written by harvest_calib3.py (harvest_calib3_G0=*.txt.hist.npz), which
holds, for every pixel of every model of that field, the fitted surface density at 13.5 arcsec, x = log10 Sigma_F, and the factor
that restores the true surface density, f = Sigma_true / Sigma_F (with two more observables that are summed over here). The factor
that the correction applies, A(x), is the factor of hicorrect2.py for G_0 = 10 held at its maximum, times the field factor g, in
three variants: without g, with g of all nine fields (GTAB of hicorrect2.py), and with g interpolated in log G_0 between the two
neighboring fields, the field tested being left out. The relative error of the corrected surface density of a pixel is A(x) / f - 1. The table gives its
median over the pixels in four ranges of x.

The histograms are binned (0.05 dex in x, about 1.4 percent in f), so the medians are accurate to about 1 percent.

Usage:

    python3 field_factor_errors.py <directory of hicorrect2.py> <directory of the histograms>
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(sys.argv[1]))
import hicorrect2 as h

TAGS = ["0p98", "2p12", "4p61", "10", "21p7", "47p1", "102p4", "222p2", "482p5"]
RANGES = [(21.6, 22.0), (22.0, 22.4), (22.4, 22.8), (22.8, 23.4)]


def load(path):
    """The histogram summed over the two other observables: counts against x and f, with the bin edges."""
    d = np.load(path)
    shape = (len(d["xedge"]) - 1, len(d["yedge"]) - 1, len(d["zedge"]) - 1, len(d["fedge"]) - 1)
    ix, iy, iz, jf = np.unravel_index(d["index"], shape)
    hist = np.zeros((shape[0], shape[3]))
    np.add.at(hist, (ix, jf), d["count"])
    return hist, d["xedge"], d["fedge"]


def leftout_factor(logs, k):
    """Field factor of the field k interpolated in log G_0 between its two neighbors, without its own column of GTAB."""
    lg = np.log10(h.GFIELD)
    a = (lg[k] - lg[k - 1]) / (lg[k + 1] - lg[k - 1])
    xs = h.GX0 + h.GDX * np.arange(h.GTAB.shape[0])
    lr = (1.0 - a) * np.log(h.GTAB[:, k - 1]) + a * np.log(h.GTAB[:, k + 1])
    return np.exp(np.interp(logs, xs, lr))


def medians(hist, xe, fe, applied):
    """Median relative error of the corrected surface density in each range of x."""
    xc = 0.5 * (xe[1:] + xe[:-1])
    fc = np.sqrt(fe[1:] * fe[:-1])
    a = applied(xc)
    out = []
    for lo, hi in RANGES:
        m = (xc >= lo - 1e-9) & (xc < hi - 1e-9)
        err = (a[m][:, None] / fc[None, :] - 1.0).ravel()
        w = hist[m].ravel()
        o = np.argsort(err)
        c = np.cumsum(w[o])
        out.append(100.0 * err[o][np.searchsorted(c, 0.5 * c[-1])])
    return out


def main():
    hdir = sys.argv[2]
    fitted, factor = h.build_relation([160.0, 250.0, 350.0, 500.0])
    lf = np.log10(fitted)
    print(" G_0      without g             with g                left out")
    for k, (tag, g0) in enumerate(zip(TAGS, h.GFIELD)):
        hist, xe, fe = load(os.path.join(hdir, "harvest_calib3_G0=%s.txt.hist.npz" % tag))
        rows = []
        for variant in ("without", "with", "leftout"):
            if variant == "leftout" and (k == 0 or k == len(TAGS) - 1 or g0 == 10.0):
                rows.append("        --         ")
                continue
            g = {"without": np.ones_like(lf), "with": h.field_factor(lf, g0),
                 "leftout": leftout_factor(lf, k) if variant == "leftout" and 0 < k < len(TAGS) - 1 else None}[variant]
            f = np.maximum.accumulate(factor) * g
            rows.append(" ".join("%+4.0f" % v for v in medians(hist, xe, fe, lambda x: np.interp(x, lf, f))))
        print(" %6.2f  " % g0 + "   ".join(rows))


if __name__ == "__main__":
    main()
