#!/usr/bin/env python3
"""simsky_numbers.py -- the numbers of Sect. 5.3 of the paper, for the simulated sky simsky6 with the enhanced RADMC-3D cloud.

The images are corrected with the factor of hicorrect2.py (G_0 = 10, held at its maximum), as fig_simsky6 of
make_hires2_figures12.py does, and compared with the true surface density at 13.5 arcsec. It prints:

  1. for the sky that contains the cloud alone: the median accuracy ratio of the corrected and of the uncorrected surface density
     over all pixels, and the median corrected ratio within +-0.1 dex of 10^21 and 10^22 cm^-2 of true surface density;
  2. for the full sky, with the cores masked (a square of half side max(2.5 R_BE, 30 arcsec) plus three pixels about each core):
     the median accuracy ratio over the field and the ratio of the sums;
  3. the relative error of the corrected surface density across each of the seven filaments, as the median over the outer rows
     (30-300 and n-300 to n-30), where the cloud is weak, and over the central rows (600-1200), where it is dense, at offsets of
     0, 15, 30, 45, 60, 90, 150 and 240 arcsec from the crest (both sides averaged), with the minimum and maximum of the profile
     within +-300 arcsec, and the same for the cloud alone at the same positions.

Usage:

    python3 simsky_numbers.py <directory of hicorrect2.py> <sky directory> <cloud-only directory> <core catalog>

The sky directory holds simsky6.w40.surfdens.r13p5.cores.fits and hi.surface.density.r13p5.fits, the cloud-only directory
simsky6.w40.surfdens.r13p5.fits and hi.surface.density.r13p5.fits, and the core catalog is simsky6.w40.cores.cores.txt.
"""

import os
import sys
import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.abspath(sys.argv[1]))
import hicorrect2 as h

PIX = 3.0
COLS = [224, 448, 672, 896, 1120, 1344, 1568]


def main():
    sky, cld, catalog = sys.argv[2], sys.argv[3], sys.argv[4]
    g = lambda d, f: fits.getdata(os.path.join(d, f)).astype(float)
    t = g(sky, "simsky6.w40.surfdens.r13p5.cores.fits")
    u = g(sky, "hi.surface.density.r13p5.fits")
    tc = g(cld, "simsky6.w40.surfdens.r13p5.fits")
    uc = g(cld, "hi.surface.density.r13p5.fits")
    fitted, factor = h.build_relation([160.0, 250.0, 350.0, 500.0])
    factor = np.maximum.accumulate(factor)
    cor = lambda x: x * np.interp(np.log10(np.maximum(x, 1.0e18)), np.log10(fitted), factor)
    c, cc = cor(u), cor(uc)

    print(" cloud alone: median corrected / true %.3f, uncorrected / true %.3f" % (np.median(cc / tc), np.median(uc / tc)))
    for s0 in (1.0e21, 1.0e22):
        m = np.abs(np.log10(tc / s0)) < 0.1
        print("   at %.0e cm^-2: median corrected / true %.3f" % (s0, np.median(cc[m] / tc[m])))

    mask = np.zeros(t.shape, bool)
    for line in open(catalog):
        if line.startswith("#") or not line.strip():
            continue
        w = line.split()
        y, x, rbe = int(w[1]), int(w[2]), float(w[4])
        r = int(max(2.5 * rbe, 30.0) / PIX) + 3
        mask[max(0, y - r):y + r + 1, max(0, x - r):x + r + 1] = True
    ok = ~mask
    print(" full sky, cores masked: median corrected / true %.3f, ratio of sums %.3f"
          % (np.median(c[ok] / t[ok]), c[ok].sum() / t[ok].sum()))

    rel = np.where(mask, np.nan, c / t - 1.0)
    relc = cc / tc - 1.0
    n = t.shape[0]
    off = np.arange(-100, 101)
    steps = [0, 5, 10, 15, 20, 30, 50, 80]
    for rows, title in ((np.r_[30:300, n - 300:n - 30], "outer rows"), (np.r_[600:1200], "central rows")):
        print(" %s, error at offsets of 0, 15, 30, 45, 60, 90, 150, 240 arcsec, then min and max" % title)
        for col in COLS:
            p = np.nanmedian(rel[rows][:, col + off], axis=0)
            v = [0.5 * (p[100 + s] + p[100 - s]) if s else p[100] for s in steps]
            print("   column %4d  " % col + " ".join("%+.3f" % x for x in v) + "   %+.3f %+.3f" % (np.nanmin(p), np.nanmax(p)))
        pc = np.nanmedian(np.array([np.median(relc[rows][:, col + off], axis=0) for col in COLS]), axis=0)
        print("   cloud alone  crest %+.3f, edge %+.3f" % (pc[100], pc[0]))


if __name__ == "__main__":
    main()
