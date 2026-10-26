#!/usr/bin/env python3
"""aquila_numbers.py -- the numbers of Sect. 6 of the paper, for the Aquila field corrected with G_0 = 33.

It reads the uncorrected surface density at 13.5 arcsec, the correction factor image and the corrected surface density written by
hicorrect2 with -g0=33, and the corrected and fitted temperature images, and prints the median and the 1, 25, 75 and 99th
percentiles of the factor, its maximum and the lowest surface density where it is reached, the ratio of the total corrected to the
total uncorrected surface density, the median and 99th percentile of the surface density before and after the correction, the
numbers of pixels above 10^22, 3 x 10^22 and 10^23 cm^-2, and the median and 1st and 99th percentiles of the corrected temperature.
It then recomputes the factor from hicorrect2.py, checks it against the factor image, and gives the median factor and the total
ratio with the calibration field G_0 = 10.

Usage:

    python3 aquila_numbers.py <directory of hicorrect2.py> <directory of the images>

The directory holds hi.surface.density.r13p5.fits, hi.surface.density.r13p5.cfactor_G0=33.fits,
hi.surface.density.r13p5.correct_G0=33.fits, hi.temperature.r13p5.correct_G0=33.fits and hi.temperature.r13p5.fits.
"""

import os
import sys
import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.abspath(sys.argv[1]))
import hicorrect2 as h


def main():
    g = lambda f: fits.getdata(os.path.join(sys.argv[2], f)).astype(float)
    u = g("hi.surface.density.r13p5.fits")
    f = g("hi.surface.density.r13p5.cfactor_G0=33.fits")
    c = g("hi.surface.density.r13p5.correct_G0=33.fits")
    tc = g("hi.temperature.r13p5.correct_G0=33.fits")
    tf = g("hi.temperature.r13p5.fits")
    ok = np.isfinite(u) & np.isfinite(f) & (u > 0.0) & (f > 0.0)
    u, f, c, tc, tf = u[ok], f[ok], c[ok], tc[ok], tf[ok]
    print(" %d pixels" % ok.sum())
    print(" factor: median %.2f, percentiles 1, 25, 75, 99: %s" % (np.median(f), np.round(np.percentile(f, [1, 25, 75, 99]), 2)))
    print(" maximum %.3f, reached from %.2e cm^-2" % (f.max(), u[f >= 0.999 * f.max()].min()))
    print(" total corrected / total uncorrected %.3f" % (c.sum() / u.sum()))
    print(" median Sigma %.2e -> %.2e, 99th percentile %.2e -> %.2e"
          % (np.median(u), np.median(c), np.percentile(u, 99), np.percentile(c, 99)))
    for th in (1.0e22, 3.0e22, 1.0e23):
        print(" pixels above %.0e cm^-2: %d -> %d" % (th, (u > th).sum(), (c > th).sum()))
    print(" corrected temperature: median %.1f K (fitted %.1f K), percentiles 1, 99: %s"
          % (np.nanmedian(tc), np.nanmedian(tf), np.round(np.nanpercentile(tc, [1, 99]), 1)))

    fitted, factor = h.build_relation([160.0, 250.0, 350.0, 500.0])
    for g0 in (33.0, 10.0):
        a = np.maximum.accumulate(factor)              # the factor for G_0 = 10 held at its maximum, then times g, Eq. (9)
        if g0 != 10.0:
            a = a * h.field_factor(np.log10(fitted), g0)
        ff = np.interp(np.log10(u), np.log10(fitted), a)
        print(" recomputed, G_0 = %g: median factor %.3f, total ratio %.3f" % (g0, np.median(ff), (u * ff).sum() / u.sum()))
        if g0 == 33.0:
            print("   largest difference from the factor image %.4f" % np.max(np.abs(ff / f - 1.0)))


if __name__ == "__main__":
    main()
