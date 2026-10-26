#!/usr/bin/env python3
"""cylinder_table.py -- Table tab:cyl and the numbers of Sect. 5.2 of the paper, for the seven isolated cylinders.

For each cylinder it takes the median profile across the cylinder over the middle 801 rows of the true surface density image
(mc.surfdens.fits) and of the uncorrected hires image at 13.5 arcsec (hi.surface.density.r13p5.fits), convolves the true profile
with the beam of 13.5 arcsec, and corrects the uncorrected profile with the factor of hicorrect2.py (G_0 = 10, held at its maximum),
as fig_cylinders of make_hires2_figures12.py does. It prints:

  crest   the accuracy ratio on the axis, uncorrected and corrected;
  Lambda  the accuracy ratio of the linear density within +-300 arcsec, uncorrected and corrected;
  H       the half-maximum width of the true, uncorrected and corrected profiles (arcsec);
  alpha   the slope of log Sigma against log r between 20 and 120 arcsec from the crest, for the same three profiles;
  wings   the largest relative error of the corrected profile beyond 20 arcsec from the crest.

Usage:

    python3 cylinder_table.py <directory of hicorrect2.py> <directory holding the seven cylinder directories>
"""

import os
import sys
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.abspath(sys.argv[1]))
import hicorrect2 as h

PIX = 3.0
CYLINDERS = [("cyl_sd3e21_isolated", "3e21"), ("cyl_sd1e22_isolated", "1e22"), ("cylinder_isolated", "3e22"),
             ("cyl_sd5e22_isolated", "5e22"), ("cyl_sd1e23_isolated", "1e23"), ("cyl_p2_isolated", "3e22 p=2"),
             ("cyl_w80_isolated", "3e22 H=80")]


def width(p, x0):
    """Half-maximum width of the profile p about its peak x0, in arcsec, by linear interpolation."""
    half = 0.5 * p[x0]
    i = x0
    while p[i] > half:
        i += 1
    right = i - 1 + (p[i - 1] - half) / (p[i - 1] - p[i])
    j = x0
    while p[j] > half:
        j -= 1
    left = j + 1 - (p[j + 1] - half) / (p[j + 1] - p[j])
    return (right - left) * PIX


def slope(p, x0):
    r = np.abs(np.arange(len(p)) - x0) * PIX
    m = (r >= 20.0) & (r <= 120.0)
    return -np.polyfit(np.log10(r[m]), np.log10(p[m]), 1)[0]


def main():
    fitted, factor = h.build_relation([160.0, 250.0, 350.0, 500.0])
    factor = np.maximum.accumulate(factor)
    for d, label in CYLINDERS:
        path = os.path.join(sys.argv[2], d)
        t0 = fits.getdata(os.path.join(path, "mc.surfdens.fits")).astype(float)
        u = fits.getdata(os.path.join(path, "hi.surface.density.r13p5.fits")).astype(float)
        rows = slice(u.shape[0] // 2 - 400, u.shape[0] // 2 + 401)
        pt0 = np.nanmedian(t0[rows], 0)
        pu = np.nanmedian(u[rows], 0)
        x0 = int(np.nanargmax(pt0))
        pt = gaussian_filter1d(pt0, 13.5 / 2.3548 / PIX, mode="nearest")
        pc = pu * np.interp(np.log10(np.maximum(pu, 1.0e18)), np.log10(fitted), factor)
        lin = lambda p: p[x0 - 100:x0 + 101].sum()
        rel = pc[x0 - 100:x0 + 101] / pt[x0 - 100:x0 + 101] - 1.0
        far = np.abs(np.arange(-100, 101)) * PIX > 20.0
        print(" %-10s crest %.2f -> %.2f   Lambda %.2f -> %.2f   H %3.0f %3.0f %3.0f   alpha %.2f %.2f %.2f   wings %.3f"
              % (label, pu[x0] / pt[x0], pc[x0] / pt[x0], lin(pu) / lin(pt), lin(pc) / lin(pt), width(pt, x0),
                 width(pu, int(np.nanargmax(pu))), width(pc, int(np.nanargmax(pc))), slope(pt, x0), slope(pu, x0),
                 slope(pc, x0), np.max(np.abs(rel[far]))))


if __name__ == "__main__":
    main()
