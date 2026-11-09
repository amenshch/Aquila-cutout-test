#!/usr/bin/env python3
"""footprint_fraction.py -- the fraction of the true core mass inside the footprint (the contour at 1% of the peak).

For every node of a grid, the core alone (nc.surfdens.bsl.r13p5x0.rs3p0as.fits, the embedding cloud subtracted) is summed inside
the footprint mask and over the whole image, and the ratio is reported: median and 16th and 84th percentiles over the nodes.

Usage: python3 footprint_fraction.py <root directory>
"""
import os
import sys

import numpy as np
from astropy.io import fits

CORE = "nc.surfdens.bsl.r13p5x0.rs3p0as.fits"
MASK = "nc.surfdens.bsl.r13p5x0.rs3p0as.mask.fits"

frac = []
for p, dirs, files in os.walk(sys.argv[1]):
    dirs[:] = [d for d in dirs if not d.startswith((".", "+"))]
    if CORE in files and MASK in files:
        core = np.nan_to_num(fits.getdata(os.path.join(p, CORE)).astype(float).squeeze())
        mask = fits.getdata(os.path.join(p, MASK)).squeeze() > 0
        total = core.sum()
        if total > 0.0:
            frac.append(core[mask].sum() / total)
frac = np.array(frac)
print(" %d nodes: fraction of the core mass inside the footprint, median %.3f, 16th and 84th percentiles %.3f and %.3f, minimum %.3f"
      % (len(frac), np.median(frac), *np.percentile(frac, [16, 84]), frac.min()))
      