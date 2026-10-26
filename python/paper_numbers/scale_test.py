#!/usr/bin/env python3
"""scale_test.py -- the analytic field factor when T_min and T_0 of the G_0 = 10 relation are scaled by the same factor.

The ratio T_0 / T_min, the temperature contrast along a line of sight, is then unchanged, and only the absolute temperatures rise,
which isolates the effect of the Planck function at the wavelengths of the fit.

Usage: python3 scale_test.py <laws_refit.txt>
"""
import sys
import numpy as np
import hicorrect2 as h

BANDS = [160.0, 250.0, 350.0, 500.0]
rows = np.loadtxt(sys.argv[1])
p10 = tuple(rows[rows[:, 0] == 10.0][0, 1:5])
XS = np.array([21.5, 21.8, 22.0, 22.2, 22.4, 22.6, 22.85])


def factor(p):
    h.TFLOOR, h.TSURF, h.SSTAR, h.GAMMA = p
    fit, f, _ = h.relation_geometry(BANDS, "mean")
    return np.interp(XS, np.log10(fit), np.log(f))


ref = factor(p10)
print(" T0 / T_min of the fields: " + " ".join("%g:%.2f" % (r[0], r[2] / r[1]) for r in rows))
print(" log10 Sigma_F          " + "".join("%7.2f" % x for x in XS))
for k in (0.8, 1.15, 1.3, 1.5, 1.75):
    p = (p10[0] * k, p10[1] * k, p10[2], p10[3])
    print(" both scaled by %.2f    " % k + "".join("%7.3f" % v for v in np.exp(factor(p) - ref)))
