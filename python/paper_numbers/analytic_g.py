#!/usr/bin/env python3
"""analytic_g.py -- the field factor of uniform bodies, computed with the temperature relation of each radiation field.

For every field, the relation T(Sigma_A) refitted on that field's grid (laws_refit.txt) is put into the analytic calculation of
hicorrect2 (sphere and layer averaged, no calibration of the cloud, no hold at the maximum), and the factor f_G0(Sigma_F) is
computed. Its ratio to the same at G_0 = 10 is the analytic field factor, compared here with the measured one (GTAB). Then the
four constants of the G_0 = 10 relation are replaced one at a time by those of another field, to see which drives the change.

Usage: python3 analytic_g.py <laws_refit.txt>
"""
import sys
import numpy as np
import hicorrect2 as h

BANDS = [160.0, 250.0, 350.0, 500.0]
rows = np.loadtxt(sys.argv[1])
G0 = rows[:, 0]
PAR = {g: tuple(r[1:5]) for g, r in zip(G0, rows)}
XS = np.array([21.5, 21.8, 22.0, 22.2, 22.4, 22.6, 22.85])


def factor(p):
    h.TFLOOR, h.TSURF, h.SSTAR, h.GAMMA = p
    fit, f, _ = h.relation_geometry(BANDS, "mean")
    return np.interp(XS, np.log10(fit), np.log(f))


def measured(g0):
    return np.log(h.field_factor(XS, g0))


ref = factor(PAR[10.0])
print("log10 Sigma_F      " + "".join("%7.2f" % x for x in XS))
for g in G0:
    print("G0 %6.2f analytic " % g + "".join("%7.3f" % v for v in np.exp(factor(PAR[g]) - ref)))
    print("          measured " + "".join("%7.3f" % v for v in np.exp(measured(g))))

print("\nOne constant at a time, from the relation of G0 = 10 to that of the field (analytic g):")
names = ["T_min", "T_0", "Sigma_star and gamma"]
for g in (47.1, 482.5, 0.98):
    print(" G0 = %g: constants %s" % (g, " ".join("%.4g" % v for v in PAR[g])))
    for k, name in enumerate(names):
        p = list(PAR[10.0])
        if k == 0:
            p[0] = PAR[g][0]
        elif k == 1:
            p[1] = PAR[g][1]
        else:
            p[2], p[3] = PAR[g][2], PAR[g][3]
        print("   %-22s" % name + "".join("%7.3f" % v for v in np.exp(factor(tuple(p)) - ref)))
