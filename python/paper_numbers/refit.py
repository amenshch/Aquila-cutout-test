#!/usr/bin/env python3
"""refit.py -- refit the relation T(Sigma_A) of each radiation field, with wide limits, from the tables of fit_attenuation_law.py.

The cost is that of fit_attenuation_law.py: the mass-weighted sum of squares of ln(T_fit / T_model), with the weights of the table.
It also prints the weighted median temperature in bins of Sigma_A and the departure of the fit from it.

Usage: python3 refit.py <directory of attenuation_law_G0=*.txt> <output file of constants>
"""
import sys
import numpy as np
from scipy.optimize import least_squares

TAGS = ["0p98", "2p12", "4p61", "10", "21p7", "47p1", "102p4", "222p2", "482p5"]
G0 = [0.98, 2.12, 4.61, 10.0, 21.7, 47.1, 102.4, 222.2, 482.5]


def law(p, s):
    tf, t0, ss, g = p
    return tf + (t0 - tf) * (1.0 + s / ss) ** (-g)


def unpack(x):
    sig = lambda v: 1.0 / (1.0 + np.exp(-np.clip(v, -40, 40)))
    return 2 + 18 * sig(x[0]), 10 + 50 * sig(x[1]), 10 ** (19 + 4 * sig(x[2])), 0.05 + 5.95 * sig(x[3])


def wmedian(v, w):
    o = np.argsort(v)
    c = np.cumsum(w[o])
    return v[o][np.searchsorted(c, 0.5 * c[-1])]


def main():
    d, out = sys.argv[1], sys.argv[2]
    rows = []
    for tag, g0 in zip(TAGS, G0):
        a = np.loadtxt("%s/attenuation_law_G0=%s.txt" % (d, tag), comments="!")
        s, t, w = a[:, 0], a[:, 1], a[:, 2]
        w = w / w.sum()
        ok = (s > 0) & (t > 0) & (w > 0)
        s, t, w = s[ok], t[ok], w[ok] / w[ok].sum()
        sw = np.sqrt(w)
        res = lambda p: sw * np.log(law(p, s) / t)
        lo, hi = [2.0, 10.0, 1e19, 0.05], [20.0, 60.0, 1e23, 6.0]
        best = None
        for p0 in ([7, 20, 2e21, 1.2], [8, 30, 4e21, 1.8], [10, 40, 8e21, 2.5], [6, 25, 1e21, 1.0]):
            r = least_squares(res, p0, bounds=(lo, hi), x_scale=[1, 1, 1e21, 0.3], xtol=1e-12, ftol=1e-14, max_nfev=20000)
            if best is None or r.cost < best.cost:
                best = r
        best.fun = 2 * best.cost
        p = tuple(best.x)
        sc = 100 * (np.exp(np.sqrt(best.fun)) - 1)
        rows.append((g0,) + tuple(p) + (sc,))
        print("G0 = %6.2f   Tmin %6.3f  T0 %6.3f  Sstar %.3e  gamma %5.3f  scatter %5.2f%%" % rows[-1])
        e = np.arange(19.5, 23.6, 0.25)
        ls = np.log10(s)
        line = "   median T / fit, bins of log Sigma_A:"
        for i in range(len(e) - 1):
            m = (ls >= e[i]) & (ls < e[i + 1])
            if w[m].sum() > 1e-4:
                line += " %.2f:%.3f" % (0.5 * (e[i] + e[i + 1]), wmedian(t[m], w[m]) / law(p, 10 ** (0.5 * (e[i] + e[i + 1]))))
        print(line)
    np.savetxt(out, np.array(rows), fmt="%.6g", header="G0 Tmin T0 Sstar gamma scatter_percent")


if __name__ == "__main__":
    main()
