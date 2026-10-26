#!/usr/bin/env python3
"""law_checks.py -- the checks of the temperature relation quoted in Sect. 3.2 of the paper.

1. With the table written by fit_attenuation_law.py (attenuation_law.txt, four columns: Sigma_A, T, mass weight, node index), the
   mass-weighted mean relative departure of the dust temperature from the fitted relation, T / T(Sigma_A) - 1, for all the spheres
   together and for each cylinder separately ("the spheres are 0.6 percent colder than the relation and the cylinders 0.1 to 1.1
   percent warmer").
2. With the file of the fits that leave out one cylinder each (attenuation_law-leave_one_cylinder_out.txt, the headers of the seven
   fits, each preceded by a line with its name), the largest change of the correction factor of hicorrect2.py between 10^20 and
   10^23.5 cm^-2 relative to the adopted relation ("leaving out any one cylinder changes the factor by at most 1.5 percent").

Usage:

    python3 law_checks.py <directory of hicorrect2.py> <attenuation_law.txt> [<leave-one-out file>]
"""

import os
import re
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(sys.argv[1]))
import hicorrect2 as h


def constants(text):
    """The four fitted constants from a header of attenuation_law.txt."""
    get = lambda key: float(re.search(r"!\s+%s\s*=\s*(\S+)" % key, text).group(1))
    return get("Tf"), get("T0"), get("Sstar"), get("gamma")


def departures(path):
    text = open(path).read()
    tf, t0, ss, ga = constants(text)
    nsph = int(re.search(r"!\s+(\d+) sphere shells", text).group(1))
    d = np.loadtxt(path, comments="!")
    s, t, w, k = d[:, 0], d[:, 1], d[:, 2], d[:, 3].astype(int)
    ok = t > 0.0
    r = t / (tf + (t0 - tf) * (1.0 + s / ss) ** (-ga)) - 1.0
    sph = np.zeros(len(s), bool)
    sph[:nsph] = True
    print(" spheres: mean departure %+.2f %%" % (100.0 * np.average(r[sph & ok], weights=w[sph & ok])))
    for kk in np.unique(k[nsph:]):
        m = (k == kk) & ok
        print(" cylinder node %d: mean departure %+.2f %%" % (kk, 100.0 * np.average(r[m], weights=w[m])))


def leave_one_out(path):
    text = open(path).read()
    x = np.arange(20.0, 23.51, 0.05)
    waves = [160.0, 250.0, 350.0, 500.0]

    def factor():
        f, a = h.build_relation(waves)
        return np.interp(x, np.log10(f), np.maximum.accumulate(a))

    ref = factor()
    for block in re.split(r"\n(?=law_)", text):
        name = block.split(":")[0].strip()
        h.TFLOOR, h.TSURF, h.SSTAR, h.GAMMA = constants(block)
        d = factor() / ref - 1.0
        i = int(np.argmax(np.abs(d)))
        print(" %-22s largest change of the factor %.2f %% at log10 Sigma_F = %.2f" % (name, 100.0 * abs(d[i]), x[i]))


if __name__ == "__main__":
    departures(sys.argv[2])
    if len(sys.argv) > 3:
        leave_one_out(sys.argv[3])
