#!/usr/bin/env python3
"""superposition_numbers.py -- Table tab:concentr and the numbers of Appendix B of the paper (concentration, two components).

It uses the functions of superposition_test.py, which must be in the same directory as hicorrect2.py:

  1. concentration(Sigma, ratio): the factor required by the central line of sight of a sphere with the density
     (1 + r^2 / a^2)^-1 and the central-to-edge density ratio 'ratio' (1 for a uniform sphere), for the four central surface
     densities of Table tab:concentr;
  2. two_component(Sigma, Sigma_cl, q): the factor required by a cloud of surface density Sigma_cl with a uniform sphere carrying
     the rest of Sigma, a fraction q of the cloud lying in front of it. The fitted surface density of that line of sight is
     Sigma / required, and the correction applies f(Sigma_F) of hicorrect2.py (G_0 = 10, held at its maximum), so the relative error
     of the corrected surface density is f(Sigma_F) / required - 1. It prints its range over q for Sigma = 2 x 10^22 and
     5 x 10^22 cm^-2 and Sigma_cl = 10^21, 3 x 10^21 and 10^22 cm^-2, and the worked example of the appendix (q = 0,
     Sigma = 2 x 10^22, Sigma_cl = 10^22 cm^-2).

Usage:

    python3 superposition_numbers.py <directory of hicorrect2.py and superposition_test.py>
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(sys.argv[1]))
import hicorrect2 as h
import superposition_test as S


def main():
    print(" concentration: factor for density ratios 1 (uniform), 10, 100, 1000")
    for sig in (5.0e21, 1.0e22, 3.0e22, 1.0e23):
        print("   %.0e  " % sig + " ".join("%.3f" % S.concentration(sig, r) for r in (1.0, 10.0, 100.0, 1000.0)), flush=True)
    fitted, factor = h.build_relation([160.0, 250.0, 350.0, 500.0])
    factor = np.maximum.accumulate(factor)
    f = lambda x: np.interp(np.log10(x), np.log10(fitted), factor)
    qs = np.linspace(0.0, 1.0, 11)
    print(" two components: range of the error over q")
    for sig in (2.0e22, 5.0e22):
        for scl in (1.0e21, 3.0e21, 1.0e22):
            req = np.array([S.two_component(sig, scl, q) for q in qs])
            err = f(sig / req) / req - 1.0
            print("   Sigma %.0e, Sigma_cl %.0e: %+.3f to %+.3f" % (sig, scl, err.min(), err.max()), flush=True)
    req = S.two_component(2.0e22, 1.0e22, 0.0)
    one = S.one_body(2.0e22)
    print(" example: composite Sigma_F %.3e, requires %.3f, applied %.3f, error %+.3f; uniform body Sigma_F %.3e"
          % (2.0e22 / req, req, f(2.0e22 / req), f(2.0e22 / req) / req - 1.0, 2.0e22 / one))


if __name__ == "__main__":
    main()
