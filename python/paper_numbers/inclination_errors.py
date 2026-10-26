#!/usr/bin/env python3
"""inclination_errors.py -- the over-correction of inclined layers and filaments (Sect. 7, subsection on inclination).

Seen at an angle i from its normal (or from the normal to its axis), a layer or a filament keeps its temperature profile along the
line of sight, while its surface density and its fitted surface density both grow by the factor 1/cos i. The factor that it
requires is that of the face-on line of sight, whereas the correction applies the factor of the stretched one. The relative error of
the corrected surface density is therefore f(s Sigma_F) / f(Sigma_F) - 1, with f the factor of hicorrect2.py (G_0 = 10, held at its
maximum) and s the stretch. It prints the error at i = 30, 45, 60 and 70 degrees, and the median, mean and 84th percentile of the
error for random orientations: s = 1/mu for a layer and s = 1/sqrt(1 - mu^2) for a filament, with mu uniform on (0, 1).

Usage:

    python3 inclination_errors.py <directory of hicorrect2.py>
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(sys.argv[1]))
import hicorrect2 as h


def main():
    fitted, factor = h.build_relation([160.0, 250.0, 350.0, 500.0])
    factor = np.maximum.accumulate(factor)
    fa = lambda s: np.interp(np.log10(s), np.log10(fitted), factor)
    xs = [21.0, 21.5, 22.0, 22.5, 23.0]
    print(" log10 Sigma_F   error at i = 30, 45, 60, 70 degrees")
    for x in xs:
        print("   %.1f          " % x + "  ".join("%+4.0f %%" % (100.0 * (fa(10.0 ** x / np.cos(np.radians(d))) / fa(10.0 ** x) - 1.0))
                                         for d in (30, 45, 60, 70)))
    mu = np.random.default_rng(1).uniform(0.0, 1.0, 200000)
    for name, s in (("layer", 1.0 / np.maximum(mu, 1.0e-3)), ("filament", 1.0 / np.sqrt(np.maximum(1.0 - mu ** 2, 1.0e-6)))):
        print(" %s, random orientations: median, mean, 84th percentile of the error" % name)
        for x in xs:
            e = fa(10.0 ** x * s) / fa(10.0 ** x) - 1.0
            print("   %.1f   %+5.1f %%  %+5.1f %%  %+5.1f %%" % (x, 100.0 * np.median(e), 100.0 * e.mean(), 100.0 * np.percentile(e, 84)))


if __name__ == "__main__":
    main()
