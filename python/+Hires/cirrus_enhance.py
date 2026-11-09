"""cirrus_enhance.py -- the denser central cloud added to the cirrus map of the simulated skies.

The cirrus map bd.fits, scaled by CIRRUS_SCALE, has a median surface density of about 3e21 cm^-2 and no large-scale gradient.
Real clouds such as Aquila have a brighter, extended background that rises toward their centers. The map is therefore
multiplied by a Plummer-like function of the distance r from the center of the field,

    E(r) = 1 + ENH_A [1 + (r / ENH_RC)^2]^(-ENH_P / 2),

which is 1 + ENH_A at the center, flat within the radius ENH_RC, and falls as a power law r^-ENH_P beyond it, approaching 1 far
from the center. The fine structure of the map is kept. make_cirrus3d.py and make_simsky6.py both call enhance(), so that the
RADMC-3D cirrus, the envelopes of the filaments and the simulated sky always use the same cirrus.
"""

import numpy as np

ENH_A = 5.0        # amplitude: the factor at the center is 1 + ENH_A (0 switches the enhancement off)
ENH_RC = 600.0     # flat radius (arcsec)
ENH_P = 1.5        # power-law index of the decline beyond ENH_RC


def factor(shape, pix):
    """E(r) on an image of the given shape with pixels of pix arcsec, centered on the middle of the image."""
    ny, nx = shape
    y, x = np.mgrid[0:ny, 0:nx]
    r = np.hypot((y - 0.5 * (ny - 1)) * pix, (x - 0.5 * (nx - 1)) * pix)
    return 1.0 + ENH_A * (1.0 + (r / ENH_RC) ** 2) ** (-ENH_P / 2.0)


def enhance(sig, pix):
    """The cirrus map multiplied by E(r)."""
    return sig * factor(sig.shape, pix)


def keywords():
    """FITS keywords recording the enhancement."""
    return [("ENH_A", ENH_A, "central factor of the cirrus is 1 + ENH_A"),
            ("ENH_RC", ENH_RC, "flat radius of the enhancement (arcsec)"),
            ("ENH_P", ENH_P, "power-law index beyond ENH_RC")]
