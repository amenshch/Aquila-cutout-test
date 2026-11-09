#!/usr/bin/env python3
"""cloud_tables.py -- the calibration on the diffuse cloud (CIRRUS_F of hicorrect2) for several sets of wavebands and weights.

It reads the emission maps of the RADMC-3D cloud (<prefix>.rt.<band>.fits, in the units of hicorrect2: the intensity per unit of
cnu) and its true surface density (<prefix>.sigma.fits), fits one modified blackbody to every pixel with the given wavebands and
weights exactly as hicorrect2 does, and gives the median of the true over the fitted surface density in bins of 0.1 dex of the
fitted surface density, centered on log10 Sigma_F = 20.55 to 21.65 (the rows CIRRUS_LOGX of hicorrect2).

Usage: python3 cloud_tables.py <directory of the maps> [<prefix of the maps: cloud (default) or cirrus3d>]
"""
import sys
import numpy as np
from astropy.io import fits
import hicorrect2 as h

LOGX = np.round(np.arange(20.55, 21.651, 0.1), 2)


def fit_sigma(obs, waves, err):
    """Fitted surface density of every pixel: one modified blackbody, the amplitude and temperature fitted with the weights err."""
    e = np.array([err[w] for w in waves])
    tg = np.geomspace(0.98 * h.TFLOOR, 1.02 * h.TSURF, 1500)
    model = np.array([[h.cnu(w) * h.planck(t, w) for w in waves] for t in tg])
    out = np.empty(len(obs))
    for i0 in range(0, len(obs), 4000):
        o = obs[i0:i0 + 4000]
        wg = 1.0 / (e[None, :] * o) ** 2
        amp = (wg[:, None, :] * o[:, None, :] * model[None]).sum(2) / (wg[:, None, :] * model[None] ** 2).sum(2)
        chi = (wg[:, None, :] * (o[:, None, :] - amp[:, :, None] * model[None]) ** 2).sum(2)
        j = np.argmin(chi, axis=1)
        out[i0:i0 + 4000] = o[:, -1] / model[j, -1]
    return out


def table(d, waves, err, prefix):
    sig = fits.getdata("%s/%s.sigma%s.fits" % (d, prefix, ".rt" if prefix == "cirrus3d" else "")).astype(float).ravel()
    obs = np.stack([fits.getdata("%s/%s.rt.%d.fits" % (d, prefix, w)).astype(float).ravel() for w in waves], axis=1)
    sf = fit_sigma(obs, waves, err)
    lx = np.log10(sf)
    r = sig / sf
    out = []
    for x in LOGX:
        m = np.abs(lx - x) < 0.05
        out.append(np.median(r[m]) if m.sum() >= 20 else np.nan)
    return np.array(out)


def main():
    d = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else "cloud"
    cal = {70.0: 0.20, 100.0: 0.20, 160.0: 0.20, 250.0: 0.10, 350.0: 0.10, 500.0: 0.10}
    equal = {70.0: 0.10, 100.0: 0.10, 160.0: 0.10, 250.0: 0.10, 350.0: 0.10, 500.0: 0.10}
    cases = [("standard 160-500, 20/10/10/10", [160.0, 250.0, 350.0, 500.0], cal),
             ("160-500, equal weights", [160.0, 250.0, 350.0, 500.0], equal),
             ("250-500, equal weights", [250.0, 350.0, 500.0], equal),
             ("100-500, 20/20/10/10/10", [100.0, 160.0, 250.0, 350.0, 500.0], cal),
             ("70-500, 20/20/20/10/10/10", [70.0, 100.0, 160.0, 250.0, 350.0, 500.0], cal)]
    print(" log Sigma_F      " + " ".join("%6.2f" % x for x in LOGX))
    print(" CIRRUS_F now     " + " ".join("%6.3f" % v for v in h.CIRRUS_F))
    for name, w, e in cases:
        print(" %-32s\n                  " % name + " ".join("%6.3f" % v for v in table(d, w, e, prefix)), flush=True)


if __name__ == "__main__":
    main()
