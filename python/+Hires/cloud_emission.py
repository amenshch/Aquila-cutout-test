#!/usr/bin/env python3
"""cloud_emission.py -- emission maps and true surface density of a RADMC-3D cloud, in the units of hicorrect2.

It reads dust_density.binp, dust_temperature.bdat and cirrus3d.json of a cloud computed by make_cirrus3d.py, integrates the
emission of every column of the field (margins excluded) along the line of sight in the wavebands 70 to 500 um, in units of the
modified blackbody per molecule of hicorrect2 (cnu times the Planck function), and writes cloud.rt.<band>.fits and cloud.sigma.fits,
as compare_cirrus3d.py does.

Usage: python3 cloud_emission.py <directory of the cloud> <output directory>
"""
import json
import os
import sys
import numpy as np
from astropy.io import fits
import hicorrect2 as H

AMU = 1.66053907e-24
PC = 3.0856775814913673e18
BANDS = [70.0, 100.0, 160.0, 250.0, 350.0, 500.0]


def open_binary(name, shape):
    with open(name, "rb") as f:
        hdr = np.fromfile(f, dtype=np.int64, count=4)
    return np.memmap(name, dtype=np.float64 if hdr[1] == 8 else np.float32, mode="r", offset=32, shape=shape)


def main():
    d, out = sys.argv[1], sys.argv[2]
    par = json.load(open(os.path.join(d, "cirrus3d.json")))
    nx, ny, nz, m = par["nx"], par["ny"], par["nz"], par["margin_cells"]
    dz = par["cell_pc"] * PC
    rho = open_binary(os.path.join(d, "dust_density.binp"), (nz, ny, nx))
    tem = open_binary(os.path.join(d, "dust_temperature.bdat"), (nz, ny, nx))
    conv = 1.0 / (par["dust2gas"] * par["muh2"] * AMU)
    fy, fx = ny - 2 * m, nx - 2 * m
    sig = np.zeros((fy, fx))
    em = np.zeros((len(BANDS), fy, fx))
    for y0 in range(0, fy, 16):
        y1 = min(y0 + 16, fy)
        col = np.array(rho[:, m + y0:m + y1, m:nx - m]) * conv * dz
        t = np.maximum(np.array(tem[:, m + y0:m + y1, m:nx - m]), 2.7)
        sig[y0:y1] = col.sum(axis=0)
        for k, w in enumerate(BANDS):
            em[k, y0:y1] = (col * H.cnu(w) * H.planck(t, w)).sum(axis=0)
    os.makedirs(out, exist_ok=True)
    fits.writeto(os.path.join(out, "cloud.sigma.fits"), sig.astype(np.float32), overwrite=True)
    for k, w in enumerate(BANDS):
        fits.writeto(os.path.join(out, "cloud.rt.%d.fits" % w), em[k].astype(np.float32), overwrite=True)
    print(" field %d x %d, surface density %.3e to %.3e, median %.3e" % (fx, fy, sig.min(), sig.max(), np.median(sig)))


if __name__ == "__main__":
    main()
