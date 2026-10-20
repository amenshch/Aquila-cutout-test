#!/usr/bin/env python3
"""compare_cirrus3d.py -- the RADMC-3D cirrus against the analytic cirrus of the simulated sky, and the correction tested on it.

Run it in the directory of the 3D cirrus after 'radmc3d mctherm'. It needs make_cirrus3d.py's cirrus3d.json, amr_grid.inp and
dust_density.binp, the dust temperatures (dust_temperature.bdat or .dat), and hicorrect2.py and make_simsky5.py for the relation
T(Sigma_a), the slab shielding of the analytic cirrus and the correction.

For every column of the field (the margins excluded) it computes:

  1. the emission at 160, 250, 350 and 500 microns of the RADMC-3D column, from its densities and temperatures;
  2. the emission of the analytic cirrus with the same surface density: a uniform slab, every layer shielded by the angular average
     of the column toward both faces, with the temperature from the relation T(Sigma_a), as in make_simsky5.py and make_simsky6.py;
  3. for both, the surface density that one modified blackbody fitted to the four intensities returns, with the weights of
     hicorrect2.py, and that surface density corrected with the factor f of hicorrect2.py.

Both emissions use the same modified blackbody per molecule (cnu times the Planck function of hicorrect2.py), so their ratio
measures only the difference of the temperatures. The results are printed in bins of surface density and written as images:
cirrus3d.ratio250.fits (RADMC-3D over analytic emission at 250 microns), cirrus3d.corr_rt.fits and cirrus3d.corr_an.fits (corrected
over true surface density for the RADMC-3D and for the analytic cirrus) and cirrus3d.tmean.fits (mass-weighted RADMC-3D
temperature of each column).
"""

import json
import os
import sys

import numpy as np
from astropy.io import fits

import hicorrect2 as H
import make_simsky5 as S5

S5.TFLOOR, S5.TSURF, S5.DSTAR, S5.GAMMA = H.TFLOOR, H.TSURF, H.SSTAR, H.GAMMA      # the relation of hicorrect2.py
BANDS = [160.0, 250.0, 350.0, 500.0]
ALLBANDS = [70.0, 100.0, 160.0, 250.0, 350.0, 500.0]       # the bands of the simulated sky, for the emission maps
AMU = 1.66053907e-24
PC = 3.0856775814913673e18


def open_binary(name, shape):
    """The binary cube as a read-only memory map (no copy in memory), shape (nz, ny, nx)."""
    with open(name, "rb") as f:
        hdr = np.fromfile(f, dtype=np.int64, count=4)
    return np.memmap(name, dtype=np.float64 if hdr[1] == 8 else np.float32, mode="r", offset=32, shape=shape)


def open_temperature(shape):
    if os.path.exists("dust_temperature.bdat"):
        return open_binary("dust_temperature.bdat", shape)
    with open("dust_temperature.dat") as f:
        f.readline(), f.readline(), f.readline()
        return np.array(f.read().split(), float)[:np.prod(shape)].reshape(shape)


def fit_mbb(obs):
    """Surface density of one modified blackbody fitted to the intensities obs (pixels x bands), with the weights of hicorrect2."""
    err = np.array([H.CAL_ERR.get(w, 0.2) for w in BANDS])
    tg = np.geomspace(0.98 * H.TFLOOR, 1.02 * H.TSURF, 1500)
    model = np.array([[H.cnu(w) * H.planck(t, w) for w in BANDS] for t in tg])          # (ntemp, nband) per H2 cm^-2
    out = np.empty(len(obs))
    for i0 in range(0, len(obs), 4000):
        o = obs[i0:i0 + 4000]
        wg = 1.0 / (err[None, :] * o) ** 2
        amp = (wg[:, None, :] * o[:, None, :] * model[None]).sum(2) / (wg[:, None, :] * model[None] ** 2).sum(2)
        chi = (wg[:, None, :] * (o[:, None, :] - amp[:, :, None] * model[None]) ** 2).sum(2)
        j = np.argmin(chi, axis=1)
        out[i0:i0 + 4000] = o[:, -1] / model[j, -1]
    return out


def analytic_emission(sig):
    """Emission of the analytic slab cirrus with the surface densities sig, at the four bands, per pixel."""
    zlayer = (np.arange(S5.NLAYER) + 0.5) / S5.NLAYER * 2.0 - 1.0
    t0c = S5.KAPHEAT * sig / 2.0
    out = np.zeros((len(sig), len(BANDS)))
    for zz in zlayer:
        T = S5.T_of_shielding(S5.slab_seff(np.full_like(sig, zz), t0c) * sig / 2.0)
        for k, w in enumerate(BANDS):
            out[:, k] += sig / S5.NLAYER * H.cnu(w) * H.planck(T, w)
    return out


def field_header(par):
    """The header of the images of the field (the margins excluded), with the usual keywords of the simulated skies.

    The coordinates are those of cirrus3d.sigma.field.fits if make_cirrus3d.py wrote it. Otherwise they are derived from the
    cirrus map itself (par["cirrus"], normally bd.fits, which must then be in the directory): its CDELT (or CD) multiplied by the
    number of map pixels per cell, and CRPIX rescaled to the cells. Without either file the images get no coordinates."""
    if os.path.exists("cirrus3d.sigma.field.fits"):
        h = fits.getheader("cirrus3d.sigma.field.fits")
    elif os.path.exists(par["cirrus"]):
        fac = int(round(par["cell"] / par["pix"]))
        h = fits.getheader(par["cirrus"]).copy()
        for key in ("NAXIS1", "NAXIS2", "NAXIS3", "BUNIT", "BEAM", "WAVE", "BZERO", "BSCALE", "DATAMIN", "DATAMAX"):
            if key in h:
                del h[key]
        for ax in (1, 2):
            if "CDELT%d" % ax in h:
                h["CDELT%d" % ax] = h["CDELT%d" % ax] * fac
            for i in (1, 2):
                if "CD%d_%d" % (ax, i) in h:
                    h["CD%d_%d" % (ax, i)] = h["CD%d_%d" % (ax, i)] * fac
            if "CRPIX%d" % ax in h:
                h["CRPIX%d" % ax] = (h["CRPIX%d" % ax] - 0.5) / fac + 0.5
    else:
        print(" neither cirrus3d.sigma.field.fits nor %s found: the images get no coordinates" % par["cirrus"])
        h = fits.Header()
    h["CIRSCALE"] = (par["scale"], "factor applied to %s" % par["cirrus"])
    h["CELLSIZE"] = (par["cell"], "cell of the 3D grid (arcsec)")
    h["DEPTHPC"] = (par["depth"], "depth of the cirrus (pc)")
    h["SIGMALN"] = (par["sigma_ln"], "log-normal amplitude along the depth")
    h["NZ"] = (par["nz"], "cells along the line of sight")
    return h


def write_image(name, arr, hdr, unit, comment, wave=None):
    """An image of the field with the usual keywords: BUNIT, BZERO = 0, BSCALE = 1 and WAVE (0 for surface density)."""
    h = hdr.copy()
    h["BUNIT"] = unit
    h["COMMENT"] = comment
    if "WAVE" in h:
        del h["WAVE"]
    S5.write_fits(name, arr, h, [("WAVE", 0.0 if wave is None else wave,
                                  "wavelength (microns); 0 for surface density" if wave is None else "wavelength (microns)")])


def main():
    par = json.load(open("cirrus3d.json"))
    nx, ny, nz, m = par["nx"], par["ny"], par["nz"], par["margin_cells"]
    dz = par["cell_pc"] * PC
    rho_all = open_binary("dust_density.binp", (nz, ny, nx))
    temp_all = open_temperature((nz, ny, nx))
    conv = 1.0 / (par["dust2gas"] * par["muh2"] * AMU)                         # dust density to H2 per cm^3
    # the field is processed in strips of rows, so that the cubes are never held in memory whole
    fy, fx = ny - 2 * m, nx - 2 * m
    sig = np.zeros((fy, fx))
    tsum = np.zeros((fy, fx))
    rt = np.zeros((fy, fx, len(BANDS)))
    rtall = np.zeros((len(ALLBANDS), fy, fx))
    for y0 in range(0, fy, 16):
        y1 = min(y0 + 16, fy)
        col = np.array(rho_all[:, m + y0:m + y1, m:nx - m]) * conv * dz          # H2 per cm^2 in each cell
        temp = np.maximum(np.array(temp_all[:, m + y0:m + y1, m:nx - m]), 2.7)
        sig[y0:y1] = col.sum(axis=0)
        tsum[y0:y1] = (col * temp).sum(axis=0)
        for k, w in enumerate(ALLBANDS):
            rtall[k, y0:y1] = (col * H.cnu(w) * H.planck(temp, w)).sum(axis=0)
        for k, w in enumerate(BANDS):
            rt[y0:y1, :, k] = rtall[ALLBANDS.index(w), y0:y1]
    shape = sig.shape
    tmean = tsum / sig
    rt = rt.reshape(-1, len(BANDS))
    an = analytic_emission(sig.ravel())
    fitted, factor = H.build_relation(BANDS)
    res = {}
    for name, em in (("rt", rt), ("an", an)):
        sf = fit_mbb(em)
        res[name] = (sf * np.interp(np.log10(sf), np.log10(fitted), factor) / sig.ravel()).reshape(shape)
    ratio = (rt[:, 1] / an[:, 1]).reshape(shape)

    print(" field of %d x %d columns (margins excluded); surface density %.2e to %.2e cm^-2" % (shape[1], shape[0],
          sig.min(), sig.max()))
    print(" mass-weighted temperature of the field: RADMC-3D %.2f K" % ((tmean * sig).sum() / sig.sum()))
    print("\n  Sigma bin (cm^-2)      columns   I250 RT/analytic   corrected/true RT   corrected/true analytic   T_RT (K)")
    edges = np.geomspace(sig.min(), sig.max() * 1.0001, 9)
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (sig >= lo) & (sig < hi)
        if k.sum() < 10:
            continue
        print("  %.2e-%.2e   %7d      %.3f              %.3f               %.3f                 %.2f"
              % (lo, hi, k.sum(), np.median(ratio[k]), np.median(res["rt"][k]), np.median(res["an"][k]), np.median(tmean[k])))
    hdr = field_header(par)
    for nm, arr, unit, com in (("cirrus3d.ratio250.fits", ratio, "", "RADMC-3D / analytic emission at 250 microns"),
                               ("cirrus3d.corr_rt.fits", res["rt"], "", "corrected / true surface density, RADMC-3D"),
                               ("cirrus3d.corr_an.fits", res["an"], "", "corrected / true surface density, analytic"),
                               ("cirrus3d.tmean.fits", tmean, "K", "mass-weighted RADMC-3D dust temperature")):
        write_image(nm, arr, hdr, unit, com)
    # the emission of the RADMC-3D cirrus in the six bands, column by column, in the units of hicorrect2.cnu times the Planck
    # function (per H2 per cm^2); make_simsky6.py divides it by the analytic emission computed in the same units
    for k, w in enumerate(ALLBANDS):
        write_image("cirrus3d.rt.%03d.fits" % int(w), rtall[k], hdr, "cnu*B per cm-2",
                    "emission of the RADMC-3D cirrus at %g microns, in the units of hicorrect2" % w, wave=w)
    write_image("cirrus3d.sigma.rt.fits", sig, hdr, "cm-2", "surface density of the field of the RADMC-3D cirrus")
    print("\n written cirrus3d.ratio250.fits, cirrus3d.corr_rt.fits, cirrus3d.corr_an.fits, cirrus3d.tmean.fits,")
    print(" the emission maps cirrus3d.rt.<band>.fits and the surface density of the field cirrus3d.sigma.rt.fits")


if __name__ == "__main__":
    sys.exit(main())
