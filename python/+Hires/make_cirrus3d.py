#!/usr/bin/env python3
"""make_cirrus3d.py -- a 3D cirrus for RADMC-3D whose projection is the cirrus map of the simulated sky.

A surface density map does not determine a 3D density: infinitely many cubes project to the same map. This script builds one of
them, with the same scale-free character along the line of sight as across the sky, in five steps.

  1. The map. The cirrus map (bd.fits times CIRRUS_SCALE, as in make_simsky6.py), with its border of zeros filled as in
     make_simsky5.py, is averaged in blocks of CELL arcsec, the size of the cells of the 3D grid. Its power spectrum gives the
     index BETA of its fluctuations (P(k) proportional to k^-BETA).

  2. The margins. The map is extended on every side by MARGIN pc, by mirror reflection. A finite cube is also illuminated through
     its sides, which a cirrus that continues beyond the field is not; the margins keep that unrealistic heating away from the
     field.

  3. The fluctuations along the line of sight. A 3D Gaussian random field G with the power spectrum k^-BETA fills a cube of depth
     DEPTH pc. Its log-normal transform g = exp(SIGMA_LN G - SIGMA_LN^2 / 2) is positive and has the mean one.

  4. The density. Every column of g is rescaled so that its integral equals the surface density of the map in that pixel:
     n(x, y, z) = Sigma(x, y) g(x, y, z) / (sum over z of g dz). The projection of the cube is therefore exactly the map.

  5. The files. amr_grid.inp (a regular Cartesian grid centered on the origin) and dust_density.binp (binary, double precision),
     with the dust-to-gas ratio DUST2GAS and the mass per H2 molecule MUH2, and cirrus3d.json with the parameters, which
     compare_cirrus3d.py reads.

Run RADMC-3D in the same directory with the cylinder module and cirrus3d = 1 in +radmc3d.inp, which makes the module read these
files instead of building a model:

    python3 make_cirrus3d.py
    radmc3d mctherm
    python3 compare_cirrus3d.py

The other RADMC-3D input files (dustopac.inp, the opacity, wavelength_micron.inp, external_source.inp, +radmc3d.inp) are those of
the cylinder models. Writing the dust temperature in binary (rto_style = 3 in radmc3d.inp) saves space; compare_cirrus3d.py reads
either form.
"""

import json
import sys

import numpy as np
from astropy.io import fits

import cirrus_enhance as E
import make_simsky5 as S5

CIRRUS = "bd.fits"
CIRRUS_SCALE = 0.1235      # as in make_simsky6.py
PIX = 3.0                  # pixel of the map (arcsec)
CELL = 12.0                # cell of the 3D grid (arcsec); a multiple of PIX
DEPTH = 2.0                # depth of the cirrus along the line of sight (pc)
MARGIN = 2.0               # margin added on every side of the map (pc)
SIGMA_LN = 0.5             # amplitude of the log-normal fluctuations along the line of sight
BETA = None                # index of the power spectrum; None measures it from the map
SEED = 20261004
DIST = 260.0               # distance (pc)
DUST2GAS = 0.01
MUH2 = 2.8
AMU = 1.66053907e-24
AU = 1.495978707e13
PC = 3.0856775814913673e18


def power_index(img):
    """The index of the azimuthally averaged power spectrum of the map, fitted between 4 and 1/4 of the image size."""
    a = img - img.mean()
    p = np.abs(np.fft.fft2(a)) ** 2
    ky = np.fft.fftfreq(a.shape[0])
    kx = np.fft.fftfreq(a.shape[1])
    k = np.hypot(*np.meshgrid(ky, kx, indexing="ij"))
    edges = np.geomspace(4.0 / min(a.shape), 0.25, 25)
    km, pm = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (k >= lo) & (k < hi)
        if m.sum() > 5:
            km.append(np.sqrt(lo * hi))
            pm.append(p[m].mean())
    return -np.polyfit(np.log(km), np.log(pm), 1)[0]


def gaussian_field(shape, beta, rng):
    """A Gaussian random field of unit variance with the power spectrum k^-beta."""
    k = np.sqrt(sum(np.meshgrid(*[np.fft.fftfreq(n) ** 2 for n in shape], indexing="ij", sparse=True)))
    amp = np.zeros(k.shape)
    amp[k > 0] = k[k > 0] ** (-beta / 2.0)
    noise = rng.standard_normal(shape)
    g = np.real(np.fft.ifftn(np.fft.fftn(noise) * amp))
    return (g - g.mean()) / g.std()


def grid_header(hdr0, fac, marg, shape):
    """The header of the map for the grid: the WCS of the cirrus map rescaled to the cells and shifted by the margins."""
    h = hdr0.copy()
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
            h["CRPIX%d" % ax] = (h["CRPIX%d" % ax] - 0.5) / fac + 0.5 + marg
    h["BUNIT"] = "cm-2"
    return h


def write_maps(hdr0, sig, fac, marg, cell_pc, nz):
    """The map on the grid, with margins (cirrus3d.sigma.fits) and without them (cirrus3d.sigma.field.fits), with the WCS and
    the keywords of the other surface density images of the simulated skies."""
    h = grid_header(hdr0, fac, marg, sig.shape)
    h["COMMENT"] = "cirrus map on the grid of the RADMC-3D cirrus, margins included"
    h["CIRSCALE"] = (CIRRUS_SCALE, "factor applied to %s" % CIRRUS)
    h["CELLSIZE"] = (CELL, "cell of the 3D grid (arcsec)")
    h["MARGIN"] = (marg, "margin added on every side (cells)")
    h["DEPTHPC"] = (DEPTH, "depth of the cirrus (pc)")
    h["SIGMALN"] = (SIGMA_LN, "log-normal amplitude along the depth")
    h["NZ"] = (nz, "cells along the line of sight")
    for kw in E.keywords():
        h[kw[0]] = (kw[1], kw[2])
    S5.write_fits("cirrus3d.sigma.fits", sig, h, [("WAVE", 0.0, "wavelength (microns); 0 for surface density")])
    hf = h.copy()
    for ax in (1, 2):
        if "CRPIX%d" % ax in hf:
            hf["CRPIX%d" % ax] = hf["CRPIX%d" % ax] - marg
    hf["COMMENT"] = "the field without its margins, as analyzed by compare_cirrus3d.py"
    S5.write_fits("cirrus3d.sigma.field.fits", sig[marg:sig.shape[0] - marg, marg:sig.shape[1] - marg], hf,
                  [("WAVE", 0.0, "wavelength (microns); 0 for surface density")])


def main():
    cell_pc = CELL * DIST / 206264.806
    fac = int(round(CELL / PIX))
    hdr0 = fits.getheader(CIRRUS)
    sig = CIRRUS_SCALE * S5.repair(np.array(fits.getdata(CIRRUS), float).squeeze(), "cirrus")     # its border of zeros filled
    beta = BETA if BETA is not None else power_index(sig)          # of the fluctuations, before the central enhancement
    if E.ENH_A > 0.0:
        sig = E.enhance(sig, PIX)                     # the denser central cloud, exactly as in make_simsky6.py
    ny, nx = (sig.shape[0] // fac) * fac, (sig.shape[1] // fac) * fac
    sig = sig[:ny, :nx].reshape(ny // fac, fac, nx // fac, fac).mean(axis=(1, 3))
    marg = int(round(MARGIN / cell_pc))
    sig = np.pad(sig, marg, mode="reflect")
    nz = max(int(round(DEPTH / cell_pc)), 2)
    nyc, nxc = sig.shape
    print(" map %d x %d pixels of %.0f arcsec averaged to %d x %d cells of %.0f arcsec (%.4f pc)" % (nx, ny, PIX, nx // fac,
          ny // fac, CELL, cell_pc))
    print(" margins of %d cells (%.2f pc) on every side: grid %d x %d x %d = %.2e cells, %.1f GB for one double array"
          % (marg, marg * cell_pc, nxc, nyc, nz, nxc * nyc * nz, nxc * nyc * nz * 8 / 1e9))
    print(" power-spectrum index %.3f%s; log-normal amplitude %.2f; depth %.2f pc" % (beta, " (measured)" if BETA is None
          else "", SIGMA_LN, DEPTH))

    rng = np.random.default_rng(SEED)
    g = np.exp(SIGMA_LN * gaussian_field((nz, nyc, nxc), beta, rng) - 0.5 * SIGMA_LN ** 2)
    dz = cell_pc * PC
    n = sig[None, :, :] * g / (g.sum(axis=0)[None, :, :] * dz)            # H2 per cm^3
    rho_dust = n * MUH2 * AMU * DUST2GAS

    side = cell_pc * PC
    with open("amr_grid.inp", "w") as f:
        f.write("1\n0\n1\n0\n1 1 1\n%d %d %d\n" % (nxc, nyc, nz))
        for nn in (nxc, nyc, nz):
            walls = (np.arange(nn + 1) - 0.5 * nn) * side
            f.write("".join("%.15e\n" % w for w in walls))
    with open("dust_density.binp", "wb") as f:
        np.array([1, 8, rho_dust.size, 1], dtype=np.int64).tofile(f)
        rho_dust.astype(np.float64).ravel().tofile(f)        # x varies fastest, then y, then z
    with open("cirrus3d.json", "w") as f:
        json.dump(dict(cirrus=CIRRUS, scale=CIRRUS_SCALE, enh_a=E.ENH_A, enh_rc=E.ENH_RC, enh_p=E.ENH_P,
                      pix=PIX, cell=CELL, cell_pc=cell_pc, depth=DEPTH, margin_cells=marg,
                       sigma_ln=SIGMA_LN, beta=beta, seed=SEED, dist=DIST, dust2gas=DUST2GAS, muh2=MUH2, nx=nxc, ny=nyc,
                       nz=nz), f, indent=1)
    write_maps(hdr0, sig, fac, marg, cell_pc, nz)
    proj = n.sum(axis=0) * dz
    print(" projection of the cube / map: %.6f to %.6f" % ((proj / sig).min(), (proj / sig).max()))
    print(" written amr_grid.inp, dust_density.binp, cirrus3d.json, cirrus3d.sigma.fits (the map on the grid) and"
          " cirrus3d.sigma.field.fits (the same without the margins)")


if __name__ == "__main__":
    sys.exit(main())
