#!/usr/bin/env python3
"""calib2d.py -- the correction factor as a function of the two quantities that a single-temperature fit returns.

Why two quantities.  The calibration of Sect. 4.4 follows one line of sight through the center of a uniform structure, so for each
trial surface density it returns one fitted surface density and one fitted temperature, and the pairs lie on a curve.  Every real
pixel gives such a pair as well, because the reconstruction writes both a surface density and a temperature, and a pixel that
falls off that curve is stating, in observables alone, that its line of sight is not the center of a uniform structure.  The
distance from the curve is a second coordinate, and the correction can be a function of both.

What this script computes.  It builds a population of lines of sight, each of a kind that occurs in a real cloud, and measures
three numbers for each: the surface density a single modified blackbody fit returns, the temperature it returns, and the factor
that restores the true column.  The population covers

  - a uniform sphere and a uniform cylinder, of every central column, sampled at every impact parameter, which represents a core
    or a filament seen anywhere from its crest to its edge;
  - each of those embedded in a plane-parallel diffuse medium of every column, at every depth within it, which represents the
    same structure superposed on a cirrus, warm foreground and all;
  - the diffuse medium by itself, which is the plane-parallel limit and reproduces the calibration of Sect. 4.4.

The three numbers of every member are binned in the plane of the two fitted quantities, and the correction adopted for a cell is
the median of the members falling in it.  That is what the surface represents: given what the fit returns, the expected correction
averaged over every configuration that could have produced it.  Configurations the two observables cannot tell apart are averaged
over rather than chosen between, and the spread within a cell is the irreducible uncertainty of the method at that point.

Geometry.  A parcel at the position x inside the structure, with the line of sight along z, is shielded along a direction n by the
path to the surface of its own structure,

  sphere:    t = -x.n + sqrt ( 1 - |x|^2 + (x.n)^2 ),
  cylinder:  t = [ -x_p.n_p + sqrt ( |n_p|^2 ( 1 - |x_p|^2 ) + (x_p.n_p)^2 ) ] / |n_p|^2,

where the subscript p denotes the projection onto the plane perpendicular to the axis, both in units of the radius, and by the
column of the diffuse medium lying along n, which is q Sigma_c / |n_z| toward the observer and (1 - q) Sigma_c / |n_z| away from
it.  The average of Eq. (5) is taken over directions spread evenly over the sphere, the temperature of each parcel follows from
Eq. (7), the emission of all the parcels is added with their masses, and one modified blackbody is fitted to the sum in the
wavebands and with the weights the reconstruction uses.

Writes calib2d.npz, holding the members, the binned surface and the spread within each cell.

Usage:  python3 calib2d.py [<output file>]
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hicorrect2 import planck, cnu, temperature, KAPHEAT, TFLOOR, TSURF, CAL_ERR

BANDS = [160.0, 250.0, 350.0, 500.0]
RELERR = np.array([CAL_ERR.get(w, 0.20) for w in BANDS])

NDIR = 600               # directions of the angular average, spread evenly over the sphere
NRAY = 90                # parcels along the part of the line of sight inside the structure
NCIR = 60                # parcels through the diffuse medium
NTFIT = 2000             # temperatures of the fitted grid

SIG0 = np.geomspace(1.0e21, 6.0e23, 26)                 # column through the center of the structure
IMPACT = np.array([0.0, 0.2, 0.35, 0.5, 0.62, 0.72, 0.8, 0.87, 0.93, 0.97])
SIGC = np.array([0.0, 1.0e21, 3.0e21, 6.0e21, 1.0e22, 2.0e22])       # column of the diffuse medium
QFRONT = np.array([0.0, 0.25, 0.5])                     # its fraction lying in front of the structure
SIGC_ALONE = np.geomspace(3.0e20, 6.0e22, 40)           # the diffuse medium by itself

TGRID = np.geomspace(0.98 * TFLOOR, 1.02 * TSURF, NTFIT)
MODEL = np.array([[cnu(w) * planck(t, w) for w in BANDS] for t in TGRID])

_k = np.arange(NDIR) + 0.5
_nz = 1.0 - 2.0 * _k / NDIR
_nr = np.sqrt(np.maximum(1.0 - _nz ** 2, 0.0))
_ph = np.pi * (1.0 + 5.0 ** 0.5) * _k
NX, NY, NZ = _nr * np.cos(_ph), _nr * np.sin(_ph), _nz


def fitted(obs):
    """Surface density and temperature that one modified blackbody fitted to obs returns."""
    wgt = 1.0 / (RELERR * obs) ** 2
    amp = (wgt * obs * MODEL).sum(axis=1) / (wgt * MODEL * MODEL).sum(axis=1)
    j = int(np.argmin((wgt * (obs - amp[:, None] * MODEL) ** 2).sum(axis=1)))
    return obs[-1] / (cnu(BANDS[-1]) * planck(TGRID[j], BANDS[-1])), TGRID[j]


def cirrus_temperatures(sigma_c):
    """Temperature of the parcels of a plane-parallel medium of column sigma_c, sampled uniformly in its own column."""
    x = (np.arange(NCIR) + 0.5) / NCIR * sigma_c
    mu = (np.arange(NDIR) + 0.5) / NDIR * 2.0 - 1.0
    d = np.where(mu > 0.0, x[:, None] / np.abs(mu), (sigma_c - x[:, None]) / np.abs(mu))
    j = np.exp(-KAPHEAT * d).mean(axis=1)
    return temperature(-np.log(np.maximum(j, 1e-300)) / KAPHEAT)


def structure_temperatures(sigma0, p, shape, sigma_c, q):
    """Temperature of the parcels of a ray at impact parameter p through a structure embedded in a diffuse medium.

    Returns the column of that ray through the structure and the temperatures of its parcels, which carry equal mass because the
    structure is uniform.
    """
    half = np.sqrt(max(1.0 - p * p, 1.0e-8))
    l = (np.arange(NRAY) + 0.5) / NRAY * 2.0 * half - half
    r2 = (p * p + l * l)[:, None]
    xn = p * NX[None, :] + l[:, None] * NZ[None, :]
    if shape == "sphere":
        t = -xn + np.sqrt(np.maximum(1.0 - r2 + xn ** 2, 0.0))
    else:
        npar = (NX ** 2 + NZ ** 2)[None, :]
        t = (-xn + np.sqrt(np.maximum(npar * (1.0 - r2) + xn ** 2, 0.0))) / np.maximum(npar, 1.0e-8)
    col = t * sigma0 / 2.0
    if sigma_c > 0.0:
        mu = np.maximum(np.abs(NZ[None, :]), 1.0e-3)
        col = col + np.where(NZ[None, :] > 0.0, q * sigma_c / mu, (1.0 - q) * sigma_c / mu)
    j = np.exp(-KAPHEAT * col).mean(axis=1)
    return sigma0 * half, temperature(-np.log(np.maximum(j, 1e-300)) / KAPHEAT)


def member(sigma0, p, shape, sigma_c, q):
    """One line of sight: its fitted surface density, its fitted temperature and the factor that restores its true column."""
    sig_s, ts = structure_temperatures(sigma0, p, shape, sigma_c, q)
    if sigma_c > 0.0:
        tc = cirrus_temperatures(sigma_c)
        wgt = np.concatenate([np.full(len(ts), sig_s / len(ts)), np.full(len(tc), sigma_c / len(tc))])
        tt = np.concatenate([ts, tc])
    else:
        wgt, tt = np.full(len(ts), sig_s / len(ts)), ts
    sig = sig_s + sigma_c
    obs = np.array([cnu(w) * np.sum(wgt * planck(tt, w)) for w in BANDS])
    sf, tf = fitted(obs)
    return np.log10(sf), tf, sig / sf


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "calib2d.npz"
    rows = []
    for sc in SIGC_ALONE:
        tc = cirrus_temperatures(sc)
        obs = np.array([cnu(w) * sc * np.mean(planck(tc, w)) for w in BANDS])
        sf, tf = fitted(obs)
        rows.append((np.log10(sf), tf, sc / sf))
    for shape in ("sphere", "cylinder"):
        for s0 in SIG0:
            for p in IMPACT:
                for sc in SIGC:
                    for q in (QFRONT if sc > 0.0 else [0.0]):
                        rows.append(member(s0, p, shape, sc, q))
        print(" %s done" % shape)
    rows = np.array(rows)

    ex = np.arange(19.4, 23.8, 0.05)
    ey = np.arange(7.0, 23.0, 0.20)
    pts, val, spread = [], [], []
    for i in range(len(ex) - 1):
        m1 = (rows[:, 0] >= ex[i]) & (rows[:, 0] < ex[i + 1])
        if not m1.any():
            continue
        sub = rows[m1]
        for j in range(len(ey) - 1):
            m2 = (sub[:, 1] >= ey[j]) & (sub[:, 1] < ey[j + 1])
            if m2.sum() >= 3:
                f = sub[m2, 2]
                pts.append((0.5 * (ex[i] + ex[i + 1]), 0.5 * (ey[j] + ey[j + 1])))
                val.append(np.median(f))
                spread.append(0.5 * (np.percentile(f, 84) - np.percentile(f, 16)) / np.median(f))
    pts, val, spread = np.array(pts), np.array(val), np.array(spread)
    np.savez(out, members=rows, points=pts, factor=val, spread=spread)
    print(" %d members, %d cells" % (len(rows), len(val)))
    print(" log Sigma_fit %.2f to %.2f, T_fit %.2f to %.2f K, factor %.3f to %.3f"
          % (rows[:, 0].min(), rows[:, 0].max(), rows[:, 1].min(), rows[:, 1].max(), val.min(), val.max()))
    print(" spread within a cell: median %.1f%%, 90th percentile %.1f%%"
          % (100.0 * np.median(spread), 100.0 * np.percentile(spread, 90)))
    print(" written %s" % out)


if __name__ == "__main__":
    main()
