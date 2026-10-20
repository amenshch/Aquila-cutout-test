#!/usr/bin/env python3
"""superposition_test.py -- how much the correction depends on the way the material of a line of sight is arranged.

The correction treats a line of sight as one uniform structure whose total column is the observed one.  A real line of sight
crosses a diffuse medium and one or more denser structures at unknown depths within it.  This script measures what that costs, in
two parts.

Part 1, the invariance.  For a plane-parallel geometry the emergent emission depends on the total column alone, whatever the
density profile.  The reason is that the material shielding a parcel is the column lying between it and the surface, which does
not depend on how that column is distributed, and that the mass per unit column is constant for any density profile.  Writing the
column in front of a parcel at physical depth s as

    Sigma(s) = integral of rho from 0 to s,

the emission integrated over the line of sight is

    integral B(T(Sigma_a(Sigma(s)))) rho(s) ds = integral B(T(Sigma_a(Sigma))) dSigma,

and rho has disappeared.  The first part of this script verifies that statement numerically on density profiles that differ by a
factor of thirty, as a check on the algebra and on the integration.

Part 2, the departure.  The invariance holds only because a plane-parallel layer offers no escape sideways.  A structure of finite
transverse size does, and then the arrangement matters, because the diffuse material in front of the structure and behind it is
crossed at different obliquities.  The second part measures this with a two-component line of sight that can be stated exactly:

  - a plane-parallel cirrus of column Sigma_c, in which the material at column depth x is shielded by the angular average of
    x / |mu| toward the observer and (Sigma_c - x) / |mu| away from it;
  - a uniform sphere whose column through the center is Sigma_s, embedded in that cirrus with the fraction q of the cirrus column
    in front of it and 1 - q behind it.  A parcel at the signed position z along the diameter, in units of the radius, is shielded
    by its own sphere over the path

        l(z, psi) = -z cos(psi) + sqrt ( 1 - z^2 sin^2(psi) ),

    in units of the radius, and by the cirrus over q Sigma_c / |cos(psi)| toward the observer and (1 - q) Sigma_c / |cos(psi)|
    away from it.

The angular average of Eq. (6) is taken over psi with the isotropic weight sin(psi), the temperature of each parcel follows from
the fitted relation, the emission of the cirrus and of the sphere are added with their masses, one modified blackbody is fitted to
the result exactly as the reconstruction fits it, and the correction factor is the ratio of the true column to the fitted one.
The same quantity computed for one uniform structure of the same total column is what the method applies.  The difference between
them, taken over the range of q and of Sigma_c, is the systematic uncertainty that the arrangement contributes.

The signed z matters.  With |z| in its place the sphere would be symmetric about its center whatever the cirrus in front of it,
and the answer would be wrong for every q except one half.

Requires hicorrect2.py in the same directory, from which the relation, the angular average and the constants are taken unchanged.

Usage:  python3 superposition_test.py [<output file>]
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hicorrect2 import (planck, cnu, temperature, attenuating_column, KAPHEAT, TFLOOR, TSURF, CAL_ERR)

BANDS = [160.0, 250.0, 350.0, 500.0]
NANG = 720               # directions of the angular average
NZ = 400                 # parcels along the diameter of the sphere
NX = 400                 # parcels through the cirrus
NTFIT = 2400             # temperatures of the fitted grid

PSI = (np.arange(NANG) + 0.5) / NANG * np.pi
WPSI = np.sin(PSI)
CPSI = np.cos(PSI)
MU = (np.arange(NANG) + 0.5) / NANG * 2.0 - 1.0

TGRID = np.geomspace(0.98 * TFLOOR, 1.02 * TSURF, NTFIT)
MODEL = np.array([[cnu(w) * planck(t, w) for w in BANDS] for t in TGRID])
RELERR = np.array([CAL_ERR.get(w, 0.20) for w in BANDS])


def fitted_column(obs):
    """Surface density that a single modified blackbody fitted to obs would report, from the longest waveband."""
    wgt = 1.0 / (RELERR * obs) ** 2
    amp = (wgt * obs * MODEL).sum(axis=1) / (wgt * MODEL * MODEL).sum(axis=1)
    j = int(np.argmin((wgt * (obs - amp[:, None] * MODEL) ** 2).sum(axis=1)))
    return obs[-1] / (cnu(BANDS[-1]) * planck(TGRID[j], BANDS[-1]))


def intensities(sigma, weights, temps):
    """Intensity in each waveband of a line of sight whose parcels carry the given columns and temperatures."""
    return np.array([cnu(w) * np.sum(weights * planck(temps, w)) for w in BANDS])


def one_body(sigma):
    """The correction factor the method applies: one uniform structure of the whole column."""
    t = temperature(attenuating_column(sigma))
    obs = np.array([cnu(w) * sigma * np.mean(planck(t, w)) for w in BANDS])
    return sigma / fitted_column(obs)


def layer_temperatures(sigma_c):
    """Temperature of each parcel of a plane-parallel cirrus, sampled uniformly in its own column."""
    x = (np.arange(NX) + 0.5) / NX * sigma_c
    d = np.where(MU > 0.0, x[:, None] / np.abs(MU), (sigma_c - x[:, None]) / np.abs(MU))
    j = np.exp(-KAPHEAT * d).mean(axis=1)
    return temperature(-np.log(np.maximum(j, 1e-300)) / KAPHEAT)


def sphere_temperatures(sigma_s, sigma_c, q):
    """Temperature of each parcel of a sphere embedded in the cirrus with the fraction q of that cirrus in front of it."""
    z = (np.arange(NZ) + 0.5) / NZ * 2.0 - 1.0
    lsph = -z[:, None] * CPSI + np.sqrt(np.maximum(1.0 - (z[:, None] * np.sin(PSI)) ** 2, 0.0))
    col = lsph * sigma_s / 2.0 + np.where(CPSI > 0.0, q * sigma_c / np.abs(CPSI), (1.0 - q) * sigma_c / np.abs(CPSI))
    j = (np.exp(-KAPHEAT * col) * WPSI).sum(axis=1) / WPSI.sum()
    return temperature(-np.log(np.maximum(j, 1e-300)) / KAPHEAT)


def two_component(sigma, sigma_c, q):
    """The correction factor a cirrus with an embedded sphere would require, for the same total column."""
    sigma_s = sigma - sigma_c
    if sigma_s <= 0.0:
        return np.nan
    tc = layer_temperatures(sigma_c)
    ts = sphere_temperatures(sigma_s, sigma_c, q)
    wc = np.full(NX, sigma_c / NX)
    ws = np.full(NZ, sigma_s / NZ)
    obs = intensities(sigma, np.concatenate([wc, ws]), np.concatenate([tc, ts]))
    return sigma / fitted_column(obs)


def invariance(sigma=2.0e22, nprof=3):
    """Part 1: the emergent intensities of a plane-parallel line of sight, for density profiles that differ by a factor of 30."""
    x = (np.arange(4 * NX) + 0.5) / (4 * NX) * sigma
    d = np.where(MU > 0.0, x[:, None] / np.abs(MU), (sigma - x[:, None]) / np.abs(MU))
    t = temperature(-np.log(np.maximum(np.exp(-KAPHEAT * d).mean(axis=1), 1e-300)) / KAPHEAT)
    ref = np.array([cnu(w) * sigma * np.mean(planck(t, w)) for w in BANDS])
    out = []
    s = (np.arange(40000) + 0.5) / 40000.0
    for pos in np.linspace(0.15, 0.85, nprof):
        rho = 1.0 + 30.0 * np.exp(-0.5 * ((s - pos) / 0.05) ** 2)
        col = np.cumsum(rho) / np.sum(rho) * sigma
        wgt = rho / np.sum(rho) * sigma
        tt = np.interp(col, x, t)
        val = intensities(sigma, wgt, tt)
        out.append((pos, float(np.max(np.abs(val / ref - 1.0)))))
    return ref, out


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "superposition_test.txt"
    lines = []
    lines.append("! superposition_test.py -- dependence of the correction on the arrangement of the line of sight")
    lines.append("!")
    lines.append("! Part 1: a plane-parallel line of sight of 2.0e22 cm^-2, with a clump thirty times the ambient density placed")
    lines.append("! at three depths.  The quantity given is the largest relative difference of the four intensities from those of")
    lines.append("! a line of sight of uniform density and the same total column.")
    ref, inv = invariance()
    lines.append("!   uniform: " + "  ".join("%.6e" % v for v in ref) + "  MJy sr^-1")
    for pos, dev in inv:
        lines.append("!   clump at %.2f of the depth: %.2e" % (pos, dev))
    lines.append("!")
    lines.append("! Part 2: the correction factor of a cirrus of column Sigma_c with a sphere embedded in it carrying the rest of")
    lines.append("! the column, against the fraction q of the cirrus that lies in front of the sphere, beside the factor the")
    lines.append("! method applies to the same total column.")
    lines.append("!")
    qs = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
    lines.append("! Sigma (cm^-2)  Sigma_c (cm^-2)   one body" + "".join("    q=%.2f" % q for q in qs) + "     spread")
    for sigma in (1.0e22, 2.0e22, 5.0e22, 1.0e23):
        for sigma_c in (1.0e21, 3.0e21, 1.0e22):
            if sigma_c >= 0.9 * sigma:
                continue
            one = one_body(sigma)
            vals = [two_component(sigma, sigma_c, q) for q in qs]
            lines.append("  %.2e       %.2e       %6.3f" % (sigma, sigma_c, one)
                         + "".join("   %6.3f" % v for v in vals)
                         + "    %5.1f %%" % (100.0 * (max(vals) / min(vals) - 1.0)))
    text = "\n".join(lines) + "\n"
    open(out, "w").write(text)
    print(text)
    print(" written %s" % out)


if __name__ == "__main__":
    main()


# ----------------------------------------------------------------------------------------------------------------------------
# Part 3: concentration.  The invariance of Part 1 removes the density profile from a plane-parallel line of sight, so whatever
# concentration does it must do through the transverse escape of a bounded structure.  This part measures it directly, on a sphere
# whose density follows rho(r) = rho_0 (1 + (r/a)^2)^-1 out to the radius unity, normalized so that the column through the center
# is the observed one.  Parcels are taken along the central line of sight and weighted by their density, which is what the
# emission of that line of sight weights them by, and the column each of them must send its radiation through is integrated along
# the ray to the boundary, direction by direction.

NRAY = 240               # steps of the integral along a ray to the boundary


def sphere_profile_temperatures(sigma, ratio):
    """Temperature of the parcels of the central line of sight of a sphere of central-to-edge density ratio 'ratio'.

    Returns their densities, used as the weights of the emission, and their temperatures.
    """
    a2 = 1.0 / max(ratio - 1.0, 1.0e-6)                 # rho(1)/rho(0) = 1 / (1 + 1/a^2) fixes a from the ratio
    z = (np.arange(NZ) + 0.5) / NZ * 2.0 - 1.0
    rho = 1.0 / (1.0 + z ** 2 / a2)
    norm = sigma / (np.sum(1.0 / (1.0 + z ** 2 / a2)) * 2.0 / NZ)     # column through the center equals sigma
    tau = np.empty(NZ)
    for i, zi in enumerate(z):
        # the ray from (0, 0, zi) in the direction psi leaves the sphere after the length l(zi, psi)
        l = -zi * CPSI + np.sqrt(np.maximum(1.0 - (zi * np.sin(PSI)) ** 2, 0.0))
        t = (np.arange(NRAY) + 0.5) / NRAY
        # the point reached after the fraction t of that length, and its distance from the center
        x = (l[:, None] * t[None, :]) * np.sin(PSI)[:, None]
        y = zi + (l[:, None] * t[None, :]) * CPSI[:, None]
        r2 = x ** 2 + y ** 2
        col = norm * np.sum(1.0 / (1.0 + r2 / a2), axis=1) * l / NRAY
        tau[i] = -np.log(max(np.sum(np.exp(-KAPHEAT * col) * WPSI) / np.sum(WPSI), 1e-300)) / KAPHEAT
    return rho, temperature(tau)


def concentration(sigma, ratio):
    """Correction factor required by the central line of sight of a sphere of that central-to-edge density ratio."""
    rho, t = sphere_profile_temperatures(sigma, ratio)
    w = rho / np.sum(rho) * sigma
    obs = intensities(sigma, w, t)
    return sigma / fitted_column(obs)
