#!/usr/bin/env python3
"""
fit_attenuation_law.py -- fit the relation between the dust temperature and the
effective attenuating column, on a grid of RADMC-3D models.

This script produces the relation T(Sigma_A) of Eq. (5) of the paper.

What it does, for each model node:

  1. reads the radius, the number density of molecular hydrogen and the dust
     temperature of every shell from mc.denstemps.tab;
  2. for every shell, computes the column Sigma(n) of material that the
     incoming radiation crosses on its way in along each of NANG directions n,
     by integrating the density along that ray through the spherical model;
  3. replaces the whole angular distribution by the single column that would
     produce the same mean intensity, the effective attenuating column of
     Eq. (4),

         exp( -sigma Sigma_A ) = < exp( -sigma Sigma(n) ) >  over 4 pi ,

     with the extinction cross section sigma gray at the wavelengths that do
     the heating;
  4. weights every shell by its mass, n r^2 dr.

The shells of all the nodes are then fitted together with

    T(Sigma_A) = T_min + ( T_0 - T_min ) ( 1 + Sigma_A / Sigma_star )^(-gamma)

by minimizing the mass-weighted sum of squares of ln(T_fit / T_model).

The mass of a shell is computed here from the density and the shell volume.
The MASSFRAC column of mc.denstemps.tab is never read, so its value does not
matter to this fit.

Usage:

    python3 fit_attenuation_law.py <root directory> [<node> ...] [-file=<name>] [-out=<file>]
                                   [-exclude=<text>] [-noclean] [-keepouter] [-cylinders=<dir>,<dir>,...] [-cylweight=<fraction>]

With no node names the whole tree below the root is searched, at any depth,
for files called mc.denstemps.tab (or the name given with -file), and every
directory holding one is taken as a node, named by its path relative to the
root. Node names given explicitly are also paths relative to the root. With
-default, the fifteen nodes of an early test fit are used instead of the
search. The fitted constants are printed and written, with the table of
shells behind them (Sigma_A, T, mass weight and the index of the node, listed in
the header), to attenuation_law.txt in the current directory, or to the
file given with -out=. The noisy temperatures of the innermost shells of dense
cores are rejected before the fit (see clean_shells), unless -noclean is given.
The outermost shell of every sphere is also rejected, unless -keepouter is
given. In the thicker embedding clouds that shell spans an extinction of several
magnitudes, and its single Monte Carlo temperature, an average over a steep
gradient, is warmer by up to 2.4 K than the dust at its center, where Sigma_A is
evaluated. The mass weights of every sphere are renormalized after the
rejections, so that every sphere keeps the same total weight.
With -exclude=<text>, repeatable, every node whose path contains the text is
left out, for example -exclude=additional for the three models that are
optically thick in the far infrared.

With -cylinders, the radial temperature profiles of RADMC-3D cylinder models (the directories holding their +radmc3d.inp and
cylinder_temperature_profile.txt) are fitted together with the spheres. The cylinders together carry the fraction of the
total weight given by -cylweight (0.5 by default), shared equally among them, the spheres the rest.

Columns of mc.denstemps.tab, one-based: 1 index, 2 radius in cm, 3 relative
zone size, 4 radius in au, 5 radius in pc, 6 radius in arcsec, 7 mass density,
8 number density of H2 in cm^-3, 9 dust mass density, 10 dust number density,
11 mass fraction, 12 dust temperature in K.
"""

import os
import sys
import glob
import numpy as np
from scipy.optimize import minimize

# Extinction cross section per hydrogen molecule at the wavelengths that heat the dust (sigma in the paper), taken gray with
# A_V = N(H2) / 0.94e21 and tau = A_V / 1.086. This is an average over the ultraviolet and optical, where the extinction in
# fact varies, and it is the only assumption of the calculation that is not a geometric one.
KAPHEAT = 1.0 / (0.94e21 * 1.086)

NANG = 240           # directions used in the angular average
NRAY = 120           # samples along each ray
CLEAN = True         # reject the noisy central shells, see clean_shells; -noclean keeps them
NOOUTER = True       # reject the outermost shell of every sphere; -keepouter keeps it
RMIN_CELL = 0        # innermost cells to discard; 0, because clean_shells
                     # rejects the noisy central shells instead

# The fifteen nodes of an early test fit, used with -default. They span a
# factor of 256 in the surface density of the embedding cloud, a factor of 12
# in the size of the core and a factor of 270 in its central to edge density
# ratio. The first eight vary the concentration, the next four the embedding
# cloud, the last three the size.
DEFAULT_NODES = ["i01j19k13", "i01j20k13", "i01j22k13", "i01j24k13",
                 "i01j25k13", "i01j26k13", "i01j27k13", "i01j28k13",
                 "i03j24k13", "i05j24k13", "i07j24k13", "i09j24k13",
                 "i01j24k09", "i01j24k16", "i01j24k19"]


TABLE = "mc.denstemps.tab"   # the file of each node, changed with -file=
OUT = "attenuation_law.txt"  # the output, in the current directory, changed with -out=


def read_node(path):
    """Radius in cm, number density of H2 in cm^-3, dust temperature in K."""
    r, n, t = [], [], []
    for line in open(os.path.join(path, TABLE)):
        if line.lstrip().startswith("!") or not line.strip():
            continue
        c = line.split()
        r.append(float(c[1]))
        n.append(float(c[7]))
        t.append(float(c[11]))
    return np.array(r), np.array(n), np.array(t)


def shell_edges(r):
    e = np.empty(len(r) + 1)
    e[1:-1] = 0.5 * (r[:-1] + r[1:])
    e[0] = max(2.0 * r[0] - e[1], 0.0)
    e[-1] = 2.0 * r[-1] - e[-2]
    return e


def node_columns(path):
    """Effective attenuating column, temperature and mass weight of each shell.

    The model is spherically symmetric, so a ray leaving a shell at radius r in
    a direction making an angle psi with the outward radial reaches the
    boundary after a distance L(psi) = -r cos psi + sqrt(R^2 - r^2 sin^2 psi).
    The density is integrated along that ray to give the attenuating column
    Sigma(n) in that direction, and the directions are averaged with the
    solid-angle weights sin psi, as in Eq. (4).
    """
    r, n, t = read_node(path)
    e = shell_edges(r)
    R = e[-1]
    mass = n * r ** 2 * np.diff(e)
    mass = mass / mass.sum()

    ps = np.linspace(0.0, np.pi, NANG + 1)
    ps = 0.5 * (ps[1:] + ps[:-1])
    wps = np.sin(ps)
    u = np.linspace(0.0, 1.0, NRAY)

    seff = np.empty(len(r))
    for i, ri in enumerate(r):
        lray = (-ri * np.cos(ps)
                + np.sqrt(np.maximum(R ** 2 - (ri * np.sin(ps)) ** 2, 0.0)))
        # points along each ray, by the law of cosines
        pts = np.sqrt(ri ** 2 + (lray[:, None] * u) ** 2
                      + 2.0 * ri * (lray[:, None] * u) * np.cos(ps)[:, None])
        dens = np.interp(np.clip(pts, r[0], R), r, n,
                         left=n[0], right=n[-1])
        col = _trapz(dens, lray[:, None] * u, axis=1) \
            if hasattr(np, "trapezoid") else _trapz(dens,
                                                      lray[:, None] * u,
                                                      axis=1)
        j = float(np.sum(wps * np.exp(-KAPHEAT * col)) / np.sum(wps))
        seff[i] = -np.log(max(j, 1e-300)) / KAPHEAT
    k = RMIN_CELL
    seff, t, mass = seff[k:], t[k:], mass[k:]
    if CLEAN:
        keep = clean_shells(t)
        seff, t, mass = seff[keep], t[keep], mass[keep]
    if NOOUTER:
        seff, t, mass = seff[:-1], t[:-1], mass[:-1]
        mass = mass / mass.sum()
    return seff, t, mass


def clean_shells(t, inner=10, tol=0.10):
    """The shells whose temperature is not Monte Carlo noise.

    Near the center of a dense core few photon packets arrive, and the temperatures of the innermost cells scatter wildly
    about the smooth profile, down to zero. A shell is rejected if its temperature is not positive, or if, among the inner
    shells, it departs by more than tol from the running median of five neighbors; every shell inside the outermost rejected
    one is rejected with it. The rejected shells carry almost no mass, so the fit hardly changes, but they no longer appear
    in the table and the figure made from it.
    """
    n = len(t)
    med = np.array([np.median(t[max(0, i - 2):i + 3]) for i in range(n)])
    bad = (t <= 0.0) | ((np.arange(n) < inner) & (np.abs(t / np.where(med > 0.0, med, 1.0) - 1.0) > tol))
    keep = np.ones(n, bool)
    if bad.any():
        keep[:np.nonzero(bad)[0].max() + 1] = False
    return keep


# ------------------------------------------------------------------------------------------------ cylinders
PC = 3.0856775814913673e18
AUCM = 1.495978707e13
_trapz = getattr(np, "trapezoid", None) or np.trapz


def read_cyl_parameters(path):
    """The parameters of a cylinder model from its +radmc3d.inp (or radmc3d.inp): surface density on the axis (cm^-2),
    half-maximum width H of the surface density (arcsec), index p, outer radius (pc) and distance (pc)."""
    import re
    for name in ("+radmc3d.inp", "radmc3d.inp"):
        f = os.path.join(path, name)
        if os.path.exists(f):
            text = open(f).read()
            val = {}
            for key in ("cylsd", "cylfwhm", "cylpind", "cylrout", "distpc"):
                m = re.search(r"^\s*%s\s*=\s*([0-9.eEdD+-]+)" % key, text, re.M)
                if m:
                    val[key] = float(m.group(1).replace("d", "e").replace("D", "e"))
            if len(val) == 5:
                return val
    raise ValueError("cylinder parameters not found in %s" % path)


def cyl_chord(b, rflat, rout, pind):
    """Surface density through the cylinder at the impact parameter b, in units of the central density times cm."""
    if b >= rout:
        return 0.0
    smax = np.sqrt(rout ** 2 - b ** 2)
    sa = min(smax, rflat)
    s = np.concatenate([np.linspace(0.0, sa, 400), np.geomspace(sa * 1.0001, max(smax, sa * 1.0002), 3000)])
    s = s[s <= smax]
    return 2.0 * _trapz((1.0 + (b ** 2 + s ** 2) / rflat ** 2) ** (-pind / 2.0), s)


def cylinder_columns(path, nrad=45, nphi=180, nmu=200):
    """Effective attenuating column, dust temperature and mass weight at radii of a cylinder model.

    The temperatures are those of cylinder_temperature_profile.txt, written by cylinder_smooth_temperature.py: the Monte Carlo
    temperatures of the middle of the cylinder averaged along its axis. The density profile is rebuilt from the parameters of
    the model exactly as the user module builds it, the flat radius being solved from the half-maximum width H of the
    surface density of the truncated profile. The effective attenuating column is the angular average of Eq. (4), the
    cylinder being taken as infinite. A direction is given by its azimuth phi in the cross-section and its polar angle theta
    to the axis, sampled uniformly in phi and in mu = cos(theta), so that every direction has the same solid angle. The column
    Sigma(n) along a direction is the column along the chord of the cross-section divided by sin(theta). Every radius is
    weighted by the mass of the cylinder at that radius, n(R) R dR per unit length.
    """
    par = read_cyl_parameters(path)
    rout = par["cylrout"] * PC
    pind = par["cylpind"]
    bhalf = 0.5 * par["cylfwhm"] * par["distpc"] * AUCM
    lo, hi = np.log(1.0e-3 * bhalf), np.log(1.0e1 * bhalf)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        rf = np.exp(mid)
        if cyl_chord(bhalf, rf, rout, pind) / cyl_chord(0.0, rf, rout, pind) > 0.5:
            hi = mid
        else:
            lo = mid
    rflat = np.exp(0.5 * (lo + hi))
    nc = par["cylsd"] / cyl_chord(0.0, rflat, rout, pind)

    prof = np.loadtxt(os.path.join(path, "cylinder_temperature_profile.txt"), comments="!")
    R, T = prof[:, 0] * AUCM, prof[:, 3]
    keep = (R > 0.05 * rflat) & (R < 0.95 * rout) & (T > 0.0)
    R, T = R[keep], T[keep]
    idx = np.unique(np.linspace(0, len(R) - 1, nrad).astype(int))
    R, T = R[idx], T[idx]

    phi = (np.arange(nphi) + 0.5) / nphi * 2.0 * np.pi
    mu = (np.arange(nmu) + 0.5) / nmu * 2.0 - 1.0
    sint = np.sqrt(1.0 - mu ** 2)
    ux, uy = np.cos(phi), np.sin(phi)
    seff = np.empty(len(R))
    for i, r0 in enumerate(R):
        b = r0 * ux
        smax = -b + np.sqrt(b ** 2 - (r0 ** 2 - rout ** 2))
        cols = np.empty(nphi)
        for k in range(nphi):
            s = np.concatenate([np.linspace(0.0, min(smax[k], 4.0 * rflat), 150),
                                np.geomspace(min(smax[k], 4.0 * rflat) * 1.0001, max(smax[k], 4.0 * rflat * 1.0002), 150)])
            s = s[s <= smax[k]]
            rr = np.sqrt((r0 + s * ux[k]) ** 2 + (s * uy[k]) ** 2)
            cols[k] = _trapz(nc * (1.0 + (rr / rflat) ** 2) ** (-pind / 2.0), s)
        j = np.mean(np.exp(-KAPHEAT * cols[:, None] / np.maximum(sint[None, :], 1.0e-3)))
        seff[i] = -np.log(max(j, 1e-300)) / KAPHEAT
    dlnr = np.gradient(np.log(R))
    wgt = nc * (1.0 + (R / rflat) ** 2) ** (-pind / 2.0) * R * R * dlnr
    return seff, T, wgt


# The limits of the four parameters. They are wide enough for every radiation field of the grids, from G_0 = 0.98 to 483, where the
# unattenuated dust reaches 40 K and more; the earlier limits (T_0 up to 30 K, gamma up to 3) stopped the fits of the strongest fields.
TFMIN, TFMAX = 2.0, 20.0          # T_min, K
T0MIN, T0MAX = 10.0, 60.0         # T_0, K
LSMIN, LSMAX = 19.0, 23.0         # log10 Sigma_star, cm^-2
GAMIN, GAMAX = 0.05, 6.0          # gamma


def unpack(x):
    """Bounded parameters, so that a small sample cannot run to a degenerate
    solution with an enormous Sigma_star compensated by an enormous gamma."""
    sig = lambda v: 1.0 / (1.0 + np.exp(-np.clip(v, -40.0, 40.0)))
    tf = TFMIN + (TFMAX - TFMIN) * sig(x[0])                           # T_min
    t0 = T0MIN + (T0MAX - T0MIN) * sig(x[1])                           # T_0
    sstar = 10.0 ** (LSMIN + (LSMAX - LSMIN) * sig(x[2]))              # Sigma_star
    gam = GAMIN + (GAMAX - GAMIN) * sig(x[3])                          # gamma
    return tf, t0, sstar, gam


def law(p, s):
    tf, t0, sstar, gam = p
    return tf + (t0 - tf) * (1.0 + s / sstar) ** (-gam)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    global TABLE, OUT, CLEAN, NOOUTER
    exclude = []
    cyldirs = []
    cylweight = 0.5
    args = []
    default = False
    for a in sys.argv[1:]:
        if a.startswith("-file="):
            TABLE = a[6:]
        elif a.startswith("-out="):
            OUT = a[5:]
        elif a == "-default":
            default = True
        elif a == "-noclean":
            CLEAN = False
        elif a == "-keepouter":
            NOOUTER = False
        elif a.startswith("-exclude="):
            exclude.append(a[9:])
        elif a.startswith("-cylinders="):
            cyldirs += [x for x in a[11:].split(",") if x]
        elif a.startswith("-cylweight="):
            cylweight = float(a[11:])
        else:
            args.append(a)
    root = args[0]
    names = args[1:] if len(args) > 1 else None
    if default:
        names = DEFAULT_NODES
    if names is None:
        names = []
        for d, dirs, files in os.walk(root, followlinks=True):
            dirs[:] = sorted(x for x in dirs if not x.startswith("."))
            if TABLE in files:
                names.append(os.path.relpath(d, root))
        names = sorted(names)
        print(" %d directories holding %s found below %s" % (len(names), TABLE, root))
        if not names:
            print(" nothing to fit: check the name of the file with -file=, or give the nodes explicitly")
            return
    if exclude:
        dropped = [x for x in names if any(e in x for e in exclude)]
        names = [x for x in names if x not in dropped]
        print(" %d nodes excluded: %s" % (len(dropped), ", ".join(dropped)))
    used, A, T, W, I = [], [], [], [], []
    for name in names:
        path = os.path.join(root, name)
        if not os.path.exists(os.path.join(path, TABLE)):
            print(" %-14s no %s, skipped" % (name, TABLE))
            continue
        a, t, w = node_columns(path)
        I.append(np.full(len(a), len(used)))
        used.append(name)
        A.append(a)
        T.append(t)
        W.append(w)
        print(" %-14s %3d shells, Sigma_a from %9.3e to %9.3e cm^-2,"
              " T from %6.2f to %6.2f K"
              % (name, len(a), a.min(), a.max(), t.min(), t.max()))
    if not used:
        print(" no usable node found in %s" % root)
        return
    A = np.concatenate(A)
    T = np.concatenate(T)
    W = np.concatenate(W)
    I = np.concatenate(I)
    good = (A > 0.0) & (T > 0.0) & (W > 0.0)
    A, T, W, I = A[good], T[good], W[good], I[good]
    W = W / W.sum()
    nsph = len(A)
    if cyldirs:
        CA, CT, CW, CI = [], [], [], []
        for cdir in cyldirs:
            a, t, w = cylinder_columns(cdir)
            CI.append(np.full(len(a), len(used)))
            CA.append(a)
            CT.append(t)
            CW.append(w / w.sum() / len(cyldirs))
            used.append("cylinder:" + os.path.basename(os.path.abspath(cdir)))
            print(" %-14s %3d radii, Sigma_a from %9.3e to %9.3e cm^-2, T from %6.2f to %6.2f K"
                  % (os.path.basename(os.path.abspath(cdir)), len(a), a.min(), a.max(), t.min(), t.max()))
        A = np.concatenate([A] + CA)
        T = np.concatenate([T] + CT)
        W = np.concatenate([W * (1.0 - cylweight)] + [w * cylweight for w in CW])
        I = np.concatenate([I] + CI)
        print(" cylinders carry %.2f of the total weight, spheres %.2f" % (cylweight, 1.0 - cylweight))

    def cost(x):
        return float(np.sum(W * np.log(law(unpack(x), A) / T) ** 2))

    best = (np.inf, None)
    for seed in ([-0.5, 0.5, 0.4, -0.6], [0.0, 0.0, 0.0, 0.0],
                 [0.5, 1.0, 0.8, 0.5], [-1.0, -0.5, -0.5, -1.0]):
        res = minimize(cost, seed, method="Nelder-Mead",
                       options={"maxiter": 60000, "maxfev": 60000,
                                "fatol": 1e-14, "xatol": 1e-8})
        if res.fun < best[0]:
            best = (res.fun, res.x)
    p = unpack(best[1])
    tf, t0, sstar, gam = p
    for name, v, lo, hi in (("T_min", tf, TFMIN, TFMAX), ("T_0", t0, T0MIN, T0MAX),
                            ("log10 Sigma_star", np.log10(sstar), LSMIN, LSMAX), ("gamma", gam, GAMIN, GAMAX)):
        if min(v - lo, hi - v) < 1e-3 * (hi - lo):
            print(" WARNING: %s = %g lies at a limit of the fit (%g to %g); widen the limits" % (name, v, lo, hi))
    pred = law(p, A)
    rms = np.sqrt(np.sum(W * np.log(pred / T) ** 2) / W.sum())

    print()
    print(" T(Sigma_a) = %.4f + (%.4f - %.4f) (1 + Sigma_a / %.5e)^(-%.5f) K"
          % (tf, t0, tf, sstar, gam))
    print(" %d nodes, %d shells, mass-weighted scatter %.1f per cent"
          % (len(used), len(A), 100.0 * (np.exp(rms) - 1.0)))

    out = os.path.abspath(OUT)
    with open(out, "w") as fh:
        fh.write("! Dust temperature against the effective attenuating"
                 " column, fitted on\n! %d RADMC-3D nodes, %d shells,"
                 " mass weighted.\n" % (len(used), len(A)))
        if cyldirs:
            fh.write("! %d sphere shells and %d cylinder radii; the cylinders carry %.2f of the total weight\n"
                     % (nsph, len(A) - nsph, cylweight))
        fh.write("! outermost shell of every sphere %s\n" % ("rejected" if NOOUTER else "kept"))
        fh.write("! kappa at the heating wavelengths: %.6e cm^2 per H2\n"
                 % KAPHEAT)
        fh.write("! T(Sigma_a) = Tf + (T0 - Tf) (1 + Sigma_a / Sstar)^(-gamma)\n")
        fh.write("!   Tf    = %14.6f K\n" % tf)
        fh.write("!   T0    = %14.6f K\n" % t0)
        fh.write("!   Sstar = %14.6e cm^-2\n" % sstar)
        fh.write("!   gamma = %14.6f\n" % gam)
        fh.write("!   mass-weighted scatter %.2f per cent\n"
                 % (100.0 * (np.exp(rms) - 1.0)))
        fh.write("! nodes, with the index used in the fourth column:\n")
        for k, name in enumerate(used):
            fh.write("!   node %6d  %s\n" % (k, name))
        fh.write("!\n!     Sigma_a (cm^-2)      T_model (K)      mass weight    node\n")
        for a, t, w, k in zip(A, T, W, I):
            fh.write("   %18.6e %16.4f %16.6e %7d\n" % (a, t, w, k))
    print(" written %s" % out)


if __name__ == "__main__":
    main()
