#!/usr/bin/env python3
"""cylinder_smooth_temperature.py -- average the Monte Carlo temperatures of a cylinder along its axis.

A RADMC-3D model of a cylinder is computed on an axisymmetric grid in spherical coordinates, whose cells near the origin and near
the axis are small. Few photon packets cross them, so their temperatures scatter widely, and some are never crossed at all and are
left at zero. The middle part of a long cylinder is, however, the same at every height: its temperature depends on the distance R
from the axis alone. This script uses that symmetry to beat the noise. It bins the cells of the middle part in R, replaces the
temperature of every cell there by the mass-weighted mean temperature of its bin (cells at zero excluded), and writes the result
back to dust_temperature.dat, keeping the original as dust_temperature.mc.dat. Cells near the ends of the cylinder, where the
temperature also depends on the height, are left unchanged.

Run it in the model directory after 'radmc3d mctherm' and before the images are computed:

    python3 cylinder_smooth_temperature.py [<fraction of the half-length>] [<bins per decade in R>]

The middle part is |z| < f L/2 - R, with R the outer radius of the cylinder or of its envelope, f = 0.7 by default, and the
bins are 50 per decade of R by default. When the envelope varies along the axis (cylenvz = 1), the middle part is |z| < L/2 - 1.5 R
and every cell is averaged within a running window of DZ = 0.2 pc along the axis. The profile of the mass-weighted
temperature against R is written to cylinder_temperature_profile.txt.
"""

import os
import re
import shutil
import sys

import numpy as np

AU = 1.495978707e13
DZ = 0.2           # height of the running window along the axis (pc) when the envelope varies along the axis
PC = 3.0856775814913673e18


def read_grid():
    with open("amr_grid.inp") as fh:
        lines = [ln.split(";")[0].strip() for ln in fh]
    nr, nt = [int(v) for v in lines[5].split()[:2]]
    walls = np.array([float(v.replace("d", "e").replace("D", "e")) for v in lines[6:6 + nr + 1 + nt + 1]])
    return nr, nt, walls[:nr + 1], walls[nr + 1:]


def read_scalar(name):
    with open(name) as fh:
        head = [fh.readline() for _ in range(3)]
        data = np.array(fh.read().split(), float)
    return head, data


def envz_flag():
    """True if the envelope varies along the axis (cylenvz = 1)."""
    for name in ("radmc3d.inp", "+radmc3d.inp"):
        if os.path.exists(name):
            m = re.search(r"^\s*cylenvz\s*=\s*([0-9]+)", open(name).read(), re.M)
            if m:
                return int(m.group(1)) == 1
    return False


def model_parameters():
    """The outer radius and the length of the cylinder, from radmc3d.inp or +radmc3d.inp, in pc."""
    for name in ("radmc3d.inp", "+radmc3d.inp"):
        if os.path.exists(name):
            text = open(name).read()
            rout = re.search(r"^\s*cylrout\s*=\s*([0-9.eEdD+-]+)", text, re.M)
            clen = re.search(r"^\s*cyllen\s*=\s*([0-9.eEdD+-]+)", text, re.M)
            embs = re.search(r"^\s*cylembsd\s*=\s*([0-9.eEdD+-]+)", text, re.M)
            renv = re.search(r"^\s*cylrenv\s*=\s*([0-9.eEdD+-]+)", text, re.M)
            if rout and clen:
                r = float(rout.group(1).replace("d", "e"))
                # with an envelope, the ends of the cylinder disturb the temperatures over the radius of the envelope
                if embs and renv and float(embs.group(1).replace("d", "e")) > 0.0:
                    r = max(r, float(renv.group(1).replace("d", "e")))
                return r, float(clen.group(1).replace("d", "e"))
    print(" cylrout and cyllen not found in radmc3d.inp")
    sys.exit(1)


def main():
    frac = float(sys.argv[1]) if len(sys.argv) > 1 else 0.7
    perdec = float(sys.argv[2]) if len(sys.argv) > 2 else 50.0
    rout, clen = model_parameters()
    nr, nt, rw, tw = read_grid()
    _, dens = read_scalar("dust_density.inp")
    source = "dust_temperature.mc.dat" if os.path.exists("dust_temperature.mc.dat") else "dust_temperature.dat"
    thead, temp = read_scalar(source)
    print(" Monte Carlo temperatures read from %s" % source)
    if dens.size != nr * nt or temp.size != nr * nt:
        print(" the files hold %d and %d cells, the grid %d; only one dust species is supported" % (dens.size, temp.size, nr * nt))
        sys.exit(1)

    ir = np.tile(np.arange(nr), nt)
    it = np.repeat(np.arange(nt), nr)
    rc = 0.5 * (rw[ir] + rw[ir + 1])
    tc = 0.5 * (tw[it] + tw[it + 1])
    R = rc * np.sin(tc)
    z = rc * np.cos(tc)
    vol = (2.0 * np.pi / 3.0) * (rw[ir + 1] ** 3 - rw[ir] ** 3) * (np.cos(tw[it]) - np.cos(tw[it + 1]))
    mass = dens * vol

    zmax = frac * 0.5 * clen * PC - rout * PC
    if envz_flag():
        # the whole cylinder except its ends, which are disturbed over about the radius of the envelope
        zmax = 0.5 * clen * PC - 1.5 * rout * PC
    if zmax <= 0.0:
        print(" the cylinder is too short for a middle part: shorten nothing, or lower the fraction")
        sys.exit(1)
    mid = (np.abs(z) < zmax) & (dens > 0.0)
    lo = np.log10(max(R[mid].min(), 1.0))
    hi = np.log10(R[mid].max())
    edges = 10.0 ** np.arange(lo, hi + 1.0 / perdec, 1.0 / perdec)
    k = np.clip(np.searchsorted(edges, R) - 1, 0, len(edges) - 2)

    # With an envelope that varies along the axis, the temperature depends on z as well. Every cell of the middle part then gets the
    # mass-weighted mean of the cells of its bin of R within DZ / 2 of its own height: a running window, which follows the
    # variation along the axis smoothly, without steps between fixed windows.
    envz = envz_flag()
    new = temp.copy()
    rows = []
    for b in range(len(edges) - 1):
        cells = np.nonzero(mid & (k == b))[0]
        if cells.size == 0:
            continue
        good = cells[temp[cells] > 0.0]
        if good.size == 0:
            continue
        if envz:
            o = np.argsort(z[good])
            zg = z[good][o]
            cm = np.concatenate([[0.0], np.cumsum(mass[good][o])])
            ct = np.concatenate([[0.0], np.cumsum(mass[good][o] * temp[good][o])])
            lo = np.searchsorted(zg, z[cells] - 0.5 * DZ * PC, side="left")
            hi = np.searchsorted(zg, z[cells] + 0.5 * DZ * PC, side="right")
            ok = hi > lo
            new[cells[ok]] = (ct[hi[ok]] - ct[lo[ok]]) / (cm[hi[ok]] - cm[lo[ok]])
            mid_cells = good[np.abs(z[good]) < 0.5 * DZ * PC]
        else:
            new[cells] = np.average(temp[good], weights=mass[good])
            mid_cells = good
        if mid_cells.size == 0:
            continue
        tmean = np.average(temp[mid_cells], weights=mass[mid_cells])
        rows.append((np.sqrt(edges[b] * edges[b + 1]) / AU, mid_cells.size, 1.0 - good.size / cells.size, tmean,
                     np.sqrt(np.average((temp[mid_cells] - tmean) ** 2, weights=mass[mid_cells]))))
    if envz:
        print(" the envelope varies along the axis: running windows of %.2f pc along the axis" % DZ)

    if not os.path.exists("dust_temperature.mc.dat"):
        shutil.copyfile("dust_temperature.dat", "dust_temperature.mc.dat")
    with open("dust_temperature.dat", "w") as fh:
        fh.writelines(thead)
        for v in new:
            fh.write("%16.8e\n" % v)
    with open("cylinder_temperature_profile.txt", "w") as fh:
        fh.write("! mass-weighted dust temperature of the middle part of the cylinder, |z| < %.3f pc, against R\n" % (zmax / PC))
        fh.write("!%13s %8s %10s %10s %10s\n" % ("R(AU)", "cells", "frac T=0", "T(K)", "rms(K)"))
        for r, n, f, t, s in rows:
            fh.write(" %13.5e %8d %10.3f %10.4f %10.4f\n" % (r, n, f, t, s))
    print(" %d cells of the middle part, |z| < %.3f pc, averaged in %d bins of R; original kept as dust_temperature.mc.dat"
          % (mid.sum(), zmax / PC, len(rows)))
    print(" profile written to cylinder_temperature_profile.txt")


if __name__ == "__main__":
    main()
