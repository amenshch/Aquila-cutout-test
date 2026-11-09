"""g0_requirement.py -- the radiation fields that the dust temperatures of a
cloud require of the models of the correction.

The correction of the single-temperature bias is calibrated on radiative
transfer models computed in one radiation field, G_0 = 10 in the units of
Habing (1968). A cloud heated by a stronger field is warmer at every surface
density than the models of the same surface density. The grid of radiation
fields repeats the models in nine fields, from G_0 = 0.98 to 482.5, and
harvest_calib3.py has tabulated, for each field, the temperature that the
reconstruction fits against the surface density of the finest image.

For each bin of surface density, this script takes the median temperature of
the models in each field, which rises monotonically with G_0, and interpolates
log G_0 against it. An observed temperature is thereby converted into the
field that the models need to reach it at the same surface density. The
scaling T ~ G_0^(1/(4+beta)) of unshielded dust is not used: at a fixed
surface density, the fitted temperature of the models rises only as about
G_0^0.1, and the script reports the exponent measured in each bin.

For each cloud and each bin, the script gives the number of pixels, the 1st,
16th, 50th, 84th and 99th percentiles of the observed temperature, and the
fields that reproduce them. A temperature below the weakest field of the grid
is marked with '<', and one above the strongest with '>'. Pixels warmer than
the strongest field allows cannot be heated from outside; they are internally
heated, by protostars or nearby massive stars, and their fraction is reported.

The surface density and the temperature images must have the resolution of
the images in the harvested tables, that is, the finest image of the
reconstruction (13.5 arcsec for hires with 160 to 500 um).

Usage:

    python3 g0_requirement.py <directory of the harvest tables> <name>:<surface density image>:<temperature image> ...
                              [-bin=0.2] [-xmin=21.0] [-xmax=23.6] [-minpix=100] [-xfield=21.8]

The directory must contain the files harvest_calib3_G0=<field>.txt, with the
field written as in 0p98, 2p12, 10, 482p5. Each cloud is given as a name, its
surface density image (H2 cm^-2) and its temperature image (K), separated by
colons. Bins with fewer than minpix pixels of the cloud are not listed.

The field of a cloud is estimated from the bins below 10^xfield cm^-2 only.
There, a line of sight crosses mostly the diffuse cloud. In denser bins, the
fitted temperature is raised also by superposition (a dense structure seen
through the warm cloud is warmer than an isolated core of the same surface
density) and by internal heating, so the fields there are upper bounds. The
summary gives both estimates.
"""

import os
import sys
import glob
import numpy as np
from astropy.io import fits


def read_harvest(path):
    """Columns x (log10 surface density), y (fitted temperature, K) and the
    number of pixels of each cell of a table of harvest_calib3.py."""
    a = []
    for line in open(path):
        c = line.split()
        if not c or c[0].startswith("!"):
            continue
        a.append((float(c[0]), float(c[1]), float(c[3])))
    a = np.array(a)
    return a[:, 0], a[:, 1], a[:, 2]


def wmedian(v, w):
    """Median of the values v with the weights w."""
    k = np.argsort(v)
    v, w = v[k], w[k]
    c = np.cumsum(w)
    c = (c - 0.5 * w) / c[-1]
    return np.interp(0.5, c, v)


def option(name, default):
    for a in sys.argv[1:]:
        if a.startswith("-" + name + "="):
            return float(a.split("=", 1)[1])
    return default


def ladder(directory):
    """The fields of the grid and their harvested tables, sorted by field."""
    out = []
    for p in glob.glob(os.path.join(directory, "harvest_calib3_G0=*.txt")):
        tag = os.path.basename(p)[len("harvest_calib3_G0="):-4]
        out.append((float(tag.replace("p", ".")), read_harvest(p)))
    out.sort(key=lambda t: t[0])
    return out


def field_of(temp, gs, tmed):
    """Field that reproduces the temperature temp, from the medians tmed of
    the models in the fields gs, and a flag: -1 below, +1 above, 0 inside."""
    if temp < tmed[0]:
        return gs[0], -1
    if temp > tmed[-1]:
        return gs[-1], 1
    return 10.0 ** np.interp(temp, tmed, np.log10(gs)), 0


def rnd(g):
    """A field rounded to two significant digits, written without an exponent."""
    if not np.isfinite(g):
        return "--"
    return "%d" % int(float("%.2g" % g)) if g >= 10.0 else "%.2g" % g


def fmt(g, flag):
    s = "%.3g" % g
    return ("<" + s) if flag < 0 else ((">" + s) if flag > 0 else s)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if len(args) < 2:
        print(__doc__)
        return
    step = option("bin", 0.2)
    xmin = option("xmin", 21.0)
    xmax = option("xmax", 23.6)
    minpix = int(option("minpix", 100))
    xfield = option("xfield", 21.8)
    edges = np.arange(xmin, xmax + 0.5 * step, step)

    lad = ladder(args[0])
    if len(lad) < 2:
        print(" fewer than two harvest tables found in %s" % args[0])
        return
    print(" fields of the grid: " + ", ".join("%g" % g for g, _ in lad))

    # median model temperature in each bin and field
    med = np.full((len(edges) - 1, len(lad)), np.nan)
    for k in range(len(edges) - 1):
        for j, (g, (hx, hy, hn)) in enumerate(lad):
            m = (hx >= edges[k]) & (hx < edges[k + 1]) & (hn > 0)
            if m.any():
                med[k, j] = wmedian(hy[m], hn[m])

    print("\n median fitted temperature of the models (K), in bins of log10 Sigma, and the exponent of T against G_0")
    print(" %-12s" % "log10 Sigma" + "".join("%7s" % ("%g" % g) for g, _ in lad) + "  exponent")
    for k in range(len(edges) - 1):
        ok = np.isfinite(med[k])
        if ok.sum() == 0:
            continue
        gs = np.array([g for g, _ in lad])[ok]
        ex = np.polyfit(np.log(gs), np.log(med[k, ok]), 1)[0] if ok.sum() >= 3 else np.nan
        print(" %4.1f -- %4.1f" % (edges[k], edges[k + 1])
              + "".join("%7.1f" % v if np.isfinite(v) else "%7s" % "--" for v in med[k]) + "%9.3f" % ex)

    q = [1.0, 16.0, 50.0, 84.0, 99.0]
    summary = []
    for spec in args[1:]:
        name, fd, ft = spec.split(":")
        d = np.nan_to_num(np.array(fits.getdata(fd), float))
        t = np.nan_to_num(np.array(fits.getdata(ft), float))
        if d.shape != t.shape:
            print(" %s: the two images differ in shape, skipped" % name)
            continue
        ok = (d > 0.0) & (t > 0.0)
        x = np.log10(d[ok])
        tt = t[ok]

        print("\n %s: %d valid pixels" % (name, ok.sum()))
        print(" %-12s %9s   %-34s   %-40s %s" % ("log10 Sigma", "pixels", "T (K) at 1, 16, 50, 84, 99%",
                                                  "G_0 at 1, 16, 50, 84, 99%", "too warm"))
        nall, nhot, g16, g50, g84 = 0, 0, [], [], []
        d16, d50, d84 = [], [], []
        latex = []
        for k in range(len(edges) - 1):
            s = (x >= edges[k]) & (x < edges[k + 1])
            npx = int(s.sum())
            okm = np.isfinite(med[k])
            if npx < minpix or okm.sum() < 2:
                continue
            gs = np.array([g for g, _ in lad])[okm]
            tm = np.maximum.accumulate(med[k, okm])
            tp = np.percentile(tt[s], q)
            fl = [field_of(v, gs, tm) for v in tp]
            hot = int((tt[s] > tm[-1]).sum())
            nall += npx
            nhot += hot
            g16.append((fl[1][0], npx))
            g50.append((fl[2][0], npx))
            g84.append((fl[3][0], npx))
            if edges[k + 1] <= xfield + 1e-6:
                d16.append((fl[1][0], npx))
                d50.append((fl[2][0], npx))
                d84.append((fl[3][0], npx))
            print(" %4.1f -- %4.1f %9d   %-34s   %-40s %5.1f%%"
                  % (edges[k], edges[k + 1], npx, " ".join("%5.1f" % v for v in tp),
                     " ".join("%6s" % fmt(*f) for f in fl), 100.0 * hot / npx))
            latex.append("%4.1f -- %4.1f & %s & %.1f -- %.1f -- %.1f & %s -- %s -- %s \\\\"
                         % (edges[k], edges[k + 1], ("{:,}".format(npx)).replace(",", "\\,"),
                            tp[1], tp[2], tp[3], fmt(*fl[1]), fmt(*fl[2]), fmt(*fl[3])))
        if nall == 0:
            continue

        def wmed(lst):
            v = np.array([a for a, _ in lst])
            w = np.array([b for _, b in lst], float)
            return wmedian(v, w)

        dif = (wmed(d16), wmed(d50), wmed(d84)) if d50 else (np.nan, np.nan, np.nan)
        summary.append((name, dif, (wmed(g16), wmed(g50), wmed(g84)), nhot / nall))
        print(" %s: fields at the 16th, 50th and 84th percentiles, pixel-weighted medians over the bins below 10^%.1f:"
              " %.3g, %.3g, %.3g; over all bins: %.3g, %.3g, %.3g; %.1f%% of the pixels warmer than the strongest"
              " field allows" % ((name, xfield) + dif + summary[-1][2] + (100.0 * nhot / nall,)))
        print("\n LaTeX rows for %s (log10 Sigma, pixels, T at 16, 50, 84%%, G_0 at 16, 50, 84%%):" % name)
        for row in latex:
            print(row)

    print("\n Summary: cloud; fields at the 16th, 50th and 84th percentiles of the temperature, as pixel-weighted"
          " medians over the bins below 10^%.1f (diffuse) and over all bins; fraction of pixels too warm for"
          " external heating" % xfield)
    for name, dif, alls, h in summary:
        print(" %-16s diffuse %7.3g %7.3g %7.3g   all %7.3g %7.3g %7.3g %7.1f%%" % ((name,) + dif + alls + (100.0 * h,)))
    print("\n LaTeX rows of the summary (diffuse fields at 16, 50, 84%, fields over all bins at 50%, too warm):")
    for name, dif, alls, h in summary:
        print("%-16s & %s & %s & %s & %s & %.1f{\\%%} \\\\" % (name, rnd(dif[0]), rnd(dif[1]), rnd(dif[2]), rnd(alls[1]),
                                                     100.0 * h))


if __name__ == "__main__":
    main()
