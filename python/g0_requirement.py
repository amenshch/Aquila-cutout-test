"""g0_requirement.py -- the range of radiation fields that a cloud requires of
the models of the correction.

The correction of the single-temperature bias is calibrated on radiative
transfer models computed in one radiation field, G_0 = 10 in the units of
Habing (1968). A cloud heated by a stronger field is warmer at every surface
density than any model of the same surface density, and its pixels then lie
outside the calibration. This script carries out the procedure of the
appendix "Estimating the field that a cloud requires" of the paper on the
correction:

 1. The pixels of the surface density image of a cloud are binned by their
    surface density.
 2. In each bin, the range of the fitted temperature of the cloud is taken,
    from its 1st to its 99th percentile.
 3. In the same bin, the same range is taken for the models, from the table
    written by harvest_calib3.py, weighting each of its cells by its pixels.
 4. The field that reproduces an observed temperature T_obs in that bin is
    G_0 = G_ref (T_obs / T_mod)^(4 + beta), evaluated for the lowest
    temperatures of the cloud and of the models, and for the highest ones.
 5. The models must cover the fields from the smallest to the largest of
    these values over all bins.

The surface density and the temperature images must have the resolution of
the images in the harvested table, that is, the finest image of the
reconstruction (13.5 arcsec for hires with 160 to 500 um).

Usage:

    python3 g0_requirement.py <harvest table> <name>:<surface density image>:<temperature image> ...
                              [-gref=10] [-beta=2] [-bin=0.2] [-xmin=21.4] [-xmax=23.4]

Each cloud is given as a name, its surface density image (H2 cm^-2) and its
temperature image (K), separated by colons. For each cloud the script prints
a table of the bins, and at the end one summary line per cloud: the fraction
of its pixels within the range of the models, and the smallest and largest
fields that the procedure requires. A LaTeX version of both is also printed.
"""

import sys
import numpy as np
from astropy.io import fits


def read_harvest(path):
    """Columns x (log10 surface density), y (fitted temperature, K) and the
    number of pixels of each cell of the table of harvest_calib3.py."""
    x, y, n = [], [], []
    for line in open(path):
        c = line.split()
        if not c or c[0].startswith("!"):
            continue
        x.append(float(c[0]))
        y.append(float(c[1]))
        n.append(float(c[3]))
    return np.array(x), np.array(y), np.array(n)


def wpercentile(v, w, q):
    """Percentiles q (in percent) of the values v with the weights w."""
    k = np.argsort(v)
    v, w = v[k], w[k]
    c = np.cumsum(w)
    c = (c - 0.5 * w) / c[-1]
    return np.interp(np.asarray(q) / 100.0, c, v)


def option(name, default):
    for a in sys.argv[1:]:
        if a.startswith("-" + name + "="):
            return float(a.split("=", 1)[1])
    return default


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if len(args) < 2:
        print(__doc__)
        return
    gref = option("gref", 10.0)
    beta = option("beta", 2.0)
    step = option("bin", 0.2)
    xmin = option("xmin", 21.4)
    xmax = option("xmax", 23.4)
    edges = np.arange(xmin, xmax + 0.5 * step, step)

    hx, hy, hn = read_harvest(args[0])
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
        print(" %-13s %9s %15s %15s %8s %17s" % ("log10 Sigma", "pixels", "models (K)", "cloud (K)",
                                                 "inside", "G_0 required"))
        rows = []
        ninside, ntotal, glo, ghi = 0, 0, np.inf, 0.0
        for k in range(len(edges) - 1):
            a, b = edges[k], edges[k + 1]
            s = (x >= a) & (x < b)
            npx = int(s.sum())
            if npx == 0:
                continue
            olo, ohi = np.percentile(tt[s], [1.0, 99.0])
            m = (hx >= a) & (hx < b)
            ntotal += npx
            if m.sum() == 0 or hn[m].sum() == 0:
                rows.append((a, b, npx, None, None, olo, ohi, None, None, None))
                print(" %5.1f -- %4.1f %9d %15s %6.1f -- %4.1f %8s %17s"
                      % (a, b, npx, "no models", olo, ohi, "--", "--"))
                continue
            mlo, mhi = wpercentile(hy[m], hn[m], [1.0, 99.0])
            inside = ((tt[s] >= mlo) & (tt[s] <= mhi)).sum()
            ninside += inside
            g1 = gref * (olo / mlo) ** (4.0 + beta)
            g2 = gref * (ohi / mhi) ** (4.0 + beta)
            glo, ghi = min(glo, g1, g2), max(ghi, g1, g2)
            rows.append((a, b, npx, mlo, mhi, olo, ohi, inside / npx, g1, g2))
            print(" %5.1f -- %4.1f %9d %6.1f -- %4.1f %6.1f -- %4.1f %7.1f%% %7.2f -- %6.1f"
                  % (a, b, npx, mlo, mhi, olo, ohi, 100.0 * inside / npx, g1, g2))
        frac = ninside / ntotal if ntotal else 0.0
        summary.append((name, frac, glo, ghi))
        print(" %s: %.1f%% of the binned pixels within the models; G_0 required from %.2f to %.1f"
              % (name, 100.0 * frac, glo, ghi))

        print("\n LaTeX rows for %s:" % name)
        for a, b, npx, mlo, mhi, olo, ohi, fr, g1, g2 in rows:
            if mlo is None:
                print("%4.1f -- %4.1f & %s & no models & %4.1f -- %4.1f & -- \\\\"
                      % (a, b, ("{:,}".format(npx)).replace(",", "\\,"), olo, ohi))
            else:
                print("%4.1f -- %4.1f & %s & %4.1f -- %4.1f & %4.1f -- %4.1f & %.1f{\\%%} \\\\"
                      % (a, b, ("{:,}".format(npx)).replace(",", "\\,"), mlo, mhi, olo, ohi, 100.0 * fr))

    print("\n Summary: cloud, fraction of pixels within the models, smallest and largest G_0 required")
    for name, frac, glo, ghi in summary:
        print(" %-16s %6.1f%%  %8.2f  %8.1f" % (name, 100.0 * frac, glo, ghi))
    print("\n LaTeX rows of the summary:")
    for name, frac, glo, ghi in summary:
        print("%-16s & %.0f{\\%%} & %.2g -- %.3g \\\\" % (name, 100.0 * frac, glo, ghi))


if __name__ == "__main__":
    main()
