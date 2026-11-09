"""make_field_factor.py -- the field factor g of the correction factor for
radiation fields other than G_0 = 10, and its test on the grid of fields.

The relation of hicorrect2 is calibrated for G_0 = 10 in Habing units. The
grid of radiation fields repeats the models in nine fields, from G_0 = 0.98 to
482.5, and harvest_calib3.py has tabulated for each field the factor f that
restores the true surface density, against the surface density that the
reconstruction reports. This script takes, for each field, the pixel-weighted
mean of log f in bins of 0.1 dex of the fitted surface density, smooths it
over three bins, and subtracts the same quantity at G_0 = 10. The result,
exponentiated, is the field factor g(Sigma_f, G_0) that hicorrect2 applies with
the option -g0. Bins with fewer than 20000 pixels take the value of the nearest
populated bin.

The script prints the table, the lines to paste into hicorrect2.py and
hicorrect2.for, and the relative error of the corrected surface density in
four ranges of the fitted surface density for every field: without the
field factor, with the field factor built from all fields, and with the
field factor built without the field being tested (leave one field out).
Above a fitted surface density of 10^22.6 cm^-2 the grid has too few pixels to
determine g, which is therefore held at its value there. The applied factor,
f times g, is made non-decreasing with the fitted surface density, as in
hicorrect2. With a third argument, the harvest table of the subdivision grid
(G_0 = 10), the errors on those independent models are printed as well.

Usage:

    python3 make_field_factor.py <directory of harvest_calib3_G0=*.txt> <directory of hicorrect2.py> [<subdivision harvest>]
"""

import os
import sys
import glob
import numpy as np

DX = 0.1
XC = np.arange(20.95, 23.30, DX)        # centers of the bins of log10 of the fitted surface density
MINPIX = 2.0e4                          # fewest pixels for a bin to be used
XHOLD = 22.6                            # above this log10 of the fitted surface density g is held constant
EDGES = [21.6, 22.0, 22.4, 22.8, 23.4]  # ranges in which the errors are reported


def read(path):
    return np.array([[float(v) for v in l.split()] for l in open(path) if l.strip() and not l.startswith("!")])


def meanlog(d):
    """Pixel-weighted mean of log f in each bin, smoothed over three bins."""
    m = np.full(len(XC), np.nan)
    for i, c in enumerate(XC):
        s = (d[:, 0] >= c - DX / 2) & (d[:, 0] < c + DX / 2)
        if d[s, 3].sum() >= MINPIX:
            m[i] = np.average(np.log(d[s, 5]), weights=d[s, 3])
    o = m.copy()
    for i in range(len(m)):
        w = m[max(0, i - 1):i + 2]
        w = w[np.isfinite(w)]
        o[i] = w.mean() if len(w) else np.nan
    return o


def build(data, fields):
    """log10 of the fields and log g in each bin, for the fields given."""
    ref = meanlog(data[10.0])
    rows = []
    for g in fields:
        r = meanlog(data[g]) - ref
        ok = np.isfinite(r)
        r[~ok] = np.interp(XC[~ok], XC[ok], r[ok])
        r[XC > XHOLD] = np.interp(XHOLD, XC, r)
        rows.append(r)
    return np.log10(np.array(fields)), np.array(rows)


def field_factor(x, g, tab):
    lgs, G = tab
    col = np.array([np.interp(x, XC, r) for r in G])
    t = np.clip(np.log10(g), lgs[0], lgs[-1])
    j = np.clip(np.searchsorted(lgs, t) - 1, 0, len(lgs) - 2)
    a = (t - lgs[j]) / (lgs[j + 1] - lgs[j])
    return np.exp((1.0 - a) * col[j] + a * col[j + 1])


def errors(d, g, fan, tab):
    out = []
    for a, b in zip(EDGES[:-1], EDGES[1:]):
        s = (d[:, 0] >= a) & (d[:, 0] < b)
        if d[s, 3].sum() == 0:
            out.append(np.nan)
            continue
        x, req, w = d[s, 0], d[s, 5], d[s, 3]
        app = fan(x, g, tab)
        out.append(np.exp(np.average(np.log(app / req), weights=w)) - 1.0)
    return out


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return
    sys.path.insert(0, sys.argv[2])
    import hicorrect2 as h
    fa, fadopt = h.build_relation([160.0, 250.0, 350.0, 500.0])
    lfa = np.log10(fa)

    def fan(x, g=10.0, tab=None):
        """Applied factor: f times g, made non-decreasing in the fitted surface density."""
        f = fadopt * (field_factor(lfa, g, tab) if tab is not None else 1.0)
        return np.interp(x, lfa, np.maximum.accumulate(f))

    data = {}
    for p in glob.glob(os.path.join(sys.argv[1], "harvest_calib3_G0=*.txt")):
        g = float(os.path.basename(p)[len("harvest_calib3_G0="):-4].replace("p", "."))
        data[g] = read(p)
    fields = sorted(data)
    tab = build(data, fields)

    print(" field factor g: rows log10 of the fitted surface density (cm^-2), columns G_0 (Habing units)")
    print(" %6s" % "" + "".join("%7g" % g for g in fields))
    for i, c in enumerate(XC):
        print(" %6.2f" % c + "".join("%7.3f" % np.exp(tab[1][j][i]) for j in range(len(fields))))

    print("\n relative error of the corrected surface density (%) in the ranges of log10 Sigma_f",
          ", ".join("%.1f-%.1f" % (a, b) for a, b in zip(EDGES[:-1], EDGES[1:])))
    print(" %6s  %-24s %-24s %s" % ("G_0", "without g", "with g", "leave one field out"))
    for g in fields:
        e0 = errors(data[g], g, fan, None)
        e1 = errors(data[g], g, fan, tab)
        inner = g not in (fields[0], fields[-1], 10.0)
        e2 = errors(data[g], g, fan, build(data, [f for f in fields if f != g])) if inner else [np.nan] * 4
        fmt = lambda e: " ".join("%+5.0f" % (100 * v) for v in e)
        print(" %6g  %-24s %-24s %s" % (g, fmt(e0), fmt(e1), fmt(e2)))

    if len(sys.argv) > 3:
        d = read(sys.argv[3])
        old = lambda x, g=10.0, t=None: np.interp(x, lfa, fadopt)
        print("\n subdivision grid (G_0 = 10): %-24s %s" % ("previous factor", "non-decreasing factor"))
        print("                              %-24s %s" % (" ".join("%+5.0f" % (100 * v) for v in errors(d, 10.0, old, None)),
                                                      " ".join("%+5.0f" % (100 * v) for v in errors(d, 10.0, fan, None))))

    print("\n rows for hicorrect2.py (GTAB):")
    for i, c in enumerate(XC):
        print("    [%s],  # %.2f" % (", ".join("%.3f" % np.exp(tab[1][j][i]) for j in range(len(fields))), c))


if __name__ == "__main__":
    main()
