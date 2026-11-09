#!/usr/bin/env python3
"""
diagnose_shielding_branches.py -- find the models that form the separate strips in the plot of the dust temperature against
the effective attenuating column (fig_shielding.pdf).

It reads the table written by fit_attenuation_law.py, which holds for every shell the effective attenuating column Sigma_A,
the dust temperature T, the mass weight and the index of the node it belongs to, with the names of the nodes in the header.
For every node it computes the mass-weighted mean departure of its shells from the fitted relation T(Sigma_A) of Eq. (5),
both in kelvin and relative to the relation. The nodes are written to a text file sorted by that departure, so that the
warmest and the coldest families can be read off directly.

When the node names contain the grid indices iNNjNNkNN (i the surface density of the embedding cloud, j the mass of the
Bonnor-Ebert sphere, k its size), the shells are also plotted in three panels colored by each index. A strip made by one
family of models then appears in a single color in one of the panels.

Usage:

    python3 diagnose_shielding_branches.py <attenuation_law.txt> [-out=<prefix>]

Output: <prefix>_nodes.txt (default prefix: shielding_branches) and <prefix>.png.
"""

import re
import sys
import numpy as np


def read_table(path):
    """The shells, the node names and the fitted constants of the table."""
    names, par = {}, {}
    for line in open(path):
        if not line.startswith("!"):
            continue
        m = re.match(r"!\s+node\s+(\d+)\s+(\S+)", line)
        if m:
            names[int(m.group(1))] = m.group(2)
        m = re.match(r"!\s+(Tf|T0|Sstar|gamma)\s*=\s*([0-9.eE+-]+)", line)
        if m:
            par[m.group(1)] = float(m.group(2))
    d = np.loadtxt(path, comments="!")
    if d.ndim != 2 or d.shape[1] < 4:
        sys.exit(" the table has no node column: rerun the updated fit_attenuation_law.py")
    if len(par) < 4:
        sys.exit(" the fitted constants were not found in the header of the table")
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3].astype(int), names, par


def law(s, par):
    return par["Tf"] + (par["T0"] - par["Tf"]) * (1.0 + s / par["Sstar"]) ** (-par["gamma"])


def grid_indices(name):
    """The grid indices i, j and k of a node name, or None when the name does not contain them."""
    m = re.search(r"i(\d+)j(\d+)k(\d+)", name)
    return tuple(int(x) for x in m.groups()) if m else None


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    table = sys.argv[1]
    prefix = "shielding_branches"
    for a in sys.argv[2:]:
        if a.startswith("-out="):
            prefix = a[5:]
    s, t, w, node, names, par = read_table(table)
    dt = t - law(s, par)

    rows = []
    for k in np.unique(node):
        m = node == k
        ww = w[m] / w[m].sum()
        rows.append((float(np.sum(ww * dt[m])), float(np.sum(ww * dt[m] / law(s[m], par))),
                     float(np.sum(ww * np.abs(dt[m]))), int(m.sum()), float(s[m].min()), float(s[m].max()),
                     k, names.get(k, "?")))
    rows.sort(key=lambda r: -r[0])
    with open(prefix + "_nodes.txt", "w") as fh:
        fh.write("! Mass-weighted mean departure of the shells of every node from the fitted relation T(Sigma_A),\n")
        fh.write("! sorted from the warmest to the coldest node. dT in K, dT/T relative, |dT| the mean absolute departure.\n")
        fh.write("!      dT      dT/T     |dT|  shells   Sigma_A min   Sigma_A max   node  name\n")
        for r in rows:
            fh.write(" %8.3f %8.4f %8.3f %7d %13.4e %13.4e %6d  %s\n" % r)
    print(" written %s_nodes.txt" % prefix)
    print(" the ten warmest nodes:")
    for r in rows[:10]:
        print("   dT = %+6.2f K  %s" % (r[0], r[7]))

    idx = np.array([grid_indices(names.get(k, "")) or (-1, -1, -1) for k in node])
    if not (idx[:, 0] >= 0).any():
        print(" the node names contain no grid indices iNNjNNkNN; no figure made")
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    good = idx[:, 0] >= 0
    order = np.argsort(w[good])
    sx = np.geomspace(1.0e19, 1.0e24, 300)
    fig, ax = plt.subplots(1, 3, figsize=(15.0, 4.6), sharey=True)
    for p, label in enumerate(("i (embedding cloud)", "j (mass)", "k (size)")):
        c = idx[good, p][order]
        sc = ax[p].scatter(s[good][order], t[good][order], c=c, s=1.5, cmap="turbo", lw=0.0, rasterized=True)
        ax[p].plot(sx, law(sx, par), "k--", lw=1.0)
        ax[p].set_xscale("log")
        ax[p].set_xlim(2.0e19, 4.0e23)
        ax[p].set_ylim(6.0, 23.0)
        ax[p].set_xlabel(r"Effective attenuating column $\Sigma_{\rm A}$ (cm$^{-2}$)")
        ax[p].set_title("colored by grid index " + label)
        fig.colorbar(sc, ax=ax[p], pad=0.01)
    ax[0].set_ylabel("Dust temperature $T$ (K)")
    fig.tight_layout()
    fig.savefig(prefix + ".png", dpi=150)
    print(" written %s.png" % prefix)


if __name__ == "__main__":
    main()
