#!/usr/bin/env python3
"""
make_letter_figures.py -- figures for the Letter on the correction of the
temperature bias of high-resolution surface density images.

Run in the directory holding the Aquila cutout and the outputs of
hires_correct_levels.py produced with --save-diagnostics --out figrun.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from astropy.io import fits
from scipy import ndimage

sys.path.insert(0, "/home/claude/proj/python")
import aa_plotstyle
aa_plotstyle.apply()

PIX = 3.0
MAP = "hi.surface.density.r13p5-nosmooth.fits"
COR = "hi.surface.density.r13p5-nosmooth.corrected.fits"
FAC = "hi.surface.density.r13p5-nosmooth.corrfactor.fits"
TCOR = "hi.temperature.r13p5-nosmooth.corrected.fits"
DIAMS = [27, 40, 61, 91, 137, 205]
from dkbluered import my_cmap as SD_CMAP     # standard scheme for images
PLASMA = plt.get_cmap("plasma", 256)          # line colours only
PLASMA = matplotlib.colors.ListedColormap(PLASMA(np.linspace(0.0, 0.88, 256)))


def get(name):
    return np.array(fits.getdata(name), float).squeeze()


# ----------------------------------------------------------------------
def fig_levels():
    """One-dimensional illustration of the decomposition into levels."""
    m = get(MAP)
    row = m[330, 120:420] / 1e21
    x = np.arange(len(row)) * PIX / 60.0
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.8))
    ax.plot(x, row, color="#d62728", lw=1.2, zorder=5, label="image")
    # a few admissible positions of the flat plate of diameter d = 91 arcsec:
    # each lies entirely beneath the surface, and the opening is the envelope
    # of all such positions rather than the set of contact points
    d_as = 91.0
    half = 0.5 * d_as / 60.0
    for xc in (3.2, 6.6, 11.6):
        k = (x >= xc - half) & (x <= xc + half)
        if k.sum() < 2:
            continue
        h = row[k].min()
        ax.plot([xc - half, xc + half], [h, h], color="#1a9641", lw=1.6,
                solid_capstyle="butt", zorder=6)
    prev = row.copy()
    cols = PLASMA(np.linspace(0.05, 0.75, 4))
    for i, d in enumerate([40, 91, 205, 400]):
        r = max(int(0.5 * d / PIX), 1)
        bg = ndimage.grey_opening(prev, size=2 * r + 1)
        bg = ndimage.uniform_filter(bg, max(int(0.4 * d / PIX), 1))
        bg = np.minimum(bg, prev)
        ax.plot(x, bg, color=cols[i], lw=1.0,
                label=r"$O_{%d^{\prime\prime}}$" % d)
        prev = bg
    ax.set_xlabel("offset along the cut (arcmin)")
    ax.set_ylabel(r"$N_{\rm H_2}$ ($10^{21}$ cm$^{-2}$)")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0, None)
    ax.plot([], [], color="#1a9641", lw=1.6, label=r"plate, $d=91''$")
    ax.legend(loc="upper left", fontsize=7, ncol=1, handlelength=1.3,
              labelspacing=0.25, borderpad=0.3, framealpha=0.85)
    fig.savefig("fig_levels_cut.pdf")
    plt.close(fig)
    print("fig_levels_cut.pdf")


def fig_maps():
    """The original image, the corrected image, and the correction factor."""
    m, c, f = get(MAP), get(COR), get(FAC)
    n = m.shape[0]
    ext = [-n * PIX / 120.0, n * PIX / 120.0] * 2
    fig = plt.figure(figsize=(7.1, 2.55))
    # a narrow spacer follows the first colour bar, so that its tick labels are
    # not covered by the panel beside it
    gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 0.045, 0.125, 1, 0.045],
                          wspace=0.02, left=0.005, right=0.945,
                          bottom=0.01, top=0.90)
    a0 = fig.add_subplot(gs[0, 0])
    a1 = fig.add_subplot(gs[0, 1])
    ca = fig.add_subplot(gs[0, 2])
    a2 = fig.add_subplot(gs[0, 4])
    cb_ax = fig.add_subplot(gs[0, 5])
    vmin, vmax = 3e21, 3e23
    for ax, d, lab in ((a0, m, "original"), (a1, c, "corrected")):
        im = ax.imshow(d, origin="lower", extent=ext, cmap=SD_CMAP,
                       norm=LogNorm(vmin=vmin, vmax=vmax))
        ax.set_title(lab, fontsize=8, pad=3)
    im2 = a2.imshow(f, origin="lower", extent=ext, cmap=SD_CMAP,
                    vmin=1.0, vmax=1.8)
    a2.set_title("correction factor", fontsize=8, pad=3)
    for ax in (a0, a1, a2):
        ax.set_xticks([])
        ax.set_yticks([])
    cb = fig.colorbar(im, cax=ca)
    # the label is placed on the left of this bar, next to the panel it
    # belongs to, since the space beyond its tick labels is otherwise empty
    cb.ax.yaxis.set_label_position("left")
    cb.set_label(r"$N_{\rm H_2}$ (cm$^{-2}$)", fontsize=8, labelpad=-24)
    cb.ax.tick_params(labelsize=7, pad=1)
    cb.ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    cb2 = fig.colorbar(im2, cax=cb_ax)
    cb2.set_label(r"$\mathcal{D}_{\rm c}\,/\,\mathcal{D}$", fontsize=8, labelpad=3)
    cb2.ax.tick_params(labelsize=7, pad=1)
    # the colour bars are made exactly as tall as the images beside them, and
    # each is placed clear of the panel that follows it
    fig.canvas.draw()
    for ax, cax in ((a1, ca), (a2, cb_ax)):
        p_ax, p_cb = ax.get_position(), cax.get_position()
        cax.set_position([p_cb.x0, p_ax.y0, p_cb.width, p_ax.height])
    fig.savefig("fig_aquila_maps.pdf")
    plt.close(fig)
    print("fig_aquila_maps.pdf")


def fig_leveltemp():
    """Temperature of each level, and the resulting factor, against width."""
    sys.path.insert(0, "/home/claude/proj/python")
    from hires_correct_levels import levels as make_levels
    tfit = get("figrun.Tfit.fits")
    m = get(MAP)
    good = np.isfinite(m) & (m > 0)
    filled = np.where(good, m, np.nanmedian(m[good]))
    lev, _ = make_levels(filled, [float(d) for d in DIAMS], PIX)
    ts, cs = [], []
    for k, d in enumerate(DIAMS):
        t = get("figrun.T%d.fits" % d)
        c = get("figrun.C%d.fits" % d)
        # only where the level actually carries material, since elsewhere the
        # temperature is undetermined and the factor is not applied
        w = lev[k]
        sel = good & np.isfinite(t) & np.isfinite(c) \
            & (w > 0.1 * np.percentile(w[good], 99.5))
        ts.append(np.percentile(t[sel], [16, 50, 84]))
        cs.append(np.percentile(c[sel], [16, 50, 84]))
    ts = np.array(ts); cs = np.array(cs)
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.7))
    ax.errorbar(DIAMS, ts[:, 1], yerr=[ts[:, 1] - ts[:, 0], ts[:, 2] - ts[:, 1]],
                fmt="o-", color="k", ms=3.5, lw=1.0, capsize=2,
                label=r"$T_{k}$, level")
    good = np.isfinite(tfit)
    ax.axhline(np.median(tfit[good]), color=PLASMA(0.6), lw=1.2, ls="--",
               label=r"$T_{\rm eff}$, single fit")
    ax.set_xscale("log")
    ax.set_xlabel(r"ball diameter $d_{k}$ (arcsec)")
    ax.set_ylabel(r"dust temperature (K)")
    ax.set_xlim(20, 300)
    aa_plotstyle.clean_log(ax, which="x")
    ax2 = ax.twinx()
    ax2.plot(DIAMS, cs[:, 1], "s:", color="k", mfc="0.6", ms=3.5, lw=1.0)
    ax2.fill_between(DIAMS, cs[:, 0], cs[:, 2], color="0.6", alpha=0.25, lw=0)
    ax2.set_ylabel(r"correction factor $C_{k}$")
    ax2.tick_params(axis="y", direction="in")
    ax.legend(loc="lower right", fontsize=7)
    fig.savefig("fig_level_temperatures.pdf")
    plt.close(fig)
    print("fig_level_temperatures.pdf")


def fig_profile():
    """Radial profile through the brightest core, before and after."""
    m, c = get(MAP), get(COR)
    ok = np.isfinite(m)
    a = np.where(ok, m, -1)
    j, i = np.unravel_index(np.argmax(a), a.shape)
    n = m.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    r = np.hypot(xx - i, yy - j) * PIX
    bins = np.arange(0, 120, PIX)
    pm, pc = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        k = ok & (r >= lo) & (r < hi)
        pm.append(np.median(m[k])); pc.append(np.median(c[k]))
    mid = 0.5 * (bins[:-1] + bins[1:])
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.7))
    ax.plot(mid, np.array(pm) / 1e22, color="k", lw=1.2, label="original")
    ax.plot(mid, np.array(pc) / 1e22, color=PLASMA(0.3), lw=1.2,
            label="corrected")
    ax.set_xlabel(r"radius (arcsec)")
    ax.set_ylabel(r"$N_{\rm H_2}$ ($10^{22}$ cm$^{-2}$)")
    ax.set_yscale("log")
    ax.set_xlim(0, 117)
    aa_plotstyle.clean_log(ax, which="y")
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig("fig_core_profile.pdf")
    plt.close(fig)
    print("fig_core_profile.pdf")


def fig_flowchart():
    """Schematic of the method."""
    fig, ax = plt.subplots(1, 1, figsize=(7.1, 2.2))
    ax.set_xlim(0, 11.9); ax.set_ylim(0, 3.1); ax.axis("off")
    boxes = [
        (0.15, 2.0, 2.1, 0.75, "surface density\n" + r"image $\mathcal{S}$"),
        (0.15, 0.35, 2.1, 0.75, "waveband images\n" +
         r"$\mathcal{I}_\lambda$, common beam"),
        (2.85, 2.0, 2.0, 0.75, "levels of the\nimage, " +
         r"$\mathcal{L}_k(\mathcal{D})$"),
        (2.85, 0.35, 2.0, 0.75, "levels of each\nwaveband, " +
         r"$\mathcal{L}_k^\lambda$"),
        (5.4, 2.0, 1.9, 0.75, "effective\ntemperature " + r"$T_{\rm eff}$"),
        (5.4, 0.35, 1.9, 0.75, "temperature of\neach level, " + r"$T_k$"),
        (7.55, 1.10, 2.05, 1.0,
         "correction factor\n" + r"$C_k=B(T_{\rm eff})\,/\,B(T_k)$"),
        (9.95, 1.10, 1.85, 1.0,
         "corrected image\n" +
         r"$\mathcal{D}_{\rm c}=\sum_k C_k \mathcal{L}_k$"),
    ]
    for x, y, w, h, txt in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                                    lw=0.8, ec="k", fc="0.96"))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center",
                fontsize=7.2)
    arrows = [(2.25, 2.38, 2.85, 2.38), (2.25, 0.72, 2.85, 0.72),
              (4.85, 2.38, 5.4, 2.38), (4.85, 0.72, 5.4, 0.72),
              (7.3, 2.38, 7.85, 2.15), (7.3, 0.72, 7.85, 1.05),
              (9.60, 1.60, 9.95, 1.60)]
    for x0, y0, x1, y1 in arrows:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), lw=0.8,
                                     arrowstyle="-|>", mutation_scale=8,
                                     color="k"))

    fig.savefig("fig_method_flowchart.pdf")
    plt.close(fig)
    print("fig_method_flowchart.pdf")


if __name__ == "__main__":
    fig_flowchart()
    fig_levels()
    fig_maps()
    fig_leveltemp()
    fig_profile()
