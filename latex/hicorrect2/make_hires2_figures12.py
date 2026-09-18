#!/usr/bin/env python3
"""
make_hires2_figures.py -- figures of the paper on the calibrated correction of
the hires surface density images.

The figures produced are:

  fig_method.pdf         a schematic of the method, from the grid of models to
                         the corrected images
  fig_superposition.pdf  what the correction assumes about a line of sight that
                         crosses more than one structure, and what that costs
  fig_deficit.pdf        the reconstruction against the truth at the coarsest
                         and the finest resolution, showing that the deficit is
                         carried by the base image alone
  fig_shielding.pdf      the dust temperature against the effective attenuating
                         column, from the grid of radiative transfer models,
                         with the relation fitted to it
  fig_geometry.pdf       the directions the heating radiation arrives from in
                         the three limiting shapes, and the attenuating column
                         and dust temperature along a line of sight
  fig_angular.pdf        how the angular average collapses a distribution of
                         directions into a single effective column, and how
                         little the shape assumed matters
  fig_mixing.pdf         why a single modified blackbody fitted to a
                         non-isothermal line of sight returns too little
                         surface density, and the residual pattern it leaves
  fig_relation.pdf       the correction factor against the surface density of
                         the uncorrected image, for the two limiting shapes and
                         their average
  fig_simsky.pdf         the simulated sky before and after correction, with
                         the true surface density, in three panels
  fig_core_recovery.pdf  the recovered fraction of the integrated excess of the
                         injected cores against their true mass, before and
                         after correction
  fig_aquila.pdf
  fig_residual.pdf
  fig_masses.pdf

fig_orientation.pdf and the routine that serves it belong to a superseded
scheme and are retained only until it is decided what replaces them.

Run in the directory holding the simulated sky, the images corrected by
hicorrect2, and attenuation_law.txt.
"""
import os
import sys

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from matplotlib import font_manager
import matplotlib.pyplot as plt

font_manager.fontManager.addfont('/Users/amenshch/Library/Fonts/Helvetica.ttf')
font_manager.fontManager.addfont('/Users/amenshch/Library/Fonts/Helvetica-Oblique.ttf')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Helvetica'
plt.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
plt.rcParams['mathtext.sf'] = 'Helvetica'   

sys.path.insert(0, "/mnt/project")
try:
    import aa_plotstyle
    from dkbluered import my_cmap as SD_CMAP
    aa_plotstyle.apply()
except ImportError:                                   # outside the project
    SD_CMAP = plt.get_cmap("magma")
    matplotlib.rcParams.update({"font.size": 8, "axes.linewidth": 0.5, "xtick.direction": "in", "ytick.direction": "in"})

plt.rcParams['lines.linewidth'] = 0.9   # default is 1.5   
plt.rcParams['axes.linewidth'] = 0.5      # default is 0.8   
plt.rcParams['xtick.major.width'] = 0.4   # default 0.8
plt.rcParams['ytick.major.width'] = 0.4
plt.rcParams['xtick.minor.width'] = 0.4
plt.rcParams['ytick.minor.width'] = 0.4

PLASMA = plt.get_cmap("plasma", 256)
PIX = 3.0
BEAM, FINE = 36.3, 13.5
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
TFLOOR, TSURF, SSTAR, GAMMA = 6.2090, 22.4463, 2.07379e21, 1.14100
KAPHEAT = 1.0 / (0.94e21 * 1.086)    # cm^2 per H2 at the wavelengths that heat the dust
NANG = 160                           # directions used in the angular average
CAL_ERR = {70.0: 0.20, 100.0: 0.20, 160.0: 0.20, 250.0: 0.10, 350.0: 0.10, 500.0: 0.10}
H, K, C = 6.62607015e-27, 1.380649e-16, 2.99792458e10
KAPPA0, DTG, MU, AMU = 10.0, 1.0e-2, 2.8, 1.66053906660e-24
BANDS = [160.0, 250.0, 350.0, 500.0]

g = lambda f: np.array(fits.getdata(f), float).squeeze()


def planck(t, w):
    nu = C / (w * 1.0e-4)
    x = np.clip(H * nu / (K * np.asarray(t, float)), 1e-8, 500.0)
    return 2.0 * H * nu ** 3 / C ** 2 / np.expm1(x)


def cnu(w):
    nu = C / (w * 1.0e-4)
    return 1.0e17 * KAPPA0 * (nu / 1.0e12) ** 2 * DTG * MU * AMU


def law(s):
    return TFLOOR + (TSURF - TFLOOR) * (1.0 + np.asarray(s, float)
                                        / SSTAR) ** (-GAMMA)


def attenuating_column(sigma, geometry="mean", ndep=512):
    """Effective attenuating column at each depth through a uniform structure.

    The angular average is taken over NANG directions, for a sphere sampled along a diameter, for a plane-parallel layer sampled
    along its normal, for the average of the two, or for a uniform cube sampled along the line joining two opposite faces.  The
    first three reproduce attenuating_column() of hicorrect2 exactly; the cube is used only for the comparison of Table 5.
    """
    t0 = KAPHEAT * sigma / 2.0
    z = (np.arange(ndep) + 0.5) / ndep * 2.0 - 1.0
    ps = (np.arange(NANG) + 0.5) / NANG * np.pi
    wps = np.sin(ps)
    lsph = (-np.abs(z)[:, None] * np.cos(ps)
            + np.sqrt(np.maximum(1.0 - (np.abs(z)[:, None] * np.sin(ps)) ** 2, 0.0)))
    jsph = (np.exp(-t0 * lsph) * wps).sum(axis=1) / wps.sum()
    mu = (np.arange(NANG) + 0.5) / NANG * 2.0 - 1.0
    dlay = np.where(mu > 0.0, (1.0 - z[:, None]) / np.abs(mu), (1.0 + z[:, None]) / np.abs(mu))
    jlay = np.exp(-t0 * dlay).mean(axis=1)
    ssph = -np.log(np.maximum(jsph, 1e-300)) / max(t0, 1e-12)
    slay = -np.log(np.maximum(jlay, 1e-300)) / max(t0, 1e-12)
    if geometry == "cube":
        # a uniform cube, sampled along the line through the centers of two opposite faces, with the directions drawn from a
        # spiral that covers the sphere evenly, since a cube is not symmetric about the line of sight
        nd = 4 * NANG
        k = np.arange(nd) + 0.5
        cz = 1.0 - 2.0 * k / nd
        cr = np.sqrt(np.maximum(1.0 - cz ** 2, 0.0))
        phi = np.pi * (1.0 + 5.0 ** 0.5) * k
        nx, ny, nz = cr * np.cos(phi), cr * np.sin(phi), cz
        far = np.where(nz > 0.0, (1.0 - z[:, None]) / np.maximum(nz, 1e-12), (1.0 + z[:, None]) / np.maximum(-nz, 1e-12))
        lcub = np.minimum(np.minimum(1.0 / np.maximum(np.abs(nx), 1e-12), 1.0 / np.maximum(np.abs(ny), 1e-12)), far)
        jcub = np.exp(-t0 * lcub).mean(axis=1)
        return -np.log(np.maximum(jcub, 1e-300)) / max(t0, 1e-12) * sigma / 2.0
    if geometry == "short":
        u = (np.arange(len(z)) + 0.5) / len(z)
        return np.minimum(u, 1.0 - u) * sigma
    use = {"sphere": ssph, "layer": slay, "mean": 0.5 * (ssph + slay)}[geometry]
    return use * sigma / 2.0


def relation(geometry="mean", n=90, ndep=512):
    """Fitted surface density and the factor that restores the true one.

    The wavebands are weighted by their calibration uncertainties, which is what the reconstruction does.
    """
    tg = np.geomspace(0.98 * TFLOOR, 1.02 * TSURF, 2000)
    model = np.array([[cnu(w) * planck(t, w) for w in BANDS] for t in tg])
    err = np.array([CAL_ERR.get(w, 0.20) for w in BANDS])
    sig = np.geomspace(1.0e20, 1.0e24, n)
    fit = np.empty(n)
    for i, s in enumerate(sig):
        t = law(attenuating_column(s, geometry, ndep))
        obs = np.array([cnu(w) * s * np.mean(planck(t, w)) for w in BANDS])
        wg = 1.0 / (err * obs) ** 2
        amp = ((wg * obs * model).sum(axis=1) / (wg * model * model).sum(axis=1))
        j = int(np.argmin((wg * (obs - amp[:, None] * model) ** 2).sum(axis=1)))
        fit[i] = (cnu(500.0) * s * np.mean(planck(t, 500.0))
                  / (cnu(500.0) * planck(tg[j], 500.0)))
    return fit, sig / fit


def fig_superposition(out="fig_superposition.pdf"):
    """What the method assumes about a pixel that receives emission from more than one structure.

    Every pixel collects the emission of everything along its line of sight, and one modified blackbody is fitted to that total,
    so what biases the fit is the whole distribution of temperatures along the line of sight.  The method replaces that
    distribution by the one a single uniform body of the same total column would have.

    Left: the dust temperature along a line of sight of 2e22 cm^-2, as the method assumes it and as it would be if half the column
    were a diffuse cirrus and half a sphere embedded in it.  The cirrus is warm because it is shielded by its own column only.
    Middle: the emission that follows, with the single modified blackbody fitted to each, and the surface density each fit
    returns.  Right: the correction factor required, the curve the method applies and the range spanned by embedding a sphere in a
    cirrus of 1e21, 3e21 and 1e22 cm^-2 at every position within it.  All of it is computed by superposition_test.py, which is
    described in the appendix; nothing here is drawn by hand.
    """
    import superposition_test as sup

    one, sep = PLASMA(120), PLASMA(30)
    fig, axx = plt.subplots(2, 2, figsize=(7.1, 5.0))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.095, top=0.930, wspace=0.30, hspace=0.42)
    ax = [axx[0, 0], axx[0, 1], axx[1, 0]]
    axq = axx[1, 1]

    sig, sig_c = 2.0e22, 1.0e22
    sig_s, q = sig - sig_c, 0.0        # the sphere against the near side of the cirrus, the arrangement that departs most

    # ---- left: the temperature along the line of sight, in the column coordinate measured from the observer's side
    u = (np.arange(400) + 0.5) / 400.0
    sc = 1.0e22
    ax[0].plot(u * sig / sc, law(attenuating_column(sig, "mean", 400)), "-", lw=1.5, color=one,
               label="one uniform body,\nas the method assumes")
    tc = sup.layer_temperatures(sig_c)
    ts = sup.sphere_temperatures(sig_s, sig_c, q)
    nc = len(tc) // 2
    xc1 = (np.arange(nc) + 0.5) / nc * q * sig_c
    xs = q * sig_c + (np.arange(len(ts)) + 0.5) / len(ts) * sig_s
    xc2 = q * sig_c + sig_s + (np.arange(nc) + 0.5) / nc * (1.0 - q) * sig_c
    ax[0].plot(xc1 / sc, tc[:nc], "-", lw=1.5, color=sep, label="cirrus with an embedded sphere")
    ax[0].plot(xs / sc, ts[::-1], "-", lw=1.5, color=sep)
    ax[0].plot(xc2 / sc, tc[nc:], "-", lw=1.5, color=sep)
    ax[0].axvspan(q * sig_c / sc, (q * sig_c + sig_s) / sc, color="0.90", zorder=0)
    ax[0].text((q * sig_c + 0.5 * sig_s) / sc, 21.2, "sphere", fontsize=6.4, color="0.45", ha="center")
    ax[0].set_xlim(0.0, sig / sc)
    ax[0].set_ylim(4.0, 23.5)
    ax[0].set_xlabel(r"column from the near side ($10^{22}$ cm$^{-2}$)")
    ax[0].set_ylabel(r"$T$ (K)")
    ax[0].legend(frameon=False, fontsize=6.2, loc="lower center", handlelength=1.6, labelspacing=0.3)
    ax[0].annotate("", xy=(0.02 * sig / sc, 24.7), xytext=(0.34 * sig / sc, 24.7), annotation_clip=False,
                   arrowprops=dict(arrowstyle="-|>", lw=1.0, color="0.25"))
    ax[0].text(0.38 * sig / sc, 24.7, "to the observer", fontsize=6.4, color="0.30", va="center")

    # ---- middle: the emission each produces, and what one blackbody fitted to it returns
    wl = np.geomspace(130.0, 620.0, 200)
    tg = np.geomspace(0.98 * TFLOOR, 1.02 * TSURF, 2000)
    for label, temps, wgts, col in (("one uniform body", law(attenuating_column(sig, "mean", 400)),
                                     np.full(400, sig / 400.0), one),
                                    ("cirrus and sphere", np.concatenate([tc, ts]),
                                     np.concatenate([np.full(len(tc), sig_c / len(tc)),
                                                     np.full(len(ts), sig_s / len(ts))]), sep)):
        obs = np.array([cnu(w) * np.sum(wgts * planck(temps, w)) for w in BANDS])
        sf = sup.fitted_column(obs)
        ax[1].plot(wl, [cnu(w) * np.sum(wgts * planck(temps, w)) for w in wl], "-", lw=1.5, color=col,
                   label=r"%s, $\Sigma_{\rm f} = %.2f$" % (label, sf / 1.0e22))
        ax[1].plot(BANDS, obs, "o", ms=3.4, mfc="none", mec=col, mew=1.0)
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlim(130.0, 620.0)
    ax[1].set_xticks(BANDS)
    ax[1].set_xticklabels(["160", "250", "350", "500"])
    ax[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax[1].set_xlabel(r"wavelength ($\mu$m)")
    ax[1].set_ylabel(r"$I_\nu$ (MJy sr$^{-1}$)")
    ax[1].legend(frameon=False, fontsize=6.0, loc="lower center", handlelength=1.6, labelspacing=0.3)
    ax[1].set_title(r"the same total column, $\Sigma = 2$", fontsize=6.8, pad=4.0)
    ax[1].text(0.02, 0.98, r"all in $10^{22}$ cm$^{-2}$", fontsize=6.0, color="0.45", ha="left", va="top",
               transform=ax[1].transAxes)

    # ---- right: the factor the method applies, and the range the arrangement spans
    sg = np.geomspace(6.0e21, 1.2e23, 22)
    base = np.array([sup.one_body(v) / 1.0 for v in sg])
    xs_f = np.array([v / sup.one_body(v) for v in sg])
    lo, hi = np.full(len(sg), np.inf), np.full(len(sg), -np.inf)
    for scir in (1.0e21, 3.0e21, 1.0e22):
        for qq in (0.0, 0.25, 0.5):
            v = np.array([sup.two_component(x, scir, qq) if x > 1.6 * scir else np.nan for x in sg])
            lo, hi = np.fmin(lo, v), np.fmax(hi, v)
    ok = np.isfinite(lo) & np.isfinite(hi)
    ax[2].fill_between(xs_f[ok], lo[ok], hi[ok], color=sep, alpha=0.22, lw=0.0,
                       label="range over the arrangement")
    ax[2].plot(xs_f, base, "-", lw=1.5, color=one, label="the factor the method applies")
    ax[2].set_xscale("log")
    ax[2].set_xlim(5.0e21, 6.0e22)
    ax[2].set_ylim(1.0, 3.0)
    ax[2].set_xlabel(r"uncorrected $N_{\rm H_2}$ (cm$^{-2}$)")
    ax[2].set_ylabel(r"correction factor $f$")
    ax[2].legend(frameon=False, fontsize=6.2, loc="upper left", handlelength=1.6, labelspacing=0.3)

    # ---- the factor against where the structure sits in the cirrus
    qq = np.linspace(0.0, 1.0, 21)
    for scir, col in ((1.0e21, PLASMA(20)), (3.0e21, PLASMA(110)), (1.0e22, PLASMA(200))):
        for tot, ls in ((2.0e22, "-"), (5.0e22, "--")):
            axq.plot(qq, [sup.two_component(tot, scir, v) for v in qq], ls, lw=1.3, color=col)
        axq.plot([], [], "-", lw=1.3, color=col,
                 label=r"$\Sigma_{\rm c} = 10^{%.2f}$" % np.log10(scir))
    for tot, ls in ((2.0e22, "-"), (5.0e22, "--")):
        axq.axhline(sup.one_body(tot), ls=":", lw=0.9, color="0.4")
    axq.set_xlim(0.0, 1.0)
    axq.set_xlabel(r"fraction $q$ of the cirrus in front of the structure")
    axq.set_ylabel(r"correction factor $f$")
    axq.legend(frameon=False, fontsize=6.2, loc="lower center", handlelength=1.6, labelspacing=0.25, ncol=3,
               columnspacing=0.8)
    axq.set_title(r"solid $\Sigma = 2\times10^{22}$, dashed $5\times10^{22}$ cm$^{-2}$;"
                  "\ndotted: the factor the method applies", fontsize=6.4, pad=3)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_deficit(stem="simsky5.w40", out="fig_deficit.pdf"):
    """Where the deficit sits: in the base image, not in the increments.

    The reconstruction is compared with the true surface density of the simulated sky, bin by bin in the true surface density, at
    the coarsest resolution and at the finest.  The two curves lie almost on top of one another, which means that the increments
    the reconstruction adds between those resolutions are very nearly right and that the whole of the deficit is carried by the
    base image.  That is what makes a single additive correction of the base image enough.
    """
    cols = (PLASMA(30), PLASMA(120), PLASMA(210))
    T13 = g("%s.surfdens.r13p5.cores.fits" % stem)
    D13 = g("hi.surface.density.r13p5.fits")
    D36 = g("hi.surface.density.r36p3.fits")
    sig = np.sqrt(BEAM ** 2 - FINE ** 2) * FWHM_TO_SIGMA / PIX
    T36 = gaussian_filter(T13, sig, mode="nearest")
    ok = (D36 > 0.0) & (D13 > 0.0) & np.isfinite(T13)

    edges = np.geomspace(3.0e21, 6.0e22, 19)
    mid = np.sqrt(edges[:-1] * edges[1:])
    r36, r13, inc = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = ok & (T36 >= lo) & (T36 < hi)
        if m.sum() < 200:
            r36.append(np.nan)
            r13.append(np.nan)
            inc.append(np.nan)
            continue
        r36.append(np.median(D36[m] / T36[m]))
        r13.append(np.median(D13[m] / T13[m]))
        inc.append(np.median((D13[m] - D36[m]) / np.maximum(T13[m] - T36[m], 1.0)))
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.7))
    ax.plot(mid, r36, "-", color=cols[0], label="%.1f″ base image" % BEAM)
    ax.plot(mid, r13, "--", color=cols[1], label="%.1f″ reconstruction" % FINE)
    ax.axhline(1.0, ls=":", lw=0.5, color="0.3")
    ax.set_xscale("log")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0.45, 1.05)
    ax.set_xlabel(r"true $N_{\rm H_2}$ (cm$^{-2}$)")
    ax.set_ylabel("reconstructed / true")
    ax.legend(frameon=False, fontsize=6.8, loc="lower center", handlelength=1.8, borderpad=0.9)
    fig.tight_layout(pad=0.4)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_shielding(table="attenuation_law.txt", out="fig_shielding.pdf"):
    """The measured relation between the effective attenuating column and the dust temperature.

    The table is the one written by fit_attenuation_law.py: the effective attenuating column of a shell in cm^-2, the dust
    temperature of that shell in K, and its mass weight.
    """
    ls, te = [], []
    for line in open(table):
        if line.lstrip().startswith("!") or not line.strip():
            continue
        c = line.split()
        ls.append(np.log10(float(c[0])))
        te.append(float(c[1]))
    ls, te = np.array(ls), np.array(te)
    s = np.geomspace(1.0e18, 1.0e23, 400)
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.7))
    ax.plot(10.0 ** ls, te, "o", ms=3.0, mfc="none", mec=PLASMA(40),
            mew=0.9, label="model grid")
    ax.plot(s, law(s), "-", lw=1.2, color=PLASMA(180), label="Eq. (6)")
    ax.axhline(TFLOOR, ls=":", lw=0.8, color="0.5")
    ax.text(1.3e18, TFLOOR + 0.4, r"$T_{\rm f}$", fontsize=7, color="0.4")
    ax.set_xscale("log")
    ax.set_xlim(1.0e18, 1.0e23)
    ax.set_ylim(5.0, 24.0)
    ax.set_xlabel(r"shielding column $s$ (cm$^{-2}$)")
    ax.set_ylabel(r"dust temperature $T$ (K)")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    fig.tight_layout(pad=0.4)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def _arrows(ax, paths, angles, r_out, r_in, col, tau_ref):
    """Draw the incoming radiation as arrows whose width follows how much of it survives the journey inward.

    paths holds the column crossed on the way in along each direction, in units of the half thickness of the structure, and
    tau_ref is the optical depth of that half thickness.  The arrow is drawn thick and opaque where little material has been
    crossed and thin and pale where the radiation is extinguished.
    """
    tr = np.exp(-tau_ref * np.asarray(paths, float))
    tr = tr / tr.max()
    for ang, t in zip(angles, tr):
        dx, dy = np.cos(ang), np.sin(ang)
        ax.annotate("", xy=(r_in * dx, r_in * dy), xytext=(r_out * dx, r_out * dy),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=0.25 + 2.4 * t,
                                    mutation_scale=3.0 + 8.0 * t, shrinkA=0.0, shrinkB=0.0,
                                    alpha=0.30 + 0.70 * t))


def fig_geometry(out="fig_geometry.pdf"):
    import matplotlib.gridspec as gridspec
    """Where the heating radiation comes from, and what it implies along a line of sight.

    Upper panels: a parcel of dust, the black point, at mid-depth in a plane-parallel layer, at the center of a sphere, and on the
    axis of a filament.  The arrows are the incoming radiation, drawn thick where it has crossed little material, so that the eye
    sees which directions actually deliver energy.  Lower panels: the effective attenuating column of Eq. (6) and the dust
    temperature that follows from it, along a line of sight through a layer, for three total surface densities.
    """
    cols = (PLASMA(30), PLASMA(120), PLASMA(210))
    fig = plt.figure(figsize=(7.1, 4.75))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.9, 1.0], hspace=0.15, wspace=0.0, left=0.1, right=0.985, bottom=0.095, top=0.945)
    ax = [fig.add_subplot(gs[0, i]) for i in range(3)]
    na, tau_ref = 36, 1.3
    ang = (np.arange(na) + 0.5) / na * 2.0 * np.pi
    sa = np.maximum(np.abs(np.sin(ang)), 1.0e-3)
    ca = np.maximum(np.abs(np.cos(ang)), 1.0e-3)

    # ---- plane-parallel layer: the body runs past both edges, the escape is short only along the normal
    ax[0].add_patch(plt.Rectangle((-1.30, -0.46), 2.60, 0.92, fc="0.8", ec="0.55", lw=0.5, zorder=0))
    _arrows(ax[0], np.minimum(1.0 / sa, 40.0), ang, 1.00, 0.13, cols[0], tau_ref)
    ax[0].set_title("Plane-parallel layer", fontsize=7.6, pad=4.0)

    # ---- sphere: every direction crosses the same column, none is favored
    ax[1].add_patch(plt.Circle((0.0, 0.0), 0.46, fc="0.8", ec="0.55", lw=0.5, zorder=0))
    _arrows(ax[1], np.ones(na), ang, 1.00, 0.13, cols[1], tau_ref)
    ax[1].set_title("Sphere", fontsize=7.6, pad=4.0)

    # ---- filament: short across the axis, long along it, so it lies between the other two
    ax[2].add_patch(plt.Rectangle((-1.30, -0.30), 2.60, 0.60, fc="0.8", ec="0.55", lw=0.5, zorder=0))
    _arrows(ax[2], np.minimum(np.minimum(1.0 / sa, 3.2 / ca), 40.0), ang, 1.00, 0.13, cols[2], tau_ref)
    ax[2].set_title("Filament", fontsize=7.6, pad=4.0)

    for a in ax:
        a.plot([0.0], [0.0], "o", ms=3.4, color="k", zorder=5)
        a.set_xlim(-1.06, 1.06)
        a.set_ylim(-1.06, 1.06)
        a.set_aspect("equal")
        a.set_xticks([])
        a.set_yticks([])
        for sp in a.spines.values():
            sp.set_visible(False)
    
    # ---- along a line of sight through a layer
    # bottom row: 2 equal-width panels spanning the full width
    gs_bottom = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, :], wspace=0.3)
    b1 = fig.add_subplot(gs_bottom[0, 0])
    b2 = fig.add_subplot(gs_bottom[0, 1])   
##    b1 = fig.add_subplot(gs[1, 0:2])
##    b2 = fig.add_subplot(gs[1, 2])
    u = (np.arange(512) + 0.5) / 512.0
    for sig, col, lab in ((1.0e21, cols[0], "10^{21}"), (1.0e22, cols[1], "10^{22}"), (1.0e23, cols[2], "10^{23}")):
        sa_u = attenuating_column(sig, "layer")
        b1.plot(u, sa_u, "-", color=col, label=r"$\Sigma = %s$ cm$^{-2}$" % lab)
        b1.plot(u, np.minimum(u, 1.0 - u) * sig, ":", color=col)
        b2.plot(u, law(sa_u), "-", color=col)
    b1.plot([], [], ":", color="0.3", label="shortest path")
    b1.set_yscale("log")
    b1.set_xlim(0.0, 1.0)
    b1.set_ylim(1.0e19, 3.0e23)
    b1.annotate("", xy=(0.02, 0.92), xytext=(0.13, 0.92), xycoords=b1.transAxes, textcoords=b1.transAxes,
                annotation_clip=False, arrowprops=dict(arrowstyle="-|>", color="0.25", lw=0.5))
    b1.text(0.14, 0.92, "to observer", fontsize=7.0, color="0.20", va="center", transform=b1.transAxes)
    b1.set_xlabel("Fractional depth along line of sight")
    b1.set_ylabel(r"$\Sigma_{\rm a}$ (cm$^{-2}$)", labelpad=0)
    b1.legend(frameon=False, fontsize=7, loc="lower center", ncol=2, handlelength=1.6, columnspacing=1.0, 
              borderpad=0.9, labelspacing=0.25)
    b2.axhline(TFLOOR, ls=":", color="0.3")
    b2.text(0.04, TFLOOR + 0.3, r"$T_{\rm f}$", fontsize=8, color="0.3")
    b2.set_xlim(0.0, 1.0)
    b2.set_ylim(4.5, 20.0)
    b2.set_xlabel("Fractional depth")
    b2.set_ylabel(r"$T$ (K)", labelpad=0)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)

def _bar(ax, y, h, txt, fc, fs=7.8, w=0.72):
    """One bar of the flowchart, centered, with its sentence inside it."""
    ax.add_patch(FancyBboxPatch((0.5 - 0.5 * w, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.010",
                                lw=0.5, ec="0.45", fc=fc, transform=ax.transAxes, clip_on=False))
    ax.text(0.5, y + 0.5 * h, txt, ha="center", va="center", fontsize=fs, linespacing=1.4, transform=ax.transAxes)


def fig_method(out="fig_method.pdf"):
    """The method step by step, in the linear form of the flowchart of the getsf paper.

    The first group is the calibration, done once, which yields the temperature relation.  The second is the calculation of
    Sect. 4.4, repeated for every trial surface density, which yields the table of the correction factor against the surface
    density that an observed image actually holds.  The third is what is then done to an image.
    """
    groups = [
        (PLASMA(205), [
            "Grid of radiative transfer models: dense cores embedded in clouds",
            r"Take the dust temperature $T$ and the attenuating column $\Sigma_{\rm a}$ of every shell",
            r"Fit $T\,(\Sigma_{\rm a}) = T_{\rm f} + (T_0 - T_{\rm f})\,(1 + \Sigma_{\rm a}\,/\,\Sigma_{\bigstar})^{-\gamma}$, "
            "Eq. (7)"]),
        (PLASMA(120), [
            r"Trial line of sight of surface density $\Sigma$, one uniform structure",
            r"Attenuating column $\Sigma_{\rm a}(u)$ at every depth, averaged over directions, Eq. (6)",
            r"Temperature $T\,(\Sigma_{\rm a}(u))$ at every depth, from the fitted relation",
            r"Emission of all the depths added, giving the intensities $I_\lambda(\Sigma)$, Eq. (8)",
            r"One modified blackbody fitted as the image was, giving $\Sigma_{\rm f}$, Eq. (9)",
            r"Tabulate $f\,(\Sigma_{\rm f}) = \Sigma\,/\,\Sigma_{\rm f}$, what an image holds, Eq. (10)"]),
        (PLASMA(35), [
            r"Interpolate $f$ onto the image $\mathcal{D}$ at the coarsest resolution",
            r"Add $(f - 1)\,\mathcal{D}$ to that image and to every finer one, Eq. (12)"])]

    nbar = sum(len(t) for _, t in groups)
    ngap, ngrp = nbar - 1, len(groups) - 1
    fig, ax = plt.subplots(1, 1, figsize=(7.1, 0.435 * nbar + 0.30))
    fig.subplots_adjust(left=0.004, right=0.996, bottom=0.008, top=0.992)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    gap, extra = 0.54, 0.30                     # the gap between bars, and the additional gap between groups, in units of h
    h = 0.75 / (nbar + gap * ngap + extra * ngrp)
    y = 0.875
    for gi, (col, texts) in enumerate(groups):
        for k, t in enumerate(texts):
            y -= h
            _bar(ax, y, h, t, (col[0], col[1], col[2], 0.50 if k == 0 else 0.24))
            last = (gi == len(groups) - 1) and (k == len(texts) - 1)
            if not last:
                g = gap * h + (extra * h if k == len(texts) - 1 else 0.0)
                ax.annotate("", xy=(0.5, y - g + 0.004), xytext=(0.5, y - 0.005), xycoords=ax.transAxes, textcoords=ax.transAxes,
                            arrowprops=dict(arrowstyle="-|>", lw=0.8, color=(col[0], col[1], col[2], 0.85),
                                            mutation_scale=9, shrinkA=0.0, shrinkB=0.0))
                y -= g

#    fig.savefig(out)
    fig.savefig(out, bbox_inches='tight', pad_inches=-0.6)
    plt.close(fig)
    print("written %s" % out)


def fig_angular(out="fig_angular.pdf"):
    """How the angular average collapses a whole distribution of directions into one column.

    Left: the fraction of the heating radiation that survives the journey inward, against the cosine of the angle to the normal,
    for a parcel at mid-depth of a plane-parallel layer of three total surface densities.  The mean of each curve is the mean
    intensity of Eq. (5), marked by the dashed line, and the single column that would produce that same mean intensity is the
    effective attenuating column of Eq. (6).  Right: that column against the total surface density of the line of sight, for the
    two limiting shapes and for the average of them that is adopted, with the shortest path for comparison.  The three lie close
    together, which is why the shape assumed matters so little.
    """
    cols = (PLASMA(30), PLASMA(120), PLASMA(210))
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.85))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.175, top=0.955, wspace=0.28)

    mu = np.linspace(0.005, 1.0, 400)
    for sig, col, lab in ((3.0e20, cols[0], r"3\times10^{20}"), (1.0e21, cols[1], "10^{21}"),
                          (3.0e21, cols[2], r"3\times10^{21}")):
        tr = np.exp(-KAPHEAT * 0.5 * sig / mu)
        ax[0].plot(mu, tr, "-", lw=1.4, color=col, label=r"$\Sigma = %s$ cm$^{-2}$" % lab)
        ax[0].axhline(np.trapezoid(tr, mu), ls="--", lw=0.8, color=col)
    ax[0].set_xlim(0.0, 1.0)
    ax[0].set_ylim(0.0, 1.0)
    ax[0].set_xlabel(r"$\cos i$, the direction the radiation comes from")
    ax[0].set_ylabel(r"$e^{-\kappa\Sigma_{\rm a}(\hat n)}$, the fraction that survives")
    ax[0].legend(frameon=False, fontsize=6.8, loc="upper left", handlelength=1.7)
    ax[0].text(0.97, 0.06, "dashed: the angular mean of each curve,\nwhich is what defines "
               r"$\Sigma_{\rm a}$", fontsize=6.6, color="0.35", ha="right", va="bottom",
               transform=ax[0].transAxes)

    sg = np.geomspace(1.0e20, 1.0e24, 60)
    mid = 512 // 2
    for geom, col, lab in (("layer", cols[0], "plane-parallel layer"), ("mean", cols[1], "adopted, the average"),
                           ("sphere", cols[2], "sphere")):
        y = np.array([attenuating_column(v, geom)[mid] for v in sg])
        ax[1].plot(sg, y, "-", lw=1.4 if geom == "mean" else 1.0, color=col, label=lab)
    ax[1].plot(sg, 0.5 * sg, ":", lw=1.0, color="0.4", label="shortest path")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlim(1.0e20, 1.0e24)
    ax[1].set_xlabel(r"$\Sigma$ of the line of sight (cm$^{-2}$)")
    ax[1].set_ylabel(r"$\Sigma_{\rm a}$ at mid-depth (cm$^{-2}$)")
    ax[1].legend(frameon=False, fontsize=6.8, loc="upper left", handlelength=1.7)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_mixing(out="fig_mixing.pdf"):
    """Why one modified blackbody fitted to a non-isothermal line of sight returns too little surface density.

    Left: the spectral energy distribution of a line of sight through a layer of 3e22 cm^-2, the sum over depth of the emission of
    material at every temperature, together with the best single modified blackbody fitted to it in the four wavebands the
    reconstruction uses, weighted by their calibration uncertainties.  The fit is pulled toward the warmer, brighter material, so
    the temperature it returns is too high and the surface density it implies is too low.  Right: the relative residuals of that
    fit, waveband by waveband, for three total surface densities.  The pattern is always the same, positive at the ends and
    negative in the middle, and it grows with the surface density; the fraction of the true value recovered is printed beside each
    curve.
    """
    cols = (PLASMA(30), PLASMA(120), PLASMA(210))
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.85))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.175, top=0.955, wspace=0.28)
    err = np.array([CAL_ERR.get(w, 0.20) for w in BANDS])
    tg = np.geomspace(0.98 * TFLOOR, 1.02 * TSURF, 2000)

    def one(sig):
        t = law(attenuating_column(sig, "mean"))
        obs = np.array([cnu(w) * sig * np.mean(planck(t, w)) for w in BANDS])
        model = np.array([[cnu(w) * planck(tt, w) for w in BANDS] for tt in tg])
        wg = 1.0 / (err * obs) ** 2
        amp = (wg * obs * model).sum(axis=1) / (wg * model * model).sum(axis=1)
        j = int(np.argmin((wg * (obs - amp[:, None] * model) ** 2).sum(axis=1)))
        fit = amp[j] * model[j]
        rec = obs[-1] / (cnu(500.0) * planck(tg[j], 500.0)) / sig
        return t, obs, fit, tg[j], rec

    sig0 = 3.0e22
    t, obs, fit, tfit, rec = one(sig0)
    wl = np.geomspace(120.0, 650.0, 300)
    smooth = np.array([cnu(w) * sig0 * np.mean(planck(t, w)) for w in wl])
    sfit = np.array([fit[-1] / (cnu(500.0) * planck(tfit, 500.0)) * cnu(w) * planck(tfit, w) for w in wl])
    ax[0].plot(wl, smooth, "-", lw=1.4, color=cols[1], label="the true, mixed emission")
    ax[0].plot(wl, sfit, "--", lw=1.2, color="0.3", label=r"best single fit, $T = %.1f$ K" % tfit)
    ax[0].plot(BANDS, obs, "o", ms=4.0, mfc="none", mec=cols[1], mew=1.1)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlim(120.0, 650.0)
    ax[0].set_xticks(BANDS)
    ax[0].set_xticklabels(["160", "250", "350", "500"])
    ax[0].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax[0].set_xlabel(r"wavelength ($\mu$m)")
    ax[0].set_ylabel(r"$I_\nu$ (MJy sr$^{-1}$)")
    ax[0].legend(frameon=False, fontsize=6.8, loc="lower left", handlelength=1.8)
    ax[0].text(0.97, 0.95, r"$\Sigma = 3\times10^{22}$ cm$^{-2}$" "\n" r"recovered $%.3f$" % rec, fontsize=6.8,
               color="0.35", ha="right", va="top", transform=ax[0].transAxes)

    for sig, col, lab in ((3.0e21, cols[0], r"$3\times10^{21}$"), (1.0e22, cols[1], r"$10^{22}$"),
                          (3.0e22, cols[2], r"$3\times10^{22}$")):
        _, o, f, _, r = one(sig)
        ax[1].plot(BANDS, 100.0 * (o - f) / o, "o-", ms=3.4, lw=1.2, color=col,
                   label=r"%s cm$^{-2}$, recovered %.3f" % (lab, r))
    ax[1].axhline(0.0, ls=":", lw=0.8, color="0.5")
    ax[1].set_xscale("log")
    ax[1].set_xlim(140.0, 570.0)
    ax[1].set_xticks(BANDS)
    ax[1].set_xticklabels(["160", "250", "350", "500"])
    ax[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax[1].set_ylim(-11.0, 27.0)
    ax[1].set_xlabel(r"wavelength ($\mu$m)")
    ax[1].set_ylabel("relative residual of the fit (%)")
    ax[1].legend(frameon=False, fontsize=6.6, loc="upper center", handlelength=1.7)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def _relation(eta, ns=36):
    """Fitted surface density and the restoring factor, for one value of the superseded geometry parameter.

    This and fig_orientation below belong to the earlier scheme, in which a single parameter eta set the deepest shielding of a
    line of sight.  That parameter no longer exists: the orientation average is now inside the angular average of
    attenuating_column.  They are left here only until it is decided what replaces fig_orientation, and they are not called by
    anything else.
    """
    u = (np.arange(256) + 0.5) / 256.0
    dep = 2.0 * np.minimum(u, 1.0 - u)
    tg = np.arange(5.0, 30.0, 0.02)
    mod = np.array([[cnu(w) * planck(t, w) for w in BANDS] for t in tg])
    sg = np.geomspace(3.0e20, 3.0e23, ns)
    rec, fs = np.empty(ns), np.empty(ns)
    for i, S in enumerate(sg):
        d = np.array([cnu(w) * S * np.mean(planck(law(eta * S * dep), w))
                      for w in BANDS])
        wg = 1.0 / (0.20 * d) ** 2
        a = (wg * d * mod).sum(1) / (wg * mod * mod).sum(1)
        j = int(np.argmin((wg * (d - a[:, None] * mod) ** 2).sum(1)))
        rec[i] = d[-1] / (cnu(500.0) * planck(tg[j], 500.0))
        fs[i] = S / rec[i]
    return rec, fs


def fig_orientation(out="fig_orientation.pdf"):
    """The effective geometry and the correction when orientations are random.

    Left: eta against the angle of the structure to the line of sight, for a
    filament, whose axis makes the angle theta, and for a sheet, whose normal
    makes the angle i.  The grey band is the isotropic probability density of
    that angle, proportional to its sine, drawn on an arbitrary scale so that
    the eye can weight the curves correctly: angles near the line of sight are
    rare.  Right: the correction factor averaged over isotropic orientations,
    against the face-on case and against the adopted geometry.
    """
    cyl, sht, ado = PLASMA(120), PLASMA(210), PLASMA(30)
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.8))
    th = np.linspace(0.0, 90.0, 400)
    r = np.radians(th)
    ax[0].fill_between(th, 0.0, 0.55 * np.sin(r), color="0.92", lw=0.0,
                       zorder=0)
    ax[0].text(72.0, 0.045, "isotropic probability\nof the angle",
               fontsize=6.4, color="0.45", ha="center")
    ax[0].plot(th, 0.5 * np.sin(r), "-", lw=1.4, color=cyl,
               label=r"filament, $\eta=\frac{1}{2}\sin\theta$")
    ax[0].plot(th, 0.5 * np.cos(r), "-", lw=1.4, color=sht,
               label=r"sheet, $\eta=\frac{1}{2}\cos i$")
    ax[0].axhline(0.393, ls="--", lw=0.8, color=cyl)
    ax[0].axhline(0.250, ls="--", lw=0.8, color=sht)
    ax[0].text(3.0, 0.405, r"$\langle\eta\rangle=0.39$", fontsize=6.8,
               color=cyl)
    ax[0].text(3.0, 0.262, r"$\langle\eta\rangle=0.25$", fontsize=6.8,
               color=sht)
    ax[0].set_xlim(0.0, 90.0)
    ax[0].set_ylim(0.0, 0.55)
    ax[0].set_xlabel(r"angle to the line of sight (degrees)")
    ax[0].set_ylabel(r"$\eta$")
    ax[0].legend(frameon=False, fontsize=7, loc="lower center")
    mu = (np.arange(24) + 0.5) / 24.0
    sth = np.sqrt(1.0 - mu ** 2)
    rec0, f0 = _relation(0.50)
    reca, fa = _relation(0.35)
    gc = np.mean([np.interp(np.log10(rec0), np.log10(_relation(max(0.5 * s, 0.02))[0]),
                            _relation(max(0.5 * s, 0.02))[1]) for s in sth], axis=0)
    gs_ = np.mean([np.interp(np.log10(rec0), np.log10(_relation(max(0.5 * m, 0.02))[0]),
                             _relation(max(0.5 * m, 0.02))[1]) for m in mu], axis=0)
    ax[1].plot(rec0, f0, ":", lw=1.0, color="0.45", label=r"face-on, $\eta=0.5$")
    ax[1].plot(rec0, gc, "-", lw=1.4, color=cyl, label="filaments, random")
    ax[1].plot(rec0, gs_, "-", lw=1.4, color=sht, label="sheets, random")
    ax[1].plot(reca, fa, "--", lw=1.2, color=ado, label=r"adopted, $\eta=0.35$")
    ax[1].set_xscale("log")
    ax[1].set_xlim(1.0e21, 1.0e23)
    ax[1].set_ylim(0.95, 4.0)
    ax[1].set_xlabel(r"$N_{\rm H_2}$ of the uncorrected image (cm$^{-2}$)")
    ax[1].set_ylabel(r"correction factor $f$")
    ax[1].legend(frameon=False, fontsize=7, loc="upper left")
    fig.tight_layout(pad=0.4)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_relation(out="fig_relation.pdf"):
    """The correction factor against the uncorrected surface density, for the shapes and for the shortest path.

    The sphere and the layer are the two extremes of a uniform structure and the cube, which has the symmetry of neither, falls
    between them; the adopted factor is the average of the first two.  The dotted curve is what is obtained by replacing the
    angular average of Eq. (6) by the column along the shortest path with its geometric factor of one quarter, which agrees below
    3e21 cm^-2 and diverges above 1e22.
    """
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.8))
    for geom, col, lw, lab in (("layer", PLASMA(20), 1.0, "plane-parallel layer"),
                               ("cube", PLASMA(90), 1.0, "cube"),
                               ("mean", PLASMA(150), 1.6, "adopted, the average"),
                               ("sphere", PLASMA(215), 1.0, "sphere")):
        fit, f = relation(geom)
        ax.plot(fit, f, "-", lw=lw, color=col, label=lab)
    fit, f = relation("short")
    ax.plot(fit, f, ":", lw=1.2, color="0.35", label="shortest path")
    ax.set_xscale("log")
    ax.set_xlim(5.0e20, 3.0e23)
    ax.set_ylim(0.95, 4.0)
    ax.set_xlabel(r"$N_{\rm H_2}$ of the uncorrected image (cm$^{-2}$)")
    ax.set_ylabel(r"correction factor $f$")
    ax.legend(frameon=False, fontsize=6.8, loc="upper left", handlelength=1.8, labelspacing=0.3)
    fig.tight_layout(pad=0.4)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_masses(stem="simsky4.w40", out="fig_masses.pdf"):
    """How accurately a mass is recovered, against the aperture over which it is measured.

    Left: the mass of a core, against the radius of the aperture in units of the Bonnor-Ebert radius, with the background taken
    between 2 and 2.5 radii.  Right: the mass per unit length of a filament, against the half-width of the window over which the
    profile is integrated, with the background level taken from the fit of A r^-alpha + B.  The core mass is recovered to a few
    per cent out to one radius and drifts upward beyond it as the ring of over-correction enters the aperture; the filament mass
    is over-corrected at every window, because the ring lies inside all of them.
    """
    from scipy.optimize import curve_fit

    tag = lambda b: ("%.1f" % b).replace(".", "p")
    T = g("%s.surfdens.r%s.cores.fits" % (stem, tag(FINE)))
    D = g("hi.surface.density.r%s.fits" % tag(FINE))
    cname = "hi.surface.density.r%s.correct.fits" % tag(FINE)
    if not os.path.exists(cname):
        cname = "hi.surface.density.r%s.corr.fits" % tag(FINE)
    C = g(cname)
    cat = [l.split() for l in open("%s.cores.cores.txt" % stem) if not l.startswith("#") and l.strip()]
    ny, nx = T.shape

    aps = np.arange(0.3, 2.01, 0.1)
    core = {k: [] for k in ("u", "c")}
    for ap in aps:
        ru, rc = [], []
        for q in cat:
            y, x, rb = int(q[1]), int(q[2]), float(q[4]) / PIX
            h = int(np.ceil(2.5 * rb)) + 2
            if min(y, x) < h or y > ny - h - 1 or x > nx - h - 1:
                continue
            sl = (slice(y - h, y + h + 1), slice(x - h, x + h + 1))
            gy, gx = np.mgrid[-h:h + 1, -h:h + 1]
            rr = np.hypot(gy, gx)
            m1, m2 = rr <= ap * rb, (rr > 2.0 * rb) & (rr <= 2.5 * rb)
            o = [float(np.sum(M[sl][m1] - np.median(M[sl][m2]))) for M in (T, D, C)]
            if o[0] > 0.0:
                ru.append(o[1] / o[0])
                rc.append(o[2] / o[0])
        core["u"].append(np.median(ru))
        core["c"].append(np.median(rc))

    prof = np.nanmedian(T, axis=0)
    pk = [i for i in range(110, len(prof) - 110)
          if prof[i] == np.nanmax(prof[i - 30:i + 31]) and prof[i] > 1.5 * np.nanmedian(prof)]
    dd = np.arange(-105, 106)
    wins = np.arange(30.0, 201.0, 10.0)

    def lev(M, x0):
        p = np.nanmedian(M[:, x0 - 105:x0 + 106], axis=0)
        r = np.array([k * PIX for k in range(1, 106) if 20.0 <= k * PIX <= 180.0])
        v = np.array([0.5 * (p[dd == -k][0] + p[dd == k][0]) for k in range(1, 106) if 20.0 <= k * PIX <= 180.0])
        po, _ = curve_fit(lambda z, A, a, B: A * z ** (-a) + B, r, v,
                          p0=[v[0] * r[0] ** 0.75, 0.75, v[-1]], maxfev=20000)
        return p - po[2]
    qs = {k: [lev(M, x) for x in pk] for k, M in (("t", T), ("u", D), ("c", C))}
    fil = {k: [] for k in ("u", "c")}
    for R in wins:
        m = np.abs(dd) * PIX <= R
        for k in ("u", "c"):
            fil[k].append(np.median([np.sum(qs[k][i][m]) / np.sum(qs["t"][i][m]) for i in range(len(pk))]))

    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.85))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.175, top=0.905, wspace=0.28)
    for a, xx, dat, xl, ti, lc in ((ax[0], aps, core, r"aperture radius in units of $R_{\rm BE}$", "cores", "lower right"),
                                   (ax[1], wins, fil, "half-width of the window (arcsec)", "filaments", "center left")):
        a.plot(xx, dat["u"], "-", lw=1.4, color=PLASMA(30), label="uncorrected")
        a.plot(xx, dat["c"], "-", lw=1.4, color=PLASMA(150), label="corrected")
        a.axhline(1.0, ls="--", lw=0.9, color="0.4")
        a.set_xlim(xx[0], xx[-1])
        a.set_ylim(0.35, 1.65)
        a.set_xlabel(xl)
        a.set_ylabel("fraction of the true mass")
        a.set_title(ti, fontsize=7.4, pad=3)
        a.legend(frameon=False, fontsize=6.8, loc="center right", handlelength=1.7)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_base(stem="simsky4.w40", out="fig_base.pdf"):
    """What the choice of the base image costs, measured on the simulated sky.

    The factor is evaluated on the reconstruction at one resolution and the difference is added to the image at the finest
    resolution, which is Eq. (12) with the base image taken in turn at each available resolution.  Left: the mass of a core
    recovered within three apertures.  Right: the shape of a filament, as ratios to the true image so that one axis serves for
    all four quantities.  The coarsest base gives the best masses and the poorest shapes, and the gain in shape is exhausted by
    18.2 arcsec.
    """
    import superposition_test as sup

    tag = lambda b: ("%.1f" % b).replace(".", "p")
    beams = [36.3, 24.9, 18.2, 13.5]
    T = g("%s.surfdens.r%s.cores.fits" % (stem, tag(FINE)))
    D = g("hi.surface.density.r%s.fits" % tag(FINE))
    fitted, factor = sup.build_relation(BANDS) if hasattr(sup, "build_relation") else (None, None)
    if fitted is None:
        from hicorrect2 import build_relation
        fitted, factor = build_relation(BANDS)
    fac = lambda x: np.interp(np.log10(np.maximum(x, 1.0)), np.log10(fitted), factor,
                              left=factor[0], right=factor[-1])
    images = {}
    for b in beams:
        B = g("hi.surface.density.r%s.fits" % tag(b))
        images[b] = D + (fac(B) - 1.0) * B

    cat = [l.split() for l in open("%s.cores.cores.txt" % stem) if not l.startswith("#") and l.strip()]
    ny, nx = T.shape

    def coremass(M, ap):
        out = []
        for q in cat:
            y, x, rb = int(q[1]), int(q[2]), float(q[4]) / PIX
            h = int(np.ceil(2.5 * rb)) + 2
            if min(y, x) < h or y > ny - h - 1 or x > nx - h - 1:
                continue
            sl = (slice(y - h, y + h + 1), slice(x - h, x + h + 1))
            gy, gx = np.mgrid[-h:h + 1, -h:h + 1]
            rr = np.hypot(gy, gx)
            ap_m, an_m = rr <= ap * rb, (rr > 2.0 * rb) & (rr <= 2.5 * rb)
            a = float(np.sum(T[sl][ap_m] - np.median(T[sl][an_m])))
            v = float(np.sum(M[sl][ap_m] - np.median(M[sl][an_m])))
            if a > 0.0:
                out.append(v / a)
        return np.median(out)

    prof = np.nanmedian(T, axis=0)
    pk = [i for i in range(110, len(prof) - 110)
          if prof[i] == np.nanmax(prof[i - 30:i + 31]) and prof[i] > 1.5 * np.nanmedian(prof)]
    dd = np.arange(-105, 106)

    def filament(M, x0):
        """Index, half-maximum width, mass per unit length and crest amplitude of one filament, with the level fitted."""
        from scipy.optimize import curve_fit
        p = np.nanmedian(M[:, x0 - 105:x0 + 106], axis=0)
        r = np.array([k * PIX for k in range(1, 106) if 20.0 <= k * PIX <= 180.0])
        v = np.array([0.5 * (p[dd == -k][0] + p[dd == k][0]) for k in range(1, 106) if 20.0 <= k * PIX <= 180.0])
        po, _ = curve_fit(lambda z, A, a, B: A * z ** (-a) + B, r, v,
                          p0=[v[0] * r[0] ** 0.75, 0.75, v[-1]], maxfev=20000)
        q = p - po[2]
        a0 = q[dd == 0][0]
        half = []
        for sgn in (-1, 1):
            k = np.arange(1, 106)
            qq = np.array([q[dd == sgn * j][0] for j in k])
            i = int(np.argmax(qq < 0.5 * a0))
            half.append(k[i - 1] + (0.5 * a0 - qq[i - 1]) * (k[i] - k[i - 1]) / (qq[i] - qq[i - 1]))
        m = np.abs(dd) * PIX <= 100.0
        return po[1], (half[0] + half[1]) * PIX, np.sum(q[m]) * PIX, a0

    ref = np.median(np.array([filament(T, x) for x in pk]), axis=0)
    fil = {b: np.median(np.array([filament(images[b], x) for x in pk]), axis=0) for b in beams}
    filu = np.median(np.array([filament(D, x) for x in pk]), axis=0)

    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.95))
    fig.subplots_adjust(left=0.098, right=0.985, bottom=0.175, top=0.905, wspace=0.30)
    cols = (PLASMA(30), PLASMA(120), PLASMA(210))
    for k, ap in enumerate((0.5, 1.0, 1.5)):
        ax[0].plot(beams, [coremass(images[b], ap) for b in beams], "o-", ms=3.6, lw=1.4, color=cols[k],
                   label=r"aperture $%.1f\,R_{\rm BE}$" % ap)
        ax[0].axhline(coremass(D, ap), ls=":", lw=0.9, color=cols[k])
    ax[0].axhline(1.0, ls="--", lw=0.9, color="0.4")
    ax[0].set_xlim(39.0, 11.0)
    ax[0].set_ylim(0.40, 1.20)
    ax[0].set_xlabel("resolution of the base image (arcsec)")
    ax[0].set_ylabel("fraction of the true core mass")
    ax[0].legend(frameon=False, fontsize=6.6, loc="center right", handlelength=1.7)
    ax[0].text(38.0, 0.60, "dotted: uncorrected", fontsize=6.4, color="0.40")
    ax[0].set_title("cores", fontsize=7.4, pad=3)

    names = [("index " + r"$\alpha$", 0), ("half-maximum width", 1), ("mass per unit length", 2), ("crest amplitude", 3)]
    styl = ["o-", "s-", "^-", "v-"]
    for (lab, j), st, col in zip(names, styl, (PLASMA(20), PLASMA(90), PLASMA(160), PLASMA(225))):
        ax[1].plot(beams, [fil[b][j] / ref[j] for b in beams], st, ms=3.4, lw=1.3, color=col, label=lab)
        ax[1].axhline(filu[j] / ref[j], ls=":", lw=0.9, color=col)
    ax[1].axhline(1.0, ls="--", lw=0.9, color="0.4")
    ax[1].set_xlim(39.0, 11.0)
    ax[1].set_xlabel("resolution of the base image (arcsec)")
    ax[1].set_ylabel("ratio to the true image")
    ax[1].legend(frameon=False, fontsize=6.4, loc="upper left", handlelength=1.7, labelspacing=0.25)
    ax[1].text(0.97, 0.03, "dotted: uncorrected", fontsize=6.4, color="0.40", ha="right", va="bottom",
               transform=ax[1].transAxes)
    ax[1].set_title("filaments", fontsize=7.4, pad=3)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_residual(stem="simsky4.w40", out=None, region=slice(150, 1650)):
    """Where the reconstruction and the correction depart from the truth, as maps and as profiles.

    Upper panels: the relative residual of the reconstruction at the finest resolution before and after the correction.  The
    uncorrected image is too low almost everywhere, most of all in the dense material.  The corrected one is right to a few per
    cent over the diffuse field, but leaves a ring of over-correction around every compact structure and along the flanks of every
    filament.  Lower panels: the same thing as profiles, averaged over the filaments and over the injected cores, which shows that
    the ring is not an artifact of a single structure.
    """
    if out is None:
        out = "fig_%s_residual.pdf" % stem.split(".")[0]
    fine = ("%.1f" % FINE).replace(".", "p")
    T = g("%s.surfdens.r%s.cores.fits" % (stem, fine))
    D = g("hi.surface.density.r%s.fits" % fine)
    cname = "hi.surface.density.r%s.correct.fits" % fine
    if not os.path.exists(cname):
        cname = "hi.surface.density.r%s.corr.fits" % fine
    C = g(cname)
    ok = (D > 0.0) & np.isfinite(T) & (T > 0.0)
    rd = np.where(ok, D / np.maximum(T, 1.0) - 1.0, np.nan)
    rc = np.where(ok, C / np.maximum(T, 1.0) - 1.0, np.nan)

    fig = plt.figure(figsize=(7.1, 5.9))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.45, 1.0], hspace=0.1, wspace=0.43,
                          left=0.0, right=0.905, bottom=0.085, top=0.955)
    for k, (x, lab) in enumerate(((rd, "$hires$ 13.5″ surface density"), 
                                  (rc, "$hires$ 13.5″ surface density, corrected"))):
        a = fig.add_subplot(gs[0, 3 * k:3 * k + 3])
        n = x[region, region].shape[0]
        ext = [-n * PIX / 120.0, n * PIX / 120.0] * 2
        im = a.imshow(x[region, region], origin="lower", extent=ext, cmap=SD_CMAP, vmin=-0.45, vmax=0.45)
        a.set_title(lab, fontsize=7, pad=3)
        a.set_xticks([])
        a.set_yticks([])
        cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.015)
        cb.ax.tick_params(labelsize=7.3)

##    cax = fig.add_axes([0.915, 0.475, 0.016, 0.480])
##    cb = fig.colorbar(im, cax=cax)
    cb.set_label("(reconstructed $-$ true) / true")

    prof = np.nanmedian(np.where(ok, T, np.nan), axis=0)
    off = np.arange(-70, 71)
    pk = [i for i in range(80, len(prof) - 80)
          if prof[i] == np.nanmax(prof[i - 30:i + 31]) and prof[i] > 1.5 * np.nanmedian(prof)]

    # ---- the profiles themselves, stacked over the filaments
    p0 = fig.add_subplot(gs[1, 0:2])
    for M, col, lab, lw in ((T, "0.25", "true", 1.0), (D, PLASMA(30), "uncorrected", 1.0),
                            (C, PLASMA(150), "corrected", 1.0)):
        st = np.array([[np.nanmedian(M[:, x + d]) for d in off] for x in pk])
        q = np.nanmedian(st, axis=0)
        q = q - np.median(np.concatenate([q[:12], q[-12:]]))
        p0.plot(off * PIX, q, "-", lw=lw, color=col, label=lab)
    p0.set_yscale("log")
    p0.set_xlim(off[0] * PIX - 20, off[-1] * PIX + 20)
    p0.set_ylim(3.0e20, 5.0e22)
    p0.set_xlabel("Offset from crest (″)")
    p0.set_ylabel(r"Background-subtracted $N_{\rm H_2}$ (cm$^{-2}$)", labelpad=-1)
    p0.legend(frameon=False, fontsize=6.6, loc="upper left", handlelength=1.6, labelspacing=0.25)
    p0.text(0.5 * FINE + 8.0, 0.04, "%.1f″ beam" % FINE, fontsize=7, color="0.40", va="bottom",
            transform=p0.get_xaxis_transform())
    for v in (-0.5 * FINE, 0.5 * FINE):
        p0.axvline(v, ls="--", lw=0.5, color="0.3")

    # ---- across the filaments, stacked
    b = fig.add_subplot(gs[1, 2:4])
    for r, col, lab in ((rd, PLASMA(30), "uncorrected"), (rc, PLASMA(150), "corrected")):
        stack = np.array([[np.nanmedian(r[:, x + d]) for d in off] for x in pk])
        b.plot(off * PIX, np.nanmedian(stack, axis=0), "-", lw=1.0, color=col, label=lab)
    b.axhline(0.0, ls=":", lw=0.8, color="0.5")
    for v in (-0.5 * FINE, 0.5 * FINE):
        b.axvline(v, ls="--", lw=0.5, color="0.3")
    b.set_xlim(off[0] * PIX - 20, off[-1] * PIX + 20)
    b.set_xlabel("Offset from crest (″)")
##    b.set_ylabel("(reconstructed $-$ true) / true")
    b.legend(frameon=False, fontsize=6.6, loc="lower right", handlelength=1.6, labelspacing=0.25)
    b.set_title("Relative residuals (%d filaments)" % len(pk), fontsize=7.3, pad=3)
    p0.set_title("Profiles of %d filaments" % len(pk), fontsize=7.3, pad=3)

    # ---- around the injected cores, stacked
    c = fig.add_subplot(gs[1, 4:6])
    cat = [l.split() for l in open("%s.cores.cores.txt" % stem) if not l.startswith("#") and l.strip()]
    edges = np.linspace(0.0, 3.0, 16)
    mid = 0.5 * (edges[:-1] + edges[1:])
    ny, nx = T.shape
    for r, col, lab in ((rd, PLASMA(30), "uncorrected"), (rc, PLASMA(150), "corrected")):
        rows = []
        for q in cat:
            y, x, rb = int(q[1]), int(q[2]), float(q[4]) / PIX
            h = int(np.ceil(3.0 * rb)) + 2
            if y - h < 0 or y + h + 1 > ny or x - h < 0 or x + h + 1 > nx:
                continue
            sl = (slice(y - h, y + h + 1), slice(x - h, x + h + 1))
            gy, gx = np.mgrid[-h:h + 1, -h:h + 1]
            rr = np.hypot(gy, gx) / rb
            v = r[sl]
            rows.append([np.nanmedian(v[(rr >= a0) & (rr < a1)]) for a0, a1 in zip(edges[:-1], edges[1:])])
        c.plot(mid, np.nanmedian(np.array(rows), axis=0), "-", lw=1.0, color=col, label=lab)
    c.axhline(0.0, ls=":", lw=0.5, color="0.3")
    rbe = np.median([float(q[4]) for q in cat])
    c.axvline(0.5 * FINE / rbe, ls="--", lw=0.5, color="0.3")
    c.text(0.5 * FINE / rbe + 0.06, 0.93, "%.1f″ beam" % FINE, fontsize=7, color="0.3", va="top",
           transform=c.get_xaxis_transform())
    c.set_xlim(0.0, 3.0)
    c.set_ylim(-0.45, 0.35)
    c.set_xlabel(r"Offset from center ($r\,R^{-1}_{\rm BE}$)")
##    c.set_yticklabels([])   
##    c.set_ylabel("(reconstructed $-$ true) / true")
    c.legend(frameon=False, fontsize=6.6, loc="lower right", handlelength=1.6, labelspacing=0.25)
    c.set_title("Relative residuals (%d cores)" % len(rows), fontsize=7.3, pad=3)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_simsky_bands(stem="simsky5.w40", out=None, region=slice(400, 1400)):
    """The simulated sky as it is observed: the five wavebands used in the figure and the true surface density.

    The panels are arranged by resolution, from the sharpest at the upper left to the coarsest at the lower right, with the true
    surface density in the place its resolution would occupy.  Each panel carries its own logarithmic scale, because the sky is
    two orders of magnitude fainter at 70 than at 250 micron and one common scale would show nothing at the ends.  The 100 micron
    image exists and is passed to hires, but it is left out here to leave room for the true surface density.
    """
    if out is None:
        out = "fig_%s_bands.pdf" % stem.split(".")[0]
    panels = [("070", "r8p4", 8.4), ("160", "r13p5", 13.5), (None, None, FINE),
              ("250", "r18p2", 18.2), ("350", "r24p9", 24.9), ("500", "r36p3", 36.3)]
    fig = plt.figure(figsize=(7.1, 4.9))
    gs = fig.add_gridspec(2, 3, wspace=0.16, hspace=0.16, left=0.006, right=0.955, bottom=0.006, top=0.945)
    for k, (band, rtag, beam) in enumerate(panels):
        a = fig.add_subplot(gs[k // 3, k % 3])
        if band is None:
            x = g("%s.surfdens.r%s.cores.fits" % (stem, ("%.1f" % FINE).replace(".", "p")))[region, region]
            cmap, lab, unit = SD_CMAP, "true surface density", r"$N_{\rm H_2}$ (cm$^{-2}$)"
        else:
            x = g("%s.%s.%s.cores.fits" % (stem, band, rtag))[region, region]
            cmap, unit = SD_CMAP, r"$I_\nu$ (MJy sr$^{-1}$)"
            lab = r"%s{\um}, $%.1f^{\prime\prime}$" % (band, beam)
            lab = "%s $\\mu$m, %.1f arcsec" % (band, beam)
        v = x[np.isfinite(x) & (x > 0.0)]
        vmin, vmax = np.percentile(v, 1.0), np.percentile(v, 99.9)
        n = x.shape[0]
        ext = [-n * PIX / 120.0, n * PIX / 120.0] * 2
        im = a.imshow(x, origin="lower", extent=ext, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
        a.set_title(lab, fontsize=7.4, pad=3)
        a.set_xticks([])
        a.set_yticks([])
        cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.015)
        cb.set_label(unit, fontsize=6.4, labelpad=1)
        cb.ax.tick_params(labelsize=6.0, pad=1)
        cb.ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_aquila(out="fig_aquila.pdf"):
    """The Aquila field before and after the correction, with the factor applied and the dust temperature that follows.

    The two surface density panels share one logarithmic scale, so the change can be read directly.  The factor is the one
    interpolated onto the base image at the coarsest resolution and added to both resolutions; the temperature is the one implied
    by the corrected surface density at the coarsest resolution.
    """
    dd = g("hi.surface.density.r%s.fits" % ("%.1f" % FINE).replace(".", "p"))
    cc = g("hi.surface.density.r%s.correct.fits" % ("%.1f" % FINE).replace(".", "p"))
    ff = g("hi.surface.density.r%s.cfactor.fits" % ("%.1f" % BEAM).replace(".", "p"))
    tt = g("hi.surface.density.r%s.temperature.corr.fits" % ("%.1f" % BEAM).replace(".", "p"))
    lo, hi = np.nanpercentile(dd[dd > 0.0], 0.5), np.nanpercentile(cc[cc > 0.0], 99.9)
    n = dd.shape[0]
    ext = [-n * PIX / 120.0, n * PIX / 120.0] * 2

    fig = plt.figure(figsize=(7.1, 6.9))
    gs = fig.add_gridspec(2, 2, wspace=0.16, hspace=0.16, left=0.010, right=0.940, bottom=0.010, top=0.955)
    panels = [(dd, "uncorrected", SD_CMAP, LogNorm(vmin=lo, vmax=hi), r"$N_{\rm H_2}$ (cm$^{-2}$)"),
              (cc, "corrected", SD_CMAP, LogNorm(vmin=lo, vmax=hi), r"$N_{\rm H_2}$ (cm$^{-2}$)"),
              (ff, "correction factor", SD_CMAP, None, r"$f$"),
              (tt, "corrected dust temperature", SD_CMAP, None, r"$T$ (K)")]
    for k, (x, lab, cmap, nor, unit) in enumerate(panels):
        a = fig.add_subplot(gs[k // 2, k % 2])
        v = np.where(np.isfinite(x) & (x > 0.0), x, np.nan)          # the blank border must not set the scale
        kw = dict(norm=nor) if nor is not None else dict(vmin=np.nanpercentile(v, 0.5),
                                                         vmax=np.nanpercentile(v, 99.5))
        im = a.imshow(v, origin="lower", extent=ext, cmap=cmap, **kw)
        a.set_title(lab, fontsize=7.6, pad=3)
        a.set_xticks([])
        a.set_yticks([])
        cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.015)
        cb.set_label(unit, fontsize=6.6, labelpad=1)
        cb.ax.tick_params(labelsize=6.2, pad=1)
        cb.ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_simsky(stem="simsky4.w40", out="fig_simsky.pdf"):
    """The simulated sky before and after correction, and the truth."""
    T = g("%s.surfdens.r13p5.cores.fits" % stem)
    D = g("hi.surface.density.r13p5.fits")
    C = g("hi.surface.density.r13p5.correct.fits")
    s = slice(500, 1300)
    d = [X[s, s] for X in (D, C, T)]
    n = d[0].shape[0]
    ext = [-n * PIX / 120.0, n * PIX / 120.0] * 2
    vmin, vmax = 6.0e21, 8.0e22
    fig = plt.figure(figsize=(7.1, 2.55))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.045], wspace=0.02, left=0.005, right=0.945, bottom=0.01, top=0.90)
    ax = [fig.add_subplot(gs[0, i]) for i in range(3)]
    cax = fig.add_subplot(gs[0, 3])
    for a, x, lab in zip(ax, d, ("hires", "corrected", "true")):
        im = a.imshow(x, origin="lower", extent=ext, cmap=SD_CMAP, norm=LogNorm(vmin=vmin, vmax=vmax))
        a.set_title(lab, fontsize=8, pad=3)
        a.set_xticks([])
        a.set_yticks([])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"$N_{\rm H_2}$ (cm$^{-2}$)", fontsize=8, labelpad=-14)
    cb.ax.tick_params(labelsize=7, pad=1)
    cb.ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.canvas.draw()
    p_ax, p_cb = ax[2].get_position(), cax.get_position()
    cax.set_position([p_cb.x0, p_ax.y0, p_cb.width, p_ax.height])
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def core_table(stem="simsky4.w40"):
    """Recovered fraction of the integrated excess of every injected core."""
    T = g("%s.surfdens.r13p5.cores.fits" % stem)
    D = g("hi.surface.density.r13p5.fits")
    C = g("hi.surface.density.r13p5.correct.fits")
    cat = [l.split() for l in open("%s.cores.cores.txt" % stem)
           if not l.startswith("#") and l.strip()]
    ny, nx = D.shape
    rows = []
    for c in cat:
        y, x, mbe, rb = int(c[1]), int(c[2]), float(c[3]), float(c[4]) / PIX
        h = int(np.ceil(2.0 * rb)) + 2
        if y - h < 0 or y + h + 1 > ny or x - h < 0 or x + h + 1 > nx:
            continue
        sl = (slice(y - h, y + h + 1), slice(x - h, x + h + 1))
        gy, gx = np.mgrid[-h:h + 1, -h:h + 1]
        rr = np.hypot(gy, gx)
        ap, an = rr <= rb, (rr > 1.5 * rb) & (rr <= 2.0 * rb)
        w = [X[sl] for X in (T, D, C)]
        if any(not np.all(np.isfinite(v[ap | an])) for v in w):
            continue
        o = [float(np.sum(v[ap] - np.median(v[an]))) for v in w]
        if o[0] <= 0:
            continue
        rows.append((mbe, o[1] / o[0], o[2] / o[0]))
    return np.array(rows)


def fig_core_recovery(out="fig_core_recovery.pdf"):
    r = core_table()
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.7))
    ax.axhline(1.0, ls=":", color="0.5")
    ax.plot(r[:, 0], r[:, 1], "o", ms=3.0, mfc="none", mec=PLASMA(40),
            mew=0.9, label="hires")
    ax.plot(r[:, 0], r[:, 2], "s", ms=3.0, mfc="none", mec=PLASMA(180),
            mew=0.9, label="corrected")
    ax.set_xscale("log")
    ax.set_xlabel(r"true mass of the core ($M_\odot$)")
    ax.set_ylabel("recovered fraction of the excess")
    ax.set_ylim(0.0, 1.6)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout(pad=0.4)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)
    print("   hires     median %.3f" % np.median(r[:, 1]))
    print("   corrected median %.3f" % np.median(r[:, 2]))


if __name__ == "__main__":
    if os.path.exists("attenuation_law.txt"):
        fig_shielding()
    fig_method()
    fig_geometry()
    fig_simsky()
    fig_deficit()
    core_table()
    fig_base()
    fig_masses()
    fig_core_recovery()
    fig_residual()
#    fig_aquila()
#    fig_angular()
#    fig_superposition()
#    fig_mixing()
#    fig_relation()
