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

# Helvetica on the author's machine, Liberation Sans (metrically identical) elsewhere. The italic face is looked up by the name
# and style that the font file itself declares, so that symbols in math mode, such as the N of N_H2, come out italic even when
# the oblique face is registered under a family name of its own.
FONTS = (("/Users/amenshch/Library/Fonts/Helvetica.ttf", "/Users/amenshch/Library/Fonts/Helvetica-Oblique.ttf"),
         ("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
          "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf"))
for _reg, _ita in FONTS:
    if os.path.exists(_reg) and os.path.exists(_ita):
        font_manager.fontManager.addfont(_reg)
        font_manager.fontManager.addfont(_ita)
        _er = [e for e in font_manager.fontManager.ttflist if e.fname == _reg][0]
        _ei = [e for e in font_manager.fontManager.ttflist if e.fname == _ita][0]
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = [_er.name]
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.rm'] = _er.name
        plt.rcParams['mathtext.sf'] = _er.name
        plt.rcParams['mathtext.it'] = "%s:%s" % (_ei.name, "italic" if _ei.style == "italic" else "oblique")
        plt.rcParams['mathtext.bf'] = "%s:bold" % _er.name
        break

sys.path.insert(0, "/mnt/project")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))     # the codes next to this script come first
try:
    import aa_plotstyle
    from dkbluered import my_cmap as SD_CMAP
    aa_plotstyle.apply()
except ImportError:                                   # outside the project
    SD_CMAP = plt.get_cmap("magma")
    matplotlib.rcParams.update({"font.size": 8, "axes.linewidth": 0.5, "xtick.direction": "in", "ytick.direction": "in"})

plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.labeltop'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['ytick.labelright'] = True
plt.rcParams['lines.linewidth'] = 0.9   # default is 1.5   
plt.rcParams['axes.linewidth'] = 0.5      # default is 0.8   
plt.rcParams['xtick.major.width'] = 0.4   # default 0.8
plt.rcParams['ytick.major.width'] = 0.4
plt.rcParams['xtick.minor.width'] = 0.4
plt.rcParams['ytick.minor.width'] = 0.4

PLASMA = plt.get_cmap("plasma", 256)

# The house style of the figures, as set by Figs. 3, 5, 6 and 9: curves at the default width, reference lines thin and gray,
# titles and legends small, legends without frames, and labels that start with a capital letter.
FS_TITLE, FS_LEG, FS_NOTE, FS_CB = 7.3, 7.0, 6.6, 7.0
LEG = dict(frameon=False, fontsize=FS_LEG, handlelength=1.6, labelspacing=0.5)
REF = dict(ls=":", lw=0.75, color="black")
COL3 = (PLASMA(30), PLASMA(120), PLASMA(210))
CCURVE = ('red', 'green', 'blue', 'magenta')


def clean_log(ax, which="both"):
    """Remove the labels of the minor ticks of logarithmic axes, keeping the ticks themselves."""
    if which in ("x", "both"):
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if which in ("y", "both"):
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
PIX = 3.0
BEAM, FINE = 36.3, 13.5
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
TFLOOR, TSURF, SSTAR, GAMMA = 7.045362, 22.830197, 2.731651e+21, 1.682779   # fitted to 2086 spheres and 7 full-grid cylinders
##TFLOOR, TSURF, SSTAR, GAMMA = 6.812424, 22.626125, 2.538614e+21, 1.606502   # fitted to 2086 spheres and 7 cylinders
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
    along its normal, for the average of the two, or for a uniform cube sampled along the line joining two opposite faces. The
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



def colorbar_beside(fig, ax, im, size="4.5%", pad=0.04, **kw):

    """A color bar exactly as tall as the image axes beside it.

    fig.colorbar(im, ax=...) takes its space from the axes and draws a bar somewhat shorter than an image with a fixed aspect
    ratio. A divider of mpl_toolkits.axes_grid1 attaches a new axes to the right of the image axes instead, and that axes always
    has the height of the image."""

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    cax = make_axes_locatable(ax).append_axes("right", size=size, pad=pad)
    cb = fig.colorbar(im, cax=cax, **kw)
    cax.tick_params(left=False, labelleft=False, right=True, labelright=True, top=False, bottom=False, labelbottom=False)

    return cb


def fig_superposition(out="fig_superposition.pdf"):
  
    """What the method assumes about a pixel that receives emission from more than one structure.

    Every pixel collects the emission of everything along its line of sight, and one modified blackbody is fitted to that total,
    so what biases the fit is the whole distribution of temperatures along the line of sight. The method replaces that
    distribution by the one a single uniform body of the same total surface density would have.

    Upper left: the dust temperature along a line of sight of 2e22 cm^-2, as the method assumes it and as it would be if half the
    surface density were a diffuse cirrus and half a sphere embedded in it. Upper right: the emission that follows, with the
    surface density that a single modified blackbody fitted to each returns. Lower left: the correction factor required, the curve
    the method applies and the range spanned by embedding a sphere in a cirrus of 1e21, 3e21 and 1e22 cm^-2 at every position
    within it. Lower right: the factor against the position of the sphere within the cirrus. All of it is computed by
    superposition_test.py, which is described in the appendix.
    """

    import superposition_test as sup

    one, sep = CCURVE[2], CCURVE[0]
    fig, axx = plt.subplots(2, 2, figsize=(7.0, 5.6))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.095, top=0.945, wspace=0.28, hspace=0.23)
    ax = [axx[0, 0], axx[0, 1], axx[1, 0]]
    axq = axx[1, 1]

    sig, sig_c = 2.0e22, 1.0e22
    sig_s, q = sig - sig_c, 0.0        # the sphere against the near side of the cirrus, the arrangement that departs most

    # ---- the temperature along the line of sight, against the surface density measured from the observer's side
    u = (np.arange(400) + 0.5) / 400.0
    sc = 1.0e22
    ax[0].plot(u * sig / sc, law(attenuating_column(sig, "mean", 400)), "-", color=one, label="one uniform body")
    tc = sup.layer_temperatures(sig_c)
    ts = sup.sphere_temperatures(sig_s, sig_c, q)
    nc = len(tc) // 2
    xc1 = (np.arange(nc) + 0.5) / nc * q * sig_c
    xs = q * sig_c + (np.arange(len(ts)) + 0.5) / len(ts) * sig_s
    xc2 = q * sig_c + sig_s + (np.arange(nc) + 0.5) / nc * (1.0 - q) * sig_c
    ax[0].plot(xc1 / sc, tc[:nc], "-", color=sep, label="cloud with an embedded sphere")
    ax[0].plot(xs / sc, ts[::-1], "-", color=sep)
    ax[0].plot(xc2 / sc, tc[nc:], "-", color=sep)
    ax[0].axvspan(q * sig_c / sc, (q * sig_c + sig_s) / sc, color="0.90", zorder=0)
    ax[0].text((q * sig_c + 0.5 * sig_s) / sc, 14.0, "Sphere", fontsize=7.75, color="black", ha="center")
    ax[0].set_xlim(0.0, sig / sc + 0.001)
    ax[0].set_ylim(7.5, 19.0)
    ax[0].tick_params(labelsize=7.7)
    ax[0].set_xlabel(r"Surface density from the near side ($10^{22}$ cm$^{-2}$)", fontsize=7.7, labelpad=1.5)
    ax[0].set_ylabel(r"Temperature $T$ (K)", fontsize=7.7, labelpad=1.5)
    ax[0].legend(loc="upper right", **LEG)
    ax[0].annotate("", xy=(0.02, 0.9), xytext=(0.15, 0.9), xycoords=ax[0].transAxes, textcoords=ax[0].transAxes,
                   annotation_clip=False, arrowprops=dict(arrowstyle="-|>", color="black", lw=0.5))
    ax[0].text(0.16, 0.9, "to observer", fontsize=7.5, color="black", va="center", transform=ax[0].transAxes)

    # ---- the emission each produces, and what one blackbody fitted to it returns
    wl = np.geomspace(130.0, 620.0, 200)
    for label, temps, wgts, col in (("one uniform body", law(attenuating_column(sig, "mean", 400)), np.full(400, sig / 400.0), one),
                                    ("cloud and sphere", np.concatenate([tc, ts]), np.concatenate([np.full(len(tc), sig_c / len(tc)),
                                                                                                    np.full(len(ts), sig_s / len(ts))]), sep)):
        obs = np.array([cnu(w) * np.sum(wgts * planck(temps, w)) for w in BANDS])
        sf = sup.fitted_column(obs)
        ax[1].plot(wl, [cnu(w) * np.sum(wgts * planck(temps, w)) for w in wl], "-", color=col,
                   label=r"%s, fitted $\Sigma_{\rm f} = %.2f{\times}10^{22}$ cm$^{-2}$" % (label, sf / 1.0e22))
        ax[1].plot(BANDS, obs, "o", ms=3.0, mfc="none", mec=col, mew=0.7)

    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlim(130.0, 620.0)
    ax[1].set_ylim(35.0, 200.0)
    ax[1].set_xticks(BANDS)
    ax[1].set_xticklabels(["160", "250", "350", "500"])
    ax[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax[1].tick_params(labelsize=7.7)
    ax[1].set_xlabel(r"Wavelength ($\mu$m)", fontsize=7.7, labelpad=1.5)
    ax[1].set_ylabel(r"Intensity $I_\nu$ (MJy sr$^{-1}$)", fontsize=7.7, labelpad=1.5)
    clean_log(ax[1], "y")
    ax[1].legend(loc="lower center", **LEG)
    ax[1].text(280.0, 60.0, r"Total $\Sigma = 2{\times}10^{22}$ cm$^{-2}$", fontsize=7.75, color="black", ha="center")

    # ---- the factor the method applies, and the factor each arrangement requires, both against the observed (fitted) surface
    # density: the method sees only the fitted surface density of a pixel, so the two must be compared at the same fitted value
    sg = np.geomspace(6.0e21, 1.2e23, 22)
    base = np.array([sup.one_body(v) / 1.0 for v in sg])
    xs_f = np.array([v / sup.one_body(v) for v in sg])
    fmeth = lambda sf: np.interp(np.log10(sf), np.log10(xs_f), base)
    lo, hi = np.full(len(sg), np.inf), np.full(len(sg), -np.inf)

    for scir in (1.0e21, 3.0e21, 1.0e22):
        for qq in (0.0, 0.25, 0.5):
            v = np.array([sup.two_component(x, scir, qq) if x > 1.6 * scir else np.nan for x in sg])
            good = np.isfinite(v)
            if good.sum() < 2:
                continue
            xf = sg[good] / v[good]
            vi = np.interp(np.log10(xs_f), np.log10(xf), v[good], left=np.nan, right=np.nan)
            lo, hi = np.fmin(lo, vi), np.fmax(hi, vi)

    ok = np.isfinite(lo) & np.isfinite(hi)
    ax[2].fill_between(xs_f[ok], lo[ok], hi[ok], color=sep, alpha=0.22, lw=0.0, label="range over the arrangement")
    ax[2].plot(xs_f, base, "-", color=one, label="factor the method applies")
    ax[2].set_xscale("log")
    ax[2].set_xlim(5.0e21, 6.0e22)
    clean_log(ax[2], "x")
    ax[2].set_ylim(1.0, 2.2)
    ax[2].tick_params(labelsize=7.7)
    ax[2].set_xlabel(r"Uncorrected $\Sigma$ (cm$^{-2}$)", fontsize=7.7, labelpad=0.5)
    ax[2].set_ylabel(r"Correction factor $f$", fontsize=7.7, labelpad=1.5)
    ax[2].legend(loc="lower right", **LEG)

    # ---- the error of the corrected surface density against where the structure sits in the cirrus
    qq = np.linspace(0.0, 1.0, 21)
    for scir, col in ((1.0e21, CCURVE[0]), (3.0e21, CCURVE[1]), (1.0e22, CCURVE[2])):
        for tot, ls in ((2.0e22, "-"), (5.0e22, "--")):
            req = np.array([sup.two_component(tot, scir, v) for v in qq])
            axq.plot(qq, fmeth(tot / req) / req - 1.0, ls, color=col)
        axq.plot([], [], "-", color=col, label=r"$\Sigma_{\rm c} = %s$" % {1.0e21: r"1{\times}10^{21}", 3.0e21: r"3{\times}10^{21}",
                                                                           1.0e22: r"1{\times}10^{22}"}[scir])
    axq.axhline(0.0, **REF)
    axq.plot([], [], "-", color="1", label=r"solid: $\Sigma = 2{\times}10^{22}$")
    axq.plot([], [], "--", color="1", label=r"dash: $\Sigma = 5{\times}10^{22}$")
    axq.set_xlim(0.0, 1.0)
    axq.set_ylim(-0.13, 0.23)
    axq.tick_params(labelsize=7.7)
    axq.set_xlabel(r"Fraction $q$ of the cloud in front of sphere", fontsize=7.7, labelpad=2)
    axq.set_ylabel(r"Relative error of the corrected $\Sigma$", fontsize=7.7, labelpad=1.5)
    axq.legend(loc="upper center", columnspacing=0.8, **LEG)

    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_deficit(stem="simsky5.w40", out="fig_deficit.pdf"):
    """Where the deficit sits: in the base image, not in the increments.

    The reconstruction is compared with the true surface density of the simulated sky, bin by bin in the true surface density, at
    the coarsest resolution and at the finest. The two curves lie almost on top of one another, which means that the increments
    the reconstruction adds between those resolutions are very nearly right and that the whole of the deficit is carried by the
    base image. That is what makes a single additive correction of the base image enough.
    """
    T13 = g("%s.surfdens.r13p5.cores.fits" % stem)
    D13 = g("hi.surface.density.r13p5.fits")
    D36 = g("hi.surface.density.r36p3.fits")
    sig = np.sqrt(BEAM ** 2 - FINE ** 2) * FWHM_TO_SIGMA / PIX
    T36 = gaussian_filter(T13, sig, mode="nearest")
    ok = (D36 > 0.0) & (D13 > 0.0) & np.isfinite(T13)

    edges = np.geomspace(3.0e21, 6.0e22, 19)
    mid = np.sqrt(edges[:-1] * edges[1:])
    r36, r13 = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = ok & (T36 >= lo) & (T36 < hi)
        r36.append(np.median(D36[m] / T36[m]) if m.sum() >= 200 else np.nan)
        r13.append(np.median(D13[m] / T13[m]) if m.sum() >= 200 else np.nan)
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.7))
    ax.plot(mid, r36, "-", color=COL3[0], label="%.1f″ base image" % BEAM)
    ax.plot(mid, r13, "--", color=COL3[1], label="%.1f″ reconstruction" % FINE)
    ax.axhline(1.0, **REF)
    ax.set_xscale("log")
    ax.set_xlim(3.0e21, 6.0e22)
#    ax.set_xticks([5.0e21, 1.0e22, 2.0e22, 5.0e22])
#    ax.set_xticklabels([r"$5\times10^{21}$", r"$10^{22}$", r"$2\times10^{22}$", r"$5\times10^{22}$"])
    clean_log(ax, "x")
    ax.set_ylim(0.45, 1.05)
    ax.set_xlabel(r"True surface density $\Sigma_{\rm T}$ (cm$^{-2}$)")
    ax.set_ylabel("Reconstructed / true", labelpad=1.5)
    ax.legend(loc="lower left", **LEG)
    fig.tight_layout(pad=0)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_shielding(table="/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260908_RT_BES_radmc3d/attenuation_law.txt", out="fig_shielding.pdf"):

    """The measured relation between the effective attenuating column and the dust temperature.

    The table is the one written by fit_attenuation_law.py: the effective attenuating column of every shell of every model in
    cm^-2, its dust temperature in K and its mass weight. With thousands of models the shells are too many to draw one by one,
    so they are shown as their mass-weighted density in the plane, with the weighted median and the 16th to 84th percentiles
    of the temperature in bins of the column, and the relation fitted to them, Eq. (7).
    """

    d = np.loadtxt(table, comments="!")
    # the cylinder radii, when fitted with the spheres, follow the sphere shells at the end of the table; the header says how many
    nsph = len(d)

    for line in open(table):
        if "sphere shells and" in line:
            nsph = int(line.split()[1])
            break

    cyl = d[nsph:]
    d = d[:nsph]
    s, te, w = d[:, 0], d[:, 1], d[:, 2]
    ok = (s > 0.0) & (te > 0.0) & (w > 0.0)
    s, te, w = s[ok], te[ok], w[ok]
    xe = np.linspace(18.0, 24.0, 121)
    ye = np.linspace(4.0, 24.0, 101)
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.8))

    if len(s) > 5000:
        h, _, _ = np.histogram2d(np.log10(s), te, bins=[xe, ye], weights=w)
        h = np.where(h > 0.0, h, np.nan)
        ax.pcolormesh(10.0 ** xe, ye, h.T, cmap="Greys", norm=LogNorm(vmin=np.nanmax(h) * 1.0e-4, vmax=np.nanmax(h)), rasterized=True, zorder=0)
    else:
        ax.plot(s, te, "o", ms=2.2, mfc="none", mec=COL3[0], mew=0.5, label="model shells")

    xb = np.arange(18.0, 24.01, 0.2)
    med, lo, hi, xm = [], [], [], []

    for a0, a1 in zip(xb[:-1], xb[1:]):
        m = (np.log10(s) >= a0) & (np.log10(s) < a1)
        if w[m].sum() <= 0.0 or m.sum() < 20:
            continue
        o = np.argsort(te[m])
        c = np.cumsum(w[m][o]) / w[m].sum()
        q = [te[m][o][np.searchsorted(c, p)] for p in (0.16, 0.50, 0.84)]
        lo.append(q[0])
        med.append(q[1])
        hi.append(q[2])
        xm.append(10.0 ** (0.5 * (a0 + a1)))

    if len(s) > 5000:
        ax.fill_between(xm, lo, hi, color=CCURVE[3], alpha=0.35, lw=0.0, label="spheres, 16th to 84th pct")
        ax.plot(xm, med, "-", color=CCURVE[1], label="spheres, median", lw=1.5)

    if len(cyl):
        ax.plot(cyl[:, 0], cyl[:, 1], "o", ms=3.5, mfc="none", mec=CCURVE[2], mew=0.6, label="cylinders")

    sx = np.geomspace(1.0e18, 1.0e24, 400)
    ax.plot(sx, law(sx), "--", color=CCURVE[0], label="fitted relation", lw=1.0)
    ax.axhline(TFLOOR, **REF)
    ax.text(3e19, TFLOOR + 0.35, r"$T_{\rm min}$", fontsize=7.7)
    ax.set_xscale("log")
    ax.set_yscale("linear")
    ax.set_xlim(2e19, 4.0e23)
    ax.set_ylim(6.0, 23.0)
    ax.tick_params(labelsize=7.7)
    ax.set_xlabel(r"Effective attenuating column $\Sigma_{\rm A}$ (cm$^{-2}$)", fontsize=7.7, labelpad=1.5)
    ax.set_ylabel(r"Dust temperature $T$ (K)", fontsize=7.7, labelpad=1.5)
    ax.legend(loc="upper right", **LEG ) ##, bbox_to_anchor=(1.0, 0.97)

    ax.set_box_aspect(1.0 / 1.25)

    fig.tight_layout(pad=0)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def _arrows(ax, paths, angles, r_out, r_in, col, tau_ref):
    """Draw the incoming radiation as arrows whose width follows how much of it survives the journey inward.

    paths holds the column crossed on the way in along each direction, in units of the half thickness of the structure, and
    tau_ref is the optical depth of that half thickness. The arrow is drawn thick and opaque where little material has been
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
    axis of a filament. The arrows are the incoming radiation, drawn thick where it has crossed little material, so that the eye
    sees which directions actually deliver energy. Lower panels: the effective attenuating column of Eq. (6) and the dust
    temperature that follows from it, along a line of sight through a layer, for three total surface densities.
    """

    cols = (CCURVE[0], CCURVE[1], CCURVE[2])
    ARROW = "#e07b00"          # one color for the radiation in all three geometries, distinct from the curves below
    fig = plt.figure(figsize=(7.1, 4.75))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.6, 1.0], hspace=0.12, wspace=-0.0, left=-0.0, right=0.8, bottom=0.195, top=0.945)
    ax = [fig.add_subplot(gs[0, i]) for i in range(3)]
    na, tau_ref = 36, 1.3
    ang = (np.arange(na) + 0.5) / na * 2.0 * np.pi
    sa = np.maximum(np.abs(np.sin(ang)), 1.0e-3)
    ca = np.maximum(np.abs(np.cos(ang)), 1.0e-3)

    # ---- plane-parallel layer: the body runs past both edges, the escape is short only along the normal

    ax[0].add_patch(plt.Rectangle((-1.30, -0.46), 2.60, 0.92, fc="0.8", ec="0.55", lw=0.5, zorder=0))
    _arrows(ax[0], np.minimum(1.0 / sa, 40.0), ang, 1.00, 0.13, ARROW, tau_ref)
    ax[0].set_title("Plane-parallel layer", fontsize=7.6, pad=4.0)

    # ---- sphere: every direction crosses the same column, none is favored

    ax[1].add_patch(plt.Circle((0.0, 0.0), 0.46, fc="0.8", ec="0.55", lw=0.5, zorder=0))
    _arrows(ax[1], np.ones(na), ang, 1.00, 0.13, ARROW, tau_ref)
    ax[1].set_title("Sphere", fontsize=7.6, pad=4.0)

    # ---- filament: short across the axis, long along it, so it lies between the other two

    ax[2].add_patch(plt.Rectangle((-1.30, -0.30), 2.60, 0.60, fc="0.8", ec="0.55", lw=0.5, zorder=0))
    _arrows(ax[2], np.minimum(np.minimum(1.0 / sa, 3.2 / ca), 40.0), ang, 1.00, 0.13, ARROW, tau_ref)
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
    gs_bottom = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, :], wspace=0.275)
    b1 = fig.add_subplot(gs_bottom[0, 0])
    b2 = fig.add_subplot(gs_bottom[0, 1])   
##    b1 = fig.add_subplot(gs[1, 0:2])
##    b2 = fig.add_subplot(gs[1, 2])
    u = (np.arange(512) + 0.5) / 512.0
    for sig, col, lab in ((1.0e23, cols[2], "10^{23}"), (1.0e22, cols[1], "10^{22}"), (1.0e21, cols[0], "10^{21}")):
        sa_u = attenuating_column(sig, "layer")
        b1.plot(u, sa_u, "-", color=col, label=r"$\Sigma = %s$ cm$^{-2}$" % lab)
        b1.plot(u, np.minimum(u, 1.0 - u) * sig, ":", color=col)
        b2.plot(u, law(sa_u), "-", color=col)
    b1.plot([], [], ":", color="1.0", label="dotted: shortest path")
    b1.set_yscale("log")
    b1.set_xlim(0.0, 1.0)
    b1.set_ylim(1.0e19, 2.0e23)
    b1.annotate("", xy=(0.02, 0.92), xytext=(0.13, 0.92), xycoords=b1.transAxes, textcoords=b1.transAxes,
                annotation_clip=False, arrowprops=dict(arrowstyle="-|>", color="black", lw=0.5))
    b1.text(0.14, 0.92, "to observer", fontsize=7.0, color="black", va="center", transform=b1.transAxes)

    for b, x, ha in ((b1, 0.97, "right"), (b2, 0.5, "center")):
        b.text(x, 0.92, "Plane-parallel layer ", fontsize=7.6, color="black", ha=ha, va="center", transform=b.transAxes)

    b1.tick_params(labelsize=7.8)
    b1.set_xlabel("Fractional depth along line of sight", fontsize=7.8)
    b1.set_ylabel(r"Attenuating column $\Sigma_{\rm A}$ (cm$^{-2}$)", fontsize=7.8, labelpad=1.5)
    b1.legend(frameon=False, fontsize=7, loc="lower center", ncol=1, handlelength=1.6, columnspacing=1.0, 
              borderpad=0.9, labelspacing=0.25)

    b2.axhline(TFLOOR, ls=":", color="black")
    b2.text(0.04, TFLOOR + 0.3, r"$T_{\rm min}$", fontsize=8, color="black")
    b2.set_xlim(0.0, 1.0)
    b2.set_ylim(6.5, 19.0)
    b2.tick_params(labelsize=7.8)
    b2.set_xlabel("Fractional depth", fontsize=7.8)
    b2.set_ylabel(r"Temperature $T$ (K)", fontsize=7.8, labelpad=1.5)

    b1.set_box_aspect(1.0 / 1.25)
    b2.set_box_aspect(1.0 / 1.25)

##    fig.tight_layout(pad=0.0, w_pad=1)
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

    The first group is the calibration, done once, which yields the temperature relation. The second is the calculation of
    Sect. 4.4, repeated for every trial surface density, which yields the table of the correction factor against the surface
    density that an observed image actually holds. The third is what is then done to an image.
    """
    import matplotlib.colors as mcolors
    
    groups = [
        (mcolors.to_rgb(CCURVE[0]), [
            "Radiative transfer models: dense cores embedded in clouds, and isolated filaments",
            r"Take the dust temperature $T$ and the attenuating column $\Sigma_{\rm A}$ of every shell and radius",
            r"Fit $T\,(\Sigma_{\rm A}) = T_{\rm min} + (T_0 - T_{\rm min})\,(1 + \Sigma_{\rm A}\,/\,\Sigma_{\bigstar})^{-\gamma}$, "
            "Eq. (5)"]),
        (mcolors.to_rgb(CCURVE[1]), [
            r"Trial line of sight of surface density $\Sigma$, a uniform sphere and a uniform layer",
            r"Attenuating column $\Sigma_{\rm A}(u)$ at every depth, averaged over directions and shapes, Eq. (4)",
            r"Temperature $T\,(\Sigma_{\rm A}(u))$ at every depth, from the fitted relation",
            r"Emission of all the depths added, giving the intensities $I_\lambda(\Sigma)$, Eq. (6)",
            r"One modified blackbody fitted as the image was, giving $\Sigma_{\rm F}$, Eq. (7)",
            r"Tabulate $f\,(\Sigma_{\rm F}) = \Sigma\,/\,\Sigma_{\rm F}$, what an image holds, Eq. (8)",
            r"Below $4{\times}10^{21}$ cm$^{-2}$, average the sphere with a radiative transfer cloud, blend to $10^{22}$, Sect. 4.3",
            r"Multiply by the field factor $g\,(\Sigma_{\rm F}, G_0)$ for the field of the cloud and hold at the maximum, Eq. (10)"]),
        (mcolors.to_rgb(CCURVE[2]), [
            r"Evaluate $f_{G_0}$ on the original image $\mathcal{D}$, create the correction factors image $\mathcal{F}$",
            r"Multiply the original image by the correction factors: $\mathcal{F}\cdot\mathcal{D}$, Eq. (11)"])]

    nbar = sum(len(t) for _, t in groups)
    ngap, ngrp = nbar - 1, len(groups) - 1
    fig, ax = plt.subplots(1, 1, figsize=(7.1, 0.435 * nbar + 0.30))
    fig.subplots_adjust(left=0.004, right=0.96, bottom=0.008, top=0.992)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    gap, extra = 0.54, 0.30                     # the gap between bars, and the additional gap between groups, in units of h
    h = 0.75 / (nbar + gap * ngap + extra * ngrp)
    y = 0.875
    for gi, (col, texts) in enumerate(groups):
        for k, t in enumerate(texts):
            y -= h
            _bar(ax, y, h, t, (col[0], col[1], col[2], 0.15 if k == 0 else 0.15))
            last = (gi == len(groups) - 1) and (k == len(texts) - 1)
            if not last:
                g = gap * h + (extra * h if k == len(texts) - 1 else 0.0)
                ax.annotate("", xy=(0.5, y - g + 0.003), xytext=(0.5, y - 0.005), xycoords=ax.transAxes, textcoords=ax.transAxes,
                            arrowprops=dict(arrowstyle="-|>", lw=0.7, color=(col[0], col[1], col[2], 1.0),
                                            mutation_scale=9, shrinkA=0.0, shrinkB=0.0))
                y -= g

#    fig.savefig(out)
    fig.savefig(out, bbox_inches='tight', pad_inches=-0.64)
    plt.close(fig)
    print("written %s" % out)


def fig_angular(out="fig_angular.pdf"):

    """How the angular average collapses a whole distribution of directions into one column.

    Left: the fraction of the heating radiation that survives the journey inward, against the cosine of the angle to the normal,
    for a parcel at mid-depth of a plane-parallel layer of three total surface densities. The dashed lines are the angular means
    of the curves, the mean intensity of Eq. (5), and the single column that gives the same mean intensity is the effective
    attenuating column of Eq. (6). Right: that column at mid-depth against the total surface density of the line of sight, for
    the two limiting shapes and their adopted average, with the shortest path for comparison.
    """

    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.85))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.165, top=0.965, wspace=0.3)
    
    mu = np.linspace(0.005, 1.0, 400)

    for sig, col, lab in ((3.0e20, CCURVE[0], r"3\times10^{20}"), (1.0e21, CCURVE[1], "10^{21}"),
                          (3.0e21, CCURVE[2], r"3\times10^{21}")):
        tr = np.exp(-KAPHEAT * 0.5 * sig / mu)
        ax[0].plot(mu, tr, "-", color=col, label=r"$\Sigma = %s$ cm$^{-2}$" % lab)
        ax[0].axhline(np.trapezoid(tr, mu), ls="--", lw=0.5, color=col)

    ax[0].set_xlim(0.0, 1.0)
    ax[0].set_ylim(0.01, 1.0)
    ax[0].tick_params(labelsize=7.7)
    ax[0].set_xlabel(r"Direction of incidence $\cos i$", fontsize=7.7, labelpad=3.5)
    ax[0].set_ylabel(r"Transmitted fraction $\exp({-\kappa\,\Sigma_{\rm A}(\hat n))}$", fontsize=7.7, labelpad=1.5)
    ax[0].legend(loc="upper left", **LEG)    ##, bbox_to_anchor=(0.02, 0.98)

    sg = np.geomspace(1.0e20, 1.0e24, 60)
    mid = 512 // 2

    for geom, col, lab in (("layer", CCURVE[0], "plane-parallel layer"), ("mean", CCURVE[1], "average (adopted)"),
                           ("sphere", CCURVE[2], "sphere")):
        y = np.array([attenuating_column(v, geom)[mid] for v in sg])
        ax[1].plot(sg, y, "-", color=col, label=lab)

    ax[1].plot(sg, 0.5 * sg, ":", color=CCURVE[3], lw=3.0, label="shortest path")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlim(1.0e20, 1.0e24)
    ax[1].set_ylim(5.0e19, 5.0e23)
    ax[1].tick_params(labelsize=7.7)
    ax[1].set_xlabel(r"Surface density $\Sigma$ of the line of sight (cm$^{-2}$)", fontsize=7.7, labelpad=1.5)
    ax[1].set_ylabel(r"Attenuating column $\Sigma_{\rm A}$ (cm$^{-2}$) at mid-depth", fontsize=7.7, labelpad=1.5)
    ax[1].legend(loc="upper left", **LEG)    ##, bbox_to_anchor=(0.02, 0.98)

    for a in ax:
        a.set_box_aspect(1.0 / 1.25)

    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_mixing(out="fig_mixing.pdf"):
  
    """Why one modified blackbody fitted to a non-isothermal line of sight returns too little surface density.

    Left: the spectral energy distribution of a line of sight through a layer of 3e22 cm^-2, the sum over depth of the emission of
    material at every temperature, with the best single modified blackbody fitted to it in the four wavebands the reconstruction
    uses, weighted by their calibration uncertainties. Right: the relative residuals of that fit, waveband by waveband, for three
    total surface densities, with the fraction of the true surface density that the fit recovers.
    """

    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.85))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.165, top=0.965, wspace=0.18)
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

    ax[0].plot(wl, smooth, "-", color=CCURVE[0], label="mixed emission")
    ax[0].plot(wl, sfit, "--", color=CCURVE[2], label=r"single fit, $T = %.1f$ K" % tfit)
    ax[0].plot(BANDS, obs, "o", ms=3.0, mfc="none", mec=CCURVE[0], mew=0.7)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlim(120.0, 650.0)
    ax[0].set_ylim(30.0, 200.0)
    ax[0].set_xticks(BANDS)
    ax[0].set_xticklabels(["160", "250", "350", "500"])
    ax[0].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax[0].tick_params(labelsize=7.5)
    ax[0].set_xlabel(r"Wavelength ($\mu$m)", fontsize=7.5)
    ax[0].set_ylabel(r"Intensity $I_\nu$ (MJy sr$^{-1}$)", fontsize=7.5, labelpad=1.5)
    clean_log(ax[0], "y")
    ax[0].legend(loc="lower center", **LEG)
##    ax[0].set_title(r"$\Sigma = 3\times10^{22}$ cm$^{-2}$, recovered %.3f" % rec, fontsize=FS_TITLE, pad=3)
    ax[0].text(280.0, 50.0, r"$\Sigma = 3{\times}10^{22}$ cm$^{-2}$, recovered: %.3f" % rec, fontsize=7.5, color="black", ha="center")

    for sig, col, lab in ((3.0e21, CCURVE[0], r"$3{\times}10^{21}$"), (1.0e22, CCURVE[1], r"$1{\times}10^{22}$"),
                          (3.0e22, CCURVE[2], r"$3{\times}10^{22}$")):
        _, o, f, _, r = one(sig)
        ax[1].plot(BANDS, 100.0 * (o - f) / o, "o--", ms=3.0, mfc="none", mew=0.7, color=col,
                   label=r"%s cm$^{-2}$, recovered: %.3f" % (lab, r))

    ax[1].axhline(0.0, **REF)
    ax[1].set_xscale("log")
    ax[1].set_xlim(140.0, 570.0)
    ax[1].set_xticks(BANDS)
    ax[1].set_xticklabels(["160", "250", "350", "500"])
    ax[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax[1].set_ylim(-9.0, 24.0)
    ax[1].tick_params(labelsize=7.5)
    ax[1].set_xlabel(r"Wavelength ($\mu$m)", fontsize=7.5)
    ax[1].set_ylabel("Relative residual of the fit (%)", fontsize=7.5, labelpad=1.8)
    ax[1].legend(loc="upper center", frameon=False, fontsize=FS_LEG, handlelength=1.6, labelspacing=0.2)

    for a in ax:
        a.set_box_aspect(1.0 / 1.25)

    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def _relation(eta, ns=36):

    """Fitted surface density and the restoring factor, for one value of the superseded geometry parameter.

    This and fig_orientation below belong to the earlier scheme, in which a single parameter eta set the deepest shielding of a
    line of sight. That parameter no longer exists: the orientation average is now inside the angular average of
    attenuating_column. They are left here only until it is decided what replaces fig_orientation, and they are not called by
    anything else.
    """

    u = (np.arange(256) + 0.5) / 256.0
    dep = 2.0 * np.minimum(u, 1.0 - u)
    tg = np.arange(5.0, 30.0, 0.02)
    mod = np.array([[cnu(w) * planck(t, w) for w in BANDS] for t in tg])
    sg = np.geomspace(3.0e20, 3.0e23, ns)
    rec, fs = np.empty(ns), np.empty(ns)

    for i, S in enumerate(sg):
        d = np.array([cnu(w) * S * np.mean(planck(law(eta * S * dep), w)) for w in BANDS])
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
    makes the angle i. The gray band is the isotropic probability density of
    that angle, proportional to its sine, drawn on an arbitrary scale so that
    the eye can weight the curves correctly: angles near the line of sight are
    rare. Right: the correction factor averaged over isotropic orientations,
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
    between them; the adopted factor is the average of the first two. The dotted curve replaces the angular average of Eq. (6)
    by the column along the shortest path, with its geometric factor of one quarter.
    """

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.8))

    for geom, col, lab in (("layer", CCURVE[0], "plane-parallel layer"), ("cube", CCURVE[1], "cube"),
                           ("sphere", CCURVE[2], "sphere"), ("mean", 'black', "average of sphere and layer")):
        fit, f = relation(geom)
        ax.plot(fit, f, "--" if geom == "mean" else "-", lw=1.0 if geom == "mean" else None, color=col, label=lab)
    fit, f = relation("short")
    ax.plot(fit, f, ":", color="orange", label="shortest path", lw=1.2)
    # the relation adopted in hicorrect2: below 10^21.6 cm^-2 the layer is replaced by the RADMC-3D cirrus

    import hicorrect2 as h2

    fit, f = h2.build_relation([160.0, 250.0, 350.0, 500.0])
    # the adopted factor never decreases with the surface density (held at its maximum), Eq. (12) with g = 1
    ax.plot(fit, np.maximum.accumulate(f), "-", color=CCURVE[3], lw=1.5, label="adopted, calibrated on a radiative transfer cloud")
    # the adopted factor for two other fields, with the field factor g of Eq. (12)
    for g0, ls in ((5.0, (0, (5, 2))), (33.0, (0, (1, 1.5)))):
        fg = np.maximum.accumulate(f * h2.field_factor(np.log10(fit), g0))
        ax.plot(fit, fg, linestyle=ls, color=CCURVE[3], lw=1.1, label=r"adopted for $G_0 = %g$" % g0)
    ax.set_xscale("log")
    ax.set_xlim(1e21, 3.0e23)
    ax.set_ylim(0.9, 2.98)
    ax.tick_params(labelsize=7.7)
    ax.set_xlabel(r"Surface density $\Sigma$ (cm$^{-2}$), uncorrected", fontsize=7.7, labelpad=1.5)
    ax.set_ylabel(r"Correction factor $f$", fontsize=7.7, labelpad=1.5)
    ax.legend(loc="upper left", **LEG)    ##, bbox_to_anchor=(0.02, 0.98)

    ax.set_box_aspect(1.0 / 1.25)

    fig.tight_layout(pad=0)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_masses(stem="simsky4.w40", out="fig_masses.pdf"):
    """How accurately a mass is recovered, against the aperture over which it is measured.

    Left: the mass of a core, against the radius of the aperture in units of the Bonnor-Ebert radius, with the background taken
    between 2 and 2.5 radii. Right: the linear density of a filament, against the half-width of the window over which the
    profile is integrated, with the background level taken from the fit of A r^-alpha + B. The core mass is recovered to a few
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
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.165, top=0.915, wspace=0.25)
    for a, xx, dat, xl, ti in ((ax[0], aps, core, r"Aperture radius ($R_{\rm BE}$)", "Cores"),
                               (ax[1], wins, fil, "Half-width of the window (″)", "Filaments")):
        a.plot(xx, dat["u"], "-", color=CCURVE[2], label="uncorrected")
        a.plot(xx, dat["c"], "-", color=CCURVE[0], label="corrected")
        a.axhline(1.0, **REF)
        a.set_xlim(xx[0], xx[-1]-0.1)
        a.set_ylim(0.35, 1.65)
        a.set_xlabel(xl)
        a.set_ylabel("Fraction of the true mass", labelpad=1.5)
##        a.set_title(ti, fontsize=FS_TITLE, pad=3)
        a.text(xx[-1]/2, 1.3, ti, fontsize=7.3)
        a.legend(loc="lower center", **LEG)

    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_base(stem="simsky4.w40", out="fig_base.pdf"):
    """What the choice of the base image costs, measured on the simulated sky.

    The factor is evaluated on the reconstruction at one resolution and the difference is added to the image at the finest
    resolution, which is Eq. (12) with the base image taken in turn at each available resolution. Left: the mass of a core
    recovered within three apertures. Right: the shape of a filament, as ratios to the true image so that one axis serves for
    all four quantities. The coarsest base gives the best masses and the poorest shapes, and the gain in shape is exhausted by
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
        """Index, half-maximum width, linear density and crest amplitude of one filament, with the level fitted."""
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
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.165, top=0.915, wspace=0.28)

    for k, ap in enumerate((0.5, 1.0, 1.5)):
        ax[0].plot(beams, [coremass(images[b], ap) for b in beams], "o-", ms=3.0, mfc="none", mew=0.7, color=CCURVE[k],
                   label=r"aperture $%.1f\,R_{\rm BE}$" % ap)
        ax[0].axhline(coremass(D, ap), ls=":", lw=0.7, color=CCURVE[k])

    ax[0].plot([], [], ls=":", lw=0.7, color="black", label="uncorrected")
    ax[0].axhline(1.0, ls="--", lw=0.5, color="black")
    ax[0].set_xlim(39.0, 11.0)
    ax[0].set_ylim(0.40, 1.20)
    ax[0].set_xlabel("Resolution of the base image (″)")
    ax[0].set_ylabel("Fraction of the true core mass", labelpad=1.5)
    ax[0].legend(loc="center right", **LEG)
##    ax[0].set_title("Cores", fontsize=FS_TITLE, pad=3)
    ax[0].text(30.0, 0.8, "Cores", fontsize=7.7, color="black", ha="center")

    names = [("index " + r"$\alpha$", 0), ("half-maximum width", 1), ("linear density", 2), ("crest amplitude", 3)]
    styl = ["o-", "s-", "^-", "v-"]

    for (lab, j), st, col in zip(names, styl, (CCURVE[0], CCURVE[1], CCURVE[2], CCURVE[3])):
        ax[1].plot(beams, [fil[b][j] / ref[j] for b in beams], st, ms=3.0, mfc="none", mew=0.7, color=col, label=lab)
        ax[1].axhline(filu[j] / ref[j], ls=":", lw=0.7, color=col)

    ax[1].plot([], [], ls=":", lw=0.7, color="black", label="uncorrected")
    ax[1].axhline(1.0, ls="--", lw=0.5, color="black")
    ax[1].set_xlim(39.0, 11.0)
    ax[1].set_ylim(0.5, 2.2)
    ax[1].set_xlabel("Resolution of the base image (″)")
    ax[1].set_ylabel("Ratio to the true image", labelpad=1.5)
    ax[1].legend(loc="upper center", ncol=2, columnspacing=1.0, **LEG)
##    ax[1].set_title("Filaments", fontsize=FS_TITLE, pad=3)
    ax[1].text(30.0, 1.35, "Filaments", fontsize=7.7, color="black", ha="center")

    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_residual(stem="simsky4.w40", out=None, region=slice(150, 1650)):
    """Where the reconstruction and the correction depart from the truth, as maps and as profiles.

    Upper panels: the relative residual of the reconstruction at the finest resolution before and after the correction. The
    uncorrected image is too low almost everywhere, most of all in the dense material. The corrected one is right to a few per
    cent over the diffuse field, but leaves a ring of over-correction around every compact structure and along the flanks of every
    filament. Lower panels: the same thing as profiles, averaged over the filaments and over the injected cores, which shows that
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
    gs = fig.add_gridspec(2, 6, height_ratios=[1.45, 1.0], hspace=0.1, wspace=0.53, left=0.0, right=0.905, bottom=0.085, top=0.955)
    for k, (x, lab) in enumerate(((rd, "$hires$ 13.5″ surface density"), (rc, "$hires$ 13.5″ surface density, corrected"))):
        a = fig.add_subplot(gs[0, 3 * k:3 * k + 3])
        n = x[region, region].shape[0]
        ext = [-n * PIX / 120.0, n * PIX / 120.0] * 2
        im = a.imshow(x[region, region], origin="lower", extent=ext, cmap=SD_CMAP, vmin=-0.45, vmax=0.45)
        a.set_title(lab, fontsize=7.3, pad=3)
        a.set_xticks([])
        a.set_yticks([])
        a.tick_params(labelsize=7.8)
        cb = colorbar_beside(fig, a, im)
        cb.ax.tick_params(labelsize=7.3, pad=1.7)

##    cax = fig.add_axes([0.915, 0.475, 0.016, 0.480])
##    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Relative accuracy of surface density", labelpad=1.8)

    prof = np.nanmedian(np.where(ok, T, np.nan), axis=0)
    off = np.arange(-70, 71)
    pk = [i for i in range(80, len(prof) - 80)
          if prof[i] == np.nanmax(prof[i - 30:i + 31]) and prof[i] > 1.5 * np.nanmedian(prof)]

    # ---- the profiles themselves, stacked over the filaments
    p0 = fig.add_subplot(gs[1, 0:2])
    for M, col, lab, lw in ((T, CCURVE[1], "true", 1.0), (D, CCURVE[2], "uncorrected", 1.0), (C, CCURVE[0], "corrected", 1.0)):
        st = np.array([[np.nanmedian(M[:, x + d]) for d in off] for x in pk])
        q = np.nanmedian(st, axis=0)
        q = q - np.median(np.concatenate([q[:12], q[-12:]]))
        p0.plot(off * PIX, q, "-", lw=lw, color=col, label=lab)
    p0.set_yscale("log")
    p0.set_xlim(off[0] * PIX - 20, off[-1] * PIX + 20)
    p0.set_ylim(3.0e20, 5.0e22)
    p0.tick_params(labeltop=False, labelright=False, labelsize=7.8)
    p0.set_xlabel("Offset from crest (″)", fontsize=7.7, labelpad=1.5)
    p0.set_ylabel(r"Surface density $\Sigma$ (cm$^{-2}$)", labelpad=1, fontsize=7.8)
    p0.legend(frameon=False, fontsize=6.6, loc="upper left", handlelength=1.6, labelspacing=0.25)
    p0.text(0.5 * FINE + 8.0, 0.04, "%.1f″ beam" % FINE, fontsize=7, color="black", va="bottom",
            transform=p0.get_xaxis_transform())
    for v in (-0.5 * FINE, 0.5 * FINE):
        p0.axvline(v, ls="--", lw=0.5, color="0.3")

    # ---- across the filaments, stacked
    b = fig.add_subplot(gs[1, 2:4])
    for r, col, lab in ((rd, CCURVE[2], "uncorrected"), (rc, CCURVE[0], "corrected")):
        stack = np.array([[np.nanmedian(r[:, x + d]) for d in off] for x in pk])
        b.plot(off * PIX, np.nanmedian(stack, axis=0), "-", lw=1.0, color=col, label=lab)
    b.axhline(0.0, ls=":", lw=0.8, color="0.5")
    for v in (-0.5 * FINE, 0.5 * FINE):
        b.axvline(v, ls="--", lw=0.5, color="0.3")
    b.set_xlim(off[0] * PIX - 20, off[-1] * PIX + 20)
    b.set_ylim(-0.55, 0.55)
    b.tick_params(labeltop=False, labelright=False, labelsize=7.8)   
    b.set_xlabel("Offset from crest (″)", fontsize=7.8, labelpad=1.5)
##    b.set_ylabel("(reconstructed $-$ true) / true")
    b.legend(frameon=False, fontsize=6.6, loc="lower right", handlelength=1.6, labelspacing=0.25)
    b.set_title("Relative residuals (%d filaments)" % len(pk), fontsize=7.3, pad=4)
    p0.set_title("Profiles of %d filaments" % len(pk), fontsize=7.3, pad=4)

    # ---- around the injected cores, stacked
    c = fig.add_subplot(gs[1, 4:6])
    cat = [l.split() for l in open("%s.cores.cores.txt" % stem) if not l.startswith("#") and l.strip()]
    edges = np.linspace(0.0, 3.0, 16)
    mid = 0.5 * (edges[:-1] + edges[1:])
    ny, nx = T.shape
    for r, col, lab in ((rd, CCURVE[2], "uncorrected"), (rc, CCURVE[0], "corrected")):
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
    c.text(0.5 * FINE / rbe + 0.06, 0.93, "%.1f″ beam" % FINE, fontsize=7, color="black", va="top",
           transform=c.get_xaxis_transform())
    c.set_xlim(0.0, 3.0)
    c.set_ylim(-0.55, 0.55)
    c.tick_params(labeltop=False, labelright=False, labelsize=7.8)   
    c.set_xlabel(r"Offset from center ($r\,R^{-1}_{\rm BE}$)", fontsize=7.8, labelpad=1.3)
##    c.set_yticklabels([])   
##    c.set_ylabel("(reconstructed $-$ true) / true")
    c.legend(frameon=False, fontsize=6.6, loc="lower right", handlelength=1.6, labelspacing=0.25)
    c.set_title("Relative residuals (%d cores)" % len(rows), fontsize=7.3, pad=4)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_simsky_bands(stem="simsky5.w40", out=None, region=slice(400, 1400)):
    """The simulated sky as it is observed: the five wavebands used in the figure and the true surface density.

    The panels are arranged by resolution, from the sharpest at the upper left to the coarsest at the lower right, with the true
    surface density in the place its resolution would occupy. Each panel carries its own logarithmic scale, because the sky is
    two orders of magnitude fainter at 70 than at 250 micron and one common scale would show nothing at the ends. The 100 micron
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
            cmap, lab, unit = SD_CMAP, "True surface density, 13.5″", r"$N_{\rm H_2}$ (cm$^{-2}$)"
        else:
            x = g("%s.%s.%s.cores.fits" % (stem, band, rtag))[region, region]
            cmap, unit = SD_CMAP, r"$I_\nu$ (MJy sr$^{-1}$)"
            lab = r"%s{\um}, $%.1f^{\prime\prime}$" % (band, beam)
            lab = "%s $\\mu$m, %.1f″" % (band, beam)
        v = x[np.isfinite(x) & (x > 0.0)]
        vmin, vmax = np.percentile(v, 1.0), np.percentile(v, 99.9)
        n = x.shape[0]
        ext = [-n * PIX / 120.0, n * PIX / 120.0] * 2
        im = a.imshow(x, origin="lower", extent=ext, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
        a.set_title(lab, fontsize=FS_TITLE, pad=3)
        a.set_xticks([])
        a.set_yticks([])
        cb = colorbar_beside(fig, a, im)
        cb.set_label(unit, fontsize=FS_CB, labelpad=1)
        cb.ax.tick_params(labelsize=FS_CB, pad=1)
        cb.ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_aquila(out="fig_aquila.pdf"):

    """The Aquila field before and after the correction, with the factor applied and the dust temperature that follows.

    The two surface density panels share one logarithmic scale, so the change can be read directly. The factor is the one
    interpolated onto the base image at the coarsest resolution and added to both resolutions; the temperature is the one implied
    by the corrected surface density at the coarsest resolution.
    """

    dd = g("/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/+hires_cutout/+results/hi.surface.density.r%s.fits" % ("%.1f" % FINE).replace(".", "p"))
    cc = g("/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/+hires_cutout/+results/hi.surface.density.r%s.correct.fits" % ("%.1f" % FINE).replace(".", "p"))
    ff = g("/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/+hires_cutout/+results/hi.surface.density.r%s.cfactor.fits" % ("%.1f" % FINE).replace(".", "p"))
    tt = g("/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/+hires_cutout/+results/hi.temperature.r%s.correct.fits" % ("%.1f" % FINE).replace(".", "p"))

    lo, hi = np.nanpercentile(dd[dd > 0.0], 0.5), np.nanpercentile(cc[cc > 0.0], 99.9)
    n = dd.shape[0]
    ext = [-n * PIX / 120.0, n * PIX / 120.0] * 2

    fig = plt.figure(figsize=(7.1, 6.9))
    gs = fig.add_gridspec(2, 2, wspace=0.075, hspace=0.0, left=0.010, right=0.940, bottom=0.015, top=0.955)

    panels = [(dd, r"Surface density $\mathcal{D}$ (cm$^{-2}$)", SD_CMAP, LogNorm(vmin=lo, vmax=hi), ""),
              (cc, r"Corrected surface density $\mathcal{D}_{\rm c}$ (cm$^{-2}$)", SD_CMAP, LogNorm(vmin=lo, vmax=hi), ""),
              (ff, r"Correction factors $\mathcal{F}$", SD_CMAP, None, ""),
              (tt, r"Corrected dust temperature $T$ (K)", SD_CMAP, None, "")]

#    ax.set_xticks([0.1, 0.2, 0.5, 1.0, 2.0])
#    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2"])

    for k, (x, lab, cmap, nor, unit) in enumerate(panels):
        a = fig.add_subplot(gs[k // 2, k % 2])
        v = np.where(np.isfinite(x) & (x > 0.0), x, np.nan)          # the blank border must not set the scale
        kw = dict(norm=nor) if nor is not None else dict(vmin=np.nanpercentile(v, 0.5), vmax=np.nanpercentile(v, 99.5))
        im = a.imshow(v, origin="lower", extent=ext, cmap=cmap, **kw)
        a.set_title(lab, fontsize=FS_TITLE, pad=4)
        a.set_xticks([])
        a.set_yticks([])
        cb = colorbar_beside(fig, a, im, "3%")
        if k == 0: cb.set_label(unit, fontsize=FS_CB, labelpad=-11)
        if k == 1: cb.set_label(unit, fontsize=FS_CB, labelpad=-11)
        if k == 2: cb.set_label(unit, fontsize=FS_CB, labelpad=-4)
        if k == 3: cb.set_label(unit, fontsize=FS_CB, labelpad=-4)
        cb.ax.tick_params(labelsize=7.5, pad=1.5)
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


def fig_grid_validation(table="/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260910_RT_BES_radmc3d_subdivision/validate3_finRTlaw.txt", 
                        out="fig_grid_validation.pdf"):

    """The correction measured on the radiative transfer models, and how much a misjudged footprint costs.

    Left: the relative residual of the reconstruction at 13.5 arcsec against the distance from the center of the core, the median
    over the 1815 nodes of the subdivision grid, before the correction and after it with the base image at 36.3 and at 13.5
    arcsec. Right: the core mass, measured over the footprint of the core less the background just outside it, divided by the true
    mass measured the same way with the correct footprint, against the error of the footprint in units of its equivalent radius.
    The lines are the medians over the nodes and the band the 16th to 84th percentiles for the correction with the 36.3 arcsec
    base. Both panels are read from the output of validate_calib3.py.
    """

    names = {"uncorrected": "uncorrected", "hicorrect2, coarse": "corrected, $f$ at 36.3″", "hicorrect2, fine": "corrected, $f$ at 13.5″"}
    cols = {"uncorrected": CCURVE[2], "hicorrect2, coarse": CCURVE[1], "hicorrect2, fine": CCURVE[0]}
    lines = open(table).read().splitlines()
    mass0, foot, prof = {}, {}, []
    svals, sect = None, None

    for ln in lines:
        t = ln[1:].strip() if ln.startswith("!") else ln.strip()

        if t.startswith("core mass over the mask"):
            sect = "mass"
        elif t.startswith("the same core mass with the footprint"):
            sect = "foot"
        elif t.startswith("ring minus center"):
            sect = None
        elif t.startswith("s ") and sect == "foot":
            svals = [float(v) for v in t.split()[1:]]
        elif not ln.startswith("!") and t:
            prof.append([float(v) for v in t.split()])
        else:
            for key in names:
                if t.startswith(key + " ") and t[len(key):len(key) + 3].strip()[:1] in "0123456789+-[":
                    rest = t[len(key):].split()
                    if sect == "mass" and not rest[0].startswith("["):
                        mass0[key] = [float(v) for v in rest[:3]]
                    elif sect == "foot":
                        v = [float(x.strip("[]")) for x in rest]
                        foot[key] = np.array(v).reshape(-1, 3)

    prof = np.array(prof)

    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.165, top=0.915, wspace=0.26)

    for k, key in enumerate(names):
        ax[0].plot(prof[:, 0], prof[:, 1 + k], "-", color=cols[key], label=names[key])

    ax[0].axhline(0.0, **REF)
    ax[0].axvline(0.5 * FINE, ls="--", lw=0.5, color="0.3")
    ax[0].text(0.5 * FINE + 2.0, 0.09, "%.1f″ beam" % FINE, fontsize=7, color="black", va="top", transform=ax[0].get_xaxis_transform())
    ax[0].set_xlim(0.0, 124.0)
    ax[0].set_ylim(-0.45, 0.12)
    ax[0].tick_params(labelsize=7.7)
    ax[0].set_xlabel("Distance from center (″)", fontsize=7.7, labelpad=1.5)
    ax[0].set_ylabel(r"Relative error of $\Sigma(r\,)$", fontsize=7.7, labelpad=1.5)
    ax[0].text(0.5 * FINE + 3.0, 0.06, "Radial profiles (1815 models)", fontsize=7.3)
    ax[0].legend(loc="lower right", **LEG)

    s = np.array(svals[:2] + [0.0] + svals[2:])

    for key in names:
        f = foot[key]
        med = np.concatenate([f[:2, 0], [mass0[key][0]], f[2:, 0]])
        lo = np.concatenate([f[:2, 1], [mass0[key][1]], f[2:, 1]])
        hi = np.concatenate([f[:2, 2], [mass0[key][2]], f[2:, 2]])
        ax[1].plot(s, med, "o-", ms=3.0, mfc="none", mew=0.7, color=cols[key], label=names[key])

        if key == "hicorrect2, coarse":
            ax[1].fill_between(s, lo, hi, color=cols[key], alpha=0.2, lw=0.0)

    ax[1].axhline(1.0, **REF)
    ax[1].axvline(0.0, ls="--", lw=0.5, color="black")
    ax[1].set_xlim(-0.55, 1.05)
    ax[1].set_ylim(0.2, 2.0)
    ax[1].tick_params(labelsize=7.7)
    ax[1].set_xlabel("Footprint error (equivalent radii)", fontsize=7.7, labelpad=1.5)
    ax[1].set_ylabel("Accuracy ratio of core mass", fontsize=7.7, labelpad=1.5)
    ax[1].text(-0.42, 1.806, "Core masses (1750 models)", fontsize=7.3)
    ax[1].legend(loc="lower right", **LEG)

    for a in ax:
        a.set_box_aspect(1.0 / 1.25)

    fig.tight_layout(pad=0.0, w_pad=2)
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
    ax.axhline(1.0, **REF)
    ax.plot(r[:, 0], r[:, 1], "o", ms=3.0, mfc="none", mec=PLASMA(30), mew=0.7, label="uncorrected")
    ax.plot(r[:, 0], r[:, 2], "s", ms=3.0, mfc="none", mec=PLASMA(150), mew=0.7, label="corrected")
    ax.set_xscale("log")
    ax.set_xticks([0.1, 1.0])
    ax.set_xticklabels(["0.1", "1.0"])
    clean_log(ax, "x")
    ax.set_xlabel(r"True core mass ($M_\odot$)")
    ax.set_ylabel("Recovered mass fraction", labelpad=1.5)
    ax.set_ylim(0.05, 1.7)
    ax.legend(loc="lower center", **LEG)
    fig.tight_layout(pad=0.4)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)
    print("   hires     median %.3f" % np.median(r[:, 1]))
    print("   corrected median %.3f" % np.median(r[:, 2]))


def fig_cylinders(dirs=None, out="fig_cylinders.pdf"):

    """The seven isolated RADMC-3D cylinders reconstructed by hires and corrected.

    Left: the profile across the cylinder with 3e22 cm^-2 on its axis, true (convolved to 13.5 arcsec), uncorrected and
    corrected. Right: the corrected over the true surface density across all seven cylinders. Each profile is the median over the
    middle rows of the image, where the cylinder is uniform along its axis.
    """

    import hicorrect2 as h2
    from scipy.ndimage import gaussian_filter1d

    if dirs is None:
        dirs = [("/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261005_RT_cyl_sd3e21_isolated", r"$3{\times}10^{21}$"), 
                ("/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261005_RT_cyl_sd1e22_isolated", r"$1{\times}10^{22}$"),
                ("/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261005_RT_cylinder_isolated"  , r"$3{\times}10^{22}$"), 
                ("/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261005_RT_cyl_sd5e22_isolated", r"$5{\times}10^{22}$"),
                ("/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261005_RT_cyl_sd1e23_isolated", r"$1{\times}10^{23}$"), 
                ("/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261005_RT_cyl_p2_isolated"    , r"$3{\times}10^{22}$, $p=2$"),
                ("/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261005_RT_cyl_w80_isolated"   , r"$3{\times}10^{22}$, $H=80″$")]

    fitted, factor = h2.build_relation([160.0, 250.0, 350.0, 500.0])
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))       ##, gridspec_kw={"wspace": 0.35}
##    cols = plt.cm.viridis(np.linspace(0.0, 0.9, len(dirs)))
    cols = ["red", "green", "blue", "cyan", "magenta", "brown", "orange"]
    off = np.arange(-100, 101)

    for k, (d, lab) in enumerate(dirs):
        t0 = g(os.path.join(d, "mc.surfdens.fits"))
        u = g(os.path.join(d, "+hires/+results/hi.surface.density.r13p5.fits"))
        ny = u.shape[0]
        rows = slice(ny // 2 - 400, ny // 2 + 401)
        pt0 = np.nanmedian(t0[rows], 0)
        pu = np.nanmedian(u[rows], 0)
        x0 = int(np.nanargmax(pt0))
        pt = gaussian_filter1d(pt0, 13.5 / 2.3548 / PIX, mode="nearest")
        pc = pu * np.interp(np.log10(np.maximum(pu, 1.0e18)), np.log10(fitted), factor)

        if d.startswith("/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261005_RT_cylinder_"):
            for p, col, l in ((pt, CCURVE[1], r"true $3{\times}10^{22}$"), (pu, CCURVE[2], "uncorrected"), (pc, CCURVE[0], "corrected")):
                ax[0].plot(off * PIX, p[x0 + off], "-", color=col, label=l)

        ax[1].plot(off * PIX, pc[x0 + off] / pt[x0 + off] - 1, "-", color=cols[k], label=lab)

##    for a in ax:
##        a.tick_params(labeltop=False, labelright=False)

    ax[0].set_yscale("log")
    ax[0].set_xlim(-249, 249)
    ax[0].set_ylim(9.0e20, 4.0e22)
    ax[0].tick_params(labelsize=7.0)
    ax[0].set_xlabel("Offset from crest (″)", fontsize=7.0, labelpad=1.5)
    ax[0].set_ylabel(r"Surface density $\Sigma$ (cm$^{-2}$)", fontsize=7.0, labelpad=1.5)
    ax[0].legend(loc="upper left", fontsize=6.2)  # **LEG

    ax[1].axhline(1.0, **REF)
    ax[1].set_xlim(-249, 249)
    ax[1].set_ylim(-0.22, 0.07)
    ax[1].tick_params(labelsize=7.0)
    ax[1].set_xlabel("Offset from crest (″)", fontsize=7.0, labelpad=1.5)
    ax[1].set_ylabel(r"Relative error of $\Sigma$", fontsize=7.0, labelpad=1.5)
    ax[1].legend(loc="lower left", ncol=1, fontsize=6.2)  #**LEG

    for a in ax:
        a.set_box_aspect(1.0 / 1.25)
    
    fig.tight_layout(pad=0.0, w_pad=0.7)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_simsky6(sky="simsky6_RT_enh_nocrater", cirrus="simsky6_RT_enh_cirronly", out="fig_simsky6.pdf"):

    """The simulated sky with RADMC-3D cirrus, filaments and cores, corrected with hicorrect2.

    Upper panels: the true surface density, and the relative error of the corrected surface density, (corrected - true) / true,
    for the full sky and for the cirrus alone. Lower panels: that error across each filament, median over the outer rows of the
    field, where the cirrus is a few 10^21 cm^-2, and over its central rows, where the enhanced cirrus reaches 2e22 cm^-2. The
    cores are masked in the profiles.
    """

    import hicorrect2 as h2

    fitted, factor = h2.build_relation([160.0, 250.0, 350.0, 500.0])
    cor = lambda u: u * np.interp(np.log10(np.maximum(u, 1.0e18)), np.log10(fitted), factor)

    T = g(os.path.join(sky, "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_simsky6_RT_enh_cirrdepth=2_nocrater/simsky6.w40.surfdens.r13p5.cores.fits"))
    U = g(os.path.join(sky, "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_simsky6_RT_enh_cirrdepth=2_nocrater/+hires/+results/hi.surface.density.r13p5.fits"))
    Tc = g(os.path.join(cirrus, "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_simsky6_RT_enh_cirrdepth=2_cirronly/simsky6.w40.surfdens.r13p5.fits"))
    Uc = g(os.path.join(cirrus, "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_simsky6_RT_enh_cirrdepth=2_cirronly/+hires/+results/hi.surface.density.r13p5.fits"))

    R = cor(U) / T - 1.0
    Rc = cor(Uc) / Tc - 1.0
    cat = [l.split() for l in open(os.path.join(sky, "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_simsky6_RT_enh_cirrdepth=2_nocrater/simsky6.w40.cores.cores.txt")) if not l.startswith("#") and l.strip()]
    mask = np.zeros(T.shape, bool)

    for c in cat:
        y, x, rbe = int(c[1]), int(c[2]), float(c[4])
        r = int(max(2.5 * rbe, 30.0) / PIX) + 3
        mask[max(0, y - r):y + r + 1, max(0, x - r):x + r + 1] = True
    Rm = np.where(mask, np.nan, R)

    fig = plt.figure(figsize=(7.1, 4.7))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.8], hspace=-0.2, wspace=0.18, left=0.075, right=0.94, bottom=0.085, top=0.965)
    n = T.shape[0]
    ext = [-n * PIX / 120.0, n * PIX / 120.0] * 2
    maps = []
    a = fig.add_subplot(gs[0, 0])
    im = a.imshow(np.log10(T), origin="lower", extent=ext, cmap=SD_CMAP, vmin=20.99, vmax=23.0)
    maps.append(a)
    a.set_title(r"True $\log_{10}\mathcal{D}$ (cm$^{-2}$)", fontsize=7.0, pad=5)
    cb = colorbar_beside(fig, a, im)
    cb.ax.tick_params(labelsize=7.0, pad=2)
    cb.ax.tick_params(which='major', length=3.5)
    cb.ax.tick_params(which='minor', length=2.3)

    for k, (x, lab) in enumerate(((R, r"Relative errors of $\mathcal{D}_{\rm C}$"), 
                                  (Rc, r"Relative errors of $\mathcal{D}_{\rm C}$ (cloud)"))):
        a = fig.add_subplot(gs[0, k + 1])
        im = a.imshow(x, origin="lower", extent=ext, cmap=SD_CMAP, vmin=-0.25, vmax=0.25)
        maps.append(a)
        a.set_title(lab, fontsize=7.0, pad=4)
        cb = colorbar_beside(fig, a, im, "4%")
        cb.ax.tick_params(labelsize=7.0, pad=2)
        cb.ax.tick_params(which='major', length=3.5)
        cb.ax.tick_params(which='minor', length=2.3)

    for a in maps:
        a.set_xticks([])
        a.set_yticks([])

    cols = [224, 448, 672, 896, 1120, 1344, 1568]
    labs = [r"$3{\times}10^{21}$", r"$1{\times}10^{22}$", r"$3{\times}10^{22}$", r"$5{\times}10^{22}$", r"$1{\times}10^{23}$", 
            r"$3{\times}10^{22}$, $p=2$", r"$3{\times}10^{22}$, $H=80″$"]
##    colors = plt.cm.viridis(np.linspace(0.0, 0.9, len(cols)))
    colors = ["red", "green", "blue", "cyan", "magenta", "brown", "orange"]
    off = np.arange(-100, 101)

    for k, (rows, title) in enumerate(((np.r_[30:300, n - 300:n - 30], "Weak cloud area (outer rows)"),
                                       (np.r_[600:1200], "Dense cloud area (central rows)"))):
        b = fig.add_subplot(gs[1, k])

        for c, l, col in zip(cols, labs, colors):
            b.plot(off * PIX, np.nanmedian(Rm[rows][:, c + off], axis=0), "-", color=col, lw=1.0, label=l)

        b.plot(off * PIX, np.nanmedian(np.array([np.median(Rc[rows][:, c + off], axis=0) for c in cols]), axis=0), "--",
               color="black", lw=1.0, label="cloud alone")
               
        b.axhline(0.0, **REF)
        b.set_xlim(-300, 300)
        b.set_ylim(-0.199, 0.199)
        b.set_title(title, fontsize=7.0, pad=4)
        b.tick_params(labeltop=False, labelright=False, labelsize=7.0)
        b.set_xlabel("Offset from crest (″)", labelpad=1.5, fontsize=7.0)
        if k == 0: b.set_ylabel("Relative error of corrected $\\Sigma$", labelpad=1, fontsize=7.0)
        if k == 1: b.tick_params(labelright=True, labelsize=7.0)
        b.set_box_aspect(1.0 / 1.25)

    b.legend(loc="upper left", bbox_to_anchor=(1.17, 1.06), **LEG)

#    fig.tight_layout(pad=0.0, w_pad=0.7)
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_cube(slice_file="/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261008_RT_cirrus3d_12_enh_8e10/cirrus3d_slice_row224.npz", out="fig_cube.pdf"):
  
    """A vertical slice through the RADMC-3D diffuse cloud, compared with the analytic slab.

    The slice contains one row of the field (horizontal) and the line of sight (vertical). Upper panels: the H2 number density, the
    dust temperature of the radiative transfer model, and its difference from the analytic slab with the same surface density, in
    which each layer is shielded by the angular average of the column toward both faces and has the temperature of Eq. (7). Lower
    panel: the temperature against the depth for three columns of different surface density, radiative transfer (solid) and
    analytic slab (dashed).
    """
  
    import make_simsky5 as S5
    import hicorrect2 as h2

    S5.TFLOOR, S5.TSURF, S5.DSTAR, S5.GAMMA = h2.TFLOOR, h2.TSURF, h2.SSTAR, h2.GAMMA
    d = np.load(slice_file)
    n, T, sig, cell = d["n"].astype(float), d["T"].astype(float), d["sigma"].astype(float), float(d["cell_pc"])
    nz, nx = n.shape
    col = n * cell * 3.0856775814913673e18
    zz = (np.cumsum(col, axis=0) - 0.5 * col) / sig[None, :] * 2.0 - 1.0          # column position of each cell, -1 to +1
    Ta = S5.T_of_shielding(S5.slab_seff(zz, S5.KAPHEAT * sig[None, :] / 2.0) * sig[None, :] / 2.0)
    ext = [0.0, nx * cell, 0.0, nz * cell]

    # all four panels have the depth along the line of sight on the horizontal axis

    fig = plt.figure(figsize=(7.1, 4.6))
    gs = fig.add_gridspec(2, 5, width_ratios=[1.0, 1.0, 1.0, 0.2, 2.6], height_ratios=[0.05, 2.0], hspace=0.02, wspace=0.1,
                          left=0.065, right=0.99, bottom=0.09, top=0.90)
    panels = ((np.log10(n), r"$\log_{10} n_{\rm H_2}$ (cm$^{-3}$)", None, None),
              (T, r"$T_{\rm RT}$ (K)", None, None),
              (T - Ta, r"$T_{\rm RT} - T_{\rm slab}$ (K)", -3.0, 3.0))
    extr = [0.0, nz * cell, 0.0, nx * cell]

    for k, (img, lab, vmin, vmax) in enumerate(panels):
        a = fig.add_subplot(gs[1, k])
        im = a.imshow(img.T, origin="lower", extent=extr, aspect="auto", cmap=SD_CMAP, vmin=vmin, vmax=vmax)
        a.tick_params(labeltop=False, labelright=False, labelleft=(k == 0), labelsize=7.5)
        a.set_xticks([0.0, 1.0, 2.0])
        a.set_xlabel("Depth (pc)", labelpad=1.5, fontsize=8.0)
        if k == 0:
            a.set_ylabel("Position along the row (pc)", labelpad=1.5, fontsize=8.0)
        cax = fig.add_subplot(gs[0, k])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal")
        cax.xaxis.set_ticks_position("top")
        cax.xaxis.set_label_position("top")
        cb.ax.tick_params(labelsize=7.3, pad=1)
        cb.set_label(lab, fontsize=8.0, labelpad=3)

    b = fig.add_subplot(gs[1, 4])
    order = np.argsort(sig)
    picks = [order[int(q * (nx - 1))] for q in (0.1, 0.5, 0.97)]
    depth = (np.arange(nz) + 0.5) * cell

    for j, col_i in zip((0, 1, 2), picks):
        c = CCURVE[j]
        e = int(np.floor(np.log10(sig[col_i])))
        b.plot(depth, T[:, col_i], "-", color=c, lw=1.0, label=r"$%.1f{\times}10^{%d}$ cm$^{-2}$, RT" % (sig[col_i] / 10 ** e, e))
        b.plot(depth, Ta[:, col_i], "--", color=c, lw=1.0, label="the same, slab")

    b.tick_params(labeltop=True, labelleft=True, labelright=True, labelsize=7.5)
##    b.yaxis.set_label_position("right")
    b.set_xlim(0.0, nz * cell + 0.004)
    b.set_ylim(np.nanmin(T[:, picks]) - 0.75, np.nanmax(T[:, picks]) + 0.3)
    b.set_xlabel("Depth (pc)", labelpad=1.5, fontsize=8.0)
    b.set_ylabel(r"Temperature $T$ (K)", labelpad=1.5, fontsize=8.0)
    b.legend(loc="upper center", ncol=1, fontsize=7.3, frameon=False, handlelength=1.6, labelspacing=0.3)

    pos = b.get_position()
    ratio = (pos.height * fig.get_figheight()) / (pos.width * fig.get_figwidth())   # current height / width
    b.set_box_aspect(0.8)
    b.set_anchor("S")        # keep the bottom edge aligned with the three image panels
    
    fig.savefig(out)
    plt.close(fig)
    print("written %s" % out)


def fig_cube3d(slice_file="/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261008_RT_cirrus3d_12_enh_8e10/cirrus3d_slice_row224.npz", 
               sigma_file="/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261008_RT_cirrus3d_12_enh_8e10/cirrus3d.sigma.rt.fits", 
                 col_file="/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261008_RT_cirrus3d_12_enh_8e10/cirrus3d_slice_col447.npz", 
                 out="fig_cube3d.pdf"):

    """A cut-away view of the RADMC-3D cloud, lying on its long side and seen from the side of the observer.

    The front face is the map of the surface density of the lower half of the field, oriented as in the maps of fig_simsky6: the
    columns of the field run to the right, the rows upward, and the line of sight points toward the reader. The depth runs into the
    page. The upper face is the cut through the middle row of the field, and the right face the cut along its last column; both
    show the dust temperature. The color bars lie horizontally below the brick.
    """

    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from mpl_toolkits.mplot3d import proj3d
    import matplotlib as mpl

    d = np.load(slice_file)
    T, cell, row = d["T"].astype(float), float(d["cell_pc"]), int(d["row"])
    sig = g(sigma_file)
    nz, nx = T.shape
    tlim, slim = (7.5, 21.0), (21.0, 22.55)
    norm = lambda v, lim: SD_CMAP(np.clip((v - lim[0]) / (lim[1] - lim[0]), 0, 1))

    # plot coordinates

    X1, D, Z1 = nx * cell, nz * cell, row * cell
    fig = plt.figure(figsize=(7.0, 4.2))
    ax = fig.add_axes([-0.04, 0.12, 1.08, 0.90], projection="3d")
    kw = dict(rstride=1, cstride=1, shade=False, linewidth=0, antialiased=False, rasterized=True)
    xx = np.arange(nx + 1) * cell
    yy = np.arange(row + 1) * cell
    dd = np.arange(nz + 1) * cell

    # the front face: the surface density map of rows 0 to row, as in fig_simsky6 (origin at the lower left)

    XX, ZZ = np.meshgrid(xx, yy)
    ax.plot_surface(XX, np.zeros_like(XX), ZZ, facecolors=norm(np.log10(sig[:row, :]), slim), **kw)

    # the upper face: the cut through the middle row, temperature against column and depth

    XX2, DD2 = np.meshgrid(xx, dd)
    ax.plot_surface(XX2, DD2, np.full_like(XX2, Z1), facecolors=norm(T, tlim), **kw)

    # the right face: the cut along the last column, temperature against depth and row

    if os.path.exists(col_file):
        Tc = np.load(col_file)["T"].astype(float)[:, :row]
        DD3, ZZ3 = np.meshgrid(dd, yy)
        ax.plot_surface(np.full_like(DD3, X1), DD3, ZZ3, facecolors=norm(Tc.T, tlim), **kw)

    for (a0, b0, c0), (a1, b1, c1) in [((0, 0, 0), (X1, 0, 0)), ((0, 0, Z1), (X1, 0, Z1)), ((0, 0, 0), (0, 0, Z1)),
                                       ((X1, 0, 0), (X1, 0, Z1)), ((0, 0, Z1), (0, D, Z1)), ((X1, 0, Z1), (X1, D, Z1)),
                                       ((0, D, Z1), (X1, D, Z1)), ((X1, 0, 0), (X1, D, 0)), ((X1, D, 0), (X1, D, Z1))]:
        ax.plot([a0, a1], [b0, b1], [c0, c1], color="0.2", lw=0.6)

    ax.text(0.47 * X1, 0.0, -0.15, "%.1f pc" % X1, fontsize=7.5, ha="center", va="top")
    ax.text(-0.08, 0.0, 0.5 * Z1, "%.1f pc" % Z1, fontsize=7.5, ha="right", va="center")
    ax.text(X1 + 0.1, 0.46 * D, 0.1, "%.0f pc" % D, fontsize=7.5, ha="left", va="top")
    ax.set_box_aspect((X1, D, Z1), zoom=1.12)
    ax.view_init(elev=32, azim=-64)
    ax.set_axis_off()

    # the line of sight points toward the reader: the usual symbol of a vector out of the page, on the front face

    fig.canvas.draw()
    px, py, _ = proj3d.proj_transform(-0.2, 1.7, Z1 - 0.1, ax.get_proj())
    fx, fy = fig.transFigure.inverted().transform(ax.transData.transform((px, py)))
    ax.annotate("", xy=(-0.015, fy - 0.625), xytext=(0.0855, fy - 0.522), xycoords=ax.transAxes, textcoords=ax.transAxes,
                   annotation_clip=False, arrowprops=dict(arrowstyle="-|>", color="black", lw=0.5))
##    fig.text(fx, fy, "$\\odot$", fontsize=13, ha="center", va="center", color="black")
    fig.text(fx + 0.008, fy - 0.52, "to observer", rotation=0, fontsize=7.5, ha="right", va="center")

    bars = [(r"Front face: $\log_{10}\Sigma$ (cm$^{-2}$)", slim), (r"Upper face: $T$ (K)", tlim)]

    if os.path.exists(col_file):
        bars.append((r"Right face: $T$ (K)", tlim))
    w = 0.84 / len(bars)

    for k, (lab, lim) in enumerate(bars):
        cax = fig.add_axes([0.23 + k * w * 0.64, 0.18, w - 0.11, 0.022])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(*lim), cmap=SD_CMAP), cax=cax, orientation="horizontal")
        cax.tick_params(labelsize=7.5, pad=1.5, top=False, labeltop=False)
        cax.tick_params(which='major', length=3.5)
        cax.tick_params(which='minor', length=2.3)
        if k == 0: 
          cb.set_label(lab, fontsize=7.5, labelpad=1)
        else: cb.set_label(lab, fontsize=7.5, labelpad=2)

    fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)
    print("written %s" % out)

if __name__ == "__main__":

#    fig_geometry()
#    fig_shielding()
    fig_method()
#    fig_relation()
#    fig_grid_validation()
#    fig_cylinders()
#    fig_simsky6()
#    fig_aquila()
#    fig_cube3d()
#    fig_cube()
#    fig_superposition()
#    fig_mixing()
#    fig_angular()
