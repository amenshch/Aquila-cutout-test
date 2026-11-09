#!/usr/bin/env python3
"""Radial profiles of the true, uncorrected and corrected surface density around a model core.

The accuracy of the correction has so far been measured at the peak of a core and over its whole footprint, and the two disagree:
the integrated surface density is over-corrected by a third at several nodes while the peak is over-corrected by a tenth. The two
can only differ if the correction is too strong away from the centre, so this script measures where in radius that happens.

Two things are computed on the same annuli. The first is the mean surface density in each annulus, for the true image, for the
uncorrected one and for the corrected one, together with the ratios of the last two to the first and the mean correction factor.
This says at which radii the correction is too strong. The second is the running sum within each radius, for the same three images,
with the ratios of the last two to the first. This says at which radii the error of the integrated surface density accumulates, and
its last value is the number quoted for the footprint as a whole.

The centre is the brightest pixel of the true image unless one is given. The pixel size is taken from the header of the true image
and may be overridden.

usage:

  radial_profiles.py --true model.fits --uncorrected hi.fits --corrected hi.corrected.fits
                     --factor hi.corrfactor.fits --rmax 120 --nbins 24
                     --out radial_profiles.pdf --table radial_profiles.txt

The correction factor is optional; without it the ratio of the corrected image to the uncorrected one is used instead, which is the
same quantity wherever both are defined.
"""
import argparse
import os
import sys

import numpy as np

try:
    from astropy.io import fits
except ImportError:
    sys.exit("ERROR: this script requires astropy (pip install astropy)")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import aa_plotstyle
except ImportError:
    aa_plotstyle = None


def read_image(name):
    """The first two-dimensional image of a FITS file, as double precision, with its header."""
    with fits.open(name) as hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                return np.array(hdu.data, dtype=np.float64).squeeze(), hdu.header
    sys.exit("ERROR: no image found in '%s'" % name)


def pixel_size(header):
    """The pixel size in arcsec, from whichever of the usual keywords the header carries."""
    for key in ("CDELT2", "CDELT1", "CD2_2", "CD1_1"):
        if key in header:
            value = abs(float(header[key]))
            if value > 0.0:
                return value * 3600.0
    return None


def main():
    p = argparse.ArgumentParser(description="radial profiles of the correction around a model core")
    p.add_argument("--true", required=True, help="the true surface density image, at the same resolution as the others")
    p.add_argument("--uncorrected", required=True, help="the hires surface density image given to hicorrect")
    p.add_argument("--corrected", required=True, help="the corrected surface density image")
    p.add_argument("--factor", default=None, help="the correction factor image; the ratio of the two above is used without it")
    p.add_argument("--centre", type=float, nargs=2, default=None, metavar=("X", "Y"),
                   help="the centre in pixels, counted from one; the brightest pixel of the true image by default")
    p.add_argument("--pix", type=float, default=None, help="pixel size in arcsec, overriding the header")
    p.add_argument("--rmax", type=float, default=None, help="largest radius in arcsec; a quarter of the shorter side by default")
    p.add_argument("--nbins", type=int, default=24, help="number of annuli (default: %(default)s)")
    p.add_argument("--out", default="radial_profiles.pdf", help="the figure to write (default: %(default)s)")
    p.add_argument("--table", default=None, help="a text file to write the profiles to")
    p.add_argument("--label", default="", help="a label for the title of the figure")
    a = p.parse_args()

    true, header = read_image(a.true)
    unc, _ = read_image(a.uncorrected)
    cor, _ = read_image(a.corrected)
    for name, image in (("uncorrected", unc), ("corrected", cor)):
        if image.shape != true.shape:
            sys.exit("ERROR: the %s image has shape %s and the true image has shape %s" % (name, image.shape, true.shape))
    if a.factor is not None:
        fac, _ = read_image(a.factor)
        if fac.shape != true.shape:
            sys.exit("ERROR: the factor image has shape %s and the true image has shape %s" % (fac.shape, true.shape))
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            fac = np.where(unc > 0.0, cor / unc, np.nan)

    pix = a.pix if a.pix is not None else pixel_size(header)
    if pix is None:
        sys.exit("ERROR: the pixel size is not in the header of the true image; give it with --pix")

    ny, nx = true.shape
    if a.centre is None:
        jc, ic = np.unravel_index(int(np.nanargmax(np.where(np.isfinite(true), true, -np.inf))), true.shape)
        xc, yc = float(ic), float(jc)
    else:
        xc, yc = a.centre[0] - 1.0, a.centre[1] - 1.0

    y, x = np.mgrid[0:ny, 0:nx]
    r = np.sqrt((x - xc) ** 2 + (y - yc) ** 2) * pix
    rmax = a.rmax if a.rmax is not None else 0.25 * min(nx, ny) * pix
    edges = np.linspace(0.0, rmax, a.nbins + 1)
    centre = 0.5 * (edges[:-1] + edges[1:])

    good = np.isfinite(true) & np.isfinite(unc) & np.isfinite(cor) & (true > 0.0)
    index = np.digitize(r, edges) - 1

    print(" centre at pixel (%.1f, %.1f), pixel size %.3f arcsec, radius out to %.1f arcsec in %d annuli"
          % (xc + 1.0, yc + 1.0, pix, rmax, a.nbins))
    print(" %s" % a.label)
    print()
    print("   radius     true       uncorrected   corrected    unc/true   cor/true   factor    within r: unc/true  cor/true  pixels")

    mt = np.full(a.nbins, np.nan)
    mu = np.full(a.nbins, np.nan)
    mc = np.full(a.nbins, np.nan)
    mf = np.full(a.nbins, np.nan)
    st = np.zeros(a.nbins)
    su = np.zeros(a.nbins)
    sc = np.zeros(a.nbins)
    num = np.zeros(a.nbins, dtype=np.int64)
    for k in range(a.nbins):
        sel = good & (index == k)
        num[k] = sel.sum()
        if num[k] > 0:
            mt[k] = true[sel].mean()
            mu[k] = unc[sel].mean()
            mc[k] = cor[sel].mean()
            ff = fac[sel]
            mf[k] = np.nanmean(ff) if np.isfinite(ff).any() else np.nan
            st[k] = true[sel].sum()
            su[k] = unc[sel].sum()
            sc[k] = cor[sel].sum()
    ct = np.cumsum(st)
    cu = np.cumsum(su)
    cc = np.cumsum(sc)
    with np.errstate(divide="ignore", invalid="ignore"):
        rt_u = np.where(ct > 0.0, cu / ct, np.nan)
        rt_c = np.where(ct > 0.0, cc / ct, np.nan)

    for k in range(a.nbins):
        if num[k] > 0:
            print("   %7.2f %11.4e %11.4e %11.4e %10.4f %10.4f %9.4f %14.4f %9.4f %8d"
                  % (centre[k], mt[k], mu[k], mc[k], mu[k] / mt[k], mc[k] / mt[k], mf[k], rt_u[k], rt_c[k], num[k]))
    print()
    print(" within the largest radius: uncorrected/true %.4f, corrected/true %.4f" % (rt_u[-1], rt_c[-1]))

    if a.table is not None:
        with open(a.table, "a") as fh:
            fh.write("# %s   centre (%.1f, %.1f), pixel %.3f arcsec\n" % (a.label, xc + 1.0, yc + 1.0, pix))
            fh.write("# radius, true, uncorrected, corrected, unc/true, cor/true, factor, cumulative unc/true, cor/true, pixels\n")
            for k in range(a.nbins):
                if num[k] > 0:
                    fh.write("%9.3f %13.6e %13.6e %13.6e %10.5f %10.5f %9.5f %10.5f %10.5f %9d\n"
                             % (centre[k], mt[k], mu[k], mc[k], mu[k] / mt[k], mc[k] / mt[k], mf[k], rt_u[k], rt_c[k], num[k]))
            fh.write("\n")

    if aa_plotstyle is not None:
        aa_plotstyle.apply()
    fig, ax = plt.subplots(2, 1, figsize=(3.5, 4.6), sharex=True)
    ax[0].plot(centre, mt, color="0.2", linewidth=1.4, label="true")
    ax[0].plot(centre, mu, color="#1f4e9c", linewidth=1.2, linestyle=(0, (4, 2)), label="uncorrected")
    ax[0].plot(centre, mc, color="#c1272d", linewidth=1.2, label="corrected")
    ax[0].set_yscale("log")
    ax[0].set_ylabel(r"$N_{\rm H_2}$ (cm$^{-2}$)")
    ax[0].legend(loc="upper right", fontsize=7.5)
    ax[1].axhline(1.0, color="0.5", linewidth=0.7, linestyle=(0, (4, 3)))
    ax[1].plot(centre, mu / mt, color="#1f4e9c", linewidth=1.2, linestyle=(0, (4, 2)), label="uncorrected / true")
    ax[1].plot(centre, mc / mt, color="#c1272d", linewidth=1.2, label="corrected / true")
    ax[1].plot(centre, rt_c, color="#1a9641", linewidth=1.2, linestyle=(0, (1, 2)), label="corrected / true within r")
    ax[1].set_xlabel("radius (arcsec)")
    ax[1].set_ylabel("ratio")
    ax[1].legend(loc="best", fontsize=7.5)
    if a.label:
        ax[0].set_title(a.label, fontsize=8)
    fig.tight_layout()
    fig.savefig(a.out)
    print(" figure written to '%s'" % a.out)
    if a.table is not None:
        print(" table written to '%s'" % a.table)


if __name__ == "__main__":
    main()
