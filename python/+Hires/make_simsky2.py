#!/usr/bin/env python3
"""
make_simsky.py -- a simulated sky in which the dust temperature follows the
total shielding column along each line of sight.

Each line of sight carries the total surface density D = D_cirrus + D_filament.
Illuminated from both sides, the material at fractional depth u in [0,1] is
shielded by the column s = D min(u, 1-u), and its temperature follows one
shielding law,

    T(s) = T_0 (1 + s / D_star)^(-gamma).

The intensity is the integral of the emission over depth,

    I_lambda = D  integral_0^1  m_lambda(T(s(u))) du ,
    m_lambda(T) = B_lambda(T) kappa_lambda eta mu_H2 m_H ,
    kappa_lambda = 10 (nu / 1e12 Hz)^2 cm^2 g^-1 .

Every line of sight therefore has a smooth temperature gradient rather than one
or two discrete values, the surface of the cloud is at T_0 everywhere, and a
dense narrow structure is cold in its interior because its own column shields
it.  No component is assigned a temperature of its own, so the construction is
self-consistent: two lines of sight with the same total column have the same
temperature profile whatever material produced that column.

The parameters are chosen to span the range Herschel observes: T_0 = 20 K at
the cloud surface, about 17 K at the mid-plane of the diffuse cirrus and about
10 K at the mid-plane of the filament crest.
"""
import sys

import numpy as np
from astropy.io import fits
from scipy import ndimage

H_PLANCK, K_BOLTZ, C_LIGHT = 6.62607015e-27, 1.380649e-16, 2.99792458e10
AMU, MU_H2, ETA = 1.6605402e-24, 2.8, 1.0e-2
KAPPA0, NU0, BETA = 10.0, 1.0e12, 2.0
PIX = 3.0
FWHM_TO_SIGMA = 1.0 / 2.354820045

def T_of_shielding(s):
    """Dust temperature behind a shielding column s, in cm-2.

    The form is the one used in the getsf paper, applied here to the column
    shielding each layer rather than to the total column of the pixel, so that
    a line of sight carries a gradient instead of a single value.  It gives
    20.0 K at the surface, 16.7 K behind 1e21 cm-2, and saturates at 10 K deep
    inside a cloud.
    """
    return 200.0 / (1.0e-20 * np.asarray(s, float) + 20.0) + 10.0

T0 = 20.0            # dust temperature at the surface, K (value of the law at s = 0)
NLAYER = 64          # layers used for the integral over depth
VERSION = "2026-08-27: WAVE and BEAM written, band 0 for surface density"

BANDS = [(70., 8.4, "070", "r8p4"), (100., 9.4, "100", "r9p4"),
         (160., 13.5, "160", "r13p5"), (250., 18.2, "250", "r18p2"),
         (350., 24.9, "350", "r24p9"), (500., 36.3, "500", "r36p3")]
WIDTH = 40           # half-maximum width of the filaments, arcsec
STEM = "simsky2.w%d.%s.%s.fits"


def m_lambda(t, wave_um):
    nu = C_LIGHT / (wave_um * 1.0e-4)
    kap = KAPPA0 * (nu / NU0) ** BETA
    x = np.clip(H_PLANCK * nu / (K_BOLTZ * np.asarray(t, float)), 1e-8, 700.0)
    b = (2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2) / np.expm1(x)
    return b * kap * ETA * MU_H2 * AMU / 1.0e-17          # MJy/sr per cm-2


def _double_card(key, value, comment):
    """A header card holding a double in fixed exponential form.

    Astropy writes a float in the shortest form that round-trips, so a
    wavelength of 160 becomes '160.0'.  That is valid FITS and is read
    correctly, but it does not look like the cards written by cfitsio through
    ftpkyd, which are of the form ' 1.6000000000000E+02'.  Writing the same
    form here removes any doubt about how the value will be read.
    """
    return fits.Card.fromstring(
        ("%-8s= %20.13E / %s" % (key, float(value), comment))[:80].ljust(80))


def write_fits(name, data, header, doubles=()):
    """Write an image, then set the double-valued keywords in a second pass.

    Astropy drops BZERO and BSCALE when it writes floating-point data, even
    when the header supplied carries them, so they have to be put back by
    reopening the file in update mode.  The keywords listed in 'doubles', each
    a triple of keyword, value and comment, are written there as well, in
    fixed exponential form.
    """
    fits.PrimaryHDU(np.asarray(data, np.float32), header).writeto(
        name, overwrite=True)
    with fits.open(name, mode="update") as hdul:
        h = hdul[0].header
        h["BZERO"] = 0.0
        h["BSCALE"] = 1.0
        for key, value, comment in doubles:
            if key in h:
                del h[key]
            h.append(_double_card(key, value, comment))


def repair(a, name=""):
    """Replace non-positive pixels by the nearest positive one.

    The cirrus map carries a border of zeros -- its first and last columns and
    its last three rows -- which a convolution would smear inward as a dark rim.
    A surface density of zero is not a measurement, so such pixels are filled
    from their nearest valid neighbour before anything else is done.
    """
    bad = ~(a > 0) | ~np.isfinite(a)
    if not bad.any():
        return a
    idx = ndimage.distance_transform_edt(bad, return_distances=False,
                                         return_indices=True)
    print("   %s: %d non-positive pixels (%.4f per cent) filled from their "
          "nearest valid neighbour" % (name, bad.sum(), 100.0 * bad.mean()))
    return a[tuple(idx)]


def smooth(a, fwhm_as):
    """Convolve with a Gaussian, padding the image first and cropping after.

    The image is expanded on every side by twice the half-maximum width of the
    beam, by replicating the edge, so that the convolution never runs off the
    array, and the result is cut back to the original size.
    """
    if fwhm_as <= 0:
        return a.copy()
    n = int(np.ceil(2.0 * fwhm_as / PIX))
    b = np.pad(a, n, mode="edge")
    b = ndimage.gaussian_filter(b, fwhm_as * FWHM_TO_SIGMA / PIX,
                                mode="nearest", truncate=4.0)
    return b[n:-n, n:-n]


def main():
    global WIDTH
    if len(sys.argv) > 1:
        WIDTH = int(sys.argv[1])
    get = lambda f: np.array(fits.getdata(f), float).squeeze()
    hdr = fits.getheader("bdm1500x1p5.fits")
    print("repairing the input maps:")
    D = (repair(get("bdm1500x1p5.fits"), "cirrus")
         + repair(get("fils.sha%dx0p75-.x3e22.fits" % WIDTH), "filament"))
    print("total surface density: min %.3g, median %.3g, max %.3g cm-2"
          % (D.min(), np.median(D), D.max()))

    u = (np.arange(NLAYER) + 0.5) / NLAYER
    depth = np.minimum(u, 1.0 - u)                        # shielded from both sides
    tgrid = np.geomspace(9.0, T0 * 1.001, 4000)

    # mass-weighted mean temperature of each line of sight, for reference only
    tmean = np.zeros_like(D)
    for dd in depth:
        tmean += T_of_shielding(dd * D)
    tmean /= NLAYER

    for wave, beam, tag, rtag in BANDS:
        table = m_lambda(tgrid, wave)
        acc = np.zeros_like(D)
        for dd in depth:
            t = T_of_shielding(dd * D)
            acc += np.interp(t, tgrid, table)
        img = smooth(D * acc / NLAYER, beam)
        h = hdr.copy()
        h["BUNIT"] = "MJy/sr"
        h["BZERO"], h["BSCALE"] = 0.0, 1.0
        h["COMMENT"] = "simulated sky, temperature from total shielding column"
        h["TLAW"] = ("200/(1e-20*s+20)+10 K", "s = shielding column, cm-2")
        name = STEM % (WIDTH, tag, rtag)
        # hires reads the wavelength of an image from the WAVE keyword
        write_fits(name, img, h,
                   [("WAVE", wave, "wavelength (microns)"),
                    ("BEAM", beam, "half-maximum beam width (arcsec)")])
        print("%-58s median %10.3f MJy/sr" % (name, np.median(img)))

    for nm, arr, unit, bm in (
            (("simsky2.w%d.surfdens.pix.fits" % WIDTH), D, "cm-2", 0.0),
            (("simsky2.w%d.surfdens.r13p5.fits" % WIDTH), smooth(D, 13.5), "cm-2", 13.5),
            (("simsky2.w%d.surfdens.r36p3.fits" % WIDTH), smooth(D, 36.3), "cm-2", 36.3),
            (("simsky2.w%d.tmean.pix.fits" % WIDTH), tmean, "K", 0.0)):
        h = hdr.copy()
        h["BUNIT"], h["BZERO"], h["BSCALE"] = unit, 0.0, 1.0
        # a surface density or temperature image is band 0 micron in the getsf
        # convention, so the wavelength is set to zero rather than inherited
        dbl = [("WAVE", 0.0, "wavelength (microns); 0 for surface density")]
        if bm > 0:
            dbl.append(("BEAM", bm, "half-maximum beam width (arcsec)"))
        elif "BEAM" in h:
            del h["BEAM"]
        write_fits(nm, arr, h, dbl)
        print("%-58s median %10.4g %s" % (nm, np.median(arr), unit))

    with open("+hires.cfg.simsky2.w%d" % WIDTH, "w") as f:
        f.write("#\n# HIRES: CONFIGURATION PARAMETERS FOR HIRES\n#\n"
                "# PM  WAVE  BEAM  OFFSET   IMAGE   [MASKIMAGE]\n#\n")
        for wave, beam, tag, rtag in BANDS:
            f.write("  %s   %s  %4.1f   0.0   %s\n"
                    % ("-" if wave < 160 else "+", tag, beam,
                       STEM % (WIDTH, tag, rtag)))
    print("written +hires.cfg.simsky2.w%d" % WIDTH)

    print()
    print("make_simsky.py version %s" % VERSION)
    print("header keywords written:")
    print("%42s %10s %10s" % ("file", "WAVE", "BEAM"))
    import glob
    for f in sorted(glob.glob("simsky2.w%d.*.fits" % WIDTH)):
        hh = fits.getheader(f)
        print("%42s %10s %10s"
              % (f, hh.get("WAVE", "absent"), hh.get("BEAM", "absent")))

    print()
    print("temperature at the mid-plane, as a function of the total column:")
    print("%16s %14s %14s"
          % ("N_H2 (cm-2)", "mid-plane (K)", "mass mean (K)"))
    for q in (5, 25, 50, 75, 95, 99.9):
        n = np.percentile(D, q)
        print("%16.3g %14.2f %14.2f"
              % (n, T_of_shielding(0.5 * n),
                 np.median(tmean[(D > 0.98 * n) & (D < 1.02 * n)])))


if __name__ == "__main__":
    main()
