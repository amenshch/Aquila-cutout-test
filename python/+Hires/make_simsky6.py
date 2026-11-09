#!/usr/bin/env python3
"""make_simsky6.py -- a simulated sky with radiative transfer filaments on the analytic cirrus.

The sky of make_simsky5.py has the same cirrus as this one, but its filaments were analytic: each parcel of a filament was
shielded by the column along the radius, from the parcel outward. That makes the flanks of a filament too warm. Here the
filaments are replaced by the RADMC-3D models of cylinders, which have the right temperatures.

The sky is built in four steps.

  1. The cirrus. Its intensities are computed as in make_simsky5.py: a plane-parallel slab, every layer shielded by the angular
     average of the column toward both faces, with the dust temperature from the final relation T(Sigma_a) of the paper.

     With RTCIRRUS = True, these intensities are multiplied, band by band, by the ratio of the emission of the RADMC-3D cirrus
     of make_cirrus3d.py to that of the same analytic slab, computed column by column on the 3D grid and interpolated to the
     pixels: the sky then has the temperatures of the radiative transfer cirrus and the fine structure of the cirrus map.

  2. The filaments. Each filament comes from one cylinder model, computed inside a uniform envelope that stands for the cirrus
     around it (cylembsd and cylrenv of the cylinder module). Two image runs are made of each model: the full one, and the
     background run with the argument doback, in which the volume of the filament is emptied while the envelope and all
     temperatures are kept. The emission of the filament in the cirrus is their difference, exactly as for the injected cores.
     In the images the cylinder axis is vertical. Near the ends of the cylinder the temperature differs, so only the middle part
     is used: the median of the middle rows gives one profile across the filament, for each waveband and for the true surface
     density. The profile is then repeated in every row of the sky, with its axis at the position of one filament of
     make_simsky5.py. The filament therefore has no ends.

  3. The sum. In every waveband, the intensities of the cirrus and of the filaments are added, and the sum is convolved to the
     beam of that waveband, exactly as in make_simsky5.py. The true surface density is the sum of the true surface densities.

  4. The cores. They are not added here. Run inject_cores_simsky.py on the output, with the stem simsky6.w40, as for simsky4.

The envelope gives each filament the shielding of the cirrus around it. Its surface density through the axis should equal that of
the cirrus beneath the filament; python3 make_simsky6.py --cirrus prints it for the seven positions.

The script needs, in the working directory: make_simsky5.py (for its functions), the cirrus map bdm1500x1p5.fits, the filament
map fils.sha40x0p75-.x3e22.fits of make_simsky5.py (only to find the positions of the seven filaments), and one directory per
cylinder model with the images that the cylinder module writes: the unconvolved intensities at 70, 100, 160, 250, 350 and 500
microns (mc.070um.fits, mc.100um.fits, ..., mc.500um.fits), in MJy/sr, and the true surface density mc.surfdens.fits, all on
pixels of 3 arcsec.

Usage:

    python3 make_simsky6.py [<cylinder directory> ...]
    python3 make_simsky6.py --cirrus           (prints the cirrus surface density along the seven filament axes)
    python3 make_simsky6.py --envelopes        (writes cylenv.inp, the cirrus along each axis, into each model directory)
    python3 make_simsky6.py --cirrus-only      (the sky with the cirrus alone, to test the correction on it separately)

With ENVZ = True, each model is computed with an envelope that follows the cirrus along its filament (cylenvz = 1 and the table
cylenv.inp), and the full images are placed, shifted only sideways, instead of one repeated middle profile. The images must then
have the size of the sky, 1793 x 1793 pixels of 3 arcsec, and the cylinder must be longer than the sky (for instance 10.3 pc), so
that its ends, disturbed over about the radius of the envelope, fall outside the images. The orientation of every image along the
filament is checked against the cirrus, using the envelope in its background image, and flipped if needed.

With no directories, the seven models are taken in this order, from left to right:

    cyl_sd3e21 cyl_sd1e22 cylinder cyl_sd5e22 cyl_sd1e23 cyl_p2 cyl_w80

which places them in the order of increasing surface density on the axis, followed by the steeper and the wider profile.
"""

import os
import sys

import numpy as np
from astropy.io import fits

import cirrus_enhance as E
import hicorrect2 as H
import make_simsky5 as S5

# the final relation T(Sigma_a) of the paper, fitted to 2086 spheres and 7 cylinders
S5.TFLOOR, S5.TSURF, S5.DSTAR, S5.GAMMA = H.TFLOOR, H.TSURF, H.SSTAR, H.GAMMA     # the relation of hicorrect2.py

PREFIX = "simsky6"
WIDTH = 40
CIRRUS = "bd.fits"                 # the unmodified cirrus, smoother along the filaments than bdm1500x1p5.fits
CIRRUS_SCALE = 0.1235              # factor applied to the cirrus, matching the median of bdm1500x1p5.fits (--cirrus)
CIRRUS_REF = "bdm1500x1p5.fits"    # the cirrus of simsky4 and simsky5, used only to suggest CIRRUS_SCALE
FILMAP = "fils.sha%dx0p75-.x3e22.fits" % WIDTH
MODELS = ["/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_RT_cyl_sd3e21_enh", 
          "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_RT_cyl_sd1e22_enh", 
          "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_RT_cylinder_enh", 
          "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_RT_cyl_sd5e22_enh", 
          "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_RT_cyl_sd1e23_enh", 
          "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_RT_cyl_p2_enh", 
          "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_RT_cyl_w80_enh"]
RTNAME = "mc.%03gum.fits"          # the intensity image of a cylinder model: mc.070um.fits, mc.100um.fits, ...
RTSURF = "mc.surfdens.fits"        # its true surface density
RTBACK = ".bgl.fits"               # the suffix of the background images of the doback run, replacing ".fits"
ENVZ = True                        # envelopes that follow the cirrus along each filament (cylenvz = 1): full images are placed
DIST = 260.0                       # distance (pc), to convert the rows of the sky into heights along the cylinders
ENVSMOOTH = 5                      # rows over which the cirrus is smoothed when its profile along a filament is tabulated
RTCIRRUS = True                    # scale the analytic cirrus to the emission of the RADMC-3D cirrus (compare_cirrus3d.py)
RTCIRRUS_DIR = "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/261007_RT_cirrus3d_12_enh"       # the directory of that cirrus, with cirrus3d.json, cirrus3d.rt.*.fits, cirrus3d.sigma.rt.fits
RTEDGE = 2                         # cells at the edges of the RADMC-3D field replaced by their inner neighbors
CRATER = False                     # remove from the cirrus the material that each filament displaces
CIRRUS_DEPTH = 2.0                 # depth of the cirrus along the line of sight (pc), which sets its volume density
ROUT = 0.5                         # outer radius of the filaments (pc), as cylrout of the models
ENVMATCH = "shielding"             # "shielding": the envelope shields the axis of the filament as the cirrus would;
                                   # "column": its surface density through the axis equals that of the cirrus
MIDDLE = 0.45                      # the middle rows used: this fraction of the image height, centered
STEM = PREFIX + ".w%d.%s.%s.fits"


def getimg(path):
    with fits.open(path) as h:
        return np.nan_to_num(np.array(h[0].data, float)).squeeze(), h[0].header


def pixel_arcsec(header):
    for key in ("CDELT2", "CDELT1"):
        if key in header:
            return abs(float(header[key])) * 3600.0
    return S5.PIX


def filament_columns(nfil=7):
    """The columns of the axes of the seven filaments of make_simsky5.py, from left to right."""
    f, _ = getimg(FILMAP)
    prof = np.median(f, axis=0)
    cols = []
    for i in range(1, len(prof) - 1):
        lo, hi = max(0, i - 30), min(len(prof), i + 31)
        if prof[i] == prof[lo:hi].max() and prof[i] > 0.5 * prof.max():
            cols.append(i)
    if len(cols) != nfil:
        raise ValueError("found %d filaments in %s instead of %d" % (len(cols), FILMAP, nfil))
    return cols


def filament_image(model, name):
    """The image of the filament alone: the full image minus the crater background of the doback run, if there is one."""
    img, h = getimg(os.path.join(model, name))
    back = os.path.join(model, name.replace(".fits", RTBACK))
    if os.path.exists(back):
        bg, _ = getimg(back)
        img = img - bg
    else:
        print(" %s: no %s, the filament is used without background subtraction" % (model, os.path.basename(back)))
    return img, h


def shell_equivalent(sig_c):
    """The surface density through the axis of a cylindrical envelope that shields the axis of the filament as a plane-parallel
    cirrus of surface density sig_c shields its mid-plane.

    For the cirrus, a direction at the angle arccos(mu) to its normal crosses the column sig_c / (2 |mu|). For the envelope, a
    direction at the angle arccos(mu) to the axis of the cylinder crosses sig_e / (2 sqrt(1 - mu^2)). Directions in the plane of
    the cirrus are shielded strongly by the cirrus but weakly by the envelope, so an envelope with the same surface density through
    the axis lets in about twice as much radiation. The envelope is therefore made denser, so that the mean intensities of
    Eq. (5), averaged over isotropic directions with the heating opacity, are equal at the axis.
    """
    from scipy.optimize import brentq
    mu = (np.arange(4000) + 0.5) / 4000.0
    out = np.empty_like(sig_c)
    for i, sc in enumerate(sig_c):
        if sc <= 0.0:
            out[i] = 0.0
            continue
        js = np.mean(np.exp(-S5.KAPHEAT * sc / 2.0 / mu))
        out[i] = brentq(lambda se: np.mean(np.exp(-S5.KAPHEAT * se / 2.0 / np.sqrt(1.0 - mu ** 2))) - js, 0.5 * sc, 20.0 * sc)
    return out


def write_envelopes(models, cols, Dc):
    """Write cylenv.inp into every model directory: the cirrus surface density along the axis of its filament, against z."""
    from scipy.ndimage import uniform_filter1d
    ny = Dc.shape[0]
    z = (np.arange(ny) - 0.5 * (ny - 1)) * S5.PIX * DIST / 206264.806
    for model, col in zip(models, cols):
        prof = uniform_filter1d(np.median(Dc[:, col - 2:col + 3], axis=1), ENVSMOOTH, mode="nearest")
        if ENVMATCH == "shielding":
            levels = np.geomspace(max(prof.min(), 1e19), prof.max(), 60)
            prof = np.interp(prof, levels, shell_equivalent(levels))
        os.makedirs(model, exist_ok=True)
        with open(os.path.join(model, "cylenv.inp"), "w") as f:
            f.write("! surface density of the envelope through the axis, from the cirrus %s x %g along column %d of the sky,\n"
                    % (CIRRUS, CIRRUS_SCALE, col))
            f.write("! matched by %s\n" % ("the shielding of the axis" if ENVMATCH == "shielding" else "the surface density"))
            f.write("! z (pc), measured from the middle of the cylinder, and surface density (cm^-2); row 0 of the sky is at"
                    " the smallest z\n")
            for a, b in zip(z, prof):
                f.write("%10.5f %12.5e\n" % (a, b))
        print("   %-12s cylenv.inp written; set cylembsd = %.3e (the median) and cylenvz = 1" % (model, np.median(prof)))


def oriented(img, model, col_axis, Dc, col):
    """The image itself, or flipped upside down if its envelope runs the other way along the filament than the cirrus."""
    back = os.path.join(model, RTSURF.replace(".fits", RTBACK))
    if not os.path.exists(back):
        return img
    bg, _ = getimg(back)
    off = int(0.75 * 1.0 * 206264.806 / DIST / S5.PIX)          # three quarters of a parsec from the axis, inside the envelope
    env = bg[:, min(col_axis + off, bg.shape[1] - 1)]
    cir = Dc[:, col]
    if np.corrcoef(env, cir)[0, 1] < np.corrcoef(env[::-1], cir)[0, 1]:
        print("   %-12s image flipped upside down to match the cirrus" % model)
        return img[::-1]
    return img


def middle_profile(img):
    """The median over the middle rows of a cylinder image, and the column of its axis."""
    ny = img.shape[0]
    half = int(0.5 * MIDDLE * ny)
    prof = np.median(img[ny // 2 - half:ny // 2 + half + 1, :], axis=0)
    return prof


def shift_columns(img, axis_col, target_col):
    """The image shifted sideways so that the axis of the filament is at the target column; rows are kept as they are."""
    out = np.zeros_like(img)
    d = target_col - axis_col
    if d >= 0:
        out[:, d:] = img[:, :img.shape[1] - d]
    else:
        out[:, :d] = img[:, -d:]
    return out


def place(prof, axis_col, target_col, shape):
    """A sky image of the given shape whose every row is the profile, shifted so that its axis is at the target column."""
    out = np.zeros(shape)
    nx = shape[1]
    j = np.arange(nx)
    src = j - target_col + axis_col
    ok = (src >= 0) & (src < len(prof))
    row = np.zeros(nx)
    row[ok] = prof[src[ok]]
    out[:, :] = row[None, :]
    return out


def crater_fraction(shape, cols):
    """The fraction of the cirrus column that the filaments displace, pixel by pixel.

    The cirrus is a slab of depth CIRRUS_DEPTH along the line of sight, so its volume density is its surface density divided by
    that depth. A filament of radius ROUT lies in its mid-plane and displaces the cirrus within its volume. A line of sight at the
    distance b from the axis crosses the filament over the chord 2 sqrt(ROUT^2 - b^2), so the fraction of the cirrus column that is
    removed is that chord divided by the depth, at most one.
    """
    nx = shape[1]
    b = np.abs(np.arange(nx)[None, :] - np.array(cols)[:, None]) * S5.PIX * DIST / 206264.806     # pc
    chord = 2.0 * np.sqrt(np.clip(ROUT ** 2 - b ** 2, 0.0, None))
    frac = np.clip(chord.max(axis=0) / CIRRUS_DEPTH, 0.0, 1.0)
    return np.repeat(frac[None, :], shape[0], axis=0)


def rt_cirrus_ratio(wave, shape):
    """The emission of the RADMC-3D cirrus divided by that of the analytic slab with the same surface density, in one band,
    interpolated from the cells of the 3D grid to the pixels of the sky.

    Both emissions are computed column by column in the same units (hicorrect2.cnu times the Planck function), the RADMC-3D one by
    compare_cirrus3d.py from the Monte Carlo temperatures, the analytic one here from the relation of hicorrect2.py. Their ratio
    varies smoothly across the field, because the temperatures do, and it multiplies the analytic emission of the sky at 3 arcsec,
    which keeps the fine structure of the cirrus map. The cells nearest the edges of the field are replaced by their inner
    neighbors, because the first version of make_cirrus3d.py did not fill the border of zeros of bd.fits.
    """
    import json
    from scipy.ndimage import map_coordinates
    par = json.load(open(os.path.join(RTCIRRUS_DIR, "cirrus3d.json")))
    fac = int(round(par["cell"] / par["pix"]))
    rt, _ = getimg(os.path.join(RTCIRRUS_DIR, "cirrus3d.rt.%03d.fits" % int(wave)))
    sig, _ = getimg(os.path.join(RTCIRRUS_DIR, "cirrus3d.sigma.rt.fits"))
    zl = (np.arange(S5.NLAYER) + 0.5) / S5.NLAYER * 2.0 - 1.0
    t0c = S5.KAPHEAT * sig / 2.0
    an = np.zeros_like(sig)
    for zz in zl:
        T = S5.T_of_shielding(S5.slab_seff(np.full_like(sig, zz), t0c) * sig / 2.0)
        an += sig / S5.NLAYER * H.cnu(wave) * H.planck(T, wave)
    ratio = rt / an
    e = RTEDGE
    inner = ratio[e:-e, e:-e]
    ratio = np.pad(inner, e, mode="edge")
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]].astype(float)
    coords = [(yy - 0.5 * (fac - 1)) / fac, (xx - 0.5 * (fac - 1)) / fac]
    return map_coordinates(ratio, coords, order=1, mode="nearest")


def cirrus_emission(Dc, wave, frac=None):
    """Intensity of the cirrus at one wavelength, in MJy/sr, as in make_simsky5.py with ANGULAR = True.

    With frac, the layers of the slab nearest its mid-plane are removed, the fraction frac of the column in each pixel: the
    material displaced by a filament. The temperatures of the remaining layers are those of the full slab; the filament that
    replaces the removed layers shields them at least as much."""
    tgrid = np.geomspace(0.999 * S5.TFLOOR, S5.TSURF * 1.001, 4000)
    table = S5.m_lambda(tgrid, wave)
    zlayer = (np.arange(S5.NLAYER) + 0.5) / S5.NLAYER * 2.0 - 1.0
    t0c = S5.KAPHEAT * Dc / 2.0
    acc = np.zeros_like(Dc)
    for zz in zlayer:
        sa = S5.slab_seff(np.full_like(Dc, zz), t0c) * Dc / 2.0
        keep = 1.0 if frac is None else (np.abs(zz) >= frac).astype(float)
        acc += keep * Dc * np.interp(S5.T_of_shielding(sa), tgrid, table)
    return acc / S5.NLAYER


def main():
    known = ("--cirrus", "--envelopes", "--cirrus-only")
    unknown = [a for a in sys.argv[1:] if a.startswith("--") and a not in known]
    if unknown:
        print(" unknown option %s; the options are %s" % (", ".join(unknown), ", ".join(known)))
        sys.exit(1)
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    models = args if args else MODELS
    hdr = fits.getheader(CIRRUS)
    Dc = CIRRUS_SCALE * S5.repair(np.array(fits.getdata(CIRRUS), float).squeeze(), "cirrus")
    Dc0 = Dc.copy()                                   # before the central enhancement
    if E.ENH_A > 0.0:
        Dc = E.enhance(Dc, S5.PIX)
        print(" cirrus enhanced toward the center: factor 1 + %.2f [1 + (r/%.0f\")^2]^(-%.2f/2); median now %.3e cm-2"
              % (E.ENH_A, E.ENH_RC, E.ENH_P, np.median(Dc)))
    cols = filament_columns(len(models))
    if "--envelopes" in sys.argv:
        print(" writing the envelope tables from %s x %g" % (CIRRUS, CIRRUS_SCALE))
        write_envelopes(models, cols, Dc)
        return
    if "--cirrus" in sys.argv:
        if os.path.exists(CIRRUS_REF):
            ref = np.array(fits.getdata(CIRRUS_REF), float).squeeze()
            print(" median of %s: %.3e cm-2; median of %s x %g before the enhancement: %.3e cm-2; CIRRUS_SCALE = %.4f"
                  " would match them" % (CIRRUS_REF, np.nanmedian(ref), CIRRUS, CIRRUS_SCALE, np.nanmedian(Dc0),
                                         CIRRUS_SCALE * np.nanmedian(ref) / np.nanmedian(Dc0)))
        print(" after the enhancement: median %.3e, 99th percentile %.3e, maximum %.3e cm-2"
              % (np.median(Dc), np.percentile(Dc, 99), Dc.max()))
        print(" surface density of the cirrus (%s x %g, enhanced) along the axis of each filament: median, 16th and 84th"
              " percentiles" % (CIRRUS, CIRRUS_SCALE))
        for model, col in zip(models, cols):
            print("   %-12s column %4d   %.3e   %.3e   %.3e cm-2" % ((model, col) + tuple(np.percentile(Dc[:, col], [50, 16, 84]))))
        return
    print(" filament axes at columns %s, from left to right" % ", ".join(str(c) for c in cols))

    # the true surface density of the filaments, and the column of the axis in each model
    Dfil = np.zeros_like(Dc)
    axes = []
    for model, col in zip(models, cols):
        img, h = filament_image(model, RTSURF)
        if abs(pixel_arcsec(h) - S5.PIX) > 0.01:
            raise ValueError("%s: pixels of %.3f arcsec, not %.1f" % (model, pixel_arcsec(h), S5.PIX))
        prof = middle_profile(img)
        ax = int(np.argmax(prof))
        axes.append(ax)
        if ENVZ:
            if img.shape != Dc.shape:
                raise ValueError("%s: image %s, sky %s; with ENVZ the images must match the sky" % (model, img.shape, Dc.shape))
            Dfil += shift_columns(oriented(img, model, ax, Dc, col), ax, col)
        else:
            Dfil += place(prof, ax, col, Dc.shape)
        print(" %-12s axis at column %4d of its image, peak %.3e cm-2, placed at column %4d" % (model, ax, prof[ax], col))
    cirrus_only = "--cirrus-only" in sys.argv
    if cirrus_only:
        print(" --cirrus-only: the sky contains the cirrus alone, without filaments and without craters")
        Dfil[:] = 0.0
    frac = crater_fraction(Dc.shape, cols) if (CRATER and not cirrus_only) else None
    if frac is not None:
        print(" crater: cirrus depth %.2f pc, filament radius %.2f pc; at most %.0f%% of the cirrus column removed on the axes"
              % (CIRRUS_DEPTH, ROUT, 100.0 * frac.max()))
    Dc_full = Dc
    if frac is not None:
        Dc = Dc * (1.0 - frac)
    D = Dc + Dfil

    for wave, beam, tag, rtag in S5.BANDS:
        sky = cirrus_emission(Dc_full, wave, frac)
        if RTCIRRUS:
            ratio = rt_cirrus_ratio(wave, sky.shape)
            print("   %g microns: RADMC-3D / analytic cirrus emission, median %.3f, 1st and 99th percentiles %.3f and %.3f"
                  % (wave, np.median(ratio), *np.percentile(ratio, [1, 99])))
            sky = sky * ratio
        for model, col, ax in ([] if cirrus_only else zip(models, cols, axes)):
            img, _ = filament_image(model, RTNAME % wave)
            if ENVZ:
                sky += shift_columns(oriented(img, model, ax, Dc, col), ax, col)
            else:
                sky += place(middle_profile(img), ax, col, Dc.shape)
        out = S5.smooth(sky, beam)
        h = hdr.copy()
        h["BUNIT"] = "MJy/sr"
        h["BZERO"], h["BSCALE"] = 0.0, 1.0
        h["COMMENT"] = "simulated sky: analytic cirrus plus RADMC-3D cylinders as filaments"
        h["CIRSCALE"] = (CIRRUS_SCALE, "factor applied to the cirrus image %s" % CIRRUS)
        for kw in E.keywords():
            h[kw[0]] = (kw[1], kw[2])
        h["TLAW_TF"] = (S5.TFLOOR, "cirrus T(s): floor temperature T_f (K)")
        h["TLAW_T0"] = (S5.TSURF, "cirrus T(s): unshielded temperature T_0 (K)")
        h["TLAW_SS"] = (S5.DSTAR, "cirrus T(s): Sigma_star (cm-2)")
        h["TLAW_G"] = (S5.GAMMA, "cirrus T(s): index gamma")
        name = STEM % (WIDTH, tag, rtag)
        S5.write_fits(name, out, h, [("WAVE", wave, "wavelength (microns)"),
                                     ("BEAM", beam, "half-maximum beam width (arcsec)")])
        print("%-46s median %10.3f MJy/sr" % (name, np.median(out)))

    for nm, arr, bm in (((PREFIX + ".w%d.surfdens.pix.fits") % WIDTH, D, 0.0),
                        ((PREFIX + ".w%d.surfdens.r13p5.fits") % WIDTH, S5.smooth(D, 13.5), 13.5),
                        ((PREFIX + ".w%d.surfdens.r36p3.fits") % WIDTH, S5.smooth(D, 36.3), 36.3),
                        ((PREFIX + ".w%d.filaments.pix.fits") % WIDTH, Dfil, 0.0)):
        h = hdr.copy()
        h["BUNIT"], h["BZERO"], h["BSCALE"] = "cm-2", 0.0, 1.0
        dbl = [("WAVE", 0.0, "wavelength (microns); 0 for surface density")]
        if bm > 0:
            dbl.append(("BEAM", bm, "half-maximum beam width (arcsec)"))
        elif "BEAM" in h:
            del h["BEAM"]
        S5.write_fits(nm, arr, h, dbl)
        print("%-46s median %10.4g cm-2" % (nm, np.median(arr)))

    cfg = ("+hires.cfg." + PREFIX + ".w%d") % WIDTH
    with open(cfg, "w") as f:
        f.write("#\n# HIRES: CONFIGURATION PARAMETERS FOR HIRES\n#\n# PM  WAVE  BEAM  OFFSET   IMAGE   [MASKIMAGE]\n#\n")
        for wave, beam, tag, rtag in S5.BANDS:
            f.write("  %s   %s  %4.1f   0.0   %s\n" % ("-" if wave < 160 else "+", tag, beam, STEM % (WIDTH, tag, rtag)))
    print(" written %s" % cfg)


if __name__ == "__main__":
    main()
