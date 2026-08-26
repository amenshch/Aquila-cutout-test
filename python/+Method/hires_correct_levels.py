#!/usr/bin/env python3
"""
hires_correct_levels.py -- correct a hires surface density map for the
single-temperature bias, level by level, without detecting anything.

The map and the waveband images are decomposed into levels by a nested set of
morphological openings.  Picture a ball rolled underneath the surface of an
image: the background it traces is everything wide enough to hold the ball up,
and what is left above it is everything narrower.  With balls of increasing
diameter this gives

    I  =  O_N(I)  +  sum over k of L_k ,      L_k = O_(k-1)(I) - O_k(I),

where O_0 is the identity.  Each level is non-negative, the levels sum exactly
to the image, and material appears in the level that matches its own width.

For every level a dust temperature is fitted to the levels of the waveband
images alone, so it describes the material of that width by itself, with the
warmer material in front of and behind it removed.  That temperature is what
decides the correction:

    Sigma_corrected  =  O_N(Sigma)  +  sum over k of C(T_k) * L_k(Sigma) ,

with C calibrated on the Bonnor-Ebert radiative transfer grid as

    log10 C = 0.7477 - 0.5394 log10(T / K),

which gives 1.96 at 7 K, 1.71 at 9 K, 1.53 at 11 K and 1.18 at 18 K, with a
scatter of 0.117 dex.  The widest level, the cirrus, is returned unchanged.

Nothing is classified as a source or as a filament.  Material is assigned to a
level by its width, which is well defined, and the temperature of that level
decides how much correction it receives, so a wide source and a filament of the
same width are treated identically, as the physics requires.  The correction is
multiplicative on non-negative levels, so it can neither create a depression nor
change the shape of a structure.

Usage
-----
    python3 hires_correct_levels.py --map hi.surface.density.r13p5.fits \\
        --config bands.txt --out corrected --save-diagnostics
"""

import argparse
import sys

import numpy as np
from scipy import ndimage

try:
    from astropy.io import fits
except ImportError:                                              # pragma: no cover
    sys.exit("ERROR: this script requires astropy (pip install astropy)")

H_PLANCK = 6.62607015e-27
K_BOLTZ = 1.380649e-16
C_LIGHT = 2.99792458e10
AMU = 1.6605402e-24
MU_H2 = 2.8
DUST_TO_GAS = 1.0e-2
KAPPA0_DUST = 10.0
NU0 = 1.0e12
BETA = 2.0
FWHM_TO_SIGMA = 1.0 / 2.354820045

# calibration of the correction against the colour temperature of a level
CAL_SLOPE, CAL_CONST = -0.5394, 0.7477


def invert_planck(intensity_mjysr, sigma, wave_um):
    """Dust temperature satisfying I = sigma B(T) kappa eta mu m_H."""
    nu = C_LIGHT / (wave_um * 1.0e-4)
    kap = KAPPA0_DUST * (nu / NU0) ** BETA
    b = 1.0e-17 * intensity_mjysr / np.maximum(sigma, 1e-30) \
        / (kap * DUST_TO_GAS * MU_H2 * AMU)
    pref = 2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2
    x = np.log1p(pref / np.maximum(b, 1e-300))
    return H_PLANCK * nu / (K_BOLTZ * np.maximum(x, 1e-12))


def planck(t, wave_um):
    nu = C_LIGHT / (np.asarray(wave_um, float) * 1.0e-4)
    x = H_PLANCK * nu / (K_BOLTZ * np.asarray(t, float))
    return (2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2) / \
        np.expm1(np.clip(x, 1e-8, 700.0))


def unit_surfden(t, wave_um):
    nu = C_LIGHT / (wave_um * 1.0e-4)
    kap = KAPPA0_DUST * (nu / NU0) ** BETA
    return planck(t, wave_um) * kap * DUST_TO_GAS * MU_H2 * AMU / 1.0e-17


def fill_invalid(image, valid):
    """Replace pixels outside the coverage by the nearest covered value, so
    that neither the openings nor the convolutions see the boundary."""
    if valid.all():
        return image
    idx = ndimage.distance_transform_edt(~valid, return_distances=False,
                                         return_indices=True)
    return image[tuple(idx)]


def smooth(a, fwhm_as, pix):
    if fwhm_as <= 0:
        return a.copy()
    return ndimage.gaussian_filter(a, fwhm_as * FWHM_TO_SIGMA / pix,
                                   mode="nearest", truncate=4.0)


def _disk(r):
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= r * r


def opening(a, diameter_as, pix, rmax=8):
    """Morphological opening with a disk, on a decimated grid when the disk is
    large.  The cost of an opening grows with the area of the disk, so a large
    ball is applied to a coarser grid, where it spans few pixels; the result is
    a smooth surface and nothing is lost by doing so."""
    r = max(0.5 * diameter_as / pix, 1.0)
    if r <= rmax:
        bg = ndimage.grey_opening(a, footprint=_disk(int(round(r))))
    else:
        step = int(np.ceil(r / rmax))
        small = a[::step, ::step]
        bs = ndimage.grey_opening(small, footprint=_disk(int(round(r / step))))
        bg = ndimage.zoom(bs, (a.shape[0] / bs.shape[0],
                               a.shape[1] / bs.shape[1]), order=1)
        if bg.shape != a.shape:
            out = np.empty_like(a)
            out[:] = np.pad(bg, ((0, max(0, a.shape[0] - bg.shape[0])),
                                 (0, max(0, a.shape[1] - bg.shape[1]))),
                            mode="edge")[:a.shape[0], :a.shape[1]]
            bg = out
    return np.minimum(smooth(bg, 0.4 * diameter_as, pix), a)


def levels(a, diameters, pix):
    """Return the list of levels and the widest background."""
    out, prev = [], a
    for d in diameters:
        bg = opening(prev, d, pix)
        out.append(np.maximum(prev - bg, 0.0))
        prev = bg
    return out, prev


def incoherent_amplitude(level_images, valid):
    """Amplitude of the part of a level that is not shared between wavebands.

    Instrumental noise is independent from waveband to waveband, whereas real
    structure is common to all of them.  Regressing every waveband on the
    longest one and taking the robust scatter of the residual therefore
    measures the noise of that level, in the units of the longest waveband,
    without any noise map being supplied.
    """
    ref = level_images[-1]
    out = 0.0
    for im in level_images[:-1]:
        m = valid & np.isfinite(im) & np.isfinite(ref)
        den = float(np.sum(ref[m] ** 2))
        if den <= 0:
            continue
        a = float(np.sum(im[m] * ref[m])) / den
        res = (im - a * ref)[m]
        mad = np.median(np.abs(res - np.median(res)))
        out = max(out, 1.4826 * float(mad) / max(abs(a), 1e-30))
    return out


def fit_level_temperature(level_images, waves, tmin, tmax, ntemp, rel=0.2):
    comps = np.array(level_images)
    best = np.full(comps.shape[1:], np.inf)
    tb = np.full(comps.shape[1:], np.nan)
    scale = max(np.nanmedian(comps[-1]), 1e-6)
    for t in np.geomspace(tmin, tmax, int(ntemp)):
        m = np.array([unit_surfden(t, l) for l in waves])[:, None, None]
        w = 1.0 / np.maximum(rel * np.maximum(comps, 1e-3 * scale), 1e-30) ** 2
        s = (w * comps * m).sum(axis=0) / np.maximum((w * m * m).sum(axis=0),
                                                     1e-300)
        chi = (w * (comps - s * m) ** 2).sum(axis=0)
        u = chi < best
        best[u] = chi[u]
        tb[u] = t
    return tb


def read_config(path):
    bands = []
    for raw in open(path):
        line = raw.split("#")[0].strip()
        if not line:
            continue
        f = line.split()
        bands.append((float(f[0]), float(f[1]), f[2],
                      None if len(f) < 4 or f[3] in ("-", "none") else f[3],
                      float(f[4]) if len(f) > 4 and f[4] != "-" else 0.0))
    return sorted(bands, key=lambda b: b[1])


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--map", required=True, help="hires surface density map")
    p.add_argument("--config", required=True, help="waveband list, as for hires_apply")
    p.add_argument("--out", default=None,
                   help="prefix of the output files.  If omitted it is the "
                        "name of the input map without its extension, with "
                        "'.corrected' appended; the temperature map is written "
                        "beside it with 'surface.density' replaced by "
                        "'tempers' where that appears in the name")
    p.add_argument("--diameters", default="auto",
                   help="ball diameters in arcsec, increasing, or 'auto' for a "
                        "geometric ladder from twice the beam to --max-ball "
                        "with a ratio of --ladder-factor (default: %(default)s)")
    p.add_argument("--ladder-factor", type=float, default=1.5,
                   help="ratio between successive ball diameters when the "
                        "ladder is automatic (default: %(default)s)")
    p.add_argument("--max-ball", type=float, default=300.0,
                   help="largest ball in arcsec; structures wider than this "
                        "are treated as cirrus and left unchanged")
    p.add_argument("--protostar-image", default=None,
                   help="a 70 micron image, used only as a flag: where it "
                        "shows a compact source the correction is switched off "
                        "smoothly, because the calibration is derived from "
                        "starless cores and does not apply to internally "
                        "heated ones.  If not given, the shortest waveband of "
                        "the configuration below --min-wave is used")
    p.add_argument("--no-protostar-flag", action="store_true",
                   help="do not use any 70 micron flag")
    p.add_argument("--protostar-contrast", type=float, default=1.0,
                   help="a compact 70 micron source must also exceed this "
                        "fraction of its own local background before the "
                        "correction is switched off")
    p.add_argument("--protostar-radius", type=float, default=45.0,
                   help="radius in arcsec over which the correction is "
                        "switched off around a compact 70 micron source; it "
                        "must cover the heated envelope, not only the peak, or "
                        "a ring is left where the correction is still applied")
    p.add_argument("--protostar-sigma", type=float, default=10.0,
                   help="significance above which a compact 70 micron source "
                        "switches the correction off")
    p.add_argument("--min-wave", type=float, default=160.0,
                   help="shortest waveband used; 70 micron is far too noisy "
                        "for this measurement and must be excluded")
    p.add_argument("--pixel", type=float, default=None)
    p.add_argument("--tmin", type=float, default=6.0)
    p.add_argument("--tmax", type=float, default=30.0)
    p.add_argument("--ntemp", type=int, default=240)
    p.add_argument("--noise-sigma", type=float, default=2.0,
                   help="a level must stand this many times above its own "
                        "incoherent part, measured as the scatter between "
                        "wavebands, before it is corrected; zero disables the "
                        "test (default: %(default)s)")
    p.add_argument("--calibration", default="physical",
                   choices=["physical", "empirical"],
                   help="'physical' converts each level with its own measured "
                        "temperature instead of the single temperature of the "
                        "pixel, which requires no calibration and returns "
                        "exactly one where there is no temperature mixing; "
                        "'empirical' uses the relation fitted to the "
                        "Bonnor-Ebert grid (default: %(default)s)")
    p.add_argument("--c-max", type=float, default=2.5,
                   help="the correction factor is never allowed above this")
    p.add_argument("--smooth-factor", type=float, default=0.0,
                   help="smooth the temperature of each level on this multiple "
                        "of the ball diameter.  With the physical calibration "
                        "it must be left at zero: the pixel temperature is a "
                        "per-pixel quantity, so the level temperature must be "
                        "too, or the two no longer cancel where there is no "
                        "temperature mixing and a spurious correction appears "
                        "wherever the temperature has a gradient, which is "
                        "what produces rings around warm sources")
    p.add_argument("--contract", type=int, default=2,
                   help="pixels by which the common coverage mask is contracted")
    p.add_argument("--edge-taper", type=float, default=1.0,
                   help="the correction is reduced near the edge of the "
                        "coverage over this multiple of the largest ball, "
                        "since a structure cut by the edge is not measured "
                        "correctly; zero disables the taper")
    p.add_argument("--beam", type=float, default=13.5,
                   help="resolution of the map in arcsec.  The correction "
                        "field is smoothed on this scale, because a "
                        "temperature has no structure below the beam and the "
                        "map has no real structure there either")
    p.add_argument("--save-diagnostics", action="store_true")
    a = p.parse_args()

    if a.out is None:
        base = a.map[:-5] if a.map.endswith(".fits") else a.map
        a.out = base + ".corrected"
        print("output prefix from the input map: %s" % a.out)
    if a.diameters.strip().lower() == "auto":
        diam, d0 = [], 2.0 * a.beam
        while d0 <= a.max_ball * 1.0001:
            diam.append(round(d0, 1))
            d0 *= a.ladder_factor
    else:
        diam = [float(x) for x in a.diameters.replace(",", " ").split()]
    print("ball diameters: " + ", ".join("%g" % d for d in diam) + " arcsec")
    with fits.open(a.map) as hdul:
        hdu = next(h for h in hdul if h.data is not None)
        smap = np.array(hdu.data, float).squeeze()
        hdr = hdu.header
    pix = a.pixel
    if pix is None:
        for key in ("CDELT2", "CD2_2", "CDELT1"):
            if key in hdr:
                pix = abs(float(hdr[key])) * 3600.0
                break
    # a pixel of the map that is not positive carries no information; the
    # outermost rows of a reconstruction are often zero, and inverting the
    # Planck function there would give the clamp of 100 K and a surface density
    # of zero in the output
    valid = np.isfinite(smap) & (smap > 0.0)

    bands = [b for b in read_config(a.config) if b[0] >= a.min_wave]
    if len(bands) < 2:
        sys.exit("ERROR: at least two wavebands at or above %.0f micron are "
                 "needed" % a.min_wave)
    coarsest = max(b[1] for b in bands)
    raw = []
    for w, beam, path, mask, off in bands:
        d = np.array(fits.getdata(path), float).squeeze() + off
        ok = np.isfinite(d)
        if mask:
            mm = np.array(fits.getdata(mask), float).squeeze()
            ok &= np.isfinite(mm) & (mm > 0)
        valid = valid & ok
        raw.append(d)
    if a.contract > 0:
        valid = ndimage.binary_erosion(valid, iterations=int(a.contract))
    print("common coverage: %.1f per cent of the map"
          % (100.0 * valid.mean()))
    imgs = []
    for (w, beam, path, mask, off), d in zip(bands, raw):
        imgs.append(smooth(fill_invalid(np.where(np.isfinite(d), d, 0.0), valid),
                           np.sqrt(max(coarsest ** 2 - beam ** 2, 0.0)), pix))
    finite = valid
    filled = fill_invalid(np.where(valid, smap, 0.0), valid)
    waves = [b[0] for b in bands]
    print("map %s, %s pixels of %.2f arcsec" % (a.map, smap.shape, pix))
    print("wavebands used: " + ", ".join("%g" % w for w in waves))

    lev_img, _ = levels(np.array(imgs), diam, pix) if False else (None, None)
    per_band = [levels(im, diam, pix)[0] for im in imgs]
    lev_map, cirrus = levels(filled, diam, pix)

    # The temperature that the conventional reconstruction used in each pixel.
    # The map was built as I / (B(T) kappa eta mu m_H) from the longest
    # waveband, so that temperature can be recovered exactly by inversion.
    long_i = imgs[int(np.argmax(waves))]
    long_w = max(waves)
    t_fit = invert_planck(long_i, filled, long_w)
    t_fit = np.clip(np.where(np.isfinite(t_fit), t_fit, 15.0), 3.0, 100.0)

    # everything above the widest level is returned unchanged, so the
    # correction is expressed as a multiplicative field
    numer = cirrus.copy()
    diag = {"Tfit": t_fit}
    for k, d in enumerate(diam):
        t = fit_level_temperature([pb[k] for pb in per_band], waves,
                                  a.tmin, a.tmax, a.ntemp)
        w = lev_map[k]
        num = smooth(np.where(np.isfinite(t), t, 0.0) * w,
                     a.smooth_factor * d, pix)
        den = smooth(w, a.smooth_factor * d, pix)
        ts = np.where(den > 0, num / np.maximum(den, 1e-300), np.nan)
        ts = np.where(np.isfinite(ts), ts, np.nanmedian(t[np.isfinite(t)]))
        if a.calibration == "physical":
            # The level was converted into surface density with the single
            # temperature of the pixel; converting it with its own temperature
            # instead multiplies it by the ratio of the Planck functions.  No
            # empirical calibration enters, and where every level shares the
            # temperature of the pixel, as in a sky with one temperature per
            # line of sight, the factor is exactly one and nothing is changed.
            c = planck(t_fit, long_w) / np.maximum(planck(ts, long_w), 1e-300)
        else:
            c = 10.0 ** (CAL_CONST + CAL_SLOPE * np.log10(
                np.clip(ts, a.tmin, a.tmax)))
        c = np.clip(c, 1.0, a.c_max)
        # A level that is dominated by noise carries no temperature
        # information, and since the factor is convex in the temperature and is
        # not allowed below one, random scatter would produce a systematic
        # excess.  The correction is therefore faded out where the level does
        # not stand above its own incoherent part.
        if a.noise_sigma > 0:
            sig = incoherent_amplitude([pb[k] for pb in per_band], valid)
            ref = per_band[-1][k]
            thr = a.noise_sigma * sig
            wgt = ref ** 2 / (ref ** 2 + thr ** 2)
            c = 1.0 + (c - 1.0) * wgt
            diag["N%.0f" % d] = wgt
        numer = numer + c * w
        diag["T%.0f" % d] = ts
        diag["C%.0f" % d] = c
        sel = w > 0.05 * np.nanmax(w)
        print("  level below %5.0f arcsec: carries %5.1f per cent of the map, "
              "temperature median %5.2f K, factor median %4.2f"
              % (d, 100.0 * np.nansum(w) / max(np.nansum(filled), 1e-30),
                 np.median(ts[sel]) if sel.any() else np.nan,
                 np.median(c[sel]) if sel.any() else np.nan))
    # a compact 70 micron source is internally heated; the calibration comes
    # from starless cores and must not be applied there
    proto = a.protostar_image
    if a.no_protostar_flag:
        proto = None
    elif proto is None:
        short = [b for b in read_config(a.config) if b[0] < a.min_wave]
        if short:
            proto = min(short, key=lambda b: b[0])[2]
            print("  70 micron flag taken from the configuration: %s" % proto)
    if proto:
        p70 = np.array(fits.getdata(proto), float).squeeze()
        p70 = fill_invalid(np.where(np.isfinite(p70), p70, 0.0), valid)
        bg70 = opening(p70, diam[0], pix)
        c70 = p70 - bg70
        med = np.median(np.abs(c70[valid] - np.median(c70[valid])))
        sig = 1.4826 * max(med, 1e-30)
        # a protostar must stand out both above the noise and above its own
        # local background; the second condition is what keeps bright cirrus
        # structure at 70 micron from switching the correction off
        thr = np.maximum(a.protostar_sigma * sig,
                         a.protostar_contrast * np.maximum(bg70, 0.0))
        warm = c70 / np.maximum(thr, 1e-30)
        keep = 1.0 / (1.0 + np.maximum(warm, 0.0) ** 2)
        # The suppression must cover the whole heated region, not only the
        # pixels where the 70 micron source peaks.  Taking the minimum over a
        # disk first spreads it across the source, so that the correction is
        # switched off over a smooth bowl rather than in a ring around a hole,
        # which is what produced the donut shapes.
        r = int(max(0.5 * a.protostar_radius / pix, 1))
        keep = ndimage.grey_erosion(keep, footprint=_disk(r))
        # smoothed over the same radius again, so that the boundary of the
        # suppressed region is a gradual slope and not a step that would show
        # as a ring in the correction field
        keep = smooth(keep, 1.5 * a.protostar_radius, pix)
        print("  70 micron flag: correction reduced below half in %.2f per cent "
              "of the covered pixels" % (100.0 * np.mean(keep[valid] < 0.5)))
        numer = filled + (numer - filled) * keep

    # The correction is a temperature effect, so it has no structure below the
    # beam, and the map has none either; smoothing the field on the beam
    # prevents any differential amplification of sub-beam noise.
    field = numer / np.maximum(filled, 1e-30)
    field = smooth(fill_invalid(np.where(valid, field, 1.0), valid), a.beam, pix)
    if a.edge_taper > 0:
        near = smooth(valid.astype(float), a.edge_taper * max(diam), pix)
        taper = np.clip((near - 0.5) / 0.5, 0.0, 1.0)
        field = 1.0 + (field - 1.0) * taper
    corrected = filled * np.maximum(field, 1.0)
    ratio = corrected / np.maximum(filled, 1e-30)
    print("corrected / original: 50th %.3f, 90th %.3f, 99th %.3f, max %.3f"
          % (*np.nanpercentile(ratio[finite], [50, 90, 99]),
             np.nanmax(ratio[finite])))

    def write(name, arr, unit, note):
        h = hdr.copy()
        h["BUNIT"] = unit
        h["BZERO"] = 0.0
        h["BSCALE"] = 1.0
        h.add_history(note)
        fits.PrimaryHDU(np.where(finite, arr, np.nan).astype(np.float32),
                        h).writeto(name, overwrite=True)
        print("written %s" % name)

    factor = np.maximum(field, 1.0)
    write("%s.fits" % a.out, corrected, "cm-2",
          "hires map corrected level by level for the single-temperature bias")
    fname = a.out
    if fname.endswith(".corrected"):
        fname = fname[:-len(".corrected")] + ".corrfactor"
    else:
        fname = fname + ".corrfactor"
    write("%s.fits" % fname, factor, "", "multiplicative correction field")
    # The dust temperature consistent with the corrected surface density.
    #
    # The conventional map is I / (B(T) kappa eta mu m_H) with the effective
    # temperature T_fit written out as the diagnostic Tfit, so multiplying the
    # map by C is the same as dividing the Planck function by C.  Defining the
    # corrected temperature that way makes it exactly consistent with the
    # corrected map and identical to T_fit wherever the correction is unity, so
    # that no feature can appear where nothing has been corrected.
    nu_ref = C_LIGHT / (long_w * 1.0e-4)
    pref = 2.0 * H_PLANCK * nu_ref ** 3 / C_LIGHT ** 2
    t_safe = np.clip(np.where(np.isfinite(t_fit), t_fit, 15.0), 3.0, 100.0)
    x_fit = np.clip(H_PLANCK * nu_ref / (K_BOLTZ * t_safe), 1e-8, 700.0)
    b_cor = (pref / np.expm1(x_fit)) / np.maximum(factor, 1e-30)
    tcor = H_PLANCK * nu_ref / (K_BOLTZ * np.maximum(
        np.log1p(pref / np.maximum(b_cor, 1e-300)), 1e-12))
    tcor = np.where(np.isfinite(tcor), tcor, t_safe)
    tname = a.out.replace("surface.density", "temperature")
    if tname == a.out:
        tname = a.out + ".temperature"
    write("%s.fits" % tname, np.clip(tcor, 3.0, 100.0), "K",
          "dust temperature consistent with the corrected surface density")
    if a.save_diagnostics:
        write("%s.cirrus.fits" % a.out, cirrus, "cm-2", "widest background")
        for key, arr in diag.items():
            write("%s.%s.fits" % (a.out, key), arr,
                  "K" if key.startswith("T") else "", key)


if __name__ == "__main__":
    main()
