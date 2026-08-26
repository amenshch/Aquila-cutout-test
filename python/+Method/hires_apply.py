#!/usr/bin/env python3
"""
hires_apply.py -- reconstruction of the dust surface density from a set of
Herschel images, with a dust temperature fitted separately for every spatial
scale.

The method keeps the telescoping structure of HIRES.  The observed intensities
are written as a sum of components that are exactly additive: a term smoothed
to the largest scale considered, the differences between successive smoothings
above the coarsest beam, and the differences between successive beams of the
resolution ladder below it.  Every component is available in all wavebands
whose beam is not coarser than the resolution at which the component is
defined, so a modified blackbody with a fixed emissivity index can be fitted to
each component separately.  The surface density is the sum of the fitted
amplitudes, which are linear in the observed intensities and may take either
sign, so that additivity is preserved exactly.

Any set of wavebands may be supplied.  Which components can have their
temperature fitted, and which must borrow it from the next coarser component,
follows from the beams that are actually present.

Usage
-----
    python3 hires_apply.py --config bands.txt --out aquila

The configuration file has one line per waveband, blank lines and lines
beginning with a hash being ignored:

    # wavelength_um  beam_arcsec  image.fits            mask.fits      offset
       70            8.4          aquila_070.fits       cov_070.fits   0.0
      160           13.5          aquila_160.fits       cov_160.fits  -12.3
      250           18.2          aquila_250.fits       cov_250.fits    4.5
      350           24.9          aquila_350.fits       cov_350.fits    2.1
      500           36.3          aquila_500.fits       cov_500.fits    0.8

The mask and the offset may be omitted or given as a hyphen.  Images must be in
MJy sr^-1, on a common pixel grid.  The coverage masks of all wavebands are
multiplied together and the result is contracted by a few pixels, so that only
pixels covered in every waveband are reconstructed, as in the original hires
script.
"""

import argparse
import os
import sys

import numpy as np
from scipy import ndimage

try:
    from astropy.io import fits
except ImportError:                                              # pragma: no cover
    sys.exit("ERROR: this script requires astropy (pip install astropy)")

# ----------------------------------------------------------------------
# constants
# ----------------------------------------------------------------------
H_PLANCK = 6.62607015e-27      # erg s
K_BOLTZ = 1.380649e-16         # erg K^-1
C_LIGHT = 2.99792458e10        # cm s^-1
AMU = 1.6605402e-24            # g
MU_H2 = 2.8                    # mean molecular weight per hydrogen molecule
DUST_TO_GAS = 1.0e-2           # dust-to-gas mass ratio
KAPPA0_DUST = 10.0             # cm^2 per gram of dust at NU0
NU0 = 1.0e12                   # Hz
BETA = 2.0                     # dust emissivity index
T_MIN, T_MAX = 3.0, 100.0      # range within which the temperature is kept
FWHM_TO_SIGMA = 1.0 / 2.354820045


def kappa_abs(wave_um):
    """Dust absorption opacity, cm^2 per gram of dust."""
    nu = C_LIGHT / (wave_um * 1.0e-4)
    return KAPPA0_DUST * (nu / NU0) ** BETA


def planck(temp_k, wave_um):
    """Planck function, erg s^-1 cm^-2 Hz^-1 sr^-1."""
    nu = C_LIGHT / (np.asarray(wave_um, float) * 1.0e-4)
    x = H_PLANCK * nu / (K_BOLTZ * np.asarray(temp_k, float))
    return (2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2) / \
        np.expm1(np.clip(x, 1e-8, 700.0))


def unit_surfden(temp_k, wave_um):
    """Intensity in MJy sr^-1 produced by a surface density of one hydrogen
    molecule per square centimetre at the given dust temperature."""
    return planck(temp_k, wave_um) * kappa_abs(wave_um) * DUST_TO_GAS \
        * MU_H2 * AMU / 1.0e-17


# ----------------------------------------------------------------------
# input
# ----------------------------------------------------------------------
class Band:
    def __init__(self, wave, beam, image_path, mask_path=None, offset=0.0):
        self.wave = float(wave)
        self.beam = float(beam)
        self.image_path = image_path
        self.mask_path = mask_path
        self.offset = float(offset)
        self.data = None
        self.header = None


def read_config(path):
    bands = []
    for raw in open(path):
        line = raw.split("#")[0].strip()
        if not line:
            continue
        f = line.split()
        if len(f) < 3:
            sys.exit("ERROR: malformed line in %s: %s" % (path, raw.rstrip()))
        mask = None if len(f) < 4 or f[3] in ("-", "none") else f[3]
        offset = float(f[4]) if len(f) > 4 and f[4] not in ("-",) else 0.0
        bands.append(Band(f[0], f[1], f[2], mask, offset))
    if len(bands) < 2:
        sys.exit("ERROR: at least two wavebands are required")
    bands.sort(key=lambda b: b.beam)
    return bands


def load(bands, contract, verbose=True):
    """Read the images, add the offsets, and build the common coverage mask."""
    shape = None
    valid = None
    for b in bands:
        with fits.open(b.image_path) as hdul:
            hdu = next(h for h in hdul if h.data is not None)
            b.data = np.array(hdu.data, dtype=float).squeeze()
            b.header = hdu.header
        if shape is None:
            shape = b.data.shape
        elif b.data.shape != shape:
            sys.exit("ERROR: %s has shape %s, expected %s"
                     % (b.image_path, b.data.shape, shape))
        unit = str(b.header.get("BUNIT", "")).strip().lower().replace(" ", "")
        if unit and unit not in ("mjy/sr", "mjysr-1", "mjy/steradian"):
            print("WARNING: BUNIT of %s is '%s', MJy/sr assumed"
                  % (b.image_path, b.header.get("BUNIT")))
        b.data = b.data + b.offset
        ok = np.isfinite(b.data)
        if b.mask_path:
            with fits.open(b.mask_path) as hdul:
                hdu = next(h for h in hdul if h.data is not None)
                m = np.array(hdu.data, dtype=float).squeeze()
            if m.shape != shape:
                sys.exit("ERROR: mask %s has shape %s, expected %s"
                         % (b.mask_path, m.shape, shape))
            ok &= np.isfinite(m) & (m > 0)
        valid = ok if valid is None else (valid & ok)
    # contract the common mask, as the hires script does with its border erosion
    if contract > 0:
        valid = ndimage.binary_erosion(valid, iterations=int(contract))
    if verbose:
        print("common coverage: %d of %d pixels (%.1f per cent)"
              % (valid.sum(), valid.size, 100.0 * valid.sum() / valid.size))
    if valid.sum() == 0:
        sys.exit("ERROR: the common coverage mask is empty")
    return valid


def pixel_size(header):
    """Pixel size in arcsec from the header."""
    for key in ("CDELT2", "CD2_2"):
        if key in header:
            return abs(float(header[key])) * 3600.0
    if "CDELT1" in header:
        return abs(float(header["CDELT1"])) * 3600.0
    sys.exit("ERROR: no CDELT2 or CD2_2 in the header; give --pixel explicitly")


# ----------------------------------------------------------------------
# image operations
# ----------------------------------------------------------------------
def fill_invalid(image, valid):
    """Replace the pixels outside the coverage by the value of the nearest
    covered pixel, so that convolutions are not biased by the boundary.  This
    reproduces the expansion by edge values used by the hires script."""
    if valid.all():
        return image
    idx = ndimage.distance_transform_edt(~valid, return_distances=False,
                                         return_indices=True)
    return image[tuple(idx)]


def convolve(image, fwhm_as, pix_as, valid=None):
    """Convolve with a circular Gaussian of the given full width at half
    maximum, after filling the pixels outside the coverage."""
    if fwhm_as <= 0.0:
        return image.copy()
    a = image if valid is None else fill_invalid(image, valid)
    sig = fwhm_as * FWHM_TO_SIGMA / pix_as
    return ndimage.gaussian_filter(a, sig, mode="nearest", truncate=4.0)


# ----------------------------------------------------------------------
# fitting one component
# ----------------------------------------------------------------------
def robust_sigma(image, valid, nsigma=3.0, iterations=6):
    """Standard deviation of an image after iterative sigma clipping, as the
    single-scale standard deviations are obtained in getsf."""
    v = image[valid]
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 0.0
    s = v.std()
    for _ in range(iterations):
        keep = np.abs(v) < nsigma * s
        if keep.sum() < 16:
            break
        s_new = v[keep].std()
        if s_new <= 0 or abs(s_new - s) < 0.01 * s:
            s = s_new
            break
        s = s_new
    return float(s)


def fit_component(comps, waves, rel_error, tgrid, chunk=200000,
                  valid=None, nsigma=3.0):
    """Fit one modified blackbody to a set of single-scale component images.

    comps : (n_bands, n, n) component intensities in MJy sr^-1, all defined at
            the same pair of resolutions and therefore referring to the same
            spatial scale.
    Returns the fitted dust temperature and a mask of the pixels where the fit
    is meaningful, that is where every waveband has the same sign and the
    smallest of them is not negligible.
    """
    comps = np.asarray(comps, float)
    nb = len(waves)
    shape = comps.shape[1:]
    a = np.abs(comps)
    same = (comps > 0).all(axis=0) | (comps < 0).all(axis=0)
    # significance is judged against a robust standard deviation of each
    # waveband of this component, not against a global maximum, so that faint
    # regions of a map with a large dynamic range are treated on their own terms
    if valid is None:
        valid = np.ones(shape, bool)
    # No mask is applied.  Every hard decision taken pixel by pixel draws a
    # contour in the reconstructed map, and those contours are visible as thin
    # curved artifacts, so the estimator must be continuous in the data
    # everywhere.  The reliability of the fit is expressed later as a weight.
    good = valid & np.ones(shape, bool)

    rel = np.asarray(rel_error, float).reshape((nb,) + (1,) * len(shape))
    temp = np.full(shape, np.nan)
    if nb < 2 or not good.any():
        return temp, np.zeros(shape, bool)

    flat = a.reshape(nb, -1)
    gflat = good.ravel()
    idx = np.flatnonzero(gflat)
    out = np.full(flat.shape[1], np.nan)
    rel1 = np.asarray(rel_error, float).reshape(nb, 1)
    for start in range(0, len(idx), chunk):
        sub = idx[start:start + chunk]
        y = flat[:, sub]
        w = 1.0 / np.maximum(rel1 * y, 1e-30) ** 2
        best = np.full(len(sub), np.inf)
        tb = np.full(len(sub), np.nan)
        for t in tgrid:
            m = np.array([unit_surfden(t, wl) for wl in waves]).reshape(nb, 1)
            s = (w * y * m).sum(axis=0) / np.maximum((w * m * m).sum(axis=0),
                                                     1e-300)
            chi = (w * (y - s * m) ** 2).sum(axis=0)
            upd = chi < best
            best[upd] = chi[upd]
            tb[upd] = t
        # parabolic refinement between the grid points: the grid spacing of
        # 0.5 per cent in temperature would otherwise appear as jumps of
        # several per cent in surface density, since the logarithmic
        # sensitivity of the Planck function reaches nine at 160 micron
        step = np.log(tgrid[1] / tgrid[0])
        for _ in range(3):
            c = []
            for fac in (np.exp(-step), 1.0, np.exp(step)):
                tt = tb * fac
                mm = np.array([unit_surfden(tt, wl) for wl in waves])
                num = (w * y * mm).sum(axis=0)
                den = (w * mm * mm).sum(axis=0)
                ss = num / np.maximum(den, 1e-300)
                c.append((w * (y - ss * mm) ** 2).sum(axis=0))
            c0, c1, c2 = c
            dd = c0 - 2.0 * c1 + c2
            sh = np.where(np.abs(dd) > 1e-300, 0.5 * (c0 - c2) / dd, 0.0)
            tb = tb * np.exp(np.clip(sh, -1.0, 1.0) * step)
            step *= 0.5
        out[sub] = tb
    temp = out.reshape(shape)
    # a fit that lands on either end of the grid is kept at that value rather
    # than rejected: clamping saturates continuously, rejection does not
    return temp, good


def amplitude(comps, waves, temp, rel_error, mode="fit"):
    """Surface density of a component, given its dust temperature.  Linear in
    the intensities, so it may be negative where the component is negative."""
    comps = np.asarray(comps, float)
    nb = len(waves)
    t = np.clip(temp, T_MIN, T_MAX)
    if mode == "longest":
        # the surface density is taken from the longest wavelength alone,
        # divided by the Planck function, as the original hires does; this is
        # less sensitive than a multi-waveband amplitude to the fact that the
        # spectral energy distribution of a component is itself a mixture of
        # temperatures, which a single modified blackbody under-predicts at the
        # longest wavelength
        k = int(np.argmax(waves))
        return comps[k] / unit_surfden(t, waves[k])
    m = np.array([unit_surfden(t, wl) for wl in waves])
    a = np.maximum(np.abs(comps), 1e-6 * max(np.nanmax(np.abs(comps)), 1e-30))
    rel = np.asarray(rel_error, float).reshape((nb,) + (1,) * (comps.ndim - 1))
    w = 1.0 / np.maximum(rel * a, 1e-30) ** 2
    return (w * comps * m).sum(axis=0) / np.maximum((w * m * m).sum(axis=0),
                                                    1e-300)


def clip_component(add, sigma, factor=1.0):
    """Limit the surface density of one component to a multiple of what the
    coarser components together have already contributed.

    The limiter saturates smoothly rather than clipping, because a hard clip
    draws a contour in the map wherever it becomes active, and such contours
    appear as thin curved artifacts.  With factor one the limit is strict, so
    the accumulated surface density can never become negative.
    """
    lim = np.maximum(factor * np.abs(sigma), 1e-30)
    out = lim * np.tanh(add / lim)
    n = 100.0 * np.mean(np.abs(add) > 0.9 * lim)
    return out, n


def incoherent_level(comps, waves, resolution_as, pix_as, valid):
    """Root-mean-square of the part of a component that is not shared between
    the wavebands and lies below the band limit.

    Each waveband of a component is regressed on the longest waveband with a
    single slope, and the residual is high-pass filtered at the resolution of
    the component.  A component cannot contain real structure finer than that
    resolution, so what remains is incoherent: noise, or any mismatch of the
    beams.  Returned in the units of the longest waveband, so that it can be
    compared directly with the amplitude of the component.
    """
    k = int(np.argmax(waves))
    ref = comps[k]
    out = 0.0
    for i in range(len(waves)):
        if i == k:
            continue
        m = valid & np.isfinite(comps[i]) & np.isfinite(ref)
        den = float(np.sum(ref[m] ** 2))
        if den <= 0:
            continue
        a = float(np.sum(comps[i][m] * ref[m])) / den
        res = comps[i] - a * ref
        fine = res - convolve(res, resolution_as, pix_as, valid)
        out = max(out, float(np.std(fine[m])) / max(abs(a), 1e-30))
    return out


def adaptive_temperature(temp, comp_ref, sigma_inc, resolution_as, pix_as,
                         valid, target_k=3.0, width_factor=4.0, dlnb=3.0):
    """Average the temperature over a neighbourhood whose size is set by the
    local amplitude of the component.

    The uncertainty of a temperature fitted from a component of amplitude C
    that carries an incoherent contribution sigma_inc is, to first order,
    sigma_T / T = (sigma_inc / |C|) / dlnb, where dlnb is the difference of the
    logarithmic sensitivities of the Planck function between the wavebands
    used.  At each pixel the smallest width of the ladder is taken at which the
    uncertainty of the averaged temperature falls below target_k, so that a
    strong component such as a dense core keeps its temperature at full
    resolution while a weak one is averaged over as much material as it needs.
    """
    t = np.where(np.isfinite(temp), temp, np.nan)
    a = np.abs(comp_ref)
    w = np.where(valid & np.isfinite(t), a, 0.0)
    filled = np.where(np.isfinite(t), t, 0.0)

    # a single smoothed temperature, weighted by the amplitude of the component
    width = width_factor * resolution_as
    num = convolve(filled * w, width, pix_as, valid)
    den = convolve(w, width, pix_as, valid)
    t_smooth = np.where(den > 0, num / np.maximum(den, 1e-300), np.nan)

    # the amplitude at which the temperature reaches the target uncertainty
    a_ref = np.abs(t_smooth) * sigma_inc / (dlnb * max(target_k, 1e-6))
    a_ref = np.where(np.isfinite(a_ref), a_ref, np.inf)

    # continuous blend: full resolution where the component is strong, the
    # smoothed temperature where it is weak, with no threshold anywhere
    f = a ** 2 / (a ** 2 + a_ref ** 2)
    out = np.where(np.isfinite(t), f * t, 0.0) + \
        np.where(np.isfinite(t_smooth), (1.0 - f) * t_smooth, 0.0)
    both = np.isfinite(t) | np.isfinite(t_smooth)
    out = np.where(both, out, np.nan)
    med = np.nanmedian(out[valid]) if np.isfinite(out[valid]).any() else 15.0
    return np.clip(np.where(np.isfinite(out), out, med), T_MIN, T_MAX)


def smooth_temperature(temp, weight, width_as, pix_as, valid):
    """Smooth a component's temperature with a weighting by its own amplitude.

    The temperature of a component varies on scales of order the scale of the
    component itself, whereas the intensities from which it is fitted are
    faint, so the fitted temperature is noisy and that noise is amplified into
    the surface density by the logarithmic sensitivity of the Planck function.
    Smoothing the temperature and leaving the amplitude at full resolution
    removes that noise without removing any structure the data support.
    """
    if width_as <= 0.0:
        return temp
    w = np.where(valid & np.isfinite(temp), np.abs(weight), 0.0)
    t = np.where(np.isfinite(temp), temp, 0.0)
    num = convolve(w * t, width_as, pix_as, valid)
    den = convolve(w, width_as, pix_as, valid)
    thr = 1e-6 * max(np.nanmax(den), 1e-30)
    return np.where(den > thr, num / np.maximum(den, 1e-300), temp)


def fill_temperature(temp, good, weight, width_as, pix_as, fallback, valid):
    """Replace the dust temperature where the fit was not meaningful, by a
    smoothing of the valid values weighted with the amplitude of the component,
    and where that has no support by the temperature of the coarser scale."""
    w = np.where(good, np.abs(weight), 0.0)
    t = np.where(good, np.nan_to_num(temp), 0.0)
    num = convolve(w * t, width_as, pix_as, valid)
    den = convolve(w, width_as, pix_as, valid)
    thr = 0.02 * max(np.nanmax(den), 1e-30)
    smooth = np.where(den > thr, num / np.maximum(den, 1e-300), np.nan)
    out = np.where(good, temp, smooth)
    out = np.where(np.isfinite(out), out, fallback)
    med = np.nanmedian(out[valid]) if np.isfinite(out[valid]).any() else 15.0
    return np.clip(np.where(np.isfinite(out), out, med), T_MIN, T_MAX)


def reliability(temp, wave_um, sigma_t):
    """Shrinkage of a component whose temperature was borrowed rather than
    fitted.  The surface density derived from a waveband responds to an error
    of the temperature with the logarithmic sensitivity x = h c / (lambda k T),
    so a borrowed temperature of uncertainty sigma_t makes the component
    uncertain by x sigma_t / T; the factor below is the estimator that
    minimizes the mean square error."""
    t = np.clip(temp, T_MIN, T_MAX)
    nu = C_LIGHT / (wave_um * 1.0e-4)
    x = H_PLANCK * nu / (K_BOLTZ * t)
    return 1.0 / (1.0 + (x * sigma_t / t) ** 2)


# ----------------------------------------------------------------------
# the reconstruction
# ----------------------------------------------------------------------
def reconstruct(bands, valid, pix_as, factor=1.2, max_scale=None,
                sed_min_wave=160.0, rel_error=0.2, sigma_t=1.0,
                fill_width=2.0, ntemp=400, nsigma=1.0, tmin=5.0, tmax=40.0,
                clip_factor=2.0, finest_beam=None, amp_mode="longest",
                use_short=False, temp_smooth=1.0,
                verbose=True):
    """Reconstruct the surface density at every resolution of the ladder.

    Returns a dictionary keyed by the beam in arcsec, each entry being the
    surface density in cm^-2 accumulated down to that resolution, and a
    dictionary of the fitted dust temperature of every component.
    """
    nb = len(bands)
    beams = np.array([b.beam for b in bands])
    waves = np.array([b.wave for b in bands])
    coarsest = nb - 1
    tgrid = np.geomspace(tmin, tmax, int(ntemp))
    if finest_beam is None:
        long_enough = [beams[j] for j in range(nb) if waves[j] >= sed_min_wave]
        finest_beam = min(long_enough) if long_enough else beams[0]
    rel = {b.wave: rel_error for b in bands}

    if verbose:
        print("wavebands: " + ", ".join("%g um / %.1f\"" % (w, b)
                                        for w, b in zip(waves, beams)))

    # ---- cross-convolution ------------------------------------------------
    conv = {}
    for t in range(nb):
        for j in range(t + 1):
            f = np.sqrt(max(beams[t] ** 2 - beams[j] ** 2, 0.0))
            conv[(t, j)] = bands[j].data if f <= 0 else \
                convolve(bands[j].data, f, pix_as, valid)

    # ---- which wavebands may enter which component ------------------------
    def component_bands(i):
        if not use_short:
            sed = [j for j in range(nb)
                   if beams[j] <= beams[i] + 1e-6 and waves[j] >= sed_min_wave]
            if len(sed) >= 1:
                return sed, sed
        """Indices of the wavebands usable for the component defined between
        the beam of band i and the next coarser beam."""
        avail = [j for j in range(nb) if beams[j] <= beams[i] + 1e-6]
        sed = [j for j in avail if waves[j] >= sed_min_wave]
        if len(sed) >= 2:
            return sed, sed
        # not enough long wavebands: add the shorter ones, longest first, to
        # constrain the temperature, but never to supply the surface density
        extra = sorted([j for j in avail if waves[j] < sed_min_wave],
                       key=lambda j: -waves[j])
        for j in extra:
            sed_plus = sed + [j]
            if len(sed_plus) >= 2:
                return sorted(sed_plus, key=lambda k: waves[k]), sed or [avail[-1]]
            sed = sed_plus
        return sorted(sed, key=lambda k: waves[k]), sed or [avail[-1]]

    # ---- scales above the coarsest beam -----------------------------------
    if max_scale is None:
        half = 0.5 * min(bands[0].data.shape) * pix_as
        max_scale = max(2.0 * beams[coarsest], 0.6 * half)
    ladder = [beams[coarsest]]
    s = beams[coarsest] * factor
    while s < max_scale:
        ladder.append(s)
        s *= factor
    if verbose:
        print("ladder above %.1f\": %d scales up to %.1f\""
              % (beams[coarsest], len(ladder) - 1, ladder[-1]))

    sed_all = [j for j in range(nb) if waves[j] >= sed_min_wave]
    if len(sed_all) < 2:
        sed_all = list(range(nb))
    smoothed = {}
    for sc in ladder:
        f = np.sqrt(max(sc ** 2 - beams[coarsest] ** 2, 0.0))
        for j in sed_all:
            smoothed[(sc, j)] = conv[(coarsest, j)] if f <= 0 else \
                convolve(conv[(coarsest, j)], f, pix_as, valid)

    sigma = np.zeros_like(bands[0].data)
    temps = {}

    top = np.array([smoothed[(ladder[-1], j)] for j in sed_all])
    wl = [waves[j] for j in sed_all]
    t, g = fit_component(top, wl, [rel[w] for w in wl], tgrid, valid=valid,
                         nsigma=nsigma)
    if temp_smooth > 0:
        si = incoherent_level(top, wl, ladder[-1], pix_as, valid)
        t = adaptive_temperature(t, top[-1], si, ladder[-1], pix_as, valid,
                                 target_k=temp_smooth)
    sigma = sigma + amplitude(top, wl, t, [rel[w] for w in wl], amp_mode)
    temps["top_%.0f" % ladder[-1]] = t
    coarser = t

    for k in range(len(ladder) - 1, 0, -1):
        comp = np.array([smoothed[(ladder[k - 1], j)] - smoothed[(ladder[k], j)]
                         for j in sed_all])
        t, g = fit_component(comp, wl, [rel[w] for w in wl], tgrid,
                             valid=valid, nsigma=nsigma)
        if temp_smooth > 0:
            si = incoherent_level(comp, wl, ladder[k], pix_as, valid)
            t = adaptive_temperature(t, comp[-1], si, ladder[k], pix_as, valid,
                                     target_k=temp_smooth)
        add = amplitude(comp, wl, t, [rel[w] for w in wl], amp_mode)
        add, nclip = clip_component(add, sigma, clip_factor)
        sigma = sigma + add
        temps["scale_%.0f" % ladder[k]] = t
        coarser = t
        if verbose:
            frac = 100.0 * g[valid].mean()
            print("  scale %7.1f\": temperature fitted in %5.1f per cent of "
                  "the covered pixels, %5.2f per cent of amplitudes limited%s"
                  % (ladder[k], frac, nclip,
                     "   <-- too few, lower --nsigma" if frac < 20.0 else ""))

    # ---- the beam ladder --------------------------------------------------
    out = {beams[coarsest]: sigma.copy()}
    for i in range(coarsest - 1, -1, -1):
        if finest_beam is not None and beams[i] < finest_beam - 1e-6:
            if verbose:
                print("  component %.1f\" -> %.1f\": skipped, the "
                      "reconstruction stops at %.1f\""
                      % (beams[i], beams[i + 1], finest_beam))
            continue
        fit_idx, amp_idx = component_bands(i)
        comp_fit = np.array([conv[(i, j)] - conv[(i + 1, j)] for j in fit_idx])
        comp_amp = np.array([conv[(i, j)] - conv[(i + 1, j)] for j in amp_idx])
        wf = [waves[j] for j in fit_idx]
        wa = [waves[j] for j in amp_idx]
        borrowed = len(fit_idx) < 2
        if borrowed:
            t = coarser.copy()
        else:
            t, g = fit_component(comp_fit, wf, [rel[w] for w in wf], tgrid,
                                 valid=valid, nsigma=nsigma)
            if temp_smooth > 0:
                si = incoherent_level(comp_fit, wf, beams[i + 1], pix_as, valid)
                t = adaptive_temperature(t, comp_fit[-1], si, beams[i + 1],
                                         pix_as, valid, target_k=temp_smooth)
        frac = 100.0 * (g[valid].mean() if not borrowed else 0.0)
        add = amplitude(comp_amp, wa, t, [rel[w] for w in wa], amp_mode)
        if borrowed:
            add = add * reliability(t, min(wa), sigma_t)
        add, nclip = clip_component(add, sigma, clip_factor)
        sigma = sigma + add
        temps["beam_%.1f" % beams[i]] = t
        out[beams[i]] = sigma.copy()
        coarser = t
        if verbose:
            print("  component %.1f\" -> %.1f\": %s from %s"
                  % (beams[i], beams[i + 1],
                     "temperature borrowed and shrunk" if borrowed
                     else "temperature fitted",
                     ", ".join("%g" % waves[j] for j in fit_idx)))
            print("      temperature fitted in %5.1f per cent of the covered "
                  "pixels, %5.2f per cent of amplitudes limited%s"
                  % (frac, nclip,
                     "   <-- too many, this component is unreliable"
                     if nclip > 2.0 else ""))
    return out, temps


# ----------------------------------------------------------------------
def classic_surfden(bands, valid, pix_as, sed_min_wave, rel_error, tgrid,
                    amp_mode="fit", verbose=True):
    """Conventional surface density at the coarsest resolution: one modified
    blackbody with a fixed emissivity index fitted to the total intensities of
    all wavebands convolved to the coarsest beam, which is what fitfluxes and
    operate do in the original hires.  Provided for comparison only."""
    nb = len(bands)
    beams = np.array([b.beam for b in bands])
    waves = np.array([b.wave for b in bands])
    coarsest = nb - 1
    use = [j for j in range(nb) if waves[j] >= sed_min_wave]
    if len(use) < 2:
        use = list(range(nb))
    stack = []
    for j in use:
        f = np.sqrt(max(beams[coarsest] ** 2 - beams[j] ** 2, 0.0))
        stack.append(bands[j].data if f <= 0 else
                     convolve(bands[j].data, f, pix_as, valid))
    stack = np.array(stack)
    wl = [waves[j] for j in use]
    rel = [rel_error] * len(use)
    temp, good = fit_component(stack, wl, rel, tgrid, valid=valid, nsigma=0.0)
    temp = np.where(np.isfinite(temp), temp,
                    np.nanmedian(temp[good]) if good.any() else 15.0)
    sd = amplitude(stack, wl, temp, rel, amp_mode)
    if verbose:
        v = sd[valid]
        print("classic surface density at %.1f\": median %.3e cm^-2, "
              "temperature median %.2f K"
              % (beams[coarsest], np.median(v), np.median(temp[valid])))
    return sd, temp


def write_fits(path, data, header, valid, bunit, comment):
    out = np.where(valid, data, np.nan)
    hdr = header.copy()
    hdr["BUNIT"] = bunit
    hdr["BZERO"] = 0.0
    hdr["BSCALE"] = 1.0
    hdr.add_history(comment)
    fits.PrimaryHDU(out.astype(np.float32), hdr).writeto(path, overwrite=True)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True,
                   help="file listing the wavebands, one per line")
    p.add_argument("--out", required=True, help="prefix of the output files")
    p.add_argument("--pixel", type=float, default=None,
                   help="pixel size in arcsec; taken from the header if absent")
    p.add_argument("--factor", type=float, default=1.2,
                   help="geometric factor of the ladder of scales above the "
                        "coarsest beam (default: %(default)s)")
    p.add_argument("--max-scale", type=float, default=None,
                   help="largest scale of that ladder, arcsec")
    p.add_argument("--sed-min-wave", type=float, default=160.0,
                   help="shortest wavelength admitted to a fit when two or "
                        "more longer wavebands are available (default: "
                        "%(default)s micron)")
    p.add_argument("--rel-error", type=float, default=0.2,
                   help="relative uncertainty assigned to every waveband")
    p.add_argument("--sigma-t", type=float, default=1.0,
                   help="assumed uncertainty in K of a borrowed temperature")
    p.add_argument("--contract", type=int, default=2,
                   help="pixels by which the common coverage mask is contracted")
    p.add_argument("--temp-smooth", type=float, default=3.0,
                   help="target uncertainty in kelvin of the temperature of "
                        "each component.  The temperature is averaged over the "
                        "smallest neighbourhood at which this is reached, so a "
                        "strong component keeps full resolution and a weak one "
                        "is averaged over more material.  Zero disables the "
                        "averaging (default: %(default)s)")
    p.add_argument("--amplitude", default="longest", choices=["fit", "longest"],
                   help="how the surface density of a component is obtained "
                        "from its intensities once the temperature is known: "
                        "'fit' uses all wavebands, 'longest' divides the "
                        "longest wavelength by the Planck function, as the "
                        "original hires does (default: %(default)s)")
    p.add_argument("--finest-beam", type=float, default=None,
                   help="stop the reconstruction at this resolution in arcsec.  "
                        "By default it is the beam of the shortest waveband at "
                        "or above --sed-min-wave, so that wavebands shorter "
                        "than that add no resolution of their own")
    p.add_argument("--use-short-for-temperature", action="store_true",
                   help="allow wavebands shorter than --sed-min-wave to enter "
                        "the temperature fit of the finest component.  On the "
                        "radiative-transfer models this helps; on Aquila the "
                        "70 micron difference image was too noisy and it did "
                        "not, so it is off by default")
    p.add_argument("--clip-factor", type=float, default=1.0,
                   help="a component may not contribute more than this "
                        "multiple of the surface density already accumulated "
                        "(default: %(default)s)")
    p.add_argument("--nsigma", type=float, default=0.0,
                   help="significance in robust standard deviations that a "
                        "component must reach in every waveband before its "
                        "temperature is fitted (default: %(default)s)")
    p.add_argument("--tmin", type=float, default=5.0,
                   help="lowest dust temperature of the grid, K")
    p.add_argument("--tmax", type=float, default=40.0,
                   help="highest dust temperature of the grid, K")
    p.add_argument("--ntemp", type=int, default=400,
                   help="number of points of the temperature grid")
    p.add_argument("--classic", action="store_true",
                   help="also compute the conventional surface density at the "
                        "coarsest resolution, from a single modified blackbody "
                        "fitted to the total intensities, for comparison")
    p.add_argument("--save-temps", action="store_true",
                   help="also write the dust temperature of every component")
    p.add_argument("--all-levels", action="store_true",
                   help="write the surface density at every resolution, not "
                        "only at the finest")
    a = p.parse_args()

    bands = read_config(a.config)
    valid = load(bands, a.contract)
    pix = a.pixel if a.pixel else pixel_size(bands[0].header)
    print("pixel size: %.3f arcsec" % pix)

    out, temps = reconstruct(bands, valid, pix, factor=a.factor,
                             max_scale=a.max_scale,
                             sed_min_wave=a.sed_min_wave,
                             rel_error=a.rel_error, sigma_t=a.sigma_t,
                             ntemp=a.ntemp, nsigma=a.nsigma,
                             tmin=a.tmin, tmax=a.tmax,
                             clip_factor=a.clip_factor,
                             finest_beam=a.finest_beam,
                             amp_mode=a.amplitude,
                             use_short=a.use_short_for_temperature,
                             temp_smooth=a.temp_smooth)

    hdr = bands[0].header
    finest = min(out)
    tag = ("%.1f" % finest).replace(".", "p")
    name = "%s.surfden.r%s.fits" % (a.out, tag)
    write_fits(name, out[finest], hdr, valid, "cm-2",
               "per-scale surface density, finest resolution %.1f arcsec" % finest)
    print("written %s" % name)
    if a.all_levels:
        for beam, data in sorted(out.items()):
            if beam == finest:
                continue
            tg = ("%.1f" % beam).replace(".", "p")
            nm = "%s.surfden.r%s.fits" % (a.out, tg)
            write_fits(nm, data, hdr, valid, "cm-2",
                       "per-scale surface density at %.1f arcsec" % beam)
            print("written %s" % nm)
    if a.classic:
        sd, tt = classic_surfden(bands, valid, pix,
                                 a.sed_min_wave, a.rel_error,
                                 np.geomspace(a.tmin, a.tmax, int(a.ntemp)),
                                 a.amplitude)
        coarse = max(out)
        tg = ("%.1f" % coarse).replace(".", "p")
        write_fits("%s.surfden.classic.r%s.fits" % (a.out, tg), sd, hdr, valid,
                   "cm-2", "conventional single-temperature surface density")
        write_fits("%s.tempers.classic.r%s.fits" % (a.out, tg), tt, hdr, valid,
                   "K", "conventional single-temperature dust temperature")
        ratio = np.where(valid, out[coarse] / np.where(sd != 0, sd, np.nan),
                         np.nan)
        r = ratio[np.isfinite(ratio)]
        print("per-scale / classic at %.1f\": 16th %.3f, median %.3f, "
              "84th %.3f" % (coarse, *np.percentile(r, [16, 50, 84])))
        print("written %s.surfden.classic.r%s.fits" % (a.out, tg))

    if a.save_temps:
        for key, data in temps.items():
            nm = "%s.tempers.%s.fits" % (a.out, key)
            write_fits(nm, data, hdr, valid, "K",
                       "dust temperature of the component %s" % key)
            print("written %s" % nm)


if __name__ == "__main__":
    main()
