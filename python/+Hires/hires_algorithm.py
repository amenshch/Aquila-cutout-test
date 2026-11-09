"""
hires_algorithm.py -- faithful Python reimplementation of the HIRES surface
density reconstruction of Men'shchikov (2021), as it is actually carried out by
the Bash script 'hires' and the Fortran utilities 'convolve', 'fitfluxes',
'operate' and 'imgstat'.

The reconstruction is

    D_H = D_500  +  sum over wavebands i  of  max over temperatures n of
          [ D(I_i, T_n) - G * D(I_i, T_n) ],

with D(I, T) = I / (B(T) kappa eta mu m_H), G the Gaussian kernel that takes
the beam of waveband i to the beam of the next longer waveband, and the
differential term clipped at zero for wavebands at and below 160 micron.

Every step can be switched off individually, so that its contribution to the
error of the reconstructed surface density can be measured separately.  The
switches are collected in the class Options.
"""

import numpy as np
from scipy import ndimage

from hires_forward import (WAVES, BEAMS, PIX_ARCSEC, FWHM_TO_SIGMA,
                           kappa_abs, planck, gaussian_convolve,
                           DUST_TO_GAS, MU_H2, AMU)

# index of the finest waveband entering the spectral fits; wavebands with a
# smaller index (70 and 100 micron) carry structure but are not fitted.
NWAVEFIN = 3                       # 0:70  1:100  2:160  3:250  4:350  5:500
TEMP_INDICES = (3, 4, 5)           # resolutions 18.2, 24.9, 36.3 arcsec
TEMP_BANDS = {3: (2, 3),           # T_2 from 160 and 250 micron
              4: (2, 3, 4),        # T_3 from 160, 250 and 350 micron
              5: (2, 3, 4, 5)}     # T_4 from 160, 250, 350 and 500 micron

T_CLAMP = (3.0, 100.0)             # the clamp imposed inside 'operate'


class Options:
    """Switches of the reconstruction.  The defaults reproduce 'hires'."""

    def __init__(self, combine="max", clip=True, smooth_temperature=True,
                 weighting="relative", convolve_ratio=True, exact_temperature=False,
                 rel_error=0.20, noise_mjysr=None):
        # combine: 'max'  pixelwise maximum over the temperature variants
        #          'mean' inverse-variance weighted mean over the same variants
        self.combine = combine
        # clip: clip the differential term at zero for wavebands <= 160 micron
        self.clip = clip
        # smooth_temperature: circular median then mean filter of the
        #          temperature maps, over a window of one beam in diameter
        self.smooth_temperature = smooth_temperature
        # weighting: 'relative' gives sigma = rel_error * I in every waveband,
        #          as 'fitfluxes' does; 'noise' adds a constant noise term
        self.weighting = weighting
        self.rel_error = rel_error
        self.noise_mjysr = noise_mjysr
        # convolve_ratio: form the differential term by convolving the ratio
        #          I / B(T), as 'hires' does; if False, convolve the intensity
        #          first and divide afterwards
        self.convolve_ratio = convolve_ratio
        # exact_temperature: replace the fitted temperature by the temperature
        #          that reproduces the true surface density exactly
        self.exact_temperature = exact_temperature

    def label(self):
        tags = []
        if self.combine != "max":
            tags.append("combine=%s" % self.combine)
        if not self.clip:
            tags.append("no clipping")
        if not self.smooth_temperature:
            tags.append("no temperature filtering")
        if self.weighting != "relative":
            tags.append("weighting=%s" % self.weighting)
        if np.ndim(self.rel_error) > 0:
            tags.append("per-band uncertainties")
        if not self.convolve_ratio:
            tags.append("convolve intensity")
        if self.exact_temperature:
            tags.append("exact temperature")
        return ", ".join(tags) if tags else "baseline (as in hires)"


# ----------------------------------------------------------------------
# conversion between intensity and surface density
# ----------------------------------------------------------------------
def _factor(wave_um):
    """N(H2) = factor * I[MJy/sr] / B_nu(T), matching 'operate'."""
    return 1.0e-17 / (kappa_abs(wave_um) * DUST_TO_GAS * MU_H2 * AMU)


def surfden_from_intensity(image_mjysr, temp_k, wave_um):
    t = np.clip(temp_k, T_CLAMP[0], T_CLAMP[1])
    return _factor(wave_um) * image_mjysr / planck(t, wave_um)


def intensity_from_surfden(sigma, temp_k, wave_um):
    t = np.clip(temp_k, T_CLAMP[0], T_CLAMP[1])
    return sigma * planck(t, wave_um) / _factor(wave_um)


# ----------------------------------------------------------------------
# pixel-by-pixel spectral fit, exploiting the circular symmetry of the model
# ----------------------------------------------------------------------
def cross_convolve(images, pix_as=PIX_ARCSEC):
    """Every waveband convolved to the beam of every longer waveband."""
    conv = {}
    for t in range(len(WAVES)):
        for j in range(t + 1):
            f = np.sqrt(max(BEAMS[t] ** 2 - BEAMS[j] ** 2, 0.0))
            conv[(t, j)] = images[j] if f <= 0 else gaussian_convolve(images[j], f, pix_as)
    return conv


def _unique_radius_index(n):
    """Index map that identifies pixels at the same distance from the center."""
    c = (n - 1) // 2
    y, x = np.mgrid[0:n, 0:n]
    d2 = (x - c) ** 2 + (y - c) ** 2
    vals, inv = np.unique(d2, return_inverse=True)
    return inv.reshape(n, n), len(vals)


def fit_temperature(images, bands, opt, radius_index=None):
    """Fit the optically thin modified blackbody pixel by pixel.

    images : (n_bands, n, n) intensities in MJy sr^-1, all at one resolution
    bands  : indices into WAVES of the wavebands supplied
    Returns the fitted dust temperature in K and the fitted surface density in
    cm^-2, both of shape (n, n).
    """
    n = images.shape[-1]
    if radius_index is None:
        radius_index, nuni = _unique_radius_index(n)
    flat = np.empty((len(bands), radius_index.max() + 1))
    for k in range(len(bands)):
        f = np.zeros(radius_index.max() + 1)
        np.put(f, radius_index.ravel(), images[k].ravel())
        flat[k] = f

    lam = np.array([WAVES[b] for b in bands])
    rel = opt.rel_error
    rel = (np.array([rel[b] for b in bands])[:, None]
           if np.ndim(rel) > 0 else float(rel))
    if opt.weighting == "relative":
        sig = rel * flat
    else:
        noise = (np.array([opt.noise_mjysr[b] for b in bands])[:, None]
                 if opt.noise_mjysr is not None else 0.0)
        sig = np.sqrt((rel * flat) ** 2 + noise ** 2)
    w = 1.0 / np.maximum(sig, 1e-30) ** 2

    tgrid = np.geomspace(3.0, 60.0, 500)
    best_chi2 = np.full(flat.shape[1], np.inf)
    best_t = np.full(flat.shape[1], np.nan)
    best_a = np.zeros(flat.shape[1])
    for t in tgrid:
        m = np.array([planck(t, l) / _factor(l) for l in lam])[:, None]
        num = (w * flat * m).sum(axis=0)
        den = (w * m * m).sum(axis=0)
        a = num / np.maximum(den, 1e-300)
        chi2 = (w * (flat - a * m) ** 2).sum(axis=0)
        upd = chi2 < best_chi2
        best_chi2[upd] = chi2[upd]
        best_t[upd] = t
        best_a[upd] = a[upd]
    # parabolic refinement on the logarithm of the temperature
    dl = np.log(tgrid[1] / tgrid[0])
    for _ in range(3):
        tl = best_t * np.exp(-dl)
        th = best_t * np.exp(dl)
        c = []
        for t in (tl, best_t, th):
            m = np.array([planck(t, l) / _factor(l) for l in lam])
            num = (w * flat * m).sum(axis=0)
            den = (w * m * m).sum(axis=0)
            a = num / np.maximum(den, 1e-300)
            c.append((w * (flat - a * m) ** 2).sum(axis=0))
        c0, c1, c2 = c
        denom = (c0 - 2 * c1 + c2)
        step = np.where(np.abs(denom) > 1e-300, 0.5 * (c0 - c2) / denom, 0.0)
        best_t = best_t * np.exp(np.clip(step, -1.0, 1.0) * dl)
        dl *= 0.5
    m = np.array([planck(best_t, l) / _factor(l) for l in lam])
    best_a = (w * flat * m).sum(axis=0) / np.maximum((w * m * m).sum(axis=0), 1e-300)

    return best_t[radius_index], best_a[radius_index]


# ----------------------------------------------------------------------
# circular median and mean filtering, as in 'imgstat'
# ----------------------------------------------------------------------
def _circular_footprint(radius_pix):
    r = int(radius_pix)
    if r < 1:
        return None
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= r * r


def filter_temperature(temp, beam_as, pix_as=PIX_ARCSEC):
    """Median filter then mean filter over a circular window whose radius in
    pixels is one half of the beam full width at half maximum, truncated to an
    integer, exactly as the HIRES script calls 'imgstat'."""
    r = int(beam_as / (2.0 * pix_as))
    fp = _circular_footprint(r)
    if fp is None:
        return temp.copy()
    out = ndimage.median_filter(temp, footprint=fp, mode="nearest")
    k = fp.astype(float)
    k /= k.sum()
    return ndimage.convolve(out, k, mode="nearest")


# ----------------------------------------------------------------------
# the reconstruction
# ----------------------------------------------------------------------
def run_hires(images, opt=None, sigma_true=None, pix_as=PIX_ARCSEC,
              return_all=False, temps_override=None):
    """Reconstruct the surface density from the six intensity maps.

    images     : (6, n, n) intensities in MJy sr^-1, each at its own beam
    sigma_true : true surface density map, required only when the switch
                 exact_temperature is set
    Returns the reconstructed surface density in cm^-2 at 13.5 arcsec, and,
    if return_all is set, a dictionary with the intermediate products.
    """
    opt = opt or Options()
    n = images.shape[-1]
    ridx, _ = _unique_radius_index(n)

    # ---- cross-convolution to each waveband resolution ----
    conv = {}
    for t in range(len(WAVES)):
        for j in range(t + 1):
            f = np.sqrt(max(BEAMS[t] ** 2 - BEAMS[j] ** 2, 0.0))
            conv[(t, j)] = images[j] if f <= 0 else gaussian_convolve(images[j], f, pix_as)

    # ---- temperature maps ----
    temps = {}
    if temps_override is not None:
        temps = {k: v.copy() for k, v in temps_override.items()}
        if opt.smooth_temperature:
            for t in TEMP_INDICES:
                temps[t] = filter_temperature(temps[t], BEAMS[t], pix_as)
    for t in ([] if temps_override is not None else TEMP_INDICES):
        bands = TEMP_BANDS[t]
        if opt.exact_temperature:
            if sigma_true is None:
                raise ValueError("exact_temperature requires sigma_true")
            # temperature that reproduces the true surface density at this
            # resolution from the 250 micron intensity
            sig_t = gaussian_convolve(sigma_true, BEAMS[t], pix_as)
            temps[t] = _invert_planck(conv[(t, 3)], sig_t, WAVES[3])
        else:
            stack = np.array([conv[(t, b)] for b in bands])
            temps[t], _ = fit_temperature(stack, bands, opt, ridx)
        if opt.smooth_temperature:
            temps[t] = filter_temperature(temps[t], BEAMS[t], pix_as)

    # ---- base map from the 500 micron intensity ----
    hires = surfden_from_intensity(conv[(5, 5)], temps[5], WAVES[5])
    out = {"base": hires.copy(), "temps": temps, "levels": {}}

    # ---- accumulate the differential terms ----
    for i in range(len(WAVES) - 1, -1, -1):
        terms = []
        for j in range(len(WAVES) - 1, i - 1, -1):
            if j < NWAVEFIN:
                continue
            if opt.convolve_ratio:
                mono = surfden_from_intensity(images[i], temps[j], WAVES[i])
                if i == len(WAVES) - 1:
                    low = mono
                else:
                    f = np.sqrt(max(BEAMS[i + 1] ** 2 - BEAMS[i] ** 2, 0.0))
                    low = gaussian_convolve(mono, f, pix_as)
                delta = mono - low
            else:
                if i == len(WAVES) - 1:
                    dint = np.zeros_like(images[i])
                else:
                    f = np.sqrt(max(BEAMS[i + 1] ** 2 - BEAMS[i] ** 2, 0.0))
                    dint = images[i] - gaussian_convolve(images[i], f, pix_as)
                delta = surfden_from_intensity(dint, temps[j], WAVES[i])
            if opt.clip and i < NWAVEFIN:
                delta = np.maximum(delta, 0.0)
            terms.append(delta)
        if terms:
            stack = np.array(terms)
            add = stack.max(axis=0) if opt.combine == "max" else stack.mean(axis=0)
        else:
            add = np.zeros_like(hires)
        hires = hires + add
        out["levels"][BEAMS[i]] = hires.copy()

    if return_all:
        return out["levels"][13.5], out
    return out["levels"][13.5]


def _invert_planck(image_mjysr, sigma, wave_um):
    """Dust temperature satisfying I = sigma * B(T) * kappa * eta * mu * m_H,
    obtained by analytic inversion of the Planck function."""
    from hires_forward import H_PLANCK, K_BOLTZ, C_LIGHT
    b = 1.0e-17 * image_mjysr / np.maximum(sigma, 1e-30) \
        / (kappa_abs(wave_um) * DUST_TO_GAS * MU_H2 * AMU)
    nu = C_LIGHT / (wave_um * 1.0e-4)
    pref = 2.0 * H_PLANCK * nu ** 3 / C_LIGHT ** 2
    x = np.log1p(pref / np.maximum(b, 1e-300))
    return np.clip(H_PLANCK * nu / (K_BOLTZ * np.maximum(x, 1e-12)),
                   T_CLAMP[0], T_CLAMP[1])
