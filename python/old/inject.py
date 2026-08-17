"""
inject.py -- inject synthetic BE-sphere cores into real Herschel maps and
build the truth table for the injection-recovery test.

This file currently implements the FIELD side and PLACEMENT logic (loading the
maps/masks, finding valid positions for a model given its Sigma_cloud and
R_BE). The per-band stamping (map + 2C - I, masked to r<=R_BE) plugs in once
the model C/I/mask images are available; see add_stamp() stub at the end.

Conventions (confirmed): all maps are on a common pixel grid, 3"/pixel, fully
aligned; placement uses pixel coordinates only (no WCS). Sigma in H2/cm^2,
band images in MJy/sr. Coverage masks (.omask) are 0/1.
"""

import numpy as np
from astropy.io import fits
from scipy.ndimage import uniform_filter, binary_erosion, distance_transform_edt

PIX_ARCSEC = 3.0
BANDS = ['000', '070', '160', '250', '350', '500']


# ----------------------------------------------------------------------
# field loading
# ----------------------------------------------------------------------
def _read(path):
    d = fits.open(path)[0].data
    return np.asarray(d, float)


def load_field(img_paths, omask_paths, sigma_path, sigma_omask_path):
    """img_paths/omask_paths: dict band->path. Returns a field dict."""
    img = {b: _read(img_paths[b]) for b in BANDS}
    omask = {b: _read(omask_paths[b]) > 0 for b in BANDS}
    sigma = _read(sigma_path)
    sig_omask = _read(sigma_omask_path) > 0
    shape = sigma.shape
    for b in BANDS:
        assert img[b].shape == shape, f"band {b} grid mismatch"
    # combined coverage: every band + Sigma covered, and Sigma finite/positive
    coverage = sig_omask & np.isfinite(sigma) & (sigma > 0)
    for b in BANDS:
        coverage &= omask[b]
    return dict(img=img, omask=omask, sigma=sigma, coverage=coverage, shape=shape)


# ----------------------------------------------------------------------
# placement
# ----------------------------------------------------------------------
def local_sigma(sigma, R_pix):
    """Mean Sigma over a box ~ the core size, as the local-cloud estimate."""
    n = max(3, int(round(2 * R_pix + 1)))
    return uniform_filter(sigma, size=n, mode='nearest')


def valid_placement_mask(field, sigma_cloud, R_BE_arcsec,
                         sqrt2_window=True, forbidden=None,
                         stamp_halfwidth=None):
    """Boolean mask of pixels where the model's R_BE disk may be centered.

    A position is valid if the whole stamp window fits inside coverage (and
    clear of `forbidden`, e.g. source footprints), away from edges, and the
    local Sigma matches sigma_cloud within a factor sqrt(2) -- but never
    *below* sigma_cloud, so the corrector always has grid coverage.

    The Sigma window is therefore one-sided:
        sigma_cloud <= local_Sigma < sigma_cloud * sqrt(2)

    stamp_halfwidth : half-width of the model stamp in pixels (= (n-1)//2,
        where n is the stamp array size).  The edge margin is the larger of
        R_BE and stamp_halfwidth, so that the full stamp always fits inside
        the field.  Pass this whenever stamp size > R_BE (which is always
        true for 73x73 stamps at typical R_BE ~ 6 pix).
    """
    R_pix = R_BE_arcsec / PIX_ARCSEC
    from scipy.ndimage import binary_erosion as _be

    # --- field-edge erosion: stamp must fit entirely inside coverage ---
    edge_margin = int(np.ceil(stamp_halfwidth if stamp_halfwidth else R_pix))
    edge_struct = np.ones((2 * edge_margin + 1, 2 * edge_margin + 1), bool)
    fit = _be(field['coverage'], structure=edge_struct, border_value=0)

    # --- footprint avoidance: only the core disk (R_BE) must clear forbidden ---
    if forbidden is not None:
        r = int(np.ceil(R_pix))
        yy, xx = np.ogrid[-r:r + 1, -r:r + 1]
        foot_struct = (xx ** 2 + yy ** 2) <= R_pix ** 2
        fit &= _be(~forbidden, structure=foot_struct, border_value=0)
    # Sigma match: one-sided window [sigma_cloud, sigma_cloud * sqrt(2))
    # Lower bound >= sigma_cloud ensures the corrector interpolant has coverage.
    # Upper bound < sigma_cloud * sqrt(2) keeps the environment representative.
    # Averaging radius = R_BE: this is the placement criterion (does this
    # location have the right ambient column?).  The truth-table records a
    # larger 2*R_BE average as the corrector input (cloud scale), computed in
    # run_inject.py after placement.
    sloc = local_sigma(field['sigma'], R_pix)
    if sqrt2_window:
        match = (sloc >= sigma_cloud) & (sloc < sigma_cloud * np.sqrt(2.0))
    else:
        match = np.ones_like(fit)
    return fit & match


def choose_positions(valid, n, min_sep_pix, rng=None, max_tries=200000):
    """Greedily pick up to n centers from `valid` with a minimum separation."""
    rng = np.random.default_rng() if rng is None else rng
    ys, xs = np.where(valid)
    if len(ys) == 0:
        return []
    order = rng.permutation(len(ys))
    chosen = []
    s2 = min_sep_pix * min_sep_pix
    for idx in order[:max_tries]:
        y, x = ys[idx], xs[idx]
        ok = all((y - cy) ** 2 + (x - cx) ** 2 >= s2 for cy, cx in chosen)
        if ok:
            chosen.append((int(y), int(x)))
            if len(chosen) >= n:
                break
    return chosen


# ----------------------------------------------------------------------
# stamping (plugs in when model C/I/mask images are available)
# ----------------------------------------------------------------------
def add_stamp(field_img_band, center, I_bs, mask):
    """Add one model into one band's map: map += I_bs * mask, centered.

    The injected signal is the interpolated-background-subtracted image I (the
    'bs' image): at a site whose observed background matches the model cloud
    level, (sky with core) - (sky without core) = raw - const = I. This already
    contains the displaced (cratered) cloud; no separate crater digging is
    needed (equivalently, dig the depression C-I and add the crater core C,
    which reduces to +I). Only the bs image is needed; the crater image C
    cancels. `mask` is the per-band core mask (0/1, zero outside R_BE).
    """
    stamp = np.nan_to_num(I_bs * mask)
    h = (I_bs.shape[0] - 1) // 2
    cy, cx = center
    field_img_band[cy - h:cy + h + 1, cx - h:cx + h + 1] += stamp
    return field_img_band


# ----------------------------------------------------------------------
# controlled-background injection (flatten cloud fluctuations under source)
# ----------------------------------------------------------------------
def planar_background(img, center, support_pix, pad_pix=6):
    """Low-order (planar) interpolated background over the stamp window,
    fitted to an annulus just outside the source support. Returns a window
    the size of (2*support_pix+1) approximated to the stamp; here we evaluate
    over the full model-stamp window passed by the caller via `win_half`."""
    raise NotImplementedError  # see inject_model, which fits inline


def inject_model(field, center, Ibs, mask, flatten=True, pad_pix=8,
                 feather_pix=1.5, nave=1, nradmf=4):
    """Inject one model (all bands) at `center`. Returns dict band->new map.

    flatten=True : replace the cloud under the source footprint with getsf's
       own background interpolation (power-9 directional scheme + median filter,
       from getsf_bg.py) before adding I.  This is the correct choice for the
       injection test: the source sees the same background estimation that getsf
       will later apply when extracting it, so the only residual is the
       separation error of getsf itself -- not a mismatch between injection and
       extraction backgrounds.  Parameters nave and nradmf must match the getsf
       run (defaults nave=1, nradmf=4 reproduce getsf to 0.02% rms).
    flatten=False : keeps the raw cloud under the source (realism /
       fluctuation-budget test).
    Ibs, mask : dict band -> full model arrays (same odd npix, centered).
    """
    from getsf_bg import getsf_background
    out = {}
    cy, cx = center
    ny, nx = field['shape']
    for b in BANDS:
        I = Ibs[b]; m = mask[b] > 0
        n = I.shape[0]; h = (n - 1) // 2
        y0, y1 = cy - h, cy + h + 1
        x0, x1 = cx - h, cx + h + 1
        if y0 < 0 or y1 > ny or x0 < 0 or x1 > nx:
            raise ValueError(
                f'Stamp window [{y0}:{y1}, {x0}:{x1}] out of bounds '
                f'for field shape ({ny}, {nx}) at center ({cy}, {cx}). '
                f'Pass stamp_halfwidth={h} to valid_placement_mask.')
        sl = (slice(y0, y1), slice(x0, x1))
        win = field['img'][b][sl].astype(float).copy()
        if flatten:
            # replace cloud under the footprint using getsf's own background
            # scheme: power-9 directional interpolation + median filter.
            # win is the stamp-sized cutout; m is the footprint mask.
            bg = getsf_background(win, m, nave=nave, nradmf=nradmf)
            win = bg  # cloud under footprint replaced by interpolated background
        win = win + np.nan_to_num(I)             # add full stamp (mask only used for bg interpolation)
        nm = field['img'][b].astype(float).copy()
        nm[sl] = win
        out[b] = nm
    return out
