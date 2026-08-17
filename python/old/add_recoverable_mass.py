#!/usr/bin/env python3
"""
add_recoverable_mass.py -- add environment-dependent recoverable-mass columns
to the BE model grid catalog.

For each grid model, the recoverable (getsf-detectable) mass fraction is
estimated from the model's own convolved surface-density stamp using an
empirical getsf detection floor calibrated on the real Aquila ok-SED catalog:

    N_src_min(N_back) = 1.97e21 * (N_back / 1e22)^1.34   [cm^-2]

(the faintest detected source peak column PEAK^SRC03 vs the local background
column PEAK^BGF03, lower-envelope fit over 136 reliable sources).

Two flux-loss mechanisms are applied to the (background-subtracted) surfdens
stamp of each model, embedded at its own SD_emb:

  1. Outskirt truncation: the detectable footprint is limited to the radius
     where the source column drops to the empirical detection floor
     N_src_min(N_back) for the local background column.

  2. Background pedestal subtraction: getsf interpolates and subtracts a
     background pinned at the footprint rim; we subtract that pedestal
     (the column value at the truncation radius) over the footprint and
     clip negatives.

The recoverable fraction is the residual flux divided by the full stamp flux.
Output columns:  M_SED3bs_rec  = M_SED3bs  * frac
                 M_SED3bsl_rec = M_SED3bsl * frac
                 frac_rec      = the fraction itself (diagnostic)

Stamp location (matches the RADMC-3D output tree):
    <stamp_root>/cSD_{i:02d}/M_{j:02d}/{k:02d}/nc.surfdens.bs.r13p5x0.rs3p0as.fits

Usage:
    python add_recoverable_mass.py \
        --grid bes_model_params_catalog \
        --stamp-root /Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260616_RT_BES_radmc3d \
        --out bes_model_params_catalog_rec \
        --col-msed3bs 20 --col-msed3bsl 23

Requires numpy and astropy.
"""

import argparse
import os
import numpy as np
from astropy.io import fits

PIX_ARCSEC = 3.0        # all stamps resampled to a common 3" grid before getsf

# --- Empirical getsf detection floor, calibrated on the real Aquila "ok" (matched
#     reliable-SED) source catalog: the faintest detected source peak column
#     (PEAK^SRC03) as a function of the local background column (PEAK^BGF03).
#     Lower-envelope (5th-percentile) fit over 136 sources:
#         N_src_min(N_back) = FLOOR_NORM * (N_back / FLOOR_REF)^FLOOR_INDEX
#     i.e. a contrast floor of ~0.20 at 1e22, ~0.13 at low column, ~0.43 at high column.
#     This REPLACES the old 5*N_rms Konyves-cirrus cut, which used getsf's
#     single-scale combination sigma (a per-decomposed-scale quantity), not the
#     detection floor in the original combined image, and was far too harsh --
#     it drove frac_rec to zero above N_back ~ 1.2e22, where real cores are in
#     fact routinely detected.
FLOOR_NORM  = 1.971e21     # cm^-2  (source-peak floor at N_back = FLOOR_REF)
FLOOR_REF   = 1.0e22       # cm^-2
FLOOR_INDEX = 1.342

STAMP_NAME = 'nc.surfdens.bs.r13p5x0.rs3p0as.fits'


def stamp_path(stamp_root, i, j, k):
    """Path in the RADMC-3D tree: cSD_{i:02d}/M_{j:02d}/{k:02d}/<STAMP_NAME>."""
    return os.path.join(stamp_root,
                        f'cSD_{i:02d}', f'M_{j:02d}', f'{k:02d}',
                        STAMP_NAME)


def N_detect_floor(N_back):
    """Empirical getsf source-peak detection floor (cm^-2) at background column
    N_back, from the real Aquila ok-SED catalog (see constants above)."""
    return FLOOR_NORM * (N_back / FLOOR_REF) ** FLOOR_INDEX


def recoverable_fraction(stamp, sd_emb, pix_arcsec=PIX_ARCSEC,
                         eta=3.0):
    """Fraction of the model's stamp flux that getsf recovers against a
    fluctuating background at column sd_emb.

    Replicates getsf's footprint behaviour (Men'shchikov 2021, Sect. 3.4.6):
      * the measurement footprint is set at eta * H_n, where H_n is the
        half-maximum radius of the (convolved) source -- NOT a 5-sigma cut;
      * getsf expands the footprint for resolved power-law sources to reduce
        the residual background pedestal, but the expansion is detection-
        limited: it cannot extend past the radius where the source profile
        sinks to the empirical detection floor N_detect_floor(N_back).
      * whatever pedestal remains at the footprint rim is subtracted.

    Effective footprint radius = min(eta * H_n, r_noise), where r_noise is
    where the azimuthal profile drops to N_detect_floor(sd_emb). In low-noise
    (low column) environments the footprint reaches eta*H_n and recovers most
    of the flux; in high-noise (high column) environments it is truncated
    earlier and more flux is lost -- giving the environment dependence.

    stamp : 2D array, background-subtracted convolved surface density (cm^-2),
            source centered, zero far from the core.
    Returns (frac, r_foot, pedestal, fwhm_rec, conc_footfwhm,
             conc_peakmean, conc_slope).
    """
    ny, nx = stamp.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(yy - cy, xx - cx) * pix_arcsec        # radius map (arcsec)

    # azimuthally-averaged radial profile (1-pixel-wide rings)
    rmax = int(np.floor(rr.max() / pix_arcsec))
    r_edges = (np.arange(rmax + 2) - 0.5) * pix_arcsec
    prof = np.zeros(rmax + 1)
    for ir in range(rmax + 1):
        m = (rr >= r_edges[ir]) & (rr < r_edges[ir + 1])
        prof[ir] = stamp[m].mean() if m.any() else 0.0
    r_centers = np.arange(rmax + 1) * pix_arcsec

    peak = prof[0]
    if peak <= 0:
        return 0.0, 0.0, 0.0, np.nan, np.nan, np.nan, np.nan

    # half-maximum radius H_n
    half = peak / 2.0
    below_half = np.where(prof < half)[0]
    H_n = r_centers[below_half[0]] if below_half.size else r_centers[-1]

    # getsf target footprint: eta * H_n
    r_eta = eta * H_n

    # detection-limited radius: where the (bg-subtracted) source profile sinks
    # to the empirical detection floor for this background column.  The footprint
    # cannot be expanded past where the source is no longer detectable.
    thr = N_detect_floor(sd_emb)
    below_thr = np.where(prof < thr)[0]
    r_noise = r_centers[below_thr[0]] if below_thr.size else r_centers[-1]

    # effective footprint = the smaller of the two (expansion is noise-limited)
    r_foot = min(r_eta, r_noise)

    # pedestal = profile value at the footprint rim (what getsf subtracts)
    irim = int(np.searchsorted(r_centers, r_foot))
    irim = min(irim, len(prof) - 1)
    pedestal = prof[irim]

    # full flux (all positive stamp signal)
    total = stamp[stamp > 0].sum()
    if total <= 0:
        return 0.0, r_foot, pedestal, np.nan, np.nan, np.nan, np.nan

    # recovered: within footprint, minus pedestal, clipped at 0
    foot = rr <= r_foot
    resid = np.clip(stamp[foot] - pedestal, 0.0, None)
    recovered = resid.sum()

    # --- concentration diagnostics on the RECOVERED (truncated, pedestal-
    #     subtracted) profile.  These distinguish a genuine compact core from
    #     a huge diffuse core whose only central part survived: the latter has
    #     a steep, over-concentrated recovered profile. ---
    prof_sub = np.clip(prof - pedestal, 0.0, None)
    prof_sub[r_centers > r_foot] = 0.0
    pk = prof_sub[0]
    if pk > 0 and r_foot > 0:
        # recovered FWHM (half-max diameter of the truncated profile)
        halfp = pk / 2.0
        bh = np.where(prof_sub < halfp)[0]
        r_half = r_centers[bh[0]] if bh.size else r_foot
        fwhm_rec = 2.0 * r_half

        # (1) footprint-to-FWHM ratio: large for a truncated giant (footprint
        #     much larger than the surviving central half-max region), ~a few
        #     for a genuine compact core.
        conc_footfwhm = (2.0 * r_foot) / fwhm_rec if fwhm_rec > 0 else np.nan

        # (2) peak-to-mean ratio inside the footprint: high for a steep
        #     (over-concentrated) recovered profile.
        mean_in = resid[resid > 0].mean() if np.any(resid > 0) else np.nan
        conc_peakmean = pk / mean_in if (mean_in and mean_in > 0) else np.nan

        # (3) logarithmic profile slope near the footprint rim: steeper
        #     (more negative) for a truncated giant.
        inner = max(1, int(0.5 * r_foot / pix_arcsec))
        outer = max(inner + 1, int(0.9 * r_foot / pix_arcsec))
        outer = min(outer, len(prof_sub) - 1)
        if (prof_sub[inner] > 0 and prof_sub[outer] > 0
                and r_centers[outer] > r_centers[inner]):
            conc_slope = ((np.log(prof_sub[outer]) - np.log(prof_sub[inner])) /
                          (np.log(r_centers[outer]) - np.log(r_centers[inner])))
        else:
            conc_slope = np.nan
    else:
        fwhm_rec = conc_footfwhm = conc_peakmean = conc_slope = np.nan

    return (recovered / total, r_foot, pedestal,
            fwhm_rec, conc_footfwhm, conc_peakmean, conc_slope)


def load_grid_lines(path):
    """Return (header_lines, data_lines, col_index_map). Assumes whitespace
    columns; header lines start with '#'. Locates i/j/k, M_SED3bs, M_SED3bsl
    by position from the first data row using the known layout."""
    header, data = [], []
    for line in open(path):
        if line.lstrip().startswith('#'):
            header.append(line.rstrip('\n'))
        elif line.strip():
            data.append(line.rstrip('\n'))
    return header, data


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--grid', required=True, help='input grid catalog')
    ap.add_argument('--stamp-root', required=True,
                    help='root of RADMC-3D tree (contains cSD_NN/M_NN/KK/)')
    ap.add_argument('--out', required=True, help='output catalog with new columns')
    # column indices (0-based) in the data rows; defaults match the standard
    # bes_model_params_catalog layout:
    #   0:n 1:i 2:j 3:k 4:SD_emb ... 7:M_BE ... plus M_SED3bs / M_SED3bsl
    ap.add_argument('--col-i',  type=int, default=1)
    ap.add_argument('--col-j',  type=int, default=2)
    ap.add_argument('--col-k',  type=int, default=3)
    ap.add_argument('--col-sd', type=int, default=4)
    ap.add_argument('--col-msed3bs',  type=int, required=True,
                    help='0-based column index of M_SED3bs in the data rows')
    ap.add_argument('--col-msed3bsl', type=int, required=True,
                    help='0-based column index of M_SED3bsl in the data rows')
    args = ap.parse_args()

    header, data = load_grid_lines(args.grid)

    out_rows = []
    n_done = n_miss = 0
    for line in data:
        p = line.split()
        i = int(float(p[args.col_i]))
        j = int(float(p[args.col_j]))
        k = int(float(p[args.col_k]))
        sd = float(p[args.col_sd])
        m3bs  = float(p[args.col_msed3bs])
        m3bsl = float(p[args.col_msed3bsl])

        fpath = stamp_path(args.stamp_root, i, j, k)
        fname = os.path.relpath(fpath, args.stamp_root)
        if os.path.exists(fpath):
            stamp = fits.getdata(fpath).astype(float)
            # all stamps are resampled to a common PIX_ARCSEC (3") grid before getsf,
            # so the fixed pixel scale is correct here (the variable RT pixel matters
            # only upstream, at the imaging stage).
            frac, r_foot, ped, fwhm_rec, cff, cpm, cslp = recoverable_fraction(stamp, sd)
            n_done += 1
        else:
            frac, r_foot, ped, fwhm_rec, cff, cpm, cslp = (np.nan,)*7
            n_miss += 1
            print(f'  MISSING stamp: {fname}')

        rec_bs  = m3bs  * frac
        rec_bsl = m3bsl * frac
        out_rows.append(f'{line}  {rec_bs:12.5e}  {rec_bsl:12.5e}  {frac:8.4f}'
                        f'  {fwhm_rec:9.3f}  {cff:8.3f}  {cpm:9.3f}  {cslp:8.3f}')

    # append column description to header
    header.append('#')
    header.append('# Added by add_recoverable_mass.py:')
    header.append('#   M_SED3bs_rec  = M_SED3bs  * frac_rec   (Msun)')
    header.append('#   M_SED3bsl_rec = M_SED3bsl * frac_rec   (Msun)')
    header.append('#   frac_rec      = recoverable flux fraction '
                  '(Konyves noise-limited footprint + rim-pedestal subtraction)')
    header.append('#   FWHM_rec      = recovered half-max size (arcsec)')
    header.append('#   conc_footfwhm = footprint diameter / FWHM_rec')
    header.append('#   conc_peakmean = peak / mean intensity in footprint')
    header.append('#   conc_slope    = log-log profile slope near rim')

    with open(args.out, 'w') as f:
        f.write('\n'.join(header) + '\n')
        f.write('\n'.join(out_rows) + '\n')

    print(f'Wrote {args.out}: {n_done} models processed, {n_miss} stamps missing')


if __name__ == '__main__':
    main()
