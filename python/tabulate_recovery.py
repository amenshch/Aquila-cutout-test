#!/usr/bin/env python3
"""
tabulate_recovery.py -- recovery quantities as a function of truncation radius,
for every model of the BE grid and every band.

Why this replaces the single-value approach
-------------------------------------------
add_recoverable_mass.py assigns ONE recovered fraction per model, obtained by
truncating the profile at a detection floor taken from an empirical relation
    N_src_min(N_back) = 1.971e21 * (N_back/1e22)**1.342
calibrated as an average over all places in Aquila with a given background.
Two measurements show this cannot represent an individual source:

  * at FIXED background column, the per-source noise (FXP_ERR) in a real getsf
    catalogue spans a factor 3.3-5.7 between the 16th and 84th percentiles;
  * the measured footprint-to-size ratio FOOA/FWHM spans 1.5-2.5, with 62% of
    sources outside the +-10% band around the eta = 2 assumed by the models.

Both distributions are collapsed to a single value by the current scheme, and
the error goes straight into frac.  Here we instead tabulate the recovery
quantities against the truncation radius itself, so that the matching variable
can be chosen afterwards -- everything below is a monotonic function of r, so a
single table supports matching on peak/noise, on the measured footprint, or on
the concentration, without recomputing anything.

No radiative transfer is repeated: this is post-processing of existing stamps.

What is tabulated
-----------------
For each node, band and truncation radius r:

    r, r/R_conv        truncation radius, and as a fraction of the full extent
    I_rim              profile value at r -- the EFFECTIVE detection threshold
    peak_over_rim      peak / I_rim; the model counterpart of an observed
                       peak-to-noise ratio, and the variable that locates a
                       real source on this curve
    frac               recovered flux fraction, normalised so that frac = 1 at
                       r = R_conv (the full footprint)
    fwhm_rec           FWHM of the truncated, rim-subtracted profile
    conc_footfwhm      2r / fwhm_rec
    conc_peakmean      peak / mean intensity inside r  (both rim-subtracted)
    conc_slope         logarithmic profile slope between 0.5r and 0.9r
    reliable           1 if frac >= --frac-min, else 0

I_rim is essential: a real source's r/R_conv is NOT measurable, because R_conv
is the full extent, which is precisely what is not seen.  peak/I_rim is what
has an observational counterpart.

Ranges
------
Model validation may use arbitrarily small radii, but heavily truncated sources
cannot be corrected robustly on real data -- below frac ~ 0.2 the leave-one-out
error exceeds 30% even for in-family models.  --r-min-frac therefore defaults to
0.05 (validation), while --frac-min defaults to 0.2 and only sets the
`reliable` flag; nothing is discarded here, the flag is applied downstream.

Output
------
A long-format table, one line per (node, band, radius), written separately from
the grid catalogue so that the catalogue stays readable and diffable.

Usage
-----
    python tabulate_recovery.py \
        --grid bes_model_grid_final3_catalog \
        --stamp-root /path/to/RT_BES_radmc3d \
        --out   bes_model_grid_final3_recovery_tables \
        --bands SD,160,250,350,500 --nr 15
"""

import argparse
import os

import numpy as np
from astropy.io import fits

PIX_ARCSEC = 3.0

# label, beam FWHM ["], stamp token
BANDS = {
    'SD':  (13.5, 'surfdens'),
    '070': (8.4,  '070um'),
    '100': (9.4,  '100um'),
    '160': (13.5, '160um'),
    '250': (18.2, '250um'),
    '350': (24.9, '350um'),
    '500': (36.3, '500um'),
}


def radial_profile(stamp, pix_arcsec=PIX_ARCSEC):
    """Azimuthally averaged profile and the radius of each pixel."""
    ny, nx = stamp.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(yy - cy, xx - cx) * pix_arcsec
    nb = int(np.floor(rr.max() / pix_arcsec))
    edges = (np.arange(nb + 2) - 0.5) * pix_arcsec
    prof = np.zeros(nb + 1)
    for i in range(nb + 1):
        m = (rr >= edges[i]) & (rr < edges[i + 1])
        prof[i] = stamp[m].mean() if m.any() else 0.0
    return rr, np.arange(nb + 1) * pix_arcsec, prof


def full_extent(r_centers, prof, rel=1e-3):
    """R_conv from the profile: fallback when no mask is available.

    Taken as the first radius where the profile drops below `rel` times the
    peak, or turns non-positive, whichever comes first.  This threshold is
    arbitrary; prefer extent_from_mask() when the source masks exist.
    """
    peak = prof[0]
    if peak <= 0:
        return np.nan
    bad = np.where((prof <= 0) | (prof < rel * peak))[0]
    return r_centers[bad[0]] if bad.size else r_centers[-1]


def extent_from_mask(path, pix_arcsec=PIX_ARCSEC):
    """R_conv from the source mask, as the equivalent-area radius.

    The masks (1|0) are built from the surface-density image convolved to each
    band's resolution, so they carry the model's own definition of the source
    footprint -- the same one over which the 'bs' background was interpolated.
    Using them removes the arbitrary threshold of full_extent() and guarantees
    that frac = 1 corresponds to exactly the footprint the model fluxes were
    measured in.

    The equivalent-area radius, sqrt(N_pix * pix^2 / pi), is used rather than
    the maximum radius, since it is insensitive to pixelation of the rim.  For
    a spherically symmetric model the two agree to within a pixel.
    """
    with fits.open(path, memmap=True) as hd:
        m = np.squeeze(hd[0].data)
    n = int(np.count_nonzero(np.nan_to_num(m) > 0))
    if n <= 0:
        return np.nan
    return np.sqrt(n * pix_arcsec ** 2 / np.pi)


def frac_curve(rr, stamp, r_grid):
    """Fast frac(r) on a grid of radii, from cumulative sums.

    For a monotonically declining profile every pixel inside r exceeds the rim
    value, so the clip in quantities_at() is a no-op and

        F(r) = C(r) - I_rim(r) * A(r)

    with C the cumulative flux and A the enclosed area.  This lets the curve be
    built once, cheaply, and inverted to place the tabulated radii at chosen
    values of frac rather than at chosen fractions of R_conv.
    """
    rf = rr.ravel(); sf = stamp.ravel()
    o = np.argsort(rf); rs, ss = rf[o], sf[o]
    csum = np.cumsum(ss)
    idx = np.searchsorted(rs, r_grid, side='right')
    idx = np.clip(idx, 1, len(rs))
    C = csum[idx - 1]
    A = idx.astype(float)
    I_rim = np.interp(r_grid, rs, ss)
    F = np.clip(C - I_rim * A, 0.0, None)
    return F


def radii_for_fracs(rr, stamp, R_conv, targets, beam):
    """Radii at which frac reaches each target value.

    Sampling in frac rather than in r/R_conv keeps every tabulated row
    informative regardless of how far the source mask extends beyond the
    flux-bearing region (a 1e16 clip typically reaches 1.25-1.47x past the
    99%-flux radius, most of all for compact sources at coarse resolution).
    """
    fine = np.geomspace(max(beam / 4.0, 1e-3), R_conv, 400)
    F = frac_curve(rr, stamp, fine)
    if F[-1] <= 0:
        return np.array([])
    f = F / F[-1]
    f = np.maximum.accumulate(f)          # enforce monotonicity
    out = []
    for t in targets:
        if t >= f[-1]:
            out.append(R_conv)
        else:
            out.append(float(np.interp(t, f, fine)))
    return np.array(out)


def half_max_radius(r_centers, prof_sub, r):
    """Interpolated half-maximum radius of a rim-subtracted, truncated profile.

    getsf obtains AFWHM/BFWHM by interpolating the source intensity
    distribution to half the peak, so the model must do the same: taking the
    first bin centre below half quantises the result to the pixel size and
    makes the model FWHM incomparable with the measured one.
    """
    pk = prof_sub[0]
    if pk <= 0:
        return np.nan
    bh = np.where(prof_sub < pk / 2.0)[0]
    if not bh.size:
        return r
    a = bh[0]
    if a == 0:
        return r_centers[0]
    y0, y1 = prof_sub[a - 1], prof_sub[a]
    r0, r1 = r_centers[a - 1], r_centers[a]
    if y0 <= y1:
        return r0
    return r0 + (r1 - r0) * (y0 - pk / 2.0) / (y0 - y1)


def cfoot_curve(r_centers, prof, r_grid):
    """Footprint-to-FWHM ratio 2r/FWHM(r) over a grid of truncation radii.

    This is the quantity the matching should use: it is measurable per source
    as FOOA/AFWHM, it is dimensionless, and it constrains how much flux lies
    outside the footprint far better than the size alone.  For a Gaussian,
    2r/FWHM = 1.82 corresponds to frac = 0.90 and 2.58 to frac = 0.99; the
    seven getsf catalogues all sit at 1.8-2.0, so that is the range the tables
    must resolve.
    """
    out = np.empty(len(r_grid))
    for t, r in enumerate(r_grid):
        irim = int(np.clip(np.searchsorted(r_centers, r), 0, len(prof) - 1))
        ps = np.clip(prof - prof[irim], 0.0, None)
        ps[r_centers > r] = 0.0
        rh = half_max_radius(r_centers, ps, r)
        out[t] = (2.0 * r) / (2.0 * rh) if (np.isfinite(rh) and rh > 0) else np.nan
    return out


def radii_for_cfoot(r_centers, prof, R_conv, targets, beam):
    """Radii at which the footprint/FWHM ratio reaches each target value."""
    fine = np.geomspace(max(beam / 4.0, 1e-3), R_conv, 250)
    cf = cfoot_curve(r_centers, prof, fine)
    good = np.isfinite(cf)
    if good.sum() < 5:
        return np.array([])
    f, x = cf[good], fine[good]
    o = np.argsort(f); f, x = f[o], x[o]
    f = np.maximum.accumulate(f)
    out = []
    for t in targets:
        if t <= f[0]:
            out.append(x[0])
        elif t >= f[-1]:
            out.append(x[-1])
        else:
            out.append(float(np.interp(t, f, x)))
    return np.array(out)


def quantities_at(rr, r_centers, prof, stamp, r, total):
    """Recovery quantities for a truncation at radius r."""
    irim = int(np.clip(np.searchsorted(r_centers, r), 0, len(prof) - 1))
    i_rim = prof[irim]
    peak = prof[0]

    foot = rr <= r
    resid = np.clip(stamp[foot] - i_rim, 0.0, None)
    recovered = resid.sum()
    frac = recovered / total if total > 0 else np.nan

    prof_sub = np.clip(prof - i_rim, 0.0, None)
    prof_sub[r_centers > r] = 0.0
    pk = prof_sub[0]

    fwhm_rec = conc_footfwhm = conc_peakmean = conc_slope = np.nan
    if pk > 0 and r > 0:
        # half-maximum radius by linear interpolation between the bracketing
        # bins.  Taking the first bin centre below half quantises fwhm_rec to
        # the pixel size, which makes it incomparable with getsf's AFWHM /
        # BFWHM -- getsf interpolates the intensity distribution to half peak,
        # so the model must do the same for the two to be matched.
        r_half = half_max_radius(r_centers, prof_sub, r)
        fwhm_rec = 2.0 * r_half if np.isfinite(r_half) else np.nan
        if fwhm_rec > 0:
            conc_footfwhm = (2.0 * r) / fwhm_rec
        pos = resid[resid > 0]
        if pos.size:
            conc_peakmean = pk / pos.mean()
        i_in = max(1, int(0.5 * r / PIX_ARCSEC))
        i_out = min(max(i_in + 1, int(0.9 * r / PIX_ARCSEC)), len(prof_sub) - 1)
        if (prof_sub[i_in] > 0 and prof_sub[i_out] > 0
                and r_centers[i_out] > r_centers[i_in]):
            conc_slope = ((np.log(prof_sub[i_out]) - np.log(prof_sub[i_in])) /
                          (np.log(r_centers[i_out]) - np.log(r_centers[i_in])))

    peak_over_rim = peak / i_rim if i_rim > 0 else np.nan
    return (i_rim, peak_over_rim, frac, fwhm_rec, conc_footfwhm,
            conc_peakmean, conc_slope)


def stamp_path(root, i, j, k, token, beam, pattern):
    res = ('r%.1f' % beam).replace('.', 'p')
    return os.path.join(root, 'cSD_%02d' % i, 'M_%02d' % j, '%02d' % k,
                        pattern.format(band=token, res=res))


def mask_path(root, i, j, k, beam, pattern):
    res = ('r%.1f' % beam).replace('.', 'p')
    return os.path.join(root, 'cSD_%02d' % i, 'M_%02d' % j, '%02d' % k,
                        pattern.format(res=res))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid', required=True)
    ap.add_argument('--stamp-root', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--stamp-pattern', default='nc.{band}.bs.{res}x0.rs3p0as.fits')
    ap.add_argument('--mask-pattern',
                    default='nc.surfdens.bsl.{res}x0.rs3p0as.mask.fits',
                    help='source mask; defines R_conv exactly. "" disables')
    ap.add_argument('--bands', default='SD,070,100,160,250,350,500')
    ap.add_argument('--nr', type=int, default=15,
                    help='number of truncation radii per band')
    ap.add_argument('--frac-lo', type=float, default=0.02,
                    help='smallest tabulated recovered fraction')
    ap.add_argument('--r-min-frac', type=float, default=0.05,
                    help='floor on r/R_conv, applied after the frac sampling')
    ap.add_argument('--frac-min', type=float, default=0.20,
                    help='frac below which a source is flagged unreliable')
    ap.add_argument('--sample', default='both', choices=['frac', 'cfoot', 'both'],
                    help="how to place the truncation radii.  'cfoot' targets "
                         'the footprint/FWHM ratio, which is what the matching '
                         "uses; 'frac' targets the recovered fraction; 'both' "
                         'takes the union (default).')
    ap.add_argument('--cfoot-lo', type=float, default=1.2,
                    help='lowest footprint/FWHM ratio to tabulate')
    ap.add_argument('--cfoot-hi', type=float, default=3.0,
                    help='highest footprint/FWHM ratio to tabulate.  The seven '
                         'getsf catalogues sit at 1.8-2.0 (16-84%%: 1.5-2.5), '
                         'so 1.2-3.0 brackets the observed range with margin '
                         'for interpolation at the edges.')
    ap.add_argument('--cfoot-min', type=float, default=None,
                    help='discard rows below this footprint/FWHM ratio; '
                         'ratios much below ~1.5 correspond to footprints that '
                         'getsf does not produce')
    ap.add_argument('--cfoot-max', type=float, default=None,
                    help='discard rows above this footprint/FWHM ratio')
    ap.add_argument('--fwhm-floor', type=float, default=1.0,
                    help='skip radii whose truncated FWHM falls below this '
                         'multiple of the beam.  Real extractions do not '
                         'return sizes below ~1.5 beams, because background '
                         'fluctuations broaden the measured profile, so rows '
                         'below the floor describe truncations that cannot '
                         'occur in practice.  1.5 imposes the observed limit; '
                         '1.0 keeps everything physically resolvable.')
    ap.add_argument('--col-i', type=int, default=1)
    ap.add_argument('--col-j', type=int, default=2)
    ap.add_argument('--col-k', type=int, default=3)
    ap.add_argument('--col-sd', type=int, default=4)
    args = ap.parse_args()

    bands = [b.strip() for b in args.bands.split(',')]
    for b in bands:
        if b not in BANDS:
            raise SystemExit('unknown band %r; known: %s' % (b, sorted(BANDS)))

    rows = [l for l in open(args.grid)
            if l.strip() and not l.lstrip().startswith('#')]

    out = open(args.out, 'w')
    out.write('# recovery quantities vs truncation radius\n')
    out.write('#   generated by tabulate_recovery.py\n')
    out.write('#   frac is normalised to 1 at r = R_conv (full footprint)\n')
    out.write('#   peak_over_rim is the model counterpart of an observed peak/noise\n')
    out.write('#   reliable = 1 where frac >= %.2f\n' % args.frac_min)
    out.write('#   absorp = 1 where the central pixel is negative (absorption\n')
    out.write('#            against a warmer background; expected at 70/100 um for\n')
    out.write('#            cold cores, and a flag for the depression models).\n')
    out.write('#            Such rows carry no frac: the profile is not a\n')
    out.write('#            declining positive source and truncation is undefined.\n')
    out.write('# %4s %3s %3s %5s %11s %9s %9s %11s %11s %11s %9s %9s %9s %11s %4s %6s\n'
              % ('i', 'j', 'k', 'band', 'SD_emb', 'r_as', 'r_Rconv', 'peak',
                 'I_rim', 'peak_rim', 'frac', 'fwhm_rec', 'cfoot', 'cpeakmean',
                 'rel', 'absorp'))

    nmiss = nnode = nmaskmiss = nabs = nfloor = ncfoot = 0
    for line in rows:
        f = line.split()
        i, j, k = (int(float(f[args.col_i])), int(float(f[args.col_j])),
                   int(float(f[args.col_k])))
        sd = float(f[args.col_sd])
        nnode += 1
        for b in bands:
            beam, token = BANDS[b]
            path = stamp_path(args.stamp_root, i, j, k, token, beam,
                              args.stamp_pattern)
            if not os.path.exists(path):
                nmiss += 1
                continue
            with fits.open(path, memmap=True) as hd:
                stamp = np.nan_to_num(np.squeeze(hd[0].data).astype(float))
            rr, r_centers, prof = radial_profile(stamp)
            R_conv = np.nan
            if args.mask_pattern:
                mpath = mask_path(args.stamp_root, i, j, k, beam,
                                  args.mask_pattern)
                if os.path.exists(mpath):
                    R_conv = extent_from_mask(mpath)
                else:
                    nmaskmiss += 1
            if not np.isfinite(R_conv):
                R_conv = full_extent(r_centers, prof)
            if not np.isfinite(R_conv) or R_conv <= 0:
                continue
            peak_val = prof[0]
            if not np.isfinite(peak_val) or peak_val <= 0:
                # absorption (or an empty stamp): record the fact, skip frac
                nabs += 1
                out.write('%6d %3d %3d %5s %11.4e %9.2f %9.4f %11.4e %11.4e '
                          '%11.4e %9s %9s %9s %11s %4d %6d\n'
                          % (i, j, k, b, sd, R_conv, 1.0, peak_val, np.nan,
                             np.nan, 'nan', 'nan', 'nan', 'nan', 0, 1))
                continue
            # total flux with the full footprint -> frac(R_conv) = 1
            i_rim_full = prof[int(np.clip(np.searchsorted(r_centers, R_conv),
                                          0, len(prof) - 1))]
            total = np.clip(stamp[rr <= R_conv] - i_rim_full, 0.0, None).sum()
            if total <= 0:
                continue
            radii = np.array([])
            if args.sample in ('frac', 'both'):
                radii = np.concatenate([radii, radii_for_fracs(
                    rr, stamp, R_conv,
                    np.geomspace(args.frac_lo, 1.0, args.nr), beam)])
            if args.sample in ('cfoot', 'both'):
                radii = np.concatenate([radii, radii_for_cfoot(
                    r_centers, prof, R_conv,
                    np.linspace(args.cfoot_lo, args.cfoot_hi, args.nr), beam)])
            if radii.size == 0:
                continue
            radii = np.unique(np.round(radii, 3))
            radii = np.clip(radii, max(args.r_min_frac * R_conv, beam / 4.0),
                            R_conv)
            for r in radii:
                (i_rim, por, frac, fwhm_rec, cff, cpm,
                 cslp) = quantities_at(rr, r_centers, prof, stamp, r, total)
                if (np.isfinite(fwhm_rec)
                        and fwhm_rec < args.fwhm_floor * beam):
                    nfloor += 1
                    continue
                if args.cfoot_min is not None and np.isfinite(cff) \
                        and cff < args.cfoot_min:
                    ncfoot += 1
                    continue
                if args.cfoot_max is not None and np.isfinite(cff) \
                        and cff > args.cfoot_max:
                    ncfoot += 1
                    continue
                rel = 1 if (np.isfinite(frac) and frac >= args.frac_min) else 0
                out.write('%6d %3d %3d %5s %11.4e %9.2f %9.4f %11.4e %11.4e '
                          '%11.4e %9.4f %9.2f %9.4f %11.4f %4d %6d\n'
                          % (i, j, k, b, sd, r, r / R_conv, peak_val, i_rim,
                             por, frac, fwhm_rec, cff, cpm, rel, 0))
    out.close()
    print('wrote %s  (%d nodes x %d bands x %d radii)'
          % (args.out, nnode, len(bands), args.nr))
    if nmiss:
        print('WARNING: %d missing stamps -- check --stamp-pattern' % nmiss)
    if ncfoot:
        print('NOTE: %d rows outside the requested footprint/FWHM window'
              % ncfoot)
    if nfloor:
        print('NOTE: %d rows skipped with truncated FWHM below %.2f beams'
              % (nfloor, args.fwhm_floor))
    if nabs:
        print('NOTE: %d (node,band) cases in absorption -- flagged, no frac'
              % nabs)
    if nmaskmiss:
        print('NOTE: %d masks not found -- R_conv fell back to the profile '
              'threshold for those; check --mask-pattern' % nmaskmiss)


if __name__ == '__main__':
    main()
