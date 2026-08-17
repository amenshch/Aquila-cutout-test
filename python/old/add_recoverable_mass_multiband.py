#!/usr/bin/env python3
"""
add_recoverable_mass_multiband.py -- per-band recoverable-flux fractions for the
BE model grid  (stage B of the multi-band flux-recovery scheme).

Background
----------
The original add_recoverable_mass.py estimates ONE recovered fraction, from the
surface-density stamp, using an empirical getsf detection floor

    N_src_min(N_back) = 1.971e21 * (N_back/1e22)**1.342     [cm^-2]

and multiplies M_SED by it.  That entangles two physically distinct effects:
flux lost outside the detectable footprint, and the temperature bias of the SED
fit.  It also assumes the loss is grey, which it is not: each band has its own
noise, its own beam and therefore its own footprint, so each band loses a
DIFFERENT fraction.  The measured SED is distorted in colour, not merely in
normalisation, and the fitted temperature is biased as a result.

This script computes frac_b per band, so that the correction can be applied to
the band fluxes BEFORE the SED fit:

    F_b(corrected) = F_b(measured) / frac_b   ->   re-fit SED   ->  T, M_SED
    M_SED -> M_BE   via the idealised (temperature-bias) correction only.

The residual model-dependent step is then the idealised correction, which is
~5x more accurate in leave-one-out than the entangled recoverable path.

The models are noiseless
------------------------
No noise is added to the models.  The noise information enters solely through
the detection floor, which is measured from real extractions -- exactly as in
the original single-band version.  Two ways to obtain the per-band floor are
supported:

  --floor-table FILE   Per-band empirical floors, as written by
                       calibrate_floors.py from a real getsf catalogue:
                           band  A  p  I_ref     floor = A*(I_bg/I_ref)**p
                       Recommended.  Scatter about the fit is 1.03-1.3x, and
                       calibrating per field removes any dependence on a single
                       reference cloud.

  --floor-physical     Derive every band from the SINGLE surface-density floor
                       by converting the column-density fluctuation through the
                       dust emission law at the background temperature,
                           dI_b = dN * mu * m_H * kappa_b * B_b(T_bg),
                       then correcting for beam smoothing of the cirrus,
                           floor_b *= (theta_b/theta_SD)**BEAM_INDEX.
                       Tested against the independent per-band calibration on
                       Aquila: 0.97x at 250 um, 1.18x at 350, 1.77x at 500
                       before the beam term; BEAM_INDEX = -0.87 absorbs the
                       trend.  Useful as a cross-check and as a fallback for
                       bands whose empirical envelope fit is poor.

Background level per band
-------------------------
Needed in BOTH modes, because the stamps are background-subtracted and only
SD_emb is stored.  Obtained from the same dust law,
    I_bg,b = SD_emb * mu * m_H * kappa_b * B_b(T_bg),
with T_bg taken per node from the catalogue's Td_emb column when available.
On Aquila a single T_bg reproduces the observed background in the three SPIRE
bands to within 0.05 K, so this conversion is well founded.

Usage
-----
    python add_recoverable_mass_multiband.py \
        --grid  bes_model_grid_final2_catalog \
        --stamp-root /path/to/260616_RT_BES_radmc3d \
        --out   bes_model_grid_final2_catalog_mb \
        --floor-table floors_aquila.txt \
        --col-msed3bs 20 --col-msed3bsl 23 --col-tdemb 12

Column arguments are 0-BASED indices into line.split().

Requires numpy and astropy.
"""

import argparse
import os

import numpy as np
from astropy.io import fits

PIX_ARCSEC = 3.0          # stamps are resampled to a common 3" grid
ETA = 2.0                 # getsf footprint target, eta * half-maximum radius

# physical constants (cgs)
H_PL, K_B, C_L = 6.62607015e-27, 1.380649e-16, 2.99792458e10
M_H, MU = 1.6737e-24, 2.8            # mean molecular weight per H2

# Ossenkopf & Henning beta=2 power law, per gram of GAS
KAPPA_100 = 0.83364                  # cm^2 g^-1 at 100 um
BEAM_INDEX = -0.87                   # cirrus beam-smoothing exponent (see above)

# Bands to process.  label, wavelength [um], beam FWHM ["], stamp token.
# The surfdens map is listed first and is the reference for the beam term.
BANDS = [
    ('SD',  None, 13.5, 'surfdens'),
    ('070',   70,  8.4, '070'),
    ('100',  100, 10.0, '100'),
    ('160',  160, 13.5, '160'),
    ('250',  250, 18.2, '250'),
    ('350',  350, 24.9, '350'),
    ('500',  500, 36.3, '500'),
]

# surface-density detection floor (the original calibration)
FLOOR_NORM, FLOOR_REF, FLOOR_INDEX = 1.971e21, 1.0e22, 1.342


def sd_floor(n_back):
    """Empirical getsf detection floor on the surface-density map [cm^-2]."""
    return FLOOR_NORM * (n_back / FLOOR_REF) ** FLOOR_INDEX


def kappa(lam_um):
    """Dust opacity per gram of gas [cm^2/g]."""
    return KAPPA_100 * (100.0 / lam_um) ** 2


def planck(lam_um, T):
    """B_nu [cgs: erg s^-1 cm^-2 sr^-1 Hz^-1]."""
    nu = C_L / (lam_um * 1e-4)
    return 2 * H_PL * nu ** 3 / C_L ** 2 / np.expm1(H_PL * nu / (K_B * T))


def column_to_intensity(n_h2, lam_um, T):
    """Convert an H2 column [cm^-2] to band intensity [MJy/sr] at temperature T."""
    sigma = n_h2 * MU * M_H
    return sigma * kappa(lam_um) * planck(lam_um, T) / 1e-17


def read_floor_table(path):
    """band -> (A, p, I_ref) from calibrate_floors.py output."""
    out = {}
    for line in open(path):
        if line.lstrip().startswith('#') or not line.strip():
            continue
        f = line.split()
        out[f[0]] = (float(f[1]), float(f[2]), float(f[3]))
    return out


def band_floor(label, lam_um, beam, sd_emb, t_bg, floor_table, band_key):
    """Return (floor, I_bg) in the band's native units.

    Surface-density bands stay in cm^-2; flux bands are in MJy/sr.
    """
    if label == 'SD':
        i_bg = sd_emb
    else:
        i_bg = column_to_intensity(sd_emb, lam_um, t_bg)

    if floor_table is not None and band_key in floor_table:
        A, p, i_ref = floor_table[band_key]
        return A * (i_bg / i_ref) ** p, i_bg

    # physical fallback: convert the surfdens floor through the dust law,
    # then correct for beam smoothing of the cirrus fluctuations
    dn = sd_floor(sd_emb)
    if label == 'SD':
        return dn, i_bg
    f = column_to_intensity(dn, lam_um, t_bg)
    return f * (beam / BANDS[0][2]) ** BEAM_INDEX, i_bg


def recoverable_fraction(stamp, floor, pix_arcsec=PIX_ARCSEC, eta=ETA):
    """Recovered flux fraction and shape diagnostics for one band.

    `stamp` is the background-subtracted model image in that band; `floor` is
    the detection threshold in the same units.  Identical in structure to the
    original single-band routine -- the only change is that the floor is passed
    in rather than derived internally from the surface density, so that each
    band can carry its own.

    Returns (frac, r_foot, pedestal, fwhm_rec, conc_footfwhm, conc_peakmean,
             conc_slope, peak_over_floor).
    """
    nan8 = (0.0, 0.0, 0.0, np.nan, np.nan, np.nan, np.nan, np.nan)
    ny, nx = stamp.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(yy - cy, xx - cx) * pix_arcsec

    rmax = int(np.floor(rr.max() / pix_arcsec))
    r_edges = (np.arange(rmax + 2) - 0.5) * pix_arcsec
    prof = np.zeros(rmax + 1)
    for ir in range(rmax + 1):
        m = (rr >= r_edges[ir]) & (rr < r_edges[ir + 1])
        prof[ir] = stamp[m].mean() if m.any() else 0.0
    r_centers = np.arange(rmax + 1) * pix_arcsec

    peak = prof[0]
    if peak <= 0:
        return nan8

    half = peak / 2.0
    below_half = np.where(prof < half)[0]
    h_n = r_centers[below_half[0]] if below_half.size else r_centers[-1]
    r_eta = eta * h_n

    below_thr = np.where(prof < floor)[0]
    r_noise = r_centers[below_thr[0]] if below_thr.size else r_centers[-1]
    r_foot = min(r_eta, r_noise)

    irim = min(int(np.searchsorted(r_centers, r_foot)), len(prof) - 1)
    pedestal = prof[irim]

    total = stamp[stamp > 0].sum()
    if total <= 0:
        return (0.0, r_foot, pedestal, np.nan, np.nan, np.nan, np.nan,
                peak / floor if floor > 0 else np.nan)

    foot = rr <= r_foot
    resid = np.clip(stamp[foot] - pedestal, 0.0, None)
    recovered = resid.sum()

    prof_sub = np.clip(prof - pedestal, 0.0, None)
    prof_sub[r_centers > r_foot] = 0.0
    pk = prof_sub[0]
    if pk > 0 and r_foot > 0:
        bh = np.where(prof_sub < pk / 2.0)[0]
        r_half = r_centers[bh[0]] if bh.size else r_foot
        fwhm_rec = 2.0 * r_half
        conc_footfwhm = (2.0 * r_foot) / fwhm_rec if fwhm_rec > 0 else np.nan
        mean_in = resid[resid > 0].mean() if np.any(resid > 0) else np.nan
        conc_peakmean = pk / mean_in if (mean_in and mean_in > 0) else np.nan
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

    return (recovered / total, r_foot, pedestal, fwhm_rec, conc_footfwhm,
            conc_peakmean, conc_slope, peak / floor if floor > 0 else np.nan)


def stamp_path(root, i, j, k, token, beam, pattern):
    res = ('r%.1f' % beam).replace('.', 'p')
    name = pattern.format(band=token, res=res)
    return os.path.join(root, 'cSD_%02d' % i, 'M_%02d' % j, '%02d' % k, name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid', required=True)
    ap.add_argument('--stamp-root', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--floor-table', default=None,
                    help='per-band floors from calibrate_floors.py (recommended)')
    ap.add_argument('--floor-physical', action='store_true',
                    help='derive all band floors from the surfdens floor + dust law')
    ap.add_argument('--stamp-pattern', default='nc.{band}.bs.{res}x0.rs3p0as.fits')
    ap.add_argument('--t-bg-default', type=float, default=15.4,
                    help='background dust temperature if Td_emb is unavailable')
    ap.add_argument('--col-i', type=int, default=1)
    ap.add_argument('--col-j', type=int, default=2)
    ap.add_argument('--col-k', type=int, default=3)
    ap.add_argument('--col-sd', type=int, default=4)
    ap.add_argument('--col-tdemb', type=int, default=None,
                    help='0-based column of the embedding dust temperature')
    ap.add_argument('--col-msed3bs', type=int, required=True)
    ap.add_argument('--col-msed3bsl', type=int, required=True)
    ap.add_argument('--bands', default='SD,160,250,350,500',
                    help='comma-separated band labels to process')
    args = ap.parse_args()

    if args.floor_table is None and not args.floor_physical:
        raise SystemExit('give either --floor-table or --floor-physical')
    table = read_floor_table(args.floor_table) if args.floor_table else None
    wanted = [b.strip() for b in args.bands.split(',')]
    bands = [b for b in BANDS if b[0] in wanted]
    # map our labels onto the two-digit band keys of the floor table, in order
    keys = {lab: '%02d' % (n + 1) for n, (lab, _, _, _) in enumerate(bands)}

    header, data = [], []
    for line in open(args.grid):
        (header if line.lstrip().startswith('#') else data).append(line.rstrip('\n'))
    data = [d for d in data if d.strip()]

    cols = []
    for lab, _, _, _ in bands:
        cols += ['frac_%s' % lab, 'peakfloor_%s' % lab, 'rfoot_%s' % lab,
                 'fwhmrec_%s' % lab, 'concfoot_%s' % lab, 'concpm_%s' % lab,
                 'concslope_%s' % lab]

    out = list(header)
    out.append('#')
    out.append('#   Per-band recoverable fractions (stage B).  Floor source: %s'
               % ('empirical table %s' % args.floor_table if table else
                  'physical conversion of the surfdens floor'))
    for c in cols:
        out.append('#   %-16s ..... per-band recovery diagnostic' % c)

    nmiss = 0
    for line in data:
        f = line.split()
        i, j, k = int(float(f[args.col_i])), int(float(f[args.col_j])), int(float(f[args.col_k]))
        sd = float(f[args.col_sd])
        t_bg = (float(f[args.col_tdemb]) if args.col_tdemb is not None
                else args.t_bg_default)
        if not np.isfinite(t_bg) or t_bg <= 0:
            t_bg = args.t_bg_default
        vals = []
        for lab, lam, beam, token in bands:
            path = stamp_path(args.stamp_root, i, j, k, token, beam, args.stamp_pattern)
            floor, _ = band_floor(lab, lam, beam, sd, t_bg, table, keys[lab])
            if not os.path.exists(path):
                nmiss += 1
                vals += ['%12.5e' % np.nan] * 7
                continue
            with fits.open(path, memmap=True) as hd:
                stamp = np.nan_to_num(np.squeeze(hd[0].data).astype(float))
            (frac, r_foot, _ped, fwhm_rec, cff, cpm, cslp,
             pf) = recoverable_fraction(stamp, floor)
            vals += ['%12.5e' % v for v in
                     (frac, pf, r_foot, fwhm_rec, cff, cpm, cslp)]
        out.append(line + ' ' + ' '.join(vals))

    with open(args.out, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('wrote %s  (%d models, %d bands)' % (args.out, len(data), len(bands)))
    if nmiss:
        print('WARNING: %d missing stamps -> nan columns; check --stamp-pattern' % nmiss)


if __name__ == '__main__':
    main()
