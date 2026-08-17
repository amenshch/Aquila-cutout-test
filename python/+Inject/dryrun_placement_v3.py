#!/usr/bin/env python
"""dryrun_placement_v3.py -- exercise the model selection and the placement of
run_inject_v3.py without injecting anything into the images.

The purpose is to predict, before any getsf run, the joint distribution of
background column density and contrast that the injected cores would have, and
to compare it with the distribution of the real sources.  Nothing is written
into the maps and no model stamp is needed, so this runs anywhere the
surface-density map, its coverage mask, the footprint map, the model grid
catalogue and the target table are available.

For each placed core the following are recorded:

    tag           grid node label
    x_pix, y_pix  position, one-based
    SD_emb        embedding column density of the node               (cm^-2)
    local_Sigma   column density of the map averaged over a disk of
                  radius equal to the observed size of the model     (cm^-2)
    contrast      1 + ICSDbs / local_Sigma, the contrast the core
                  would show at the position it was given            (-)
    M_BE          true mass of the node                             (M_sun)
    M_SED         reported four-band SED mass of the node           (M_sun)
    FWHM_SD       observed major axis at half maximum of the node   (arcsec)
    R_BE_as       truncation radius of the node                     (arcsec)

Usage:
    python dryrun_placement_v3.py [config_module] [output_file]
"""
import os
import sys
import math

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import inject                                     # noqa: E402


def main(cfg_name='inject_config_dryrun', out_path='dryrun_placement_v3.txt'):
    sys.argv = [sys.argv[0], cfg_name]
    import run_inject_v3 as R                      # noqa: E402

    rng = np.random.default_rng(R.RANDOM_SEED)

    # ---------------- 1. selection ----------------------------------
    models = R._models_from_population(
        R.CMF_N_CORES, R.CMF_SLOPE,
        R.CMF_MASS_RANGE[0], R.CMF_MASS_RANGE[1],
        rng,
        R._load_joint_target(R.JOINT_TARGET_FILE),
        fwhm_min=R.FWHM_MIN_AS, fwhm_max=R.FWHM_MAX_AS)
    print('\nselection produced %d distinct grid nodes for %d cores'
          % (len(models), sum(m['n_cores'] for m in models)))

    # ---------------- 2. field --------------------------------------
    sigma = fits.getdata(R.SIGMA_PATH).astype(float)
    omask = fits.getdata(R.SIGMA_OMASK_PATH) > 0
    coverage = omask & np.isfinite(sigma) & (sigma > 0)
    field = dict(sigma=sigma, coverage=coverage, shape=sigma.shape)
    foots = fits.getdata(R.FOOTS_PATH).astype(float)
    forbidden = (foots > 0) if R.AVOID_REAL_FOOTPRINTS \
        else np.zeros(foots.shape, bool)
    print('field %s, coverage %d pixels, footprints cover %.1f%%'
          % (sigma.shape, coverage.sum(), 100.0 * np.mean(foots > 0)))

    # node properties, looked up once
    g = R._grid_arrays()
    idx = {t: i for i, t in enumerate(g['tag'])}

    # ---------------- 3. placement ----------------------------------
    md_list = []
    for m in models:
        i = idx[m['tag']]
        md_list.append(dict(tag=m['tag'], n_cores=m['n_cores'],
                            sd_emb=float(g['SD_emb'][i]),
                            r_be_arcsec=float(g['R_BE_as'][i]),
                            fwhm_arcsec=float(g['FWHM_as'][i]),
                            icsd=float(g['ICSDbs'][i]),
                            m_be=float(g['M_BE'][i]),
                            m_sed=float(g['M_SED'][i]),
                            targets=list(m['targets'])))
    # same ordering as run_inject_v3.main(): by exclusion radius, descending
    md_list.sort(key=lambda d: (-R.model_radius_as(d['fwhm_arcsec'],
                                                   d['r_be_arcsec'],
                                                   d['tag']),
                                -d['sd_emb']))

    # the stamp half-width is not known without the stamps; use the radius the
    # model occupies on the sky, so that the edge margin is realistic
    placed = []
    used = []
    n_want_total = sum(d['n_cores'] for d in md_list)
    for md in md_list:
        half = int(math.ceil(max(md['fwhm_arcsec'], 36.3) / inject.PIX_ARCSEC))
        valid = inject.valid_placement_mask(
            field, sigma_cloud=md['sd_emb'],
            R_BE_arcsec=md['r_be_arcsec'],
            sqrt2_window=True, forbidden=forbidden,
            stamp_halfwidth=half,
            window=R.SIGMA_WINDOW,
            forbid_margin_as=(R.real_source_margin_as(md['fwhm_arcsec'],
                                                      md['r_be_arcsec'],
                                                      md['tag'])
                              if R.AVOID_REAL_FOOTPRINTS else None),
            local_average_radius_as=md['fwhm_arcsec'])
        ys, xs = np.where(valid)
        if len(ys) == 0:
            print('  %-12s no valid position' % md['tag'])
            continue
        lsig_map = inject.local_sigma(field['sigma'],
                                      md['fwhm_arcsec'] / inject.PIX_ARCSEC)
        lsig = lsig_map[ys, xs]
        jitter = 1e-6 * rng.random(len(lsig))
        taken = np.zeros(len(ys), bool)
        rad = R.model_radius_as(md['fwhm_arcsec'], md['r_be_arcsec'], md['tag'])

        found = 0
        for (tgt_sig, tgt_con, tgt_siz) in md['targets'][:md['n_cores']]:
            key = np.abs(np.log10(np.maximum(lsig, 1e-30) / tgt_sig)) + jitter
            key = np.where(taken, np.inf, key)
            order = np.argsort(key)
            for k in order:
                if not np.isfinite(key[k]):
                    break
                y, x = int(ys[k]), int(xs[k])
                if used:
                    uc = np.asarray(used, float)
                    smin = np.array([R.min_separation_as(rad, r)
                                     for r in uc[:, 2]]) / inject.PIX_ARCSEC + 1.0
                    d2 = (y - uc[:, 0]) ** 2 + (x - uc[:, 1]) ** 2
                    apply = np.array([R.separation_applies(md['fwhm_arcsec'], f)
                                      for f in uc[:, 3]])
                    ok = bool(np.all(d2[apply] >= (smin * smin)[apply]))
                else:
                    ok = True
                taken[k] = True
                if not ok:
                    continue
                used.append((y, x, rad, md['fwhm_arcsec']))
                ls = float(lsig_map[y, x])
                placed.append(dict(tag=md['tag'], x=x + 1, y=y + 1,
                                   sd_emb=md['sd_emb'], local_sigma=ls,
                                   contrast=1.0 + md['icsd'] / ls,
                                   m_be=md['m_be'], m_sed=md['m_sed'],
                                   fwhm=md['fwhm_arcsec'], r_be=md['r_be_arcsec'],
                                   tgt_sigma=tgt_sig, tgt_contrast=tgt_con,
                                   tgt_size=tgt_siz))
                found += 1
                if (R.MAX_FOOTPRINT_OVERLAP is not None
                        or R.OVERLAP_FREE_FWHM_AS is None
                        or md['fwhm_arcsec'] <= R.OVERLAP_FREE_FWHM_AS):
                    r = int(round(0.5 * R.min_separation_as(rad, rad)
                                  / inject.PIX_ARCSEC))
                    if r > 0:
                        y0, y1 = max(0, y - r), min(field['shape'][0], y + r + 1)
                        x0, x1 = max(0, x - r), min(field['shape'][1], x + r + 1)
                        yy, xx = np.ogrid[y0 - y:y1 - y, x0 - x:x1 - x]
                        forbidden[y0:y1, x0:x1] |= (xx ** 2 + yy ** 2 < r ** 2)
                break

    print('\nplaced %d of %d requested cores' % (len(placed), n_want_total))

    with open(out_path, 'w') as f:
        f.write('# Predicted properties of the injected cores, from the '
                'selection and placement of run_inject_v3.py.\n')
        f.write('# No image was modified.  Seed %d.\n' % R.RANDOM_SEED)
        f.write('#\n# Columns\n')
        f.write('#   1 model        grid node label\n')
        f.write('#   2 x_pix        column of the center, one-based\n')
        f.write('#   3 y_pix        row of the center, one-based\n')
        f.write('#   4 SD_emb       embedding column of the node        (cm^-2)\n')
        f.write('#   5 local_Sigma  map column averaged over the        (cm^-2)\n')
        f.write('#                  observed size of the model\n')
        f.write('#   6 contrast     1 + ICSDbs / local_Sigma            (-)\n')
        f.write('#   7 M_BE         true mass of the node               (M_sun)\n')
        f.write('#   8 M_SED        reported four-band SED mass         (M_sun)\n')
        f.write('#   9 FWHM_SD      observed size of the node           (arcsec)\n')
        f.write('#  10 R_BE_as      truncation radius of the node       (arcsec)\n')
        f.write('#  11 tgt_Sigma    background column of the real       (cm^-2)\n')
        f.write('#                  source this core imitates\n')
        f.write('#  12 tgt_contrast contrast of that real source        (-)\n')
        f.write('#  13 tgt_FWHM     observed size of that real source   (arcsec)\n')
        f.write('# %-12s %6s %6s %12s %12s %10s %10s %10s %8s %9s %12s %12s %9s\n'
                % ('model', 'x_pix', 'y_pix', 'SD_emb', 'local_Sigma',
                   'contrast', 'M_BE', 'M_SED', 'FWHM_SD', 'R_BE_as',
                   'tgt_Sigma', 'tgt_contrast', 'tgt_FWHM'))
        for p in placed:
            f.write('  %-12s %6d %6d %12.4e %12.4e %10.4f %10.5g %10.5g '
                    '%8.2f %9.2f %12.4e %12.4f %9.2f\n'
                    % (p['tag'], p['x'], p['y'], p['sd_emb'], p['local_sigma'],
                       p['contrast'], p['m_be'], p['m_sed'], p['fwhm'],
                       p['r_be'], p['tgt_sigma'], p['tgt_contrast'],
                       p['tgt_size']))
    print('wrote %s' % out_path)

    # ---------------- 4. drawn versus realised ----------------------
    if placed:
        ts = np.array([p['tgt_sigma'] for p in placed])
        rs = np.array([p['local_sigma'] for p in placed])
        tc = np.array([p['tgt_contrast'] for p in placed])
        rc = np.array([p['contrast'] for p in placed])
        print('\nDrawn versus realised, over the %d placed cores:' % len(placed))
        print('  background column, median realised / drawn = %.3f, '
              'median |log10 ratio| = %.3f dex'
              % (np.median(rs / ts), np.median(np.abs(np.log10(rs / ts)))))
        print('  contrast,          median realised / drawn = %.3f, '
              'median |log10 ratio| = %.3f dex'
              % (np.median(rc / tc), np.median(np.abs(np.log10(rc / tc)))))
        from scipy.stats import spearmanr as _sp
        print('  Spearman rank correlation of contrast on background column:')
        print('    among the real sources drawn      %+.3f'
              % _sp(np.log10(ts), np.log10(tc)).statistic)
        print('    among the injected cores          %+.3f'
              % _sp(np.log10(rs), np.log10(rc)).statistic)
    return out_path


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0] if len(a) > 0 else 'inject_config_dryrun',
         a[1] if len(a) > 1 else 'dryrun_placement_v3.txt')
