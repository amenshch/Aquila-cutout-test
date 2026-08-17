#!/usr/bin/env python
"""smoke_test_v3.py -- run run_inject_v3.main() end to end on synthetic inputs.

Purpose.  dryrun_placement_v3.py exercises the model selection and the
placement, but it builds the field itself and never calls main(), so a name
that main() needs and does not have goes undetected.  This script runs main()
in full, on a small synthetic field with synthetic stamps, and therefore
catches that class of error.  It verifies that the script runs, not that the
numbers are right; the numbers are the business of dryrun_placement_v3.py and
check_injection_match_v3.py.

Everything is written under a temporary directory and removed afterwards, so
nothing in the working tree is touched.

Usage:
    python smoke_test_v3.py [--keep]
"""
import os
import shutil
import sys
import tempfile

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
BANDS = ['000', '070', '160', '250', '350', '500']
NPIX = 300
PIX_ARCSEC = 3.0


def _wcs_header():
    h = fits.Header()
    h['CTYPE1'] = 'RA---TAN'; h['CTYPE2'] = 'DEC--TAN'
    h['CRVAL1'] = 277.5;      h['CRVAL2'] = -2.0
    h['CRPIX1'] = NPIX / 2.0; h['CRPIX2'] = NPIX / 2.0
    h['CDELT1'] = -PIX_ARCSEC / 3600.0
    h['CDELT2'] = PIX_ARCSEC / 3600.0
    return h


def build_field(d, real_sigma_path, real_omask_path):
    """A synthetic field whose column density spans the grid's rungs."""
    sig = fits.getdata(real_sigma_path).astype(float)
    om = fits.getdata(real_omask_path) > 0
    # take a corner of the real map so the column distribution is realistic
    sub = sig[:NPIX, :NPIX].copy()
    sub[~np.isfinite(sub)] = 3.0e21
    sub[sub <= 0] = 3.0e21
    hdr = _wcs_header()
    paths = {}
    for b in BANDS:
        arr = sub if b == '000' else sub / 1.0e21   # arbitrary intensity units
        p = os.path.join(d, 'field.%s.fits' % b)
        fits.writeto(p, arr.astype(np.float32), hdr, overwrite=True)
        po = os.path.join(d, 'field.%s.omask.fits' % b)
        fits.writeto(po, np.ones_like(arr, dtype=np.float32), hdr,
                     overwrite=True)
        paths[b] = (p, po)
    pf = os.path.join(d, 'field.foots.fits')
    foot = np.zeros((NPIX, NPIX), np.float32)
    foot[40:60, 40:60] = 1.0
    fits.writeto(pf, foot, hdr, overwrite=True)
    return paths, pf


def build_stamps(paths_of_tag, fwhm_of):
    """One Gaussian stamp per band per model, at the exact paths the script
    will look for.  The layout is produced by run_inject_v3._stamp_files_from_tag,
    so it is taken from there rather than reconstructed here."""
    n = 73
    c = n // 2
    yy, xx = np.mgrid[0:n, 0:n]
    for tag, sf in paths_of_tag.items():
        s = max(fwhm_of[tag], 15.0) / PIX_ARCSEC / 2.3548
        g = np.exp(-((xx - c) ** 2 + (yy - c) ** 2) / (2 * s * s))
        for b, path in sf.items():
            if path is None:
                continue
            os.makedirs(os.path.dirname(path), exist_ok=True)
            arr = g * (5.0e21 if b == '000' else 10.0)
            fits.writeto(path, arr.astype(np.float32), overwrite=True)


def write_config(d, paths, foots, stamp_root, out_dir):
    txt = ['SET_TAG = "smoke"',
           'CMF_SLOPE = 2.0',
           'CMF_N_CORES = 60',
           'CMF_MASS_RANGE = (0.05, 6.0)',
           'JOINT_TARGET_FILE = %r' % os.path.join(HERE,
                                                   'joint_target_Aquila_v3.txt'),
           'W_MASS = 0.1', 'W_CONTRAST = 1.0', 'W_SIZE = 1.0',
           'MODEL_PICK_TOLERANCE = 0.005',
           'FWHM_MIN_AS = 15.0', 'FWHM_MAX_AS = 80.0',
           'OVERLAP_FREE_FWHM_AS = 50.0',
           'SIGMA_WINDOW = "symmetric"', 'SEP_FACTOR = 1.3',
           'AVOID_REAL_FOOTPRINTS = True', 'FOOTPRINT_MARGIN_AS = 18.2',
           'FLATTEN_BACKGROUND = False',
           'RANDOM_SEED = 7', 'N_CORES = 5',
           'CONTROL_ISOLATED = %s' % ('--isolated' in sys.argv),
           'OUT_DIR = %r' % out_dir,
           'IMG_PATHS = {%s}' % ', '.join(
               '%r: %r' % (b, paths[b][0]) for b in BANDS if b != '000'),
           'OMASK_PATHS = {%s}' % ', '.join(
               '%r: %r' % (b, paths[b][1]) for b in BANDS if b != '000'),
           'SIGMA_PATH = %r' % paths['000'][0],
           'SIGMA_OMASK_PATH = %r' % paths['000'][1],
           'HIRES_SURFDENS_PATH = %r' % paths['000'][0],
           'FOOTS_PATH = %r' % foots,
           'STAMP_ROOT = %r' % stamp_root,
           'BES_CATALOG = %r' % os.path.join(ROOT, 'cats',
                                             'bes_model_grid_final2_catalog'),
           ]
    p = os.path.join(d, 'smoke_config.py')
    open(p, 'w').write('\n'.join(txt) + '\n')
    return p


def main():
    keep = '--keep' in sys.argv
    _cwd0 = os.getcwd()
    d = tempfile.mkdtemp(prefix='smoke_v3_')
    try:
        out_dir = os.path.join(d, 'out')
        os.makedirs(out_dir)
        paths, foots = build_field(
            d,
            os.path.join(ROOT, 'fits', 'hi.surface.density.r13p5.fits'),
            os.path.join(ROOT, 'fits', 'hi.surface.density.r13p5.omask.fits'))
        stamp_root = os.path.join(d, 'stamps')
        os.makedirs(stamp_root, exist_ok=True)

        cfg = write_config(d, paths, foots, stamp_root, out_dir)
        # run_inject_v3.py resolves the configuration relative to the current
        # working directory, so move there before importing it
        _cwd0 = os.getcwd()
        os.chdir(d)
        sys.path.insert(0, d)
        sys.argv = [sys.argv[0], 'smoke_config']
        sys.path.insert(0, HERE)
        import run_inject_v3 as R

        # the sampler chooses the tags, so the stamps must be built after it
        models = R._models_from_population(
            R.CMF_N_CORES, R.CMF_SLOPE, R.CMF_MASS_RANGE[0],
            R.CMF_MASS_RANGE[1], np.random.default_rng(R.RANDOM_SEED),
            R._load_joint_target(R.JOINT_TARGET_FILE),
            fwhm_min=R.FWHM_MIN_AS, fwhm_max=R.FWHM_MAX_AS)
        build_stamps(
            {m['tag']: R._stamp_files_from_tag(m['tag'], stamp_root)
             for m in models},
            {m['tag']: m['fwhm_as'] for m in models})

        print('\n' + '=' * 70)
        print('running run_inject_v3.main()')
        print('=' * 70)
        R.main()

        produced = sorted(os.listdir(out_dir))
        print('\n' + '=' * 70)
        print('SMOKE TEST PASSED: main() completed, %d files written'
              % len(produced))
        for f in produced:
            print('   %s' % f)
        truth = [f for f in produced if f.endswith('_truth.txt')]
        if not truth:
            print('WARNING: no truth table was written')
            return 1
        n = sum(1 for ln in open(os.path.join(out_dir, truth[0]))
                if ln.strip() and ln.lstrip()[0] not in '#!')
        print('   truth table lists %d injected cores' % n)
        return 0
    finally:
        try:
            os.chdir(_cwd0)
        except Exception:
            pass
        if keep:
            print('\ntemporary directory kept at %s' % d)
        else:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main())
