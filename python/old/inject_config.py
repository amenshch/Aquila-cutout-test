#!/usr/bin/env python
"""
inject_config.py -- user settings for run_inject.py

This is the ONLY file you need to edit between injection runs.
run_inject.py reads it automatically -- never edit run_inject.py itself.

To inject a different model:
  1. Change MODEL_TAG, SD_EMB, R_BE_ARCSEC, and STAMP_FILES.
  2. Run:  python run_inject.py
"""
# -----------------------------------------------------------------------
# Model parameters (from bes_model_params_catalog for this model)
# -----------------------------------------------------------------------

MODEL_TAG    = 'i1j2k4'   # label used in output filenames
SD_EMB       = 3.0e21     # cloud column density matching this model (H2/cm^2)
R_BE_ARCSEC  = 17.97      # outer BE sphere radius (arcsec)

# -----------------------------------------------------------------------
# Injection settings
# -----------------------------------------------------------------------

CMF_SLOPE            = 2.0
CMF_N_CORES          = 600
CMF_MASS_RANGE       = (0.15, 3.0)
##SIGMA_TARGET_FILE    = 'sigma_target_Aquila.txt'
##CONTRAST_TARGET_FILE = 'contrast_target_Aquila.txt'
JOINT_TARGET_FILE = 'joint_target_Aquila.txt'
R_BE_MIN_AS          = 15.0
R_BE_MAX_AS          = 70.0
SEP_FACTOR           = 1.3
OVERLAP_FREE_R_BE_AS = 50.0
FLATTEN_MAX_R_BE_AS  = 50.0
FLATTEN_BACKGROUND   = False
AVOID_REAL_FOOTPRINTS = False

RANDOM_SEED          = 2028     # set to any integer; same seed -> same positions

N_CORES        = 5        # number of positions to inject
MIN_SEP_ARCSEC = 90.0     # minimum centre-to-centre separation (arcsec)

# -----------------------------------------------------------------------
# Output directory
# '.' = current directory open in the terminal when you run the script
# -----------------------------------------------------------------------

OUT_DIR = '.'

# -----------------------------------------------------------------------
# Herschel band images (full paths on your machine)
# -----------------------------------------------------------------------

IMG_PATHS = {
    '070': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/aquilaM2-070.image.resamp.zoom2.fits',
    '160': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/aquilaM2-160.image.resamp.zoom2+159p9.fits',
    '250': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/aquilaM2-250.image.resamp.zoom2+169p7.fits',
    '350': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/aquilaM2-350.image.resamp.zoom2+93p0.fits',
    '500': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/aquilaM2-500.image.resamp.zoom2+37p0.fits',
}

# -----------------------------------------------------------------------
# Coverage masks (1 = observed, 0 = no data)
# -----------------------------------------------------------------------

OMASK_PATHS = {
    '070': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/aquilaM2-070.image.resamp.zoom2.omask.fits',
    '160': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/aquilaM2-160.image.resamp.zoom2+159p9.omask.fits',
    '250': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/aquilaM2-250.image.resamp.zoom2+169p7.omask.fits',
    '350': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/aquilaM2-350.image.resamp.zoom2+93p0.omask.fits',
    '500': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/aquilaM2-500.image.resamp.zoom2+37p0.omask.fits',
}

# -----------------------------------------------------------------------
# Hires surface-density map and its coverage mask
# -----------------------------------------------------------------------

SIGMA_PATH       = '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/+hires/+results/hi.surface.density.r13p5.fits'
SIGMA_OMASK_PATH = '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/+hires/+results/hi.surface.density.r13p5.omask.fits'

# -----------------------------------------------------------------------
# Source footprint avoidance mask (getsf output: nonzero = real source)
# -----------------------------------------------------------------------

FOOTS_PATH = '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/260625_Herschel-zoom2/runs/results/+sources/+visuals/sm.Aquila.s.161.obs.foots.fits'

# -----------------------------------------------------------------------
# Model stamp files (background-subtracted, MJy/sr, centered, odd size)
# Set a band to None if no stamp exists for it (zeros will be used).
# -----------------------------------------------------------------------

STAMP_ROOT = '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260721_RT_BES_radmc3d'
HIRES_SURFDENS_PATH = SIGMA_PATH
STAMP_FILES = {
    '070': None,                                    # no 70um stamp for this model
    '160': '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260721_RT_BES_radmc3d/cSD_01/M_02/04/nc.160um.bs.r13p5x0.rs3p0as.fits',
    '250': '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260721_RT_BES_radmc3d/cSD_01/M_02/04/nc.250um.bs.r18p2x0.rs3p0as.fits',
    '350': '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260721_RT_BES_radmc3d/cSD_01/M_02/04/nc.350um.bs.r24p9x0.rs3p0as.fits',
    '500': '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260721_RT_BES_radmc3d/cSD_01/M_02/04/nc.500um.bs.r36p3x0.rs3p0as.fits',
}
