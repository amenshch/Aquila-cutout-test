#!/usr/bin/env python
"""
inject_config_field.py -- shared field paths used by all injection runs.

Place this file anywhere on your PYTHONPATH (e.g. your scripts directory).
Each model config imports from it, so field paths are maintained in one place.
"""

# -----------------------------------------------------------------------
# Path to the BE sphere model grid catalog
# -----------------------------------------------------------------------
BES_CATALOG = '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260616_RT_BES_radmc3d/bes_model_params_catalog'

# -----------------------------------------------------------------------
# Herschel band images (full absolute paths)
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
# Hires surface-density map (used for Sigma matching during placement)
# and its binary coverage mask
# -----------------------------------------------------------------------
SIGMA_PATH       = '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260616_RT_BES_radmc3d/cSD_01/tBE_02/inj_surfdens.fits'
SIGMA_OMASK_PATH = '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/+hires/+results/hi.surface.density.r13p5.omask.fits'

# -----------------------------------------------------------------------
# Hires surfdens image to inject column stamps into.
# run_inject.py will add each model's column stamp and write
# inj_SETTAG.surfdens.fits alongside the band images.
# hires is then run on the injected band images; this injected surfdens
# can be used as the starting reference or for comparison.
# Set to None to skip surfdens injection.
# -----------------------------------------------------------------------
HIRES_SURFDENS_PATH = '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/+hires/+results/hi.surface.density.r13p5.fits'

# -----------------------------------------------------------------------
# Source footprint avoidance mask
# Use the injected version (with 18 original cores masked out)
# -----------------------------------------------------------------------
FOOTS_PATH = '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260616_RT_BES_radmc3d/cSD_01/tBE_02/sm.Aquila.s.161.obs.foots_injected_hiresSD.fits'

# -----------------------------------------------------------------------
# Injection settings shared across all models
# -----------------------------------------------------------------------
N_CORES        = 5
MIN_SEP_ARCSEC = 90.0
RANDOM_SEED    = 2026
OUT_DIR        = '.'     # current directory where you run the script
