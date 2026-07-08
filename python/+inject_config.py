#!/usr/bin/env python
"""
inject_config_all24.py -- all settings for injecting 23 unique Tier-2 models.

This is the single config file for one injection run. All paths and parameters
are here; inject_config_field.py is no longer needed or imported.

To run a new realization: change only RANDOM_SEED (SET_TAG updates automatically).

  python -B run_inject.py inject_config_all24
"""
# -----------------------------------------------------------------------
# Change only this line for each new run
# -----------------------------------------------------------------------
RANDOM_SEED = 1111
SET_TAG     = f'all24_s{RANDOM_SEED}'

# -----------------------------------------------------------------------
# N_CORES per model (field capacity limits):
#   Sigma=3e21  (~44k valid px): fine
#   Sigma=6e21  (~87k valid px): fine
#   Sigma=1.2e22 (~8k valid px): tight
#   Sigma=4.8e22  (~775 valid px): i5j2k1 fails (acceptable)
# N_CORES=2: places 44/46
# -----------------------------------------------------------------------
N_CORES = 2
R_BE_MAX_AS = 250.0   # arcsec; increase to e.g. 1400 for the full 3720x3720 field

# -----------------------------------------------------------------------
# Field scaling (1.0 = original; e.g. 8.0 shifts 3e21->2.4e22)
# -----------------------------------------------------------------------
SCALE_FACTOR = 8.0

# -----------------------------------------------------------------------
# Output directory ('.' = current directory in the terminal)
# -----------------------------------------------------------------------
OUT_DIR = '.'

# -----------------------------------------------------------------------
# Path to the BE sphere model grid catalog
# -----------------------------------------------------------------------
BES_CATALOG = '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260616_RT_BES_radmc3d/bes_model_params_catalog'

# -----------------------------------------------------------------------
# Herschel band images (full absolute paths)
# -----------------------------------------------------------------------
IMG_PATHS = {
    '070': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/260625_Herschel-zoom2/runs/results/+sources/+visuals/Aquila.s.070.obs.cb.fits',
    '160': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/260625_Herschel-zoom2/runs/results/+sources/+visuals/Aquila.s.160.obs.cb.fits',
    '250': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/260625_Herschel-zoom2/runs/results/+sources/+visuals/Aquila.s.250.obs.cb.fits',
    '350': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/260625_Herschel-zoom2/runs/results/+sources/+visuals/Aquila.s.350.obs.cb.fits',
    '500': '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/260625_Herschel-zoom2/runs/results/+sources/+visuals/Aquila.s.500.obs.cb.fits',
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
SIGMA_PATH       = '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/260625_Herschel-zoom2/runs/results/+sources/+visuals/Aquila.s.161.obs.cb.fits'
SIGMA_OMASK_PATH = '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/+Images/260625_Herschel-zoom2/+hires/+results/hi.surface.density.r13p5.omask.fits'

# -----------------------------------------------------------------------
# Hires surfdens image to inject column stamps into.
# run_inject.py will add each model's column stamp and write
# inj_SETTAG.surfdens.fits alongside the band images.
# hires is then run on the injected band images; this injected surfdens
# can be used as the starting reference or for comparison.
# Set to None to skip surfdens injection.
# -----------------------------------------------------------------------
HIRES_SURFDENS_PATH = '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/260625_Herschel-zoom2/runs/results/+sources/+visuals/Aquila.s.161.obs.cb.fits'

# -----------------------------------------------------------------------
# Source footprint avoidance mask
# Use the injected version (with 18 original cores masked out)
# -----------------------------------------------------------------------
FOOTS_PATH = '/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc/260625_Herschel-zoom2/runs/results/+sources/+visuals/sm.Aquila.s.161.obs.foots.fits'

# Models ordered largest-sep first, then by Sigma bin alternating, so the
# greedy placement never blocks a large model with a smaller one.
# Values from bes_model_params_catalog (SD_emb != 2.4e22).

STAMP_ROOT = '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260616_RT_BES_radmc3d'

MODELS = [
          'i06j01k01',
          'i06j01k02',
          'i06j01k03',
          'i06j02k01',
          'i06j02k02',
          'i06j02k03',
          'i06j03k01',
          'i06j03k02',
          'i06j03k03',
          'i06j04k01',
          'i06j04k02',
          'i06j04k03',
          'i06j05k01',
          'i06j05k02',
          'i06j05k03',
          'i06j06k01',
          'i06j06k02',
          'i06j06k03',
          'i06j07k01',
          'i06j07k02',
          'i06j07k03',
          'i06j07k04',
          'i06j07k05',
          'i06j07k06',
]
