#!/usr/bin/env python
"""inject_config_v3.py -- settings for run_inject_v3.py.

Run with:
    python -B run_inject_v3.py inject_config_v3

This is the population-matched run: the injected cores are chosen to reproduce
the joint distribution of background column density, contrast and observed size
of the real Aquila sources, as tabulated in joint_target_Aquila_v3.txt.

For the companion run with the background flattened under the footprints of the
injected models, use inject_config_v3_flat.py.
"""

# SET_TAG is NOT used when the models are drawn from the core mass function,
# which is the case here.  The output names are built from the sampling itself,
# as cmf<slope>_n<cores>_s<seed>, with _flat appended automatically when
# FLATTEN_BACKGROUND is on.  This run therefore writes files named
# cmf2.00_n600_s2029, and its flattened companion cmf2.00_n600_s2029_flat, so
# the two never overwrite each other and the tag needs no editing between them.

# ---------------------------------------------------------------------
# Population of injected cores
# ---------------------------------------------------------------------
# Reported masses are drawn from dN/dM proportional to M^-2 over the range
# below.  The mass that is actually injected follows from the grid node the
# observables select, and the truth table records the true mass M_BE, which is
# the only mass any later comparison should use.
CMF_SLOPE         = 2.0
CMF_N_CORES       = 600
CMF_MASS_RANGE    = (0.05, 6.0)      # M_sun; above 6 the grid has no node of
                                     # realistic observed size
JOINT_TARGET_FILE = 'joint_target_Aquila_v3.txt'
JOINT_NEIGHBOURS  = 40

# ---------------------------------------------------------------------
# Node score
# ---------------------------------------------------------------------
# The contrast and the observed size are measured directly by getsf and carry
# unit weight.  The reported mass comes from a spectral energy distribution fit
# whose bias this work exists to correct, so it is weighted down: matching it
# tightly would be fitting to a quantity known to be wrong.
W_MASS               = 0.1
W_CONTRAST           = 1.0
W_SIZE               = 1.0
MODEL_PICK_TOLERANCE = 0.005

# ---------------------------------------------------------------------
# Size limits and placement
# ---------------------------------------------------------------------
# Limits are on the observed size FWHMSDbs, not on the truncation radius R_BE.
# The lower limit is above the 13.5 arcsec surface-density beam.  The upper
# limit is near the 95th percentile of the real sources, 75.5 arcsec.  Do not
# raise the lower limit or lower the upper one without re-measuring: the
# faintest models on this grid are the large diffuse ones, so restricting the
# range to 15 to 55 arcsec forces the sampler into models that are too
# prominent and degrades the contrast match from a Kolmogorov-Smirnov
# statistic of 0.063 to 0.109.
# Which background-subtracted stamp is injected.  'bsl' removes the true
# crater background, the floor of the depression the model digs in the cloud it
# displaces, so the stamp is the model's own excess and is non-negative
# everywhere.  'bs' removes a flat background at the cloud level instead, and
# is negative in an annulus inside R_BE for every model whose outer envelope is
# rarer than the cloud, which is 15% of the grid and 56% of the highest
# embedding rung.  Earlier versions injected 'bs' on the argument that getsf
# subtracts an interpolated background by the same construction and the choice
# would cancel; it does not, because getsf interpolates over the real map, not
# over the model's embedding cloud.  The setting also selects which catalogue
# columns describe the injected core, ICSD<suffix> and M_SED4<suffix>, so that
# the model chosen is the model injected.
STAMP_SUFFIX         = 'bsl'

FWHM_MIN_AS          = 15.0
FWHM_MAX_AS          = 80.0
TARGET_SIZE_RANGE    = None          # use all real sources as the target
# Separation between injected cores.  MAX_FOOTPRINT_OVERLAP applies to EVERY
# pair: two footprints may not interpenetrate at all when it is zero,
# or by at most this fraction of the sum of their radii when it is positive.
#
# The value is set to reproduce the clustering of the real sources rather than
# to eliminate blending.  Real Aquila sub-field sources overlap the footprint
# of their nearest neighbour in 44% of cases, by more than 30% in 22% and by
# more than half in 10%; their nearest-neighbour separation has a median of
# 1.72 in units of the pair's mean half-maximum width and a 5th percentile of
# 0.78.  Forbidding overlap entirely gives injected cores a 5th percentile of
# 1.91, further apart than real cores typically are, so the test would be
# cleaner than the data it represents.  At 0.30 the injected distribution has a
# median of 1.63 and a 5th percentile of 1.36: the median clustering is
# reproduced, the extreme close pairs are not.  Blending is measured instead by
# the isolated control, inject_config_v3_iso.py, which is exactly paired.  It replaces the earlier size-class exemption, under which a
# compact model could be dropped onto an extended one with no separation
# requirement at all; measured on the first population-matched injection, that
# exemption left 33% of injected cores with a companion closer than the sum of
# their two half-maximum radii, and the error in the recovered mass correlated
# with crowding at +0.42.  Setting MAX_FOOTPRINT_OVERLAP to None restores the
# older behaviour, controlled by OVERLAP_FREE_FWHM_AS.
MAX_FOOTPRINT_OVERLAP = 0.30
FOOTPRINT_FACTOR      = 1.84
# Cores are separated by the radius enclosing FLUX_FRACTION of each model's own
# flux, measured from its azimuthally averaged profile in BES_PROFILES, and
# never less than the observed footprint radius.  Neither geometric radius
# serves: the observed footprint ignores the envelope of a concentrated model
# entirely, while R_BE, where the signal formally stops, is up to fifteen times
# the radius containing 95% of the flux, and enforcing it would exclude such
# models from a field containing real sources altogether.  Without the profiles
# catalogue the code falls back to the larger of R_BE and the observed
# footprint.
FLUX_FRACTION         = 0.95
# The same criterion is applied against the footprints of the REAL sources, by
# scaling the clearance with the footprint radius of the model being placed
# instead of using a fixed margin.  With a fixed 18.2 arcsec margin, 14% of
# injected cores overlapped a real footprint and the worst overlap was 44%;
# scaled, those become 2% and 4%.  A real source inside an injected footprint
# contributes flux the extraction cannot attribute correctly, which is the
# same contamination as two overlapping models.  Set to 'fixed' to restore the
# older behaviour.
FOOTPRINT_MARGIN_MODE = 'scaled'
OVERLAP_FREE_FWHM_AS = 50.0
SIGMA_WINDOW         = 'symmetric'
SEP_FACTOR           = 1.3
AVOID_REAL_FOOTPRINTS = True
FOOTPRINT_MARGIN_AS  = 18.2          # 500 micron beam radius
FLATTEN_BACKGROUND   = False

RANDOM_SEED = 2029
N_CORES     = 5
OUT_DIR     = '.'

# ---------------------------------------------------------------------
# Field paths -- edit to match the machine
# ---------------------------------------------------------------------
_IMG = ('/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc'
        '/+Images/260625_Herschel-zoom2')

IMG_PATHS = {
    '070': _IMG + '/aquilaM2-070.image.resamp.zoom2.fits',
    '160': _IMG + '/aquilaM2-160.image.resamp.zoom2+159p9.fits',
    '250': _IMG + '/aquilaM2-250.image.resamp.zoom2+169p7.fits',
    '350': _IMG + '/aquilaM2-350.image.resamp.zoom2+93p0.fits',
    '500': _IMG + '/aquilaM2-500.image.resamp.zoom2+37p0.fits',
}
OMASK_PATHS = {b: p.replace('.fits', '.omask.fits')
               for b, p in IMG_PATHS.items()}

SIGMA_PATH       = _IMG + '/+hires/+results/hi.surface.density.r13p5.fits'
SIGMA_OMASK_PATH = _IMG + '/+hires/+results/hi.surface.density.r13p5.omask.fits'
HIRES_SURFDENS_PATH = SIGMA_PATH

FOOTS_PATH = ('/Users/amenshch/Astronomy/+HERSCHEL_EXTRACTIONS/+AQUILA~260pc'
              '/260625_Herschel-zoom2/runs/results/+sources/+visuals'
              '/sm.Aquila.s.161.obs.foots.fits')

BES_PROFILES = None   # defaults to bes_model_grid_final2_profiles beside BES_CATALOG

STAMP_ROOT = '/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260721_RT_BES_radmc3d'
