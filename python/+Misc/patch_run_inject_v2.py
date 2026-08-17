#!/usr/bin/env python3
"""
patch_run_inject_v2.py -- write run_inject_v2.py from run_inject.py, adding
core-mass-function sampling, contrast-matched placement, and a switch for
background flattening.

    python3 patch_run_inject_v2.py run_inject.py

WHAT IS ADDED
-------------
Three configuration options, all optional.  With none of them set, v2 behaves
exactly as v1.

  CMF_SLOPE, CMF_MASS_RANGE, CMF_N_CORES
      Build the MODELS list by drawing CMF_N_CORES masses from
      dN/dM proportional to M^(-alpha) with alpha = CMF_SLOPE, between the two
      limits of CMF_MASS_RANGE in solar masses, and mapping each drawn mass onto
      the nearest model of the Bonnor-Ebert grid.  In the logarithmic convention
      dN/dlogM proportional to M^(-x) the exponent is x = alpha - 1, so the
      Salpeter value x = 1.35 corresponds to CMF_SLOPE = 2.35, and the shallower
      values discussed for core mass functions, x = 0.7 and x = 1.0, correspond
      to CMF_SLOPE = 1.7 and 2.0.

  SIGMA_TARGET_FILE
      A plain text file of background column densities in cm^-2, one per line,
      giving the distribution the injected cores should sample.  The natural
      choice is the PEAK^BGF values of the real sources in the field being
      imitated.  Each drawn mass is paired with a background column drawn from
      this distribution, and the grid model is then chosen to match BOTH the
      mass and that column.  Because valid_placement_mask() already restricts a
      model to pixels whose cloud column matches its own SD_emb, matching the
      model's SD_emb to the target distribution is what makes the placement
      contrast-matched; no change to the placement loop is needed.

      This matters.  In the unscaled sub-field run used so far, the injected
      cores spanned background columns of 1.1 to 1.9e22 cm^-2, whereas the real
      Aquila sources span 2.6e21 to 1.9e22, so 83% of real sources lay below
      anything the injections sampled and the injected median was 2.6 times the
      real one.

  FLATTEN_BACKGROUND
      True reproduces the earlier behaviour, interpolating the cloud background
      away under each core before adding the model.  False, the new default,
      adds the model on top of the cloud with its own fluctuations.

      False is the more informative choice.  Background fluctuations scatter the
      measured masses, and scatter in log-mass convolved with a falling mass
      function moves cores preferentially upward from the populous low-mass end,
      flattening the recovered slope.  The slope changes being tested are only
      about 0.1, so a fluctuation-induced flattening could be comparable to the
      signal.  Flattening also removes the very structural noise that the
      footprint definition, sigma(O_lambda; k_low = 1/D), exists to handle, so
      the recovered footprints would not be the ones real cores get.  Running
      one flattened realisation with the same seed and the same core list
      isolates the fluctuation contribution as a number.

MIN_SEP_ARCSEC now defaults to 36.3 arcsec, the 500 micron beam, when the config
does not set it, instead of being derived from the largest model radius.

CAVEAT
------
I could not execute this: the container has neither the field FITS maps nor the
inject module.  It is written to be reviewed before use, and the first
realisation should be checked before further runs are launched.
"""
import shutil, sys, os

BLOCK = '''
# Everything the v2 setup prints is mirrored here and written into the header of
# the truth table, so that a realisation documents the sampling that produced it.
_V2_LOG = []


def _vprint(msg):
    print(msg)
    _V2_LOG.append(msg)


# ===========================================================================
# v2: core-mass-function sampling and contrast-matched model selection
# ===========================================================================

def _draw_cmf_masses(n, alpha, m_lo, m_hi, rng):
    """Draw n masses from dN/dM proportional to M^(-alpha) on [m_lo, m_hi].

    Inverse transform of the truncated power law.  alpha = 1 is handled
    separately because the integral becomes logarithmic there.
    """
    u = rng.random(n)
    if abs(alpha - 1.0) < 1e-8:
        return m_lo * (m_hi / m_lo) ** u
    p = 1.0 - alpha
    return (m_lo ** p + u * (m_hi ** p - m_lo ** p)) ** (1.0 / p)


def _grid_column(name, n_expected):
    """Read one column of the Bonnor-Ebert grid catalogue by its header name.

    _GRID as returned by fit_mass.load_catalog carries only the basic model
    parameters, so quantities such as ICSDbs, the beam-convolved surface density
    at the model centre, are read from the catalogue file itself.  The column is
    located by the position of its name on the commented header line, so this
    does not depend on the column count.
    """
    try:
        lines = open(_bes_catalog_path).read().splitlines()
    except Exception:
        return np.full(n_expected, np.nan)
    hdr = None
    for ln in lines:
        if ln.lstrip()[:1] == '#' and name in ln.split():
            hdr = ln.lstrip('# ').split()
    if hdr is None or name not in hdr:
        print('  note: column %s not found in the grid catalogue; '
              'contrast matching disabled' % name)
        return np.full(n_expected, np.nan)
    j = hdr.index(name)
    out = []
    for ln in lines:
        if not ln.strip() or ln.lstrip()[:1] == '#':
            continue
        f = ln.split()
        out.append(float(f[j]) if j < len(f) else np.nan)
    out = np.array(out, dtype=float)
    return out if len(out) == n_expected else np.full(n_expected, np.nan)


def _grid_arrays():
    """(tags, M_BE, SD_emb, R_BE_as) for every node of the Bonnor-Ebert grid.

    _GRID is column-oriented: a dict of equal-length arrays keyed 'i', 'j', 'k',
    'SD_emb', 'T_BE', 'rho_BE', 'M_BE', 'R_BE_as', exactly as main() uses it.
    """
    i = np.asarray(_GRID['i'], dtype=int)
    j = np.asarray(_GRID['j'], dtype=int)
    k = np.asarray(_GRID['k'], dtype=int)
    tags = np.array(['i%02dj%02dk%02d' % (p, q, r) for p, q, r in zip(i, j, k)])
    mbe = np.asarray(_GRID['M_BE'], dtype=float)
    sde = np.asarray(_GRID['SD_emb'], dtype=float)
    rbe = np.asarray(_GRID['R_BE_as'], dtype=float)
    icsd = _grid_column('ICSDbs', len(mbe))
    if len(tags) == 0:
        raise ValueError('the Bonnor-Ebert grid catalogue produced no nodes; '
                         'check BES_CATALOG in the config file')
    return tags, mbe, sde, rbe, icsd


def _models_from_cmf(n_cores, alpha, m_lo, m_hi, sigma_target, rng,
                     r_be_max=None, r_be_min=None, contrast_target=None,
                     joint_target=None):
    """Build a MODELS-style list by drawing masses from a core mass function and
    background columns from a target distribution, then choosing for each pair
    the grid node closest in both.

    Distance is measured in log10 of each quantity, with the two given equal
    weight, so neither dominates by virtue of its dynamic range.
    """
    tags, mbe, sde, rbe, icsd = _grid_arrays()
    ok = np.isfinite(mbe) & (mbe > 0) & np.isfinite(sde) & (sde > 0)
    if r_be_min is not None:
        n_before = int(ok.sum())
        ok = ok & (rbe >= r_be_min)
        _vprint('  R_BE_MIN_AS = %.1f": %d of %d grid nodes retained'
              % (r_be_min, int(ok.sum()), n_before))
    if r_be_max is not None:
        n_before = int(ok.sum())
        ok = ok & (rbe <= r_be_max)
        _vprint('  R_BE_MAX_AS = %.1f": %d of %d grid nodes retained'
              % (r_be_max, int(ok.sum()), n_before))
        if not ok.any():
            raise ValueError('R_BE_MAX_AS excludes every grid node')
    tags, mbe, sde, icsd = tags[ok], mbe[ok], sde[ok], icsd[ok]
    lm, ls = np.log10(mbe), np.log10(sde)

    masses = _draw_cmf_masses(n_cores, alpha, m_lo, m_hi, rng)
    if sigma_target is None:
        sig = np.full(n_cores, np.nan)
    else:
        sig = rng.choice(sigma_target, size=n_cores, replace=True)

    # The contrast of a real core is not independent of its background: in
    # Aquila the median contrast rises from 1.18 below 4e21 cm^-2 to 1.64 above
    # 2e22, while the median footprint radius falls from 31 to 19 arcsec.
    # Drawing the two independently would impose the population-median contrast
    # at every column, which is wrong at both ends.  When the two target files
    # have equal length they are assumed to list the same sources in the same
    # order, and a single index is drawn for both, so the contrast is
    # automatically conditioned on the column.
    if joint_target is not None:
        # for each drawn mass, take the background and contrast of a real source
        # of similar mass, so that the coupling between the three is preserved
        jM, jS, jC = joint_target
        ljM = np.log10(jM)
        sig = np.empty(n_cores); con = np.empty(n_cores)
        for t in range(n_cores):
            near = np.argsort(np.abs(ljM - np.log10(masses[t])))[:40]
            pick = near[rng.integers(0, len(near))]
            sig[t] = jS[pick]; con[t] = jC[pick]
        _vprint('  background and contrast drawn from real sources of similar '
                'mass (%d rows in the joint target)' % len(jM))
    elif contrast_target is not None:
        if sigma_target is not None and len(contrast_target) == len(sigma_target):
            k = rng.integers(0, len(sigma_target), size=n_cores)
            sig = sigma_target[k]
            con = contrast_target[k]
            _vprint('  contrast drawn jointly with the background column '
                    '(the two target files are paired)')
        else:
            con = rng.choice(contrast_target, size=n_cores, replace=True)
            _vprint('  WARNING: target files differ in length, so contrast is '
                    'drawn independently of the background column')
    else:
        con = np.full(n_cores, np.nan)

    counts = {}
    picked = []
    tgt_of = {}
    chosen_m, chosen_s, chosen_c = [], [], []
    for M, Sg, Cn in zip(masses, sig, con):
        # A model may only be placed where the local column lies in
        # [SD_emb, sqrt(2)*SD_emb), because that is the column it was computed
        # embedded in.  So the background is a CONSTRAINT, not a term to
        # minimise, and the contrast must be evaluated at the column the core
        # will actually sit at rather than at the model's own SD_emb.  Scoring
        # the contrast at SD_emb, as earlier versions did, matched on a value
        # the core never has, and the error varied across the band, which made
        # the contrast fall with column within every model (Spearman -1) and
        # produced an overall correlation of -0.32 where the real cloud has
        # +0.51.  Constraining the band instead raises the two-dimensional
        # overlap with the real distribution from 0.32 to 0.61 and halves the
        # Kolmogorov-Smirnov statistic of the background marginal.
        cand = np.where((sde <= Sg) & (Sg < sde * np.sqrt(2.0)))[0] \
            if np.isfinite(Sg) else np.arange(len(mbe))
        if len(cand) == 0:
            cand = np.array([int(np.argmin(np.abs(np.log10(sde / Sg))))])
        d = (np.log10(M) - lm[cand]) ** 2
        if np.isfinite(Cn) and np.isfinite(Sg):
            # contrast a model of central column icsd would show on background Sg
            # getsf defines the source contrast as PEAK_SBF / PEAK_BGF, that is
            # (source + background) / background, so it exceeds unity.  The
            # model-side counterpart is (ICSDbs + SD_emb) / SD_emb, since
            # ICSDbs is the background-subtracted central column, the direct
            # analogue of PEAK_SRC.
            with np.errstate(all='ignore'):
                cgrid = 1.0 + icsd[cand] / Sg          # contrast AT the target column
                d = d + (np.log10(cgrid) - np.log10(max(Cn, 1.0 + 1e-6))) ** 2
        if MODEL_PICK_TOLERANCE > 0:
            w = np.exp(-(d - d.min()) / MODEL_PICK_TOLERANCE)
            ssum = w.sum()
            b = int(cand[rng.choice(len(d), p=w / ssum)]) if (np.isfinite(ssum)
                and ssum > 0) else int(cand[np.argmin(d)])
        else:
            b = int(cand[np.argmin(d)])
        counts[tags[b]] = counts.get(tags[b], 0) + 1
        picked.append(b)
        chosen_m.append(mbe[b]); chosen_s.append(sde[b])
        chosen_c.append(1.0 + icsd[b] / Sg if np.isfinite(Sg) else np.nan)
        tgt_of.setdefault(tags[b], []).append(Sg)
    chosen_m = np.array(chosen_m); chosen_s = np.array(chosen_s)
    chosen_c = np.array(chosen_c)

    _vprint('Core mass function sampling:')
    _vprint('  dN/dM proportional to M^(-%.2f)  '
          '(dN/dlogM proportional to M^(-%.2f))' % (alpha, alpha - 1.0))
    _vprint('  %d masses drawn on [%.3g, %.3g] Msun; median %.3g'
          % (n_cores, m_lo, m_hi, np.median(masses)))
    if sigma_target is not None:
        q = np.percentile(sigma_target, [5, 50, 95])
        _vprint('  background columns drawn from %d target values: '
              '5th %.2e, median %.2e, 95th %.2e cm^-2'
              % (len(sigma_target), q[0], q[1], q[2]))
    _vprint('  mapped onto %d distinct grid models' % len(counts))
    _vprint('  drawn vs mapped mass: median %.3g vs %.3g Msun, '
          'median |log10 ratio| %.3f'
          % (np.median(masses), np.median(chosen_m),
             np.median(np.abs(np.log10(chosen_m / masses)))))
    if sigma_target is not None:
        _vprint('  drawn vs mapped background column: median %.2e vs %.2e cm^-2, '
              'median |log10 ratio| %.3f'
              % (np.median(sig), np.median(chosen_s),
                 np.median(np.abs(np.log10(chosen_s / sig)))))
    if contrast_target is not None:
        cc = chosen_c[np.isfinite(chosen_c) & (chosen_c > 1.0)]
        _vprint('  target vs mapped contrast: median %.3f vs %.3f'
              % (np.median(contrast_target), np.median(cc)))
    # The mapping onto discrete grid nodes changes the mass distribution, so the
    # slope actually placed is not the nominal input.  Report it.
    for mmin in (0.1, 0.2):
        v = chosen_m[chosen_m >= mmin]
        if len(v) > 10:
            x = len(v) / np.sum(np.log(v / mmin))
            _vprint('  slope of the MAPPED masses above %.2f Msun: '
                  'x = %.2f +/- %.2f  (nominal input x = %.2f)'
                  % (mmin, x, x / np.sqrt(len(v)), alpha - 1.0))
    # the median column drawn for each model, so that placement can prefer
    # positions where the local column actually matches it
    return [dict(tag=t, n_cores=c,
                 sigma_target=float(np.nanmedian(tgt_of.get(t, [np.nan]))),
                 sigma_targets=list(tgt_of.get(t, [])))
            for t, c in sorted(counts.items())]


def _load_joint_target(path):
    """Read the three-column joint target: mass, background column, contrast."""
    if not path:
        return None
    M, S, C = [], [], []
    for line in open(path):
        line = line.strip()
        if not line or line[0] in '#!':
            continue
        f = line.split()
        if len(f) < 3:
            continue
        try:
            a, b, c = float(f[0]), float(f[1]), float(f[2])
        except ValueError:
            continue
        if a > 0 and b > 0 and c > 1:
            M.append(a); S.append(b); C.append(c)
    if not M:
        raise ValueError('no usable rows in %s' % path)
    return np.array(M), np.array(S), np.array(C)


def _load_sigma_target(path):
    """Read one background column density per line, in cm^-2."""
    if not path:
        return None
    v = []
    for line in open(path):
        line = line.strip()
        if not line or line[0] in '#!':
            continue
        try:
            x = float(line.split()[0])
        except ValueError:
            continue
        if np.isfinite(x) and x > 0:
            v.append(x)
    if not v:
        raise ValueError('no usable values in %s' % path)
    return np.array(v)

'''

CFG = """
# Everything the v2 setup prints is mirrored here and written into the header
# of the truth table, so that a realisation documents the sampling that
# produced it.
_V2_LOG = []


def _vprint(msg):
    print(msg)
    _V2_LOG.append(msg)


# --- v2 options -----------------------------------------------------------
# These are read from the config file, NOT from the environment.  Add them to
# inject_config.py; omit them and run_inject_v2.py behaves exactly as v1.
#
#   CMF_SLOPE = 2.0            # dN/dM ~ M^-alpha; alpha = x + 1, so x = 1.0
#   CMF_N_CORES = 500          # number of cores to draw
#   CMF_MASS_RANGE = (0.05, 20.0)   # solar masses
#   SIGMA_TARGET_FILE = 'sigma_target_Aquila.txt'
#   FLATTEN_BACKGROUND = False
#
# Defined here, after both the multi-model and single-model config branches, so
# that they exist whichever format the config file uses.
STAMP_ROOT         = globals().get('STAMP_ROOT', getattr(_cfg, 'STAMP_ROOT', None))
CMF_SLOPE          = getattr(_cfg, 'CMF_SLOPE', None)
CMF_MASS_RANGE     = getattr(_cfg, 'CMF_MASS_RANGE', (0.05, 20.0))
CMF_N_CORES        = getattr(_cfg, 'CMF_N_CORES', None)
SIGMA_TARGET_FILE  = getattr(_cfg, 'SIGMA_TARGET_FILE', None)
# CONTRAST_TARGET_FILE holds one value per line of the getsf source contrast,
# PEAK_SBF / PEAK_BGF = (source + background) / background, which exceeds unity,
# measured for real sources in the field being imitated.  When given, models are chosen
# to reproduce that distribution of contrast rather than only the distribution
# of background column: matching the background alone does not match the
# contrast, because the contrast also depends on the model's own peak.  In the
# first core-mass-function injection the injected cores had a median contrast of
# 1.731 against 1.235 for real Aquila sources, so their source component was
# three times too prominent relative to the background.
CONTRAST_TARGET_FILE = getattr(_cfg, 'CONTRAST_TARGET_FILE', None)
# R_BE_MIN_AS drops grid models smaller than this from the draw.  Cores below
# the 13.5 arcsec surface-density beam cannot be measured: in the first
# injection their recovered mass scattered by a factor of 8.8 and was biased
# high by up to a factor of 260, because the beam measures the surrounding
# cloud rather than the core.
R_BE_MIN_AS        = getattr(_cfg, 'R_BE_MIN_AS', None)
# SIGMA_TARGET_MIN restricts the drawn background columns to those at or above
# it.  Below about 5e21 cm^-2 the grid cannot reach the contrast that real cores
# have there: a core of 0.15 Msun and 15-50 arcsec would need a central column of
# only ~7e20 cm^-2 to show a contrast of 1.18 on a 4e21 background, and no
# Bonnor-Ebert node is that faint.  The sampler then accepts a large contrast
# error, giving injected cores of contrast ~2 wherever the column is low.  The
# same limit, 5e21, was used for model placement in the original Aquila
# completeness simulations.
SIGMA_TARGET_MIN   = getattr(_cfg, 'SIGMA_TARGET_MIN', None)
# JOINT_TARGET_FILE holds three columns per real source: reported mass,
# background column, and contrast.  When given, it supersedes
# SIGMA_TARGET_FILE and CONTRAST_TARGET_FILE: for each drawn mass the
# background and contrast are taken from a real source of SIMILAR mass, so the
# coupling between the three is preserved.  In Aquila that coupling is strong,
# with Spearman r(log M, log Sigma) = +0.73, and drawing them independently
# places cores of a given mass in regions far less dense than they occupy in
# reality, which makes them stand out too clearly and be measured too easily.
JOINT_TARGET_FILE  = getattr(_cfg, 'JOINT_TARGET_FILE', None)
# MODEL_PICK_TOLERANCE spreads the choice of grid node over the nodes that fit
# nearly as well as the closest one, instead of always taking the best match.
# Each node is scored by a squared distance in base-10 logarithms,
#     d = (log M - log M_node)^2 + (log Sigma - log Sigma_node)^2
#         + (log C - log C_node)^2,
# and is then chosen with probability proportional to exp(-d / TOL).  TOL
# therefore has the units of d, squared dex; it is a selection tolerance and has
# nothing to do with any physical temperature.  sqrt(TOL) is the mismatch, in
# dex, that halves a node's weight roughly: TOL = 0.005 corresponds to 0.07 dex,
# about 17% spread over the three quantities, while TOL = 0.15 corresponds to
# 0.39 dex, a factor of 2.4, which destroys the match.  TOL = 0 always takes the
# single best node.
#
# Some spreading is needed because with mass, background column and contrast all
# constrained at once the best-match rule collapses onto a few nodes: the s2028
# run drew 253 cores from only 34 distinct models, five of which supplied 63% of
# them.  On this grid the second-nearest node is materially worse, however, so
# the fidelity falls as TOL rises; 0.005 buys a third more diversity at
# negligible cost, and larger values are not worth it.
MODEL_PICK_TOLERANCE = float(getattr(_cfg, 'MODEL_PICK_TOLERANCE', 0.005))
# R_BE_MAX_AS was previously read only in the multi-model branch, so a
# single-model config silently ignored it.  Bound here for both formats.
R_BE_MAX_AS        = globals().get('R_BE_MAX_AS',
                                   getattr(_cfg, 'R_BE_MAX_AS', None))
STAMP_J_PREFIX     = getattr(_cfg, 'STAMP_J_PREFIX', STAMP_J_PREFIX)
# Separation between two placed cores is SEP_FACTOR * (r_i + r_j), with
# r = sqrt(R_BE^2 + BEAM500^2) the convolved radius.  SEP_FACTOR = 1.5 is the
# value calibrated empirically against the 500 micron footprints; lower it only
# after measuring the realised gaps with check_separation.py, since the gap
# scales with (r_i + r_j) and a global change affects large pairs most.
SEP_FACTOR         = float(getattr(_cfg, 'SEP_FACTOR', 1.5))
# AVOID_REAL_FOOTPRINTS controls whether the footprints of the real sources are
# forbidden to injected cores.  The injection field is the background image, in
# which every real core has already been interpolated away, so there is nothing
# to collide with and the exclusion only removes area: the footprints of the
# 1814 real Aquila sources cover 87% of the test sub-field, preferentially the
# filaments, which is where compact cores at high column density belong.  The
# footprint map itself is still read, because inject_model uses it for the
# background interpolation.
AVOID_REAL_FOOTPRINTS = bool(getattr(_cfg, 'AVOID_REAL_FOOTPRINTS', False))
# OVERLAP_FREE_R_BE_AS separates compact models from large ones for the purpose
# of the non-overlap criterion.  Requiring every pair to be separated is too
# strong: a large, faint, flat-topped model is barely detectable against the
# cloud and does not disturb the measurement of a small core lying on top of it,
# yet its exclusion zone denies that core the dense filament where it belongs.
# The criterion is therefore applied only within a class: two compact models may
# not overlap, and two large models may not overlap either, since intersecting
# flat tops would double the intensity over an extended area and create an
# artifact.  A compact model may be placed on a large one.  Set to None to
# restore the criterion for every pair.
OVERLAP_FREE_R_BE_AS = getattr(_cfg, 'OVERLAP_FREE_R_BE_AS', 50.0)
# FLATTEN_MAX_R_BE_AS restricts background flattening to models no larger than
# this, defaulting to the same threshold that governs overlap.  Interpolating
# the background under a large footprint is wrong on three counts: it removes a
# large area of genuine cloud structure, and so removes the very fluctuations
# that smaller cores overlapping that model are supposed to experience; the
# interpolation itself becomes progressively less accurate as the footprint
# grows, with an error of either sign; and it is very slow.  The same loss of
# accuracy affects getsf's own background subtraction for large footprints, and
# is a bias to keep in mind when interpreting recovered masses.
FLATTEN_MAX_R_BE_AS = getattr(_cfg, 'FLATTEN_MAX_R_BE_AS', OVERLAP_FREE_R_BE_AS)
FLATTEN_BACKGROUND = getattr(_cfg, 'FLATTEN_BACKGROUND', False)
# Restrict the band list to the bands the configuration actually provides.
# inject.BANDS is ['000', '070', '160', '250', '350', '500'], where '000' is the
# high-resolution surface-density map, supplied through HIRES_SURFDENS_PATH.  A
# configuration that omits a band previously failed with a bare KeyError deep in
# the file inventory; now the band is dropped with a printed note.  Set
# HIRES_SURFDENS_PATH in the config if the surface-density band is wanted.
_bands_have = [b for b in inject.BANDS
               if IMG_PATHS.get(b) and OMASK_PATHS.get(b)]
_bands_miss = [b for b in inject.BANDS if b not in _bands_have]
if _bands_miss:
    _vprint('  bands omitted (no image or observation mask in the config): %s'
          % ', '.join(_bands_miss))
    if '000' in _bands_miss:
        _vprint('    note: band 000 is the surface-density map that the correction')
        _vprint('    itself works in -- Sigma_cloud, the concentration zeta, the')
        _vprint('    footprint factor phi and the footprint radius are all measured')
        _vprint('    there.  Without it the injected cores are absent from that map')
        _vprint('    and the extraction will not exercise the correction.')
        _vprint('    SIGMA_PATH is used read-only, for placement validity and for')
        _vprint('    the local_Sigma of each core from the PRE-injection map.')
        _vprint('    To also inject into it and write an injected copy, add')
        _vprint('    to the config:')
        _vprint('        HIRES_SURFDENS_PATH = SIGMA_PATH')
        _vprint('    Loading the same file twice is intended: one copy stays')
        _vprint('    pristine as the placement and truth reference, the other')
        _vprint('    receives the cores.')
if not _bands_have:
    raise SystemExit('no usable bands: check IMG_PATHS and OMASK_PATHS')
inject.BANDS = _bands_have

if CMF_SLOPE is not None and CMF_N_CORES:
    # In core-mass-function mode the models are drawn, so a single MODEL_TAG
    # from the config would be misleading in every output filename.  The tag is
    # rebuilt from what actually defines the run, including the seed so that
    # realisations do not overwrite one another.
    SET_TAG = 'cmf%.2f_n%d_s%d' % (CMF_SLOPE, CMF_N_CORES, RANDOM_SEED)
    if FLATTEN_BACKGROUND:
        SET_TAG = SET_TAG + '_flat'
    if SCALE_FACTOR != 1.0:
        SET_TAG = SET_TAG + '_x%.4g' % SCALE_FACTOR
    _vprint('  v2: core mass function mode, dN/dM ~ M^-%.2f, %d cores, '
          'flatten=%s' % (CMF_SLOPE, CMF_N_CORES, FLATTEN_BACKGROUND))
    _vprint('  v2: output set tag "%s"' % SET_TAG)
"""


def main():
    src_path = sys.argv[1] if len(sys.argv) > 1 else 'run_inject.py'
    if not os.path.exists(src_path):
        sys.exit('file not found: %s' % src_path)
    out_path = os.path.join(os.path.dirname(os.path.abspath(src_path)),
                            'run_inject_v2.py')
    s = open(src_path).read()
    if 'CMF_SLOPE' in s:
        sys.exit('source already contains v2 additions; nothing done')

    # 1. helper block, before main()
    anchor = 'def main():'
    if anchor not in s:
        sys.exit('could not locate main(); this does not look like run_inject.py')
    s = s.replace(anchor, BLOCK.strip('\n') + '\n\n\n' + anchor, 1)

    # 2. configuration defaults, AFTER both config branches close, so that the
    #    names exist for the single-model format as well as the multi-model one
    marker = "    print(f'  Single-model: {SET_TAG}')"
    if marker not in s:
        sys.exit('could not locate the end of the config section; aborting')
    s = s.replace(marker, marker + '\n' + CFG.strip('\n') + '\n', 1)

    # 3. build MODELS from the core mass function, just before the inventory
    inv = "    print('File inventory:')"
    build = '''    _rng_cmf = np.random.default_rng(RANDOM_SEED)
    if CMF_SLOPE is not None and CMF_N_CORES:
        _sigma_target = _load_sigma_target(SIGMA_TARGET_FILE)
        _contrast_target = _load_sigma_target(CONTRAST_TARGET_FILE)
        if (SIGMA_TARGET_MIN is not None and _sigma_target is not None):
            _keep = _sigma_target >= SIGMA_TARGET_MIN
            _vprint('  SIGMA_TARGET_MIN = %.2e cm^-2: %d of %d target sources kept'
                    % (SIGMA_TARGET_MIN, int(_keep.sum()), len(_keep)))
            _sigma_target = _sigma_target[_keep]
            if (_contrast_target is not None
                    and len(_contrast_target) == len(_keep)):
                _contrast_target = _contrast_target[_keep]
        MODELS = [_normalise_model_entry(m, STAMP_ROOT)
                  for m in _models_from_cmf(CMF_N_CORES, CMF_SLOPE,
                                            CMF_MASS_RANGE[0], CMF_MASS_RANGE[1],
                                            _sigma_target, _rng_cmf,
                                            R_BE_MAX_AS, R_BE_MIN_AS,
                                            _contrast_target,
                                            _load_joint_target(JOINT_TARGET_FILE))]

'''
    s = s.replace(inv, build + inv, 1)

    # 2b. stamp subdirectory prefix: tBE_ renamed to M_, kept configurable
    old_sub = "subdir = os.path.join(root, f'cSD_{i:02d}', f'tBE_{j:02d}', f'{k:02d}')"
    new_sub = "subdir = os.path.join(root, f'cSD_{i:02d}', f'{STAMP_J_PREFIX}{j:02d}', f'{k:02d}')"
    if old_sub not in s:
        sys.exit('could not locate the stamp subdirectory expression; aborting')
    s = s.replace(old_sub, new_sub, 1)
    s = s.replace("        <root>/cSD_{i:02d}/tBE_{j:02d}/{k:02d}/nc.<band>um.bs.<res>x0.rs3p0as.fits",
                  "        <root>/cSD_{i:02d}/{STAMP_J_PREFIX}{j:02d}/{k:02d}/"
                  "nc.<band>um.bs.<res>x0.rs3p0as.fits\n"
                  "    STAMP_J_PREFIX defaults to 'M_' (formerly 'tBE_') and may be\n"
                  "    overridden in the config file.", 1)
    # define the constant before the function that uses it
    anchor_c = "def _stamp_files_from_tag(tag, root):"
    s = s.replace(anchor_c,
                  "STAMP_J_PREFIX = 'M_'   # middle level of the stamp tree; was 'tBE_'\n\n\n"
                  + anchor_c, 1)

    # 3. default grid catalogue name
    s = s.replace("_find_catalog('bes_model_params_catalog')",
                  "_find_catalog('bes_model_grid_final2_catalog')", 1)

    # 4. separation default: the 500 micron beam
    old_sep = "        R_max = max(md['r_be_arcsec'] for md in model_data)"
    new_sep = ("        # v2: the 500 micron beam is sufficient to keep convolved\n"
               "        # models from overlapping, and is far less restrictive than\n"
               "        # a criterion scaled to the largest model radius.\n"
               "        _min_sep = BEAM500\n"
               "        print('  separation criterion = %.1f\" (500 um beam)' % _min_sep)\n"
               "    if False:\n"
               "        R_max = max(md['r_be_arcsec'] for md in model_data)")
    if old_sep in s:
        s = s.replace(old_sep, new_sep, 1)

    # 4b. pairwise separation.  The original test compared every candidate with
    #     every placed centre using the separation of the model being placed,
    #     1.5 * 2 * sqrt(R_BE^2 + BEAM500^2).  Because models are placed largest
    #     first, a small core placed later was tested against its own small
    #     separation, so it could sit on top of a much larger neighbour.  The
    #     criterion is replaced by the symmetric 1.5 * (r_i + r_j), with
    #     r = sqrt(R_BE^2 + BEAM500^2) the convolved radius.  For two models of
    #     the same size this is identical to the original expression, so the
    #     calibration against the 500 um footprints is unchanged; only pairs of
    #     different size are affected.
    old_sep_calc = """        import math
        _sep_this_pix = int(math.ceil(
            1.5 * 2.0 * math.sqrt(md['r_be_arcsec']**2 + BEAM500**2) / inject.PIX_ARCSEC)) + 1"""
    new_sep_calc = """        import math
        _r_conv = math.sqrt(md['r_be_arcsec']**2 + BEAM500**2)   # arcsec"""
    if old_sep_calc not in s:
        sys.exit('could not locate the per-model separation expression; aborting')
    s = s.replace(old_sep_calc, new_sep_calc, 1)

    old_test = """            s2 = _sep_this_pix ** 2
            if all((y-cy)**2 + (x-cx)**2 >= s2 for cy,cx in used_centers):
                used_centers.append((y, x))"""
    new_test = """            if len(used_centers):
                _uc = np.asarray(used_centers, dtype=float)
                _smin = (SEP_FACTOR * (_r_conv + _uc[:, 2])
                         / inject.PIX_ARCSEC) + 1.0
                _ok = np.all((y - _uc[:, 0])**2 + (x - _uc[:, 1])**2
                             >= _smin * _smin)
            else:
                _ok = True
            if _ok:
                used_centers.append((y, x, _r_conv))"""
    if old_test not in s:
        sys.exit('could not locate the separation test; aborting')
    s = s.replace(old_test, new_test, 1)
    s = s.replace("    used_centers  = []   # all placed centers for global separation check",
                  "    used_centers  = []   # (y, x, convolved radius in arcsec) per placed core",
                  1)

    # 4c. record the sampling diagnostics in the truth-table header, so that
    #     a realisation documents the sampling that produced it
    key = 'Injection truth table'
    src_lines = s.split(chr(10))
    hit = [k for k, ln in enumerate(src_lines) if key in ln and 'f.write' in ln]
    if not hit:
        sys.exit('could not locate the truth-table header line; aborting')
    k = hit[0]
    src_lines.insert(k + 1, '        for _ln in _V2_LOG:')
    src_lines.insert(k + 2,
                     '            f.write(' + chr(39) + '# %s' + chr(92) +
                     'n' + chr(39) + ' % _ln.strip())')
    s = chr(10).join(src_lines)
    # 4d. the footprint map still loads, but forbids placement only if asked
    old_forb = "    forbidden = foots > 0"
    new_forb = ("    # the injection field has every real core interpolated away, so the\n"
                "    # footprints of the former sources need not be forbidden; see\n"
                "    # AVOID_REAL_FOOTPRINTS above\n"
                "    forbidden = (foots > 0) if AVOID_REAL_FOOTPRINTS \\\n"
                "        else np.zeros(foots.shape, dtype=bool)\n"
                "    if not AVOID_REAL_FOOTPRINTS:\n"
                "        _vprint('  real-source footprints are NOT excluded from placement '\n"
                "                '(%.0f%% of the field would otherwise be lost)'\n"
                "                % (100.0 * np.mean(foots > 0)))")
    if old_forb not in s:
        sys.exit('could not locate the forbidden-mask line; aborting')
    s = s.replace(old_forb, new_forb, 1)

    # 4e. size-class overlap rule
    old_test2 = """            if len(used_centers):
                _uc = np.asarray(used_centers, dtype=float)
                _smin = (SEP_FACTOR * (_r_conv + _uc[:, 2])
                         / inject.PIX_ARCSEC) + 1.0
                _ok = np.all((y - _uc[:, 0])**2 + (x - _uc[:, 1])**2
                             >= _smin * _smin)
            else:
                _ok = True
            if _ok:
                used_centers.append((y, x, _r_conv))"""
    new_test2 = """            if len(used_centers):
                _uc = np.asarray(used_centers, dtype=float)
                _smin = (SEP_FACTOR * (_r_conv + _uc[:, 2])
                         / inject.PIX_ARCSEC) + 1.0
                _d2 = (y - _uc[:, 0])**2 + (x - _uc[:, 1])**2
                if OVERLAP_FREE_R_BE_AS is None:
                    _apply = np.ones(len(_uc), dtype=bool)
                else:
                    # the criterion binds only within a size class
                    _this_small = md['r_be_arcsec'] <= OVERLAP_FREE_R_BE_AS
                    _other_small = _uc[:, 3] <= OVERLAP_FREE_R_BE_AS
                    _apply = (_other_small == _this_small)
                _ok = np.all(_d2[_apply] >= (_smin * _smin)[_apply])
            else:
                _ok = True
            if _ok:
                used_centers.append((y, x, _r_conv, md['r_be_arcsec']))"""
    if old_test2 not in s:
        sys.exit('could not locate the pairwise separation test; aborting')
    s = s.replace(old_test2, new_test2, 1)
    s = s.replace("    used_centers  = []   # (y, x, convolved radius in arcsec) per placed core",
                  "    used_centers  = []   # (y, x, convolved radius, R_BE) per placed core", 1)

    # 4f. the forbidden disk must obey the same size-class rule as the pairwise
    #     test.  Previously every placed core stamped a disk of radius
    #     stamp_halfwidth, the model image half-width, which encloses emission
    #     far below anything detectable and which valid_placement_mask then
    #     refused unconditionally.  That second mechanism contradicted the
    #     pairwise criterion: a compact core allowed to sit on a large flat model
    #     could still only reach its rim.  Now a large model stamps no disk at
    #     all, relying on the pairwise test to keep large models apart from one
    #     another, and a compact model stamps SEP_FACTOR * r, the same radius the
    #     pairwise test uses, rather than the full stamp.
    old_disk = """                # also mark as forbidden for subsequent models
                r = md['stamp_halfwidth']"""
    new_disk = """                # mark as forbidden for subsequent models, on the same
                # size-class rule as the pairwise test above
                if (OVERLAP_FREE_R_BE_AS is not None
                        and md['r_be_arcsec'] > OVERLAP_FREE_R_BE_AS):
                    r = 0
                else:
                    r = int(round(SEP_FACTOR * _r_conv / inject.PIX_ARCSEC))"""
    if old_disk not in s:
        sys.exit('could not locate the forbidden-disk radius; aborting')
    s = s.replace(old_disk, new_disk, 1)
    old_stamp = """                y0,y1 = max(0,y-r), min(field['shape'][0],y+r+1)"""
    new_stamp = """                if r > 0:
                    y0,y1 = max(0,y-r), min(field['shape'][0],y+r+1)"""
    if old_stamp not in s:
        sys.exit('could not locate the disk-stamping block; aborting')
    s = s.replace(old_stamp, new_stamp, 1)
    for ln in ("                x0,x1 = max(0,x-r), min(field['shape'][1],x+r+1)",
               "                yy,xx  = np.ogrid[y0-y:y1-y, x0-x:x1-x]",
               "                forbidden[y0:y1, x0:x1] |= (xx**2+yy**2 < r**2)"):
        s = s.replace(ln, "    " + ln, 1)

    # 4g. flatten only models small enough for the interpolation to be sound,
    #     and record local_Sigma from the map as it stands at placement, which
    #     is what the core actually sits on once earlier models and any
    #     flattening have been applied
    old_flat = "            flatten=True,"
    new_flat = ("            flatten=(FLATTEN_BACKGROUND\n"
                "                     and (FLATTEN_MAX_R_BE_AS is None\n"
                "                          or md['r_be_arcsec'] <= FLATTEN_MAX_R_BE_AS)),")
    if old_flat not in s:
        sys.exit('could not locate the flatten argument; aborting')
    s = s.replace(old_flat, new_flat, 1)

    old_sig = "        sigma_patch = field['sigma'][y0s:y1s, x0s:x1s]"
    new_sig = ("        # the surface-density map as it stands now, if that band is being\n"
               "        # injected into, so that local_Sigma describes the background this\n"
               "        # core really sits on rather than the pristine field\n"
               "        _sig_map = current.get('000', field['sigma'])\n"
               "        if _sig_map.shape != field['sigma'].shape:\n"
               "            _sig_map = field['sigma']\n"
               "        sigma_patch = _sig_map[y0s:y1s, x0s:x1s]")
    if old_sig not in s:
        sys.exit('could not locate the local_Sigma patch; aborting')
    s = s.replace(old_sig, new_sig, 1)
    s = s.replace("        local_sig   = float(sigma_patch[valid_px].mean()) "
                  "if valid_px.any() else field['sigma'][cy, cx]",
                  "        local_sig   = float(sigma_patch[valid_px].mean()) "
                  "if valid_px.any() else _sig_map[cy, cx]", 1)

    # report the model peak, not the net change to the map
    old_pk = "        peak_added = patch.max()"
    new_pk = ("        # with flattening, patch.max() is the NET change to the map, which\n"
              "        # goes negative where the interpolation removed more than the model\n"
              "        # added; report the model's own peak instead\n"
              "        peak_added = float(np.nanmax(md['Ibs']['250'])) \\\n"
              "            if isinstance(md.get('Ibs'), dict) and '250' in md['Ibs'] \\\n"
              "            else patch.max()")
    if old_pk in s:
        s = s.replace(old_pk, new_pk, 1)

    # 4h. order the candidate positions by how well the local column matches the
    #     column that was drawn for this model, instead of at random.
    #     valid_placement_mask accepts any position whose local column lies
    #     between SD_emb and sqrt(2)*SD_emb.  Taking positions at random from
    #     that band means a model's realised contrast, 1 + IC/local_Sigma with IC
    #     fixed, varies inversely with where it happened to land: within each of
    #     the five most-used models of the s2028 run the Spearman correlation
    #     between local column and contrast was exactly -1.  Pooled over models
    #     that manufactured an overall correlation of -0.32, where the real cloud
    #     has +0.51 and the between-model value is +0.02.  Preferring the
    #     positions closest to the drawn column removes the artifact and also
    #     tightens the background marginal, whose median ran 1.2 to 1.3 times
    #     the drawn value for the same reason.
    old_order = "        order  = rng.permutation(len(ys))"
    new_order = ("        # prefer positions whose local column is closest to the one drawn\n"
                 "        _tgt = md.get('sigma_target', md.get('sd_emb'))\n"
                 "        if _tgt and _tgt > 0:\n"
                 "            _lsig = field['sigma'][ys, xs]\n"
                 "            _key = np.where(_lsig > 0,\n"
                 "                            np.abs(np.log10(np.maximum(_lsig, 1e-30)\n"
                 "                                            / _tgt)), 9.0)\n"
                 "            _key = _key + 1e-6 * rng.random(len(_key))\n"
                 "            order = np.argsort(_key)\n"
                 "        else:\n"
                 "            order = rng.permutation(len(ys))")
    if old_order not in s:
        sys.exit('could not locate the candidate ordering; aborting')
    s = s.replace(old_order, new_order, 1)

    # 5. background flattening switch
    old_fl = '            flatten=True,'
    if old_fl in s:
        s = s.replace(old_fl, '            flatten=FLATTEN_BACKGROUND,', 1)
    elif 'FLATTEN_BACKGROUND' not in s:
        sys.exit('could not locate the flatten=True argument; aborting')

    open(out_path, 'w').write(s)
    print('wrote     : %s' % out_path)
    print('added     : CMF sampling, contrast-matched model choice, '
          'FLATTEN_BACKGROUND')
    print('changed   : default separation criterion to the 500 um beam (36.3")')
    print('unchanged : behaviour when CMF_SLOPE is not set in the config')
    import py_compile
    py_compile.compile(out_path, doraise=True)
    print('syntax    : OK')


if __name__ == '__main__':
    main()
