#!/usr/bin/env python3
"""
mass_correction_pipeline_4d.py  (v4 -- self-contained, internal SED fitting;
adds the four-observable, physical-size correction alongside the original
three-observable one; see NOTATION and the cm_size_interp() docstring below)
=============================================================================
Corrects Herschel SED-derived prestellar core masses for the systematic
under-estimate measured against a radiative-transfer (radmc-3d),
pressure-truncated Bonnor-Ebert model grid, and reproduces the paper's
leave-one-out validation and seven-cloud Table 1.

NOTATION: this script's internal variable/column names (conc, R_foot, FWHM,
AFWHM, FOOA, FOOB, RT, ...) predate a later notation pass in the paper and
have NOT been renamed here, to avoid breaking file-format compatibility
with existing catalogues and downstream scripts. The mapping to the paper's
symbols, used throughout its text and figures, is:
  conc / C_pm      ->  zeta        (peak-to-mean concentration)
  R_foot           ->  phi         (footprint factor, A_F/H)
  FWHM, AFWHM      ->  H           (half-maximum size)
  FOOA, FOOB       ->  A_F, B_F    (footprint major/minor axis)
  SD (surface-density proxy band)  ->  Sigma
  RT, RADMC-3D     ->  radiative transfer, radmc-3d
Read code comments below with this table in hand; a future pass may rename
the identifiers themselves once the paper's notation has converged.

WHAT THIS SCRIPT IS
-------------------
ONE file.  Depends only on numpy + scipy.  It imports NO other project module
(not fit_mass.py, getsf_catalog.py, tabulate_recovery.py, ...): the grid loader,
the getsf reader, the SED fitter, the correction-factor grid, the mass grid,
the leave-one-out validation, and the seven-cloud table are all here.

WHAT YOU NEED (inputs)
----------------------
Two model-grid files (produced once, distance-independent -- the same for every
cloud):
  1. GRID CATALOGUE      cats/bes_model_grid_final2_catalog
                         (518 nodes, 67 columns: per-band model fluxes + M_BE)
  2. RECOVERY TABLES     cats/bes_model_grid_final2_recovery_tables_v2
                         (per-node, per-band, per-truncation-radius flux recovery)
  (A third file, ..._catalog_rec, is used only by the `loo` validation step.)

Per cloud, ONE getsf catalogue.  Either form works (auto-detected per cloud):
  * the concatenated getsf+SED file  <Cloud>.sw.sources.ok.cat=...=thin...00.cat,
    OR
  * the two raw getsf files supplied together, the extraction catalogue
    <Cloud>.sw.sources.ok.cat and the additional-quantities catalogue
    <Cloud>.sw.sources.ok.add.cat  (same sources, same order, complementary
    columns; '!' comments and '_' column names are handled).  Pass both to
    load_getsf([cat, add]).
In BOTH cases the SED mass is computed INSIDE this script (MASS_SOURCE='internal')
by refitting the same modified-blackbody ("thinbody") SED that getsf/fitfluxes
uses -- so the thin.*.00 SED catalogue is NOT needed.  For matched sources the
internal mass equals the catalogue TOTL^MASS exactly.  Source selection differs
by input: the concatenated file carries the thin-catalogue SED-fit QUALITY flag
and is gated on QUALITY=='ok'; raw .cat/.add.cat have no such flag, so the getsf
header's recommended per-band criteria (|GOODM|>1, |SIGNM|>1, FXP/FXP_ERR>2,
FXT/FXT_ERR>2, AFWHM/BFWHM<2, FOOA/AFWHM>1.15, on the surface-density band) are
applied instead -- a different but each-appropriate quality gate, so the two
inputs give C_M within ~0.1 but not an identical source count.  Set
MASS_SOURCE='catalog' to read TOTL^MASS instead (concatenated file only).

HOW TO RUN
----------
Single cloud (writes one per-source corrected catalogue):
  python mass_correction_pipeline.py correct <catalogue> --distance D [--out F]
  python mass_correction_pipeline.py correct <cat> <add> --distance D [--out F]
Model-grid validation (no cloud needed):
  python mass_correction_pipeline.py loo       # mass-corrector node-LOO
  python mass_correction_pipeline.py cmloo     # correction-factor node-LOO
  python mass_correction_pipeline.py test <catalogue> <recovery_tables>
                                               # apply the main-grid correction to
                                               # an INDEPENDENT set of grid models
                                               # (e.g. a subdivision grid) and
                                               # compare corrected mass to true M_BE
Many clouds at once: use the separate driver  run_clouds.py, which calls
correct_cloud() per field and writes a run-summary table.

OUTPUT (per cloud): an aligned, getsf-style catalogue with a '!' header, one row
per corrected QUALITY==ok source in the model hull, carrying the measured
observables, the INITIAL (reported) SED fit (T, M, errors, chi2, dof, QUALITY),
and the correction factor C_M decomposed into three multiplicative mass factors
  C_M = F_FLUX (footprint flux recovery)
      x F_BG   (background subtraction, FXSDbsl/FXSDbs: interpolated -> true crater)
      x F_TEMP (temperature-gradient bias),
followed by the corrected true mass M_CORR = M_INIT * C_M.

METHOD SUMMARY
--------------
The correction factor C_M = M_BE/M_reported is predicted on three
distance-invariant getsf observables measured on the 13.5" H2 (surface-density,
Sigma in the paper) band: the background column Sigma = PEAK^BGF, the
peak-to-mean concentration (zeta in the paper) C_pm = PEAK^SRC/(FXT_BST/
footprint_area), and the footprint factor (phi in the paper)
R_foot = sqrt(FOOA*FOOB)/sqrt(AFWHM*BFWHM). The reported mass is not an axis, so
the same grid corrects every cloud regardless of distance. The grid is built
from the model recovery tables: each SD-band non-absorption row gives one
(Sigma, C_pm, R_foot) -> C_M sample, where the reported mass M_rec is obtained
by reducing each band's model flux by its own recovered fraction f_b(r), refitting
the SED, and taking the amplitude ratio (per-band, colour-distorting loss -- not
a grey rescaling).

NOT YET IMPLEMENTED HERE: a later notation/validation pass (see the paper)
found that these three observables leave a genuine, non-shrinking degeneracy
in C_M -- sources sharing all three can still need correction factors
differing by 2x or more -- and that adding the model's physical truncation
radius R_BE (parsecs) as a fourth axis, matched to a real source's footprint
radius sqrt(FOOA*FOOB)/2 at its KNOWN distance, removes most of it. This is
validated in three independent tests (injection, an independent subdivision
grid, and node-level leave-one-out) but is not yet built into cm_interp() /
build_grids() / correct_cloud() below, which remain the three-axis version;
a fourth-axis corrector was only exercised via standalone test scripts.
Note this reintroduces a genuine distance dependence (unlike the pure
three-axis design), which the paper discusses explicitly.

SED FITTER (validated against getsf/fitfluxes 'thinbody')
---------------------------------------------------------
Model  F_nu = B_nu(T) * kappa(lambda) * M / D^2,  kappa ~ lambda^(-beta), beta=2.
Fit (T, amplitude) to the getsf total fluxes FXT_BST at 160/250/350/500 um with
per-band errors sqrt(FXT_ERR^2 + (adderr*F)^2), adderr = (0.2,0.1,0.1,0.1),
reproducing fitfluxes' per-band selection (snrmin=1): a band enters the fit only
if FXP_BST, FXP_ERR, FXT_BST, FXT_ERR are all > 0, |SIGNM| > 1,
FXP_BST/FXP_ERR >= 1, and FXT_BST/FXT_ERR > 1 (two bands minimum).  The mass
scale is one cloud-independent constant OPACITY_K (kappa_0, reference wavelength,
beta=2, dust-to-gas).  Matches the catalogue DUST^TEMP to <0.1 K and TOTL^MASS to
<5% for 97-99.8% of 'ok'-quality sources across all seven clouds; N_corr and the
C_M column are identical whether the internal or the catalogue mass is used.
=============================================================================
"""
import argparse, glob, os, re, sys, time
from collections import defaultdict
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree

# ------------------------------- CONFIG -------------------------------------
CAT67  = 'cats/bes_model_grid_final2_catalog'          # grid, 67-col (model fluxes + M_BE)
CATREC = 'cats/bes_model_grid_final2_catalog_rec'      # grid, recoverable cols (loo step)
RECTAB = 'cats/bes_model_grid_final2_recovery_tables_v2'
CLOUD_DIR = '{cloud}-Guoyin'    # per-cloud getsf directory

# ------------------------------- VERBOSITY -----------------------------------
# VERBOSE=0: main steps only (what's running, when it's done, headline counts).
# VERBOSE=1: + deeper comments (per-stage sample/node counts, filtering stats).
# VERBOSE=2: + most detailed (per-node/per-source progress, timing of substeps).
# Set via the --verbose N CLI flag, or directly: mc.VERBOSE = 2 when used as a
# library. vlog() is the single choke point all pipeline stages print through.
VERBOSE = 0

def vlog(level, *args):
    if VERBOSE >= level:
        print(*args, flush=True)

class _Timer:
    """Small helper for verbose=2 timing messages: `with _Timer(2,'step'):`."""
    def __init__(self, level, label):
        self.level, self.label = level, label
    def __enter__(self):
        self.t0 = time.time()
        vlog(self.level, "  [%s] starting..." % self.label)
        return self
    def __exit__(self, *exc):
        vlog(self.level, "  [%s] done in %.1fs" % (self.label, time.time() - self.t0))


def cloud_input(cloud):
    """Return the getsf input for a cloud: the concatenated catalogue if present,
    else the raw [.cat, .add.cat] pair."""
    d = CLOUD_DIR.format(cloud=cloud)
    concat = glob.glob(d + '/*=*=*')
    if concat:
        return concat[0]
    cat = [f for f in glob.glob(d + '/*.sw.sources.ok.cat') if '=' not in f]
    add = glob.glob(d + '/*.sw.sources.ok.add.cat')
    if cat and add:
        return [cat[0], add[0]]
    raise FileNotFoundError('no getsf catalogue found in ' + d)
CLOUDS = [('Scorpius', 130), ('Ophiuchus', 144), ('Aquila', 260),
          ('OrionA', 432), ('California', 470), ('CygnusX', 1150),
          ('W3W4W5', 1700)]
FRAC_FLOOR = 1e-4
MATCHED_RANGE = (0.1, 2.0)
CMLOO_STRIDE = 1                 # 1 = full node set (~5 min); >1 subsamples
CLAMP_CM = False                  # physical floor: every mechanism only ADDS mass, so
                                 # C_M and each factor (F_flux,F_bg,F_temp) >= 1.  Stops
                                 # the corrector from reducing a mass where little/no
                                 # correction is needed (LinearND edge overshoot).
CAP_FFLUX = False                # output cap on the flux-recovery factor: removed
MIN_FRAC = 0.25                  # floor on the implied recovered fraction: the flux-
                                 # recovery factor is capped at F_flux <= 1/MIN_FRAC.
                                 # Real getsf cores bottom out near frac ~0.3-0.5, and
                                 # below that C_M is huge and very uncertain; capping
                                 # under-corrects that tail but cuts scatter (net +5-7pt
                                 # within-2x in the realistic regime).  Set None to disable.
CUMULATIVE_FRAC = True           # per-band recovered flux fraction = cumulative-source
                                 # integral int_0^R_lambda prof_lambda dA / total, evaluated
                                 # inside a self-consistent, beam-enlarged per-band footprint
                                 # R_lambda (see build_recovery_samples).  Set False for the
                                 # legacy rim-subtracted frac at a single common footprint.
#
# Per-band footprint model (physical structural-noise derivation).
# ----------------------------------------------------------------
# R_lambda is the radius at which the source signal-to-noise ratio at the
# footprint edge equals the (band-independent) edge value, with the structural
# noise convolved to that band's beam.  That constant edge S/N cancels between
# the band and surface-density footprints; using the dust-emissivity cancellation
# (source and noise share e_lambda(T)=kappa(lambda)B_lambda(T)) the footprint
# condition reduces to
#   prof_lambda(R_lambda)/prof_lambda(0) = [prof_SD(r)/prof_SD(0)] * c_lambda,
#   c_lambda = sigma(beam_lambda; 1/2R_lambda) / sigma(beam_SD; 1/2r),
# where sigma(beam;k_low) is the RMS of a structural field with 2D power
# spectrum P(k) ~ k^GAMMA, convolved to a Gaussian beam and cut off below k_low
# (= inverse footprint DIAMETER).  R_lambda is solved self-consistently (its own
# footprint sets its noise cut-off), which enforces 2 R_lambda > beam_lambda for
# resolved bands and needs no ad hoc exclusion.  A band is used only where it is
# resolved (2 R_lambda > beam_lambda); an unresolved band is dropped from the SED
# refit (not the whole core).  GAMMA from the Aquila 250um power spectrum (Roy et
# al. 2019, -2.26); c_lambda varies <5% over -3 <= GAMMA <= -2 so -2.5 is used.
GAMMA    = -2.5
BAND_BEAM = {'160': 13.5, '250': 18.2, '350': 24.9, '500': 36.3}
SD_BEAM = 13.5

# getsf R_foot: symmetric, area-equivalent analogue of the model 2r/FWHM (adopted).
RFOOT_MODE = 'symmetric'         # 'symmetric' | 'major' (chat-9 production variant)

# SED fitting (getsf 'thinbody' equivalent)
MASS_SOURCE = 'internal'         # 'internal' (fit here) | 'catalog' (read TOTL^MASS)
SED_WAVES   = (160.0, 250.0, 350.0, 500.0)
SED_ADDERR  = (0.20, 0.10, 0.10, 0.10)
SED_BETA    = 2.0
OPACITY_K   = 5.326076e-24       # M[Msun] = A * OPACITY_K * D_pc^2 (calibrated to
                                 # getsf; encodes kappa_0, lambda_ref, beta, d2g)
SYS_OPACITY_RELERR = 0.20        # relative errors folded into MASS_ERRT (as fitfluxes)
SYS_DUSTGAS_RELERR = 0.20

H_PL, C_L, K_B = 6.62607015e-27, 2.99792458e10, 1.380649e-16


# =====================  MODEL-GRID SED refit (for M_rec)  ====================
_LAM3 = np.array([250., 350., 500.]); _KAP3 = (100.0 / _LAM3) ** 2
_TG = np.linspace(3.0, 60.0, 571)
def _planck(lam, T):
    nu = C_L / (np.asarray(lam)[:, None] * 1e-4)
    return 2 * H_PL * nu**3 / C_L**2 / np.expm1(H_PL * nu / (K_B * T[None, :]))
_G3 = _KAP3[:, None] * _planck(_LAM3, _TG); _GG3 = np.sum(_G3 * _G3, axis=0)

def _mbb_amp3(F):
    A_T = (F @ _G3) / _GG3
    return A_T[np.argmin(np.sum((F[:, None] - A_T[None, :] * _G3) ** 2, axis=0))]

def _mbb_amp_sub(F, mask):
    """MBB amplitude (proportional to mass) fitted to the subset of the three
    continuum bands selected by boolean 'mask' (>=2 bands)."""
    G = _G3[mask]; GG = np.sum(G * G, axis=0); Fs = F[mask]
    A_T = (Fs @ G) / GG
    return A_T[np.argmin(np.sum((Fs[:, None] - A_T[None, :] * G) ** 2, axis=0))]

# structural-noise RMS after convolution to a Gaussian beam (FWHM in arcsec),
# for a 2D power spectrum P(k) ~ k^GAMMA integrated from k_low upward.
_A_BEAM = 4 * np.pi**2 / (8 * np.log(2))          # Gaussian-beam transfer-function coeff
_SIGB = {}
def _sigma_beam(beam, klow):
    """RMS(beam, k_low)/sqrt(2 pi A); only ratios are used so the constant drops out.
    Cached per beam on a log grid of k_low (arcsec^-1)."""
    if beam not in _SIGB:
        from scipy import integrate
        klg = np.logspace(np.log10(1.0 / 2000.0), np.log10(1.0 / 8.0), 90)
        sv = np.array([np.sqrt(max(integrate.quad(
                lambda k: k**(GAMMA + 1) * np.exp(-_A_BEAM * beam**2 * k**2),
                kl, np.inf, limit=200)[0], 0.0)) for kl in klg])
        _SIGB[beam] = (klg, sv)
    klg, sv = _SIGB[beam]
    return float(np.interp(klow, klg, sv))


# =====================  cloud SED fitter (thinbody)  ========================
_LAMs = np.array(SED_WAVES); _ADD = np.array(SED_ADDERR)
_KAPs = (1.0 / _LAMs) ** SED_BETA
_Gall = _KAPs[:, None] * _planck(_LAMs, _TG)          # (nband, nT)

SNR_MIN_PEAK = 2.0     # FXP_BST / FXP_ERR, per band, to admit a flux
SNR_MIN_TOTAL = 2.0    # FXT_BST / FXT_ERR, per band, to admit a flux
# Minimum number of bands that must pass those criteria for a spectral energy
# distribution to be fitted at all.  A modified blackbody has two free
# parameters, so a two-band fit has no degrees of freedom: it passes exactly
# through both points however badly they are measured, chi-squared is zero by
# construction, and no goodness-of-fit test can reject it.  Measured on the
# population-matched injections, every core whose reported mass was wrong by
# more than a factor of five was fitted on two bands, and no core fitted on
# three or four was.  The minimum is therefore three.
MIN_SED_BANDS = 3

def sed_fit(F, E, FP, FPE, S, dist_pc):
    """Vectorised thinbody fit reproducing getsf/fitfluxes.  Per-band arrays
    (n, nband): F,E = FXT_BST/FXT_ERR; FP,FPE = FXP_BST/FXP_ERR; S = SIGNM.
    A band enters the fit iff (fitfluxes, snrmin=1) FXP,FXP_ERR,FXT,FXT_ERR>0,
    FXP/FXP_ERR>SNR_MIN_PEAK, FXT/FXT_ERR>SNR_MIN_TOTAL (>=2 bands);
    SIGNM is NOT used.  Returns dict of arrays
    mass [Msun], temp [K], mass_err, temp_err [1-sigma fit errors], chi2, dof
    (=nbands-2), quality ('ok' iff chi2<dof+1, matching fitfluxes to ~99.7%)."""
    F = np.asarray(F, float); E = np.asarray(E, float)
    FP = np.asarray(FP, float); FPE = np.asarray(FPE, float); S = np.asarray(S, float)
    n = len(F)
    mass = np.full(n, np.nan); temp = np.full(n, np.nan)
    m_err = np.full(n, np.nan); t_err = np.full(n, np.nan)
    chi2 = np.full(n, np.nan); dof = np.full(n, -1, int)
    quality = np.array(['bad'] * n, dtype=object)
    with np.errstate(divide='ignore', invalid='ignore'):
        # A band enters the fit on the strength of its own measurement,
        # not on getsf's per-band detection significance SIGNM.  SIGNM
        # carries the sentinel 9.999e-31 wherever the source was not
        # detected in that band's clean single scales; since every source
        # is detected in at least one band and then measured in all of
        # them, the sentinel records where the detection happened rather
        # than whether the fluxes are usable.  The thresholds below are
        # getsf's own recommended optimal selection criteria, quoted in
        # every catalogue header.  See patch_sed_signm.py.
        keepm = ((FP > 0) & (FPE > 0) & (F > 0) & (E > 0)
                 & (FP / FPE > SNR_MIN_PEAK) & (F / E > SNR_MIN_TOTAL))
    Tg = _TG
    for i in range(n):
        keep = keepm[i]
        nb = int(keep.sum())
        if nb < MIN_SED_BANDS:
            continue
        g = _Gall[keep]                                   # (nk, nT)
        f = F[i][keep]; sig2 = E[i][keep]**2 + (_ADD[keep] * f)**2
        A_T = (f / sig2) @ g / ((1.0 / sig2) @ (g * g))   # (nT,)
        resid = np.sum(((f[:, None] - A_T[None, :] * g) ** 2) / sig2[:, None], axis=0)
        j = np.argmin(resid)
        T = Tg[j]; A = A_T[j]; c2 = resid[j]; d = nb - 2
        M = A * OPACITY_K * dist_pc**2
        mass[i] = M; temp[i] = T; chi2[i] = c2; dof[i] = d
        quality[i] = 'ok' if c2 < d + 1 else 'bad'
        # 1-sigma fit errors from the local chi2 curvature (delta-chi2 = 1).
        # Valid for 2-band fits (dof=0) too: profiling the amplitude out leaves a
        # residual constraint, so chi2(T) still has curvature about the minimum.
        within = (resid - c2) <= 1.0
        Tw = Tg[within]
        if len(Tw):
            t_err[i] = 0.5 * (Tw.max() - Tw.min())
            Mw = A_T[within] * OPACITY_K * dist_pc**2
            m_fit = 0.5 * (Mw.max() - Mw.min())
            m_err[i] = np.hypot.reduce([m_fit, SYS_OPACITY_RELERR * M,
                                        SYS_DUSTGAS_RELERR * M])
    return dict(mass=mass, temp=temp, mass_err=m_err, temp_err=t_err,
                chi2=chi2, dof=dof, quality=quality)


# ==========================  grid catalog loaders  ==========================
P = dict(i=1, j=2, k=3, SD=4, MBE=7, RBEpc=10, M3bs=20, M3bsl=23,
         FXSDbs=32, FXSDbsl=53, FX250=36, FX350=37, FX500=38)

def load_cat67(path):
    nodes = {}
    for l in open(path):
        if not l.strip() or l.lstrip().startswith('#'):
            continue
        r = l.split()
        nodes[(int(r[P['i']]), int(r[P['j']]), int(r[P['k']]))] = dict(
            SD=float(r[P['SD']]), MBE=float(r[P['MBE']]), RBEpc=float(r[P['RBEpc']]),
            M3bs=float(r[P['M3bs']]),
            M3bsl=float(r[P['M3bsl']]), FXSDbs=float(r[P['FXSDbs']]),
            FXSDbsl=float(r[P['FXSDbsl']]),
            F=np.array([float(r[P['FX250']]), float(r[P['FX350']]), float(r[P['FX500']])]))
    return nodes

COLS_REC = ['n','i','j','k','SD_emb','T_BE','rho_BE','M_BE','R_BE_as','R_BE_au',
    'R_BE_pc','Td_avg','Td_emb','T_SED4bs','M_SED4bs','Chi2_4bs','T_SED4bsl',
    'M_SED4bsl','Chi2_4bsl','T_SED3bs','M_SED3bs','Chi2_3bs','T_SED3bsl',
    'M_SED3bsl','Chi2_3bsl','ICSDbs','IC070bs','IC100bs','IC160bs','IC250bs',
    'IC350bs','IC500bs','FWHMSDbs','FWHM070bs','FWHM100bs','FWHM160bs',
    'FWHM250bs','FWHM350bs','FWHM500bs','ICSDbsl','IC070bsl','IC100bsl',
    'IC160bsl','IC250bsl','IC350bsl','IC500bsl','FWHMSDbsl','FWHM070bsl',
    'FWHM100bsl','FWHM160bsl','FWHM250bsl','FWHM350bsl','FWHM500bsl']
RECX = ['M_SED3bs_rec','M_SED3bsl_rec','frac_rec','FWHM_rec','conc_footfwhm',
        'conc_peakmean','conc_slope']

def load_catrec(path):
    rows = [l.split() for l in open(path)
            if l.strip() and not l.lstrip().startswith('#')]
    out = {c: np.array([r[i] for r in rows], float) for i, c in enumerate(COLS_REC)}
    nc = len(rows[0])
    for off, name in enumerate(RECX):
        out[name] = np.array([r[nc - len(RECX) + off] for r in rows], float)
    return out


# ============  per-band M_rec samples from recovery tables  ==================

# Outcome counts of the per-band footprint solver _band_R(), accumulated over a
# call to build_recovery_samples().  Reset it before the call if you want counts
# for that call alone.
BANDR_DIAG = dict(converged=0, table_too_short=0, root_below_r=0, no_profile=0)

def build_recovery_samples(cat_path=CAT67, rectab_path=RECTAB):
    vlog(0, "Building recovery samples from %s + %s ..." % (cat_path, rectab_path))
    t0 = time.time()
    cat = load_cat67(cat_path)
    vlog(1, "  loaded %d model-grid nodes from %s" % (len(cat), cat_path))
    band_fr = defaultdict(lambda: defaultdict(lambda: ([], [])))   # rim-subtracted frac
    band_prof = defaultdict(lambda: defaultdict(list))             # (r, I_rim) profile
    band_peak = defaultdict(lambda: defaultdict(float))            # central intensity
    sd_rows = defaultdict(list)
    n_rectab_lines = 0
    for l in open(rectab_path):
        if l.startswith('#') or not l.strip():
            continue
        n_rectab_lines += 1
        f = l.split(); ijk = (int(f[0]), int(f[1]), int(f[2])); b = f[3]
        try:
            r, peak, irim, frac, cf, cpm, ab = (float(f[5]), float(f[7]), float(f[8]),
                                                float(f[10]), float(f[12]),
                                                float(f[13]), int(f[15]))
        except ValueError:
            continue
        if b in ('160', '250', '350', '500'):
            if np.isfinite(frac):
                band_fr[ijk][b][0].append(r); band_fr[ijk][b][1].append(frac)
        if b in ('SD', '160', '250', '350', '500') and np.isfinite(irim):
            band_prof[ijk][b].append((r, irim))    # incl. SD, for the footprint solver
            band_peak[ijk][b] = max(band_peak[ijk][b], peak)
        if b == 'SD' and np.isfinite(frac) and ab == 0:
            sd_rows[ijk].append((r, frac, cf, cpm))

    vlog(1, "  parsed %d recovery-table lines -> %d nodes with usable SD-band rows"
         % (n_rectab_lines, len(sd_rows)))

    def frac_at(node, b, r):                       # rim-subtracted frac (legacy)
        rs, fr = band_fr[node][b]
        if len(rs) < 2:
            return np.nan
        rs = np.asarray(rs); fr = np.asarray(fr); o = np.argsort(rs)
        return float(np.clip(np.interp(r, rs[o], fr[o]), FRAC_FLOOR, 1.0))

    _cum = {}
    def cumfrac_at(node, b, r):                     # cumulative-source frac (default)
        key = (node, b)
        if key not in _cum:
            pr = band_prof[node].get(b, [])
            if len(pr) < 2:
                _cum[key] = None
            else:
                a = np.array(sorted(pr))
                rr = np.linspace(0.0, a[-1, 0], 400)
                pp = np.clip(np.interp(rr, np.concatenate([[0.0], a[:, 0]]),
                                       np.concatenate([[band_peak[node][b]], a[:, 1]])),
                             0, None)
                w = 2 * np.pi * rr
                cum = np.concatenate([[0.0], np.cumsum(
                    0.5 * (pp[1:] * w[1:] + pp[:-1] * w[:-1]) * np.diff(rr))])
                _cum[key] = (rr, cum / cum[-1]) if cum[-1] > 0 else None
        c = _cum[key]
        if c is None:
            return np.nan
        return float(np.clip(np.interp(r, c[0], c[1]), FRAC_FLOOR, 1.0))

    band_frac = cumfrac_at if CUMULATIVE_FRAC else frac_at

    # ---- self-consistent, beam-enlarged per-band footprints R_lambda ----
    _pn = {}
    def _prof_norm(node, b):                        # (r_grid, prof/prof(0)), monotone decr.
        key = (node, b)
        if key not in _pn:
            pr = band_prof[node].get(b, [])
            if len(pr) < 3:
                _pn[key] = None
            else:
                a = np.array(sorted(pr))
                r = np.concatenate([[0.0], a[:, 0]])
                p = np.concatenate([[band_peak[node][b]], a[:, 1]])
                _pn[key] = (r, p / p[0]) if p[0] > 0 else None
        return _pn[key]

    def _R_of_level(node, b, level):                # radius where prof/prof(0) = level
        pa = _prof_norm(node, b)
        if pa is None:
            return np.nan
        r, pn = pa; below = np.where(pn <= level)[0]
        if len(below) == 0:
            return r[-1]
        i = below[0]
        return 0.0 if i == 0 else float(np.interp(level, [pn[i], pn[i - 1]],
                                                  [r[i], r[i - 1]]))

    def _level_SD(node, r):                         # prof_SD(r)/prof_SD(0)
        pa = _prof_norm(node, 'SD')
        if pa is None:
            return np.nan
        rr, pn = pa; return float(np.interp(r, rr, pn))

    def _band_R(node, b, r):
        """Per-band footprint radius R_lambda, from the condition of equal source
        signal-to-noise ratio at the footprint edge (paper Sect. 3.2):

            p_b(R) = p_SD(r) * c_b(R),
            c_b(R) = sigma(O_b; k_low = 1/(2R)) / sigma(O_SD; k_low = 1/(2r)),

        with p_x(s) the band-x intensity profile normalised to its own centre.
        p_b(R) decreases with R while c_b(R) increases with R, so
        g(R) = p_b(R) - p_SD(r) c_b(R) is strictly decreasing and its root is
        unique.  Solved by bracketed root finding (Brent) on [r, R_tab], where
        R_tab is the last tabulated radius of band b.

        This replaces an unguarded fixed-point iteration whose map has negative
        derivative of magnitude greater than one over much of the grid.  That
        iteration oscillated in 72% of calls, failed its own convergence test in
        23%, and silently returned the last iterate or clamped to the end of the
        profile table, making R_lambda(r) non-monotone in r in 80-91% of nodes.
        Where no root exists inside the tabulated profile this version returns
        NaN, so the sample is dropped rather than assigned a fabricated radius.
        Outcome counts accumulate in the module-level dict BANDR_DIAG.
        """
        if not CUMULATIVE_FRAC:
            return r
        levSD = _level_SD(node, r)
        if not (levSD > 0):
            BANDR_DIAG['no_profile'] += 1
            return np.nan
        pa = _prof_norm(node, b)
        if pa is None:
            BANDR_DIAG['no_profile'] += 1
            return np.nan
        rtab, ptab = pa
        Rhi = float(rtab[-1])
        sSD = _sigma_beam(SD_BEAM, 1.0 / (2 * r))

        def _g(R):
            pb = float(np.interp(R, rtab, ptab))
            cb = _sigma_beam(BAND_BEAM[b], 1.0 / (2 * R)) / sSD
            return pb - levSD * cb

        if not (Rhi > r):
            BANDR_DIAG['table_too_short'] += 1
            return np.nan
        if _g(r) < 0.0:                     # root would lie below r: the band
            BANDR_DIAG['root_below_r'] += 1  # footprint cannot be smaller than
            return r                         # the Sigma-band one
        if _g(Rhi) > 0.0:                   # profile table does not reach the root
            BANDR_DIAG['table_too_short'] += 1
            return np.nan
        from scipy.optimize import brentq
        R = float(brentq(_g, r, Rhi, xtol=1e-3, rtol=1e-8, maxiter=200))
        BANDR_DIAG['converged'] += 1
        return R

    node, sig, conc, cfoot, mrec, mbe, cm = [], [], [], [], [], [], []
    rbe_pc = []                                   # physical truncation radius (parsecs);
                                                   # the 4th ("4D") corrector axis
    frac = []                                     # SD recovered fraction (for validation)
    f_flux, f_bg = [], []                         # decomposed mass factors
    f160, f250, f350, f500 = [], [], [], []       # per-band recovery, for the refit
    n_nodes_total = len(sd_rows)
    n_nodes_done = 0
    n_drop_notincat = n_drop_badF = n_drop_norbe = 0
    n_drop_unresolved = n_drop_nonfinite = n_drop_unreliable = n_kept = 0
    for ijk, rows in sd_rows.items():
        n_nodes_done += 1
        if VERBOSE >= 2:
            vlog(2, "  node %d/%d  ijk=%s  (%d truncation-radius rows)"
                 % (n_nodes_done, n_nodes_total, ijk, len(rows)))
        elif VERBOSE >= 1 and n_nodes_done % 50 == 0:
            vlog(1, "  ...%d/%d nodes processed" % (n_nodes_done, n_nodes_total))
        if ijk not in cat:
            n_drop_notincat += 1; continue
        c = cat[ijk]
        if not (c['M3bs'] > 0 and np.all(np.isfinite(c['F'])) and np.all(c['F'] > 0)):
            n_drop_badF += 1; continue
        if not (np.isfinite(c['RBEpc']) and c['RBEpc'] > 0):
            n_drop_norbe += 1; continue
        fbg = (c['FXSDbsl'] / c['FXSDbs']
               if c['FXSDbs'] > 0 and c['FXSDbsl'] > 0 else np.nan)
        for (r, frSD, cf, cpm) in rows:
            bands = ('250', '350', '500')
            R = {b: _band_R(ijk, b, r) for b in bands}
            # a band enters the SED refit only where resolved: 2 R_lambda > beam_lambda
            resolved = np.array([np.isfinite(R[b]) and 2 * R[b] > BAND_BEAM[b]
                                 for b in bands])
            if resolved.sum() < 2:                 # need >= 2 resolved bands for the MBB fit
                n_drop_unresolved += 1; continue
            fb = np.array([band_frac(ijk, b, R[b]) if resolved[i] else 1.0
                           for i, b in enumerate(bands)])
            if not np.all(np.isfinite(fb[resolved])):
                n_drop_nonfinite += 1; continue
            # A band recovering less than MIN_FRAC of its flux is not a reliable
            # constraint on the SED refit: at low recovered fraction, differential
            # per-band footprint noise can distort the reduced-flux SED shape into
            # one no physical modified blackbody can reproduce (Sect.~sec:pathology),
            # driving the fit to an unphysical temperature and an amplitude wrong by
            # up to ~1000x. Require >=3 bands clearing MIN_FRAC; otherwise the sample
            # is too unreliable to train on and is discarded (not merely down-weighted).
            reliable = resolved & (fb >= MIN_FRAC)
            if reliable.sum() < 3:
                n_drop_unreliable += 1; continue
            n_kept += 1
            mask = reliable
            A0 = _mbb_amp_sub(c['F'], mask)
            Mrec = c['M3bs'] * _mbb_amp_sub(c['F'] * fb, mask) / A0
            if not Mrec > 0:
                continue
            f16 = band_frac(ijk, '160', r)        # 160 beam = SD beam; may be nan
            node.append(ijk); sig.append(c['SD']); conc.append(cpm); cfoot.append(cf)
            mrec.append(Mrec); mbe.append(c['MBE']); cm.append(c['MBE'] / Mrec)
            rbe_pc.append(c['RBEpc'])
            frac.append(frSD)
            f_flux.append(c['M3bs'] / Mrec)       # flux recovery (M_SED3bs / M_reported)
            f_bg.append(fbg)                      # bg subtraction (FXSDbsl / FXSDbs)
            f160.append(f16)
            f250.append(fb[0] if mask[0] else np.nan)
            f350.append(fb[1] if mask[1] else np.nan)
            f500.append(fb[2] if mask[2] else np.nan)

    vlog(1, "  filtering summary: kept=%d  dropped: not-in-catalog=%d  bad-flux=%d  "
         "no-R_BE=%d  unresolved(<2 bands)=%d  non-finite=%d  unreliable(<3 "
         "bands>=MIN_FRAC)=%d" % (n_kept, n_drop_notincat, n_drop_badF, n_drop_norbe,
                                   n_drop_unresolved, n_drop_nonfinite, n_drop_unreliable))
    vlog(0, "Recovery samples built: %d samples from %d/%d nodes, in %.1fs"
         % (len(cm), len(set(node)), n_nodes_total, time.time() - t0))

    return dict(node=np.array(node, dtype=object), sig=np.array(sig),
                conc=np.array(conc), cfoot=np.array(cfoot), mrec=np.array(mrec),
                mbe=np.array(mbe), cm=np.array(cm), rbe_pc=np.array(rbe_pc),
                frac=np.array(frac),
                f_flux=np.array(f_flux), f_bg=np.array(f_bg),
                f160=np.array(f160), f250=np.array(f250),
                f350=np.array(f350), f500=np.array(f500))

# ----------------------- local linear corrector estimator -------------------
GRID_DISTANCE_PC = 260.0   # distance at which every model stamp of the grid was
                           # rendered.  Used ONLY to express the grid's physical
                           # truncation radius R_BE as the angle it subtends, so
                           # that the fourth corrector axis is angular on both
                           # the model side and the real-source side.

LL_BOUND = True          # bound each prediction by the range of the training
                         # values used for it; see the note in LocalLinear
LL_NEIGHBOURS = 320      # training samples used per prediction; the leave-one-out
                         # residual is flat over 80-640, so this is not critical

class LocalLinear:
    """Tricube-weighted local linear regression, as a drop-in replacement for
    scipy.interpolate.LinearNDInterpolator.

    Constructed from (points, values) and called with an (n, ndim) array of
    query points, returning an (n,) array of predicted values.  Each axis is
    divided by its own standard deviation before distances are computed, so the
    metric does not depend on the units of the axes.  Unlike a Delaunay
    interpolator it extrapolates outside the convex hull of the training points
    rather than returning NaN.
    """

    def __init__(self, points, values, k=None):
        X = np.asarray(points, float)
        v = np.asarray(values, float)
        ok = np.all(np.isfinite(X), 1) & np.isfinite(v)
        self.X0 = X[ok]
        self.v = v[ok]
        s = self.X0.std(0)
        s[~(s > 0)] = 1.0
        self.scale = s
        self.Xs = self.X0 / s
        self.ndim = self.Xs.shape[1]
        self.tree = cKDTree(self.Xs)
        self.k = int(min(k or LL_NEIGHBOURS, len(self.v)))
        self.n_bounded = 0
        self._flag = np.zeros(0, bool)
        self._support = np.zeros(0)

    def __call__(self, q):
        Q = np.atleast_2d(np.asarray(q, float)) / self.scale
        out = np.full(len(Q), np.nan)
        self._flag = np.zeros(len(Q), bool)
        self._support = np.full(len(Q), np.nan)
        if self.k < self.ndim + 2 or len(self.v) < self.ndim + 2:
            return out
        dd, ii = self.tree.query(Q, k=self.k)
        dd = np.atleast_2d(dd)
        ii = np.atleast_2d(ii)
        for m in range(len(Q)):
            d, idx = dd[m], ii[m]
            g = np.isfinite(d)
            d, idx = d[g], idx[g]
            if len(idx) < self.ndim + 2:
                continue
            h = d[-1] if d[-1] > 0 else 1.0
            w = np.clip(1.0 - (d / h) ** 3, 0.0, None) ** 3
            w = np.maximum(w, 1e-6)
            A = np.column_stack([np.ones(len(idx)), self.Xs[idx] - Q[m]])
            sw = np.sqrt(w)
            try:
                b, *_ = np.linalg.lstsq(A * sw[:, None], self.v[idx] * sw,
                                        rcond=None)
            except Exception:
                continue
            v = b[0]
            lo, hi = self.v[idx].min(), self.v[idx].max()
            if LL_BOUND and not (lo <= v <= hi):
                # The plane has been evaluated outside the range spanned by the
                # training samples it was fitted to, i.e. this is extrapolation
                # of a first-order fit, not interpolation.  Bound it and record
                # the fact, rather than returning an unbounded value.
                self.n_bounded += 1
                v = min(max(v, lo), hi)
                self._flag[m] = True
            self._support[m] = d[-1]
            out[m] = v
        return out

    def last_call_bounded(self):
        """Boolean array, one entry per query of the most recent call, True where
        the prediction had to be bounded by the local training range."""
        return self._flag.copy()

    def last_call_support(self):
        """Distance to the furthest training sample used, per query of the most
        recent call, in units of the standard deviation of each axis.  Large
        values mark queries far from any training data."""
        return self._support.copy()


# CORRECTOR AXES.  The corrector predicts log10(C_M) from three distance- AND
# cloud-invariant geometric observables measured on the 13.5" H2 band: Sigma
# (PEAK^BGF), conc [zeta] (PEAK^SRC/footprint-mean) and R_foot [phi]
# (footprint/FWHM [H]).
#
# Dust temperature was evaluated as extra axes (peak T from FXP, total T from FXT,
# rim/ambient T) and REJECTED: it improves the model self-test (within-2x 58->83%)
# but does not transfer to real clouds.  The model grid's SED temperatures run
# cooler and narrower than observed (T_SED3bs ~9 K vs Aquila T_INIT ~13 K), because
# the absolute dust temperature is set by the local ISRF (G0) -- a cloud property.
# Any temperature axis therefore makes the corrector cloud-specific and breaks the
# universality that is the whole point of the geometric method.  Kept geometry-only.
#
# LATER FINDING (see paper, not implemented in cm_interp below): these three axes
# leave a real degeneracy in C_M -- near-identical (Sigma, conc, R_foot) can still
# imply true C_M differing by 2x+ -- traced to missing absolute physical scale, NOT
# to missing temperature information (that was the rejected axis above). A fourth
# axis, the model's physical truncation radius R_BE [parsecs], matched to a real
# source's footprint radius at its known distance, tested well in three independent
# checks. Unlike the temperature axis above, this does NOT make the corrector
# cloud-specific, but it does require a known distance, unlike the pure
# three-axis design.
def cm_interp(S):
    """Three-observable corrector: predicts log10(C_M) from
    (Sigma, zeta, phi).  Uses a local linear regression (LocalLinear above)
    rather than a Delaunay interpolator; see patch_estimator.py for the
    leave-one-out comparison that motivated the change."""
    return LocalLinear(np.column_stack(
        [np.log10(S['sig']), np.log10(S['conc']), np.log10(S['cfoot'])]),
        np.log10(S['cm']))

def mass_idx_interp(S):
    """Predicts M_BE from (M_rec, Sigma, conc, R_foot) -- an experimental,
    mass-INDEXED 4-axis approach kept for diagnostic comparison only. This is
    the approach the paper's Appendix B rejects (reported mass scales as
    distance^2 and is itself the biased quantity under correction, so using
    it as a matching axis reintroduces a distance-confounded failure mode).
    NOT the paper's adopted fourth observable -- see cm_size_interp() below,
    which uses the model's physical size instead and is unrelated to this
    function despite the superficially similar name.
    """
    return LinearNDInterpolator(np.column_stack(
        [np.log10(S['mrec']), np.log10(S['sig']), np.log10(S['conc']),
         np.log10(S['cfoot'])]), np.log10(S['mbe']))

def cm_size_interp(S):
    """Four-observable corrector: predicts log10(C_M) from
    (Sigma, zeta, phi, theta_BE), where theta_BE is the ANGLE subtended by the
    model's truncation radius, in arcsec, matched on the real-source side by the
    footprint radius sqrt(A_F*B_F)/2 measured directly in arcsec.

    The axis is angular, not physical.  The correction factor is a ratio of two
    masses, C_M = M_BE / M_reported, and each is proportional to the square of
    the distance, so the distance cancels exactly.  Surface brightness does not
    depend on distance either, and both the beam and the structural noise that
    fixes the footprint are properties of the map in its angular frame.  Every
    quantity C_M depends on is therefore an angle, and C_M cannot depend on how
    far away the source lies.

    Because the whole grid was rendered at one distance, the model-side axis is
    the same variable whether expressed in parsecs or arcsec, and no grid-based
    test can tell the two apart.  The difference appears only when a real source
    is looked up.  Converting its angular footprint to parsecs at its own
    distance, as earlier versions did, introduced a spurious distance trend:
    across seven clouds spanning 139 to 2000 pc the median correction factor
    ranged over a ratio of 1.70, correlating with log distance at 0.958.  With
    the angular axis the same seven clouds range over 1.12, with a correlation
    of 0.464.
    """
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    theta_be = S['rbe_pc'][ok] * 206264.806 / GRID_DISTANCE_PC     # arcsec
    return LocalLinear(np.column_stack(
        [np.log10(S['sig'][ok]), np.log10(S['conc'][ok]), np.log10(S['cfoot'][ok]),
         np.log10(theta_be)]), np.log10(S['cm'][ok]))


def factor_interp(S, key):
    """Interpolator for a decomposed mass factor on (Sigma, conc, R_foot)."""
    v = S[key]; ok = np.isfinite(v) & (v > 0)
    X = np.column_stack([np.log10(S['sig'][ok]), np.log10(S['conc'][ok]),
                         np.log10(S['cfoot'][ok])])
    return LinearNDInterpolator(X, np.log10(v[ok]))

def recovery_interp(S):
    """Per-band flux-recovery fraction f_b(Sigma, conc, R_foot), for the refit.
    Keyed by wavelength; a dict of LinearNDInterpolators on log axes."""
    out = {}
    for w in (160, 250, 350, 500):
        f = S['f%d' % w]; ok = np.isfinite(f) & (f > 0)
        X = np.column_stack([np.log10(S['sig'][ok]), np.log10(S['conc'][ok]),
                             np.log10(S['cfoot'][ok])])
        out[w] = LinearNDInterpolator(X, np.log10(f[ok]))
    return out


# ===============================  getsf reader  =============================
def _num(s):
    try:
        float(s); return True
    except ValueError:
        return False

# ===========================================================================
# Column location by header position
#
# getsf catalogues are fixed-width tables whose column COUNT varies with the
# waveband set, with which auxiliary blocks have been appended, and with the
# getsf version, so addressing columns by number is correct only for the file a
# script was written against.  The code below reproduces the scheme used inside
# getsf itself: find the line that names the columns, record the character span
# of each name, and read each value from the data token lying under that span.
# Names are normalised by replacing '^' with '_', so TOTL^MASS and TOTL_MASS are
# a single lookup.
# ===========================================================================
_COMMENT = ('#', '!')


def _normalise(name):
    return name.replace('^', '_')


def _tokens_with_spans(line):
    """[(token, start, end)] for every whitespace-delimited token."""
    return [(m.group(0), m.start(), m.end())
            for m in re.finditer(r'\S+', line)]


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class GetsfTable(object):

    def __init__(self, path=None, _empty=False):
        self.path = path
        self.names = []
        self.header = []
        self._cols = {}          # normalised name -> list of raw strings
        self.unmatched_names = []
        self.nrows = 0
        if not _empty:
            self._parse(path)

    # ------------------------------------------------------------------ parse
    def _parse(self, path):
        lines = open(path).read().split('\n')
        is_com = lambda l: l.lstrip()[:1] in _COMMENT if l.strip() else False
        self.header = [l for l in lines if is_com(l)]

        data_idx = [i for i, l in enumerate(lines)
                    if l.strip() and not is_com(l)]
        if not data_idx:
            raise ValueError('no data rows found in %s' % path)
        first_data = data_idx[0]

        # The column-name line is the last comment line above the data that
        # contains several plausible names.  Scanning upward rather than
        # matching fixed content keeps this independent of the waveband set and
        # of any appended blocks.
        hdr_line = None
        for i in range(first_data - 1, -1, -1):
            if not is_com(lines[i]):
                continue
            stripped = lines[i].lstrip('#! ').rstrip()
            toks = stripped.split()
            if len(toks) < 4:
                continue
            good = sum(1 for t in toks
                       if re.match(r'^[A-Za-z][A-Za-z0-9_^]*$', t))
            if good >= max(4, int(0.8 * len(toks))):
                hdr_line = lines[i]
                break
        if hdr_line is None:
            raise ValueError('could not find a column-name line in %s' % path)

        # strip only the leading comment marker, preserving character positions
        lead = len(hdr_line) - len(hdr_line.lstrip('#! '))
        hdr_spans = [(_normalise(t), s, e)
                     for (t, s, e) in _tokens_with_spans(hdr_line)
                     if s >= lead - 1]
        hdr_spans = [(n, s, e) for (n, s, e) in hdr_spans
                     if re.match(r'^[A-Za-z][A-Za-z0-9_]*$', n)]
        self.names = [n for (n, _, _) in hdr_spans]

        cols = dict((n, []) for n in self.names)
        seen = set()
        rows = [lines[i] for i in data_idx]
        for row in rows:
            toks = _tokens_with_spans(row)
            used = [False] * len(toks)
            vals = {}
            # assign by closest right edge, then by closest centre
            for (n, hs, he) in hdr_spans:
                best, bestd = None, None
                for k, (t, s, e) in enumerate(toks):
                    if used[k]:
                        continue
                    d = abs(e - he)
                    if bestd is None or d < bestd:
                        best, bestd = k, d
                if best is None:
                    continue
                hc, tc = 0.5 * (hs + he), 0.5 * (toks[best][1] + toks[best][2])
                for k, (t, s, e) in enumerate(toks):
                    if used[k]:
                        continue
                    if abs(0.5 * (s + e) - hc) < abs(tc - hc):
                        best = k
                        tc = 0.5 * (s + e)
                used[best] = True
                vals[n] = toks[best][0]
                seen.add(n)
            for n in self.names:
                cols[n].append(vals.get(n, ''))
        self._cols = cols
        self.nrows = len(rows)
        self.unmatched_names = [n for n in self.names if n not in seen]

    # ----------------------------------------------------------------- access
    def has(self, name):
        return _normalise(name) in self._cols

    def raw(self, name):
        return self._cols[_normalise(name)]

    def col(self, name):
        return np.array([float(v) if _is_number(v) else np.nan
                         for v in self._cols[_normalise(name)]])

    def find(self, pattern):
        """Names matching a regular expression, e.g. r'^AFWHM\\d\\d$'."""
        rx = re.compile(pattern)
        return [n for n in self.names if rx.match(n)]

    # ------------------------------------------------------------------ merge
    @classmethod
    def merge(cls, tables):
        """Join tables row by row.  Duplicated names keep the first occurrence."""
        n = min(t.nrows for t in tables)
        out = cls(_empty=True)
        out.path = [t.path for t in tables]
        for t in tables:
            out.header += t.header
            for nm in t.names:
                if nm in out._cols:
                    continue
                out.names.append(nm)
                out._cols[nm] = t._cols[nm][:n]
        out.nrows = n
        return out


def summarise(path):
    t = GetsfTable(path)
    print('%s: %d columns, %d rows' % (path, len(t.names), t.nrows))
    if t.unmatched_names:
        print('  WARNING: no data token found for %s' % t.unmatched_names)
    return t


if __name__ == '__main__':
    import sys
    for p in sys.argv[1:]:
        t = summarise(p)
        print('  first names: %s' % ' '.join(t.names[:8]))
        for probe in ('TOTL_MASS', 'DUST_TEMP', 'QUALITY', 'FOOA03', 'AFWHM03'):
            if t.has(probe):
                v = t.col(probe)
                print('  %-10s median %.4g over %d finite values'
                      % (probe, np.nanmedian(v), int(np.sum(np.isfinite(v)))))


def _read_getsf_file(path):
    """Return (col_names, data_rows, header_text) for one getsf/thin catalogue.

    Columns are located by the character position of their header names, not by
    column number, so this works unchanged on the raw .cat, the .add.cat, and
    any concatenated catalogue regardless of how many columns it carries.
    """
    t = GetsfTable(path)
    if t.unmatched_names:
        vlog(1, "    %s: no data token found for %s"
             % (path, t.unmatched_names))
    data = [[t.raw(n)[i] for n in t.names] for i in range(t.nrows)]
    return t.names, data, t.header


def _wave_config(header):
    waves = []
    for l in header:
        s = l.lstrip('#!').strip()
        m = re.match(r'^(\d{2,3})\s+([\d.]+)\s+-?\d+\s+-?\d+\s+[yn]\s+\w{3}\b', s)
        if m:
            w = int(m.group(1))
            if w in waves:            # a second config block (merged headers) -> stop
                break
            waves.append(w)
    return waves

def load_getsf(path):
    """path: a single concatenated catalogue, or a list [cat, add] (+optional thin).
    Returns per-source arrays for the correction; SED mass fitted internally
    unless MASS_SOURCE='catalog'."""
    paths = [path] if isinstance(path, str) else list(path)
    vlog(1, "  load_getsf: reading %d file(s): %s" % (len(paths), paths))
    names, data, header = [], None, []
    for p in paths:
        nm, dt, hd = _read_getsf_file(p)
        vlog(2, "    %s: %d columns, %d rows" % (p, len(nm), len(dt)))
        header += hd
        if names == []:
            names, data = nm, dt
        else:
            names += nm
            data = [a + b for a, b in zip(data, dt)] if data else dt
    idx = {n: i for i, n in enumerate(names)}
    vlog(1, "  load_getsf: %d sources, %d total columns after merge" % (len(data), len(names)))

    def col(nm):
        i = idx[nm]
        return np.array([float(r[i]) if _num(r[i]) else np.nan for r in data])

    # surface-density 13.5" band: the band measured in H2 column-density units
    # (PEAK_BGF in the 1e21-5e22 range, unlike the continuum bands which are in
    # flux units) with the finest resolution.  Identifying by units is robust; the
    # 13.5" beam is only a sanity bound, since the smallest MEASURED source FWHM can
    # sit slightly above the beam (e.g. 14.6" here) and a strict 13.5" test misses it.
    bands = sorted({m.group(1) for n in idx
                    for m in [re.match(r'AFWHM(\d\d)', n)] if m})
    cand = []
    for c in bands:
        fw = np.sqrt(col('AFWHM'+c) * col('BFWHM'+c)); fw = fw[np.isfinite(fw) & (fw > 0)]
        pb = col('PEAK_BGF'+c); pb = pb[np.isfinite(pb) & (pb > 0)]
        ok = len(fw) and len(pb) and 1e21 < np.median(pb) < 2e23 and fw.min() < 16.0
        vlog(2, "    band %s: median PEAK_BGF=%.3e  median AFWHM~%.2f\"  candidate=%s"
             % (c, np.median(pb) if len(pb) else np.nan,
                np.median(fw) if len(fw) else np.nan, ok))
        if ok:
            cand.append((np.median(fw), c))
    sd = min(cand)[1] if cand else None
    vlog(1, "  load_getsf: SD (surface-density) band identified as code '%s'" % sd)
    if sd is None:
        raise ValueError('surface-density (13.5" H2-column) band not found in ' + str(paths))

    area = np.pi * (col('FOOA'+sd)/2) * (col('FOOB'+sd)/2)
    with np.errstate(divide='ignore', invalid='ignore'):
        conc = col('PEAK_SRC'+sd) / (col('FXT_BST'+sd) / area)
        if RFOOT_MODE == 'symmetric':
            rfoot = (np.sqrt(col('FOOA'+sd) * col('FOOB'+sd))
                     / np.sqrt(col('AFWHM'+sd) * col('BFWHM'+sd)))
        else:
            rfoot = col('FOOA'+sd) / np.sqrt(col('AFWHM'+sd) * col('BFWHM'+sd))
    cat_quality = np.array([r[idx['QUALITY']] for r in data]) if 'QUALITY' in idx \
        else None

    # position (for output catalogues), if present
    def maybe(nm):
        return col(nm) if nm in idx else np.full(len(data), np.nan)
    NO = np.array([int(float(r[0])) if _num(r[0]) else -1 for r in data])
    ra = maybe('WCS_ACOOR'); dec = maybe('WCS_DCOOR')

    # reported mass + SED fit (temperature, errors, chi2, dof, computed QUALITY)
    temp = mass_err = temp_err = chi2 = dof = None
    sed_flux = None
    if MASS_SOURCE == 'catalog' and 'TOTL_MASS' in idx:
        mass = col('TOTL_MASS'); temp = maybe('DUST_TEMP')
        mass_err = maybe('MASS_ERRO'); temp_err = maybe('TEMP_ERRO')
        quality = cat_quality if cat_quality is not None else np.array(['ok']*len(data))
    else:
        waves = _wave_config(header)
        sedidx = [k + 1 for k, w in enumerate(waves) if float(w) in SED_WAVES]
        if len(sedidx) != len(SED_WAVES):
            raise ValueError('could not map SED bands %s from waves %s'
                             % (list(SED_WAVES), waves))
        F = np.array([col('FXT_BST%02d' % b) for b in sedidx]).T
        E = np.array([col('FXT_ERR%02d' % b) for b in sedidx]).T
        FP = np.array([col('FXP_BST%02d' % b) for b in sedidx]).T
        FPE = np.array([col('FXP_ERR%02d' % b) for b in sedidx]).T
        Sg = np.array([col('SIGNM%02d' % b) for b in sedidx]).T
        if 'DISTANCE' in idx and np.isfinite(np.nanmedian(col('DISTANCE'))):
            dpc = np.nanmedian(col('DISTANCE'))
        else:
            dpc = load_getsf._dist
        fit = sed_fit(F, E, FP, FPE, Sg, dpc)
        mass = fit['mass']; temp = fit['temp']
        mass_err = fit['mass_err']; temp_err = fit['temp_err']
        chi2 = fit['chi2']; dof = fit['dof']; quality = fit['quality']
        sed_flux = dict(F=F, E=E, FP=FP, FPE=FPE, S=Sg, dpc=dpc)
    vlog(1, "  load_getsf: %d sources parsed, %d with QUALITY=='ok'"
         % (len(quality), int(np.sum(np.asarray(quality) == 'ok'))))
    return dict(mass=mass, Nbg=col('PEAK_BGF'+sd), peak=col('PEAK_SRC'+sd),
                conc=conc, rfoot=rfoot, quality=quality, temp=temp,
                mass_err=mass_err, temp_err=temp_err, chi2=chi2, dof=dof,
                NO=NO, ra=ra, dec=dec, sed=sed_flux,
                foot_radius_as=np.sqrt(area/np.pi))   # footprint radius (arcsec):
                                                       # sqrt(FOOA*FOOB)/2, the
                                                       # paper's 4th-observable
                                                       # real-source proxy, needs
                                                       # x distance/206265 for pc
load_getsf._dist = None          # fallback distance if DISTANCE column absent

def source_mask(g):
    base = ((g['mass'] > 0) & (g['Nbg'] > 0) & (g['peak'] > 0)
            & np.isfinite(g['conc']) & (g['conc'] > 0)
            & np.isfinite(g['rfoot']) & (g['rfoot'] > 0))
    return base & (g['quality'] == 'ok')       # SED-fit QUALITY (chi2 < dof+1)


# ===================  single-cloud correction + output  =====================
_GRIDS = None
class _HullTest(object):
    """Is a query point inside the convex hull of the training samples?

    Built from the (n, ndim) array of training coordinates.  Calling it with an
    (m, ndim) array of queries returns a boolean array, True where the query
    lies inside the hull.

    The hull is represented by its facet inequalities, from
    scipy.spatial.ConvexHull, so a query is inside when every facet gives
    A.x + b <= 0.  That is far cheaper to build than the Delaunay
    triangulation the pipeline used before version 6, and it answers exactly
    the same question: whether the prediction interpolates between training
    models or extends beyond them.  Axes are standardised first so the
    tolerance means the same on each.

    If the hull cannot be built (degenerate input, or scipy missing), every
    query is reported as inside and a warning is logged, so that the
    diagnostic degrades to uninformative rather than to wrong.
    """

    TOL = 1e-12

    def __init__(self, points):
        self.eq = None
        self.mu = None
        self.sd = None
        if points is None or len(points) == 0:
            return
        X = np.asarray(points, float)
        X = X[np.all(np.isfinite(X), 1)]
        if len(X) < X.shape[1] + 1:
            return
        self.mu = X.mean(0)
        sd = X.std(0)
        sd[~(sd > 0)] = 1.0
        self.sd = sd
        try:
            from scipy.spatial import ConvexHull as _CH
            self.eq = _CH((X - self.mu) / sd).equations
        except Exception as e:
            vlog(0, "  warning: convex-hull test unavailable (%s); "
                    "IN_HULL will be reported as 1 for every source" % e)

    def __call__(self, q):
        Q = np.atleast_2d(np.asarray(q, float))
        if self.eq is None:
            return np.ones(len(Q), bool)
        Qs = (Q - self.mu) / self.sd
        d = Qs @ self.eq[:, :-1].T + self.eq[:, -1]
        return np.all(d <= self.TOL, axis=1)


def build_grids():
    """Build the (cloud-independent) correction grids once and cache them."""
    global _GRIDS
    if _GRIDS is None:
        vlog(0, "Building correctors (3-axis, 4-axis, and diagnostic grids)...")
        S = build_recovery_samples()
        with _Timer(2, 'cm_interp (3-axis)'):
            cm = cm_interp(S)
        with _Timer(2, 'cm_size_interp (4-axis)'):
            cm4 = cm_size_interp(S)
        with _Timer(2, 'mass_idx_interp (diagnostic)'):
            midx = mass_idx_interp(S)
        with _Timer(2, 'factor_interp (F_flux, F_bg)'):
            ff = factor_interp(S, 'f_flux'); fbg = factor_interp(S, 'f_bg')
        with _Timer(2, 'recovery_interp'):
            rec = recovery_interp(S)
        # Explicit domain test.  The correctors extrapolate outside the
        # training set; these two objects say whether a query lies inside the
        # convex hull of the training samples, i.e. whether the prediction is
        # an interpolation between real models or an extension beyond them.
        # This is the criterion the Delaunay interpolator applied implicitly,
        # by returning not-a-number outside the hull, before the estimator was
        # replaced.  It is now reported rather than silently applied, so the
        # caller decides whether to trust an extrapolated correction.
        with _Timer(2, 'hull tests (3-axis, 4-axis)'):
            hull3 = _HullTest(cm.X0 if hasattr(cm, 'X0') else None)
            hull4 = _HullTest(cm4.X0 if hasattr(cm4, 'X0') else None)
        _GRIDS = dict(cm=cm, cm4=cm4, midx=midx, ff=ff, fbg=fbg, rec=rec,
                      hull3=hull3, hull4=hull4, nsamp=len(S['cm']))
        vlog(0, "Correctors built (%d training samples)." % _GRIDS['nsamp'])
    else:
        vlog(1, "Correctors already built, reusing cached grids.")
    return _GRIDS


def _fit_profile(keep, f, sig2, dpc):
    """Profile-likelihood thinbody fit over the T grid on the kept bands.
    Returns (T, M, T_err, M_fit_err)."""
    g = _Gall[keep]
    A_T = (f / sig2) @ g / ((1.0 / sig2) @ (g * g))
    resid = np.sum(((f[:, None] - A_T[None, :] * g) ** 2) / sig2[:, None], axis=0)
    j = np.argmin(resid); T = _TG[j]; A = A_T[j]; c2 = resid[j]
    M = A * OPACITY_K * dpc**2
    within = (resid - c2) <= 1.0
    Terr = 0.5 * (_TG[within].max() - _TG[within].min()) if within.any() else np.nan
    Mw = A_T[within] * OPACITY_K * dpc**2
    Mfit = 0.5 * (Mw.max() - Mw.min()) if within.any() else np.nan
    return T, M, Terr, Mfit


def correct_cloud(inp, distance):
    """Correct one cloud.  inp: concatenated catalogue path, or [cat, add] pair.
    Returns (table, summary): the per-source table for QUALITY=='ok' sources
    inside the model hull (initial reported SED fit; the correction factor C_M
    decomposed into flux-recovery x background-subtraction x temperature-bias
    mass factors; the corrected true mass; and the recovery-corrected refit),
    and summary medians (R_foot over all ok sources).

    Also computes, alongside the primary (three-observable) C_M/M_CORR, an
    additional four-observable correction C_M_4D/M_CORR_4D (paper's adopted
    fourth observable, physical size -- see cm_size_interp()), using this
    function's own `distance` argument to convert each source's footprint
    radius to physical units. Both are reported; the three-observable C_M
    remains the primary output column for backward compatibility, pending a
    decision on which to treat as default -- see the paper's validation
    comparison (Table~tab:validation) before choosing.
    """
    G = build_grids(); cmf = G['cm']; cm4f = G['cm4']; midx = G['midx']
    ff = G['ff']; fbg = G['fbg']; rec = G['rec']
    hull3 = G['hull3']; hull4 = G['hull4']
    vlog(0, "Correcting cloud: %s  (distance=%.1f pc)" % (inp, distance))
    load_getsf._dist = distance
    g = load_getsf(inp)
    ok = source_mask(g)
    vlog(1, "  source_mask: %d/%d sources pass QUALITY+hull filtering"
         % (int(np.sum(ok)), len(ok)))
    sed = g['sed']
    lg = lambda x: np.log10(x)

    out = defaultdict(list)
    n_no4d = 0
    for i in np.where(ok)[0]:
        if VERBOSE >= 2:
            vlog(2, "  source NO=%s: Sigma=%.3e conc=%.3f rfoot=%.3f"
                 % (g['NO'][i], g['Nbg'][i], g['conc'][i], g['rfoot'][i]))
        Sig, cc, rf = g['Nbg'][i], g['conc'][i], g['rfoot'][i]
        q = [[lg(Sig), lg(cc), lg(rf)]]
        cM = cmf(q)[0]
        # Domain diagnostics.  LocalLinear always returns a value, because it
        # fits a plane to the nearest training samples and extrapolates outside
        # their convex hull, bounding the result to the range those samples
        # span.  That is deliberate, but it means a source far outside the
        # training set is no longer distinguishable from one inside it by the
        # value alone.  Until version 6 the pipeline used a Delaunay
        # interpolator, which returned not-a-number outside the hull, and the
        # validation tests silently discarded those sources; the estimator
        # change therefore removed an implicit reliability filter without
        # replacing it.  These two diagnostics restore the distinction
        # explicitly:
        #   EXTRAP  the local plane had to be bounded to the range of its own
        #           neighbours, i.e. the query lies outside them
        #   SUPPORT the standardised distance to the furthest of the
        #           LL_NEIGHBOURS training samples used, so larger means the
        #           prediction rests on more distant analogues
        _in3 = bool(hull3(q)[0])
        _ex3 = bool(np.asarray(getattr(cmf, '_flag', [False]))[0])
        _sp3 = float(np.asarray(getattr(cmf, '_support', [np.nan]))[0])
        if not np.isfinite(cM):
            continue                                   # outside the model hull
        C_M = 10 ** cM
        Fflux = 10 ** ff(q)[0]; Fbg = 10 ** fbg(q)[0]
        Ftemp = C_M / (Fflux * Fbg) if np.isfinite(Fflux) and np.isfinite(Fbg) else np.nan
        if CAP_FFLUX and MIN_FRAC and np.isfinite(Fflux) and Fflux > 1.0 / MIN_FRAC:
            C_M *= (1.0 / MIN_FRAC) / Fflux          # cap flux recovery at 1/MIN_FRAC:
            Fflux = 1.0 / MIN_FRAC                    # scale C_M down consistently
        if CLAMP_CM:                                   # physical floor: every mechanism
            C_M = max(C_M, 1.0)                        # only adds mass, so C_M >= 1 and
            Fflux = max(Fflux, 1.0); Fbg = max(Fbg, 1.0)   # each factor >= 1; never let
            if np.isfinite(Ftemp):                     # the corrector reduce a mass
                Ftemp = max(Ftemp, 1.0)
        Minit = g['mass'][i]; Merr = g['mass_err'][i]
        Mcorr = Minit * C_M
        wv = midx([[lg(Minit), lg(Sig), lg(cc), lg(rf)]])[0]
        Cmass = 10 ** wv / Minit if np.isfinite(wv) else np.nan

        # four-observable correction (paper's adopted fourth observable):
        # real source's footprint radius, converted to physical units (pc)
        # at this cloud's distance, matched against the model grid's own
        # physical truncation radius R_BE. NaN if the footprint radius is
        # not finite/positive (e.g. a malformed FOOA/FOOB pair) -- falls
        # back gracefully, the primary C_M/M_CORR above is unaffected.
        # fourth axis is ANGULAR: no conversion by distance (see cm_size_interp)
        foot_as = g['foot_radius_as'][i]
        C_M_4D = M_CORR_4D = np.nan
        _ex4, _sp4, _in4 = True, np.nan, False
        if np.isfinite(foot_as) and foot_as > 0:
            q4 = [[lg(Sig), lg(cc), lg(rf), lg(foot_as)]]
            cM4 = cm4f(q4)[0]
            _in4 = bool(hull4(q4)[0])
            _ex4 = bool(np.asarray(getattr(cm4f, '_flag', [False]))[0])
            _sp4 = float(np.asarray(getattr(cm4f, '_support', [np.nan]))[0])
            if np.isfinite(cM4):
                C_M_4D = 10 ** cM4
                if CLAMP_CM:
                    C_M_4D = max(C_M_4D, 1.0)
                M_CORR_4D = Minit * C_M_4D
        if not np.isfinite(C_M_4D):
            n_no4d += 1
        # recovery-corrected refit: boost each band by 1/f_b(Sigma,conc,R_foot), refit.
        # Use only bands that both pass selection AND have a valid recovery fraction
        # (the 160um hull is smaller due to absorption); >=2 such bands needed.
        sel = ((sed['FP'][i] > 0) & (sed['FPE'][i] > 0) & (sed['F'][i] > 0)
               & (sed['E'][i] > 0) & (np.abs(sed['S'][i]) > 1)
               & (sed['FP'][i]/sed['FPE'][i] >= 1) & (sed['F'][i]/sed['E'][i] > 1))
        fbnd = np.array([10 ** rec[w](q)[0] for w in SED_WAVES])
        keep = sel & np.isfinite(fbnd) & (fbnd > 0)
        Trec = Mrec = Trec_err = Mrec_err = np.nan
        if keep.sum() >= 2:
            fboost = sed['F'][i][keep] / fbnd[keep]
            sig2 = sed['E'][i][keep]**2 + (_ADD[keep] * fboost)**2
            Trec, Mrec, Trec_err, Mfit = _fit_profile(keep, fboost, sig2, distance)
            Mrec_err = np.hypot.reduce([Mfit, SYS_OPACITY_RELERR * Mrec,
                                        SYS_DUSTGAS_RELERR * Mrec])
        out['NO'].append(g['NO'][i]); out['ra'].append(g['ra'][i]); out['dec'].append(g['dec'][i])
        out['sigma'].append(Sig); out['conc'].append(cc); out['rfoot'].append(rf)
        out['T_init'].append(g['temp'][i]); out['T_init_err'].append(g['temp_err'][i])
        out['M_init'].append(Minit); out['M_init_err'].append(Merr)
        out['chi2'].append(g['chi2'][i]); out['dof'].append(g['dof'][i])
        out['F_flux'].append(Fflux); out['F_bg'].append(Fbg); out['F_temp'].append(Ftemp)
        out['C_M'].append(C_M); out['M_corr'].append(Mcorr)
        out['M_corr_err'].append(Merr * C_M); out['C_Mmass'].append(Cmass)
        out['C_M_4D'].append(C_M_4D); out['M_corr_4D'].append(M_CORR_4D)
        # domain diagnostics; see the comment where they are computed
        out['EXTRAP_3D'].append(_ex3); out['EXTRAP_4D'].append(_ex4)
        out['SUPPORT_3D'].append(_sp3); out['SUPPORT_4D'].append(_sp4)
        out['IN_HULL_3D'].append(_in3); out['IN_HULL_4D'].append(_in4)
        out['IN_DOMAIN'].append(_in3 and _in4)
        out['foot_as'].append(foot_as)
        out['T_recov'].append(Trec); out['T_recov_err'].append(Trec_err)
        out['M_recov'].append(Mrec); out['M_recov_err'].append(Mrec_err)
    table = {k: np.array(v) for k, v in out.items()}
    med = lambda a: np.nanmedian(a) if len(a) and np.isfinite(a).any() else np.nan
    inr = (table['M_init'] > MATCHED_RANGE[0]) & (table['M_init'] < MATCHED_RANGE[1]) \
        if len(table) else np.array([], bool)
    summary = dict(N=len(table.get('C_M', [])),
                   Rfoot_ok=med(g['rfoot'][ok]),          # R_foot over ALL ok sources
                   C_M=med(table.get('C_M', [])),
                   C_M_0p1_2=med(table['C_M'][inr]) if len(table) and inr.any() else np.nan,
                   C_Mmass=med(table.get('C_Mmass', [])),
                   F_flux=med(table.get('F_flux', [])), F_bg=med(table.get('F_bg', [])),
                   F_temp=med(table.get('F_temp', [])))
    vlog(1, "  4-axis correction unavailable (no valid footprint radius) for %d/%d sources"
         % (n_no4d, summary['N']))
    vlog(0, "Cloud corrected: %d sources, median C_M=%.3f, median R_foot=%.3f"
         % (summary['N'], summary['C_M'], summary['Rfoot_ok']))
    return table, summary


# getsf-style catalogue columns: (name, fmt, unit, description)
_OUTCOLS = [
    ('NO',         '%7d',    '',        'Source running number (from getsf)'),
    ('WCS_ACOOR',  '%13.7f', 'deg',     'WCS (J2000) alpha-coordinate'),
    ('WCS_DCOOR',  '%13.7f', 'deg',     'WCS (J2000) delta-coordinate'),
    ('SIGMA',      '%11.4e', 'cm^-2',   'Background column PEAK^BGF on the 13.5" H2 band'),
    ('CONC',       '%8.3f',  '',        'Peak-to-mean concentration PEAK^SRC/(FXT/area)'),
    ('RFOOT',      '%8.3f',  '',        'Footprint-to-size ratio sqrt(FOOA*FOOB)/sqrt(AFWHM*BFWHM)'),
    ('T_INIT',     '%8.3f',  'K',       'Initial (reported) SED dust temperature'),
    ('T_INIT_ERR', '%8.3f',  'K',       'Uncertainty of T_INIT'),
    ('M_INIT',     '%12.5e', 'Msun',    'Initial (reported) SED mass'),
    ('M_INIT_ERR', '%12.5e', 'Msun',    'Total uncertainty of M_INIT (fit + opacity + dust/gas)'),
    ('CHI2',       '%9.3f',  '',        'chi^2 of the initial SED fit'),
    ('DOF',        '%4d',    '',        'Degrees of freedom (N_bands - 2)'),
    ('QUALITY',    '%8s',    '',        'Fit quality: ok if CHI2 < DOF+1, else bad'),
    ('F_FLUX',     '%8.3f',  '',        'Flux-recovery mass factor (footprint flux loss)'),
    ('F_BG',       '%8.3f',  '',        'Background-subtraction mass factor FXSDbsl/FXSDbs (interp->true crater)'),
    ('F_TEMP',     '%8.3f',  '',        'Temperature-bias mass factor (residual; internal T gradient)'),
    ('C_M',        '%8.3f',  '',        'Total correction factor = F_FLUX * F_BG * F_TEMP = M_BE/M_reported'),
    ('M_CORR',     '%12.5e', 'Msun',    'Corrected true mass = M_INIT * C_M'),
    ('M_CORR_ERR', '%12.5e', 'Msun',    'Uncertainty of M_CORR'),
    ('C_M_4D',     '%8.3f',  '',        'Four-observable correction factor (+physical size; needs DISTANCE); paper-adopted'),
    ('M_CORR_4D',  '%12.5e', 'Msun',    'Corrected mass from C_M_4D = M_INIT * C_M_4D'),
    ('T_RECOV',    '%8.3f',  'K',       'Recovery-corrected SED temperature (flux-boosted refit)'),
    ('T_RECOV_ERR','%8.3f',  'K',       'Uncertainty of T_RECOV'),
    ('M_RECOV',    '%12.5e', 'Msun',    'Recovery-corrected SED mass (flux-boosted refit; = M_INIT*F_FLUX)'),
    ('M_RECOV_ERR','%12.5e', 'Msun',    'Uncertainty of M_RECOV'),
    ('C_MMASS',    '%8.3f',  '',        'Mass-indexed factor (distance-dependent; diagnostic only)'),
    ('FOOT_AS',    '%9.2f',  'arcsec',  'Footprint radius sqrt(FOOA*FOOB)/2, the fourth observable'),
    ('EXTRAP_3D',  '%10d',   '',        '1 if the three-observable prediction was bounded to its neighbours range, i.e. extrapolated'),
    ('EXTRAP_4D',  '%10d',   '',        '1 if the four-observable prediction was extrapolated, or no fourth observable existed'),
    ('SUPPORT_3D', '%11.4f', '',        'Standardised distance to the furthest training sample used by the three-observable fit'),
    ('SUPPORT_4D', '%11.4f', '',        'Standardised distance to the furthest training sample used by the four-observable fit'),
    ('IN_HULL_3D', '%11d',   '',        '1 if the three observables lie inside the convex hull of the training set'),
    ('IN_HULL_4D', '%11d',   '',        '1 if the four observables lie inside the convex hull of the training set'),
    ('IN_DOMAIN',  '%10d',   '',        '1 if both hull tests pass; the recommended reliability filter'),
]
_KEYMAP = dict(NO='NO', WCS_ACOOR='ra', WCS_DCOOR='dec', SIGMA='sigma', CONC='conc',
    RFOOT='rfoot', T_INIT='T_init', T_INIT_ERR='T_init_err', M_INIT='M_init',
    M_INIT_ERR='M_init_err', CHI2='chi2', DOF='dof', F_FLUX='F_flux', F_BG='F_bg',
    F_TEMP='F_temp', C_M='C_M', M_CORR='M_corr', M_CORR_ERR='M_corr_err',
    C_M_4D='C_M_4D', M_CORR_4D='M_corr_4D',
    T_RECOV='T_recov', T_RECOV_ERR='T_recov_err', M_RECOV='M_recov',
    M_RECOV_ERR='M_recov_err', C_MMASS='C_Mmass', FOOT_AS='foot_as',
    EXTRAP_3D='EXTRAP_3D', EXTRAP_4D='EXTRAP_4D', SUPPORT_3D='SUPPORT_3D',
    SUPPORT_4D='SUPPORT_4D', IN_HULL_3D='IN_HULL_3D', IN_HULL_4D='IN_HULL_4D',
    IN_DOMAIN='IN_DOMAIN')


def write_catalog(table, path, cloud, distance):
    """Write the per-source corrected catalogue in getsf style (aligned, '!' header)."""
    n = len(table['NO']) if len(table) else 0
    lines = []
    bar = '!' + '_' * 100
    lines += [bar, '!', '! MASS-CORRECTED PRESTELLAR CORE CATALOGUE',
              '! cloud: %s   distance: %g pc   sources: %d' % (cloud, distance, n),
              '! Correction factor C_M = M_BE/M_reported on distance-invariant',
              '! observables (SIGMA, CONC, RFOOT), decomposed into three',
              '! multiplicative mass factors: C_M = F_FLUX * F_BG * F_TEMP',
              '! (flux recovery x background subtraction x temperature bias).',
              '! INIT = reported SED fit; CORR = M_INIT * C_M (true mass).', '!', bar,
              '!', '! TABULATED QUANTITIES:', '!']
    for k, (nm, fmt, unit, desc) in enumerate(_OUTCOLS, 1):
        lines.append('! %2d %-12s %-8s %s' % (k, nm, unit, desc))
    lines += [bar, '!', '! NOTES ON nan VALUES:',
              '!   T_RECOV, M_RECOV (+ errors): the recovery-corrected refit boosts each',
              '!     SED band by its model flux-recovery fraction and refits; it needs at',
              '!     least two bands with a valid recovery fraction.  A few sources sit',
              '!     outside the per-band recovery-interpolation hull (mainly at 160 um,',
              '!     which can be in absorption) in more than one band and cannot be',
              '!     refit -> nan.  M_CORR (= M_INIT * C_M) is unaffected and always finite.',
              '!   C_MMASS: distance-dependent, mass-indexed diagnostic from a 4-D',
              '!     interpolation (mass + the three observables); its hull is smaller',
              '!     than the 3-D C_M corrector, so sources just outside it are nan.',
              '!     It is diagnostic only and not used in the correction.',
              '!   All other columns (SIGMA..M_CORR_ERR) are finite for every listed',
              '!     source: only QUALITY==ok sources inside the C_M hull are written,',
              '!     and fit errors are defined down to 2-band (DOF=0) fits.', '!', bar, '!']
    # column-name header, aligned to the widths implied by the formats
    widths = [max(len(nm), _fmt_width(fmt)) for nm, fmt, _, _ in _OUTCOLS]
    hdr = '! ' + ' '.join('%*s' % (w, nm) for (nm, *_), w in zip(_OUTCOLS, widths))
    lines.append(hdr)
    for i in range(n):
        row = []
        for (nm, fmt, _, _), w in zip(_OUTCOLS, widths):
            v = table[_KEYMAP[nm]][i] if nm != 'QUALITY' else \
                ('ok' if table['chi2'][i] < table['dof'][i] + 1 else 'bad')
            if nm == 'QUALITY':
                s = fmt % v
            elif fmt.endswith('d'):
                s = fmt % int(v)
            else:
                s = (fmt % v) if np.isfinite(v) else 'nan'
            row.append('%*s' % (w, s))
        lines.append('  ' + ' '.join(row))
    open(path, 'w').write('\n'.join(lines) + '\n')
    return n


_SUMMARY_COLS = [
    ('cloud',      '%-16s'),
    ('distance_pc','%10.1f'),
    ('N',          '%6d'),
    ('N_no4D',     '%8d'),
    ('Rfoot_med',  '%10.3f'),
    ('C_M_med',    '%9.3f'),
    ('C_M_0.1_2',  '%10.3f'),
    ('F_flux_med', '%10.3f'),
    ('F_bg_med',   '%9.3f'),
    ('F_temp_med', '%10.3f'),
    ('C_M4D_med',  '%10.3f'),
]

def write_summary_row(path, cloud, distance, table, summary, n_no4d=None, append=True):
    """Append (or start) a running, space-separated summary-statistics table:
    one row per cloud/run, with the median correction factors and headline
    counts from `summary` (as returned by correct_cloud). Safe to call
    repeatedly across many clouds -- writes the header only if the file
    doesn't exist yet or append=False.

    This is the table to use for a quick look across a whole survey (are
    correction factors behaving sensibly cloud to cloud?), as opposed to
    write_catalog()'s per-source detail for one cloud at a time.
    """
    import os
    if n_no4d is None:
        n_no4d = int(np.sum(~np.isfinite(table.get('C_M_4D', []))))
    c_m4d = table.get('C_M_4D', [])
    c_m4d_med = float(np.nanmedian(c_m4d)) if len(c_m4d) and np.isfinite(c_m4d).any() else np.nan
    row = dict(cloud=cloud, distance_pc=distance, N=summary['N'], N_no4D=n_no4d,
               Rfoot_med=summary['Rfoot_ok'], C_M_med=summary['C_M'],
               **{'C_M_0.1_2': summary['C_M_0p1_2']},
               F_flux_med=summary['F_flux'], F_bg_med=summary['F_bg'],
               F_temp_med=summary['F_temp'], C_M4D_med=c_m4d_med)

    write_header = append and not os.path.exists(path)
    mode = 'a' if append else 'w'
    _widths = dict(cloud=16, distance_pc=10, N=6, N_no4D=8, Rfoot_med=10,
                   C_M_med=9, **{'C_M_0.1_2': 10}, F_flux_med=10, F_bg_med=9,
                   F_temp_med=10, C_M4D_med=10)
    with open(path, mode) as fh:
        if write_header or not append:
            fh.write("# Per-cloud summary statistics from mass_correction_pipeline_4d.py\n")
            fh.write("# N: sources corrected (QUALITY=ok, inside 3-axis hull)\n")
            fh.write("# N_no4D: of those, how many lack a valid 4-axis correction\n")
            fh.write("# Rfoot_med: median footprint factor phi over all ok sources\n")
            fh.write("# C_M_med, C_M4D_med: median 3-axis / 4-axis correction factor\n")
            fh.write("# C_M_0.1_2: median 3-axis C_M restricted to M_init in 0.1-2 Msun\n")
            fh.write("# F_flux_med, F_bg_med, F_temp_med: median decomposed factors\n")
            fh.write('#' + ' '.join(nm.rjust(_widths[nm]) for nm, _ in _SUMMARY_COLS)[1:] + '\n')
        fh.write(' '.join(fmt % row[nm] for nm, fmt in _SUMMARY_COLS) + '\n')
    vlog(0, "Summary row written to %s (cloud=%s)" % (path, cloud))


def _fmt_width(fmt):
    import re as _re
    m = _re.search(r'%(\d+)', fmt)
    return int(m.group(1)) if m else 8


# ==============  correction-factor node-LOO (chat-9 C1)  ====================
def run_cmloo():
    S = build_recovery_samples()
    nodes = np.array([str(n) for n in S['node']])
    logX = np.column_stack([np.log10(S['sig']), np.log10(S['conc']), np.log10(S['cfoot'])])
    logY = np.log10(S['cm']); truth = S['mrec']; mbe = S['mbe']
    errs = []
    for u in np.unique(nodes)[::CMLOO_STRIDE]:
        held = nodes == u
        if (~held).sum() < 10:
            continue
        f = LinearNDInterpolator(logX[~held], logY[~held])
        pc = f(logX[held]); good = np.isfinite(pc)
        if good.any():
            errs.append(np.median(np.abs(truth[held][good] * 10**pc[good] / mbe[held][good] - 1)))
    errs = np.array(errs)
    print("=== correction-factor node-LOO (stride=%d) ===" % CMLOO_STRIDE)
    print("  nodes %d   median-of-node-medians |dM|/M: %.0f%%  (chat-9: 55%%)"
          % (len(errs), 100 * np.median(errs)))


# ===============  mass-corrector node-LOO on catalog_rec  ===================
def run_loo():
    """LEGACY / SUPERSEDED: predates the self-consistent per-band footprint
    model (uses FWHMSDbs/conc_peakmean axes directly from a recovery-table
    dump, not the cumulative-fraction footprint solver of build_recovery_
    samples above). Not used for any of the paper's reported validation
    numbers; kept for reference only. Its internal 'm4'/'4D' local variables
    are a THIRD, unrelated use of that name (a mass+concentration "hybrid"
    refinement), distinct from both mass_idx_interp (the rejected
    mass-indexed approach) and cm_size_interp (the paper's adopted fourth
    observable) above -- do not confuse the three.
    """
    grid = load_catrec(CATREC)
    n = len(grid['M_BE']); Mt = grid['M_BE']
    mrec = grid['M_SED3bs_rec']; sig = grid['SD_emb']
    fw = grid['FWHMSDbs']; conc = grid['conc_peakmean']; frac = grid['frac_rec']

    def build_hybrid(mask):
        m, s, f, c, fr = mrec[mask], sig[mask], fw[mask], conc[mask], frac[mask]
        ok = (np.isfinite(m) & (m > FRAC_FLOOR) & np.isfinite(f) & (f > 0)
              & np.isfinite(s) & (s > 0) & (fr > FRAC_FLOOR))
        y = np.log10(Mt[mask][ok])
        base = np.column_stack([np.log10(m[ok]), np.log10(s[ok]), np.log10(f[ok])])
        f3 = LinearNDInterpolator(base, y); f4 = None
        good = np.isfinite(c[ok]) & (c[ok] > 0)
        if good.sum() > 20:
            f4 = LinearNDInterpolator(np.column_stack([base[good], np.log10(c[ok][good])]), y[good])
        return f3, f4

    def build_inv(mask):
        m, s, c, fr = mrec[mask], sig[mask], conc[mask], frac[mask]
        ok = (np.isfinite(m) & (m > FRAC_FLOOR) & np.isfinite(s) & (s > 0)
              & np.isfinite(c) & (c > 0) & (fr > FRAC_FLOOR))
        return LinearNDInterpolator(
            np.column_stack([np.log10(m[ok]), np.log10(s[ok]), np.log10(c[ok])]),
            np.log10(Mt[mask][ok]))

    m3, m4, mh, idxh = [], [], [], []
    for x in range(n):
        keep = np.ones(n, bool); keep[x] = False
        f3, f4 = build_hybrid(keep)
        if mrec[x] <= 0:
            continue
        q = [np.log10(mrec[x]), np.log10(sig[x]), np.log10(fw[x])]
        v = f3([q])[0]
        if np.isfinite(v):
            m3.append(10**v / Mt[x])
        if f4 is not None and conc[x] > 0:
            v4 = f4([q + [np.log10(conc[x])]])[0]
            if np.isfinite(v4):
                m4.append(10**v4 / Mt[x])
        v = (f4([q + [np.log10(conc[x])]])[0] if (f4 is not None and conc[x] > 0) else np.nan)
        if not np.isfinite(v):
            v = f3([q])[0]
        if np.isfinite(v):
            mh.append(10**v / Mt[x]); idxh.append(x)
    idxh = np.array(idxh)

    def rep(name, r):
        r = np.asarray(r, float)
        print("  %-22s median|err| %5.1f%%  within2x %4.1f%%  scored %3d"
              % (name, 100*np.median(np.abs(r-1)), 100*np.mean((r > 0.5) & (r < 2.0)), len(r)))
    unc = mrec[idxh] / Mt[idxh]
    print("=== Hybrid family node-LOO (sec:accuracy) ===")
    rep("uncorrected M_rec", unc); rep("3D (M_rec,Sig,FWHM)", m3)
    rep("4D (+conc_peakmean)", m4); rep("hybrid", mh)
    print("  hybrid makes worse: %.1f%%" % (100*np.mean(np.abs(np.array(mh)-1) > np.abs(unc-1))))
    inv, fr_ = [], []
    for x in range(n):
        keep = np.ones(n, bool); keep[x] = False
        f = build_inv(keep)
        if mrec[x] > 0 and conc[x] > 0:
            v = f([[np.log10(mrec[x]), np.log10(sig[x]), np.log10(conc[x])]])[0]
            if np.isfinite(v):
                inv.append(10**v / Mt[x]); fr_.append(frac[x])
    inv = np.array(inv); fr_ = np.array(fr_)
    print("\n=== InvariantCorrector node-LOO (sec:reliability) ===")
    rep("invariant overall", inv)
    for lo, hi in [(0, .05), (.05, .10), (.10, .20), (.20, .40), (.40, np.inf)]:
        m = (fr_ >= lo) & (fr_ < hi)
        if m.sum():
            lbl = "<%.2f" % hi if lo == 0 else (">%.2f" % lo if np.isinf(hi) else "%.2f-%.2f" % (lo, hi))
            print("  %-12s %4d nodes  median|dM|/M %3.0f%%" % (lbl, m.sum(), 100*np.median(np.abs(inv[m]-1))))


# =====================  test on an independent model grid  ==================
def _apply_caps(cm, fflux):
    """Apply the MIN_FRAC flux-recovery cap and the CLAMP_CM physical floor to an
    array of predicted C_M, given the (predicted) flux-recovery factors."""
    cm = np.asarray(cm, float).copy(); fflux = np.asarray(fflux, float)
    if CAP_FFLUX and MIN_FRAC:
        over = np.isfinite(fflux) & (fflux > 1.0 / MIN_FRAC)
        cm[over] *= (1.0 / MIN_FRAC) / fflux[over]
    if CLAMP_CM:
        cm = np.maximum(cm, 1.0)
    return cm


def run_test(cat_path, rectab_path):
    """Validate the method on an independent set of grid models, in line with the
    LOO tests: build the correction from the MAIN grid, then apply it to the test
    models (processed identically -- same recoverable-mass and observable
    construction) and compare the corrected mass to the true M_BE.  The test
    catalogue must be a grid-model catalogue (with per-band model fluxes) plus its
    recovery tables, i.e. the same inputs as the main grid."""
    G = build_grids(); cmf = G['cm']; ff = G['ff']         # corrector from main grid
    St = build_recovery_samples(cat_path, rectab_path)     # test-set samples
    lg = np.log10
    X = np.column_stack([lg(St['sig']), lg(St['conc']), lg(St['cfoot'])])
    pc = cmf(X); good = np.isfinite(pc)
    Cpred = _apply_caps(10 ** pc[good], 10 ** ff(X[good]))
    corr = St['mrec'][good] * Cpred / St['mbe'][good]
    unc = St['mrec'][good] / St['mbe'][good]
    n_nodes = len(np.unique([str(x) for x in St['node'][good]]))
    print("=== test on model grid: %s ===" % os.path.basename(cat_path))
    print("  test samples: %d (%d in corrector hull, %d nodes)   main-grid samples: %d"
          % (len(St['cm']), good.sum(), n_nodes, G['nsamp']))

    def rep(name, r):
        print("  %-14s median M_out/M_BE %.3f   median|err| %4.0f%%   within 30%% %3.0f%%   within 2x %3.0f%%"
              % (name, np.median(r), 100*np.median(np.abs(r-1)),
                 100*np.mean(np.abs(r-1) < 0.3), 100*np.mean((r > 0.5) & (r < 2.0))))
    if good.any():
        rep("uncorrected", unc)
        rep("corrected", corr)


def run_validate(cat_path, rectab_path):
    """Error-budget validation: apply the main-grid corrector to independent test
    models and tabulate the corrected-mass accuracy vs the recovered flux fraction
    'frac' (= truncation fraction).  frac->1 is the idealized full-footprint case
    where the flux-recovery factor is 1, so the residual error there is the
    background+temperature floor; the growth toward small frac isolates the
    flux-recovery contribution."""
    G = build_grids(); cmf = G['cm']; ff = G['ff']
    St = build_recovery_samples(cat_path, rectab_path)
    lg = np.log10
    X = np.column_stack([lg(St['sig']), lg(St['conc']), lg(St['cfoot'])])
    pc = cmf(X); good = np.isfinite(pc)
    Cpred = _apply_caps(10 ** pc[good], 10 ** ff(X[good]))
    r = St['mrec'][good] * Cpred / St['mbe'][good]
    fr = St['frac'][good]; ct = St['cm'][good]
    print("=== error-budget validation: %s   (CLAMP_CM=%s, MIN_FRAC=%s) ===" %
          (os.path.basename(cat_path), CLAMP_CM, MIN_FRAC))
    print("  %-13s %5s %9s %14s %9s %9s" %
          ("frac bin", "n", "C_M_true", "med(Mout/MBE)", "med|err|", "within2x"))
    for a, b in [(0.98, 1.0001), (0.90, 0.98), (0.70, 0.90), (0.50, 0.70),
                 (0.30, 0.50), (0.0, 0.30)]:
        m = (fr >= a) & (fr < b)
        if m.sum() >= 8:
            print("  %.2f-%-6.2f %5d %9.2f %14.3f %8.0f%% %8.0f%%" %
                  (a, b, m.sum(), np.median(ct[m]), np.median(r[m]),
                   100*np.median(np.abs(r[m]-1)), 100*np.mean((r[m] > .5) & (r[m] < 2))))
    m = fr >= 0.98
    if m.any():
        print("  --> floor (frac>=0.98, F_flux=1, bg+temp only): "
              "median M_out/M_BE %.3f  median|err| %.0f%%  within2x %.0f%%"
              % (np.median(r[m]), 100*np.median(np.abs(r[m]-1)),
                 100*np.mean((r[m] > .5) & (r[m] < 2))))


def main():
    import os
    global VERBOSE
    ap = argparse.ArgumentParser(
        description="Mass-correction pipeline (single cloud).  Corrects one "
                    "getsf cloud extraction and writes a per-source catalogue; "
                    "loop over clouds with run_clouds.py.")
    ap.add_argument('--verbose', type=int, default=0, choices=(0, 1, 2),
                    help='0=main steps only (default), 1=+deeper comments, '
                         '2=+most detailed (per-node/per-source progress)')
    sub = ap.add_subparsers(dest='cmd', required=True)
    pc = sub.add_parser('correct', help='correct one cloud and write its catalogue')
    pc.add_argument('inputs', nargs='+',
                    help='concatenated getsf catalogue, OR the raw pair "cat add"')
    pc.add_argument('--distance', type=float, required=True, help='cloud distance [pc]')
    pc.add_argument('--out', default=None, help='output catalogue path')
    pc.add_argument('--cloud', default=None, help='cloud name for the header')
    pc.add_argument('--summary-table', default=None,
                    help='append a summary-statistics row to this file (space-'
                         'separated; header written automatically if the file '
                         'is new). Useful when looping over many clouds to get '
                         'one running table instead of digging through per-'
                         'source catalogues.')
    sub.add_parser('loo',   help='mass-corrector node-LOO on the model grid')
    sub.add_parser('cmloo', help='correction-factor node-LOO on the model grid')
    pt = sub.add_parser('test', help='validate on an independent grid-model catalogue (in line with LOO)')
    pt.add_argument('catalog', help='test grid-model catalogue (with per-band model fluxes), like the main grid catalogue')
    pt.add_argument('rectab', help='recovery tables for the test catalogue')
    pv = sub.add_parser('validate', help='error budget vs truncation fraction (frac=1 floor, flux-recovery growth)')
    pv.add_argument('catalog', help='test grid-model catalogue')
    pv.add_argument('rectab', help='recovery tables for the test catalogue')
    a = ap.parse_args()
    VERBOSE = a.verbose

    if a.cmd == 'loo':
        run_loo()
    elif a.cmd == 'cmloo':
        run_cmloo()
    elif a.cmd == 'test':
        run_test(a.catalog, a.rectab)
    elif a.cmd == 'validate':
        run_validate(a.catalog, a.rectab)
    elif a.cmd == 'correct':
        inp = a.inputs[0] if len(a.inputs) == 1 else a.inputs
        cloud = a.cloud or os.path.basename(a.inputs[0]).split('.')[0]
        out = a.out or (cloud + '.corrected.cat')
        table, s = correct_cloud(inp, a.distance)
        n = write_catalog(table, out, cloud, a.distance)
        print("%s: corrected %d sources -> %s" % (cloud, n, out))
        print("  medians  R_foot(ok) %.2f  C_M %.2f = F_flux %.2f x F_bg %.3f x F_temp %.3f"
              % (s['Rfoot_ok'], s['C_M'], s['F_flux'], s['F_bg'], s['F_temp']))
        if a.summary_table:
            write_summary_row(a.summary_table, cloud, a.distance, table, s)


if __name__ == '__main__':
    main()
