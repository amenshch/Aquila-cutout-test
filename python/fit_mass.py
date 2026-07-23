"""
fit_mass.py -- mass-correction pipeline for Herschel SED-derived prestellar
core masses, built on the RADMC-3D critical-BE-sphere model grid.

Two ways to use the grid to recover the true mass M_BE from an observed core:

  (1) MassCorrector  -- the production method.  The grid yields a smooth
      bias map  f(M_SED, Sigma) = M_BE / M_SED  (crater Phi3).  Given an
      observed crater SED mass and the local cloud surface density, it
      interpolates the correction factor and returns M_true = f * M_SED.
      The bias is a slowly varying ratio (~1), so it interpolates to ~1%
      on the factor-2 grid -- consistent with the grid-spacing convergence
      test.

  (2) ForwardFitter  -- a Levenberg-Marquardt forward fit of the raw image
      observables (central intensities + FWHMs) for (T_BE, rho_c) at the
      supplied Sigma, with the mass taken analytically from the fitted
      (T, rho).  PROVIDED FOR COMPARISON ONLY: the raw observables vary
      steeply across a cell and interpolate at the several-percent level,
      and inverting them amplifies that error, so on the factor-2 grid this
      route recovers mass only to ~15-25% rms.  It is *not* the production
      path.  Documented here so the choice is evidence-based.

Validation (`validate`) is a hold-out test: the subdivision (cell-centre)
models are fed as fake observations, the fitter uses only the production
grid, and the recovered mass is compared with the known true M_BE.
"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import least_squares

COLS = ['n','i','j','k','SD_emb','T_BE','rho_BE','M_BE','R_BE_as','R_BE_au',
        'R_BE_pc','Td_avg','Td_emb','T_SED4bs','M_SED4bs','Chi2_4bs','T_SED4bsl',
        'M_SED4bsl','Chi2_4bsl','T_SED3bs','M_SED3bs','Chi2_3bs','T_SED3bsl',
        'M_SED3bsl','Chi2_3bsl','ICSDbs','IC070bs','IC100bs','IC160bs','IC250bs',
        'IC350bs','IC500bs','FWHMSDbs','FWHM070bs','FWHM100bs','FWHM160bs',
        'FWHM250bs','FWHM350bs','FWHM500bs','ICSDbsl','IC070bsl','IC100bsl',
        'IC160bsl','IC250bsl','IC350bsl','IC500bsl','FWHMSDbsl','FWHM070bsl',
        'FWHM100bsl','FWHM160bsl','FWHM250bsl','FWHM350bsl','FWHM500bsl']


def load_catalog(path):
    """Read a collect-format catalog into a dict of float arrays.

    Reads the base COLS by position.  If the file also carries the appended
    recoverable-mass columns (added by add_recoverable_mass.py):
        M_SED3bs_rec  M_SED3bsl_rec  frac_rec  FWHM_rec
        conc_footfwhm  conc_peakmean  conc_slope
    those are read from the trailing columns as well.  Files without them load
    exactly as before.
    """
    rows = [ln.split() for ln in open(path)
            if not ln.lstrip().startswith('#') and ln.strip()]
    out = {c: np.array([r[i] for r in rows], float) for i, c in enumerate(COLS)}

    ncol = len(rows[0]) if rows else len(COLS)
    if ncol >= len(COLS) + 7:
        # trailing recoverable-mass block (last 7 columns, fixed order)
        rec_cols = ['M_SED3bs_rec', 'M_SED3bsl_rec', 'frac_rec', 'FWHM_rec',
                    'conc_footfwhm', 'conc_peakmean', 'conc_slope']
        for off, name in enumerate(rec_cols):
            idx = ncol - 7 + off
            out[name] = np.array([r[idx] for r in rows], float)
    return out


class MassCorrector:
    """Production correction: M_true = M_SED * f(M_SED, Sigma)."""

    def __init__(self, cat, sed_mass='M_SED3bsl'):
        self.sed_mass = sed_mass
        x = np.column_stack([np.log10(cat[sed_mass]), np.log10(cat['SD_emb'])])
        self._logf = LinearNDInterpolator(x, np.log10(cat['M_BE'] / cat[sed_mass]))

    def correct(self, m_sed_obs, sigma_obs):
        """Return the bias-corrected (true) mass for an observed core."""
        lf = self._logf([[np.log10(m_sed_obs), np.log10(sigma_obs)]])[0]
        return np.nan if not np.isfinite(lf) else m_sed_obs * 10 ** lf


class ObservableCorrector:
    """Operational correction from interpolated-background quantities only.

    The crater background cannot be measured on real data, so the correction is
    keyed on what an interpolated-background extraction actually provides: the
    interp-bg SED mass M_SED3bs, the cloud surface density, and the source size
    (FWHM).  The size is essential -- the interpolated background over-subtracts
    by (edge - crater) x footprint, which grows with angular extent, so two
    cores with the same (M_SED3bs, Sigma) but different sizes have different
    true masses.

    IMPORTANT: the FWHM coordinate must be measured on the *same* (interpolated)
    background as the observation; build the map from the model's interp-bg FWHM
    (FWHM_int), not the crater FWHM.

    Note on axis scaling: all three coordinates enter as log10, matching
    RecoverableCorrector.  An earlier version used a linear FWHM axis; a
    leave-one-out test on the model grid favors log10 on every metric (median
    error 2.5% vs 3.4%, worst case 276% vs 529%, and marginally wider hull
    coverage), which is expected since the grid samples size geometrically.
    Any external hull test must use the same log space.

    The FWHM key must match the resolution the grid was built at.  The default
    FWHMSDbs is the 161 um surface-density band (13.5 arcsec), which is the band
    the extraction and the model grid share; passing a single-band FWHM such as
    FWHM250bs silently mixes in that band's coarser beam (18.2 arcsec) and
    biases the size coordinate.

    mode='direct'  : f(M_SED3bs, Sigma, FWHM) -> M_BE                  (~3% rms)
    mode='twostep' : (M_SED3bs, Sigma, FWHM) -> M_SED3bsl, then the crater
                     map M_BE/M_SED3bsl                                 (~2.5% rms)
    """

    def __init__(self, cat, fwhm_key='FWHMSDbs', mode='direct'):
        self.mode, self.fwhm_key = mode, fwhm_key
        fw = np.asarray(cat[fwhm_key], float)
        ok = (np.isfinite(fw) & (fw > 0)
              & np.isfinite(cat['M_SED3bs']) & (cat['M_SED3bs'] > 0)
              & np.isfinite(cat['SD_emb']) & (cat['SD_emb'] > 0))
        if not ok.all():
            cat = {k: np.asarray(v)[ok] for k, v in cat.items()}
            fw = fw[ok]
        X = np.column_stack([np.log10(cat['M_SED3bs']),
                             np.log10(cat['SD_emb']), np.log10(fw)])
        if mode == 'direct':
            self._f = LinearNDInterpolator(X, np.log10(cat['M_BE']))
        elif mode == 'twostep':
            self._g = LinearNDInterpolator(X, np.log10(cat['M_SED3bsl']))
            self._crater = MassCorrector(cat, 'M_SED3bsl')
        else:
            raise ValueError("mode must be 'direct' or 'twostep'")

    def correct(self, m_sed_obs, sigma_obs, fwhm_obs):
        if not (m_sed_obs > 0 and sigma_obs > 0 and fwhm_obs > 0):
            return np.nan
        p = [[np.log10(m_sed_obs), np.log10(sigma_obs), np.log10(fwhm_obs)]]
        if self.mode == 'direct':
            v = self._f(p)[0]
            return np.nan if not np.isfinite(v) else 10 ** v
        lm = self._g(p)[0]
        if not np.isfinite(lm):
            return np.nan
        return self._crater.correct(10 ** lm, sigma_obs)


class RecoverableCorrector:
    """Environment-aware correction: f(M_SED_rec, Sigma, concentration) -> M_BE.

    This is the production corrector for real fields.  It differs from
    ObservableCorrector in two ways that together remove the environment
    dependence of the getsf flux loss:

      1. The mass coordinate is the *recoverable* SED mass, M_SED3bs_rec =
         frac_rec * M_SED3bs, where frac_rec is the fraction of a model's flux
         that getsf recovers against the Konyves cirrus fluctuations at the
         model's own column (computed by add_recoverable_mass.py).  In dense
         fields getsf under-measures flux; folding that loss into the grid means
         the observed (already flux-reduced) M_SED maps back to M_BE correctly.

      2. A concentration axis and the source FWHM together break the
         recoverable-mass fold: at fixed column
         the coldest cores lose so much flux that a higher-mass cold core
         recovers less mass than a lower-mass warm one, so (mass, Sigma) alone
         is degenerate.  The footprint-to-FWHM ratio distinguishes a genuine
         compact core from the surviving centre of an over-resolved diffuse
         core, and it is measured identically in the grid (conc_footfwhm) and in
         a getsf catalogue (FOOA/AFWHM at the 161 um column-density band), so it
         enters the interpolation without a unit/definition mismatch.

    The observed concentration for a real source is FOOA/AFWHM at the 161 um
    (surface-density) band from its getsf catalogue row.

    The grid must carry the recoverable columns (load a catalog produced by
    add_recoverable_mass.py).  Models with frac_rec below `frac_floor` (nearly
    undetectable) and non-finite concentrations are dropped from the interpolant.
    """

    def __init__(self, cat, mass_key='M_SED3bs_rec', conc_key='conc_footfwhm',
                 fwhm_key='FWHMSDbs', use_fwhm=True, frac_floor=1e-4):
        if mass_key not in cat or conc_key not in cat:
            raise KeyError(
                "catalog lacks recoverable columns; load a catalog written by "
                "add_recoverable_mass.py (mass_key=%r, conc_key=%r)"
                % (mass_key, conc_key))
        self.mass_key, self.conc_key = mass_key, conc_key
        self.use_fwhm = use_fwhm and (fwhm_key in cat)
        self.fwhm_key = fwhm_key
        mrec = cat[mass_key]
        conc = cat[conc_key]
        keep = (mrec > frac_floor) & np.isfinite(conc) & (conc > 0) & (cat['SD_emb'] > 0)
        if self.use_fwhm:
            fw = cat[fwhm_key]
            keep = keep & np.isfinite(fw) & (fw > 0)
        cols = [np.log10(mrec[keep]), np.log10(cat['SD_emb'][keep]),
                np.log10(conc[keep])]
        if self.use_fwhm:
            cols.append(np.log10(cat[fwhm_key][keep]))
        X = np.column_stack(cols)
        self._f = LinearNDInterpolator(X, np.log10(cat['M_BE'][keep]))

    def correct(self, m_sed_obs, sigma_obs, conc_obs, fwhm_obs=None):
        """Return the corrected (true) mass.

        m_sed_obs : observed interp-bg SED mass (getsf).
        sigma_obs : local cloud surface density (cm^-2); for a filament source
                    this legitimately includes the filament background.
        conc_obs  : observed concentration = FOOA/AFWHM at the 161 um band.
        fwhm_obs  : observed FWHM (arcsec) at the 161 um band (AFWHM).  Required
                    when the corrector was built with use_fwhm=True; it pins the
                    physical size so the recoverable-mass fold (a compact core
                    vs the surviving centre of a diffuse giant) is broken.
        """
        if not (m_sed_obs > 0 and sigma_obs > 0 and conc_obs > 0):
            return np.nan
        q = [np.log10(m_sed_obs), np.log10(sigma_obs), np.log10(conc_obs)]
        if self.use_fwhm:
            if not (fwhm_obs and fwhm_obs > 0):
                return np.nan
            q.append(np.log10(fwhm_obs))
        v = self._f([q])[0]
        return np.nan if not np.isfinite(v) else 10 ** v


class ForwardFitter:
    """LM forward fit of raw observables -> (T,rho) -> M.  Comparison only."""

    OBS = ['ICSDbs', 'FWHMSDbs', 'IC250bs', 'IC350bs', 'IC500bs',
           'FWHM250bs', 'FWHM350bs', 'FWHM500bs']

    def __init__(self, cat):
        self.cat = cat
        self._P = np.column_stack([np.log10(cat['SD_emb']),
                                   np.log10(cat['T_BE']),
                                   np.log10(cat['rho_BE'])])
        self._I = {o: LinearNDInterpolator(self._P, cat[o]) for o in self.OBS}
        self._C = np.median(cat['M_BE'] / (cat['T_BE'] ** 1.5 / np.sqrt(cat['rho_BE'])))
        self._loT = (np.log10(cat['T_BE'].min()), np.log10(cat['T_BE'].max()))
        self._loR = (np.log10(cat['rho_BE'].min()), np.log10(cat['rho_BE'].max()))

    def _sigma(self, t, o):
        if o.startswith('IC'):
            return 0.05 * abs(t[o]) + 0.02 * max(abs(t['IC250bs']), 1e-3)
        return 0.05 * abs(t[o]) + 0.5  # FWHM, arcsec

    def fit(self, sigma_obs, observables):
        t, logS = observables, np.log10(sigma_obs)

        def resid(x):
            r = []
            for o in self.OBS:
                v = self._I[o]([[logS, x[0], x[1]]])[0]
                if not np.isfinite(v):
                    return np.full(len(self.OBS), 1e3)
                r.append((v - t[o]) / self._sigma(t, o))
            return np.array(r)

        d = sum(((self.cat[o] - t[o]) / self._sigma(t, o)) ** 2 for o in self.OBS)
        j0 = int(np.argmin(d))
        x0 = [np.clip(np.log10(self.cat['T_BE'][j0]), *self._loT),
              np.clip(np.log10(self.cat['rho_BE'][j0]), *self._loR)]
        s = least_squares(resid, x0, method='trf',
                          bounds=([self._loT[0], self._loR[0]],
                                  [self._loT[1], self._loR[1]]),
                          xtol=1e-12, ftol=1e-12)
        m = self._C * (10 ** s.x[0]) ** 1.5 / np.sqrt(10 ** s.x[1])
        return m, 10 ** s.x[0]


class InvariantCorrector:
    """(M_SED_rec, Sigma, conc_peakmean) -> M_BE, fully distance-invariant.

    The size axis of the other correctors is an angular FWHM, which requires a
    distance and a Gaussian size deconvolution to interpret.  Men'shchikov
    (2023) shows that deconvolution errs by factors up to ~20 for unresolved
    structures and up to ~6 for power-law profiles, i.e. exactly the
    marginally-resolved and non-Gaussian regime that distant clouds occupy.

    This corrector replaces the size axis with the peak-to-mean surface
    brightness over the source footprint (conc_peakmean).  Being a ratio of two
    intensities measured in the same map, it is distance-invariant by
    construction: no distance, no deconvolution.  The same operation applies
    unchanged to any cloud at any distance.

    Leave-one-out on the 517-node refined grid, recoverable path:

        FWHM (angular)      median |err| 10.2%   within 2x  92%   313 scored
        conc_peakmean       median |err| 14.0%   within 2x  82%   348 scored

    The invariant axis is nominally ~4 points less accurate in LOO, but that
    baseline uses the grid's own exactly-measured FWHM; on real data the FWHM
    path additionally carries the deconvolution error above, which is absent
    here.  It also reaches more sources (348 vs 313), since it does not lose
    marginally resolved cores.
    """

    def __init__(self, cat, conc_key='conc_peakmean',
                 mass_key='M_SED3bs_rec', frac_floor=1e-4):
        self.conc_key, self.mass_key = conc_key, mass_key
        m = np.asarray(cat[mass_key], float)
        sd = np.asarray(cat['SD_emb'], float)
        c = np.asarray(cat[conc_key], float)
        fr = np.asarray(cat.get('frac_rec', np.ones_like(m)), float)
        ok = (np.isfinite(m) & (m > frac_floor) & np.isfinite(sd) & (sd > 0)
              & np.isfinite(c) & (c > 0) & (fr > frac_floor))
        X = np.column_stack([np.log10(m[ok]), np.log10(sd[ok]), np.log10(c[ok])])
        self.f = LinearNDInterpolator(X, np.log10(np.asarray(cat['M_BE'], float)[ok]))

    def correct(self, m_rec_obs, sigma_obs, conc_obs):
        """Return corrected M_BE, or nan if outside the hull. No distance needed."""
        if not (m_rec_obs > 0 and sigma_obs > 0 and conc_obs > 0):
            return np.nan
        v = self.f([[np.log10(m_rec_obs), np.log10(sigma_obs), np.log10(conc_obs)]])[0]
        return 10 ** v if np.isfinite(v) else np.nan


class HybridRecoverableCorrector:
    """(M_rec, Sigma, FWHM [, concentration]) -> M_BE, with a 4D->3D fallback.

    Adding a concentration axis sharpens the interpolation but shrinks the
    convex hull, so a pure 4D map leaves a third of the queries unreachable.
    This class builds both and uses the 4D value where the query lies inside
    the 4D hull, falling back to 3D elsewhere.  Leave-one-out on the 517-node
    refined grid (313 scored):

        uncorrected M_rec   median |err| 52.7%   within 2x of truth  45%
        3D                  median |err| 10.2%   within 2x           92%
        4D (conc_peakmean)  median |err|  9.2%   within 2x           99%  (only 209 reachable)
        hybrid              median |err|  8.4%   within 2x           96%  (all 313)

    The hybrid keeps the coverage of the 3D map while recovering most of the
    tail improvement of the 4D one: the 90th-percentile error falls from 55%
    to 46%, and the fraction of nodes made WORSE than the uncorrected mass
    falls from 8% to 6%.  conc_peakmean works; conc_footfwhm does not
    (median 10.4%, no better than 3D).
    """

    def __init__(self, cat, fwhm_key='FWHMSDbs', conc_key='conc_peakmean',
                 mass_key='M_SED3bs_rec', frac_floor=1e-4):
        self.fwhm_key, self.conc_key, self.mass_key = fwhm_key, conc_key, mass_key
        m = np.asarray(cat[mass_key], float)
        fw = np.asarray(cat[fwhm_key], float)
        sd = np.asarray(cat['SD_emb'], float)
        fr = np.asarray(cat.get('frac_rec', np.ones_like(m)), float)
        ok = (np.isfinite(m) & (m > frac_floor) & np.isfinite(fw) & (fw > 0)
              & np.isfinite(sd) & (sd > 0) & (fr > frac_floor))
        base = [np.log10(m[ok]), np.log10(sd[ok]), np.log10(fw[ok])]
        y = np.log10(np.asarray(cat['M_BE'], float)[ok])
        self.f3 = LinearNDInterpolator(np.column_stack(base), y)

        self.f4 = None
        if conc_key in cat:
            c = np.asarray(cat[conc_key], float)[ok]
            good = np.isfinite(c) & (c > 0)
            if good.sum() > 20:
                cols = [b[good] for b in base] + [np.log10(c[good])]
                self.f4 = LinearNDInterpolator(np.column_stack(cols), y[good])

    def correct(self, m_rec_obs, sigma_obs, fwhm_obs, conc_obs=None):
        """Return the corrected M_BE, or nan if outside both hulls."""
        if not (m_rec_obs > 0 and sigma_obs > 0 and fwhm_obs > 0):
            return np.nan
        q = [np.log10(m_rec_obs), np.log10(sigma_obs), np.log10(fwhm_obs)]
        if self.f4 is not None and conc_obs is not None and conc_obs > 0:
            v = self.f4([q + [np.log10(conc_obs)]])[0]
            if np.isfinite(v):
                return 10 ** v
        v = self.f3([q])[0]
        return 10 ** v if np.isfinite(v) else np.nan

    def used_4d(self, m_rec_obs, sigma_obs, fwhm_obs, conc_obs=None):
        """True if the 4D map supplied the answer (useful for diagnostics)."""
        if self.f4 is None or conc_obs is None or conc_obs <= 0:
            return False
        if not (m_rec_obs > 0 and sigma_obs > 0 and fwhm_obs > 0):
            return False
        return bool(np.isfinite(self.f4([[np.log10(m_rec_obs), np.log10(sigma_obs),
                                          np.log10(fwhm_obs), np.log10(conc_obs)]])[0]))


def validate(grid_cat, holdout_cat, sed_mass='M_SED3bsl', fwhm_key='FWHMSDbs'):
    """Hold-out test on the subdivision models. Returns a dict of error arrays (%)."""
    crater = MassCorrector(grid_cat, sed_mass)
    obs_d = ObservableCorrector(grid_cat, fwhm_key, 'direct')
    obs_2 = ObservableCorrector(grid_cat, fwhm_key, 'twostep')
    ff = ForwardFitter(grid_cat)
    out = {'crater': [], 'obs_direct': [], 'obs_twostep': [], 'forward': []}
    for x in range(len(holdout_cat['n'])):
        Mt, Sig = holdout_cat['M_BE'][x], holdout_cat['SD_emb'][x]
        m = crater.correct(holdout_cat[sed_mass][x], Sig)
        if np.isfinite(m): out['crater'].append(100 * (m / Mt - 1))
        f = holdout_cat[fwhm_key][x]
        m = obs_d.correct(holdout_cat['M_SED3bs'][x], Sig, f)
        if np.isfinite(m): out['obs_direct'].append(100 * (m / Mt - 1))
        m = obs_2.correct(holdout_cat['M_SED3bs'][x], Sig, f)
        if np.isfinite(m): out['obs_twostep'].append(100 * (m / Mt - 1))
        m, _ = ff.fit(Sig, {o: holdout_cat[o][x] for o in ForwardFitter.OBS})
        if np.isfinite(m): out['forward'].append(100 * (m / Mt - 1))
    return {k: np.array(v) for k, v in out.items()}


if __name__ == '__main__':
    grid = load_catalog('cat4.txt')
    sub = load_catalog('sub4.txt')
    res = validate(grid, sub)

    def stat(e):
        return (f"N={len(e):2d}  mean {e.mean():+5.1f}%  rms {np.sqrt((e**2).mean()):5.1f}%  "
                f"max|.| {np.max(np.abs(e)):5.1f}%")

    print("Hold-out on 36 subdivision (cell-centre) models:\n")
    print("  crater map  M_BE/M_SED3bsl,Sigma   (ideal, model-only) :", stat(res['crater']))
    print("  observable  direct  (M_SED3bs,Sigma,FWHM)             :", stat(res['obs_direct']))
    print("  observable  two-step via crater mass                  :", stat(res['obs_twostep']))
    print("  LM forward fit of raw images (cross-check only)       :", stat(res['forward']))

