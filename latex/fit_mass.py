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
    """Read a collect-format catalog into a dict of float arrays."""
    rows = [ln.split() for ln in open(path)
             if not ln.lstrip().startswith('#') and ln.strip()]
    return {c: np.array([r[i] for r in rows], float) for i, c in enumerate(COLS)}


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

    mode='direct'  : f(M_SED3bs, Sigma, FWHM) -> M_BE                  (~3% rms)
    mode='twostep' : (M_SED3bs, Sigma, FWHM) -> M_SED3bsl, then the crater
                     map M_BE/M_SED3bsl                                 (~2.5% rms)
    """

    def __init__(self, cat, fwhm_key='FWHM250bs', mode='direct'):
        self.mode, self.fwhm_key = mode, fwhm_key
        X = np.column_stack([np.log10(cat['M_SED3bs']),
                             np.log10(cat['SD_emb']), cat[fwhm_key]])
        if mode == 'direct':
            self._f = LinearNDInterpolator(X, np.log10(cat['M_BE']))
        elif mode == 'twostep':
            self._g = LinearNDInterpolator(X, np.log10(cat['M_SED3bsl']))
            self._crater = MassCorrector(cat, 'M_SED3bsl')
        else:
            raise ValueError("mode must be 'direct' or 'twostep'")

    def correct(self, m_sed_obs, sigma_obs, fwhm_obs):
        p = [[np.log10(m_sed_obs), np.log10(sigma_obs), fwhm_obs]]
        if self.mode == 'direct':
            v = self._f(p)[0]
            return np.nan if not np.isfinite(v) else 10 ** v
        lm = self._g(p)[0]
        if not np.isfinite(lm):
            return np.nan
        return self._crater.correct(10 ** lm, sigma_obs)


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


def validate(grid_cat, holdout_cat, sed_mass='M_SED3bsl', fwhm_key='FWHM250bs'):
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

