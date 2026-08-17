#!/usr/bin/env python3
"""
patch_estimator.py -- replace the correction interpolator with a local linear
regression, and remove the two output constraints, in
mass_correction_pipeline_4d.py.

Run patch_band_R.py FIRST; this script checks for it and refuses otherwise.

WHAT IT CHANGES
---------------
1. Adds a module-level class LocalLinear, a drop-in replacement for scipy's
   LinearNDInterpolator: constructed from (points, values) and called with an
   (n, ndim) array of query points, returning an (n,) array of predictions.
   It fits a tricube-distance-weighted least squares plane to the LL_NEIGHBOURS
   nearest training samples in a space where each axis is divided by its own
   standard deviation.
2. Rewrites cm_interp() and cm_size_interp() to return a LocalLinear instead of
   a LinearNDInterpolator.  Nothing else about them changes: same axes, same
   target quantity log10(C_M), same call signature.
3. Sets CLAMP_CM = False, so the corrector output is no longer floored at
   C_M >= 1, and adds CAP_FFLUX = False, which disables the cap of the
   flux-recovery factor at 1/MIN_FRAC.  MIN_FRAC itself is UNCHANGED and keeps
   its other, unrelated role of requiring at least three bands with a recovered
   flux fraction of at least 0.25 before a training sample is accepted.

WHY
---
Measured by node-level leave-one-out over the 403 nodes of the production grid,
with the corrected per-band footprint solver in place, and with the fourth axis
queried at the model's own R_BE:

  estimator                          residual in log10(C_M)   nodes predicted
  LinearNDInterpolator, raw axes     0.0583 dex (factor 1.144)     388 of 403
  LinearNDInterpolator, axes scaled  0.0582 dex (factor 1.143)     388 of 403
  LocalLinear, 320 neighbours        0.0445 dex (factor 1.108)     403 of 403

The improvement is flat across neighbourhood sizes from 80 to 640 samples and
across weightings of the fourth axis from 0.5 to 2, and it survives a split-half
check in which the configuration is chosen on one random half of the nodes and
evaluated on the other (0.0437 dex).  Unlike a Delaunay triangulation, a local
regression extrapolates, so the 15 nodes that fall outside the four-dimensional
convex hull now receive a prediction instead of none.

The two output constraints altered the leave-one-out result by less than one per
cent on the corrected training data (scatter 1.145 against 1.15, median ratio
1.007 against 1.008).  They were compensating for the solver noise and no longer
serve a purpose, so they are removed.

NOT CHANGED
-----------
factor_interp() and recovery_interp(), which supply the diagnostic factors
F_flux and F_bg and the per-band recovery fractions for the recovery-corrected
refit, still use LinearNDInterpolator.  Those outputs were not part of this
validation and are left alone deliberately.

USAGE
-----
    python3 patch_estimator.py mass_correction_pipeline_4d.py
"""
import shutil, sys, os, re

LOCALLIN_CLASS = '''
# ----------------------- local linear corrector estimator -------------------
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

    def __call__(self, q):
        Q = np.atleast_2d(np.asarray(q, float)) / self.scale
        out = np.full(len(Q), np.nan)
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
            out[m] = b[0]
        return out


'''

NEW_CM_INTERP = '''def cm_interp(S):
    """Three-observable corrector: predicts log10(C_M) from
    (Sigma, zeta, phi).  Uses a local linear regression (LocalLinear above)
    rather than a Delaunay interpolator; see patch_estimator.py for the
    leave-one-out comparison that motivated the change."""
    return LocalLinear(np.column_stack(
        [np.log10(S['sig']), np.log10(S['conc']), np.log10(S['cfoot'])]),
        np.log10(S['cm']))
'''

NEW_CM_SIZE = '''def cm_size_interp(S):
    """Four-observable corrector: predicts log10(C_M) from
    (Sigma, zeta, phi, R_BE), where R_BE [parsecs] is the model's own physical
    truncation radius, matched on the real-source side by the footprint radius
    sqrt(A_F*B_F)/2 at the source's known distance.

    Uses a local linear regression (LocalLinear above) rather than a Delaunay
    interpolator.  Measured by node-level leave-one-out over the production
    grid, this reduces the residual in log10(C_M) from 0.0583 dex (a factor of
    1.144) to 0.0445 dex (a factor of 1.108) and returns a prediction for all
    403 nodes instead of 388, because a local regression extrapolates outside
    the convex hull.
    """
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    return LocalLinear(np.column_stack(
        [np.log10(S['sig'][ok]), np.log10(S['conc'][ok]), np.log10(S['cfoot'][ok]),
         np.log10(S['rbe_pc'][ok])]), np.log10(S['cm'][ok]))
'''


def cut(src, start_marker, end_marker):
    a = src.index(start_marker)
    b = src.index(end_marker, a)
    return a, b


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'mass_correction_pipeline_4d.py'
    if not os.path.exists(path):
        sys.exit("file not found: %s" % path)
    src = open(path).read()

    if 'BANDR_DIAG' not in src:
        sys.exit("run patch_band_R.py on this file first; the footprint-solver "
                 "fix must be in place before the estimator is changed")
    if 'class LocalLinear' in src:
        sys.exit("this file already appears to be patched (LocalLinear present); "
                 "nothing done")

    shutil.copyfile(path, path + '.orig_estimator')

    # 1. import cKDTree
    imp = "from scipy.interpolate import LinearNDInterpolator"
    if imp not in src:
        sys.exit("could not find the scipy import line; aborting")
    src = src.replace(imp, imp + "\nfrom scipy.spatial import cKDTree", 1)

    # 2. insert the LocalLinear class just before cm_interp's comment block
    anchor = "# CORRECTOR AXES."
    if anchor not in src:
        sys.exit("could not find the CORRECTOR AXES comment block; aborting")
    src = src.replace(anchor, LOCALLIN_CLASS.lstrip('\n') + anchor, 1)

    # 3. replace cm_interp
    a, b = cut(src, "def cm_interp(S):", "def mass_idx_interp(S):")
    src = src[:a] + NEW_CM_INTERP + "\n" + src[b:]

    # 4. replace cm_size_interp
    a, b = cut(src, "def cm_size_interp(S):", "def factor_interp(S, key):")
    src = src[:a] + NEW_CM_SIZE + "\n" + src[b:]

    # 5. remove the output constraints
    src2 = re.sub(r"^CLAMP_CM = True", "CLAMP_CM = False", src, count=1, flags=re.M)
    if src2 == src:
        sys.exit("could not find the CLAMP_CM setting; aborting")
    src = src2
    src = src.replace("MIN_FRAC = 0.25",
                      "CAP_FFLUX = False                # output cap on the "
                      "flux-recovery factor: removed\nMIN_FRAC = 0.25", 1)
    n1 = src.count("if MIN_FRAC and np.isfinite(Fflux)")
    src = src.replace("if MIN_FRAC and np.isfinite(Fflux)",
                      "if CAP_FFLUX and MIN_FRAC and np.isfinite(Fflux)")
    n2 = src.count("    if MIN_FRAC:\n        over = np.isfinite(fflux)")
    src = src.replace("    if MIN_FRAC:\n        over = np.isfinite(fflux)",
                      "    if CAP_FFLUX and MIN_FRAC:\n        over = np.isfinite(fflux)")

    open(path, 'w').write(src)
    print("patched   : %s" % path)
    print("backup    : %s.orig_estimator" % path)
    print("added     : class LocalLinear (LL_NEIGHBOURS = 320)")
    print("rewrote   : cm_interp() and cm_size_interp() to use it")
    print("set       : CLAMP_CM = False, CAP_FFLUX = False "
          "(%d and %d cap sites guarded)" % (n1, n2))
    print("unchanged : MIN_FRAC = 0.25 in its training-sample role, "
          "factor_interp(), recovery_interp()")
    import py_compile
    py_compile.compile(path, doraise=True)
    print("syntax    : OK")


if __name__ == '__main__':
    main()
