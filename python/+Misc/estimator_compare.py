#!/usr/bin/env python3
"""
estimator_compare.py -- compare estimators for the four-observable mass
corrector, under the same node-level leave-one-out protocol as
loo_solverfix.py.

MOTIVATION
----------
With the per-band footprint solver corrected (patch_band_R.py), the residual of
the four-observable corrector is a factor of 1.145 in the predicted correction
factor C_M, of which only a factor of 1.072 is the irreducible degeneracy of the
observables and a factor of 1.025 is pipeline numerical noise.  The remainder, a
factor of about 1.12, is interpolation error.  This script tests whether a
better estimator removes it.

ESTIMATORS TESTED
-----------------
  delaunay_raw      scipy LinearNDInterpolator on the raw log axes.  This is
                    what the production pipeline uses.
  delaunay_std      the same, with each log axis divided by its standard
                    deviation before triangulation, so the metric is isotropic.
  locallin_k<K>_w<W>  local linear regression: an ordinary least squares plane
                    fitted to the K nearest training samples, with tricube
                    distance weights, on standardized axes in which the R_BE
                    axis is additionally multiplied by W.  Unlike the Delaunay
                    estimators this extrapolates outside the convex hull, so it
                    also returns a value for the nodes the current method loses.
  rbf_n<K>_s<S>     scipy RBFInterpolator, thin-plate spline, restricted to the
                    K nearest samples, with smoothing parameter S.

PROTOCOL
--------
For each node, all of its samples are removed from the training set, the
estimator is rebuilt, and the prediction is made on that node's representative
sample -- the one whose footprint factor phi is closest to 1.84.  The Delaunay
estimators use the local-neighbourhood shortcut validated in the previous step:
piecewise linear interpolation uses only the simplex enclosing the query, so a
triangulation of the nearest NLOCAL samples gives an identical answer, and the
neighbourhood is escalated whenever a non-finite value is returned.

HONEST HYPERPARAMETER CHECK
---------------------------
Choosing K and W by inspecting the same leave-one-out numbers is a mild form of
selection on the test set.  The script therefore also reports a split-half
result: the best (K, W) is chosen on a random half of the nodes and evaluated on
the other half, with the split seeded for reproducibility.

QUANTITIES REPORTED
-------------------
All ratios are the corrected mass divided by the true Bonnor-Ebert mass M_BE;
the ideal value is 1.  "Scatter" is 10 raised to the standard deviation of
log10(ratio).  "Residual" is the root-mean-square of
log10(predicted C_M / true C_M).

USAGE
-----
    python3 estimator_compare.py --module mass_correction_pipeline_4d.py

Options:
    --cat, --rectab   catalog paths (defaults as in the pipeline)
    --module          path to the PATCHED mass_correction_pipeline_4d.py
    --nlocal          neighbourhood size for the Delaunay shortcut (default 1200)
    --skip-delaunay   omit the two Delaunay estimators (they dominate the runtime)
    --out             output CSV of per-node predictions
                      (default estimator_compare_results.csv)

Runtime: the two Delaunay estimators take a few minutes each; the local linear
and radial basis estimators together take well under a minute.
"""
import argparse, importlib.util, sys, time
import numpy as np
from scipy.interpolate import LinearNDInterpolator, RBFInterpolator
from scipy.spatial import cKDTree

PHI_REF = 1.84


def _load_module(path):
    spec = importlib.util.spec_from_file_location("mc_patched", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ----------------------------- estimators ---------------------------------
def predict_delaunay(X, y, tree, keep, i, nlocal):
    n = len(X)
    for N in (nlocal, 4 * nlocal, n):
        if N >= n:
            jj = np.where(keep)[0]
        else:
            _, cand = tree.query(X[i], k=min(N + 200, n))
            cand = cand[keep[cand]]
            jj = cand[:N]
        if len(jj) < X.shape[1] + 2:
            continue
        try:
            v = LinearNDInterpolator(X[jj], y[jj])(X[i][None, :])[0]
        except Exception:
            v = np.nan
        if np.isfinite(v):
            return float(v)
        if N >= n:
            break
    return np.nan


def predict_locallin(X, y, tree, keep, i, K):
    """Tricube-weighted local linear least squares on the K nearest samples."""
    d, cand = tree.query(X[i], k=min(K + 300, len(X)))
    m = keep[cand]
    cand, d = cand[m][:K], d[m][:K]
    if len(cand) < X.shape[1] + 3:
        return np.nan
    h = d[-1] if d[-1] > 0 else 1.0
    w = np.clip(1.0 - (d / h) ** 3, 0, None) ** 3
    w = np.maximum(w, 1e-6)
    A = np.column_stack([np.ones(len(cand)), X[cand] - X[i]])
    sw = np.sqrt(w)
    try:
        b, *_ = np.linalg.lstsq(A * sw[:, None], y[cand] * sw, rcond=None)
    except Exception:
        return np.nan
    return float(b[0])


def predict_rbf(X, y, tree, keep, i, K, smooth):
    _, cand = tree.query(X[i], k=min(K + 300, len(X)))
    cand = cand[keep[cand]][:K]
    if len(cand) < X.shape[1] + 3:
        return np.nan
    try:
        f = RBFInterpolator(X[cand], y[cand], kernel='thin_plate_spline',
                            smoothing=smooth, degree=1)
        return float(f(X[i][None, :])[0])
    except Exception:
        return np.nan


# ----------------------------- reporting ----------------------------------
def summarize(label, ratio, resid=None):
    r = np.asarray(ratio, float)
    m = np.isfinite(r) & (r > 0)
    r = r[m]
    lg = np.log10(r)
    if resid is None:
        tail = "residual        --"
    else:
        e = np.asarray(resid, float); e = e[np.isfinite(e)]
        tail = "residual %.4f dex (factor %.3f)" % (e.std(), 10 ** e.std())
    print("  %-22s n=%3d  median %.3f  scatter %.3f  within2x %4.0f%%  %s"
          % (label, len(r), np.median(r), 10 ** lg.std(),
             100 * np.mean((r > 0.5) & (r < 2.0)), tail))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cat', default='cats/bes_model_grid_final2_catalog')
    ap.add_argument('--rectab',
                    default='cats/bes_model_grid_final2_recovery_tables_v2')
    ap.add_argument('--module', default='mass_correction_pipeline_4d.py')
    ap.add_argument('--nlocal', type=int, default=1200)
    ap.add_argument('--stride', type=int, default=1)
    ap.add_argument('--skip-delaunay', action='store_true')
    ap.add_argument('--out', default='estimator_compare_results.csv')
    a = ap.parse_args()

    mc = _load_module(a.module)
    if not hasattr(mc, 'BANDR_DIAG'):
        sys.exit("the module at %s is NOT patched: run patch_band_R.py first"
                 % a.module)
    mc.VERBOSE = 1
    t0 = time.time()
    S = mc.build_recovery_samples(a.cat, a.rectab)
    print("\nband-footprint solver outcomes: %s" % mc.BANDR_DIAG)

    nk = np.array([str(tuple(n)) for n in S['node']])
    uniq, nid = np.unique(nk, return_inverse=True)
    X = np.column_stack([np.log10(S['sig']), np.log10(S['conc']),
                         np.log10(S['cfoot']), np.log10(S['rbe_pc'])])
    y = np.log10(S['cm'])
    good = np.all(np.isfinite(X), 1) & np.isfinite(y)
    X, y, nid = X[good], y[good], nid[good]
    MBE, Mrep, cmt = S['mbe'][good], S['mrec'][good], S['cm'][good]
    cfoot = S['cfoot'][good]

    tasks = []
    for u in np.unique(nid):
        idx = np.where(nid == u)[0]
        tasks.append((u, int(idx[np.argmin(np.abs(cfoot[idx] - PHI_REF))])))
    tasks = tasks[::a.stride]
    print("nodes to evaluate: %d   training samples: %d" % (len(tasks), len(y)))

    sd = X.std(0)
    spaces = {}                                    # name -> (scaled X, KD tree)
    spaces['raw'] = (X, cKDTree(X))
    for W in (0.5, 1.0, 2.0, 4.0):
        Xs = (X / sd) * np.array([1.0, 1.0, 1.0, W])
        spaces['std_w%.1f' % W] = (Xs, cKDTree(Xs))

    preds = {}

    def run(name, fn):
        t = time.time()
        out = np.full(len(tasks), np.nan)
        for k, (u, i) in enumerate(tasks):
            out[k] = fn(u, i)
        preds[name] = out
        print("   %-22s done in %5.1f s" % (name, time.time() - t), flush=True)

    if not a.skip_delaunay:
        Xr, Tr = spaces['raw']
        run('delaunay_raw',
            lambda u, i: predict_delaunay(Xr, y, Tr, nid != u, i, a.nlocal))
        Xs, Ts = spaces['std_w1.0']
        run('delaunay_std',
            lambda u, i: predict_delaunay(Xs, y, Ts, nid != u, i, a.nlocal))

    for W in (0.5, 1.0, 2.0, 4.0):
        Xs, Ts = spaces['std_w%.1f' % W]
        for K in (20, 40, 80, 160, 320, 640):
            run('locallin_k%d_w%.1f' % (K, W),
                lambda u, i, Xs=Xs, Ts=Ts, K=K: predict_locallin(Xs, y, Ts,
                                                                 nid != u, i, K))
    Xs, Ts = spaces['std_w1.0']
    for K in (60, 150):
        for smooth in (0.0, 0.1, 1.0):
            run('rbf_n%d_s%.1f' % (K, smooth),
                lambda u, i, K=K, s=smooth: predict_rbf(Xs, y, Ts, nid != u,
                                                        i, K, s))

    ir = np.array([i for (u, i) in tasks])
    truth = np.log10(cmt[ir])
    print("\n=== estimator comparison, node-level leave-one-out ===")
    print("All ratios are corrected mass divided by the true Bonnor-Ebert mass")
    print("M_BE; ideal value 1.  Scatter is 10**std(log10 ratio).  Residual is")
    print("the rms of log10(predicted C_M / true C_M).")
    summarize("uncorrected", Mrep[ir] / MBE[ir])
    order = sorted(preds, key=lambda n: np.nanstd(preds[n] - truth))
    for nm in order:
        p = preds[nm]
        summarize(nm, Mrep[ir] * 10 ** p / MBE[ir], p - truth)
        print("        %-16s finite predictions: %d of %d"
              % ('', np.isfinite(p).sum(), len(p)))

    # split-half honest check over the local linear family
    rng = np.random.default_rng(20260806)
    half = rng.random(len(ir)) < 0.5
    fam = [n for n in preds if n.startswith('locallin')]
    if fam:
        best = min(fam, key=lambda n: np.nanstd((preds[n] - truth)[half]))
        e = (preds[best] - truth)[~half]
        e = e[np.isfinite(e)]
        print("\nsplit-half check: hyperparameters chosen on a random half of the")
        print("nodes, evaluated on the other half (seed 20260806)")
        print("  chosen: %s   residual on the held-out half: %.4f dex (factor %.3f)"
              % (best, e.std(), 10 ** e.std()))

    with open(a.out, 'w') as fh:
        fh.write("# estimator comparison, node-level leave-one-out\n")
        fh.write("# M_BE (Msun), M_reported (Msun), C_M_true = M_BE/M_reported\n")
        fh.write("# remaining columns: predicted C_M from each estimator\n")
        names = list(preds)
        fh.write("node_key,M_BE,M_reported,C_M_true," + ",".join(names) + "\n")
        for k, (u, i) in enumerate(tasks):
            fh.write("%s,%.6g,%.6g,%.6g," % (uniq[u].replace(',', ' '),
                                             MBE[i], Mrep[i], cmt[i]))
            fh.write(",".join("%.6g" % 10 ** preds[n][k] for n in names) + "\n")
    print("\nwrote %s (%d rows)" % (a.out, len(tasks)))
    print("total elapsed: %.1f min" % ((time.time() - t0) / 60))


if __name__ == '__main__':
    main()
