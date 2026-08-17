#!/usr/bin/env python3
"""
loo_4axis_full_v2.py -- full node-level leave-one-out comparison of the
three-observable corrector against the four-observable corrector, using the
estimator and the axis definitions of the current pipeline.

WHAT CHANGED SINCE loo_4axis_full.py (chat 10)
----------------------------------------------
Three things, of which the first is the one that matters.

1. THE ESTIMATOR.  The old script built both correctors with
   scipy.interpolate.LinearNDInterpolator, a Delaunay linear interpolator.
   The pipeline no longer uses it: mass_correction_pipeline_5.cm_interp and
   cm_size_interp both build a LocalLinear, a tricube-weighted local linear
   regression that scales each axis by its own standard deviation and
   extrapolates outside the convex hull instead of returning not-a-number.
   A leave-one-out test built on the old interpolator therefore measures an
   estimator the pipeline has stopped using, and its numbers are not
   comparable with the paper's table.  This script uses mc.LocalLinear.

   The consequence is not cosmetic: with the Delaunay interpolator a withheld
   node lying outside the convex hull of the remaining nodes returns
   not-a-number and drops out of the test, so the old script scored only the
   nodes that were easy to reach.  LocalLinear returns a value for every node,
   bounded to the range of its own neighbours when the local plane
   extrapolates, and reports how often that bound was applied.  The count of
   bounded predictions is printed below, because a high value means the test
   is being carried by extrapolation.

2. THE FOURTH AXIS IS ANGULAR.  The old script used the physical truncation
   radius in parsecs.  The pipeline uses the angle it subtends, in arcsec,
   because the correction factor is a ratio of two masses and the distance
   cancels exactly.  On the model grid the distinction cannot show up: every
   model was rendered at 260 pc, so the angular and the physical axis differ
   by a constant in the logarithm, and translating one coordinate of the
   training set changes neither the neighbour ranking nor the fitted plane.
   The two therefore give identical leave-one-out numbers, and this script
   uses the angular form only so that the axis matches the pipeline's.  The
   difference is real when a real source is queried, not here.

3. THE RECOVERY TABLE.  Default is now the v3 table, which is the v2 table
   without the 70 and 100 micron bands.  Use --rectab to select the other.

USAGE
    Put in the same directory as this script:
        mass_correction_pipeline_5.py
        getsf_columns.py
        bes_model_grid_final2_catalog
        bes_model_grid_final2_recovery_tables_v3
    then:
        python3 loo_4axis_full_v2.py

    Options:
        --rectab FILE      recovery table to use
        --estimator {locallinear,delaunay}
                           locallinear (default) matches the pipeline;
                           delaunay reproduces the old script, for comparison
        --out FILE         output CSV (default loo4d_full_results_v2.csv)

    The run checkpoints every 20 nodes, so it can be interrupted and resumed
    by re-running the same command.  Expect roughly the same runtime as the
    old script; LocalLinear is faster to build than a 4-D Delaunay
    triangulation but is evaluated per query.

OUTPUT
    One row per node: node id, true M_BE, uncorrected reported mass, and the
    corrected mass from each corrector, plus flags recording whether each
    prediction required the local-range bound.  A console summary gives, for
    uncorrected and for both correctors, the median ratio of corrected to true
    mass, the scatter as a multiplicative factor, and the fraction within a
    factor of two.
"""
import argparse
import csv
import os
import pickle
import sys
import time

sys.path.insert(0, '.')
import numpy as np                                     # noqa: E402
import mass_correction_pipeline_5 as mc                # noqa: E402

lg = np.log10
CATALOG = 'bes_model_grid_final2_catalog'
RECTAB_DEFAULT = 'bes_model_grid_final2_recovery_tables_v3'
CKPT = 'loo4d_checkpoint_v2.pkl'
REP_RFOOT = 1.84          # grid-typical footprint factor, used to pick one
                          # representative sample per node
CHECKPOINT_EVERY = 20
GRID_DISTANCE_PC = mc.GRID_DISTANCE_PC


def check_inputs(rectab):
    missing = [f for f in (CATALOG, rectab, 'mass_correction_pipeline_5.py')
               if not os.path.exists(f)]
    if missing:
        print('ERROR: required file(s) not found in the current directory:')
        for f in missing:
            print('  - %s' % f)
        print('\nCopy them into this directory (no cats/ or python/ '
              'subdirectory needed) and re-run.')
        sys.exit(1)


def build_estimator(kind, X, y):
    if kind == 'locallinear':
        return mc.LocalLinear(X, y)
    from scipy.interpolate import LinearNDInterpolator
    return LinearNDInterpolator(X, y)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rectab', default=RECTAB_DEFAULT)
    ap.add_argument('--estimator', default='locallinear',
                    choices=('locallinear', 'delaunay'))
    ap.add_argument('--out', default='loo4d_full_results_v2.csv')
    a = ap.parse_args()
    check_inputs(a.rectab)

    print('Estimator: %s%s' % (a.estimator,
                               '  (matches the pipeline)'
                               if a.estimator == 'locallinear'
                               else '  (old script; for comparison only)'))
    print('Recovery table: %s' % a.rectab)
    print('Building training samples from the production grid ...')
    S = mc.build_recovery_samples(CATALOG, a.rectab)

    # The pipeline's own axis definition: the angle subtended by the model's
    # truncation radius, in arcsec.  See note 2 in the module docstring.
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    node = np.asarray(S['node'])[ok]
    theta = S['rbe_pc'][ok] * 206264.806 / GRID_DISTANCE_PC        # arcsec
    sig, conc, cfoot = S['sig'][ok], S['conc'][ok], S['cfoot'][ok]
    cm, mrec, mbe = S['cm'][ok], S['mrec'][ok], S['mbe'][ok]

    X3 = np.column_stack([lg(sig), lg(conc), lg(cfoot)])
    X4 = np.column_stack([lg(sig), lg(conc), lg(cfoot), lg(theta)])
    y = lg(cm)
    key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]
    uniq = np.unique(key)
    print('n samples=%d  n unique nodes=%d' % (len(key), len(uniq)))

    if os.path.exists(CKPT):
        ck = pickle.load(open(CKPT, 'rb'))
        results, done = ck['results'], ck['done']
        print('resuming from checkpoint: %d/%d nodes already done'
              % (done, len(uniq)))
    else:
        results, done = [], 0

    done0, t0 = done, time.time()
    for i in range(done, len(uniq)):
        kk = uniq[i]
        out_mask = key == kk
        in_mask = ~out_mask

        e3 = build_estimator(a.estimator, X3[in_mask], y[in_mask])
        e4 = build_estimator(a.estimator, X4[in_mask], y[in_mask])

        idxs = np.where(out_mask)[0]
        j = idxs[np.argmin(np.abs(cfoot[idxs] - REP_RFOOT))]

        p3, p4 = e3(X3[j:j + 1]), e4(X4[j:j + 1])
        b3 = bool(getattr(e3, '_flag', np.zeros(1, bool))[0]) \
            if a.estimator == 'locallinear' else False
        b4 = bool(getattr(e4, '_flag', np.zeros(1, bool))[0]) \
            if a.estimator == 'locallinear' else False
        p3 = float(p3[0]) if np.isfinite(p3[0]) else np.nan
        p4 = float(p4[0]) if np.isfinite(p4[0]) else np.nan

        results.append(dict(node=int(kk), mrec=float(mrec[j]),
                            mbe=float(mbe[j]), theta=float(theta[j]),
                            log_cm_pred_3axis=p3, log_cm_pred_4axis=p4,
                            bounded3=b3, bounded4=b4))
        done = i + 1
        if done % CHECKPOINT_EVERY == 0 or done == len(uniq):
            el = time.time() - t0
            n = max(done - done0, 1)
            print('  %4d/%4d nodes done (%.0fs elapsed, ~%.0fs left)'
                  % (done, len(uniq), el, (len(uniq) - done) * el / n))
            pickle.dump(dict(results=results, done=done), open(CKPT, 'wb'))

    print('\nAll %d nodes done. Writing %s ...' % (len(uniq), a.out))
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['node_key', 'M_BE_Msun', 'M_reported_Msun',
                    'M_corrected_3axis_Msun', 'M_corrected_4axis_Msun',
                    'theta_BE_arcsec', 'bounded_3axis', 'bounded_4axis'])
        for r in results:
            f3 = (r['mrec'] * 10 ** r['log_cm_pred_3axis']
                  if np.isfinite(r['log_cm_pred_3axis']) else '')
            f4 = (r['mrec'] * 10 ** r['log_cm_pred_4axis']
                  if np.isfinite(r['log_cm_pred_4axis']) else '')
            w.writerow([r['node'], r['mbe'], r['mrec'], f3, f4,
                        '%.2f' % r['theta'], int(r['bounded3']),
                        int(r['bounded4'])])

    mbe_a = np.array([r['mbe'] for r in results])
    mrec_a = np.array([r['mrec'] for r in results])
    p3a = np.array([r['log_cm_pred_3axis'] for r in results])
    p4a = np.array([r['log_cm_pred_4axis'] for r in results])
    both = np.isfinite(p3a) & np.isfinite(p4a)
    n3 = int(np.isfinite(p3a).sum()); n4 = int(np.isfinite(p4a).sum())
    print('\npredictions returned: 3-axis %d, 4-axis %d, of %d nodes'
          % (n3, n4, len(results)))
    if a.estimator == 'locallinear':
        print('predictions requiring the local-range bound: 3-axis %d, '
              '4-axis %d' % (sum(r['bounded3'] for r in results),
                             sum(r['bounded4'] for r in results)))

    r_unc = (mrec_a / mbe_a)[both]
    r3 = (mrec_a[both] * 10 ** p3a[both]) / mbe_a[both]
    r4 = (mrec_a[both] * 10 ** p4a[both]) / mbe_a[both]

    print('\n=== Ratio of the mass obtained to the true mass M_BE, over the '
          '%d nodes\n    scored by both correctors ===' % both.sum())
    print('  %-40s %8s %10s %10s' % ('', 'median', 'scatter', 'within 2x'))
    for label, r in (('reported mass, uncorrected', r_unc),
                     ('corrected, three observables', r3),
                     ('corrected, four observables', r4)):
        lr = np.log10(r)
        print('  %-40s %8.3f %9.2f%s %9.0f%%'
              % (label, np.median(r), 10 ** lr.std(), '', 
                 100 * np.mean((r > 0.5) & (r < 2))))
    print('  (scatter is the standard deviation of log10(ratio) expressed as '
          'a multiplicative factor)')

    d3, d4 = np.abs(np.log10(r3)), np.abs(np.log10(r4))
    print('\n  median |log10 ratio|: three observables %.3f, four %.3f'
          % (np.median(d3), np.median(d4)))
    print('  paired per node: four observables better for %d, worse for %d, '
          'of %d' % (int(np.sum(d4 < d3)), int(np.sum(d4 > d3)), int(both.sum())))

    if os.path.exists(CKPT):
        os.remove(CKPT)
    print('\nDone. Per-node results in %s' % a.out)


if __name__ == '__main__':
    main()
