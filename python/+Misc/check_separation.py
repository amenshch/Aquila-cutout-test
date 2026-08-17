#!/usr/bin/env python3
"""
check_separation.py -- measure the realised gaps between injected cores, to
decide whether the separation factor can be reduced.

    python3 check_separation.py inj_<tag>_truth.txt bes_model_grid_final2_catalog \
            [footprint_radii.txt] [--sep-factor=1.3]

Pass --sep-factor with the value the run actually used; it defaults to 1.5.

WHAT IT MEASURES
----------------
The placement loop enforces a centre-to-centre distance of

    d_ij >= SEP_FACTOR * (r_i + r_j),     r = sqrt(R_BE^2 + BEAM500^2)

with SEP_FACTOR = 1.5 and BEAM500 = 36.3 arcsec.  Whether that factor can be
lowered depends not on the typical gap but on the TIGHTEST pairs, because in a
sparse packing most realised distances sit far above the enforced minimum.

For every pair of injected cores the script reports

    slack_ij = d_ij / (r_i + r_j)

the realised separation in units of the sum of convolved radii.  A pair sitting
exactly on the constraint has slack = SEP_FACTOR.  The smallest slack over all
pairs is the largest factor the current placement would still satisfy, so it is
the value below which nothing is gained by lowering SEP_FACTOR, and the value
above which cores would have had to be rejected.

It also reports the gap in pixels between the 500 micron footprint edges, taking
the footprint radius as SEP_FACTOR * r, which is what the calibration against
those footprints implies.
"""
import sys
import numpy as np

BEAM500 = 36.3
PIX = 3.0
SEP_FACTOR = 1.5


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    truth_path, cat_path = sys.argv[1], sys.argv[2]
    # SEP_FACTOR used for the run being checked; overrides the module default
    global SEP_FACTOR
    args = [a for a in sys.argv[3:] if not a.startswith('--sep-factor=')]
    for a in sys.argv[3:]:
        if a.startswith('--sep-factor='):
            SEP_FACTOR = float(a.split('=', 1)[1])
    sys.argv = sys.argv[:3] + args

    rbe_of = {}
    for line in open(cat_path):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        f = line.split()
        rbe_of[(int(f[1]), int(f[2]), int(f[3]))] = float(f[8])

    x, y, r = [], [], []
    for line in open(truth_path):
        if line.lstrip().startswith('#') or not line.strip():
            continue
        f = line.split()
        x.append(float(f[2])); y.append(float(f[3])); r.append(float(f[9]))
    x, y, rbe = np.array(x), np.array(y), np.array(r)
    n = len(x)
    if n < 2:
        sys.exit('fewer than two cores in %s' % truth_path)
    rc = np.sqrt(rbe ** 2 + BEAM500 ** 2)

    d = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :]) * PIX
    rsum = rc[:, None] + rc[None, :]
    iu = np.triu_indices(n, 1)
    d, rsum = d[iu], rsum[iu]
    slack = d / rsum

    print('injected cores: %d   pairs: %d' % (n, len(d)))
    print('R_BE of injected cores: median %.1f", range %.1f-%.1f"'
          % (np.median(rbe), rbe.min(), rbe.max()))
    print()
    print('slack = centre distance / (r_i + r_j); the constraint is '
          'slack >= %.2f' % SEP_FACTOR)
    q = np.percentile(slack, [0, 0.1, 1, 5, 50])
    print('   minimum %.3f   0.1st pct %.3f   1st %.3f   5th %.3f   median %.3f'
          % tuple(q))
    print()
    tight = slack < SEP_FACTOR * 1.05
    print('pairs within 5%% of the constraint: %d of %d (%.2f%%)'
          % (tight.sum(), len(d), 100 * tight.mean()))
    print()
    # Solve for k, the ratio of the true 500 micron footprint radius to the
    # convolved radius r = sqrt(R_BE^2 + BEAM500^2).  Two cores just touch when
    # d_ij = k * (r_i + r_j), so for every pair k_ij = d_ij / (r_i + r_j) is an
    # UPPER limit on k if that pair does not overlap.  The minimum over pairs is
    # therefore the tightest constraint the placement provides, and equals the
    # smallest factor that would have kept every pair separated.
    print('implied footprint ratio k = footprint radius / convolved radius r:')
    print('   the placement guarantees k <= %.3f' % slack.min())
    print('   (a pair overlaps once SEP_FACTOR falls below k)')
    print()
    if len(sys.argv) > 3:
        # optional: measured footprint radii, one per core, same order as truth
        f = np.array([float(v) for v in open(sys.argv[3]).read().split()])
        if len(f) == n:
            k = f / rc
            print('measured footprint radii supplied: k = %.3f +/- %.3f '
                  '(median %.3f, 95th pct %.3f)'
                  % (k.mean(), k.std(), np.median(k), np.percentile(k, 95)))
            fs = f[:, None] + f[None, :]
            clear = (d - fs[iu]) / PIX
            print('true clearance between footprint edges, pixels: '
                  'minimum %.1f, 1st pct %.1f, median %.1f'
                  % tuple(np.percentile(clear, [0, 1, 50])))
            print()
            print('recommended SEP_FACTOR: %.2f  '
                  '(95th percentile of k, plus 0.05 margin)'
                  % (np.percentile(k, 95) + 0.05))
            return
    print('To pin k down, supply a third argument: a file of measured 500 micron')
    print('footprint radii in arcsec, one per injected core, in the order of the')
    print('truth table.  The script then reports k directly and recommends a')
    print('SEP_FACTOR.  Without it, only the upper limit above is available.')


if __name__ == '__main__':
    main()
