#!/usr/bin/env python
"""
run_validation.py -- run every validation test of the mass-correction method
and report the scores in the form used by the paper's validation table.

WHY THIS SCRIPT EXISTS
======================
Earlier validation drivers hard-coded their own interpolator instead of using
the one the pipeline defines.  When the pipeline's estimator was replaced, the
drivers went on measuring the old one, and the published numbers could no
longer be reproduced by any script in the repository.  This script never builds
an estimator of its own: every correction comes from
mass_correction_pipeline_6.correct_cloud or from the pipeline's own corrector
objects, so it cannot drift away from the method it is validating.

THE DOMAIN QUESTION, AND WHY EVERY TEST IS SCORED TWICE
=======================================================
The correctors are built by interpolating a grid of radiative-transfer models.
A real source whose observables fall outside the region those models span
cannot be corrected by interpolation; the answer is an extrapolation.

Until pipeline version 6 the corrector was a Delaunay interpolator, which
returned not-a-number outside the convex hull of the training samples.  The
validation drivers discarded those sources, so the published scores were
implicitly restricted to sources inside the domain -- but nothing said so, and
no flag recorded it.  The estimator was then replaced by a local linear
regression, which extrapolates and always returns a value.  That silently
removed the filter: sources far outside the training set began to be scored,
and they are the ones the extraction handles worst.

Version 6 restores the distinction explicitly.  correct_cloud now reports
IN_HULL_3D, IN_HULL_4D and IN_DOMAIN per source, from a convex-hull test on
the training samples.  This script therefore reports every test twice:

    all matched sources     every source matched to a known truth
    inside the domain       IN_DOMAIN = 1, i.e. both corrections are
                            interpolations rather than extrapolations

The second is the honest comparison for judging the method; the first shows
what happens when it is applied outside its domain, which is also worth
knowing.  Quoting either alone, without saying which, is what caused the
earlier confusion.

THE FOUR TESTS
==============
loo           Node-level leave-one-out on the model grid.  Every sample of one
              node is removed, both correctors are rebuilt from the rest, and
              the withheld node's representative sample is predicted.  Tests
              interpolation between models, with no observational chain.
subdivision   The independent subdivision grid, whose nodes lie between the
              production nodes.  Tests interpolation at points the training
              set does not contain.
injection     Models injected into real Herschel maps and re-extracted with
              getsf, then corrected.  Tests the full observational chain.
bgscan        Seven fixed models placed at many different real background
              positions.  Holds the model constant and varies the environment.
              By design it cannot exercise the degeneracy the fourth
              observable exists to resolve -- different models sharing the
              same three observables -- so it is reported separately and its
              scores should not be pooled with the others.

SCORES REPORTED
===============
For each test and each of the three treatments (uncorrected reported mass,
three-observable correction, four-observable correction), over the ratio
r = M_obtained / M_true:

    n                number of sources scored
    median r         ideal 1; below 1 means mass is still under-estimated
    scatter          standard deviation of log10(r), as a multiplicative
                     factor; ideal 1
    within 2x        percentage with 0.5 < r < 2
    median |log10 r| typical absolute error in dex; the most robust single
                     number, being insensitive to the tails
    improved/worse   per-source comparison of four observables against three

USAGE
=====
    python run_validation.py all
    python run_validation.py loo --rectab bes_model_grid_final2_recovery_tables_v3

Required in the working directory (or give --cats-dir / --tests-dir):
    mass_correction_pipeline_6.py, getsf_columns.py
    cats/bes_model_grid_final2_catalog
    cats/bes_model_grid_final2_recovery_tables_v3
    cats/bes_model_grid_subdivision2_catalog
    cats/bes_model_grid_subdivision2_catalog_recovery_tables
    tests/  (injection and background-scan catalogs and truth tables)

OUTPUT
======
    validation_<test>.txt   one row per scored source, with a header giving
                            every column's meaning and units, the provenance
                            of the run, and the scores
    validation_summary.txt  all tests together, in the layout of the paper's
                            validation table
"""
import argparse
import datetime
import os
import platform
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '.')
import mass_correction_pipeline_6 as mc                      # noqa: E402

DIST_AQUILA = 260.0
REP_PHI = 1.84            # grid-typical footprint factor; picks one
                          # representative sample per node
MATCH_PIX = 10.0          # positional match radius, pixels
CSD_BLOCKS = [('01', 'x0.25'), ('02', 'x0.5'), ('03', ''),
              ('04', 'x2'), ('05', 'x4')]
BGSCAN_MODELS = ['bgscan_i1j02k05', 'bgscan_i1j07k06', 'bgscan_i1j08k08',
                 'bgscan_i1j12k08', 'bgscan_i1j13k06', 'bgscan_i2j07k05',
                 'bgscan_i3j15k07']


# ----------------------------------------------------------------------
# scoring
# ----------------------------------------------------------------------
def score(mtrue, mobt):
    """Scores of one treatment: see SCORES REPORTED in the module docstring."""
    r = np.asarray(mobt, float) / np.asarray(mtrue, float)
    g = np.isfinite(r) & (r > 0)
    if g.sum() == 0:
        return dict(n=0, median=np.nan, scatter=np.nan, within2=np.nan,
                    medabs=np.nan)
    lr = np.log10(r[g])
    return dict(n=int(g.sum()), median=float(np.median(r[g])),
                scatter=float(10 ** lr.std()),
                within2=float(100 * np.mean((r[g] > 0.5) & (r[g] < 2))),
                medabs=float(np.median(np.abs(lr))))


def score_block(rows, label):
    """rows: list of (M_true, M_uncorrected, M_3axis, M_4axis, in_domain)."""
    a = np.asarray([r[:4] for r in rows], float)
    dom = np.asarray([r[4] for r in rows], bool)
    out = {}
    for subset, m in (('all', np.ones(len(a), bool)), ('domain', dom)):
        if m.sum() == 0:
            out[subset] = None
            continue
        s = {}
        for nm, j in (('uncorrected', 1), ('3-axis', 2), ('4-axis', 3)):
            s[nm] = score(a[m, 0], a[m, j])
        d3 = np.abs(np.log10(a[m, 2] / a[m, 0]))
        d4 = np.abs(np.log10(a[m, 3] / a[m, 0]))
        ok = np.isfinite(d3) & np.isfinite(d4)
        s['paired'] = (int(np.sum(d4[ok] < d3[ok])), int(np.sum(d4[ok] > d3[ok])),
                       int(ok.sum()))
        out[subset] = s
    out['label'] = label
    return out


def print_block(res):
    print('\n=== %s ===' % res['label'])
    for subset, title in (('all', 'all matched sources'),
                          ('domain', 'inside the correction domain '
                                     '(IN_DOMAIN = 1)')):
        s = res[subset]
        if s is None:
            print('  %s: no sources' % title)
            continue
        print('  %s' % title)
        print('    %-14s %5s %9s %9s %10s %14s'
              % ('treatment', 'n', 'median', 'scatter', 'within 2x',
                 'median |log10 r|'))
        for nm in ('uncorrected', '3-axis', '4-axis'):
            v = s[nm]
            print('    %-14s %5d %9.3f %9.2f %9.0f%% %14.3f'
                  % (nm, v['n'], v['median'], v['scatter'], v['within2'],
                     v['medabs']))
        b, w, n = s['paired']
        print('    four observables better for %d, worse for %d, of %d'
              % (b, w, n))


# ----------------------------------------------------------------------
# output files
# ----------------------------------------------------------------------
def write_table(basename, rows, res, description, provenance):
    path = '%s.txt' % basename
    a = np.asarray([r[:4] for r in rows], float)
    dom = np.asarray([r[4] for r in rows], bool)
    with open(path, 'w') as f:
        f.write('# %s\n#\n' % res['label'])
        for line in description.strip().split('\n'):
            f.write('# %s\n' % line.strip())
        f.write('#\n# PROVENANCE\n')
        for line in provenance:
            f.write('#   %s\n' % line)
        f.write('#\n# COLUMNS\n')
        f.write('#   1 M_true       true mass of the source            '
                '(M_sun)\n')
        f.write('#   2 M_uncorr     reported SED mass, uncorrected     '
                '(M_sun)\n')
        f.write('#   3 M_3axis      mass after the three-observable    '
                '(M_sun)\n'
                '#                  correction\n')
        f.write('#   4 M_4axis      mass after the four-observable     '
                '(M_sun)\n'
                '#                  correction\n')
        f.write('#   5 r_uncorr     M_uncorr / M_true                  (-)\n')
        f.write('#   6 r_3axis      M_3axis / M_true                   (-)\n')
        f.write('#   7 r_4axis      M_4axis / M_true                   (-)\n')
        f.write('#   8 in_domain    1 if both corrections are interpolations\n'
                '#                  rather than extrapolations; see the\n'
                '#                  header of run_validation.py\n')
        f.write('#\n# SCORES\n')
        for subset, title in (('all', 'all matched sources'),
                              ('domain', 'inside the correction domain')):
            s = res[subset]
            if s is None:
                continue
            f.write('#   %s\n' % title)
            f.write('#     %-14s %5s %9s %9s %10s %16s\n'
                    % ('treatment', 'n', 'median', 'scatter', 'within_2x',
                       'median_abs_log10'))
            for nm in ('uncorrected', '3-axis', '4-axis'):
                v = s[nm]
                f.write('#     %-14s %5d %9.3f %9.2f %9.0f%% %16.3f\n'
                        % (nm, v['n'], v['median'], v['scatter'],
                           v['within2'], v['medabs']))
            b, w, n = s['paired']
            f.write('#     four observables better for %d, worse for %d, '
                    'of %d\n' % (b, w, n))
        f.write('#\n# %10s %12s %12s %12s %10s %10s %10s %10s\n'
                % ('M_true', 'M_uncorr', 'M_3axis', 'M_4axis',
                   'r_uncorr', 'r_3axis', 'r_4axis', 'in_domain'))
        for i in range(len(a)):
            f.write('  %10.5g %12.5e %12.5e %12.5e %10.4f %10.4f %10.4f '
                    '%10d\n'
                    % (a[i, 0], a[i, 1], a[i, 2], a[i, 3],
                       a[i, 1] / a[i, 0], a[i, 2] / a[i, 0],
                       a[i, 3] / a[i, 0], int(dom[i])))
    print('  wrote %s' % path)


# ----------------------------------------------------------------------
# truth-table reading and matching
# ----------------------------------------------------------------------
def load_truth(path):
    """Injection truth table: x_pix, y_pix, M_BE.

    Column 9 (one-based) is M_BE, the true Bonnor-Ebert mass of the injected
    model, in solar masses.  The last column, scale, records SCALE_FACTOR,
    which multiplied the FIELD images and the column-density map before
    injection so that models built for dense clouds had somewhere to sit.  It
    did not scale the models, so M_BE needs no conversion.
    """
    rows = []
    for line in open(path):
        if not line.strip() or line.startswith('#'):
            continue
        f = line.split()
        rows.append((float(f[2]), float(f[3]), float(f[8])))
    return np.array(rows)


def match_injection(ext, truth_path, label, dedupe=True):
    """Correct one injected field and match its sources to the truth table.

    Matching is one-to-one when dedupe is set: several extracted sources can
    fall within the match radius of the same injected core, and counting each
    of them would weight that core several times over.  The closest is kept.
    """
    truth = load_truth(truth_path)
    table, _ = mc.correct_cloud(ext, DIST_AQUILA)
    main_cat = ext if isinstance(ext, str) else ext[0]
    pix = {}
    for line in open(main_cat):
        s = line.strip()
        if not s or s[0] in '#!':
            continue
        f = s.split()
        try:
            pix[int(float(f[0]))] = (float(f[4]), float(f[5]))
        except (ValueError, IndexError):
            pass

    best = {}
    for i, no in enumerate(table['NO']):
        no = int(no)
        if no not in pix:
            continue
        xp, yp = pix[no]
        d = np.hypot(truth[:, 0] - xp, truth[:, 1] - yp)
        j = int(np.argmin(d))
        if d[j] >= MATCH_PIX:
            continue
        if not np.isfinite(table['M_corr_4D'][i]):
            continue
        rec = (d[j], truth[j, 2], table['M_init'][i], table['M_corr'][i],
               table['M_corr_4D'][i], bool(table['IN_DOMAIN'][i]))
        if dedupe:
            if j in best and best[j][0] <= d[j]:
                continue
            best[j] = rec
        else:
            best[len(best) + 1000000] = rec
    rows = [(r[1], r[2], r[3], r[4], r[5]) for r in best.values()]
    n_dom = sum(1 for r in rows if r[4])
    mc.vlog(0, '  %-18s %3d truth rows, %3d matched, %3d inside the domain'
            % (label + ':', len(truth), len(rows), n_dom))
    return rows


# ----------------------------------------------------------------------
# the four tests
# ----------------------------------------------------------------------
def test_injection(tests_dir, cats_dir, dedupe=True):
    rows = []
    for n, x in CSD_BLOCKS:
        ext = ('%s/Aquila.s.sources.ok.cat=Aquila.s.sources.ok.add.cat='
               'thin.Aquila.s.sources.ok.00.cat_cSD_%s' % (tests_dir, n))
        tr = ('%s/inj_all24_s1111_%s_truth.txt_cSD_%s' % (tests_dir, x, n)) \
            if x else ('%s/inj_all24_s1111_truth.txt_cSD_%s' % (tests_dir, n))
        rows += match_injection(ext, tr, 'cSD_%s' % n, dedupe)
    return rows, ('Injection into real Herschel maps (5 cSD blocks)',
                  """Bonnor-Ebert models injected into the real Aquila sub-field and
                  re-extracted with getsf, then corrected.  The five blocks differ
                  in SCALE_FACTOR, which multiplied the field images and the
                  column-density map by 0.25, 0.5, 1, 2 and 4 so that models built
                  for denser clouds had positions available.  The models themselves
                  were not scaled.  Blocks 04 and 05 contribute almost nothing,
                  because at 2x and 4x the extraction itself breaks down.
                  This test exercises the full observational chain.""")


def test_bgscan(tests_dir, cats_dir, dedupe=True):
    rows = []
    for d in BGSCAN_MODELS:
        ext = ['%s/%s/Aquila.s.sources.ok.cat' % (tests_dir, d),
               '%s/%s/Aquila.s.sources.ok.add.cat' % (tests_dir, d)]
        tr = '%s/%s/inj_%s_truth.txt' % (tests_dir, d, d)
        rows += match_injection(ext, tr, d, dedupe)
    return rows, ('Background scan (7 fixed models, varying real background)',
                  """Seven fixed models placed at many different real background
                  positions each.  Holding the model constant and varying only the
                  environment, this test cannot exercise the degeneracy the fourth
                  observable exists to resolve, namely different models sharing the
                  same three observables.  Its scores should therefore be reported
                  separately and not pooled with the other three tests.""")


def _node_test(cat_path, rectab_path, label, description):
    """Leave-one-out over nodes of a model grid, using the pipeline's own
    estimator so that the test follows the method rather than a copy of it."""
    lg = np.log10
    S = mc.build_recovery_samples(cat_path, rectab_path)
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    node = np.asarray(S['node'])[ok]
    theta = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
    sig, conc, phi = S['sig'][ok], S['conc'][ok], S['cfoot'][ok]
    cm, mrec, mbe = S['cm'][ok], S['mrec'][ok], S['mbe'][ok]
    X3 = np.column_stack([lg(sig), lg(conc), lg(phi)])
    X4 = np.column_stack([lg(sig), lg(conc), lg(phi), lg(theta)])
    y = lg(cm)
    key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]
    uniq = np.unique(key)
    mc.vlog(0, '  %d samples from %d nodes' % (len(key), len(uniq)))

    rows = []
    for k in uniq:
        out_m = key == k
        in_m = ~out_m
        e3 = mc.LocalLinear(X3[in_m], y[in_m])
        e4 = mc.LocalLinear(X4[in_m], y[in_m])
        h3 = mc._HullTest(X3[in_m])
        h4 = mc._HullTest(X4[in_m])
        idx = np.where(out_m)[0]
        j = idx[np.argmin(np.abs(phi[idx] - REP_PHI))]
        p3 = float(e3(X3[j:j + 1])[0])
        p4 = float(e4(X4[j:j + 1])[0])
        dom = bool(h3(X3[j:j + 1])[0] and h4(X4[j:j + 1])[0])
        rows.append((mbe[j], mrec[j], mrec[j] * 10 ** p3,
                     mrec[j] * 10 ** p4, dom))
    return rows, (label, description)


def test_loo(tests_dir, cats_dir, rectab=None):
    return _node_test(
        '%s/bes_model_grid_final2_catalog' % cats_dir,
        rectab or '%s/bes_model_grid_final2_recovery_tables_v3' % cats_dir,
        'Node-level leave-one-out on the production grid',
        """Every sample of one node is removed from the training set, both
        correctors are rebuilt from the remaining nodes, and the correction is
        evaluated on that node's representative sample, the one whose footprint
        factor is closest to the grid median of 1.84.  No observational chain is
        involved, so this measures interpolation between models alone.""")


def test_subdivision(tests_dir, cats_dir):
    return _node_test(
        '%s/bes_model_grid_subdivision2_catalog' % cats_dir,
        '%s/bes_model_grid_subdivision2_catalog_recovery_tables' % cats_dir,
        'Subdivision grid',
        """The subdivision grid's nodes lie between those of the production
        grid, so predicting them tests interpolation at points the training set
        does not contain, independently of the leave-one-out protocol.""")


TESTS = dict(injection=test_injection, bgscan=test_bgscan,
             loo=test_loo, subdivision=test_subdivision)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('test', choices=sorted(TESTS) + ['all'])
    ap.add_argument('--cats-dir', default='cats')
    ap.add_argument('--tests-dir', default='tests')
    ap.add_argument('--rectab', default=None)
    ap.add_argument('--no-dedupe', action='store_true',
                    help='count every extracted source that falls within the '
                         'match radius, even when several claim the same '
                         'injected core (not recommended)')
    ap.add_argument('--verbose', type=int, default=0)
    a = ap.parse_args()
    mc.VERBOSE = a.verbose

    provenance = [
        'run_validation.py, %s' % datetime.datetime.now().isoformat(' ',
                                                                    'seconds'),
        'pipeline: %s' % os.path.abspath(mc.__file__),
        'estimator: LocalLinear, LL_NEIGHBOURS = %d, LL_BOUND = %s'
        % (mc.LL_NEIGHBOURS, mc.LL_BOUND),
        'grid catalog: %s' % (a.rectab or mc.CAT67),
        'python %s on %s' % (platform.python_version(), platform.platform()),
    ]

    which = sorted(TESTS) if a.test == 'all' else [a.test]
    summary = []
    for name in which:
        fn = TESTS[name]
        kw = {}
        if name == 'loo' and a.rectab:
            kw['rectab'] = a.rectab
        if name in ('injection', 'bgscan'):
            kw['dedupe'] = not a.no_dedupe
        rows, (label, desc) = fn(a.tests_dir, a.cats_dir, **kw)
        if not rows:
            print('\n=== %s: no sources scored ===' % label)
            continue
        res = score_block(rows, label)
        print_block(res)
        write_table('validation_%s' % name, rows, res, desc, provenance)
        summary.append((name, res))

    if summary:
        path = 'validation_summary.txt'
        with open(path, 'w') as f:
            f.write('# Validation of the mass-correction method: all tests\n#\n')
            for line in provenance:
                f.write('# %s\n' % line)
            f.write('#\n# Ratio r = M_obtained / M_true.  Scatter is the '
                    'standard deviation of\n# log10(r) as a multiplicative '
                    'factor.  Each test is scored twice: over\n# every matched '
                    'source, and over those inside the correction domain,\n'
                    '# meaning both corrections are interpolations rather than '
                    'extrapolations.\n#\n')
            for subset, title in (('all', 'ALL MATCHED SOURCES'),
                                  ('domain', 'INSIDE THE CORRECTION DOMAIN')):
                f.write('# %s\n' % title)
                f.write('# %-34s %6s %8s %9s %10s %16s\n'
                        % ('test / treatment', 'n', 'median', 'scatter',
                           'within_2x', 'median_abs_log10'))
                for name, res in summary:
                    s = res[subset]
                    if s is None:
                        continue
                    for nm in ('uncorrected', '3-axis', '4-axis'):
                        v = s[nm]
                        f.write('  %-34s %6d %8.3f %9.2f %9.0f%% %16.3f\n'
                                % ('%s / %s' % (name, nm), v['n'], v['median'],
                                   v['scatter'], v['within2'], v['medabs']))
                f.write('#\n')
        print('\nwrote %s' % path)


if __name__ == '__main__':
    main()
