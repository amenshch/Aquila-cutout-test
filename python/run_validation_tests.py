#!/usr/bin/env python3
"""
run_validation_tests.py -- run all four validation tests from the paper
(injection/cSD, subdivision, node-LOO, background-scan) using ONLY
mass_correction_pipeline_4d.py's own public functions.

WHY THIS SCRIPT EXISTS
-----------------------
Earlier, ad hoc versions of these tests re-implemented file-reading and
position-matching logic instead of reusing load_getsf()/correct_cloud()
directly, and that re-implementation had two real bugs: (1) merging
column names from only the first of several catalogue files, silently
hiding columns (like PEAK_BGF) that only exist in add.cat; (2) reading
pixel position from the wrong column offset when a second file was
present, because add.cat's own column layout differs from the main
catalogue's. Both bugs were in the *test-runner* re-implementation, not
in the pipeline itself -- load_getsf() already merges files correctly.
This script fixes that by never re-implementing file reading: every test
below calls correct_cloud()/load_getsf() exactly as the pipeline's own
'correct' CLI command does, and only adds the truth-matching logic on
top.

USAGE
-----
    cd <directory with mass_correction_pipeline_4d.py, cats/, tests/>
    python3 run_validation_tests.py injection      # cSD blocks (5), all scales
    python3 run_validation_tests.py subdivision    # needs a subdivision2 catalog
    python3 run_validation_tests.py loo            # node-level leave-one-out
    python3 run_validation_tests.py bgscan         # 7 fixed models, varying bg
    python3 run_validation_tests.py all            # all four, one after another
    python3 run_validation_tests.py all --verbose 1   # with pipeline progress messages

Each test prints a summary table (uncorrected / 3-axis / 4-axis: median,
scatter as a multiplicative factor, fraction within a factor of two,
median |log10 ratio|, and the paired per-source win count) and saves its
raw per-source results to a .pkl file for further analysis or plotting.
"""
import argparse, pickle, sys, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '.')
import mass_correction_pipeline_4d as mc

DIST_AQUILA = 260.0


# =========================================================================
# shared helpers
# =========================================================================
def _loadtruth(path):
    """Truth table columns: id, model, x_pix, y_pix, ... M_BE ... (col 8)."""
    rows = []
    for l in open(path):
        if not l.strip() or l.startswith('#'):
            continue
        f = l.split()
        rows.append((float(f[2]), float(f[3]), float(f[8])))
    return np.array(rows)


def _write_table(basename, mt, mi, m3, m4):
    """Space-separated table: ref_mass, uncorrected, corrected_3axis,
    corrected_4axis, and their ratios -- one row per matched source."""
    path = basename + '.txt'
    with open(path, 'w') as fh:
        fh.write("# %s -- per-source validation results\n" % basename)
        fh.write("# M_ref: true/reference mass (Msun)\n")
        fh.write("# M_unc: uncorrected (reported) mass\n")
        fh.write("# M_3ax, M_4ax: corrected mass, 3-axis and 4-axis correctors\n")
        fh.write("# ratio_unc, ratio_3ax, ratio_4ax: M / M_ref\n")
        fh.write("#%9s %12s %12s %12s %10s %10s %10s\n" % (
            "M_ref", "M_unc", "M_3ax", "M_4ax", "r_unc", "r_3ax", "r_4ax"))
        for a, b, c, d in zip(mt, mi, m3, m4):
            fh.write("%10.5f %12.5e %12.5e %12.5e %10.4f %10.4f %10.4f\n"
                      % (a, b, c, d, b / a, c / a, d / a))
    print("  wrote table: %s" % path)


def _plot_validation(basename, title, mt, mi, m3, m4):
    """Two-panel PDF matching the paper's fig_validation.pdf style:
    (a) corrected ratio vs reference mass, 3-axis vs 4-axis;
    (b) cumulative distribution of |log10(ratio)|."""
    r3, r4 = m3 / mt, m4 / mt
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.6, 3.6))

    axL.scatter(mt, r3, s=14, facecolors='none', edgecolors='0.45', marker='o',
                linewidths=0.8, label='3-axis', alpha=0.7)
    axL.scatter(mt, r4, s=9, c='crimson', marker='o', label='4-axis',
                alpha=0.7, linewidths=0)
    axL.axhline(1.0, color='k', ls='--', lw=0.7)
    axL.set_xscale('log'); axL.set_yscale('log')
    axL.set_xlabel('reference mass [Msun]'); axL.set_ylabel('M_corrected / M_ref')
    axL.set_ylim(0.1, 12)
    axL.legend(fontsize=8, loc='upper right', frameon=False)
    axL.set_title('(a) corrected ratio vs reference mass', fontsize=9.5)

    for r, c, lab in [(r3, '0.45', '3-axis'), (r4, 'crimson', '4-axis')]:
        acc = np.sort(np.abs(np.log10(r))); frac = np.arange(1, len(acc) + 1) / len(acc)
        axR.plot(acc, frac, color=c, lw=1.7, label=lab)
    axR.axvline(np.log10(2), color='k', ls=':', lw=0.8)
    axR.set_xlim(0, 0.9); axR.set_ylim(0, 1.02)
    axR.set_xlabel('|log10(M_corrected / M_ref)|'); axR.set_ylabel('cumulative fraction')
    axR.legend(fontsize=8, loc='lower right', frameon=False)
    axR.set_title('(b) accuracy CDF (dotted: factor of 2)', fontsize=9.5)
    axR.grid(alpha=0.25, lw=0.5)

    fig.suptitle(title, fontsize=10.5)
    plt.tight_layout()
    path = basename + '.pdf'
    plt.savefig(path)
    plt.close(fig)
    print("  wrote plot: %s" % path)


def _report(label, rows, basename=None):
    """rows: list of (true_mass, uncorrected_mass, m_corrected_3axis,
    m_corrected_4axis). Prints the standard summary table, and if
    `basename` is given, also writes a PDF plot and a space-separated
    table of the underlying per-source values."""
    if not rows:
        print("%s: no matched sources -- nothing to report" % label)
        return None
    mt, mi, m3, m4 = map(np.array, zip(*rows))
    r_unc, r3, r4 = mi / mt, m3 / mt, m4 / mt
    print("\n=== %s (n=%d) ===" % (label, len(rows)))
    for lab, r in (('uncorrected', r_unc), ('3-axis', r3), ('4-axis', r4)):
        lr = np.log10(r)
        print("  %-14s median=%.3f  scatter=factor %.2f  within2x=%.0f%%  "
              "median|logratio|=%.3f" % (
                  lab, np.median(r), 10 ** lr.std(),
                  100 * np.mean((r > 0.5) & (r < 2)), np.median(np.abs(lr))))
    d3, d4 = np.abs(np.log10(r3)), np.abs(np.log10(r4))
    print("  paired: 4-axis improved %d, worsened %d (of %d)"
          % (np.sum(d4 < d3), np.sum(d4 > d3), len(r3)))
    if basename:
        _write_table(basename, mt, mi, m3, m4)
        _plot_validation(basename, "%s (n=%d)" % (label, len(rows)), mt, mi, m3, m4)
    return dict(mtrue=mt, minit=mi, mcorr3=m3, mcorr4=m4)


# =========================================================================
# test 1: injection (cSD blocks)
# =========================================================================
CSD_BLOCKS = [
    ('01', 'x0.25'), ('02', 'x0.5'), ('03', ''), ('04', 'x2'), ('05', 'x4'),
]

def test_injection(tests_dir='tests'):
    """The 5 cSD injection blocks (0.25x-4x mass/density scaling).
    Uses correct_cloud() directly -- no re-implemented file reading."""
    all_rows = []
    for n, x in CSD_BLOCKS:
        ext = ("%s/Aquila.s.sources.ok.cat=Aquila.s.sources.ok.add.cat="
               "thin.Aquila.s.sources.ok.00.cat_cSD_%s" % (tests_dir, n))
        tr = (("%s/inj_all24_s1111_%s_truth.txt_cSD_%s" % (tests_dir, x, n))
              if x else ("%s/inj_all24_s1111_truth.txt_cSD_%s" % (tests_dir, n)))
        rows = _match_one_injection(ext, tr, "cSD_%s" % n)
        all_rows += rows
    return _report("Injection (5 cSD blocks, diverse models)", all_rows,
                   basename="validation_injection")


def _match_one_injection(ext, truth_path, label):
    """Run the pipeline's own correct_cloud() on one injected field, then
    match its output sources to the truth table by position (<10 px)."""
    truth = _loadtruth(truth_path)
    table, _ = mc.correct_cloud(ext, DIST_AQUILA)
    # positions: read from the main catalogue only. For the cSD-style single
    # concatenated filename (e.g. "cat=add.cat=thin.cat_cSD_01"), the '=' is
    # part of the literal filename on disk, not a path separator -- read
    # that one file directly. For a genuine list [cat, add.cat] (bgscan),
    # read only the FIRST element: add.cat has a different column layout at
    # the same positions, and reading it the same way silently corrupts
    # positions -- this was the second bug found in the original ad hoc code.
    main_cat = ext if isinstance(ext, str) else ext[0]
    pix = {}
    for l in open(main_cat):
        s = l.strip()
        if not s or s[0] in '#!':
            continue
        f = s.split()
        try:
            pix[int(float(f[0]))] = (float(f[4]), float(f[5]))
        except (ValueError, IndexError):
            pass
    rows = []
    for i, no in enumerate(table['NO']):
        no = int(no)
        if no not in pix:
            continue
        xp, yp = pix[no]
        d = np.hypot(truth[:, 0] - xp, truth[:, 1] - yp)
        j = int(np.argmin(d))
        if d[j] >= 10:
            continue
        if not np.isfinite(table['C_M_4D'][i]):
            continue
        rows.append((truth[j, 2], table['M_init'][i], table['M_corr'][i],
                     table['M_corr_4D'][i]))
    mc.vlog(0, "  %s: %d truth rows, %d matched" % (label, len(truth), len(rows)))
    return rows


# =========================================================================
# test 2: background-scan (7 fixed models, varying real background)
# =========================================================================
BGSCAN_MODELS = ['bgscan_i1j02k05', 'bgscan_i1j07k06', 'bgscan_i1j08k08',
                  'bgscan_i1j12k08', 'bgscan_i1j13k06', 'bgscan_i2j07k05',
                  'bgscan_i3j15k07']

def test_bgscan(tests_dir='tests'):
    all_rows = []
    for d in BGSCAN_MODELS:
        ext = ["%s/%s/Aquila.s.sources.ok.cat" % (tests_dir, d),
               "%s/%s/Aquila.s.sources.ok.add.cat" % (tests_dir, d)]
        tr = "%s/%s/inj_%s_truth.txt" % (tests_dir, d, d)
        rows = _match_one_injection(ext, tr, d)
        all_rows += rows
    return _report("Background-scan (7 fixed models, varying real bg)", all_rows,
                   basename="validation_bgscan")


# =========================================================================
# test 3: independent subdivision grid
# =========================================================================
def test_subdivision(cat_path='cats/bes_model_grid_subdivision2_catalog',
                      rectab_path='cats/bes_model_grid_subdivision2_catalog_recovery_tables',
                      one_per_node=True):
    """Independent, diverse-concentration subdivision grid. Uses the
    pipeline's own build_grids()/build_recovery_samples() -- no
    re-implemented interpolation."""
    G = mc.build_grids()
    St = mc.build_recovery_samples(cat_path, rectab_path)
    good = np.isfinite(St['rbe_pc']) & (St['rbe_pc'] > 0)
    node = St['node'][good]; sig = St['sig'][good]; conc = St['conc'][good]
    cfoot = St['cfoot'][good]; rbe = St['rbe_pc'][good]
    mrec = St['mrec'][good]; mbe = St['mbe'][good]

    lg = np.log10
    if one_per_node:
        key = np.array([n[0] * 100000 + n[1] * 1000 + n[2] for n in node])
        uniq = np.unique(key)
        REP = 1.84
        idxs = []
        for k in uniq:
            j = np.where(key == k)[0]
            idxs.append(j[np.argmin(np.abs(cfoot[j] - REP))])
        idxs = np.array(idxs)
        sig, conc, cfoot, rbe, mrec, mbe = (a[idxs] for a in
            (sig, conc, cfoot, rbe, mrec, mbe))

    X3 = np.column_stack([lg(sig), lg(conc), lg(cfoot)])
    X4 = np.column_stack([lg(sig), lg(conc), lg(cfoot), lg(rbe)])
    pc3 = G['cm'](X3); pc4 = G['cm4'](X4)
    ok = np.isfinite(pc3) & np.isfinite(pc4)
    rows = list(zip(mbe[ok], mrec[ok], mrec[ok] * 10 ** pc3[ok], mrec[ok] * 10 ** pc4[ok]))
    return _report("Subdivision grid (%s)" % ('one/node' if one_per_node else 'all samples'),
                   rows, basename="validation_subdivision")


# =========================================================================
# test 4: node-level leave-one-out
# =========================================================================
def test_loo(cat_path=mc.CAT67, rectab_path=mc.RECTAB, checkpoint='loo_checkpoint.pkl'):
    """Full node-level LOO: remove each node's samples, rebuild both
    correctors from the rest, predict the withheld node's own
    representative sample. Slow (3-axis + 4-axis rebuild per node);
    checkpoints every 20 nodes so it can be safely interrupted and
    resumed by re-running this function."""
    import os
    from scipy.interpolate import LinearNDInterpolator
    lg = np.log10
    S = mc.build_recovery_samples(cat_path, rectab_path)
    good = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    node = S['node'][good]; sig = S['sig'][good]; conc = S['conc'][good]
    cfoot = S['cfoot'][good]; rbe = S['rbe_pc'][good]
    mrec = S['mrec'][good]; mbe = S['mbe'][good]
    X3 = np.column_stack([lg(sig), lg(conc), lg(cfoot)])
    X4 = np.column_stack([lg(sig), lg(conc), lg(cfoot), lg(rbe)])
    y = lg(S['cm'][good])
    key = np.array([n[0] * 100000 + n[1] * 1000 + n[2] for n in node])
    uniq = np.unique(key)

    if os.path.exists(checkpoint):
        results, done = pickle.load(open(checkpoint, 'rb')).values()
        print("resuming LOO from checkpoint: %d/%d nodes done" % (done, len(uniq)))
    else:
        results, done = [], 0

    REP = 1.84
    t0 = time.time()
    i = done
    while i < len(uniq):       # and time.time() - t0 < 240
        k = uniq[i]
        mask_out = key == k
        mask_in = ~mask_out
        itp3 = LinearNDInterpolator(X3[mask_in], y[mask_in])
        itp4 = LinearNDInterpolator(X4[mask_in], y[mask_in])
        idxs = np.where(mask_out)[0]
        j = idxs[np.argmin(np.abs(cfoot[idxs] - REP))]
        pc3 = itp3(X3[j]); pc4 = itp4(X4[j])
        results.append((mbe[j], mrec[j],
                         float(pc3[0]) if np.isfinite(pc3[0]) else np.nan,
                         float(pc4[0]) if np.isfinite(pc4[0]) else np.nan))
        i += 1
        if i % 20 == 0:
            pickle.dump(dict(results=results, done=i), open(checkpoint, 'wb'))
            print("  ...%d/%d nodes done (%.0fs elapsed)" % (i, len(uniq), time.time() - t0))

    pickle.dump(dict(results=results, done=i), open(checkpoint, 'wb'))
    if i < len(uniq):
        print("LOO incomplete (%d/%d) -- re-run this command to continue" % (i, len(uniq)))
        return None

    rows = [(mb, mr, mr * 10 ** p3, mr * 10 ** p4) for mb, mr, p3, p4 in results
            if np.isfinite(p3) and np.isfinite(p4)]
    os.remove(checkpoint)
    return _report("Node-level leave-one-out", rows, basename="validation_loo")


# =========================================================================
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('test', choices=['injection', 'subdivision', 'loo', 'bgscan', 'all'])
    ap.add_argument('--verbose', type=int, default=0, choices=(0, 1, 2),
                    help='pipeline verbosity: 0=main steps, 1=+detail, 2=+most detail')
    ap.add_argument('--tests-dir', default='tests')
    a = ap.parse_args()
    mc.VERBOSE = a.verbose

    results = {}
    if a.test in ('injection', 'all'):
        results['injection'] = test_injection(a.tests_dir)
    if a.test in ('bgscan', 'all'):
        results['bgscan'] = test_bgscan(a.tests_dir)
    if a.test in ('subdivision', 'all'):
        results['subdivision'] = test_subdivision()
    if a.test in ('loo', 'all'):
        results['loo'] = test_loo()

    pickle.dump(results, open('validation_results.pkl', 'wb'))
    print("\nAll results saved to validation_results.pkl")
