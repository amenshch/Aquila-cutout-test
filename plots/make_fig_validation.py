#!/usr/bin/env python3
"""
make_fig_validation.py -- generate the paper's actual fig_validation.pdf:
a combined 4-row x 2-column figure (injection, subdivision, node-LOO,
background-scan), each row showing corrected-ratio-vs-reference-mass and
the accuracy CDF, 3-axis vs 4-axis corrector.

This reuses run_validation_tests.py directly rather than re-implementing
any file reading or matching logic -- see that script's own docstring for
why that matters (two real bugs were found and fixed there earlier by
re-implementing what load_getsf()/correct_cloud() already do correctly).

REQUIRES, in the working directory:
  mass_correction_pipeline_4d.py, run_validation_tests.py
  cats/  (main grid + bes_model_grid_subdivision2_catalog(+_recovery_tables))
  tests/ (the 5 cSD injection blocks + 7 bgscan model directories)

Runtime: a few minutes for injection+bgscan+subdivision; the LOO test
alone is ~45-50 min (see run_validation_tests.py's own checkpointing).
Run once with `python3 run_validation_tests.py loo` beforehand if you
don't already have loo results cached, then re-run this script -- it will
reuse validation_results.pkl if present rather than re-running LOO.
"""
import pickle, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import run_validation_tests as rv
import mass_correction_pipeline_4d as mc

CACHE = 'validation_results.pkl'


def get_results():
    """Reuse a previous full run if available (saves ~50 min by not
    re-running LOO); otherwise run all four tests now."""
    if os.path.exists(CACHE):
        d = pickle.load(open(CACHE, 'rb'))
        if all(k in d and d[k] is not None for k in
               ('injection', 'bgscan', 'subdivision', 'loo')):
            print("reusing cached results from %s" % CACHE)
            return d
        print("cache incomplete, re-running missing tests...")
    else:
        d = {}
    if 'injection' not in d or d['injection'] is None:
        d['injection'] = rv.test_injection()
    if 'bgscan' not in d or d['bgscan'] is None:
        d['bgscan'] = rv.test_bgscan()
    if 'subdivision' not in d or d['subdivision'] is None:
        d['subdivision'] = rv.test_subdivision()
    if 'loo' not in d or d['loo'] is None:
        d['loo'] = rv.test_loo()
    pickle.dump(d, open(CACHE, 'wb'))
    return d


def scatter_panel(ax, xtrue, r3, r4, xlabel, ylabel, title, alpha, size):
    ax.scatter(xtrue, r3, s=size, facecolors='none', edgecolors='0.45', marker='o',
               linewidths=0.8, label='3-axis', alpha=alpha)
    ax.scatter(xtrue, r4, s=size * 0.6, c='crimson', marker='o', label='4-axis',
               alpha=alpha, linewidths=0)
    ax.axhline(1.0, color='k', ls='--', lw=0.7)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_ylim(0.1, 12)
    ax.legend(fontsize=6.5, loc='upper right', frameon=False, handletextpad=0.3)
    ax.set_title(title, fontsize=8.8)


def cdf_panel(ax, r3, r4, title):
    for r, c, lab in [(r3, '0.45', '3-axis'), (r4, 'crimson', '4-axis')]:
        acc = np.sort(np.abs(np.log10(r))); frac = np.arange(1, len(acc) + 1) / len(acc)
        ax.plot(acc, frac, color=c, lw=1.6, label=lab)
    ax.axvline(np.log10(2), color='k', ls=':', lw=0.8)
    ax.set_xlim(0, 0.9); ax.set_ylim(0, 1.02)
    ax.set_xlabel(r'$|\log_{10}(M/M_{\rm ref})|$'); ax.set_ylabel('cumulative fraction')
    ax.legend(fontsize=6.5, loc='lower right', frameon=False)
    ax.set_title(title, fontsize=8.8); ax.grid(alpha=0.25, lw=0.5)


def main():
    d = get_results()
    inj, sub, loo, bg = d['injection'], d['subdivision'], d['loo'], d['bgscan']

    plt.rcParams.update({'font.size': 9.1, 'axes.linewidth': 0.8})
    fig, axes = plt.subplots(4, 2, figsize=(7.2, 12.4))

    r3i, r4i = inj['mcorr3'] / inj['mtrue'], inj['mcorr4'] / inj['mtrue']
    scatter_panel(axes[0, 0], inj['mtrue'], r3i, r4i, r'$M_{\rm true}$ [$M_\odot$]',
                  r'$M/M_{\rm true}$', '(a) Injection, diverse models (n=%d)' % len(r3i),
                  0.6, 16)
    cdf_panel(axes[0, 1], r3i, r4i, '(b) Injection: accuracy CDF')

    r3s, r4s = sub['mcorr3'] / sub['mtrue'], sub['mcorr4'] / sub['mtrue']
    scatter_panel(axes[1, 0], sub['mtrue'], r3s, r4s, r'$M_{\rm BE}$ [$M_\odot$]',
                  r'$M/M_{\rm BE}$', '(c) Subdivision (n=%d)' % len(r3s), 0.5, 10)
    cdf_panel(axes[1, 1], r3s, r4s, '(d) Subdivision: accuracy CDF')

    r3l, r4l = loo['mcorr3'] / loo['mtrue'], loo['mcorr4'] / loo['mtrue']
    scatter_panel(axes[2, 0], loo['mtrue'], r3l, r4l, r'$M_{\rm BE}$ [$M_\odot$]',
                  r'$M/M_{\rm BE}$', '(e) Node-LOO (n=%d)' % len(r3l), 0.25, 8)
    cdf_panel(axes[2, 1], r3l, r4l, '(f) Node-LOO: accuracy CDF')

    r3b, r4b = bg['mcorr3'] / bg['mtrue'], bg['mcorr4'] / bg['mtrue']
    scatter_panel(axes[3, 0], bg['mtrue'], r3b, r4b, r'$M_{\rm true}$ [$M_\odot$]',
                  r'$M/M_{\rm true}$', '(g) bg-scan, 7 fixed models (n=%d)' % len(r3b),
                  0.6, 16)
    cdf_panel(axes[3, 1], r3b, r4b, '(h) bg-scan: accuracy CDF')

    plt.tight_layout()
    plt.savefig('fig_validation.pdf')
    print('saved fig_validation.pdf')


if __name__ == '__main__':
    main()
