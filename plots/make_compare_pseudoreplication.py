#!/usr/bin/env python3
"""
make_compare_pseudoreplication.py -- generate compare_pseudoreplication.pdf.

STATUS: EXPLORATORY / NEVER PART OF THE PAPER. This figure compared
counting every truncation-radius sample as an independent test case
("pseudo-replicated") against counting one representative sample per
physical node ("independent"), for both the subdivision and LOO tests.
The conclusion of that comparison -- use one-sample-per-node counting for
headline validation statistics -- WAS adopted and is now baked directly
into run_validation_tests.py (test_subdivision(one_per_node=True) and
test_loo(), which only ever report one sample per node). The comparison
figure itself was exploratory and was never included in the paper.

NOT directly runnable as-is: depends on four cached files
(subdiv_data.pkl, subdiv_onepernode.pkl, loo_data.pkl,
loo_onepernode.pkl) from the pre-pathology-fix, old-subdivision-grid
state of the project (same caveat as make_fig_capping_tradeoff.py).
Reproducing the comparison on current data would mean running
test_subdivision() with one_per_node=True and False, and similarly
comparing full vs. one-per-node LOO counting, using
run_validation_tests.py's own functions rather than these old pickles.

Original logic, unmodified below.
"""
import numpy as np, pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sub_pseudo = pickle.load(open('/tmp/subdiv_data.pkl', 'rb'))
sub_one = pickle.load(open('/tmp/subdiv_onepernode.pkl', 'rb'))
loo_pseudo = pickle.load(open('/tmp/loo_data.pkl', 'rb'))
loo_one = pickle.load(open('/tmp/loo_onepernode.pkl', 'rb'))

plt.rcParams.update({'font.size': 9.3, 'axes.linewidth': 0.8})
fig, axes = plt.subplots(4, 2, figsize=(7.4, 12.2))

def scatter_panel(ax, xtrue, r_unc, r_cor, title, alpha_pts=0.5, size=14):
    ax.scatter(xtrue, r_unc, s=size, facecolors='none', edgecolors='0.4', marker='o',
               linewidths=0.9, label='uncorrected', alpha=alpha_pts)
    ax.scatter(xtrue, r_cor, s=size * 0.55, c='crimson', marker='o', label='corrected',
               alpha=alpha_pts, linewidths=0)
    ax.axhline(1.0, color='k', ls='--', lw=0.7)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'$M_{\rm BE}$  [$M_\odot$]'); ax.set_ylabel(r'$M/M_{\rm BE}$')
    ax.set_ylim(0.1, 12)
    ax.legend(fontsize=6.5, loc='upper right', frameon=False, handletextpad=0.3)
    ax.set_title(title, fontsize=9)

def cdf_panel(ax, r_unc, r_cor, title):
    for r, c, lab in [(r_unc, '0.45', 'uncorrected'), (r_cor, 'crimson', 'corrected')]:
        acc = np.sort(np.abs(np.log10(r)))
        frac = np.arange(1, len(acc) + 1) / len(acc)
        ax.plot(acc, frac, color=c, lw=1.6, label=lab)
    ax.axvline(np.log10(2), color='k', ls=':', lw=0.8)
    ax.set_xlim(0, 0.9); ax.set_ylim(0, 1.02)
    ax.set_xlabel(r'$|\log_{10}(M/M_{\rm BE})|$')
    ax.set_ylabel('cumulative fraction')
    ax.legend(fontsize=6.5, loc='lower right', frameon=False)
    ax.set_title(title, fontsize=9)
    ax.grid(alpha=0.25, lw=0.5)

rows = [
    (sub_pseudo, '(a) Subdivision, ALL truncations (n=%d, pseudo-replicated)'),
    (sub_one,    '(b) Subdivision, ONE sample/node (n=%d, independent)'),
    (loo_pseudo, '(c) LOO, ALL truncations (n=%d, pseudo-replicated)'),
    (loo_one,    '(d) LOO, ONE sample/node (n=%d, independent)'),
]
for i, (d, title_fmt) in enumerate(rows):
    mbe, mrec, mcorr = d['mbe'], d['mrec'], d['mcorr']
    r_unc, r_cor = mrec / mbe, mcorr / mbe
    n = len(mbe)
    alpha_pts = 0.5 if n < 1000 else 0.12
    size = 14 if n < 1000 else 5
    scatter_panel(axes[i, 0], mbe, r_unc, r_cor, title_fmt % n, alpha_pts, size)
    cdf_panel(axes[i, 1], r_unc, r_cor, title_fmt % n)

plt.tight_layout()
plt.savefig('compare_pseudoreplication.pdf')
plt.savefig('compare_pseudoreplication.png', dpi=140)
print('saved')
