#!/usr/bin/env python3
"""
make_fig_capping_tradeoff.py -- generate fig_capping_tradeoff.pdf.

STATUS: HISTORICAL / NOT PART OF THE CURRENT PAPER. The capping analysis
this figure supported was superseded by the four-observable corrector
(Sect. "A fourth observable" in the current paper): the paper now
describes this only briefly, as "tested and superseded"
(Sect. "An output cap was tested and superseded"), without including the
figure. This script is provided for historical reference and reproduces
exactly what was run at the time, but it is NOT directly runnable as-is
today: it depends on three cached files (inj_data.pkl, subdiv_factors.pkl,
loo_data.pkl) built from the three-axis-only corrector on the PRE-
pathology-fix production grid and the OLD 36-model, critical-only
subdivision grid -- both since replaced (see the paper's Sect.
"Interpolation and output constraints" for the pathology fix, and
"Validation" for the new subdivision grid). Running this against current
data would require rebuilding those three inputs from
run_validation_tests.py's test_injection()/test_subdivision()/test_loo()
results first (rows = list of (true_mass, uncorrected, corrected) tuples;
adapt cm_pred/cm_true below to whatever the current 3-axis-only
equivalent would be), not a drop-in replacement.

Original logic, unmodified below.
"""
import numpy as np, pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

inj = pickle.load(open('/tmp/inj_data.pkl', 'rb'))
sub = pickle.load(open('/tmp/subdiv_factors.pkl', 'rb'))
loo = pickle.load(open('/tmp/loo_data.pkl', 'rb'))

r_inj = inj['mcorr'] / inj['mtrue']
r_sub = sub['cm_pred'] / sub['cm_true']
r_loo = loo['mcorr'] / loo['mbe']

caps = np.array([1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0])

plt.rcParams.update({'font.size': 9.5, 'axes.linewidth': 0.8})
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))

for r, c, lab in [(r_inj, 'crimson', 'Injection (n=108)'),
                   (r_sub, 'tab:blue', 'Subdivision (n=873)'),
                   (r_loo, '0.4', 'Node-LOO (n=8629)')]:
    # typical error, as a multiplicative factor equivalent to delta_rms
    rms_factor = [10**np.sqrt(np.mean(np.log10(np.minimum(r, cap))**2)) for cap in caps]
    f2x = [100 * np.mean((np.minimum(r, cap) > 0.5) & (np.minimum(r, cap) < 2)) for cap in caps]
    axes[0].plot(caps, rms_factor, 'o-', ms=3.5, color=c, label=lab)
    axes[1].plot(caps, f2x, 'o-', ms=3.5, color=c, label=lab)

for r, c in [(r_inj, 'crimson'), (r_sub, 'tab:blue'), (r_loo, '0.4')]:
    axes[0].axhline(10**np.sqrt(np.mean(np.log10(r)**2)), color=c, ls=':', lw=0.8)

axes[0].set_xscale('log')
axes[0].set_xlabel(r'output cap on $M_{\rm corrected}/M_{\rm ref}$')
axes[0].set_ylabel('typical error, as a factor')
axes[0].legend(fontsize=7, frameon=False)
axes[0].set_title('(a) Total error vs. cap', fontsize=9.5)

axes[1].set_xscale('log')
axes[1].set_xlabel(r'output cap on $M_{\rm corrected}/M_{\rm ref}$')
axes[1].set_ylabel('fraction within factor 2  [%]')
axes[1].set_title('(b) Accuracy vs. cap', fontsize=9.5)
axes[1].legend(fontsize=7, frameon=False)

plt.tight_layout()
plt.savefig('fig_capping_tradeoff.pdf')
print('saved')
