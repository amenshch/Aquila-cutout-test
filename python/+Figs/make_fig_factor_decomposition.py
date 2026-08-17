"""Regenerate fig_factor_decomposition.pdf on the corrected training set.

C_M decomposes as C_M = F_flux * F_bg * F_temp, where F_flux is the mass lost to
footprint truncation, F_bg the effect of getsf's interpolated background
relative to an uncontaminated annulus, and F_temp whatever remains, which
isolates the bias from fitting a single temperature to a source with an
internal temperature gradient.  Each factor is predicted by node-level
leave-one-out on the adopted four observables and compared with its true value.
"""
import sys, numpy as np
sys.path.insert(0, 'python')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import aa_plotstyle; aa_plotstyle.apply()
from aa_plotstyle import clean_log, panel_label
from scipy.spatial import cKDTree
import mass_correction_pipeline_5 as mc

mc.VERBOSE = 0
S = mc.build_recovery_samples('cats/bes_model_grid_final2_catalog',
                              'cats/bes_model_grid_final2_recovery_tables_v3')
lg = np.log10
theta = S['rbe_pc'] * 206264.806 / mc.GRID_DISTANCE_PC
X = np.column_stack([lg(S['sig']), lg(S['conc']), lg(S['cfoot']), lg(theta)])
Ftemp = S['cm'] / (S['f_flux'] * S['f_bg'])
targets = [('F_{\\rm flux}', S['f_flux'], 'crimson'),
           ('F_{\\rm bg}',   S['f_bg'],   'darkorange'),
           ('F_{\\rm temp}', Ftemp,       'steelblue')]
nk = np.array([str(tuple(k)) for k in S['node']])
uniq, nid = np.unique(nk, return_inverse=True)
K = 320

def loo(y):
    ok = np.all(np.isfinite(X), 1) & np.isfinite(y) & (y > 0)
    Xo, yo, ni = X[ok], lg(y[ok]), nid[ok]
    mb, cf = S['mbe'][ok], S['cfoot'][ok]
    Xs = Xo / Xo.std(0); tree = cKDTree(Xs)
    out_p, out_t, out_m = [], [], []
    for u in np.unique(ni):
        i = np.where(ni == u)[0]
        i = int(i[np.argmin(np.abs(cf[i] - 1.84))])
        keep = ni != u
        dd, cd = tree.query(Xs[i], k=min(K + 600, len(Xs)))
        m = keep[cd]; cd, dd = cd[m][:K], dd[m][:K]
        h = dd[-1] if dd[-1] > 0 else 1.0
        w = np.maximum(np.clip(1 - (dd / h) ** 3, 0, None) ** 3, 1e-6)
        A = np.column_stack([np.ones(len(cd)), Xs[cd] - Xs[i]]); sw = np.sqrt(w)
        b, *_ = np.linalg.lstsq(A * sw[:, None], yo[cd] * sw, rcond=None)
        out_p.append(10 ** b[0]); out_t.append(10 ** yo[i]); out_m.append(mb[i])
    return map(np.array, (out_p, out_t, out_m))

fig, axes = plt.subplots(3, 2, figsize=(7.0, 5.6), sharex='col')
lab = 'abcdef'
for row, (name, y, col) in enumerate(targets):
    p, t, m = loo(y)
    r = p / t
    good = np.isfinite(r) & (r > 0)
    sc = 10 ** np.std(lg(r[good]))
    sym = r'$%s^{\rm pred}/%s^{\rm true}$' % (name, name)

    ax = axes[row, 0]
    ax.scatter(m, r, s=6, color=col, alpha=0.45, lw=0)
    ax.axhline(1.0, ls='--', color='k', lw=0.8)
    ax.set_xscale('log'); ax.set_yscale('log')
    clean_log(ax)
    if row == 2:
        ax.set_xlabel(r'$M_{\rm BE}$  ($M_{\odot}$)', fontsize=8)
    ax.set_ylabel(sym, fontsize=8.5)
    ax.set_ylim(0.5, 2.0)
    panel_label(ax, '(%s) prediction accuracy' % lab[2 * row],
                loc='upper left', pad=0.045)
    ax.tick_params(labelsize=7)

    ax = axes[row, 1]
    v = np.sort(r[good])
    ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=col, lw=1.5)
    ax.axvline(1.0, ls='--', color='k', lw=0.8)
    ax.set_xscale('log')
    ax.set_xlim(0.5, 2.0); ax.set_ylim(0, 1.20)
    clean_log(ax, 'x')
    if row == 2:
        # the axis is shared down the column, so the label must be generic;
        # each row is identified by the symbol on its left-hand panel
        ax.set_xlabel('predicted / true', fontsize=8.5)
    ax.set_ylabel('cumulative fraction', fontsize=8)
    panel_label(ax, '(%s) scatter %.2f' % (lab[2 * row + 1], sc),
                loc='upper left', pad=0.045)
    ax.tick_params(labelsize=7)
    out = int(((r[good] > 2.0) | (r[good] < 0.5)).sum())
    print('%-14s scatter factor %.2f  16-84 pct %.3f-%.3f  '
          'range %.3f-%.2f  outside the plotted range: %d of %d'
          % (name, sc, np.percentile(r[good], 16), np.percentile(r[good], 84),
             r[good].min(), r[good].max(), out, len(r[good])))
fig.tight_layout(h_pad=0.25, w_pad=0.9)
fig.subplots_adjust(hspace=0.06)
fig.savefig('plots/fig_factor_decomposition.pdf')
