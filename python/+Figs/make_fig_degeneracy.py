"""Regenerate fig_degeneracy.pdf on the corrected training set.

For every training sample the nearest neighbour from a DIFFERENT physical node
is found in the normalized observable space; pairs are binned by separation and
the root-mean-square difference in log10(C_M,true) is plotted against the median
separation of the bin.  Interpolation error falls as the separation goes to
zero; a degeneracy does not.
"""
import sys, numpy as np
sys.path.insert(0, 'python')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import aa_plotstyle
aa_plotstyle.apply()
from aa_plotstyle import clean_log
from scipy.spatial import cKDTree
import mass_correction_pipeline_5 as mc

mc.VERBOSE = 0
S = mc.build_recovery_samples('cats/bes_model_grid_final2_catalog',
                              'cats/bes_model_grid_final2_recovery_tables_v3')
lg = np.log10
GRID_D = mc.GRID_DISTANCE_PC
theta = S['rbe_pc'] * 206264.806 / GRID_D
y = lg(S['cm'])
node = np.array([hash(tuple(k)) for k in S['node']])

def curve(cols, nbin=10):
    X = np.column_stack([lg(c) for c in cols])
    ok = np.all(np.isfinite(X), 1) & np.isfinite(y)
    Xn = (X[ok] - X[ok].mean(0)) / X[ok].std(0)
    yy, nn = y[ok], node[ok]
    K = 25
    dd, ii = cKDTree(Xn).query(Xn, k=K)
    sep, dif = [], []
    for m in range(len(Xn)):
        for q in range(1, K):
            if nn[ii[m, q]] != nn[m]:
                sep.append(dd[m, q]); dif.append(yy[m] - yy[ii[m, q]]); break
    sep, dif = np.array(sep), np.array(dif)
    e = np.percentile(sep, np.linspace(0, 100, nbin + 1))
    xs, ys, lo, hi = [], [], [], []
    for i in range(nbin):
        m = (sep >= e[i]) & (sep <= e[i + 1])
        if m.sum() < 30: continue
        xs.append(np.median(sep[m])); ys.append(dif[m].std())
        b = [np.random.default_rng(i).choice(dif[m], m.sum()).std() for _ in range(60)]
        lo.append(np.percentile(b, 16)); hi.append(np.percentile(b, 84))
    return map(np.array, (xs, ys, lo, hi))

base = [S['sig'], S['conc'], S['cfoot']]
x3, y3, l3, h3 = curve(base)
x4, y4, l4, h4 = curve(base + [theta])

# pipeline numerical floor: nearest neighbour in the model's own parameters
kB, amu, G, PC, mu = 1.380649e-16, 1.66053906660e-24, 6.67430e-8, 3.0856775814913673e18, 2.33
phys = {}
for l in open('cats/bes_model_grid_final2_catalog'):
    if not l.strip() or l.lstrip().startswith('#'): continue
    f = l.split()
    T, rho, Rpc = float(f[5]), float(f[6]), float(f[10])
    r0 = np.sqrt(kB * T / (mu * amu)) / np.sqrt(4 * np.pi * G * rho)
    phys[(int(f[1]), int(f[2]), int(f[3]))] = (float(f[4]), float(f[7]),
                                               Rpc * PC / r0)
K = [tuple(k) for k in S['node']]
# The numerical floor is the difference in C_M between samples that are
# nearest neighbours in the model's OWN four physical degrees of freedom
# (embedding column, mass, truncation parameter, and applied truncation
# radius over R_BE), which by construction must have nearly identical C_M.
# The applied truncation radius is not returned by the production pipeline,
# so it is taken from the diagnostic build.
import mass_correction_pipeline_4d as mcd
mcd.VERBOSE = 0
Sd = mcd.build_recovery_samples('cats/bes_model_grid_final2_catalog',
                                'cats/bes_model_grid_final2_recovery_tables_v3')
Kd = [tuple(k) for k in Sd['node']]
rrel = Sd['rtrunc'] / (Sd['rbe_pc'] * 206264.806 / GRID_D)
P = np.column_stack([lg([phys[k][0] for k in Kd]), lg([phys[k][1] for k in Kd]),
                     lg([phys[k][2] for k in Kd]), lg(rrel)])
y = lg(Sd['cm'])
ok = np.all(np.isfinite(P), 1) & np.isfinite(y)
Pn = (P[ok] - P[ok].mean(0)) / P[ok].std(0)
dd, ii = cKDTree(Pn).query(Pn, k=2)
floor = (y[ok] - y[ok][ii[:, 1]]).std()

fig, ax = plt.subplots(figsize=(3.4, 2.9))
ax.fill_between(x3, l3, h3, color='0.55', alpha=0.25, lw=0)
ax.plot(x3, y3, 'o-', color='0.35', ms=3.5, lw=1.4,
        label=r'$\Sigma_{\rm cloud},\ \zeta,\ \phi$')
ax.fill_between(x4, l4, h4, color='crimson', alpha=0.22, lw=0)
ax.plot(x4, y4, 'o-', color='crimson', ms=3.5, lw=1.4,
        label=r'$+\ \theta_{\rm BE}$')
ax.axhline(floor, ls=':', color='k', lw=1.0)
ax.text(x3[-1] * 0.95, floor * 1.12, 'pipeline numerical floor',
        fontsize=6.5, color='k', ha='right')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('separation to nearest sample from a different\nmodel '
              '(normalized units)', fontsize=7.5)
ax.set_ylabel(r'rms $\Delta\log_{10}C_M$  (dex)', fontsize=8)
ax.tick_params(labelsize=7)
ax.legend(fontsize=7, frameon=False, loc='lower left',
          bbox_to_anchor=(0.02, 0.16))
ax.set_ylim(0.01, 0.3)
ax.set_xlim(1e-3, 0.3)
clean_log(ax)
fig.tight_layout()
fig.savefig('plots/fig_degeneracy.pdf')
print("three observables : %.4f dex closest bin, %.4f farthest" % (y3[0], y3[-1]))
print("four observables  : %.4f dex closest bin, %.4f farthest" % (y4[0], y4[-1]))
print("numerical floor   : %.4f dex" % floor)
