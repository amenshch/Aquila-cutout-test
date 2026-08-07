#!/usr/bin/env python3
"""
make_fig_grid_coverage.py -- generate fig_grid_coverage.pdf: the model
grid's coverage in the (H_intr, M_BE) plane, colored by concentration
(xi_max), with analytic critical curves overlaid per Sigma_cloud.

FULLY SELF-CONTAINED: unlike the version originally run in chat (which
built on several cached .npz files from earlier, separate steps), this
script reconstructs everything from scratch in one pass:
  1. Solve the isothermal Lane-Emden equation once (scipy solve_ivp).
  2. Tabulate phi(xi_max) -- the projected half-maximum width in units of
     r0 -- via a numerical Abel projection of the Lane-Emden density.
  3. For every node in the production catalog, invert its known
     (M_BE, rho_BE, T_BE) for xi_max via the tabulated C(xi_max).
  4. Convert to H_intr [pc] using the catalog's own R_BE_pc column
     (H_intr = R_BE_pc * phi(xi_max)/xi_max).
  5. Compute the analytic critical (xi_max=6.451) curve for each of the
     grid's six Sigma_cloud columns, by iterating the self-consistent
     T_BE(N_tot) relation.
  6. Plot.

Requires: numpy, scipy, matplotlib, and cats/bes_model_grid_final2_catalog
in the working directory (or edit CAT67 below).

Runtime: ~1-2 minutes (the critical-curve iteration is the slow part).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

CAT67 = 'cats/bes_model_grid_final2_catalog'
trapz = getattr(np, 'trapezoid', None) or np.trapz

# physical constants (cgs)
G_CGS = 6.67430e-8; KB = 1.380649e-16; MH = 1.6726219e-24
MSUN = 1.98892e33; MU_P = 2.33; MU_H2 = 2.8; PC_TO_CM = 3.0856775814913673e18
PC_TO_AU = 206265.0

# ---------------------------------------------------------------------
# 1. Lane-Emden equation, solved once
# ---------------------------------------------------------------------
print("Solving the isothermal Lane-Emden equation...")
def _rhs(xi, y):
    psi, dpsi = y
    return [dpsi, np.exp(-psi) - 2.0 * dpsi / xi]

xi0 = 1e-6
y0 = [xi0**2 / 6.0, xi0 / 3.0]      # small-xi series expansion
sol = solve_ivp(_rhs, [xi0, 200.0], y0, max_step=0.01, rtol=1e-10, atol=1e-12,
                 dense_output=True)
xi_grid = np.concatenate([[0.0], np.geomspace(xi0, 200.0, 4000)])
psi = np.concatenate([[0.0], sol.sol(xi_grid[1:])[0]])
dpsi = np.concatenate([[0.0], sol.sol(xi_grid[1:])[1]])
C = xi_grid**2 * dpsi                # C(xi): dimensionless enclosed mass, Eq. (Cxi)

i_crit = np.argmin(np.abs(xi_grid - 6.451))
print("  sanity check: density contrast at xi=6.451 is %.2f (expect ~14.0)"
      % np.exp(psi[i_crit]))
print("  sanity check: C(6.451) = %.4f (expect ~15.7086)" % np.interp(6.451, xi_grid, C))

Cinterp = interp1d(xi_grid, C, kind='cubic')
xi_of_C = interp1d(C, xi_grid, kind='cubic', bounds_error=False,
                    fill_value=(0.0, xi_grid[-1]))
psi_of_xi = interp1d(xi_grid, psi, kind='cubic')

# ---------------------------------------------------------------------
# 2. phi(xi_max): projected half-maximum width, via Abel projection
# ---------------------------------------------------------------------
print("Tabulating phi(xi_max) (projected half-maximum width)...")
def _phi_of_ximax(ximax, n_p=80, n_z=800):
    p = np.linspace(1e-4, ximax * 0.999, n_p)
    Sigma = np.empty(n_p)
    for i, pi in enumerate(p):
        zmax = np.sqrt(max(ximax**2 - pi**2, 0))
        z = np.linspace(0, zmax, n_z)
        r = np.clip(np.sqrt(pi**2 + z**2), 0, xi_grid[-1])
        Sigma[i] = 2 * trapz(np.exp(-psi_of_xi(r)), z)
    Sigma0 = Sigma[0]
    below = np.where(Sigma <= Sigma0 / 2)[0]
    if len(below) == 0:
        return np.nan
    i = below[0]
    if i == 0:
        return 0.0
    p_half = np.interp(Sigma0 / 2, [Sigma[i], Sigma[i - 1]], [p[i], p[i - 1]])
    return 2 * p_half     # FWHM in units of r0

xi_tab = np.geomspace(0.3, 65, 45)
phi_tab = np.array([_phi_of_ximax(x) for x in xi_tab])
print("  phi saturates at xi_max=%.0f: phi=%.3f (paper reports ~6.6)"
      % (xi_tab[-1], phi_tab[-1]))
phi_of_ximax = interp1d(xi_tab, phi_tab, kind='cubic', bounds_error=False,
                         fill_value=(phi_tab[0], phi_tab[-1]))

# ---------------------------------------------------------------------
# 3-4. Invert every catalog node for xi_max, then compute H_intr [pc]
# ---------------------------------------------------------------------
print("Loading grid catalog and inverting xi_max per node...")
rows = []
for l in open(CAT67):
    if l.startswith('#') or not l.strip():
        continue
    f = l.split()
    if len(f) < 67:
        continue
    i, j, k = int(f[1]), int(f[2]), int(f[3])
    SD_emb, T_BE, rho_BE, M_BE = float(f[4]), float(f[5]), float(f[6]), float(f[7])
    R_BE_pc = float(f[10])           # column 11 (1-based): physical truncation radius
    rows.append((SD_emb, T_BE, rho_BE, M_BE, R_BE_pc))
rows = np.array(rows)
SD_emb, T_BE, rho_BE, M_BE, R_BE_pc = rows.T
print("  n nodes:", len(rows))

cs2 = KB * T_BE / (MU_P * MH)
r0 = np.sqrt(cs2 / (4 * np.pi * G_CGS * rho_BE))
Ctarget = (M_BE * MSUN) / (4 * np.pi * rho_BE * r0**3)
ok = (Ctarget >= C.min()) & (Ctarget <= C.max())
print("  nodes with C(xi_max) in the solvable range:", ok.sum(), "/", len(Ctarget))
xi_max = np.full(len(Ctarget), np.nan)
xi_max[ok] = xi_of_C(Ctarget[ok])

H_intr_pc = R_BE_pc * phi_of_ximax(xi_max) / xi_max
good = np.isfinite(H_intr_pc) & np.isfinite(xi_max)
M_BE, H_intr_pc, xi_max, SD_emb = M_BE[good], H_intr_pc[good], xi_max[good], SD_emb[good]

# ---------------------------------------------------------------------
# 5. Analytic critical (xi_max=6.451) curve per Sigma_cloud column
# ---------------------------------------------------------------------
print("Computing analytic critical curves per Sigma_cloud...")
e_psi = np.exp(-psi)
I_cum = np.concatenate([[0.0], np.cumsum(0.5 * (e_psi[1:] + e_psi[:-1]) * np.diff(xi_grid))])
I_of_xi = interp1d(xi_grid, I_cum, kind='cubic')
XI_CRIT, C_CRIT = 6.451, 15.7086
I_CRIT = float(I_of_xi(XI_CRIT))
PHI_CRIT = float(phi_of_ximax(XI_CRIT))

a0, a1, a2 = 1021.4, -86.47, 1.842    # T_BE(N_tot) fit coefficients (Eq. TofN)
def T_of_Ntot(Ntot):
    x = np.clip(np.log10(Ntot), 21.5, 23.4)
    return a0 + a1 * x + a2 * x**2

def critical_model_H(M_msun, Sigma_cloud, T_guess=10.0, niter=30):
    """Self-consistent H_intr [pc] of a CRITICAL (xi_max=6.451) sphere of
    mass M_msun, embedded in a cloud of column Sigma_cloud."""
    M = M_msun * MSUN; T = T_guess
    for _ in range(niter):
        cs2 = KB * T / (MU_P * MH)
        rho_c = (4 * np.pi * C_CRIT * cs2**1.5 / ((4 * np.pi * G_CGS)**1.5 * M))**2
        r0 = np.sqrt(cs2 / (4 * np.pi * G_CGS * rho_c))
        N_core = rho_c * r0 / (MU_H2 * MH) * I_CRIT
        T_new = T_of_Ntot(Sigma_cloud + N_core)
        if abs(T_new - T) < 1e-5:
            T = T_new; break
        T = T_new
    return PHI_CRIT * r0 / PC_TO_CM

SDvals = [3e21, 6e21, 1.2e22, 2.4e22, 4.8e22, 9.6e22]   # the grid's own 6 columns
Mgrid = np.geomspace(0.008, 28, 60)
curves = {}
for SD in SDvals:
    curves[SD] = np.array([critical_model_H(M, SD) for M in Mgrid])
    print("  Sigma_cloud=%.1e: H_intr range %.4f-%.4f pc" % (SD, curves[SD].min(), curves[SD].max()))

# ---------------------------------------------------------------------
# 6. Plot
# ---------------------------------------------------------------------
print("Plotting...")
plt.rcParams.update({'font.size': 10, 'axes.linewidth': 0.9})
fig, ax = plt.subplots(figsize=(4.6, 4.2))

norm = Normalize(vmin=np.log10(xi_max.min()), vmax=np.log10(xi_max.max()))
sc = ax.scatter(H_intr_pc, M_BE, c=np.log10(xi_max), cmap='viridis', norm=norm,
                 s=10, alpha=0.85, linewidths=0, rasterized=True, zorder=2)

reds = plt.cm.Reds(np.linspace(0.45, 0.95, len(SDvals)))
for idx, SD in enumerate(SDvals):
    F = curves[SD]
    m = (F >= 1e-3) & (F <= 0.5)
    ax.plot(F[m], Mgrid[m], color=reds[idx], lw=1.3, zorder=4)
leg_line = Line2D([0], [0], color=reds[len(reds) // 2], lw=1.3)
ax.legend([leg_line], [r'critical ($\xi_{\max}=6.451$), analytic, per $\Sigma_{\rm cloud}$'],
          loc='upper left', fontsize=7.5, frameon=False)

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'$H_{\rm intr}$ (intrinsic, unconvolved)  [pc]')
ax.set_ylabel(r'$M_{\rm BE}$  [$M_\odot$]')
ax.set_xlim(1e-3, 0.5)
ax.set_ylim(3e-3, 60)

def pc2au(x): return x * PC_TO_AU
def au2pc(x): return x / PC_TO_AU
secax = ax.secondary_xaxis('top', functions=(pc2au, au2pc))
secax.set_xlabel(r'$H_{\rm intr}$  [AU]')

cb = fig.colorbar(sc, ax=ax, pad=0.02)
cb.set_label(r'$\log_{10}\,\xi_{\max}$ (concentration)')

plt.tight_layout()
plt.savefig('fig_grid_coverage.pdf')
print("Saved fig_grid_coverage.pdf")
