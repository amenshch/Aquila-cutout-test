#!/usr/bin/env python3
"""
make_bes_grid.py -- generate the pressure-truncated Bonnor-Ebert model grid
                    catalogue `bes_grid_final.txt`.

The grid is defined on three OBSERVABLE axes -- the embedding-cloud column
Sigma_emb, the core mass M_BE, and the intrinsic FWHM -- because those are the
quantities a source extraction measures.  Everything internal to the sphere
(xi_max, rho_c, r0, R_BE, T_BE) is DERIVED per node, so the grid is in physical
units and is region-independent.

Per node the construction is:

  1. Invert (M, FWHM) -> xi_max.  For a truncated BE sphere,
         FWHM * c_s^2 / (G M) = phi(xi_max) / Cm(xi_max),
     where phi(xi) = projected (column) FWHM in units of r0 and
     Cm(xi) = xi^2 dpsi/dxi is the dimensionless mass.  Both are tabulated
     from the Lane-Emden solution and the ratio is monotonic, so it inverts.
  2. From xi_max and M, get r0, rho_c, R_BE = xi_max * r0.
  3. T_BE is solved self-consistently against the calibration
         T_BE = a0 + a1 x + a2 x^2 ,  x = log10(N_tot),
     with N_tot = Sigma_emb + N_core the total RADIAL shielding column
     (N_core uses I(xi) = int_0^xi exp(-psi) dxi, a radial path integral).
     The loop iterates T -> xi_max -> rho_c -> N_tot -> T until fixed.
  4. Detectability: the predicted background-subtracted convolved peak is
     computed exactly (see peak_excess below) and the contrast
         contrast = 1 + peak_exc / Sigma_emb
     is flagged against FLOOR.

NOTE ON THE TWO CONSTRAINTS
  * The continuity ("hole") constraint rho_edge >= rho_cloud is NOT imposed.
    It was needed only for an older reduction that convolved before
    subtracting the background; subtracting the constant pedestal at R_BE on
    the ORIGINAL image removes the depression artefact, verified over
    edge/cloud ratios 1.13 -> 0.08.  Models with rho_edge < rho_cloud are
    retained and marked with the `hole` flag (they are the compact, massive
    cores at high column, and resemble the rarefaction of an inside-out
    collapsing core).
  * Detectability IS computed, but as a FLAG (`detect`), not a cut: the full
    regular (i,j,k) lattice is kept so the directory tree and node indices are
    stable.  `run_grid_mc3d.sh` skips detect==0 nodes unless ALL_NODES=1.

PEAK PREDICTION (validated to 0.02% against a directly imaged model)
  The model background is a clean constant-Sigma_emb pedestal, so the
  background-subtracted peak is analytic.  Two subtleties, both of which were
  bugs at some point and are the reason no fudge factor is needed:
    - it is a 2-D convolution, evaluated via the Hankel form
          C(rho) = (1/s^2) Int Sigma(r) exp(-(r^2+rho^2)/2s^2) I0(r rho/s^2) r dr
      NOT a 1-D beam-weighted radial average;
    - the integration domain must cover the BEAM (>= 5 sigma), not merely the
      source radius; and the angular scale is r0, NOT R_out = xi_max * r0.
  The result is then averaged over the central PIX_ARCSEC map pixel, matching
  the observing chain (convolve to BEAM, resample to PIX_ARCSEC).

Usage
    python3 make_bes_grid.py [-o bes_grid_final.txt]

Output columns (23):
    ID i j k SD_emb T_BE M_BE FWHM_pc xi_max contr_rho peak_exc contrast
    detect hole rho_c rho_edge rho_emb r0_AU R_out_AU R_cloud_AU N_tot a_BE stab

`run_grid_mc3d.sh` reads these BY POSITION -- if columns are added or
reordered here, its `while read -r ...` list must be updated in lockstep.
"""

import argparse
import itertools

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import iv

trapz = getattr(np, 'trapz', None) or np.trapezoid

# ----------------------------------------------------------------------
# physical constants (CGS) and pipeline conventions
# ----------------------------------------------------------------------
kB     = 1.380649e-16       # erg/K
G      = 6.674e-8           # cgs
amu    = 1.6605e-24         # g
Msun   = 1.989e33           # g
xAU    = 1.495979e13        # cm
pc     = 3.0857e18          # cm

muH2      = 2.8             # mean molecular weight per H2 (for column density)
mu_p      = 2.33            # mean molecular weight per particle (for c_s)
DIST_PC   = 260.0           # distance to Aquila
BEAM      = 13.5            # arcsec, hires surfdens beam (FWHM)
PIX_ARCSEC = 3.0            # arcsec, common map pixel after resampling
XI_CRIT   = 6.451           # critical BE truncation
FLOOR     = 1.1             # detectability floor on contrast = 1 + peak/SD_emb
                            # (Aquila's measured floor is 1.16; 1.1 is adopted
                            #  so that regions detecting deeper are included)
R_CLOUD_FACTOR = 30.0       # R_cloud = 30 * R_BE

# T_BE(N_tot) calibration:  T = a0 + a1*x + a2*x^2 ,  x = log10(N_tot), clamped
TBE_A2, TBE_A1, TBE_A0 = 1.84189, -86.4699, 1021.39
TBE_XMIN, TBE_XMAX     = 21.487, 23.398

# grid axes
SD_VALUES  = [3e21, 6e21, 1.2e22, 2.4e22, 4.8e22, 9.6e22]     # i
M_MIN, M_MAX, N_M    = 0.008, 100.0, 12                        # j
FW_MIN, FW_MAX, N_FW = 0.003, 0.25, 10                         # k  [pc]

SIG_BEAM = BEAM / 2.3548    # Gaussian sigma of the beam, arcsec


# ----------------------------------------------------------------------
# Lane-Emden solution and the derived tabulations
# ----------------------------------------------------------------------
def _lane_emden(xi_max_tab=62.0, n=600000):
    """Isothermal Lane-Emden: psi'' + 2 psi'/xi = exp(-psi), psi(0)=psi'(0)=0."""
    def rhs(xi, y):
        p, dp = y
        return [dp, (np.exp(-p) - 2.0 * dp / xi) if xi > 1e-9 else 0.0]

    grid = np.linspace(1e-6, xi_max_tab, n)
    sol = solve_ivp(rhs, [grid[0], grid[-1]], [0.0, 0.0],
                    t_eval=grid, rtol=1e-11, atol=1e-13)
    psi  = interp1d(sol.t, sol.y[0], kind='cubic')
    dpsi = interp1d(sol.t, sol.y[1], kind='cubic')
    return psi, dpsi


PSI, DPSI = _lane_emden()
dens = lambda x: np.exp(-PSI(np.clip(x, 1e-6, None)))     # rho/rho_c at xi


def _radial_column_integral():
    """I(xi) = int_0^xi exp(-psi) dxi -- the RADIAL shielding path integral."""
    xg = np.linspace(1e-4, 60.0, 150000)
    ig = np.concatenate([[0.0], np.cumsum(dens(xg[1:]) * np.diff(xg))])
    f = interp1d(xg, ig, kind='cubic')
    return lambda x: float(f(min(x, 60.0)))


I_RADIAL = _radial_column_integral()


def _tabulate_profile_functions():
    """
    Tabulate, versus xi_max:
        phi(xi)  projected (column-density) FWHM in units of r0
        Cm(xi)   = xi^2 dpsi/dxi, the dimensionless mass
        con(xi)  = exp(psi(xi)), the central-to-edge density contrast
    and the monotonic inversion key phi/Cm used to get xi_max from (M, FWHM).
    """
    xitab = np.concatenate([np.linspace(0.5, 12.0, 90),
                            np.linspace(12.4, 60.0, 50)])
    phi, cm, con = [], [], []
    for xm in xitab:
        bb = np.linspace(0.0, xm * 0.9997, 120)
        nb = []
        for b in bb:
            half = xm * xm - b * b
            if half <= 0:
                nb.append(0.0)
                continue
            t = np.linspace(0.0, np.sqrt(half), 90)
            nb.append(2.0 * trapz(dens(np.sqrt(b * b + t * t)), t))
        nb = np.asarray(nb)
        ih = np.where(nb <= nb[0] / 2.0)[0]
        phi.append(2.0 * (bb[ih[0]] if len(ih) else xm))
        cm.append(xm * xm * float(DPSI(xm)))
        con.append(float(np.exp(PSI(xm))))
    phi, cm, con = map(np.asarray, (phi, cm, con))
    key = phi / cm                                    # monotonically decreasing
    return (xitab,
            interp1d(xitab, cm,  kind='cubic', fill_value="extrapolate", bounds_error=False),
            interp1d(xitab, con, kind='cubic', fill_value="extrapolate", bounds_error=False),
            key,
            interp1d(key[::-1], xitab[::-1], bounds_error=False))


XITAB, CM_OF_XI, CON_OF_XI, KEY, XI_OF_KEY = _tabulate_profile_functions()


def T_BE_of_Ntot(n_tot):
    """Self-consistent construction temperature from the total shielding column."""
    x = np.clip(np.log10(n_tot), TBE_XMIN, TBE_XMAX)
    return TBE_A0 + TBE_A1 * x + TBE_A2 * x * x


# ----------------------------------------------------------------------
# convolved, background-subtracted peak (the detectability observable)
# ----------------------------------------------------------------------
_PIX_SAMPLES = list(itertools.product(np.linspace(-PIX_ARCSEC / 2.0,
                                                  PIX_ARCSEC / 2.0, 7), repeat=2))


def peak_excess(rho_c, r0_AU, xi_max):
    """
    Background-subtracted peak column [cm^-2] after convolution to BEAM and
    averaging over one PIX_ARCSEC pixel.

    The profile is the projected excess of the truncated BE sphere above the
    constant cloud pedestal; the pedestal cancels in the subtraction, so only
    the excess is convolved.
    """
    r0_as = r0_AU / DIST_PC
    r0_cm = r0_AU * xAU

    def sigma_at(theta_as):
        b = theta_as / r0_as
        if b >= xi_max:
            return 0.0
        t = np.linspace(0.0, np.sqrt(xi_max * xi_max - b * b), 160)
        return rho_c / (muH2 * amu) * r0_cm * 2.0 * trapz(dens(np.sqrt(b * b + t * t)), t)

    # domain must cover the BEAM, not just the source
    th = np.linspace(0.0, max(5.0 * SIG_BEAM, 1.5 * xi_max * r0_as), 1500)
    prof = np.array([sigma_at(t) for t in th])

    def convolved_at(rho):
        arg = np.clip(th * rho / SIG_BEAM ** 2, 0.0, 700.0)
        return trapz(prof * np.exp(-(th ** 2 + rho ** 2) / (2.0 * SIG_BEAM ** 2))
                     * iv(0, arg) * th, th) / SIG_BEAM ** 2

    return float(np.mean([convolved_at(np.hypot(a, b)) for a, b in _PIX_SAMPLES]))


# ----------------------------------------------------------------------
# one grid node
# ----------------------------------------------------------------------
def build_node(sd_emb, mass, fwhm_pc, max_iter=60, tol=1e-5):
    """Return the derived parameters for one (Sigma_emb, M, FWHM) node, or None
    if no truncated BE sphere with that mass and size exists."""
    fwhm_cm = fwhm_pc * pc
    T = 10.0
    for _ in range(max_iter):
        c_s = np.sqrt(kB * T / (mu_p * amu))
        ratio = fwhm_cm * c_s ** 2 / (mass * Msun * G)
        if not (KEY[-1] <= ratio <= KEY[0]):
            return None                     # (M, FWHM) outside the BE family
        xi_max = float(XI_OF_KEY(ratio))
        r0 = mass * Msun * G / (c_s ** 2 * CM_OF_XI(xi_max))
        rho_c = c_s ** 2 / (4.0 * np.pi * G * r0 ** 2)
        n_core = rho_c * r0 * I_RADIAL(xi_max) / (muH2 * amu)
        T_new = T_BE_of_Ntot(sd_emb + n_core)
        if abs(T_new - T) < tol:
            T = T_new
            break
        T = 0.5 * (T + T_new)               # damped iteration

    c_s = np.sqrt(kB * T / (mu_p * amu))
    ratio = fwhm_cm * c_s ** 2 / (mass * Msun * G)
    if not (KEY[-1] <= ratio <= KEY[0]):
        return None
    xi_max = float(XI_OF_KEY(ratio))
    r0 = mass * Msun * G / (c_s ** 2 * CM_OF_XI(xi_max))
    rho_c = c_s ** 2 / (4.0 * np.pi * G * r0 ** 2)

    r0_AU    = r0 / xAU
    R_out_AU = xi_max * r0_AU
    rho_edge = rho_c / CON_OF_XI(xi_max)
    n_tot    = sd_emb + rho_c * r0 * I_RADIAL(xi_max) / (muH2 * amu)
    R_cloud_AU = R_CLOUD_FACTOR * R_out_AU
    # uniform cloud density that reproduces Sigma_emb through the model box
    rho_emb = sd_emb * muH2 * amu / (2.0 * xAU
                                     * np.sqrt(R_cloud_AU ** 2 - R_out_AU ** 2))

    peak = peak_excess(rho_c, r0_AU, xi_max)
    contrast = 1.0 + peak / sd_emb

    # Bonnor-Ebert stability parameter (a_BE > 2 -> unbound / supercritical)
    a_BE = 1.18 * c_s ** 3 / (G ** 1.5 * np.sqrt(rho_edge)) / Msun / mass

    if xi_max < XI_CRIT - 0.02:
        stab = 'sub'
    elif xi_max < XI_CRIT + 0.05:
        stab = 'crit'
    elif xi_max <= 12.0:
        stab = 'super'
    else:
        stab = 'deep'

    return dict(SD=sd_emb, T=T, M=mass, fw=fwhm_pc, xi=xi_max,
                con=float(CON_OF_XI(xi_max)), rc=rho_c, re=rho_edge, remb=rho_emb,
                r0=r0_AU, Rout=R_out_AU, Rcl=R_cloud_AU, Nt=n_tot, aBE=a_BE,
                peak=peak, contrast=contrast,
                detect=int(contrast >= FLOOR),
                hole=int(rho_edge < rho_emb),
                stab=stab)


# ----------------------------------------------------------------------
# catalogue
# ----------------------------------------------------------------------
HEADER_FIELDS = ["ID", "i", "j", "k", "SD_emb", "T_BE", "M_BE", "FWHM_pc",
                 "xi_max", "contr_rho", "peak_exc", "contrast", "detect", "hole",
                 "rho_c", "rho_edge", "rho_emb", "r0_AU", "R_out_AU",
                 "R_cloud_AU", "N_tot", "a_BE", "stab"]

HEADER_FMT = ("{:>4} {:>3} {:>3} {:>3} {:>9} {:>7} {:>10} {:>8} {:>8} {:>9} "
              "{:>10} {:>9} {:>6} {:>4} {:>11} {:>11} {:>11} {:>10} {:>10} "
              "{:>11} {:>10} {:>7} {:>5}")

ROW_FMT = ("{:>4} {:>3} {:>3} {:>3} {:>9} {:>7.3f} {:>10.4f} {:>8.4f} {:>8.4f} "
           "{:>9.3f} {:>10.3e} {:>9.3f} {:>6d} {:>4d} {:>11.4e} {:>11.4e} "
           "{:>11.4e} {:>10.2f} {:>10.2f} {:>11.2f} {:>10.4e} {:>7.2f} {:>5}")


def build_grid():
    masses = np.logspace(np.log10(M_MIN), np.log10(M_MAX), N_M)
    fwhms  = np.logspace(np.log10(FW_MIN), np.log10(FW_MAX), N_FW)
    rows = []
    for i, sd in enumerate(SD_VALUES, 1):
        for j, m in enumerate(masses, 1):
            for k, fw in enumerate(fwhms, 1):
                node = build_node(sd, m, fw)
                if node is not None:
                    node.update(i=i, j=j, k=k)
                    rows.append(node)
    return rows


def write_catalog(rows, path):
    lines = [HEADER_FMT.format(*HEADER_FIELDS)]
    for n, r in enumerate(rows, 1):
        lines.append(ROW_FMT.format(
            n, r['i'], r['j'], r['k'], "%.2e" % r['SD'], r['T'], r['M'], r['fw'],
            r['xi'], r['con'], r['peak'], r['contrast'], r['detect'], r['hole'],
            r['rc'], r['re'], r['remb'], r['r0'], r['Rout'], r['Rcl'],
            r['Nt'], r['aBE'], r['stab']))
    with open(path, 'w') as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-o', '--out', default='bes_grid_final.txt',
                    help='output catalogue (default: bes_grid_final.txt)')
    args = ap.parse_args()

    rows = build_grid()
    write_catalog(rows, args.out)

    n_det  = sum(r['detect'] for r in rows)
    n_hole = sum(r['hole'] for r in rows)
    print("wrote %s : %d models  (detectable %d, hole %d)"
          % (args.out, len(rows), n_det, n_hole))
    print("  beam %.1f\"  pixel %.1f\"  contrast floor %.2f" % (BEAM, PIX_ARCSEC, FLOOR))
    for sd in SD_VALUES:
        sel = [r for r in rows if r['SD'] == sd]
        print("  SD_emb=%.1e : %3d models, %3d detectable, %3d hole, M<=%.1f"
              % (sd, len(sel), sum(r['detect'] for r in sel),
                 sum(r['hole'] for r in sel), max(r['M'] for r in sel)))


if __name__ == '__main__':
    main()
