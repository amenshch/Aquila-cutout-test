#!/usr/bin/env python
"""flag_hole_models.py -- identify the grid nodes whose Bonnor-Ebert outer
boundary density falls below the density of the cloud they are embedded in.

WHY
===
In the radiative-transfer model the sphere occupies r < R_BE and the cloud
begins exactly at R_BE, so along a line of sight the sphere REPLACES cloud
material rather than adding to it.  At an impact parameter just inside R_BE the
chord samples only the outermost sphere material, at rho_edge.  If
rho_edge < rho_emb that chord carries less column than the neighbouring
line of sight just outside R_BE, which passes through cloud alone, and the
model image dips BELOW the surrounding cloud in an annulus inside R_BE.

The dip is a genuine column deficit, present before any convolution, and no
constant pedestal subtraction can remove it: the far field and the annulus
cannot both be brought to zero by subtracting one number.  Subtracting the
value at R_BE, which is the cloud level, sets the far field to zero and leaves
the annulus negative, which is what appears as a depression in the
source-only images and, when the model is injected, removes flux from the map
around the core and from any neighbour that reaches into that annulus.

make_bes_grid.py computes exactly this condition and records it as the `hole`
flag, but deliberately does not impose it as a cut, and the flag is written to
the grid builder's own table rather than to
bes_model_grid_final2_catalog, which run_inject_v3.py reads.  This script
recomputes the flag from the production catalogue so that the affected nodes
can be excluded from injection.

HOW
===
Per node the catalogue gives the central density rho_c = rho_BE, the
truncation radius R_BE, and the temperature T_BE.  From these:

    c_s^2   = kB * T_BE / (mu_p * amu)          isothermal sound speed
    r0      = sqrt(c_s^2 / (4 pi G rho_c))      Bonnor-Ebert scale radius
    xi_max  = R_BE / r0                         dimensionless truncation
    rho_edge = rho_c / exp(psi(xi_max))         from the Lane-Emden solution

and the embedding density is the uniform cloud density that reproduces the
embedding column through the model box, as make_bes_grid.py defines it,

    rho_emb = SD_emb * mu_H2 * amu
              / (2 * sqrt(R_cloud^2 - R_BE^2)),   R_cloud = 30 * R_BE.

A node is flagged when rho_edge < rho_emb.

Usage:
    python flag_hole_models.py [catalogue] [-o hole_flags.txt]
"""
import argparse
import sys

import numpy as np
from scipy.integrate import solve_ivp

kB = 1.380649e-16
G = 6.674e-8
amu = 1.6605e-24
xAU = 1.495979e13
PC = 3.0856776e18
muH2 = 2.8              # per H2, for column density
mu_p = 2.33             # per free particle, for the sound speed
R_CLOUD_FACTOR = 30.0


def lane_emden(xi_max=200.0, n=400000):
    """psi(xi) for the isothermal sphere: psi'' + 2 psi'/xi = exp(-psi)."""
    def rhs(x, y):
        return [y[1], np.exp(-y[0]) - 2.0 * y[1] / x]
    x0 = 1e-6
    y0 = [x0 ** 2 / 6.0, x0 / 3.0]
    grid = np.linspace(x0, xi_max, n)
    sol = solve_ivp(rhs, (x0, xi_max), y0, t_eval=grid, rtol=1e-10, atol=1e-12)
    return sol.t, sol.y[0]


def load(path):
    hdr = None
    for ln in open(path):
        if ln.lstrip()[:1] == '#':
            w = ln.lstrip('# ').split()
            if 'ICSDbs' in w and 'M_BE' in w:
                hdr = w
                break
    if hdr is None:
        raise SystemExit('no column header found in %s' % path)
    rows = [l.split() for l in open(path)
            if l.strip() and l.lstrip()[:1] != '#']
    col = lambda n: np.array([float(r[hdr.index(n)]) for r in rows])
    tag = np.array(['i%02dj%02dk%02d' % (int(a), int(b), int(c))
                    for a, b, c in zip(col('i'), col('j'), col('k'))])
    return hdr, tag, col


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('catalogue', nargs='?',
                    default='cats/bes_model_grid_final2_catalog')
    ap.add_argument('-o', default='hole_flags.txt')
    a = ap.parse_args()

    hdr, tag, col = load(a.catalogue)
    rho_c = col('rho_BE')
    T = col('T_BE')
    R_pc = col('R_BE_pc')
    sd = col('SD_emb')
    fw = col('FWHMSDbs')
    mbe = col('M_BE')

    xi_t, psi_t = lane_emden()
    cs2 = kB * T / (mu_p * amu)
    r0 = np.sqrt(cs2 / (4.0 * np.pi * G * rho_c))          # cm
    R_cm = R_pc * PC
    xi_max = R_cm / r0
    xi_max = np.clip(xi_max, xi_t[0], xi_t[-1])
    rho_edge = rho_c / np.exp(np.interp(xi_max, xi_t, psi_t))

    R_out_AU = R_cm / xAU
    R_cloud_AU = R_CLOUD_FACTOR * R_out_AU
    rho_emb = (sd * muH2 * amu
               / (2.0 * xAU * np.sqrt(R_cloud_AU ** 2 - R_out_AU ** 2)))

    ratio = rho_edge / rho_emb
    hole = ratio < 1.0
    print('nodes in %s: %d' % (a.catalogue, len(tag)))
    print('edge-to-cloud density ratio: 5th %.3f, median %.3f, 95th %.3f'
          % tuple(np.percentile(ratio, [5, 50, 95])))
    print('nodes with rho_edge < rho_emb (a depression, "hole"): %d (%.0f%%)'
          % (hole.sum(), 100 * hole.mean()))
    sel = (fw >= 15) & (fw <= 80)
    print('  among the %d nodes usable for injection '
          '(FWHMSDbs 15-80 arcsec): %d (%.0f%%)'
          % (sel.sum(), int((hole & sel).sum()),
             100 * hole[sel].mean()))
    print('\nby embedding column:')
    print('  %12s %8s %8s %10s' % ('SD_emb', 'nodes', 'holes', 'fraction'))
    for r in np.unique(sd):
        m = sd == r
        print('  %12.2e %8d %8d %9.0f%%'
              % (r, m.sum(), int(hole[m].sum()), 100 * hole[m].mean()))
    for t in ('i05j15k06',):
        k = np.where(tag == t)[0]
        if len(k):
            k = k[0]
            print('\n%s: xi_max %.2f, rho_edge %.3e, rho_emb %.3e, '
                  'ratio %.3f -> %s'
                  % (t, xi_max[k], rho_edge[k], rho_emb[k], ratio[k],
                     'HOLE' if hole[k] else 'ok'))

    with open(a.o, 'w') as f:
        f.write('# Depression ("hole") flag per node of %s\n' % a.catalogue)
        f.write('# A node is flagged when the Bonnor-Ebert outer boundary\n')
        f.write('# density falls below the density of the embedding cloud, so\n')
        f.write('# that the model removes column from the map in an annulus\n')
        f.write('# inside R_BE instead of adding it.  See the module docstring\n')
        f.write('# of flag_hole_models.py.\n#\n')
        f.write('# %-11s %10s %12s %12s %10s %10s %6s\n'
                % ('node', 'xi_max', 'rho_edge', 'rho_emb', 'ratio',
                   'FWHMSDbs', 'hole'))
        f.write('# %-11s %10s %12s %12s %10s %10s %6s\n'
                % ('', '(-)', '(g/cm3)', '(g/cm3)', '(-)', '(arcsec)', '(0/1)'))
        for i in range(len(tag)):
            f.write('  %-11s %10.3f %12.4e %12.4e %10.4f %10.2f %6d\n'
                    % (tag[i], xi_max[i], rho_edge[i], rho_emb[i], ratio[i],
                       fw[i], int(hole[i])))
    print('\nwrote %s' % a.o)


if __name__ == '__main__':
    main()
