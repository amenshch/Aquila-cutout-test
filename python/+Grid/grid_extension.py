#!/usr/bin/env python
"""grid_extension.py -- define the additional Bonnor-Ebert models needed to
close the gaps found in the present grid, and report how many nodes each
extension costs.

The grid of make_bes_grid.py is defined on three observable axes,

    SD_VALUES              embedding column density of the cloud   (cm^-2)
    M_MIN, M_MAX, N_M      core mass, logarithmically spaced       (M_sun)
    FW_MIN, FW_MAX, N_FW   core FWHM, logarithmically spaced       (pc)

everything else (xi_max, rho_c, r0, R_BE, T_BE) being derived per node.  The
extensions below are therefore stated as changes to those three definitions.

THE THREE GAPS, AND WHAT CLOSES EACH

1.  The corrector cannot reach faint backgrounds.  Of 8106 sources in the
    seven clouds, 3820 (47%) fall outside the convex hull of the training set,
    and 3287 of those fail on the background column alone, lying below the
    lowest embedding rung of 3.0e21 cm^-2.  No other axis contributes more
    than 99 sources.  Three further rungs continuing the factor-of-two ladder
    downward, to 3.75e20 cm^-2, cover the column densities at which the
    nearby, diffuse clouds actually have their cores: Scorpius currently
    retains 24 of 187 sources and Ophiuchus 232 of 965.

2.  The grid has no faint compact models.  Once mass and FWHM are fixed the
    contrast follows, so a low-contrast core of a given size requires a low
    mass.  At the 6.0e21 rung the faintest model between 26 and 35 arcsec has
    contrast 1.87 where real Aquila sources of that size have 1.26, and at the
    3.0e21 rung no size class reaches the observed contrast at all.  The mass
    ladder starts at 0.008 M_sun; extending it down to 0.001 M_sun adds the
    three steps that bring the contrast down by the required factor of about
    three.

3.  The ladders are coarse for interpolation.  Successive masses differ by a
    factor of 2.36 and successive sizes by 1.63.  The correction factor varies
    smoothly, so this is not fatal, but the leave-one-out residual is measured
    on exactly these steps and halving them is the direct way to reduce the
    interpolation component of the error.

Usage:
    python grid_extension.py
"""
import numpy as np

# present definitions, copied from make_bes_grid.py
SD_NOW = [3e21, 6e21, 1.2e22, 2.4e22, 4.8e22, 9.6e22]
M_NOW = (0.008, 100.0, 12)
FW_NOW = (0.003, 0.25, 10)
DIST_PC = 260.0


def ladder(lo, hi, n):
    return np.logspace(np.log10(lo), np.log10(hi), n)


def describe(sd, m, fw, label):
    M = ladder(*m)
    F = ladder(*fw)
    n = len(sd) * len(M) * len(F)
    print('  %-46s %6d nodes  (%d x %d x %d)'
          % (label, n, len(sd), len(M), len(F)))
    return n


def main():
    print('Present grid')
    n0 = describe(SD_NOW, M_NOW, FW_NOW, 'as built')
    print('    embedding columns (cm^-2): %s'
          % ' '.join('%.3g' % v for v in SD_NOW))
    print('    masses (M_sun): %.3g to %.3g in %d steps, factor %.2f each'
          % (M_NOW[0], M_NOW[1], M_NOW[2],
             (M_NOW[1] / M_NOW[0]) ** (1.0 / (M_NOW[2] - 1))))
    print('    FWHM (pc): %.3g to %.3g in %d steps, factor %.2f each'
          % (FW_NOW[0], FW_NOW[1], FW_NOW[2],
             (FW_NOW[1] / FW_NOW[0]) ** (1.0 / (FW_NOW[2] - 1))))
    print('    FWHM in arcsec at %.0f pc: %s'
          % (DIST_PC, ' '.join('%.0f' % (v / DIST_PC * 206265)
                               for v in ladder(*FW_NOW))))
    print('    (518 of these %d nodes were actually produced; the rest fail\n'
          '     the detectability floor or leave the tabulated xi_max range)'
          % n0)

    SD_EXT = [3.75e20, 7.5e20, 1.5e21] + SD_NOW
    M_EXT = (0.001, 100.0, 15)          # same factor 2.36, three steps lower
    M_FINE = (0.001, 100.0, 29)         # factor 1.54
    FW_FINE = (0.003, 0.25, 19)         # factor 1.28

    print('\nProposed extensions, cumulative')
    describe(SD_EXT, M_NOW, FW_NOW,
             '1. three lower embedding rungs')
    describe(SD_EXT, M_EXT, FW_NOW,
             '1 + 2. lower rungs and masses to 0.001 M_sun')
    describe(SD_EXT, M_FINE, FW_NOW,
             '1 + 2 + finer mass ladder (factor 1.54)')
    describe(SD_EXT, M_FINE, FW_FINE,
             '1 + 2 + finer mass and size ladders')

    print('\nThe new embedding rungs (cm^-2): %s'
          % ' '.join('%.3g' % v for v in SD_EXT[:3]))
    print('The new low-mass steps (M_sun): %s'
          % ' '.join('%.4g' % v for v in ladder(*M_EXT)[:4]))

    print("""
Changes to make_bes_grid.py
---------------------------
Minimum, closing gaps 1 and 2:

    SD_VALUES  = [3.75e20, 7.5e20, 1.5e21,
                  3e21, 6e21, 1.2e22, 2.4e22, 4.8e22, 9.6e22]   # i
    M_MIN, M_MAX, N_M    = 0.001, 100.0, 15                     # j
    FW_MIN, FW_MAX, N_FW = 0.003, 0.25, 10                      # k  [pc]

Adding gap 3 as well, halving both ladder steps:

    M_MIN, M_MAX, N_M    = 0.001, 100.0, 29                     # j
    FW_MIN, FW_MAX, N_FW = 0.003, 0.25, 19                      # k  [pc]

Two cautions.  The indices i, j and k are written into every node tag, so the
new rungs shift the numbering of every existing node: either prepend the new
rungs and accept that i changes meaning, or append them and keep the ladder
unsorted.  Whichever is chosen, the recovery tables and the injection truth
tables refer to nodes by tag and must be regenerated together with the
catalogue.

The detectability floor FLOOR = 1.1 will reject a large part of the new
low-mass nodes at the low embedding rungs, since a faint core on a faint
background has low contrast by construction.  That is the intended behaviour,
but it means the produced node count will fall well short of the products
quoted above; the numbers here are upper bounds on the computation.
""")


if __name__ == '__main__':
    main()
