#!/usr/bin/env python3
"""make_cds_tables.py -- the machine-readable tables of the paper for the CDS, with their ReadMe.

Three tables are written, in the fixed-width format of the CDS, together with a ReadMe that describes them byte by byte:

    table1.dat  the correction factor f for G_0 = 10 (the adopted relation of Fig. 4), against the fitted surface density
                Sigma_F, with the true surface density, the fitted temperature, and f held at its maximum as in Eq. (9);
    table2.dat  the field factor g(Sigma_F, G_0) for the nine fields of the grid (Fig. 5), the table GTAB of hicorrect2.py;
    table3.dat  the constants of the temperature relation T(Sigma_A) of Eq. (5) for each of the nine fields (Table C.2), and,
                in the first row, the adopted relation for G_0 = 10 fitted to the main grid and the filaments.

Table 1 is computed with hicorrect2.py for the wavebands 160, 250, 350, and 500 um with the weights of hires (relative
calibration errors of 20% at 160 um and 10% at longer wavelengths), on a finer grid than the code uses (NTABLE). Table 3 reads
the refitted constants from laws_refit.txt (refit.py).

Usage:

    python3 make_cds_tables.py [<laws_refit.txt>] [<output directory>]

Requires hicorrect2.py in the same directory.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hicorrect2 as H

BANDS = [160.0, 250.0, 350.0, 500.0]
NFINE = 400                      # points of Table 1, from 10^20 to 10^24 cm^-2 in the true surface density


def table1(out):
    """The adopted relation for G_0 = 10."""
    H.NTABLE = NFINE
    fitted, factor, tfit = H.build_relation_t(BANDS)
    sigma = fitted * factor
    held = np.maximum.accumulate(factor)
    rows = []
    for s, sf, f, fh, t in zip(sigma, fitted, factor, held, tfit):
        rows.append("%6.3f %6.3f %5.3f %5.3f %5.2f" % (np.log10(sf), np.log10(s), f, fh, t))
    open(out, "w").write("\n".join(rows) + "\n")
    return len(rows)


def table2(out):
    """The field factor g for the nine fields of the grid."""
    xs = H.GX0 + H.GDX * np.arange(H.GTAB.shape[0])
    rows = []
    for x, g in zip(xs, H.GTAB):
        rows.append("%5.2f" % x + "".join(" %5.3f" % v for v in g))
    open(out, "w").write("\n".join(rows) + "\n")
    return len(rows)


def table3(laws, out):
    """The constants of Eq. (5) for each field, after the adopted relation for G_0 = 10."""
    rows = ["%6.2f %6.3f %6.3f %5.3f %5.3f %4.2f %1d"
            % (10.0, H.TFLOOR, H.TSURF, H.SSTAR / 1e21, H.GAMMA, 2.4, 0)]
    for ln in open(laws):
        if ln.startswith("#") or not ln.strip():
            continue
        g0, tmin, tmax, sstar, gam, scat = [float(v) for v in ln.split()]
        rows.append("%6.2f %6.3f %6.3f %5.3f %5.3f %4.2f %1d" % (g0, tmin, tmax, sstar / 1e21, gam, scat, 1))
    open(out, "w").write("\n".join(rows) + "\n")
    return len(rows)


README = """J/A+A/???/???      Temperatures and surface densities of externally heated clouds  (Men'shchikov, 2026)
================================================================================
Correcting underestimated surface densities of molecular clouds derived from
Herschel images by single-temperature fitting
    Men'shchikov A.
   <Astron. Astrophys. ???, ??? (2026)>
   =2026A&A...???..???M
================================================================================
ADC_Keywords: Interstellar medium ; Molecular clouds ; Infrared sources
Keywords: radiative transfer - dust, extinction - ISM: clouds - stars: formation

Abstract:
    The tables give the correction factor that restores the true surface
    density of an externally heated structure from the surface density
    returned by a single-temperature modified blackbody fit to Herschel
    images at 160-500um, its dependence on the strength of the interstellar
    radiation field, and the constants of the relation between the dust
    temperature and the effective attenuating column that underlies them.

Description:
    All quantities are computed for the dust model of the paper, the
    interstellar radiation field of Black (1994) scaled by G_0 (in Habing
    units), an isotropic external field, and no internal heating. The
    surface densities are those of molecular hydrogen, in cm-2.

    Table 1 gives the adopted relation for G_0=10 (Fig. 4 of the paper), for
    images fitted to 160, 250, 350, and 500um with relative calibration
    errors of 20% at 160um and 10% at longer wavelengths (the weights of
    hires). Below Sigma_F=4x10^21cm-2 it includes the calibration on a
    radiative transfer model of a diffuse cloud. The column fhold is f held
    at its maximum, max(f(Sigma') for Sigma'<=Sigma_F), which the correction
    applies.

    Table 2 gives the field factor g(Sigma_F,G_0) for the nine fields of the
    grid (Fig. 5 of the paper). For another field, log(g) is interpolated
    linearly in log(G_0) between the columns; above the last row, g is held
    at its values there. The correction factor for the field G_0 is
    g(Sigma_F,G_0) x fhold(Sigma_F) (Eq. 9 of the paper): f is held at its
    maximum first and then multiplied by g, which is not held.

    Table 3 gives the constants of the relation
        T(Sigma_A) = Tmin + (Tmax-Tmin) (1+Sigma_A/Sigma_*)^(-gamma)
    (Eq. 5 of the paper) between the dust temperature T and the effective
    attenuating column Sigma_A (Eq. 4 of the paper). The first row (Flag=0)
    is the adopted relation for G_0=10, fitted to 2086 models of embedded
    Bonnor-Ebert spheres and 7 cylindrical filaments; the others (Flag=1)
    were fitted separately to the grids of models in each field (Table C.2).

    The tables were computed with hicorrect2.py and make_cds_tables.py, which
    are available with the getsf distribution.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
table1.dat        31     {n1}   Correction factor f for G_0=10
table2.dat        59     {n2}   Field factor g(Sigma_F,G_0) for nine fields
table3.dat        39     {n3}   Constants of the temperature relation (Eq. 5)
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table1.dat
--------------------------------------------------------------------------------
   Bytes Format Units  Label   Explanations
--------------------------------------------------------------------------------
   1-  6  F6.3  [cm-2] logSigF Log of the fitted surface density Sigma_F
   8- 13  F6.3  [cm-2] logSig  Log of the true surface density Sigma
  15- 19  F5.3  ---    f       Correction factor f=Sigma/Sigma_F
  21- 25  F5.3  ---    fhold   Factor f held at its maximum (Eq. 9)
  27- 31  F5.2  K      TF      Temperature returned by the fit
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table2.dat
--------------------------------------------------------------------------------
   Bytes Format Units  Label   Explanations
--------------------------------------------------------------------------------
   1-  5  F5.2  [cm-2] logSigF Log of the fitted surface density Sigma_F
   7- 11  F5.3  ---    g0.98   Field factor g for G_0=0.98
  13- 17  F5.3  ---    g2.12   Field factor g for G_0=2.12
  19- 23  F5.3  ---    g4.61   Field factor g for G_0=4.61
  25- 29  F5.3  ---    g10.0   Field factor g for G_0=10.0 (unity)
  31- 35  F5.3  ---    g21.7   Field factor g for G_0=21.7
  37- 41  F5.3  ---    g47.1   Field factor g for G_0=47.1
  43- 47  F5.3  ---    g102.4  Field factor g for G_0=102.4
  49- 53  F5.3  ---    g222.2  Field factor g for G_0=222.2
  55- 59  F5.3  ---    g482.5  Field factor g for G_0=482.5
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table3.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label   Explanations
--------------------------------------------------------------------------------
   1-  6  F6.2  ---       G0      Strength of the radiation field (Habing units)
   8- 13  F6.3  K         Tmin    Temperature of fully shielded dust
  15- 20  F6.3  K         Tmax    Temperature of unattenuated dust
  22- 26  F5.3  10+21cm-2 Sig*    Attenuating column of the fall, Sigma_*
  28- 32  F5.3  ---       gamma   Index of the fall
  34- 37  F4.2  %         Scat    Scatter of the model temperatures about
                                   the relation
      39  I1    ---       Flag    [0/1] 0: adopted relation for G_0=10;
                                   1: fit to the grid of that field
--------------------------------------------------------------------------------

Acknowledgements:
    Alexander Men'shchikov, alexander.menshchikov(at)cea.fr

================================================================================
(End)
"""


def main():
    laws = sys.argv[1] if len(sys.argv) > 1 else "laws_refit.txt"
    out = sys.argv[2] if len(sys.argv) > 2 else "."
    os.makedirs(out, exist_ok=True)
    n1 = table1(os.path.join(out, "table1.dat"))
    n2 = table2(os.path.join(out, "table2.dat"))
    n3 = table3(laws, os.path.join(out, "table3.dat"))
    open(os.path.join(out, "ReadMe"), "w").write(README.replace("{n1}", "%d" % n1).replace("{n2}", "%d" % n2)
                                                 .replace("{n3}", "%d" % n3))
    print(" written table1.dat (%d rows), table2.dat (%d rows), table3.dat (%d rows) and ReadMe to %s" % (n1, n2, n3, out))


if __name__ == "__main__":
    main()
