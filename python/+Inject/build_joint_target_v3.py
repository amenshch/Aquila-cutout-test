#!/usr/bin/env python
"""build_joint_target_v3.py -- write the table of real source properties that
the injected sample is required to imitate.

One row per real source, four columns:

    M_reported   reported SED mass of the source                    (M_sun)
    Sigma_cloud  background column density under the source,
                 the getsf quantity PEAK^BGF                        (cm^-2)
    contrast     PEAK^SBF / PEAK^BGF = (source + background)
                 / background, always greater than unity            (-)
    AFWHM_SD     major axis at half maximum measured in the
                 high-resolution surface-density image at 13.5
                 arcsec, waveband 02 of the getsf extraction        (arcsec)

The size column is new in version 3.  It is the observed size of the source,
and it is the quantity that has to be compared with the model column
FWHMSDbs of the Bonnor-Ebert grid.  The truncation radius R_BE of a model is
not an observable and is far larger than the size the model presents on the
sky, so selecting and separating models by R_BE, as run_inject_v2.py did,
both discarded usable models and gave the models that survived exclusion
zones much larger than their apparent extent.

Usage:
    python build_joint_target_v3.py [output_file]
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mass_correction_pipeline_5 as mc          # noqa: E402
from getsf_columns import GetsfTable             # noqa: E402

# The extraction of the whole Aquila field by G. Zhang.  Its sources define
# what a real core looks like at a given background column.  Split by column,
# the whole-field population and the population inside the injection sub-field
# have the same median contrast and the same median reported mass in every
# octave, so the larger sample may be used without importing a bias; only the
# marginal distribution of the column differs between the two, and that
# marginal is set by the map the cores are placed in, not by this table.
CATALOGS = ['Aquila-Guoyin/Aquila.sw.sources.ok.cat',
            'Aquila-Guoyin/Aquila.sw.sources.ok.add.cat']
MAIN_CATALOG = CATALOGS[0]
SIZE_COLUMN = 'AFWHM02'        # waveband 02 = hi.surface.density.r13p5
DISTANCE_PC = 260.0
OUT_DEFAULT = 'joint_target_Aquila_v3.txt'


def main(out_path=OUT_DEFAULT):
    mc.VERBOSE = 0
    mc.load_getsf._dist = DISTANCE_PC
    g = mc.load_getsf(CATALOGS)

    size = np.asarray(GetsfTable(MAIN_CATALOG).col(SIZE_COLUMN), float)
    if len(size) != len(g['mass']):
        raise ValueError('size column has %d rows, catalogue has %d'
                         % (len(size), len(g['mass'])))

    ok = (mc.source_mask(g)
          & (g['Nbg'] > 0) & (g['peak'] > 0) & (g['conc'] > 0)
          & (g['mass'] > 0)
          & np.isfinite(size) & (size > 0) & (size < 1.0e3))

    M = g['mass'][ok]
    S = g['Nbg'][ok]
    C = 1.0 + g['peak'][ok] / g['Nbg'][ok]
    A = size[ok]

    hdr = []
    hdr.append('# Real Aquila sources used as the target distribution for the')
    hdr.append('# injection test.  Written by build_joint_target_v3.py from')
    for c in CATALOGS:
        hdr.append('#   %s' % c)
    hdr.append('# Assumed distance %.0f pc.  %d sources.' % (DISTANCE_PC, len(M)))
    hdr.append('#')
    hdr.append('# Columns')
    hdr.append('#   1 M_reported   reported SED mass of the source          (M_sun)')
    hdr.append('#   2 Sigma_cloud  background column density, PEAK^BGF      (cm^-2)')
    hdr.append('#   3 contrast     PEAK^SBF / PEAK^BGF, exceeds unity       (-)')
    hdr.append('#   4 AFWHM_SD     major axis at half maximum in the        (arcsec)')
    hdr.append('#                  surface-density image at 13.5 arcsec')
    hdr.append('#')
    for nm, v, f in (('M_reported ', M, '%10.4g'), ('Sigma_cloud', S, '%10.4g'),
                     ('contrast   ', C, '%10.4g'), ('AFWHM_SD   ', A, '%10.4g')):
        q = np.percentile(v, [5, 25, 50, 75, 95])
        hdr.append(('# %s  5th ' + f + '  25th ' + f + '  median ' + f
                    + '  75th ' + f + '  95th ' + f) % ((nm,) + tuple(q)))
    hdr.append('#')
    hdr.append('# %14s %14s %12s %12s'
               % ('M_reported', 'Sigma_cloud', 'contrast', 'AFWHM_SD'))

    with open(out_path, 'w') as fh:
        fh.write('\n'.join(hdr) + '\n')
        for m, s, c, a in zip(M, S, C, A):
            fh.write('  %14.6e %14.6e %12.5f %12.3f\n' % (m, s, c, a))

    print('wrote %s with %d real sources' % (out_path, len(M)))
    for nm, v in (('M_reported (M_sun)', M), ('Sigma_cloud (cm^-2)', S),
                  ('contrast (-)', C), ('AFWHM_SD (arcsec)', A)):
        q = np.percentile(v, [5, 50, 95])
        print('  %-22s 5th %10.4g   median %10.4g   95th %10.4g'
              % (nm, q[0], q[1], q[2]))
    r = np.corrcoef(np.log10(S), np.log10(C - 1.0))[0, 1]
    print('  correlation of log10(Sigma_cloud) with log10(contrast - 1): %+.3f' % r)
    r = np.corrcoef(np.log10(M), np.log10(S))[0, 1]
    print('  correlation of log10(M_reported) with log10(Sigma_cloud):   %+.3f' % r)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else OUT_DEFAULT)
