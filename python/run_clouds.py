#!/usr/bin/env python3
"""
run_clouds.py -- multi-cloud driver for mass_correction_pipeline.

Loops over a set of clouds, runs the single-cloud pipeline on each (writing one
per-source corrected catalogue per cloud), and writes a run-summary table (the
seven-cloud Table 1) in the same aligned getsf style.

Usage:
    python run_clouds.py [--outdir DIR]

Edit CLOUDS below for your own fields.  Each entry is (name, distance_pc); the
getsf input for a cloud is resolved from its <name>-Guoyin/ directory
(concatenated catalogue if present, else the raw .cat + .add.cat pair).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mass_correction_pipeline as mc

CLOUDS = [('Scorpius', 130), ('Ophiuchus', 144), ('Aquila', 260),
          ('OrionA', 432), ('California', 470), ('CygnusX', 1150),
          ('W3W4W5', 1700)]

SUMCOLS = [('FIELD', '%-11s', '', 'Cloud / field name'),
           ('DIST', '%6.0f', 'pc', 'Adopted distance'),
           ('N_CORR', '%7d', '', 'Number of corrected QUALITY==ok sources in the model hull'),
           ('R_FOOT', '%7.2f', '', 'Median footprint-to-size ratio over ok sources'),
           ('C_M', '%7.2f', '', 'Median total correction factor M_BE/M_reported'),
           ('F_FLUX', '%7.2f', '', 'Median flux-recovery mass factor'),
           ('F_BG', '%7.3f', '', 'Median background-subtraction mass factor'),
           ('F_TEMP', '%7.3f', '', 'Median temperature-bias mass factor'),
           ('C_M_0.1_2', '%9.2f', '', 'Median C_M over reported mass 0.1-2 Msun'),
           ('C_M_MASS', '%9.2f', '', 'Median mass-indexed factor (distance-dependent; diagnostic)')]


def write_summary(rows, path):
    bar = '!' + '_' * 90
    L = [bar, '!', '! MASS-CORRECTION RUN SUMMARY  (one line per field)',
         '! R_foot = %s   opacity const K = %.4e   SED bands = %s'
         % (mc.RFOOT_MODE, mc.OPACITY_K, ','.join('%d' % w for w in mc.SED_WAVES)),
         '!', bar, '!', '! TABULATED QUANTITIES:', '!']
    for k, (nm, _, unit, desc) in enumerate(SUMCOLS, 1):
        L.append('! %d %-11s %-4s %s' % (k, nm, unit, desc))
    L += [bar, '!']
    widths = [max(len(nm), mc._fmt_width(fmt)) for nm, fmt, _, _ in SUMCOLS]
    L.append('! ' + ' '.join('%*s' % (w, nm) for (nm, *_), w in zip(SUMCOLS, widths)))
    for r in rows:
        cells = []
        for (nm, fmt, _, _), w in zip(SUMCOLS, widths):
            v = r[nm]
            s = fmt % v if not (isinstance(v, float) and not np.isfinite(v)) else 'nan'
            cells.append('%*s' % (w, s.strip()) if fmt.strip().startswith('%-') else '%*s' % (w, s))
        L.append('  ' + ' '.join(cells))
    open(path, 'w').write('\n'.join(L) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='corrected')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rows = []
    print("%-11s %6s %7s %7s %7s %7s %7s %7s %9s %9s"
          % ('field', 'dist', 'N_corr', 'R_foot', 'C_M', 'F_flux', 'F_bg', 'F_temp',
             'C_M0.1-2', 'C_M^mass'))
    for name, dist in CLOUDS:
        inp = mc.cloud_input(name)
        table, s = mc.correct_cloud(inp, dist)
        out = os.path.join(a.outdir, '%s.corrected.cat' % name)
        n = mc.write_catalog(table, out, name, dist)
        rows.append(dict(FIELD=name, DIST=dist, N_CORR=n,
                         R_FOOT=s['Rfoot_ok'], C_M=s['C_M'], F_FLUX=s['F_flux'],
                         F_BG=s['F_bg'], F_TEMP=s['F_temp'],
                         **{'C_M_0.1_2': s['C_M_0p1_2'], 'C_M_MASS': s['C_Mmass']}))
        print("%-11s %6.0f %7d %7.2f %7.2f %7.2f %7.3f %7.3f %9.2f %9.2f"
              % (name, dist, n, s['Rfoot_ok'], s['C_M'], s['F_flux'], s['F_bg'],
                 s['F_temp'], s['C_M_0p1_2'], s['C_Mmass']))
    write_summary(rows, os.path.join(a.outdir, 'run_summary.cat'))
    print("\nwrote %d catalogues + run_summary.cat to %s/" % (len(rows), a.outdir))


if __name__ == '__main__':
    main()
