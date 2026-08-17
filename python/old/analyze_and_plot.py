#!/usr/bin/env python3
"""
analyze_and_plot.py -- injection-recovery analysis and plots for one cSD run.

Matches injected models (truth table) to a getsf extraction catalog by pixel
position, applies the RT mass correction (ObservableCorrector, direct mode,
keyed on observed M_SED, local cloud Sigma, and observed FWHM_250), and writes
a 3-page PDF:
    1. mass bias        M_SED / M_BE   vs M_BE
    2. corrected acc.   M_corr / M_BE  vs M_BE
    3. applied factor   M_corr / M_SED vs M_BE

Usage:
    python analyze_and_plot.py TRUTH.txt CATALOG.cat OUT.pdf [--label cSD_03] \
        [--scale x1] [--grid bes_model_params_catalog] [--match-rad 4.5]

analyze_and_plot.py "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260616_RT_BES_radmc3d/Aquila/+Images/260704_cSD_01/inj_all24_s1111_x0.25_truth.txt"
                    "/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260616_RT_BES_radmc3d/Aquila/260704_cSD_01/runs/results/+sources/+catalogs/Aquila.s.sources.ok.cat=Aquila.s.sources.ok.add.cat=thin.Aquila.s.sources.ok.00.cat"
                    cSD_01_recovery_plots.pdf --grid "../bes_model_params_catalog" --label cSD_01 --scale x0.25 --match-rad 4.5

Requires fit_mass.py (with load_catalog + ObservableCorrector) on the path.
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from fit_mass import load_catalog, ObservableCorrector

# --- getsf catalog column layout -------------------------------------------
# ok.cat: 11 global cols + 6 bands * 17 per-band cols.
# AFWHM is the 10th per-band col (0-based 9); 250 um is band index 3.
AFWHM_250_IDX = 11 + 3 * 17 + 9          # = 71
THIN_NCOL     = 24                        # thin block width (NO..QUALITY)


def parse_catalog(path):
    """Return list of dicts: x, y (pixels), M_SED (Msun), afwhm (250um, "), qual."""
    out = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith('#') or s.startswith('!'):
            continue
        p = line.split()
        # locate the thin block by its trailing quality token
        for qi, tok in enumerate(p):
            if tok in ('ok', 'bad'):
                ts = qi - (THIN_NCOL - 1)     # thin block start index
                if ts < 0:
                    break
                try:
                    out.append(dict(
                        x=float(p[ts + 1]), y=float(p[ts + 2]),
                        M_SED=float(p[ts + 7]), qual=tok,
                        afwhm=float(p[AFWHM_250_IDX])))
                except (ValueError, IndexError):
                    pass
                break
    return out


def parse_truth(path):
    """Return dict id -> {x, y, SD_emb, M_BE}."""
    out = {}
    for line in open(path):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        p = line.split()
        try:
            out[int(p[0])] = dict(x=float(p[2]), y=float(p[3]),
                                  SD_emb=float(p[5]), M_BE=float(p[8]))
        except (ValueError, IndexError):
            pass
    return out


def match_and_correct(truth, cat, corrector, match_rad):
    """Nearest-position match (<= match_rad px), ok fits only. Returns Nx3 array
    of (M_BE, M_SED, M_corr)."""
    rows = []
    for t in truth.values():
        best, bd = None, 1e30
        for s in cat:
            d = np.hypot(s['x'] - t['x'], s['y'] - t['y'])
            if d < bd:
                bd, best = d, s
        if best is None or bd > match_rad:
            continue
        if best['qual'] != 'ok' or best['M_SED'] <= 0 or best['afwhm'] <= 0:
            continue
        mc = corrector.correct(best['M_SED'], t['SD_emb'], best['afwhm'])
        if np.isfinite(mc):
            rows.append((t['M_BE'], best['M_SED'], mc))
    return np.array(rows)


def bin_stats(x, y, edges):
    mid, med, p16, p84 = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sub = y[(x >= lo) & (x < hi)]
        if len(sub) >= 2:
            mid.append(np.sqrt(lo * hi))
            med.append(np.median(sub))
            p16.append(np.percentile(sub, 16))
            p84.append(np.percentile(sub, 84))
    return map(np.array, (mid, med, p16, p84))


def make_plots(rows, out_pdf, label, scale):
    M_BE, M_SED, M_corr = rows[:, 0], rows[:, 1], rows[:, 2]
    edges = np.logspace(np.log10(0.008), np.log10(50), 10)
    lim = [0.005, 50]
    panels = [
        (M_SED / M_BE,  r'$M_\mathrm{SED}/M_\mathrm{BE}$',
         f'Mass bias -- {label} ({scale})',            'steelblue',  'firebrick'),
        (M_corr / M_BE, r'$M_\mathrm{corr}/M_\mathrm{BE}$',
         f'Corrected accuracy -- {label} ({scale})',   'darkorange', 'darkgreen'),
        (M_corr / M_SED, r'$M_\mathrm{corr}/M_\mathrm{SED}$',
         f'Applied correction factor -- {label} ({scale})', 'purple', 'firebrick'),
    ]
    with PdfPages(out_pdf) as pdf:
        for yv, yl, title, col, mcol in panels:
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.scatter(M_BE, yv, s=22, c=col, alpha=0.75, zorder=3)
            ax.axhline(1.0, ls='--', lw=1.2, c='gray')
            mid, med, p16, p84 = bin_stats(M_BE, yv, edges)
            if len(mid):
                ax.plot(mid, med, 'o-', c=mcol, ms=6, lw=1.5, label='bin median')
                ax.fill_between(mid, p16, p84, color=mcol, alpha=0.15, label='16-84%')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel(r'$M_\mathrm{BE}\ (M_\odot)$', fontsize=12)
            ax.set_ylabel(yl, fontsize=12)
            ax.set_title(title, fontsize=11)
            ax.set_xlim(lim)
            ax.legend(fontsize=9)
            ax.grid(True, which='both', ls=':', alpha=0.4)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('truth', help='injection truth table')
    ap.add_argument('catalog', help='getsf extraction catalog')
    ap.add_argument('out_pdf', help='output PDF path')
    ap.add_argument('--label', default='cSD', help='dataset label for titles')
    ap.add_argument('--scale', default='x1', help='scale label for titles')
    ap.add_argument('--grid', default='bes_model_params_catalog',
                    help='RT model grid catalog')
    ap.add_argument('--fwhm-key', default='FWHM250bs',
                    help='grid FWHM column for the corrector')
    ap.add_argument('--match-rad', type=float, default=4.5,
                    help='max match distance in pixels (default 4.5 = 13.5")')
    args = ap.parse_args()

    grid = load_catalog(args.grid)
    corrector = ObservableCorrector(grid, fwhm_key=args.fwhm_key, mode='direct')

    truth = parse_truth(args.truth)
    cat = parse_catalog(args.catalog)
    rows = match_and_correct(truth, cat, corrector, args.match_rad)

    if not len(rows):
        print('No finite matches -- nothing to plot.')
        return

    M_BE, M_SED, M_corr = rows[:, 0], rows[:, 1], rows[:, 2]
    print(f'{args.label}: truth={len(truth)}  matched(ok,finite)={len(rows)}')
    print(f'  M_SED/M_BE   median={np.median(M_SED / M_BE):.3f}  '
          f'rms={np.std(M_SED / M_BE):.3f}')
    print(f'  M_corr/M_BE  median={np.median(M_corr / M_BE):.3f}  '
          f'rms={np.std(M_corr / M_BE):.3f}')
    print(f'  M_corr/M_SED median={np.median(M_corr / M_SED):.3f}')

    make_plots(rows, args.out_pdf, args.label, args.scale)
    print(f'  -> {args.out_pdf}')


if __name__ == '__main__':
    main()
