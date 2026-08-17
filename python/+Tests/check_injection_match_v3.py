#!/usr/bin/env python
"""check_injection_match_v3.py -- compare the predicted properties of the
injected cores with those of the real sources.

Input is a truth table written by run_inject_v3.py, or the table written by
dryrun_placement_v3.py, together with the target table written by
build_joint_target_v3.py.  The two input formats have different columns and are
distinguished automatically by their header.  A truth table does not record the
contrast nor the reported SED mass of a core, so those are recovered from the
grid catalogue through the node tag; give its path with -g if it is not in the
working directory under its usual name.

Three quantities are compared, each on a base-10 logarithmic scale:

    Sigma_cloud  background column density under the source.  For a real
                 source this is the getsf quantity PEAK^BGF; for an injected
                 core it is the column of the map averaged over a disk of
                 radius equal to the observed size of the model.
    contrast     PEAK^SBF / PEAK^BGF = (source + background) / background,
                 greater than unity.  For an injected core it is
                 1 + ICSDbs / local_Sigma.
    M_reported   reported SED mass.  For a real source this is the mass in the
                 getsf catalogue; for an injected core it is the four-band SED
                 mass M_SED4bs of the grid node, which is the like-for-like
                 counterpart, not the true mass M_BE.

Reported for each quantity are the median in both samples and the
two-sample Kolmogorov-Smirnov statistic with its probability.  Reported for
the pair (Sigma_cloud, contrast) are the Spearman rank correlation in each
sample and the Bhattacharyya overlap of the two-dimensional histograms, which
equals unity for identical distributions.

Usage:
    python check_injection_match_v3.py dryrun_placement_v3.txt \\
        [-t joint_target_Aquila_v3.txt] [-o fig_injection_match_v3.pdf]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                     # noqa: E402
from scipy.stats import ks_2samp, spearmanr, gaussian_kde   # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from aa_plotstyle import clean_log
except Exception:                                    # pragma: no cover
    def clean_log(ax):
        return


def _load_grid(path, suffix='bsl'):
    """Per-node ICSDbs, M_SED4bs and FWHMSDbs, keyed by node tag."""
    hdr = None
    for ln in open(path):
        if ln.lstrip()[:1] == '#' and 'ICSDbs' in ln.split():
            hdr = ln.lstrip('# ').split()
    if hdr is None:
        raise ValueError('no column header found in %s' % path)
    out = {}
    for ln in open(path):
        if not ln.strip() or ln.lstrip()[:1] == '#':
            continue
        f = ln.split()
        tag = 'i%02dj%02dk%02d' % (int(float(f[hdr.index('i')])),
                                   int(float(f[hdr.index('j')])),
                                   int(float(f[hdr.index('k')])))
        out[tag] = (float(f[hdr.index('ICSD' + suffix)]),
                    float(f[hdr.index('M_SED4' + suffix)]),
                    float(f[hdr.index('FWHMSD' + suffix)]))
    return out


def _is_truth_table(path):
    """A truth table written by run_inject_v3.py begins with '# Injection
    truth table'; a table written by dryrun_placement_v3.py does not."""
    for ln in open(path):
        if ln.lstrip()[:1] in '#!':
            if 'truth table' in ln.lower():
                return True
            continue
        return False
    return False


def read_injected(path, grid_path, suffix='bsl'):
    """Background column, contrast, reported SED mass and observed size of the
    injected cores, from either kind of table.

    dryrun_placement_v3.py writes all four directly.  The truth table written
    by run_inject_v3.py records the position, the local column and the true
    properties of the node, but not the contrast nor the reported SED mass, so
    those two are recovered from the grid catalogue via the node tag:

        contrast = 1 + ICSDbs / local_Sigma
        M_SED    = M_SED4bs of the node

    Truth-table columns are id, model, x_pix, y_pix, local_Sigma, SD_emb,
    T_BE, rho_BE, M_BE, R_BE_as, FWHM_SD, scale.
    """
    S, C, M, F = [], [], [], []
    if not _is_truth_table(path):
        for ln in open(path):
            if not ln.strip() or ln.lstrip()[0] in '#!':
                continue
            f = ln.split()
            S.append(float(f[4])); C.append(float(f[5]))
            M.append(float(f[7])); F.append(float(f[8]))
        return (np.array(S), np.array(C), np.array(M), np.array(F))

    g = _load_grid(grid_path, suffix)
    missing = set()
    for ln in open(path):
        if not ln.strip() or ln.lstrip()[0] in '#!':
            continue
        f = ln.split()
        tag = f[1]
        if tag not in g:
            missing.add(tag)
            continue
        sig = float(f[4])
        icsd, msed, fw_grid = g[tag]
        # FWHM_SD is written by run_inject_v3.py but not by version 2; fall
        # back to the grid value when the column is absent
        fw = float(f[10]) if len(f) > 11 else fw_grid
        S.append(sig); C.append(1.0 + icsd / sig)
        M.append(msed); F.append(fw)
    if missing:
        print('WARNING: %d node tags in the truth table are absent from the '
              'grid catalogue and were skipped: %s'
              % (len(missing), ', '.join(sorted(missing)[:5])))
    return (np.array(S), np.array(C), np.array(M), np.array(F))


def read_target(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]


def overlap2d(x1, y1, x2, y2, nb=24):
    lo = [min(x1.min(), x2.min()), min(y1.min(), y2.min())]
    hi = [max(x1.max(), x2.max()), max(y1.max(), y2.max())]
    h1, _, _ = np.histogram2d(x1, y1, bins=nb,
                              range=[[lo[0], hi[0]], [lo[1], hi[1]]])
    h2, _, _ = np.histogram2d(x2, y2, bins=nb,
                              range=[[lo[0], hi[0]], [lo[1], hi[1]]])
    return float(np.sum(np.sqrt((h1 / h1.sum()) * (h2 / h2.sum()))))


def contours(ax, x, y, color):
    k = gaussian_kde(np.vstack([x, y]))
    gx = np.linspace(x.min(), x.max(), 60)
    gy = np.linspace(y.min(), y.max(), 60)
    X, Y = np.meshgrid(gx, gy)
    Z = k(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    Zs = np.sort(Z.ravel())[::-1]
    cum = np.cumsum(Zs) / Zs.sum()
    lev = [Zs[np.searchsorted(cum, f)] for f in (0.39, 0.68, 0.86)][::-1]
    ax.contour(10 ** X, 10 ** Y, Z, levels=lev, colors=color, linewidths=1.2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('predicted',
                    help='truth table from run_inject_v3.py, or the table '
                         'written by dryrun_placement_v3.py')
    ap.add_argument('-t', '--target', default='joint_target_Aquila_v3.txt')
    ap.add_argument('-o', default='fig_injection_match_v3.pdf')
    ap.add_argument('--stamp-suffix', default='bsl', choices=('bs', 'bsl'),
                    help='background convention of the injected stamps; must '
                         'match STAMP_SUFFIX of the run being analysed')
    ap.add_argument('-g', '--grid',
                    default='/Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260721_RT_BES_radmc3d/bes_model_grid_final2_catalog',
                    help='grid catalogue, needed only for a truth table')
    a = ap.parse_args()

    Si, Ci, Mi, Fi = read_injected(a.predicted, a.grid, a.stamp_suffix)
    Mr, Sr, Cr, Fr = read_target(a.target)

    print('injected cores: %d    real sources: %d\n' % (len(Si), len(Sr)))
    print('%-42s %12s %12s %14s' % ('quantity', 'injected', 'real',
                                    'KS / probability'))
    for nm, vi, vr in (('median Sigma_cloud (cm^-2)', Si, Sr),
                       ('median contrast PEAK^SBF/PEAK^BGF (-)', Ci, Cr),
                       ('median reported SED mass (M_sun)', Mi, Mr),
                       ('median observed size FWHM (arcsec)', Fi, Fr)):
        ks = ks_2samp(np.log10(vi), np.log10(vr))
        print('%-42s %12.4g %12.4g %8.3f / %.1e'
              % (nm, np.median(vi), np.median(vr), ks.statistic, ks.pvalue))

    xi, yi = np.log10(Si), np.log10(Ci)
    xr, yr = np.log10(Sr), np.log10(Cr)
    print('%-42s %12.3f %12.3f'
          % ('Spearman rank correlation of contrast',
             spearmanr(xi, yi).statistic, spearmanr(xr, yr).statistic))
    print('%-42s %12.3f'
          % ('on Sigma_cloud; two-dimensional overlap', overlap2d(xi, yi, xr, yr)))

    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    ax.scatter(Sr, Cr, s=3, color='0.55', alpha=0.35, lw=0,
               label='real, %d' % len(Sr))
    ax.scatter(Si, Ci, s=6, color='crimson', alpha=0.55, lw=0,
               label='injected, %d' % len(Si))
    contours(ax, xr, yr, '0.25')
    contours(ax, xi, yi, 'darkred')
    ax.set_xscale('log'); ax.set_yscale('log'); clean_log(ax)
    ax.set_xlabel(r'$\Sigma_{\rm cloud}$  (cm$^{-2}$)', fontsize=8)
    ax.set_ylabel(r'contrast  PEAK$^{\rm SBF}$/PEAK$^{\rm BGF}$', fontsize=8)
    ax.legend(fontsize=7, loc='upper right')
    ax.tick_params(labelsize=7)
    fig.tight_layout(); fig.savefig(a.o)
    print('\nwrote %s' % a.o)


if __name__ == '__main__':
    main()
