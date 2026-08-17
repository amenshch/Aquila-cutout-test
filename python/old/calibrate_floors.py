#!/usr/bin/env python3
"""
Calibrate a per-band detection floor from a getsf source catalogue.

The floor is the faintest source peak that the extraction still recovers at a
given local background level,

    floor_b(I_bg) = A_b * (I_bg / I_ref_b) ** p_b ,

expressed in the native units of band b.  It is the per-band generalisation of
the single surface-density floor used by add_recoverable_mass.py
(N_detect_floor = 1.971e21 * (N_back/1e22)**1.342).

Why this exists
---------------
The RT models are noiseless, so a recovered-flux fraction cannot be obtained
from them directly.  The floor supplies the missing information: it is measured
from where getsf actually stops detecting sources in a real map, and is then
applied to the noiseless model profile to decide how far out the source would
have been measured.  Calibrating it per band, and per field, removes the
dependence on any single reference cloud.

Two independent estimators are produced:

  envelope : a low percentile of PEAK^SRC_b in bins of PEAK^BGF_b.  This is the
             direct analogue of the existing surface-density floor and is the
             recommended one.
  sigma    : a multiple of getsf's own peak-flux error FXP_ERR_b, fitted the
             same way.  Used as a cross-check; where the two disagree strongly
             the envelope is preferred, since FXP_ERR is unreliable in some
             bands (in Aquila the 160 um and 255 um fits scatter by 30-170x).

Usage
-----
    python calibrate_floors.py <getsf_catalogue> [--pct 5] [--out floors.txt]
"""

import argparse
import re
import sys

import numpy as np

MIN_SOURCES = 50          # per band, below this the fit is not attempted
MIN_PER_BIN = 15
N_BINS = 6


def read_catalogue(path):
    """Return (columns dict-of-arrays, quality array, band list)."""
    lines = open(path).read().split('\n')
    hdr = None
    for l in lines:
        if l.startswith('#') and 'PEAK^SRC03' in l and 'TOTL^MASS' in l:
            hdr = l
            break
    if hdr is None:
        raise SystemExit('no getsf column header found in %s' % path)
    names = hdr.lstrip('#').split()
    idx = {n: i for i, n in enumerate(names)}
    rows = [l.split() for l in lines
            if l.strip() and l.lstrip()[0].isdigit() and len(l.split()) >= len(names)]

    def col(nm):
        i = idx[nm]
        out = np.empty(len(rows))
        for t, r in enumerate(rows):
            try:
                out[t] = float(r[i])
            except ValueError:
                out[t] = np.nan
        return out

    quality = np.array([r[idx['QUALITY']] for r in rows]) if 'QUALITY' in idx \
        else np.array(['ok'] * len(rows))
    # bands that carry the photometry we need
    bands = sorted({m.group(1) for n in names
                    for m in [re.match(r'PEAK\^BGF(\d\d)$', n)] if m})
    bands = [b for b in bands if ('PEAK^SRC' + b) in idx and ('FXP_ERR' + b) in idx]
    return col, quality, bands, idx


def band_units(path):
    """Map band index -> wavelength/units string from the obs.fits header comments."""
    lines = open(path).read().split('\n')[:40]
    wl = []
    for l in lines:
        m = re.search(r"\.(\d+)\.obs\.fits'\s*\(units:\s*([^)]+)\)", l)
        if m:
            wl.append((m.group(1), m.group(2).strip()))
    return wl


def fit_powerlaw(x, y):
    """Fit log y = log A + p log(x/xref); return (A, p, xref, scatter)."""
    lx, ly = np.log10(x), np.log10(y)
    xref = np.median(lx)
    p = np.polyfit(lx - xref, ly, 1)
    resid = ly - np.polyval(p, lx - xref)
    return 10 ** p[1], p[0], 10 ** xref, 10 ** np.std(resid)


def envelope_floor(bg, peak, pct, nbins=N_BINS):
    """Low-percentile peak in bins of background -> power-law floor."""
    q = np.percentile(bg, np.linspace(5, 95, nbins + 1))
    xs, ys = [], []
    for lo, hi in zip(q[:-1], q[1:]):
        m = (bg >= lo) & (bg < hi)
        if m.sum() >= MIN_PER_BIN:
            xs.append(np.median(bg[m]))
            ys.append(np.percentile(peak[m], pct))
    if len(xs) < 3:
        return None
    xs, ys = np.array(xs), np.array(ys)
    good = np.isfinite(xs) & np.isfinite(ys) & (xs > 0) & (ys > 0)
    if good.sum() < 3:
        return None
    return fit_powerlaw(xs[good], ys[good])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('catalogue')
    ap.add_argument('--pct', type=float, default=5.0,
                    help='percentile defining the detection envelope (default 5)')
    ap.add_argument('--sigma-k', type=float, default=2.0,
                    help='floor = k * FXP_ERR for the cross-check (default 2)')
    ap.add_argument('--out', default=None, help='write the fitted parameters here')
    args = ap.parse_args()

    col, quality, bands, idx = read_catalogue(args.catalogue)
    wl = band_units(args.catalogue)
    ok0 = (quality == 'ok')
    print('# %s' % args.catalogue.split('/')[-1])
    print('# bands in header comments: %s' % ', '.join('%s(%s)' % w for w in wl))
    print('#')
    print('# %-5s %5s  %12s %8s %9s   %12s %8s %9s'
          % ('band', 'n', 'A_env', 'p_env', 'scat', 'A_sig', 'p_sig', 'scat'))

    results = {}
    for b in bands:
        bg = col('PEAK^BGF' + b)
        pk = col('PEAK^SRC' + b)
        sg = col('FXP_ERR' + b)
        m = ok0 & np.isfinite(bg) & (bg > 0) & np.isfinite(pk) & (pk > 0)
        if m.sum() < MIN_SOURCES:
            continue
        env = envelope_floor(bg[m], pk[m], args.pct)
        ms = m & np.isfinite(sg) & (sg > 0)
        sig = fit_powerlaw(bg[ms], args.sigma_k * sg[ms]) if ms.sum() >= MIN_SOURCES else None
        if env is None:
            continue
        results[b] = dict(envelope=env, sigma=sig, n=int(m.sum()))
        e = env
        s = sig if sig else (np.nan,) * 4
        print('  %-5s %5d  %12.5g %8.3f %8.2fx   %12.5g %8.3f %8.2fx'
              % (b, m.sum(), e[0], e[1], e[3], s[0], s[1], s[3]))

    # flag bands where the two estimators disagree badly
    print('#')
    for b, r in results.items():
        if r['sigma'] is None:
            continue
        if r['sigma'][3] > 3.0 * r['envelope'][3]:
            print('# WARNING band %s: FXP_ERR fit scatters %.0fx vs envelope %.1fx '
                  '-- use the envelope' % (b, r['sigma'][3], r['envelope'][3]))

    if args.out:
        with open(args.out, 'w') as fh:
            fh.write('# band  A  p  I_ref   (floor = A*(I_bg/I_ref)**p, envelope estimator)\n')
            for b, r in results.items():
                A, p, xref, sc = r['envelope']
                fh.write('%s %.6e %.6f %.6e %.4f\n' % (b, A, p, xref, sc))
        print('# wrote %s' % args.out)
    return results


if __name__ == '__main__':
    main()
