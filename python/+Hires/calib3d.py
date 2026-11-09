#!/usr/bin/env python3
"""calib3d.py -- the correction factor as a function of three quantities the reconstruction provides.

Two observables are not enough.  The surface density and the temperature that a single-temperature fit returns do not determine
the configuration of a line of sight: within one cell of the plane they span, the correction that is actually required varies by
about a fifth, because the crest of a structure superposed on a warm diffuse medium and the flank of a structure without one give
the same pair.  This script adds the observable that separates them.

The third observable is the contrast between two resolutions of the reconstruction, the ratio of the surface density at the finest
beam to that at the coarsest.  At the crest of a structure the fine image stands above the coarse one, on a flank the two are
equal, and just outside a structure the fine one falls below.  That is morphological information, and the reconstruction supplies
it without any decomposition of the image.

Because a contrast depends on the size of a structure relative to the beams, the calibration acquires an angular size as a new
parameter and its members become images rather than single lines of sight.  For each member the script

  - builds the structure, a uniform sphere or a uniform cylinder of a given central column and angular radius, embedded in a
    plane-parallel diffuse medium of a given column with a given fraction of it in front;
  - computes, for every impact parameter, the temperatures along that line of sight exactly as calib2d.py does, and the
    intensity in each waveband that the line of sight emits;
  - convolves the intensities and the true surface density to the coarsest and to the finest beam;
  - fits one modified blackbody per pixel at each beam, in the wavebands and with the weights the reconstruction uses, which
    gives the surface density and the temperature that an observer would measure at each resolution;
  - records, for every pixel, the three observables and the factor that restores the true surface density at the fine beam.

The members are then binned in the three observables and the correction adopted for a cell is the median of the members in it,
exactly as in two dimensions.  The reduction in the spread within a cell, relative to calib2d.py, measures how much the contrast
adds.

Writes calib3d.npz.

Usage:  python3 calib3d.py [<output file>]
"""

import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hicorrect2 import planck, cnu, temperature, KAPHEAT, TFLOOR, TSURF, CAL_ERR

BANDS = [160.0, 250.0, 350.0, 500.0]
RELERR = np.array([CAL_ERR.get(w, 0.20) for w in BANDS])
PIX = 3.0
BEAM, FINE = 36.3, 13.5
F2S = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

NDIR = 400
NRAY = 60
NCIR = 40
NTFIT = 1200
NGRID = 160              # half width of the model image, in pixels

SIG0 = np.geomspace(2.0e21, 4.0e23, 12)
RADIUS = np.array([12.0, 22.0, 40.0, 75.0, 140.0])       # angular radius of the structure, arcsec
SIGC = np.array([0.0, 1.5e21, 5.0e21, 1.2e22])
QFRONT = np.array([0.0, 0.5])

TGRID = np.geomspace(0.98 * TFLOOR, 1.02 * TSURF, NTFIT)
MODEL = np.array([[cnu(w) * planck(t, w) for w in BANDS] for t in TGRID])
MODEL_N = MODEL / MODEL.max()

_k = np.arange(NDIR) + 0.5
_nz = 1.0 - 2.0 * _k / NDIR
_nr = np.sqrt(np.maximum(1.0 - _nz ** 2, 0.0))
_ph = np.pi * (1.0 + 5.0 ** 0.5) * _k
NX, NZ = _nr * np.cos(_ph), _nz


def fit_many(obs):
    """Surface density and temperature of a single modified blackbody fitted to each row of obs."""
    o = np.maximum(obs, 1.0e-30)
    u = 1.0 / (o * RELERR ** 2)
    num = u @ MODEL_N.T
    den = (u / o) @ (MODEL_N ** 2).T
    j = np.argmax(num * num / np.maximum(den, 1.0e-300), axis=1)
    t = TGRID[j]
    return o[:, -1] / (cnu(BANDS[-1]) * planck(t, BANDS[-1])), t


def cirrus_profile(sigma_c):
    """Temperatures and masses of the parcels of a plane-parallel medium of column sigma_c."""
    if sigma_c <= 0.0:
        return np.zeros(0), np.zeros(0)
    x = (np.arange(NCIR) + 0.5) / NCIR * sigma_c
    mu = (np.arange(NDIR) + 0.5) / NDIR * 2.0 - 1.0
    d = np.where(mu > 0.0, x[:, None] / np.abs(mu), (sigma_c - x[:, None]) / np.abs(mu))
    j = np.exp(-KAPHEAT * d).mean(axis=1)
    return temperature(-np.log(np.maximum(j, 1e-300)) / KAPHEAT), np.full(NCIR, sigma_c / NCIR)


def line_of_sight(sigma0, p, shape, sigma_c, q, tc, wc):
    """Intensities and true column of the line of sight at impact parameter p, in units of the structure radius."""
    if p < 1.0:
        half = np.sqrt(1.0 - p * p)
        l = (np.arange(NRAY) + 0.5) / NRAY * 2.0 * half - half
        r2 = (p * p + l * l)[:, None]
        xn = p * NX[None, :] + l[:, None] * NZ[None, :]
        if shape == "sphere":
            t = -xn + np.sqrt(np.maximum(1.0 - r2 + xn ** 2, 0.0))
        else:
            npar = (NX ** 2 + NZ ** 2)[None, :]
            t = (-xn + np.sqrt(np.maximum(npar * (1.0 - r2) + xn ** 2, 0.0))) / np.maximum(npar, 1.0e-8)
        col = t * sigma0 / 2.0
        if sigma_c > 0.0:
            mu = np.maximum(np.abs(NZ[None, :]), 1.0e-3)
            col = col + np.where(NZ[None, :] > 0.0, q * sigma_c / mu, (1.0 - q) * sigma_c / mu)
        js = np.exp(-KAPHEAT * col).mean(axis=1)
        ts = temperature(-np.log(np.maximum(js, 1e-300)) / KAPHEAT)
        ws = np.full(NRAY, sigma0 * half / NRAY)
        sig = sigma0 * half + sigma_c
    else:
        ts, ws, sig = np.zeros(0), np.zeros(0), sigma_c
    tt = np.concatenate([ts, tc])
    ww = np.concatenate([ws, wc])
    if len(tt) == 0:
        return np.zeros(len(BANDS)), 0.0
    return np.array([cnu(w) * np.sum(ww * planck(tt, w)) for w in BANDS]), sig


def member(sigma0, radius, shape, sigma_c, q):
    """One model image: the three observables and the required factor at every radius, at the fine beam."""
    tc, wc = cirrus_profile(sigma_c)
    npix = min(NGRID, int(3.0 * radius / PIX) + 40)
    r = (np.arange(npix) + 0.5) * PIX
    prof = np.array([line_of_sight(sigma0, rr / radius, shape, sigma_c, q, tc, wc) for rr in r], dtype=object)
    I = np.array([p[0] for p in prof])
    S = np.array([p[1] for p in prof], dtype=float)
    if shape == "cylinder":
        full_I = np.concatenate([I[::-1], I], axis=0)
        full_S = np.concatenate([S[::-1], S])
        out = []
        for beam in (BEAM, FINE):
            s = beam * F2S / PIX
            cI = np.stack([gaussian_filter1d(full_I[:, k], s, mode="nearest") for k in range(len(BANDS))], axis=-1)
            cS = gaussian_filter1d(full_S, s, mode="nearest")
            out.append((cI[npix:], cS[npix:]))
    else:
        y, x = np.mgrid[-npix:npix, -npix:npix]
        rr = np.hypot(y + 0.5, x + 0.5) * PIX
        idx = np.minimum((rr / PIX).astype(int), npix - 1)
        out = []
        for beam in (BEAM, FINE):
            s = beam * F2S / PIX
            cI = np.stack([gaussian_filter(I[idx, k], s, mode="nearest") for k in range(len(BANDS))], axis=-1)
            cS = gaussian_filter(S[idx], s, mode="nearest")
            out.append((cI[npix, npix:], cS[npix, npix:]))
    (I36, S36), (I13, S13) = out
    good = np.all(I36 > 0.0, axis=1) & np.all(I13 > 0.0, axis=1) & (S13 > 0.0)
    if not good.any():
        return np.zeros((0, 4))
    f36, t36 = fit_many(I36[good])
    f13, t13 = fit_many(I13[good])
    ok = (f36 > 0.0) & (f13 > 0.0)
    return np.column_stack([np.log10(f13[ok]), t13[ok], np.log10(f13[ok] / f36[ok]), S13[good][ok] / f13[ok]])


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "calib3d.npz"
    rows = []
    for shape in ("sphere", "cylinder"):
        for s0 in SIG0:
            for rad in RADIUS:
                for sc in SIGC:
                    for q in (QFRONT if sc > 0.0 else [0.0]):
                        rows.append(member(s0, rad, shape, sc, q))
        print(" %s done" % shape)
    rows = np.vstack([r for r in rows if len(r)])
    ex = np.arange(19.4, 23.8, 0.06)
    ey = np.arange(7.0, 23.0, 0.25)
    ez = np.array([-0.30, -0.06, -0.02, 0.02, 0.06, 0.30])
    pts, val, spread = [], [], []
    for i in range(len(ex) - 1):
        m1 = (rows[:, 0] >= ex[i]) & (rows[:, 0] < ex[i + 1])
        if not m1.any():
            continue
        s1 = rows[m1]
        for j in range(len(ey) - 1):
            m2 = (s1[:, 1] >= ey[j]) & (s1[:, 1] < ey[j + 1])
            if not m2.any():
                continue
            s2 = s1[m2]
            for k in range(len(ez) - 1):
                m3 = (s2[:, 2] >= ez[k]) & (s2[:, 2] < ez[k + 1])
                if m3.sum() >= 3:
                    f = s2[m3, 3]
                    pts.append((0.5 * (ex[i] + ex[i + 1]), 0.5 * (ey[j] + ey[j + 1]), 0.5 * (ez[k] + ez[k + 1])))
                    val.append(np.median(f))
                    spread.append(0.5 * (np.percentile(f, 84) - np.percentile(f, 16)) / np.median(f))
    pts, val, spread = np.array(pts), np.array(val), np.array(spread)
    np.savez(out, members=rows, points=pts, factor=val, spread=spread)
    print(" %d members, %d cells" % (len(rows), len(val)))
    print(" spread within a cell: median %.1f%%, 90th percentile %.1f%%"
          % (100.0 * np.median(spread), 100.0 * np.percentile(spread, 90)))
    print(" written %s" % out)


if __name__ == "__main__":
    main()
