#!/usr/bin/env python3
"""hicorrect3.py -- correct the temperature bias of a surface density image, with a factor measured on radiative transfer models.

The correction of hicorrect2 derives its factor analytically, from a uniform structure whose shielding is an angular average.
This one reads the factor from a table measured directly on the grid of models, by harvest_calib3.py, and it can use either one
observable or three.

  One observable. The factor is a function of the surface density alone, exactly as in hicorrect2, but measured rather than
  derived. Validated on a second grid of models whose parameters lie between those of the first, it leaves 9.2% of the deficit at
  the center of a core against 16.3% for the analytic factor, and recovers 0.990 of the true mass integrated over the model
  against 0.947. This form is meaningful only for a table harvested at one radiation field: a table combining grids computed for
  several fields averages over them, and the one-observable curve it gives is not smooth.

  Three observables. The factor is a function of the surface density, of the temperature fitted at the same resolution, and of the
  contrast between that resolution and the coarsest, which the reconstruction also provides. Those carry the information that one
  number cannot: whether the material of a pixel lies near the surface of a structure, where it is warm, or in the body of one,
  where it is cold, and, through the temperature, the radiation field that heats it. It removes most of the ring of
  over-correction that the one-observable forms leave on the flanks of structures.

The interpolation is bilinear in the surface density and the temperature on the grid of the table and linear between the levels
of the contrast. Cells absent from the table carry no weight, and their weight is given to a two-observable fallback: the same
table averaged over the contrast, as a function of the surface density and the temperature, completed by interpolation along both
coordinates and by the nearest cell beyond them, so that it is defined everywhere. The fraction of pixels not fully covered by the
three-observable table is reported. Tables written by harvest_calib3.py --fill have their gaps filled already, and cells with a
pixel count of 0 are those filled cells.

An image that already carries a correction is refused, since applying one twice is meaningless.

usage: hicorrect3.py <table> <image> [<coarse image> <temperature>] [<image> <temperature> ...]

  <table>          the calibration written by harvest_calib3.py
  <image>          a surface density image to correct
  <coarse image>   the surface density at the coarsest resolution, which sets the contrast
  <temperature>    the temperature image of the same resolution as the image being corrected

With one argument after the table the correction uses the surface density alone. With three or more it uses all three
quantities, and every further pair adds an image at another resolution.

output: <image>.corr3.fits, the corrected image, and <image>.cfactor3.fits, the factor applied to it
"""

import os
import sys

import numpy as np
from astropy.io import fits

KEYWORD = "HICORR3"


class Calibration(object):
    """The table of harvest_calib3.py, as a regular grid in the surface density and the temperature at each level of contrast."""

    def __init__(self, path):
        d = np.loadtxt(path, comments="!")
        if d.ndim != 2 or d.shape[1] < 7:
            raise ValueError("%s does not look like the output of harvest_calib3.py" % path)
        self.raw = d
        self.x = np.unique(d[:, 0])
        self.y = np.unique(d[:, 1])
        self.z = np.unique(d[:, 2])
        self.f = np.full((len(self.x), len(self.y), len(self.z)), np.nan)
        self.s = np.full(self.f.shape, np.nan)
        self.n = np.zeros(self.f.shape)
        ix = np.searchsorted(self.x, d[:, 0])
        iy = np.searchsorted(self.y, d[:, 1])
        iz = np.searchsorted(self.z, d[:, 2])
        self.f[ix, iy, iz] = d[:, 5]
        self.s[ix, iy, iz] = 0.5 * (d[:, 6] - d[:, 4]) / d[:, 5]
        self.n[ix, iy, iz] = d[:, 3]
        # the one-observable curve: the same table averaged over the other two, weighted by the pixels of the models
        w = np.where(np.isfinite(self.f), self.n, 0.0)
        v = np.where(np.isfinite(self.f), self.f, 0.0)
        tot = w.sum(axis=(1, 2))
        self.f1 = np.where(tot > 0.0, (v * w).sum(axis=(1, 2)) / np.maximum(tot, 1.0), np.nan)
        keep = np.isfinite(self.f1)
        self.x1, self.f1 = self.x[keep], self.f1[keep]
        self._fallback()

    def _fallback(self):
        """The two-observable table, the factor against the surface density and the temperature, which replaces a missing cell.

        Each cell is the average of the factor over the levels of contrast that exist at that surface density and temperature,
        weighted by the pixels of the models; where only filled cells (with no pixels) exist, it is their plain average. The table
        is then completed so that it is defined everywhere: gaps along the temperature are interpolated at fixed surface density,
        then gaps along the surface density at fixed temperature, and any cell still empty takes the value of the nearest cell,
        with distances counted in steps of the table and ties resolved in favor of the lower surface density, then the lower
        temperature. The temperature carries the radiation field, so the fallback distinguishes fields where the one-observable
        curve would average over them.
        """
        ok = np.isfinite(self.f)
        n = np.where(ok, self.n, 0.0)
        wn = n.sum(axis=2)
        cnt = ok.sum(axis=2)
        fv = np.where(ok, self.f, 0.0)
        sv = np.where(ok & np.isfinite(self.s), self.s, 0.0)
        f2 = np.where(wn > 0.0, (fv * n).sum(axis=2) / np.maximum(wn, 1.0),
                      np.where(cnt > 0, fv.sum(axis=2) / np.maximum(cnt, 1), np.nan))
        s2 = np.where(wn > 0.0, (sv * n).sum(axis=2) / np.maximum(wn, 1.0),
                      np.where(cnt > 0, sv.sum(axis=2) / np.maximum(cnt, 1), np.nan))
        for axis, c in ((1, self.y), (0, self.x)):
            for a in (f2, s2):
                g = np.moveaxis(a, axis, -1)
                for idx in np.ndindex(g.shape[:-1]):
                    js = np.nonzero(np.isfinite(g[idx]))[0]
                    for j1, j2 in zip(js[:-1], js[1:]):
                        if j2 - j1 > 1:
                            t = (c[j1 + 1:j2] - c[j1]) / (c[j2] - c[j1])
                            g[idx][j1 + 1:j2] = (1.0 - t) * g[idx][j1] + t * g[idx][j2]
        have = np.isfinite(f2)
        hi, hj = np.nonzero(have)
        for i, j in zip(*np.nonzero(~have)):
            k = int(np.argmin((hi - i) ** 2 + (hj - j) ** 2))
            f2[i, j], s2[i, j] = f2[hi[k], hj[k]], s2[hi[k], hj[k]]
        self.f2, self.s2 = f2, s2

    def two(self, xi, yi):
        """The factor and its spread from the surface density and the temperature, bilinear in the completed fallback table."""
        i0 = np.clip(np.searchsorted(self.x, xi) - 1, 0, len(self.x) - 2)
        j0 = np.clip(np.searchsorted(self.y, yi) - 1, 0, len(self.y) - 2)
        tx = np.clip((xi - self.x[i0]) / (self.x[i0 + 1] - self.x[i0]), 0.0, 1.0)
        ty = np.clip((yi - self.y[j0]) / (self.y[j0 + 1] - self.y[j0]), 0.0, 1.0)
        f = np.zeros(xi.shape)
        s = np.zeros(xi.shape)
        for di, wx in ((0, 1.0 - tx), (1, tx)):
            for dj, wy in ((0, 1.0 - ty), (1, ty)):
                f += wx * wy * self.f2[i0 + di, j0 + dj]
                s += wx * wy * self.s2[i0 + di, j0 + dj]
        return f, s

    def one(self, sigma):
        """The factor from the surface density alone."""
        return np.interp(np.log10(np.maximum(sigma, 1.0)), self.x1, self.f1)

    def _plane(self, k, xi, yi):
        """Bilinear interpolation of one level of contrast, with the weight of the cells that exist.

        The coordinates of the table are not a complete regular sequence: a cell holding too few pixels of the models is absent
        from the file altogether, so the distinct values have gaps in them. The bracketing pair is therefore found by searching
        the values that are present rather than by arithmetic on a step, and the weights use the spacing of that pair, which is
        exact where the sequence is regular and interpolates across a gap where it is not.
        """
        i0 = np.clip(np.searchsorted(self.x, xi) - 1, 0, len(self.x) - 2)
        j0 = np.clip(np.searchsorted(self.y, yi) - 1, 0, len(self.y) - 2)
        tx = np.clip((xi - self.x[i0]) / (self.x[i0 + 1] - self.x[i0]), 0.0, 1.0)
        ty = np.clip((yi - self.y[j0]) / (self.y[j0 + 1] - self.y[j0]), 0.0, 1.0)
        num = np.zeros(xi.shape)
        den = np.zeros(xi.shape)
        spr = np.zeros(xi.shape)
        for di, wx in ((0, 1.0 - tx), (1, tx)):
            for dj, wy in ((0, 1.0 - ty), (1, ty)):
                v = self.f[i0 + di, j0 + dj, k]
                u = self.s[i0 + di, j0 + dj, k]
                w = np.where(np.isfinite(v), wx * wy, 0.0)
                num += w * np.where(np.isfinite(v), v, 0.0)
                spr += w * np.where(np.isfinite(u), u, 0.0)
                den += w
        return num, spr, den

    def three(self, sigma, temper, coarse):
        """The factor from all three quantities, falling back smoothly to the two-observable table where cells are missing."""
        xi = np.log10(np.maximum(sigma, 1.0))
        yi = temper
        zi = np.log10(np.maximum(sigma, 1.0) / np.maximum(coarse, 1.0))
        zc = np.clip(zi, self.z[0], self.z[-1])
        k0 = np.clip(np.searchsorted(self.z, zc) - 1, 0, len(self.z) - 2)
        tz = (zc - self.z[k0]) / (self.z[k0 + 1] - self.z[k0])
        num = np.zeros(xi.shape)
        spr = np.zeros(xi.shape)
        den = np.zeros(xi.shape)
        for dk in (0, 1):
            wz = (1.0 - tz) if dk == 0 else tz
            for k in np.unique(k0 + dk):
                m = (k0 + dk) == k
                if not m.any():
                    continue
                a, b, c = self._plane(k, xi[m], yi[m])
                num[m] += wz[m] * a
                spr[m] += wz[m] * b
                den[m] += wz[m] * c
        base, bspr = self.two(xi, yi)
        w = np.clip(den, 0.0, 1.0)
        return np.where(w > 0.0, w * num / np.maximum(den, 1.0e-12) + (1.0 - w) * base, base), \
            np.where(w > 0.0, w * spr / np.maximum(den, 1.0e-12) + (1.0 - w) * bspr, bspr), w


def write(name, data, header, table, unit=None):
    h = fits.PrimaryHDU(np.float32(data), header=header.copy())
    h.header[KEYWORD] = (os.path.basename(table), "calibration of hicorrect3 applied")
    if unit is not None:
        h.header["BUNIT"] = (unit, "dimensionless")
    h.writeto(name, overwrite=True)
    print("   written %s" % name)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return
    table = sys.argv[1]
    cal = Calibration(table)
    print(" calibration %s: %d cells, %d levels of contrast, %.3e pixels of the models"
          % (table, len(cal.raw), len(cal.z), cal.raw[:, 3].sum()))
    print(" covering log10 Sigma %.2f to %.2f, T %.1f to %.1f K, contrast %.3f to %.3f"
          % (cal.x.min(), cal.x.max(), cal.y.min(), cal.y.max(), cal.z.min(), cal.z.max()))

    rest = sys.argv[2:]
    if len(rest) == 1:
        jobs, coarse = [(rest[0], None)], None
        print(" using the surface density alone")
    elif len(rest) >= 3 and (len(rest) - 1) % 2 == 0:
        coarse = rest[1]
        jobs = [(rest[0], rest[2])] + [(rest[i], rest[i + 1]) for i in range(3, len(rest), 2)]
        print(" using the surface density, the temperature and the contrast against %s" % coarse)
    else:
        print(__doc__)
        return

    sc = np.array(fits.getdata(coarse), float).squeeze() if coarse else None
    for image, temper in jobs:
        hi = fits.open(image)[0]
        s = np.array(hi.data, float).squeeze()
        if KEYWORD in hi.header or "HICORR" in hi.header:
            print(" %s already carries a correction; skipped" % image)
            continue
        good = np.isfinite(s) & (s > 0.0)
        if temper is None:
            f = np.where(good, cal.one(s), 1.0)
            sp = np.full(s.shape, np.nan)
            cov = None
        else:
            t = np.array(fits.getdata(temper), float).squeeze()
            if s.shape != t.shape or s.shape != sc.shape:
                print(" %s, %s and %s differ in shape; skipped" % (image, temper, coarse))
                continue
            good &= np.isfinite(t) & (t > 0.0) & np.isfinite(sc) & (sc > 0.0)
            f, sp, w = cal.three(np.where(good, s, 1.0e21), np.where(good, t, 14.0), np.where(good, sc, 1.0e21))
            f = np.where(good, f, 1.0)
            cov = w[good]
        print(" %s: %d valid pixels" % (image, int(good.sum())))
        if cov is not None:
            print("   covered by the table: %.1f%% of the pixels fully, %.1f%% not at all"
                  % (100.0 * np.mean(cov > 0.999), 100.0 * np.mean(cov < 0.001)))
            print("   spread within the cells used, median %.1f%%" % (100.0 * np.nanmedian(sp[good])))
        print("   factor  1%%: %.3f  25%%: %.3f  50%%: %.3f  75%%: %.3f  99%%: %.3f"
              % tuple(np.percentile(f[good], [1, 25, 50, 75, 99])))
        print("   ratio of the sums, corrected to uncorrected: %.4f" % ((f * s)[good].sum() / s[good].sum()))
        stem = image[:-5] if image.endswith(".fits") else image
        write(stem + ".corr3.fits", s * f, hi.header, table)
        write(stem + ".cfactor3.fits", np.where(good, f, np.nan), hi.header, table, unit="")


if __name__ == "__main__":
    main()
