#!/usr/bin/env python3
"""harvest_calib3.py -- build the calibration of hicorrect3 from the reconstructions of the model grid.

The correction of hicorrect2 is a function of one number, the surface density that a single-temperature fit returns, and it leaves
a residual that follows the structures: a ring of over-correction of 30 to 40 percent around every core and along the flanks of
every filament, which distorts shapes even where it leaves integrated quantities nearly right. The residual exists because one
number cannot distinguish the flank of a structure, where the material lies near the surface and is warm, from the body of a
uniform structure of the same column, where it is cold.

The reconstruction reports more than one number. At each resolution it writes a surface density and a temperature, and it writes
them at several resolutions. Three of those quantities carry the missing information:

    x = log10 ( Sigma_13.5 )                      the surface density of the finest image
    y = T_13.5                                    the temperature fitted at that resolution
    z = log10 ( Sigma_13.5 / Sigma_36.3 )         the contrast between the finest and the coarsest image

The first is what hicorrect2 already uses. The second says whether the material of the pixel is warmer or colder than a uniform
structure of that column would be, which is the signature of lying near a surface. The third says whether the pixel sits on a
crest, where the fine image stands above the coarse one, on a flank, where the two agree, or just outside a structure, where the
fine image falls below.

This script measures, over the whole grid of radiative transfer models, the factor that restores the true surface density as a
function of those three quantities. For every pixel of every node it forms x, y, z and

    f = Sigma_true,13.5 / Sigma_13.5 ,

accumulates them in a fixed four-dimensional histogram, over the three observables and the factor, and writes for each cell the
number of pixels, the median factor and the 16th and 84th percentiles. The percentiles are the spread that the three observables
cannot resolve, and they are the uncertainty of the correction at that point.

Nothing is fitted: the calibration is a direct measurement on radiative transfer models whose answer is known, and it is
transferable to any field reconstructed in the same wavebands with the same beams.

The ranges of the histogram must contain every pixel of every grid. The grids computed for a range of interstellar radiation
fields reach fitted temperatures above 40 K at G_0 = 482.5 and factors below 0.9 just outside structures, so the ranges are wide,
and the pixels falling outside each axis are counted and reported, so that a truncation can never pass unnoticed.

Each run writes, next to the table, the histogram itself as <output>.hist.npz. The calibration is a histogram, so the histograms
of several grids add exactly, and the combined table is made from their sum with --combine. Tables written by earlier runs hold
only percentiles and cannot be combined exactly.

The output is a few hundred kilobytes, so the grid is harvested where it lives and only the table and histogram are transferred.

Usage:

    python3 harvest_calib3.py <root directory> <output file> [<node> <node> ...]
    python3 harvest_calib3.py --combine [--fill[=<ygap>,<xgap>,<zext>]] <output file> <histogram file> ...
    python3 harvest_calib3.py --inspect <node directory>

With no node names the tree below the root is walked and every directory holding a hires configuration is taken to be a node.
The histogram files are the <output>.hist.npz written by the first form. With --fill, the empty cells of the table are
interpolated in the three coordinates in turn: along the temperature for gaps up to <ygap> times the temperature (0.25 by
default), which bridges the bands of a ladder of radiation fields; along log10 Sigma for gaps up to <xgap> dex (0.35), which
bridges the embedding surface densities of the grid; and along the contrast inside the populated range and <zext> levels beyond it
(1). The filled cells are written with a pixel count of 0.
"""


import os
import sys

import numpy as np
from astropy.io import fits

CFG = ["+hires/+hires.cfg", "+hires.cfg", "+hires/+results/+hires.cfg"]
RESULTS = ["+hires/+results", "+hires", "+results", ""]
TRUE_FINE = "nc.surfdens.r%sx0.rs3p0as.fits"
MASK = "nc.surfdens.bsl.r%sx0.rs3p0as.mask.fits"
SD = "hi.surface.density.r%s.fits"
TEMP = "hi.temperature.r%s.fits"

# The cells of the calibration. The steps are those of the first version, so the cells of new and old tables coincide, and only
# the ranges are wider: the grids at high radiation fields reach fitted temperatures above 40 K, and factors below 0.9 occur
# just outside structures.
XEDGE = np.round(np.arange(20.00, 24.5001, 0.05), 3)      # log10 of the surface density of the finest image
YEDGE = np.round(np.arange(3.00, 50.0001, 0.20), 2)       # fitted temperature at that resolution, K
ZEDGE = np.array([-1.00, -0.20, -0.10, -0.05, -0.02, 0.0, 0.02, 0.05, 0.10, 0.20, 1.00])
MINPIX = 20                                               # pixels below which a cell is not reported

# the factor, binned so that percentiles can be accumulated; the logarithmic step is that of the first version, 0.01355
FSTEP = np.log(6.0 / 0.9) / 140.0
FEDGE = np.exp(np.arange(np.log(0.5), np.log(8.0) + FSTEP, FSTEP))

AXES = ("x", "y", "z", "f")

# the filling done by --fill: widest gap along the temperature, as a fraction of it (the ladder of fields steps by 13.8 percent);
# widest gap along log10 Sigma (the embedding surface densities step by about 0.25 dex); levels added beyond the populated contrasts
YGAP = 0.25
XGAP = 0.35
ZEXT = 1


def tag(beam):
    return ("%.1f" % beam).replace(".", "p")


def read_cfg(path):
    rows = []
    for line in open(path):
        c = line.split()
        if len(c) < 5 or c[0] not in ("+", "-"):
            continue
        rows.append({"use": c[0] == "+", "wave": float(c[1]), "beam": float(c[2])})
    return rows


def find(node, name):
    for sub in RESULTS:
        p = os.path.join(node, sub, name) if sub else os.path.join(node, name)
        if os.path.isfile(p):
            return p
    return None


def resolve(node):
    """The files one node must supply, or the reason it cannot."""
    cfg = None
    for rel in CFG:
        p = os.path.join(node, *rel.split("/"))
        if os.path.isfile(p):
            cfg = p
            break
    if cfg is None:
        return {"ok": False, "why": "no hires configuration"}
    rows = [r for r in read_cfg(cfg) if r["use"]]
    if len(rows) < 3:
        return {"ok": False, "why": "fewer than three wavebands marked '+'"}
    rows.sort(key=lambda r: r["wave"])
    coarse, fine = rows[-1]["beam"], rows[0]["beam"]
    out = {"ok": True, "coarse": coarse, "fine": fine}
    want = {"sd_fine": SD % tag(fine), "sd_coarse": SD % tag(coarse),
            "t_fine": TEMP % tag(fine), "t_coarse": TEMP % tag(coarse),
            "true_fine": TRUE_FINE % tag(fine), "mask": MASK % tag(fine)}
    for key, name in want.items():
        p = find(node, name)
        if p is None:
            if key == "mask":
                out[key] = None
                continue
            return {"ok": False, "why": "missing %s" % name}
        out[key] = p
    return out


def getdata(path):
    return np.array(fits.getdata(path), float).squeeze()


def accumulate(info, counts, outside):
    """Add the pixels of one node to the histogram.

    Returns the number of valid pixels of the node. The pixels that fall outside the range of each axis are added to outside,
    one count per axis, so a pixel outside two ranges is counted under both axes but is excluded from the histogram once.
    """
    sd_f = getdata(info["sd_fine"])
    sd_c = getdata(info["sd_coarse"])
    t_f = getdata(info["t_fine"])
    tru = getdata(info["true_fine"])
    good = (sd_f > 0.0) & (sd_c > 0.0) & (t_f > 0.0) & (tru > 0.0)
    good &= np.isfinite(sd_f) & np.isfinite(sd_c) & np.isfinite(t_f) & np.isfinite(tru)
    if info["mask"] is not None:
        good &= getdata(info["mask"]) > 0.0
    if not good.any():
        return 0
    idx = []
    for v, edge in ((np.log10(sd_f[good]), XEDGE), (t_f[good], YEDGE),
                    (np.log10(sd_f[good] / sd_c[good]), ZEDGE), (tru[good] / sd_f[good], FEDGE)):
        idx.append(np.digitize(v, edge) - 1)
    ok = np.ones(idx[0].shape, bool)
    for a, (i, edge) in enumerate(zip(idx, (XEDGE, YEDGE, ZEDGE, FEDGE))):
        inside = (i >= 0) & (i < len(edge) - 1)
        outside[a, 0] += int((i < 0).sum())
        outside[a, 1] += int((i >= len(edge) - 1).sum())
        ok &= inside
    if ok.any():
        flat = np.ravel_multi_index(tuple(i[ok] for i in idx), counts.shape)
        np.add.at(counts.reshape(-1), flat, 1)
    return int(good.sum())


def percentiles(hist):
    """The 16th, 50th and 84th percentiles of the factor in one cell, from its histogram."""
    n = hist.sum()
    c = np.cumsum(hist)
    mid = 0.5 * (FEDGE[:-1] + FEDGE[1:])
    return [mid[int(np.searchsorted(c, q * n))] for q in (0.16, 0.50, 0.84)]


def empty():
    return np.zeros((len(XEDGE) - 1, len(YEDGE) - 1, len(ZEDGE) - 1, len(FEDGE) - 1), dtype=np.int32)


def save_hist(path, counts, outside, used, failed, npix, sources):
    """The histogram in sparse form, with its edges, so that histograms of several grids can be added exactly."""
    flat = np.flatnonzero(counts)
    np.savez_compressed(path, index=flat.astype(np.int64), count=counts.reshape(-1)[flat].astype(np.int64),
                        xedge=XEDGE, yedge=YEDGE, zedge=ZEDGE, fedge=FEDGE, outside=outside,
                        used=used, failed=failed, npix=npix, sources=np.array(sources))


def load_hist(path, counts, outside):
    """Add one histogram file to counts and outside, after checking that its cells are those of this script."""
    h = np.load(path)
    for key, edge in (("xedge", XEDGE), ("yedge", YEDGE), ("zedge", ZEDGE), ("fedge", FEDGE)):
        if h[key].shape != edge.shape or not np.allclose(h[key], edge, rtol=0.0, atol=1.0e-9):
            raise ValueError("%s: the %s cells differ from those of this script" % (path, key[0]))
    np.add.at(counts.reshape(-1), h["index"], h["count"])
    outside += h["outside"]
    return int(h["used"]), int(h["failed"]), int(h["npix"]), [str(s) for s in h["sources"]]


def measured(counts):
    """The cells holding at least MINPIX pixels: their indices, pixel counts and the 16th, 50th and 84th percentiles."""
    tot = counts.sum(axis=3)
    rows = {}
    for i, j, k in zip(*np.nonzero(tot >= MINPIX)):
        rows[(i, j, k)] = (int(tot[i, j, k]),) + tuple(percentiles(counts[i, j, k]))
    return rows


def _line(F, axis, coord, gap, relative):
    """Interpolate linearly across the interior gaps of F along one axis, for gaps up to gap (a fraction of the coordinate of the
    lower side when relative, an absolute width otherwise). F holds the three percentiles in its last dimension; NaN is empty."""
    G = np.moveaxis(F, axis, -2)
    ok = np.isfinite(G[..., 1])
    for idx in zip(*np.nonzero(ok.sum(axis=-1) >= 2)):
        js = np.nonzero(ok[idx])[0]
        for j1, j2 in zip(js[:-1], js[1:]):
            width = coord[j2] - coord[j1]
            if j2 - j1 < 2 or width > (gap * coord[j1] if relative else gap):
                continue
            t = ((coord[j1 + 1:j2] - coord[j1]) / width)[:, None]
            G[idx][j1 + 1:j2] = (1.0 - t) * G[idx][j1] + t * G[idx][j2]


def fill(rows, ygap, xgap, zext):
    """Interpolate the empty cells of the table between measured cells, in the three coordinates in turn.

    First along the temperature, at fixed surface density and contrast, for gaps up to ygap times the temperature. Grids computed
    for a ladder of radiation fields populate, at a given surface density, narrow bands of fitted temperature, one per field, and
    at fixed surface density and contrast the fitted temperature increases monotonically with the field, so this interpolates
    between neighboring fields. Then along the surface density, at fixed temperature and contrast, for gaps up to xgap in log10
    Sigma: the embedding surface densities of the grid are spaced by about 0.25 dex, and the uniform background of each forms a
    narrow band of cells. Then along the contrast, across interior gaps of any width, and beyond the populated levels by zext
    levels with the value of the nearest one, since a field has pixels a little below or above the contrasts the models produce
    at the same surface density and temperature.

    The filled cells are returned with a pixel count of zero, which marks them as interpolated: they enter the three-observable
    factor and the two-observable fallback, but no pixel-weighted average.
    """
    xc = 0.5 * (XEDGE[:-1] + XEDGE[1:])
    yc = 0.5 * (YEDGE[:-1] + YEDGE[1:])
    zc = 0.5 * (ZEDGE[:-1] + ZEDGE[1:])
    F = np.full((len(xc), len(yc), len(zc), 3), np.nan)
    for (i, j, k), r in rows.items():
        F[i, j, k] = r[1:]
    _line(F, 1, yc, ygap, True)
    _line(F, 0, xc, xgap, False)
    _line(F, 2, zc, 1.0e9, False)
    if zext > 0:
        has = np.isfinite(F[..., 1])
        for i, j in zip(*np.nonzero(has.any(axis=2))):
            ks = np.nonzero(has[i, j])[0]
            for q in range(1, zext + 1):
                if ks[0] - q >= 0:
                    F[i, j, ks[0] - q] = F[i, j, ks[0]]
                if ks[-1] + q < len(zc):
                    F[i, j, ks[-1] + q] = F[i, j, ks[-1]]
    new = {}
    for i, j, k in zip(*np.nonzero(np.isfinite(F[..., 1]))):
        if (i, j, k) not in rows:
            new[(i, j, k)] = (0,) + tuple(F[i, j, k])
    return new


def write_table(out, counts, outside, used, nfailed, npix, sources, failed=(), gaps=None):
    kept = int(counts.sum())
    rows = measured(counts)
    filled = fill(rows, *gaps) if gaps else {}
    with open(out, "w") as fh:
        fh.write("! harvest_calib3.py -- the correction factor against the three quantities the reconstruction reports\n")
        fh.write("!\n! %d nodes measured, %d not measured, %d pixels\n" % (used, nfailed, npix))
        fh.write("! %d pixels in the histogram; pixels outside the range of each axis (below, above):\n" % kept)
        for a, edge in enumerate((XEDGE, YEDGE, ZEDGE, FEDGE)):
            fh.write("!   %s from %9.4f to %9.4f: %d below, %d above\n"
                     % (AXES[a], edge[0], edge[-1], outside[a, 0], outside[a, 1]))
        if sources:
            fh.write("! combined from: %s\n" % " ".join(os.path.basename(s) for s in sources))
        if filled:
            fh.write("! %d measured cells and %d filled cells: along the temperature for gaps up to %.2f of it,\n"
                     % (len(rows), len(filled), gaps[0]))
            fh.write("! along log10 Sigma for gaps up to %.2f dex, along the contrast inside and %d level(s) beyond the populated\n"
                     % (gaps[1], gaps[2]))
            fh.write("! range;\n")
            fh.write("! the interpolated cells have a pixel count of 0\n")
        fh.write("! x = log10 of the surface density of the finest image, cm^-2\n")
        fh.write("! y = the temperature fitted at that resolution, K\n")
        fh.write("! z = log10 of the ratio of the finest to the coarsest surface density\n")
        fh.write("! f = the factor that restores the true surface density; f16 and f84 bracket its spread in the cell\n")
        fh.write("!\n!%9s %9s %9s %10s %9s %9s %9s\n" % ("x", "y", "z", "pixels", "f16", "f", "f84"))
        allrows = dict(rows)
        allrows.update(filled)
        for (i, j, k) in sorted(allrows):
            n, lo, md, hi = allrows[(i, j, k)]
            fh.write(" %9.3f %9.2f %9.4f %10d %9.4f %9.4f %9.4f\n"
                     % (0.5 * (XEDGE[i] + XEDGE[i + 1]), 0.5 * (YEDGE[j] + YEDGE[j + 1]),
                        0.5 * (ZEDGE[k] + ZEDGE[k + 1]), n, lo, md, hi))
        if failed:
            fh.write("\n! not measured:\n")
            for nd, why in failed[:50]:
                fh.write("!   %-46s %s\n" % (nd, why))
    print("\n %d nodes, %d pixels, %d in the histogram, %d measured and %d interpolated cells written to %s"
          % (used, npix, kept, len(rows), len(filled), out))
    for a, edge in enumerate((XEDGE, YEDGE, ZEDGE, FEDGE)):
        if outside[a].sum():
            print(" WARNING: %d pixels below and %d above the range of %s (%.4f to %.4f)"
                  % (outside[a, 0], outside[a, 1], AXES[a], edge[0], edge[-1]))


def inspect(node):
    print(" listing of %s" % node)
    for root, dirs, files in os.walk(node):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        print("   %s/" % os.path.relpath(root, node))
        for f in sorted(files):
            print("       %s" % f)
    info = resolve(node)
    print()
    if not info["ok"]:
        print(" NOT RESOLVED: %s" % info["why"])
        return
    for k in ("sd_fine", "sd_coarse", "t_fine", "t_coarse", "true_fine", "mask"):
        print("   %-10s %s" % (k, info[k] if info[k] else "not found"))


def combine(out, paths, gaps=None):
    counts, outside = empty(), np.zeros((4, 2), np.int64)
    used = nfailed = npix = 0
    for p in paths:
        u, f, n, _ = load_hist(p, counts, outside)
        used, nfailed, npix = used + u, nfailed + f, npix + n
        print(" %s: %d nodes, %d pixels" % (p, u, n))
    write_table(out, counts, outside, used, nfailed, npix, paths, gaps=gaps)


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "--inspect":
        inspect(sys.argv[2])
        return
    if len(sys.argv) >= 4 and sys.argv[1] == "--combine":
        args, gaps = sys.argv[2:], None
        if args[0].startswith("--fill"):
            gaps = [YGAP, XGAP, ZEXT]
            if "=" in args[0]:
                for n, v in enumerate(args[0].split("=")[1].split(",")[:3]):
                    gaps[n] = int(v) if n == 2 else float(v)
            args = args[1:]
        combine(args[0], args[1:], gaps)
        return
    if len(sys.argv) < 3:
        print(__doc__)
        return
    root, out, names = sys.argv[1], sys.argv[2], sys.argv[3:]
    if names:
        nodes = [os.path.join(root, s) for s in names]
    else:
        nodes = []
        for d, dirs, files in os.walk(root):
            dirs[:] = [x for x in dirs if not x.startswith(".") and not x.startswith("+")]
            if os.path.isfile(os.path.join(d, *CFG[0].split("/"))) or os.path.isfile(os.path.join(d, CFG[1])):
                nodes.append(d)
        nodes = sorted(set(nodes))

    counts, outside = empty(), np.zeros((4, 2), np.int64)
    used, failed, npix = 0, [], 0
    for nd in nodes:
        info = resolve(nd)
        if not info["ok"]:
            failed.append((nd, info["why"]))
            continue
        try:
            n = accumulate(info, counts, outside)
        except Exception as exc:
            failed.append((nd, "%s: %s" % (type(exc).__name__, exc)))
            continue
        used += 1
        npix += n
        if used % 100 == 0:
            print(" %d nodes, %d pixels" % (used, npix))

    write_table(out, counts, outside, used, len(failed), npix, [], failed)
    save_hist(out + ".hist.npz", counts, outside, used, len(failed), npix, [os.path.abspath(root)])
    print(" histogram written to %s.hist.npz" % out)


if __name__ == "__main__":
    main()
