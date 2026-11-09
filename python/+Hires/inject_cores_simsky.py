#!/usr/bin/env python3
"""
inject_cores_simsky.py -- put radiative transfer models of Bonnor-Ebert spheres
into the stratified simulated sky made by make_simsky.py.

The sky already has a temperature gradient along every line of sight, set by the
total shielding column.  A model core carries a different stratification, one
computed by radmc-3d for its own density structure, so a field containing both
tests whether the correction reproduces two kinds of stratification at once, and
whether the column-only relation that describes the smooth sky alone still holds
when cores are present.

What is injected
----------------
The background-subtracted stamps of the grid, 'nc.<band>um.bsl.<res>x0.rs3p0as'
and 'nc.surfdens.bsl.r13p5x0.rs3p0as', which are the excess a model adds to the
cloud it is embedded in, including the depression it digs.  Injection is
therefore a plain addition, image += stamp, exactly as in inject.py.  The same
addition is applied to the true surface density, using the '000' stamp, so the
truth stays exact by construction.

The 070 stamp does not exist in the grid and the 100 micron waveband has no
stamp at all, so those two sky images are left unchanged.  Neither enters the
correction, which fits 160 to 500 micron, but the 100 micron image of the output
field is then not consistent with the rest and must not be used for anything.

Placement
---------
A model computed for an embedding column SD_emb may only be placed where the
local column of the sky matches it, within the symmetric window

    SD_emb / sqrt(2)  <=  local column  <  SD_emb * sqrt(2),

which is one rung of the grid wide and tiles the column axis exactly.  The local
column is the mean over a box of the model's own observed size.  Cores are kept
apart by a minimum separation and away from the edge by the stamp half-width.

Usage
-----
    python3 inject_cores_simsky.py \
        --width 40 \
        --catalog /Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260908_RT_BES_radmc3d/bes_model_grid_final3_catalog \
        --stamp-root /Users/amenshch/Astronomy/+SIMULATIONS_IMAGES/260908_RT_BES_radmc3d \
        --n-per-model 8 --out-tag cores

It writes, for each waveband, 'simsky6.w<W>.<band>.<res>.<tag>.fits'; the new
truth 'simsky2.w<W>.surfdens.r13p5.<tag>.fits'; a matching '+hires.cfg'; and a
catalogue of what was placed where.  Run hires on the waveband images exactly as
for the sky without cores.
"""
import argparse
import os
import sys

import numpy as np
from astropy.io import fits
from scipy import ndimage

PIX = 3.0
# waveband, beam token of the sky image, beam token of the model stamp
BANDS = [("070", "r8p4", None), ("100", "r9p4", None),
         ("160", "r13p5", "r13p5"), ("250", "r18p2", "r18p2"),
         ("350", "r24p9", "r24p9"), ("500", "r36p3", "r36p3")]
# 1-based columns of the grid catalogue, from its own header
COL = dict(n=1, i=2, j=3, k=4, SD_emb=5, T_BE=6, rho_BE=7, M_BE=8, R_BE_as=9,
           Td_avg=12, Td_emb=13, ICSDbsl=47, FWHMSDbsl=61)


def read_grid(path):
    """Read the catalogue of radiative transfer models."""
    rows = []
    for line in open(path):
        if line.lstrip().startswith("#") or not line.strip():
            continue
        f = line.split()
        try:
            rows.append(dict(tag="i%02dj%02dk%02d" % (int(f[COL["i"] - 1]),
                                                      int(f[COL["j"] - 1]),
                                                      int(f[COL["k"] - 1])),
                             i=int(f[COL["i"] - 1]),
                             SD_emb=float(f[COL["SD_emb"] - 1]),
                             M_BE=float(f[COL["M_BE"] - 1]),
                             R_BE_as=float(f[COL["R_BE_as"] - 1]),
                             Td_avg=float(f[COL["Td_avg"] - 1]),
                             Td_emb=float(f[COL["Td_emb"] - 1]),
                             FWHM=float(f[COL["FWHMSDbsl"] - 1])))
        except (ValueError, IndexError):
            continue
    return rows


def stamp_paths(tag, root):
    """File names of the stamps of one model, following the grid's tree."""
    i, j, k = (int(tag[1:3]), int(tag[4:6]), int(tag[7:9]))
    sub = os.path.join(root, "cSD_%02d" % i, "M_%02d" % j, "%02d" % k)
    out = {"000": os.path.join(sub, "nc.surfdens.bsl.r13p5x0.rs3p0as.fits")}
    for band, _, res in BANDS:
        out[band] = (None if res is None else
                     os.path.join(sub, "nc.%sum.bsl.%sx0.rs3p0as.fits"
                                  % (band, res)))
    return out


def read_stamp(path):
    d = np.array(fits.getdata(path), float).squeeze()
    if d.shape[0] % 2 == 0 or d.shape[1] % 2 == 0:
        sys.exit("ERROR: stamp %s has even size %s; the centre is not a pixel"
                 % (path, d.shape))
    return np.nan_to_num(d)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--width", type=int, default=40,
                   help="half-maximum width of the filaments of the sky to use")
    p.add_argument("--catalog", required=True,
                   help="bes_model_grid_final3_catalog")
    p.add_argument("--stamp-root", required=True,
                   help="root of the model tree, holding cSD_NN/M_NN/NN/")
    p.add_argument("--out-tag", default="cores",
                   help="tag appended to the output file names")
    p.add_argument("--n-per-model", type=int, default=8,
                   help="cores to place per selected model")
    p.add_argument("--n-models", type=int, default=12,
                   help="models to select, spread over mass and embedding column")
    p.add_argument("--mass-range", type=float, nargs=2, default=(0.1, 5.0),
                   help="range of true model mass M_BE in solar masses")
    p.add_argument("--fwhm-range", type=float, nargs=2, default=(15.0, 80.0),
                   help="range of observed size FWHMSDbsl in arcsec")
    p.add_argument("--min-sep", type=float, default=200.0,
                   help="minimum separation between cores, arcsec")
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)

    stem = "simsky6.w%d" % a.width
    sky = {}
    for band, res, _ in BANDS:
        f = "%s.%s.%s.fits" % (stem, band, res)
        with fits.open(f) as hd:
            sky[band] = (np.array(hd[0].data, float), hd[0].header.copy())
    with fits.open("%s.surfdens.r13p5.fits" % stem) as hd:
        truth = np.array(hd[0].data, float)
        thdr = hd[0].header.copy()
    ny, nx = truth.shape
    print("sky %s: %d by %d pixels, column %.3g to %.3g cm-2"
          % (stem, ny, nx, truth.min(), truth.max()))

    grid = [g for g in read_grid(a.catalog)
            if a.mass_range[0] <= g["M_BE"] <= a.mass_range[1]
            and a.fwhm_range[0] <= g["FWHM"] <= a.fwhm_range[1]]
    if not grid:
        sys.exit("ERROR: no model in the catalogue passes the mass and size cuts")
    # keep only models whose embedding column occurs in this sky
    lo, hi = np.percentile(truth, [2, 98])
    grid = [g for g in grid
            if g["SD_emb"] * np.sqrt(2) > lo and g["SD_emb"] / np.sqrt(2) < hi]
    print("%d grid models pass the cuts and match a column present in the sky"
          % len(grid))
    # spread the selection over mass
    grid.sort(key=lambda g: g["M_BE"])
    pick = [grid[int(round(t))] for t in
            np.linspace(0, len(grid) - 1, min(a.n_models, len(grid)))]

    placed = []
    occupied = np.zeros((ny, nx), bool)
    sep_pix = a.min_sep / PIX
    for g in pick:
        sp = stamp_paths(g["tag"], a.stamp_root)
        if not os.path.exists(sp["000"]):
            print("   %s: no stamp at %s, skipped" % (g["tag"], sp["000"]))
            continue
        st = {b: read_stamp(sp[b]) for b in ["000"] + [x[0] for x in BANDS]
              if sp.get(b)}
        half = (st["000"].shape[0] - 1) // 2
        # local column on the scale of the model's own observed size
        box = max(3, int(round(g["FWHM"] / PIX)))
        sloc = ndimage.uniform_filter(truth, size=box, mode="nearest")
        ok = ((sloc >= g["SD_emb"] / np.sqrt(2))
              & (sloc < g["SD_emb"] * np.sqrt(2)))
        ok[:half + 1, :] = ok[-half - 1:, :] = False
        ok[:, :half + 1] = ok[:, -half - 1:] = False
        ok &= ndimage.distance_transform_edt(~occupied) > sep_pix \
            if occupied.any() else ok
        ys, xs = np.nonzero(ok)
        if len(ys) == 0:
            print("   %s: no position with a matching column, skipped" % g["tag"])
            continue
        order = rng.permutation(len(ys))
        n = 0
        for idx in order:
            if n >= a.n_per_model:
                break
            y, x = int(ys[idx]), int(xs[idx])
            if occupied[max(0, y - int(sep_pix)):y + int(sep_pix) + 1,
                        max(0, x - int(sep_pix)):x + int(sep_pix) + 1].any():
                continue
            sl = (slice(y - half, y + half + 1), slice(x - half, x + half + 1))
            for band, _, res in BANDS:
                if res is None or band not in st:
                    continue
                sky[band][0][sl] += st[band]
            truth[sl] += st["000"]
            occupied[y, x] = True
            placed.append((g["tag"], y, x, g["M_BE"], g["R_BE_as"],
                           g["SD_emb"], g["FWHM"], g["Td_avg"],
                           float(sloc[y, x])))
            n += 1
        print("   %s: M_BE %6.3f Msun, R_BE %6.1f\", FWHM %5.1f\", "
              "SD_emb %.3g cm-2, placed %d" % (g["tag"], g["M_BE"],
                                               g["R_BE_as"], g["FWHM"],
                                               g["SD_emb"], n))

    def write(name, data, header):
        fits.PrimaryHDU(np.asarray(data, np.float32),
                        header).writeto(name, overwrite=True)
        with fits.open(name, mode="update") as hd:
            hd[0].header["BZERO"] = 0.0
            hd[0].header["BSCALE"] = 1.0

    cfg = "+hires.cfg.%s.w%d.%s" % ("simsky6", a.width, a.out_tag)
    with open(cfg, "w") as f:
        f.write("#\n# HIRES: CONFIGURATION PARAMETERS FOR HIRES\n#\n"
                "# PM  WAVE  BEAM  OFFSET   IMAGE   [MASKIMAGE]\n#\n")
        for band, res, _ in BANDS:
            nm = "%s.%s.%s.%s.fits" % (stem, band, res, a.out_tag)
            write(nm, sky[band][0], sky[band][1])
            beam = float(fits.getheader(nm).get("BEAM", 0.0))
            f.write("  %s   %s  %4.1f   0.0   %s\n"
                    % ("-" if int(band) < 160 else "+", band, beam, nm))
    write("%s.surfdens.r13p5.%s.fits" % (stem, a.out_tag), truth, thdr)
    print("written %s and the waveband images" % cfg)

    cat = "%s.%s.cores.txt" % (stem, a.out_tag)
    with open(cat, "w") as f:
        f.write("# radiative transfer cores injected into %s\n" % stem)
        f.write("# tag  y_pix  x_pix  M_BE(Msun)  R_BE(arcsec)  "
                "SD_emb(cm-2)  FWHM(arcsec)  Td_avg(K)  local_column(cm-2)\n")
        for r in placed:
            f.write("%-12s %5d %5d %10.4f %10.2f %12.4g %10.2f %8.2f %12.4g\n"
                    % r)
    print("written %s: %d cores placed" % (cat, len(placed)))


if __name__ == "__main__":
    main()
