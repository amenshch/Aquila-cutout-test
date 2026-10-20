#!/usr/bin/env python3
"""validate_calib3.py -- compare the corrections on the models, where the answer is known.

The simulated sky cannot decide between the corrections. It is built from the analytic prescription that hicorrect2 embodies, so it
agrees with that correction by construction wherever the prescription is right, and it carries its own superposition bias
everywhere else. The grid of radiative transfer models is the only ground truth. This script applies the corrections to each node
and measures what is left. Run it on grids that did not enter the calibration: the subdivision grid, whose parameters lie between
those of the main grid, or one rung of the ladder of radiation fields with the calibration and the temperature term fitted without
it.

For every node it forms

    uncorrected:            C = Sigma_fine ,
    hicorrect2, coarse:     C = Sigma_fine + ( f( Sigma_coarse ) - 1 ) Sigma_coarse ,
    hicorrect2, fine:       C = f( Sigma_fine ) Sigma_fine ,
    hicorrect2, reshaped:   C = Sigma_fine + D ,  D the coarse correction redistributed with the shape of the fine one ,
    hicorrect3:             C = f( Sigma_fine, T_fine, log10 ( Sigma_fine / Sigma_coarse ) ) Sigma_fine ,

with f the analytic relation of hicorrect2.py, and the last factor the harvested
table applied exactly as hicorrect3.py applies it. It reports the relative residual against the true surface density at the fine
resolution in three forms.

  Integrated over the mask, the corrected divided by the true surface density, its median and percentiles over the nodes.

  Radial profiles about the peak of the true image, stacked over the nodes by taking the median in each annulus. A stacked profile
  can hide rings, because the radius of a ring scales with the size of the core and differs from node to node.

  The core mass of each node: the sum over the mask, which covers the footprint of the core in the convolved true image, of the
  surface density less the median in a ring just outside the mask, divided by the same measurement on the true image. The
  reshaped correction D starts from ( f( Sigma_fine ) - 1 ) Sigma_fine and is rescaled three times so that, smoothed to the
  coarse resolution, it equals the coarse correction.

  The same core mass with the footprint misjudged, as it is in a real image whose background has to be estimated: the mask
  shrunk or grown by a fraction of its own equivalent radius, so that concentrated and nearly uniform models are misjudged
  alike, with the background taken just outside the misjudged footprint. The change of the mass from
  one footprint to another is an error that no correction removes, and it sets the scale against which the differences between
  corrections should be judged.

  The ring of each node, measured in units of its own size. The half radius r_h is where the true profile, less its outer level,
  falls to half its central value; the ring residual is the median residual between 0.5 and 1.5 r_h, the central one the median
  within 3 arcsec, and their difference is positive where the correction leaves a donut. Its median and percentiles over the nodes
  are reported, with the fraction of nodes where it exceeds 0.1.

Requires hicorrect2.py and hicorrect3.py in the same directory, and the harvested calibration table.

Usage:

    python3 validate_calib3.py <root directory> <calibration table> <output file> [<node> ...]
"""

import os
import sys

import numpy as np
from astropy.io import fits
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hicorrect2 import build_relation_t
from hicorrect3 import Calibration

BANDS = [160.0, 250.0, 350.0, 500.0]
CFG = ["+hires/+hires.cfg", "+hires.cfg", "+hires/+results/+hires.cfg"]
RESULTS = ["+hires/+results", "+hires", "+results", ""]
TRUE_FINE = "nc.surfdens.r%sx0.rs3p0as.fits"
MASK = "nc.surfdens.bsl.r%sx0.rs3p0as.mask.fits"
SD = "hi.surface.density.r%s.fits"
TEMP = "hi.temperature.r%s.fits"

PIX = 3.0
REDGE = np.arange(0.0, 240.0, 6.0)        # annuli about the center of the core, arcsec
KEYS = ("u", "h2c", "h2f", "h2r", "h3")
LABELS = ("uncorrected", "hicorrect2, coarse", "hicorrect2, fine", "hicorrect2, reshaped", "hicorrect3")
NRESHAPE = 3                              # iterations of the reshaping of the coarse correction
BGWIDTH = 5                               # width in pixels of the ring outside the mask where the background is measured
FOOTSHIFT = (-0.5, -0.25, 0.25, 0.5, 1.0) # footprint errors tested, in units of the equivalent radius of the footprint:
                                          # the mask shrunk (negative) or grown (positive) by that fraction of its radius


def tag(beam):
    return ("%.1f" % beam).replace(".", "p")


def read_cfg(path):
    rows = []
    for line in open(path):
        c = line.split()
        if len(c) >= 5 and c[0] in ("+", "-"):
            rows.append({"use": c[0] == "+", "wave": float(c[1]), "beam": float(c[2])})
    return rows


def find(node, name):
    for sub in RESULTS:
        p = os.path.join(node, sub, name) if sub else os.path.join(node, name)
        if os.path.isfile(p):
            return p
    return None


def resolve(node):
    cfg = None
    for rel in CFG:
        p = os.path.join(node, *rel.split("/"))
        if os.path.isfile(p):
            cfg = p
            break
    if cfg is None:
        return None
    rows = [r for r in read_cfg(cfg) if r["use"]]
    if len(rows) < 3:
        return None
    rows.sort(key=lambda r: r["wave"])
    coarse, fine = rows[-1]["beam"], rows[0]["beam"]
    out = {"coarse": coarse, "fine": fine}
    for key, name in (("sd_f", SD % tag(fine)), ("sd_c", SD % tag(coarse)), ("t_f", TEMP % tag(fine)),
                      ("true", TRUE_FINE % tag(fine)), ("mask", MASK % tag(fine))):
        p = find(node, name)
        if p is None and key != "mask":
            return None
        out[key] = p
    return out


def getdata(path):
    return np.array(fits.getdata(path), float).squeeze()


def reshape(template, dc, fine, coarse):
    """The coarse correction dc redistributed with the shape of a template at the fine resolution.

    The result D is positive, follows the template, and returns dc when smoothed from the fine to the coarse resolution:
    D <- D dc / (G * D), repeated NRESHAPE times from D = template, with G the Gaussian of FWHM (coarse^2 - fine^2)^1/2. The amount
    of correction in every coarse beam is therefore the calibrated one, and only its distribution within the beam changes.
    """
    sigma = np.sqrt(coarse ** 2 - fine ** 2) / 2.3548 / PIX
    d = np.maximum(template, 0.0)
    floor = 1.0e-6 * max(dc.max(), 1.0)
    for _ in range(NRESHAPE):
        d = d * dc / np.maximum(ndimage.gaussian_filter(d, sigma, mode="nearest"), floor)
    return d


def half_radius(tr, rr, good):
    """The radius at which the true profile, less its outer level, falls to half its central value, in arcsec."""
    edges = np.arange(0.0, rr[good].max() + PIX, PIX)
    prof = np.array([np.median(tr[good & (rr >= a) & (rr < a + PIX)]) if (good & (rr >= a) & (rr < a + PIX)).any()
                     else np.nan for a in edges])
    ok = np.isfinite(prof)
    if ok.sum() < 4:
        return np.nan
    level = np.nanmedian(prof[ok][-max(3, ok.sum() // 5):])
    rel = (prof - level) / max(prof[0] - level, 1.0e-30)
    below = np.nonzero(ok & (rel <= 0.5))[0]
    return edges[below[0]] if len(below) else np.nan


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return
    args = sys.argv[1:]
    root, table, out, names = args[0], args[1], args[2], args[3:]
    cal = Calibration(table)
    fitted, factor, tfit = build_relation_t(BANDS)

    if names:
        nodes = [os.path.join(root, s) for s in names]
    else:
        nodes = []
        for p, dirs, files in os.walk(root):
            dirs[:] = [x for x in dirs if not x.startswith(".") and not x.startswith("+")]
            if os.path.isfile(os.path.join(p, *CFG[0].split("/"))) or os.path.isfile(os.path.join(p, CFG[1])):
                nodes.append(p)
        nodes = sorted(set(nodes))

    prof = {k: [] for k in KEYS}
    integ, rings, cmass = [], [], []
    shifted = {k: [] for k in FOOTSHIFT}
    rfoot = []
    used = 0
    for nd in nodes:
        info = resolve(nd)
        if info is None:
            continue
        try:
            sf, sc = getdata(info["sd_f"]), getdata(info["sd_c"])
            tf, tr = getdata(info["t_f"]), getdata(info["true"])
        except Exception:
            continue
        valid = (sf > 0.0) & (sc > 0.0) & (tf > 0.0) & (tr > 0.0)
        valid &= np.isfinite(sf) & np.isfinite(sc) & np.isfinite(tf) & np.isfinite(tr)
        inmask = getdata(info["mask"]) > 0.0 if info["mask"] is not None else valid.copy()
        good = valid & inmask
        if good.sum() < 100:
            continue
        sf1, sc1, tf1 = np.where(valid, sf, 1.0e21), np.where(valid, sc, 1.0e21), np.where(valid, tf, 14.0)
        fc = np.interp(np.log10(sc1), np.log10(fitted), factor, left=factor[0], right=factor[-1])
        ff = np.interp(np.log10(sf1), np.log10(fitted), factor, left=factor[0], right=factor[-1])
        f3 = cal.three(sf1, tf1, sc1)[0]
        dc = np.where(valid, (fc - 1.0) * sc, 0.0)
        cur = {"u": sf, "h2c": sf + dc, "h2f": ff * sf,
               "h2r": sf + reshape(np.where(valid, (ff - 1.0) * sf, 0.0), dc, info["fine"], info["coarse"]),
               "h3": f3 * sf}
        iy, ix = np.unravel_index(int(np.argmax(np.where(good, tr, 0.0))), tr.shape)
        ny, nx = tr.shape
        gy, gx = np.mgrid[0:ny, 0:nx]
        rr = np.hypot(gy - iy, gx - ix) * PIX
        rows = {k: [] for k in cur}
        for a, b in zip(REDGE[:-1], REDGE[1:]):
            m = good & (rr >= a) & (rr < b)
            for k in cur:
                rows[k].append(np.median(cur[k][m] / tr[m]) - 1.0 if m.sum() >= 8 else np.nan)
        for k in cur:
            prof[k].append(rows[k])
        rh = half_radius(tr, rr, good)
        mc = good & (rr < 3.0)
        mr = good & (rr >= 0.5 * rh) & (rr < 1.5 * rh) if np.isfinite(rh) else np.zeros_like(good)
        if mc.sum() >= 1 and mr.sum() >= 8:
            rings.append([np.median(cur[k][mr] / tr[mr]) - np.median(cur[k][mc] / tr[mc]) for k in KEYS])
        integ.append([cur[k][good].sum() / tr[good].sum() for k in KEYS])
        ring = ndimage.binary_dilation(inmask, iterations=BGWIDTH) & ~inmask & valid
        if ring.sum() >= 20:
            mtrue = (tr[good] - np.median(tr[ring])).sum()
            if mtrue > 0.0:
                cmass.append([(cur[k][good] - np.median(cur[k][ring])).sum() / mtrue for k in KEYS])
                req = np.sqrt(inmask.sum() / np.pi)
                rfoot.append(req * PIX)
                dout = ndimage.distance_transform_edt(~inmask)
                din = ndimage.distance_transform_edt(inmask)
                for n in FOOTSHIFT:
                    foot = ((din > -n * req) if n < 0 else (dout <= n * req)) & valid
                    if foot.sum() < 4:
                        continue
                    dfoot = ndimage.distance_transform_edt(~foot)
                    bg = (dfoot > 0.0) & (dfoot <= BGWIDTH) & valid
                    if bg.sum() >= 20:
                        shifted[n].append([(cur[k][foot] - np.median(cur[k][bg])).sum() / mtrue for k in KEYS])
        used += 1
        if used % 100 == 0:
            print(" %d nodes" % used)

    for k in prof:
        prof[k] = np.array(prof[k])
    integ, rings, cmass = np.array(integ), np.array(rings), np.array(cmass)
    with open(out, "w") as fh:
        fh.write("! validate_calib3.py -- the corrections measured on the models\n!\n")
        fh.write("! %d nodes, table %s\n" % (used, os.path.basename(table)))
        fh.write("!\n")
        fh.write("! integrated over the mask, corrected divided by true: median, 16th and 84th percentiles over the nodes\n")
        for k, lab in enumerate(LABELS):
            q = np.percentile(integ[:, k], [16, 50, 84])
            fh.write("!   %-24s %.4f  %.4f %.4f\n" % (lab, q[1], q[0], q[2]))
        fh.write("!\n! core mass over the mask, less the median of a ring of %d pixels just outside it, divided by the same\n"
                 % BGWIDTH)
        fh.write("! measurement on the true image, %d nodes: median, 16th and 84th percentiles over the nodes\n" % len(cmass))
        for k, lab in enumerate(LABELS):
            q = np.percentile(cmass[:, k], [16, 50, 84])
            fh.write("!   %-24s %.4f  %.4f %.4f\n" % (lab, q[1], q[0], q[2]))
        fh.write("!\n! the same core mass with the footprint misjudged: the mask shrunk (negative) or grown (positive) by the\n")
        fh.write("! fraction s of its equivalent radius (median %.0f arcsec, 16th and 84th percentiles %.0f and %.0f), the\n"
                 % tuple(np.percentile(rfoot, [50, 16, 84])))
        fh.write("! background taken just outside the misjudged footprint, divided by the true mass measured with the correct\n")
        fh.write("! mask; median, then 16th and 84th percentiles, over the nodes\n")
        fh.write("!   %-24s" % "s" + "".join(" %22.2f" % n for n in FOOTSHIFT) + "\n")
        for k, lab in enumerate(LABELS):
            row = ""
            for n in FOOTSHIFT:
                if shifted[n]:
                    q = np.percentile(np.array(shifted[n])[:, k], [50, 16, 84])
                    row += "   %6.3f [%5.3f %5.3f]" % tuple(q)
                else:
                    row += " %22s" % "-"
            fh.write("!   %-24s" % lab + row + "\n")
        fh.write("!\n! ring minus center of each node (residual between 0.5 and 1.5 half radii minus residual within 3 arcsec),\n")
        fh.write("! %d nodes: median, 16th and 84th percentiles, and the fraction of nodes above 0.1\n" % len(rings))
        for k, lab in enumerate(LABELS):
            q = np.percentile(rings[:, k], [16, 50, 84])
            fh.write("!   %-24s %+.4f  %+.4f %+.4f   %.3f\n" % (lab, q[1], q[0], q[2], np.mean(rings[:, k] > 0.1)))
        fh.write("!\n! relative residual against the radius from the center of the core, median over the nodes\n")
        fh.write("!%9s" % "radius" + "".join(" %12s" % k for k in ("uncorrect", "h2 coarse", "h2 fine", "h2 reshaped",
                                                              "hicorrect3"))
                 + "\n")
        for i in range(len(REDGE) - 1):
            fh.write(" %9.1f" % (0.5 * (REDGE[i] + REDGE[i + 1]))
                     + "".join(" %12.4f" % np.nanmedian(prof[k][:, i]) for k in KEYS) + "\n")
    print("\n %d nodes written to %s" % (used, out))


if __name__ == "__main__":
    main()
