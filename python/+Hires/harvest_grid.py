#!/usr/bin/env python3
"""harvest_grid.py -- measure on the grid of radiative transfer models every quantity of the paper on the correction of the hires
surface density images that is measured on the models rather than derived analytically.

The output is a few kilobytes of text, so the grid is harvested where it lives and only the table is transferred.

What is measured, for each node:

  1. The decomposition of the deficit by scale, summed over the mask: the base image at the coarsest resolution, the true and the
     reconstructed increments, and the total at the finest resolution, each as a fraction of the true surface density; and the
     ratio that would be obtained by setting the base image to its true value and leaving the increments untouched.
  2. The signature of the temperature mixing at the central pixel: the relative residuals of the best single-component fit to the
     noiseless intensities in each waveband, their root-mean-square, the fraction of the true surface density that the
     single-component fit recovers, and the fraction that a free two-component fit recovers.
  3. The number of pixels at which a second component is accepted by a significance test, with the adopted calibration
     uncertainties.  This is expected to be zero everywhere and is the basis of the statement that the test never fires.
  4. The value of the base image at the central pixel, which is the last column of the table of the grid in the appendix.
  5. The change in the integrated corrected surface density when Gaussian noise is added to the waveband images.

It also reports, as a check, how closely its own single-component fit reproduces the base image written by hires.  The two are the
same quantity by construction and should agree to a few tenths of a percent.  A larger disagreement means that the assumptions
behind the correction do not describe the reconstruction, and it should be reported before the numbers are used.

Usage:

    python3 harvest_grid.py <root directory> <output file> [<node> <node> ...]
    python3 harvest_grid.py --inspect <node directory>

With no node names the tree below the root directory is walked and every directory holding a hires configuration is taken to be a
node.  The second form lists what a single node directory holds and what the script resolves in it, and writes nothing; run it
first if a node fails to resolve.

Requires hicorrect2.py in the same directory, from which the temperature relation, the effective attenuating column and the
correction are imported unchanged.

Surface densities in H2 cm^-2, intensities in MJy sr^-1, temperatures in K, beam sizes in arcsec.
"""

import io
import os
import sys

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hicorrect2 import (planck, cnu, temperature, attenuating_column, build_relation, TFLOOR, TSURF, SSTAR, GAMMA)

# ------------------------------------------------------------------ what to measure and how

CAL_ERR = {70.0: 0.20, 100.0: 0.20, 160.0: 0.20, 250.0: 0.10, 350.0: 0.10, 500.0: 0.10}
FIT_ERR = 0.20                  # relative weight of each waveband in the fit, as in hicorrect2
DCHI2 = 9.0                     # improvement in chi squared required to accept a second component
NOISE = {160.0: 1.5, 250.0: 0.50, 350.0: 0.30, 500.0: 0.15}      # MJy sr^-1, added at the native resolution
NREAL = 5                       # realizations of the noise
NT1 = 400                       # points of the temperature grid of the single-component fit
NT2 = 90                        # points of each temperature grid of the two-component fit
TLO, THI = 5.0, 40.0            # range of the fitted temperatures

OUTER = 0.90                    # fraction of the field radius beyond which the embedding cloud is measured

FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

CFG = ["+hires/+hires.cfg", "+hires.cfg", "+hires/+results/+hires.cfg"]
RESULTS = ["+hires/+results", "+hires", "+results", ""]
BAND = "nc.%03dum.r%sx0.rs3p0as.fits"               # waveband image at its own beam, in the node directory
BANDCONV = "nc.%03dum.r%sx0.rs3p0as.r%03d.fits"     # the same, already convolved by hires to the beam of another waveband
TRUE_FINE = "nc.surfdens.r%sx0.rs3p0as.fits"        # true surface density at the finest beam
TRUE_COARSE = "nc.surfdens.r%sx0.rs3p0as.r%sx0.fits"   # the same, further convolved to the coarsest beam
MASK = "nc.surfdens.bsl.r%sx0.rs3p0as.mask.fits"    # unity over the model field, zero outside it
BASE = "hi.surface.density.r%s.fits"                # reconstruction at a given beam


def tag(beam):
    """The string by which a beam size appears in a file name, for instance 36.3 -> 'r36p3' without the leading letter."""
    return ("%.1f" % beam).replace(".", "p")


def node_name(node, root):
    """'iIIjJJkKK' for the grid layout cSD_ii/M_jj/kk, the relative path otherwise.

    This is the convention of harvest_calibration.py and is reproduced here unchanged.
    """
    rel = os.path.relpath(node, root)
    parts = [q for q in rel.split(os.sep) if q not in (".", "")]
    if len(parts) >= 3:
        a, b, c = parts[-3], parts[-2], parts[-1]
        if a.lower().startswith("csd_") and b.lower().startswith("m_"):
            try:
                return "i%02dj%02dk%02d" % (int(a.split("_")[1]), int(b.split("_")[1]), int(c))
            except ValueError:
                pass
    if len(parts) >= 2 and parts[-2].lower().startswith("add"):
        try:
            return "add%02d" % int(parts[-1])
        except ValueError:
            pass
    return rel.replace(os.sep, "/")


def read_cfg(path):
    """Wavebands of a hires configuration: wavelength, beam and the flag that marks a waveband as used."""
    rows = []
    for line in open(path):
        c = line.split()
        if len(c) < 5 or c[0] not in ("+", "-"):
            continue
        rows.append({"use": c[0] == "+", "wave": float(c[1]), "beam": float(c[2]), "image": c[4]})
    return rows


def find(node, name):
    """An existing file of the given name, looked for in the node and in the places hires writes its results."""
    for sub in RESULTS:
        p = os.path.join(node, sub, name) if sub else os.path.join(node, name)
        if os.path.isfile(p):
            return p
    return None


def resolve(node):
    """Everything the script needs from one node, or a dictionary carrying the reason it could not be assembled.

    Every file is addressed by its exact name.  The node holds several variants of each image, distinguished by the qualifiers
    bg, bs, bgl and bsl, and a pattern loose enough to tolerate an unknown name would match those variants instead of the plain
    one, so no pattern is used.
    """
    cfg = None
    for rel in CFG:
        p = os.path.join(node, *rel.split("/"))
        if os.path.isfile(p):
            cfg = p
            break
    if cfg is None:
        return {"ok": False, "why": "no hires configuration", "node": node}
    rows = [r for r in read_cfg(cfg) if r["use"]]
    if len(rows) < 3:
        return {"ok": False, "why": "fewer than three wavebands marked '+' in %s" % cfg, "node": node}
    rows.sort(key=lambda r: r["wave"])
    coarse, fine = rows[-1]["beam"], rows[0]["beam"]
    out = {"ok": True, "node": node, "cfg": cfg, "bands": rows, "coarse": coarse, "fine": fine}
    for r in rows:
        p = find(node, BAND % (int(round(r["wave"])), tag(r["beam"])))
        if p is None:
            return {"ok": False, "why": "missing %s" % (BAND % (int(round(r["wave"])), tag(r["beam"]))), "node": node}
        r["path"] = p
    ref = int(round(rows[-1]["wave"]))
    conv = []
    for r in rows:
        p = find(node, BANDCONV % (int(round(r["wave"])), tag(r["beam"]), ref))
        conv.append(p)
    out["conv"] = conv if all(c is not None for c in conv) else None
    want = {"base": BASE % tag(coarse), "fine_image": BASE % tag(fine), "true_fine": TRUE_FINE % tag(fine),
            "true_coarse": TRUE_COARSE % (tag(fine), tag(coarse)), "mask_coarse": MASK % tag(coarse),
            "mask_fine": MASK % tag(fine)}
    for key, name in want.items():
        p = find(node, name)
        if p is None:
            if key.startswith("mask"):
                out[key] = None
                continue
            return {"ok": False, "why": "missing %s" % name, "node": node}
        out[key] = p
    return out


def getdata(path):
    return np.array(fits.getdata(path), float).squeeze()


def pixel_scale(path):
    """Pixel scale in arcsec, from the header if it carries one."""
    h = fits.getheader(path)
    for key in ("CDELT2", "CD2_2"):
        if key in h and float(h[key]) != 0.0:
            return abs(float(h[key])) * 3600.0
    return 3.0


def to_beam(image, beam_from, beam_to, pix):
    """An image convolved from one Gaussian beam to a larger one."""
    if beam_to <= beam_from:
        return image
    s = np.sqrt(beam_to ** 2 - beam_from ** 2) * FWHM_TO_SIGMA / pix
    return gaussian_filter(image, s, mode="nearest")


def single_fit(obs, waves, tgrid, model, relerr=None):
    """Single-component fit of a stack of spectral energy distributions.

    obs has the shape (npix, nband).  relerr is the relative uncertainty of each waveband, which sets the relative weight of the
    wavebands in the fit; a single value weights them all alike, which is what hicorrect2 assumes, and None weights the absolute
    residuals alike instead.  Returns the fitted temperature, the surface density implied by the longest waveband, the relative
    residuals, and chi squared with the adopted calibration uncertainties.
    """
    if relerr is None:
        wgt = np.ones_like(obs)
    else:
        wgt = 1.0 / (np.asarray(relerr, float) * np.maximum(obs, 1e-30)) ** 2
    num = np.einsum("pb,pb,tb->pt", wgt, obs, model)
    den = np.einsum("pb,tb,tb->pt", wgt, model, model)
    amp = num / np.maximum(den, 1e-300)
    chi = np.einsum("pb,ptb->pt", wgt, (obs[:, None, :] - amp[:, :, None] * model[None, :, :]) ** 2)
    j = np.argmin(chi, axis=1)
    p = np.arange(obs.shape[0])
    tfit = tgrid[j]
    fit = amp[p, j][:, None] * model[j, :]
    res = (obs - fit) / np.maximum(obs, 1e-30)
    cal = np.array([CAL_ERR.get(w, 0.20) for w in waves])
    chi2 = (((obs - fit) / (cal * np.maximum(obs, 1e-30))) ** 2).sum(axis=1)
    sigma = obs[:, -1] / (cnu(waves[-1]) * planck(tfit, waves[-1]))
    return tfit, sigma, res, chi2


def double_fit(obs, waves, tgrid2, model2):
    """Two-component fit of one spectral energy distribution, amplitudes linear and both temperatures on a grid.

    The columns of the model matrix are the intensity per unit surface density, so the fitted amplitudes are themselves the
    surface densities of the two components and their sum is the total.  Returns that total and chi squared with the adopted
    calibration uncertainties.
    """
    cal = np.array([CAL_ERR.get(w, 0.20) for w in waves])
    sig2 = (cal * np.maximum(obs, 1e-30)) ** 2
    best, bchi = None, np.inf
    for a in range(len(tgrid2)):
        for b in range(a + 1, len(tgrid2)):
            m = np.vstack([model2[a], model2[b]]).T
            w = 1.0 / sig2
            ata = (m * w[:, None]).T @ m
            atb = (m * w[:, None]).T @ obs
            try:
                x = np.linalg.solve(ata, atb)
            except np.linalg.LinAlgError:
                continue
            if np.any(x < 0.0):
                continue
            r = obs - m @ x
            c = float((r * r / sig2).sum())
            if c < bchi:
                bchi, best = c, (x, a, b)
    if best is None:
        return np.nan, np.inf
    x, a, b = best
    return float(x[0] + x[1]), bchi


def build_mask(true_coarse, node_info):
    """The area over which the sums are taken, and a description of it.

    This is the mask image written beside the model, which is unity wherever the signal of the model exceeds a hundredth of its
    peak at that resolution and zero elsewhere.
    """
    p = node_info.get("mask_coarse")
    if p is not None:
        return getdata(p) > 0.0, "mask image %s" % os.path.basename(p)
    return (np.isfinite(true_coarse) & (true_coarse > 0.0), "every positive pixel, no mask image found")


def decompose(mask, D36, D13, T36, T13):
    """The decomposition of the deficit by scale over one area, as the rows of Table 2 of the paper."""
    s_d36, s_d13 = D36[mask].sum(), D13[mask].sum()
    s_t36, s_t13 = T36[mask].sum(), T13[mask].sum()
    return {"base": s_d36 / s_t36, "base13": s_d36 / s_t13, "trueinc": (s_t13 - s_t36) / s_t13,
            "recinc": (s_d13 - s_d36) / s_t13, "total13": s_d13 / s_t13,
            "repaired": (s_t36 + (s_d13 - s_d36)) / s_t13, "npix": int(mask.sum())}


def harvest(node_info, fh, root):
    """Measure one node and write its block to the open output file.  Returns the summary line."""
    n = node_info
    name = node_name(n["node"], root)
    waves = [r["wave"] for r in n["bands"]]
    beams = [r["beam"] for r in n["bands"]]
    pix = pixel_scale(n["base"])

    D36, D13 = getdata(n["base"]), getdata(n["fine_image"])
    T36, T13 = getdata(n["true_coarse"]), getdata(n["true_fine"])
    mask, mdesc = build_mask(T36, n)
    d = decompose(mask, D36, D13, T36, T13)
    ng = d["npix"]

    # ---- the central pixel, taken as the peak of the true surface density at the finest resolution
    iy, ix = np.unravel_index(int(np.argmax(T13)), T13.shape)

    # ---- the waveband images brought to the coarsest beam
    cube_own = np.stack([to_beam(getdata(r["path"]), r["beam"], n["coarse"], pix) for r in n["bands"]], axis=-1)
    cube_hi = np.stack([getdata(p) for p in n["conv"]], axis=-1) if n["conv"] else None
    cube = cube_hi if cube_hi is not None else cube_own
    tgrid = np.linspace(TLO, THI, NT1)
    model = np.array([[cnu(w) * planck(t, w) for w in waves] for t in tgrid])
    flat = cube.reshape(-1, len(waves))
    ok = np.all(flat > 0.0, axis=1)
    tfit = np.full(flat.shape[0], np.nan)
    sfit = np.full(flat.shape[0], np.nan)
    chi2 = np.full(flat.shape[0], 0.0)
    tf, sf, res, c2 = single_fit(flat[ok], waves, tgrid, model, FIT_ERR)
    tfit[ok], sfit[ok], chi2[ok] = tf, sf, c2
    own = sfit.reshape(cube.shape[:2])

    # ---- how closely a plain single-component fit reproduces the base image written by hires
    # The two are the same quantity by construction.  The fit is repeated with each weighting of the wavebands and with each of
    # the two ways of bringing them to the common beam, so that a disagreement can be attributed.
    schemes = [("uniform relative", FIT_ERR), ("calibration relative", np.array([CAL_ERR.get(w, 0.20) for w in waves])),
               ("unweighted", None)]
    sources = [("convolved by hires", cube_hi), ("convolved here", cube_own)]
    checks = {}
    for sname, cb in sources:
        if cb is None:
            continue
        fl = cb.reshape(-1, len(waves))
        o = np.all(fl > 0.0, axis=1)
        for wname, re_ in schemes:
            v = np.full(fl.shape[0], np.nan)
            _, sv, _, _ = single_fit(fl[o], waves, tgrid, model, re_)
            v[o] = sv
            q = v.reshape(cb.shape[:2])
            g = np.isfinite(q) & (D36 > 0.0) & mask
            checks[(sname, wname)] = q[g].sum() / max(float(D36[g].sum()), 1e-300)
    check = checks.get(("convolved by hires", "uniform relative"), checks.get(("convolved here", "uniform relative"), np.nan))

    # ---- the signature at the central pixel
    obs = cube[iy, ix, :]
    t1, s1, r1, _ = single_fit(obs[None, :], waves, tgrid, model, FIT_ERR)
    tgrid2 = np.linspace(TLO, THI, NT2)
    model2 = np.array([[cnu(w) * planck(t, w) for w in waves] for t in tgrid2])
    s2, chi2_2 = double_fit(obs, waves, tgrid2, model2)
    rec1 = float(s1[0]) / T36[iy, ix]
    rec2 = float(s2) / T36[iy, ix]
    rms = float(np.sqrt(np.mean((100.0 * r1[0]) ** 2)))

    # ---- how often a second component would be accepted
    nfire = int(np.sum(chi2[ok] > DCHI2))

    # ---- the noise test
    fitted, factor = build_relation(waves)
    def corrected_sum(images):
        c = np.stack([to_beam(images[i], beams[i], n["coarse"], pix) for i in range(len(waves))], axis=-1)
        f = c.reshape(-1, len(waves))
        o = np.all(f > 0.0, axis=1)
        s = np.full(f.shape[0], np.nan)
        _, sv, _, _ = single_fit(f[o], waves, tgrid, model, FIT_ERR)
        s[o] = sv
        b = s.reshape(cube.shape[:2])
        m = np.isfinite(b) & (b > 0.0) & mask
        fa = np.interp(np.log10(np.maximum(b[m], 1.0)), np.log10(fitted), factor, left=factor[0], right=factor[-1])
        return float((fa * b[m]).sum())
    raw = [getdata(r["path"]) for r in n["bands"]]
    ref = corrected_sum(raw)
    rng = np.random.default_rng(1)
    dev = []
    for _ in range(NREAL):
        noisy = [raw[i] + rng.normal(0.0, NOISE.get(waves[i], 0.0), raw[i].shape) for i in range(len(waves))]
        dev.append(abs(corrected_sum(noisy) / ref - 1.0))
    noise_dev = 100.0 * max(dev)

    fh.write("\n! node %s\n" % name)
    fh.write("!   area: %s, %d pixels\n" % (mdesc, ng))
    fh.write("!   central pixel (y, x): %d %d, true surface density at %.1f arcsec %.4e cm^-2\n"
             % (iy, ix, n["coarse"], T36[iy, ix]))
    fh.write("!   plain single-component fit of the base image, relative to the one written by hires:\n")
    fh.write("!     %-22s %18s %18s\n" % ("weighting", "convolved by hires", "convolved here"))
    for wname, _ in schemes:
        a = checks.get(("convolved by hires", wname))
        b = checks.get(("convolved here", wname))
        fh.write("!     %-22s %18s %18s\n" % (wname, "%.4f" % a if a is not None else "not available",
                                               "%.4f" % b if b is not None else "not available"))
    fh.write("!   relative residuals of the single-component fit at the central pixel, percent:\n")
    fh.write("!   %s\n" % "  ".join("%.0f um %+6.2f" % (w, 100.0 * v) for w, v in zip(waves, r1[0])))
    return ("%-12s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %8.2f %9d %9.3f %9.4f"
            % (name, d["base"], d["base13"], d["trueinc"], d["recinc"], d["total13"], d["repaired"], rec1, rec2, rms, nfire,
               noise_dev, np.log10(max(D36[iy, ix], 1.0))))


HEADER = """! harvest_grid.py -- quantities measured on the grid of radiative transfer models
!
! T(Sigma_a) = %.4f + (%.4f - %.4f) (1 + Sigma_a / %.5e)^(-%.5f) K
!
! node        base image at the coarsest resolution, summed over the mask, as a fraction of the true surface density at the
!             same resolution
! base13      the same sum, as a fraction of the true surface density at the finest resolution
! trueinc     true increments, as a fraction of the true surface density at the finest resolution
! recinc      reconstructed increments, in the same units
! total13     reconstruction at the finest resolution, as a fraction of the true surface density there
! repaired    the same, with the base image set to its true value and the increments left untouched
! rec1        fraction of the true surface density recovered by a single-component fit at the central pixel
! rec2        the same, for a free two-component fit
! rms         root-mean-square of the relative residuals of the single-component fit at the central pixel, percent
! nfire       pixels at which a second component improves chi squared by more than %.1f with the adopted calibration
!             uncertainties
! noise       largest change of the integrated corrected surface density over %d noise realizations, percent
! logbase     base-ten logarithm of the base image at the central pixel, cm^-2
!
!%-11s %9s %9s %9s %9s %9s %9s %9s %9s %8s %9s %9s %9s
""" % (TFLOOR, TSURF, TFLOOR, SSTAR, GAMMA, DCHI2, NREAL, " node", "base", "base13", "trueinc", "recinc", "total13",
       "repaired", "rec1", "rec2", "rms", "nfire", "noise", "logbase")


def inspect(node):
    """List what a node directory holds and what the script resolves in it."""
    print(" listing of %s" % node)
    for root, dirs, files in os.walk(node):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        rel = os.path.relpath(root, node)
        print("   %s/" % ("." if rel == "." else rel))
        for f in sorted(files):
            print("       %s" % f)
    info = resolve(node)
    print()
    if not info["ok"]:
        print(" NOT RESOLVED: %s" % info["why"])
        return
    print(" resolved:")
    print("   configuration      %s" % info["cfg"])
    for r in info["bands"]:
        print("   %6.0f micron at %5.1f arcsec   %s" % (r["wave"], r["beam"], r["path"]))
    for k in ("base", "fine_image", "true_coarse", "true_fine", "mask_coarse", "mask_fine"):
        print("   %-18s %s" % (k, info[k] if info[k] else "not found"))


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "--inspect":
        inspect(sys.argv[2])
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
            # the directories hires itself writes begin with '+' and must not be walked into or taken for nodes
            dirs[:] = [x for x in dirs if not x.startswith(".") and not x.startswith("+")]
            if os.path.isfile(os.path.join(d, *CFG[0].split("/"))) or os.path.isfile(os.path.join(d, CFG[1])):
                nodes.append(d)
        nodes = sorted(set(nodes))
    lines, failed, blocks = [], [], io.StringIO()
    for nd in nodes:
        info = resolve(nd)
        if not info["ok"]:
            failed.append((nd, info["why"]))
            print(" %-40s skipped: %s" % (nd, info["why"]))
            continue
        try:
            line = harvest(info, blocks, root)
        except Exception as exc:                              # a node that cannot be measured must not stop the run
            failed.append((nd, "%s: %s" % (type(exc).__name__, exc)))
            print(" %-40s failed: %s" % (nd, exc))
            continue
        lines.append(line)
        print(line)
    with open(out, "w") as fh:
        fh.write(HEADER)
        for line in lines:
            fh.write(" %s\n" % line)
        if failed:
            fh.write("\n! not measured:\n")
            for nd, why in failed:
                fh.write("!   %-40s %s\n" % (nd, why))
        fh.write("\n! the detail of each node follows\n")
        fh.write(blocks.getvalue())
    print("\n %d nodes measured, %d not measured, written %s" % (len(lines), len(failed), out))


if __name__ == "__main__":
    main()
