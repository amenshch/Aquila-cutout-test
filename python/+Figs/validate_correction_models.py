"""
validate_correction_models.py -- end-to-end validation of the level correction
on the radiative transfer model grid.

For each node the synthetic Herschel images are built from the density and dust
temperature profiles of the model, the conventional reconstruction is run to
give the surface density image an observer would have, the level correction is
applied to it, and both are compared with the true surface density convolved to
the same resolution.
"""
import sys
import numpy as np

sys.path.insert(0, sys.path[0] or ".")
from hires_forward import (Node, load_profile_archive, build_maps, PIX_ARCSEC,
                           WAVES, BEAMS)
from hires_algorithm import Options, run_hires
from hires_correct_levels import (levels, fit_level_temperature, planck,
                                  invert_planck, smooth, incoherent_amplitude)

BEAM = 13.5
FIT = (2, 3, 4, 5)                       # 160, 250, 350, 500 micron
LAM = [WAVES[k] for k in FIT]


def correct(smap, images, pix=PIX_ARCSEC, factor=1.5, max_ball=205.0,
            tmin=6.0, tmax=30.0, ntemp=240, noise_sigma=2.0, c_max=2.5):
    """Level correction applied in memory, as hires_correct_levels does."""
    diam, d0 = [], 2.0 * BEAM
    while d0 <= max_ball * 1.0001:
        diam.append(round(d0, 1))
        d0 *= factor
    valid = np.isfinite(smap) & (smap > 0)
    imgs = [smooth(images[k], np.sqrt(max(BEAMS[-1] ** 2 - BEAMS[k] ** 2, 0.0)),
                   pix) for k in FIT]
    per_band = [levels(im, diam, pix)[0] for im in imgs]
    lev_map, cirrus = levels(smap, diam, pix)
    long_i, long_w = imgs[-1], LAM[-1]
    t_fit = np.clip(invert_planck(long_i, smap, long_w), 3.0, 100.0)
    numer = cirrus.copy()
    tmed = []
    for k, d in enumerate(diam):
        t = fit_level_temperature([pb[k] for pb in per_band], LAM,
                                  tmin, tmax, ntemp)
        c = planck(t_fit, long_w) / np.maximum(planck(t, long_w), 1e-300)
        c = np.where(np.isfinite(c), c, 1.0)
        c = np.clip(c, 1.0, c_max)
        if noise_sigma > 0:
            sig = incoherent_amplitude([pb[k] for pb in per_band], valid)
            ref = per_band[-1][k]
            thr = noise_sigma * sig
            if np.isfinite(thr) and thr > 0:
                w = ref ** 2 / np.maximum(ref ** 2 + thr ** 2, 1e-300)
                c = 1.0 + (c - 1.0) * w
        numer = numer + c * lev_map[k]
        tmed.append(np.median(t[valid]))
    field = numer / np.maximum(smap, 1e-30)
    field = smooth(np.where(np.isfinite(field), field, 1.0), BEAM, pix)
    return smap * np.maximum(field, 1.0), np.array(tmed), diam


def main():
    d = np.genfromtxt("hires/ladder_scan.csv", delimiter=",", names=True)
    H = d["FWHM_pc"] * 206265.0 / 260.0
    sel = (H > BEAM) & (d["contrast"] >= 1.2)
    rng = np.random.default_rng(31)
    ids = rng.choice(d["ID"][sel].astype(int), 30, replace=False)
    arc = load_profile_archive("hires/profiles_final3.txt.gz")
    rows = []
    for n, nid in enumerate(ids):
        nd = Node(int(nid), "cats/bes_grid_final3.txt", prof_table=arc[int(nid)])
        half = min(max(4.5 * nd.r_out_as, 110.0, 1.15 * nd.r_cloud_as), 620.0)
        mp = build_maps(nd, half_as=half, subpix=1)
        c = mp["n"] // 2
        sh = run_hires(mp["images"], Options())
        sc, tmed, diam = correct(sh, mp["images"])
        tru = mp["sigma_13p5"]
        rows.append((nid, nd.contrast, H[d["ID"].astype(int) == nid][0],
                     sh[c, c] / tru[c, c], sc[c, c] / tru[c, c]))
        if (n + 1) % 10 == 0:
            print("  %d / %d" % (n + 1, len(ids)), flush=True)
    a = np.array(rows)
    np.save("model_validation.npy", a)
    print("n = %d resolved nodes with contrast >= 1.2" % len(a))
    for j, lab in ((3, "conventional"), (4, "corrected")):
        v = a[:, j]
        print("  %-13s peak recovered / true: median %.3f, "
              "16th %.3f, 84th %.3f" % (lab, np.median(v),
                                        *np.percentile(v, [16, 84])))


if __name__ == "__main__":
    main()
