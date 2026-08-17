#!/usr/bin/env python3
"""
patch_band_R.py -- replace the per-band footprint solver in
mass_correction_pipeline_4d.py.

WHAT IT CHANGES
---------------
Only the nested function _band_R() inside build_recovery_samples(), plus one
added module-level diagnostic dictionary BANDR_DIAG.  Nothing else in the
pipeline is touched.  The original file is copied to
mass_correction_pipeline_4d.py.orig_bandR before editing.

WHY
---
_band_R() solved

    p_b(R) = p_SD(r) * c_b(R),
    c_b(R) = sigma(O_b; k_low = 1/(2R)) / sigma(O_SD; k_low = 1/(2r)),

for the band footprint radius R_lambda by fixed-point iteration.  p_b(R)
decreases with R and c_b(R) increases with R, so the fixed-point map is a
decreasing function of R; wherever the magnitude of its derivative exceeds
unity the iteration necessarily oscillates.  Measured over the 26649 solver
calls made while building the production training set, the old iteration
oscillated in 71.9% of calls and failed its own convergence test
(|R_new - R_old| < 0.5 arcsec within 6 iterations) in 23.4%, after which it
silently returned the last iterate; in a further 5.5% the target level fell
below the whole tabulated profile and the radius was clamped to the last
tabulated point.  The result was that R_lambda(r) was NOT monotone in r in
80-91% of nodes, which is physically impossible, and the recovered flux
fractions f_rec,lambda inherited the same defect.

Because g(R) = p_b(R) - p_SD(r) c_b(R) is strictly decreasing, its root is
unique and a bracketed method is guaranteed to find it.  The replacement uses
Brent's method on [r, R_tab], with R_tab the last tabulated radius of that
band, and returns NaN where no root exists inside the tabulated profile
instead of fabricating one.

USAGE
-----
    python3 patch_band_R.py path/to/mass_correction_pipeline_4d.py

Run with no argument to patch ./mass_correction_pipeline_4d.py.
"""
import shutil, sys, os

NEW_BAND_R = '''    def _band_R(node, b, r):
        """Per-band footprint radius R_lambda, from the condition of equal source
        signal-to-noise ratio at the footprint edge (paper Sect. 3.2):

            p_b(R) = p_SD(r) * c_b(R),
            c_b(R) = sigma(O_b; k_low = 1/(2R)) / sigma(O_SD; k_low = 1/(2r)),

        with p_x(s) the band-x intensity profile normalised to its own centre.
        p_b(R) decreases with R while c_b(R) increases with R, so
        g(R) = p_b(R) - p_SD(r) c_b(R) is strictly decreasing and its root is
        unique.  Solved by bracketed root finding (Brent) on [r, R_tab], where
        R_tab is the last tabulated radius of band b.

        This replaces an unguarded fixed-point iteration whose map has negative
        derivative of magnitude greater than one over much of the grid.  That
        iteration oscillated in 72% of calls, failed its own convergence test in
        23%, and silently returned the last iterate or clamped to the end of the
        profile table, making R_lambda(r) non-monotone in r in 80-91% of nodes.
        Where no root exists inside the tabulated profile this version returns
        NaN, so the sample is dropped rather than assigned a fabricated radius.
        Outcome counts accumulate in the module-level dict BANDR_DIAG.
        """
        if not CUMULATIVE_FRAC:
            return r
        levSD = _level_SD(node, r)
        if not (levSD > 0):
            BANDR_DIAG['no_profile'] += 1
            return np.nan
        pa = _prof_norm(node, b)
        if pa is None:
            BANDR_DIAG['no_profile'] += 1
            return np.nan
        rtab, ptab = pa
        Rhi = float(rtab[-1])
        sSD = _sigma_beam(SD_BEAM, 1.0 / (2 * r))

        def _g(R):
            pb = float(np.interp(R, rtab, ptab))
            cb = _sigma_beam(BAND_BEAM[b], 1.0 / (2 * R)) / sSD
            return pb - levSD * cb

        if not (Rhi > r):
            BANDR_DIAG['table_too_short'] += 1
            return np.nan
        if _g(r) < 0.0:                     # root would lie below r: the band
            BANDR_DIAG['root_below_r'] += 1  # footprint cannot be smaller than
            return r                         # the Sigma-band one
        if _g(Rhi) > 0.0:                   # profile table does not reach the root
            BANDR_DIAG['table_too_short'] += 1
            return np.nan
        from scipy.optimize import brentq
        R = float(brentq(_g, r, Rhi, xtol=1e-3, rtol=1e-8, maxiter=200))
        BANDR_DIAG['converged'] += 1
        return R

'''

DIAG_DECL = '''
# Outcome counts of the per-band footprint solver _band_R(), accumulated over a
# call to build_recovery_samples().  Reset it before the call if you want counts
# for that call alone.
BANDR_DIAG = dict(converged=0, table_too_short=0, root_below_r=0, no_profile=0)

'''

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'mass_correction_pipeline_4d.py'
    if not os.path.exists(path):
        sys.exit("file not found: %s" % path)
    src = open(path).read()

    if 'BANDR_DIAG' in src:
        sys.exit("this file already appears to be patched (BANDR_DIAG present); "
                 "nothing done")

    marker_a = "    def _band_R(node, b, r):"
    marker_b = "    node, sig, conc, cfoot, mrec"
    if marker_a not in src or marker_b not in src:
        sys.exit("could not locate _band_R() or the end-of-block marker; the file "
                 "does not look like the expected mass_correction_pipeline_4d.py")
    a = src.index(marker_a)
    b = src.index(marker_b)
    if b < a:
        sys.exit("markers found in unexpected order; aborting")

    anchor = "def build_recovery_samples("
    if anchor not in src:
        sys.exit("could not locate build_recovery_samples(); aborting")

    shutil.copyfile(path, path + '.orig_bandR')
    out = src[:a] + NEW_BAND_R + src[b:]
    ai = out.index(anchor)
    out = out[:ai] + DIAG_DECL.lstrip('\\n') + out[ai:]
    open(path, 'w').write(out)

    print("patched   : %s" % path)
    print("backup    : %s.orig_bandR" % path)
    print("replaced  : %d characters of the old _band_R with %d of the new one"
          % (b - a, len(NEW_BAND_R)))
    print("added     : module-level BANDR_DIAG counter dictionary")
    import py_compile
    py_compile.compile(path, doraise=True)
    print("syntax    : OK")

if __name__ == '__main__':
    main()
