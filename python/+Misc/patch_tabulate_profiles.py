#!/usr/bin/env python3
"""
patch_tabulate_profiles.py -- add a dense radial-profile output to
tabulate_recovery.py.

WHAT IT ADDS
------------
A new option --profile-out PATH.  When given, the script writes a second,
companion file containing the azimuthally averaged radial profile of every
node and band on the stamp's own pixel grid (3 arcsec bins, the same binning
getsf uses), out to --profile-rmax times R_conv.  Nothing about the existing
recovery table changes, and the profile costs no extra computation: the array
is already built by radial_profile() and is currently only sampled at the
truncation radii before being discarded.

WHY
---
Fitting a Bonnor-Ebert profile to a source is the most promising route to the
corrector's fourth axis: on the production grid, knowing the true R_BE gives a
node-level leave-one-out residual of 0.0445 dex against 0.1178 for the footprint
radius that the pipeline currently substitutes for it.

Measured on the grid, the accuracy of that fit is controlled by r/R_BE, the
fraction of the Bonnor-Ebert sphere lying inside the footprint.  For resolved
models (true R_BE above the 13.5 arcsec beam) the ratio of fitted to true R_BE
has a 16th-to-84th percentile range of

    r/R_BE  0.0-0.4 : 0.468 - 1.963
    r/R_BE  0.4-0.7 : 0.835 - 3.528
    r/R_BE  0.7-1.0 : 0.949 - 1.359
    r/R_BE  1.0-1.5 : 0.859 - 1.143
    r/R_BE  > 1.5   : 0.904 - 1.039

The fit is excellent once the fitted range reaches past the Bonnor-Ebert edge and
poor when it stops short.  The production grid has a median r/R_BE of 0.700, so
most training samples are fitted over too short a range.  The stamps contain the
profile further out; it is simply not written to disk.  Real getsf profiles have
the same reach: the non-deblended surface-density profile INTBS0161 extends to a
median of 105 arcsec against 27 arcsec for the deblended INTBD0161.

Sampling density is NOT the problem and is not changed: the grid already
provides a median of 12 tabulated radii inside the truncation radius, against
about 9 for a real getsf profile at 3 arcsec spacing.

OUTPUT FORMAT
-------------
    #      i   j   k  band     R_conv        r_as         I_r
where R_conv is the model's convolved extent in arcsec (repeated on every row of
that node and band for convenience), r_as is the bin-centre radius in arcsec,
and I_r is the azimuthally averaged intensity in the stamp's own units
(MJy/sr, or H2/cm^2 for the surface-density band).  The profile is NOT
background-subtracted; subtract the value at whatever radius you treat as the
rim, or fit a pedestal.

USAGE
-----
    python3 patch_tabulate_profiles.py tabulate_recovery.py

then, when running the tabulation, add for example

    --profile-out bes_model_grid_final2_profiles --profile-rmax 3.0
"""
import shutil, sys, os

ARGS = '''    ap.add_argument('--profile-out', default=None,
                    help='if given, also write the dense azimuthally averaged '
                         'radial profile of every node and band to this file, '
                         'on the stamp pixel grid.  Costs no extra computation; '
                         'needed so that a Bonnor-Ebert profile fit can be '
                         'carried past the footprint edge, which is what '
                         'determines whether the fitted truncation radius is '
                         'usable (see patch_tabulate_profiles.py).')
    ap.add_argument('--profile-rmax', type=float, default=3.0,
                    help='outer radius of the dense profile output, in units '
                         'of R_conv.  3.0 covers the range a real getsf '
                         'INTBS profile reaches; use a larger value to write '
                         'everything the stamp contains.')
'''

OPEN_BLOCK = '''    prof_out = None
    if args.profile_out:
        prof_out = open(args.profile_out, 'w')
        prof_out.write(
            '# dense azimuthally averaged radial profiles, one row per node, '
            'band and radial bin\\n'
            '#   R_conv : convolved extent of the model (arcsec), repeated on '
            'every row of that node and band\\n'
            '#   r_as   : bin-centre radius (arcsec), on the stamp pixel grid\\n'
            '#   I_r    : azimuthally averaged intensity, NOT background '
            'subtracted (MJy/sr, or H2/cm^2 for the SD band)\\n'
            '# %4s %3s %3s %5s %11s %11s %13s\\n'
            % ('i', 'j', 'k', 'band', 'R_conv', 'r_as', 'I_r'))

'''

WRITE_BLOCK = '''            if prof_out is not None:
                sel = r_centers <= args.profile_rmax * R_conv
                for rc, pv in zip(r_centers[sel], prof[sel]):
                    prof_out.write('%6d %3d %3d %5s %11.4e %11.4e %13.5e\\n'
                                   % (i, j, k, b, R_conv, rc, pv))
'''


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'tabulate_recovery.py'
    if not os.path.exists(path):
        sys.exit("file not found: %s" % path)
    src = open(path).read()
    if '--profile-out' in src:
        sys.exit("this file already appears to be patched; nothing done")

    m_args = "    ap.add_argument('--col-i', type=int, default=1)"
    m_open = "    nmiss = nnode = nmaskmiss = nabs = nfloor = ncfoot = 0"
    m_write = "            peak_val = prof[0]"
    m_close = "    out.close()"
    for mk in (m_args, m_open, m_write, m_close):
        if mk not in src:
            sys.exit("could not locate the anchor %r; the file does not look "
                     "like the expected tabulate_recovery.py" % mk[:40])

    shutil.copyfile(path, path + '.orig_profiles')
    src = src.replace(m_args, ARGS + m_args, 1)
    src = src.replace(m_open, OPEN_BLOCK + m_open, 1)
    src = src.replace(m_write, WRITE_BLOCK + m_write, 1)
    src = src.replace(m_close,
                      "    if prof_out is not None:\n"
                      "        prof_out.close()\n"
                      "        print('wrote %s  (dense radial profiles)'\n"
                      "              % args.profile_out)\n" + m_close, 1)
    open(path, 'w').write(src)
    print("patched   : %s" % path)
    print("backup    : %s.orig_profiles" % path)
    print("added     : --profile-out and --profile-rmax; dense profile writer")
    print("unchanged : the recovery table itself, and all existing options")
    import py_compile
    py_compile.compile(path, doraise=True)
    print("syntax    : OK")


if __name__ == '__main__':
    main()
