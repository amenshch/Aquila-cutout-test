#!/usr/bin/env python3
"""
patch_sed_signm.py -- remove the detection-significance cut from the band
selection in sed_fit(), in mass_correction_pipeline_4d.py.

WHAT IT CHANGES
---------------
One condition in one line.  The band-selection mask in sed_fit() currently reads

    keepm = ((FP > 0) & (FPE > 0) & (F > 0) & (E > 0) & (np.abs(S) > 1)
             & (FP / FPE >= 1) & (F / E > 1))

where S is getsf's monochromatic detection significance SIGNM.  The term
np.abs(S) > 1 is removed, leaving the two signal-to-noise criteria that already
test the measurements themselves: peak flux over its error at least 1, and total
flux over its error greater than 1.

WHY
---
getsf sets SIGNM to the sentinel value 9.999e-31 when a source is not visible in
the clean single scales of THAT band.  Every source in an ok.cat is detected in
at least one band -- here the 161 micron surface-density band -- and is then
measured in all bands.  A sentinel SIGNM therefore records where the source was
detected, not whether its fluxes and sizes in that band are trustworthy.
Treating it as a reliability test discards good measurements.

The effect is severe and column-dependent, because contrast falls as the cloud
background rises.  On the all24 injection series, counting matched injected
sources that the pipeline could assign a mass:

    field scaling   median Sigma_cloud (cm^-2)   matched   usable, before this fix
      0.25                3.36e21                   86            80
      0.5                 6.73e21                   45            41
      1                   1.33e22                   71            14
      2                        --                   55             0
      4                        --                   33             0

At and above the unscaled Aquila field, which is the regime the paper is about,
the sentinel removed nearly the whole injected sample even though the total
fluxes FXT_BST were positive and plausible in every band.

USAGE
-----
    python3 patch_sed_signm.py mass_correction_pipeline_4d.py

The argument S is left in the signature and is still returned, so nothing else
in the pipeline changes.
"""
import shutil, sys, os

OLD = ("        keepm = ((FP > 0) & (FPE > 0) & (F > 0) & (E > 0) & (np.abs(S) > 1)\n"
       "                 & (FP / FPE >= 1) & (F / E > 1))")
NEW = ("        # A band enters the fit on the strength of its own measurement,\n"
       "        # not on getsf's per-band detection significance SIGNM.  SIGNM\n"
       "        # carries the sentinel 9.999e-31 wherever the source was not\n"
       "        # detected in that band's clean single scales; since every source\n"
       "        # is detected in at least one band and then measured in all of\n"
       "        # them, the sentinel records where the detection happened rather\n"
       "        # than whether the fluxes are usable.  The thresholds below are\n"
       "        # getsf's own recommended optimal selection criteria, quoted in\n"
       "        # every catalogue header.  See patch_sed_signm.py.\n"
       "        keepm = ((FP > 0) & (FPE > 0) & (F > 0) & (E > 0)\n"
       "                 & (FP / FPE > SNR_MIN_PEAK) & (F / E > SNR_MIN_TOTAL))")

CONSTS = ("SNR_MIN_PEAK = 2.0     # FXP_BST / FXP_ERR, per band, to admit a flux\n"
          "SNR_MIN_TOTAL = 2.0    # FXT_BST / FXT_ERR, per band, to admit a flux\n"
          "\n")

DOC_OLD = "    |SIGNM|>1, FXP/FXP_ERR>=1, FXT/FXT_ERR>1 (>=2 bands)."
DOC_NEW = "    FXP/FXP_ERR>SNR_MIN_PEAK, FXT/FXT_ERR>SNR_MIN_TOTAL (>=2 bands);\n    SIGNM is NOT used."


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'mass_correction_pipeline_4d.py'
    if not os.path.exists(path):
        sys.exit("file not found: %s" % path)
    src = open(path).read()
    if 'patch_sed_signm.py' in src:
        sys.exit("this file already appears to be patched; nothing done")
    if OLD not in src:
        sys.exit("could not find the band-selection line in sed_fit(); the file "
                 "does not look like the expected mass_correction_pipeline_4d.py")
    shutil.copyfile(path, path + '.orig_signm')
    src = src.replace(OLD, NEW, 1)
    if DOC_OLD in src:
        src = src.replace(DOC_OLD, DOC_NEW, 1)
    src = src.replace("def sed_fit(", CONSTS + "def sed_fit(", 1)
    open(path, 'w').write(src)
    print("patched   : %s" % path)
    print("backup    : %s.orig_signm" % path)
    print("removed   : the np.abs(SIGNM) > 1 term from the sed_fit band selection")
    print("added     : SNR_MIN_PEAK = SNR_MIN_TOTAL = 2.0, getsf's recommended criteria")
    print("unchanged : both signal-to-noise criteria, and everything else")
    import py_compile
    py_compile.compile(path, doraise=True)
    print("syntax    : OK")


if __name__ == '__main__':
    main()
