# Raw-background injection run (flatten=False) — manifest

Model: i1 j2 k4 (M_BE = 0.0800 M_sun), 18 copies at the SAME 18 positions as the
flattened run (truth.csv), so truth.csv stays valid and the only change vs the
first run is flatten ON -> OFF (raw Herschel cloud kept under each core).

Injection: map_injected = map_observed + I_bs * mask  (no inpainting/flattening).
I_bs = nc.<band>um.bs.<rRES>x0.rs3p0as.fits ; mask = nc.surfdens.bsl.<rRES>x0.rs3p0as.mask.fits
Positions = truth.csv (y_pix, x_pix) as 0-based numpy centers (= getsf 1-based - 1).
All outputs float32 with BZERO=0.0, BSCALE=1.0; field WCS header preserved.

## Files
  inj_070_raw.fits   (070 um, absorption for cold cores — correct)
  inj_160_raw.fits   (160 um)
  inj_250_raw.fits   (250 um)
  inj_350_raw.fits   (350 um)
  inj_500_raw.fits   (500 um)
  inj_surfdens_clean_raw.fits   (raw Sigma + clean model column; 13.5", H2/cm^2)

## getsf runs to do (completes the 2x2 with the flattened run)
  RAW + hires-from-bands : run your hires-Sigma reconstruction on the 5
        inj_*_raw bands, then getsf with detection on that 161 channel.
        -> the realistic/intrinsic bias (the paper number).
  RAW + clean-column     : getsf with detection/measurement using
        inj_surfdens_clean_raw.fits as the 161 channel.
        -> isolates over-subtraction on raw bg with NO temperature bias.

## What each comparison isolates (vs the flattened run already done / running)
  clean-column 161:  flattened (running now) vs raw (new) -> flattening artifact
  hires-from-bands:  flattened (done)        vs raw (new) -> flattening artifact
  at fixed bg, clean vs hires                              -> temperature contribution
  RAW + hires-from-bands                                   -> intrinsic surfdens bias

## Sanity (passed)
  injected-excess peak = model column peak = 4.317e21 (exact)
  sum(inj - raw) = 18 x sum(model column) (exact) -> additive injection clean.

After getsf: I match detections to truth.csv (+1 to truth coords), pull SED mass
(col121), dust temp (col123), surfdens mass (col53), apply ObservableCorrector
(M_SED3bs, Sigma, FWHM250bs) -> M_BE, and report recovered-vs-true for each cell.
