# Injection–recovery test — design notes

Running notes for the synthetic-core injection test that validates the mass
correction end-to-end (inject cores into real Herschel maps → getsf → correct →
compare to truth). Implementation detail kept here, out of the paper body.
Paper text: §3.8 "Planned end-to-end test" (label `sec:injection`).

## Purpose
- Hold-out validation only tests the *interpolation* (clean model observables in,
  true mass out, ~3%). It assumes the observables would be measured correctly.
- Injection test puts the *full extraction chain* in the loop: getsf detection,
  its interpolated-background construction, its SED fit, its size measurement.
- Specifically checks: (a) does getsf's interpolated background == our model's
  interpolated background? (b) does getsf's FWHM == our FWHM_bs? (c) coverage/
  flagging behaves. These can only be checked with a real getsf run.
- Scope caveat (in draft): validates method vs our own forward model through
  getsf; does NOT validate the forward model vs nature (uniform-cloud
  idealization, §"uniform-cloud idealization").

## Injection recipe (FINAL, corrected)
Per model, per band (160/250/350/500; +070), full cloud-scale, convolved,
resampled to common 3" grid:
- bs  = interpolated-bg-subtracted core image = I = raw - const(rim level)
- bsl = crater-bg-subtracted core image       = C = raw - crater_bg
- (raw = full RT core+cloud image; const = interpolated/rim cloud level;
  crater_bg = true depressed cloud under the core)

CORRECT injection (add the bs / interpolated-bg image):
    map_injected = map_observed + I_bs        (masked to r<=R_BE)
Derivation: at a site where observed background ~ model cloud level (matched
within sqrt2), additive = (sky with core) - (sky without core)
           = raw - const = I_bs.
The crater is included automatically (obs+I ~ raw = core + cratered cloud).
Equivalent "dig + add": dig depression (C - I), add crater core C →
  -(C-I) + C = +I. Same result. ONLY the bs image is needed; C cancels.

NOTE — earlier 2C - I was WRONG (sign error: depression is C - I, not I - C;
and the core-only image C must not be added alone). Do NOT use 2C - I.
Do NOT inject the bsl/crater image. Inject bs (interpolated) only.

Why this is right / does not pre-cancel the effect: getsf subtracts its own
interpolated background from the rim, recovering ~ the bs flux = M_SED3bs
(biased low by over-subtraction) → correction → M_BE. The over-subtraction is
produced by getsf's background step on the realistic (obs+I) map, exactly as
on real data. bsl/crater image is only the reference for the IDEAL correction
(M_SED3bsl), not used in injection.

## Masking
- Crater + core are strictly bounded by r <= R_BE (BE sphere truncated there;
  no core or crater signal beyond R_BE).
- Per-band core masks already produced by user: Sigma convolved to that band,
  clipped > 1e16 cm^-2, divided by itself → 1.0 inside convolved core, NaN→0
  outside. Binary 0/1 support at each band's resolution.
- Mask is ~1 pixel wider than the signal → no boundary step, no taper needed.
- Apply each band's mask to that band's (2C - I) stamp. Outside R_BE everything
  is defined zero; this removes the ±1e15 convolution-residual speckle (harmless
  in amplitude, ~7 dex below cloud, but masked for cleanliness).

## Placement rule (per candidate position, per model)
Accept a position only if ALL hold (use model's own R_BE disk for support):
1. stamp's R_BE disk fully inside the per-band COVERAGE mask + edge margin
   (no overhang). Use largest band (500 um) mask to be safe.
2. no stamp pixel overlaps any extracted-source FOOTPRINT (+ small margin)
   — so cores sit on genuine cloud, not on bg-subtraction residuals.
   [Note: original image == source-subtracted image OUTSIDE footprints, so with
    no-overlap enforced the two image types are equivalent; we use original.]
3. local observed Sigma_obs within sqrt(2) (±~0.15 dex) of the model's
   Sigma_cloud, measured over an aperture ~ core size.
4. minimum separation from already-injected cores (no blending).

- Iterate over models/(mass,Σ,size); for each, find all valid sites, place as
  many as fit. Report sites-found per model; flag models with too few.
- High-Σ models will have few sites (little matching sky; dense filaments are
  crowded with real sources) — expected, not an error. More realizations there
  if needed. Low/mid-Σ bulk (where correction matters) is well-sampled.

## Registration / units
- All Herschel bands resampled to identical 3" pixel grid, same Npix, fully
  aligned. Pixel coordinates only — NO astrometry/WCS needed.
- Model images already in MJy/sr, already convolved to each band's resolution,
  background-free additive stamps. Literal sum into the map.
- Model images have variable odd Npix, source on the central pixel; stamp at
  map[i-h:i+h+1, j-h:j+h+1], h=(Npix-1)/2.

## Outputs of injection module
- Injected maps (one per band).
- Truth table (one row per injected core): ID, pixel (i,j), Σ_obs at site,
  true M_BE, Σ_cloud, R_BE, FWHM (model). Keyed by ID.
- This truth table is consumed by the getsf-output reader (same reader reused
  for the real Aquila run) to compute recovered-vs-true corrected mass.

## First test (agreed)
- Reuse ONE subdivision-catalog model (known true M_BE, sits between grid nodes
  → genuine test). Inject once, run getsf, correct, compare to truth AND to the
  ~3% from feeding catalog observables directly. Agreement confirms getsf
  observables line up with model observables; any gap quantifies what the real
  chain adds.

## Files needed from user (cut-out development set, same cut-out region)
Field side (FITS, cut-out):
  - five band maps (160/250/350/500; +70 optional), MJy/sr
  - per-band coverage mask(s)
  - surface-density (Σ) map
  - extracted-source footprint image (note encoding)
Model side (one subdivision model to start; per band):
  - C (crater-bg-subtracted), I (interpolated-bg-subtracted), band mask
  - raw core+cloud optional (sanity checks)
User also has dust-temperature images (from Σ + 250 μm) — not needed for
injection; possible later diagnostic / truth-table column.

## Build order (once files arrive)
1. valid-placement map for the test model over the cut-out (visual).
2. before/after image at one position (injected core visible).
3. truth-table row.
4. then scale to many positions / many models → injected maps → getsf.
5. getsf-output reader → recovered-vs-true figure + statistics.

## Data received (GitHub repo amenshch/Aquila-cutout-test, fetched via codeload)
Field cut-out, all 975x975, 3"/pix, MJy/sr (Sigma in H2/cm^2), aligned:
  - bands: aquilaM2-<bbb>.image.resamp.zoom2[+suffix].fits  (070 no suffix;
    160 +159p9; 250 +169p7; 350 +93p0; 500 +37p0) + each .omask.fits
  - Sigma: hi.surface.density.r13p5.fits + .omask.fits
Model i4j2k5 (subdivision row 9; Sigma_cloud=3.39e22, R_BE=35.9", M_BE=0.226):
  - nc.<bbb>um.bs.<rRES>x0.rs3p0as.fits   = I (interpolated-bg) -> INJECT THIS
  - nc.<bbb>um.bsl.<rRES>x0.rs3p0as.fits  = C (crater-bg)       -> not needed for injection
  - nc.surfdens.bsl.<rRES>x0.rs3p0as.mask.fits = per-band core mask (0/1)
  band->resolution token: 070->r8p4 160->r13p5 250->r18p2 350->r24p9 500->r36p3
  (also 100->r9p4, 850->r36p3 present; not used here)
  model npix=77 (odd, centered), mask support radius ~32px=97" (beam-broadened).

## Status
- Field side + placement: DONE and validated on real data. Placement =
  coverage ∩ Sigma-omask ∩ edge ∩ R_BE/support-disk-fit ∩ |Δlog Sigma|<0.5log2.
  For model i4j2k5 at Sigma=3.39e22: ~2.5% of cut-out valid (traces the dense
  filament), tens of non-overlapping sites. inject.py: load_field,
  valid_placement_mask, choose_positions, add_stamp(map += I_bs*mask).
- Injection formula CORRECTED to map + I_bs (see recipe above); single-position
  before/after demo looks right (clean core, no spurious dip).
- PENDING: getsf source-footprint masks (getsf still running) → wire into
  `forbidden=` for no-overlap. Then multi-position / multi-model injection →
  getsf → getsf-output reader → recovered-vs-true figure + stats.
- Pipeline correctors validated on hold-out: crater 0.9% / observable-direct
  3.1% / two-step 2.4% / forward-fit 20% (cross-check). Catalog = cat4 schema
  (53 cols). fit_mass.py keyed on FWHM250bs; I_C,bs tested, NOT an axis.
