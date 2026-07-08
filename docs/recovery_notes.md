# Injection-recovery: first result (model i1 j2 k4, 18 copies)

Catalog: Aquila.s.sources.ok.injected.03=SED.00.cat (getsf 260628 + FITFLUXES,
detection on the 161 / hires-Sigma channel, SED fit over 070-500).
Truth: truth.csv (18 identical copies of i1 j2 k4). truth.csv coords are
0-based; getsf is 1-based -> add 1. All 18 matched within 0.6 px.

## Truth (model i1 j2 k4, row 014 of bes_model_params_catalog)
  M_BE        = 0.0800 Msun   (true mass)
  M_SED3bs    = 0.0716 Msun   (idealised interp-bg SED mass, 250-500; model bias 1.12)
  T_SED3bs    = 12.08 K   Td_avg = 11.98 K   Td_emb = 13.44 K
  R_BE = 18"   FWHM250bs = 26.9"   Sigma_cloud = 3.0e21

## Recovery (means over 18)
  getsf SED mass      = 0.0584   bias M_BE/<SED> = 1.37   (94% below M_BE)
  getsf surfdens mass = 0.0392   bias M_BE/<SD>  = 2.04
  getsf dust temp     = 12.6 K   (truth 12.0 K; slightly warm)

Realistic getsf SED bias (1.37) > idealised model bias (1.12):
  <SED>/M_SED3bs = 0.82  -> extra ~18% loss = getsf rim-background over-subtraction
  on the real cloud. Warm-T bias and crater/over-subtraction stack, as designed.

## Surfdens deficit decomposition (answers: aperture vs over-subtraction?)
  model bg-sub column (nc.surfdens.bs.r13p5) integrates to 0.0795 ~= M_BE
        -> the clean column has ~no intrinsic mass bias (1.01).
  getsf 161 footprint (FOOA ~ 50", r~25") captures 99% of that column
        -> aperture truncation is NEGLIGIBLE.
  col53 / (model-in-footprint) = 0.50
        -> the entire factor-2 deficit is BACKGROUND OVER-SUBTRACTION on the
           smooth/extended Sigma channel, not aperture loss.

Both mass types fail through the same mechanism (getsf rim background); the
Sigma image is hit far harder (~50% lost) than the SED (~18% beyond model bias).

## Setup (clean-model run) and what the 50% means -- CONFIRMED with getsf images
This first run injected CLEAN models onto a FLATTENED (flat) baseline (3.0e21)
with fluctuations inside the model support removed. Using getsf's own 161-band
products -- bg image Aquila.s.161.obs.cb.fits and bg-subtracted source
Aquila.s.161.obs.bs.fits -- the deficit is fully resolved:
  - integral of bs in the 161 footprint = 0.0391 == col53 (loop closed).
  - core mass in getsf input (cb+bs) above the TRUE 3.0e21 baseline = 0.086
    ~= model 0.0795 ~= M_BE -> the full core IS present in the data.
  - getsf interpolated bg under footprint = 3.65e21 = 0.65e21 (22%) ABOVE the
    true baseline; it traces the RAW cloud just outside the footprint (ring-out
    median 3.60e21), not the flattened under-core level.
  - that 22%-high bg removes 0.047 Msun = 54% of the core -> the WHOLE deficit.
=> The ~50% surfdens loss is a FLATTENING ARTIFACT (flat under core, raw cloud
   around it; getsf's 4-direction linear interpolation reads the raw rim). NOT
   the low-contrast/buried-envelope effect (that earlier hypothesis was WRONG),
   and NOT aperture (footprint holds 99%). The realistic run (raw bg everywhere)
   should largely remove it; surfdens bias 2.04 here is mostly setup-induced.
   getsf bg (per A.M.): linear interpolation of pixel values around a source in
   4 directions (2 axes + 2 diagonals), averaged, on the original image.

## Next rung: realistic extraction
Inject onto RAW Herschel background (keep fluctuations, no flat bg under cores),
then RECONSTRUCT the hires surfdens from the injected band images (not inject a
clean Sigma), then run getsf. This re-introduces (i) the -14% temperature-
structure bias into the surfdens channel and (ii) site-dependent over-subtraction
scatter from real cloud structure, and exercises the SED-fit step end to end.
The clean run isolates over-subtraction (SED 1.37, surfdens 2.04, both pure bg);
the realistic run layers temp bias + fluctuation scatter on top.

## Column map (merged catalog, 1-based)
  XCO_P=5 YCO_P=6 ; SED mass=121 ; dust temp=123 ; surfdens mass=53(=118)
  161 footprint FOOA=61 FOOB=62 ; per-band FXT_BST: 070=17 160=34 161=51
  250=68 350=85 500=102 ; FXP_BST: 070=15 160=32 161=49 250=66 350=83 500=100

## SED-fit (FITFLUXES) convention -- SETTLED, matches Konyves+15 (HGBS Aquila)
Fit bands: 160, 250, 350, 500 um (4 bands). 70 um is MEASURED and carried in the
catalog but IGNORED in the fit (command line: -1:070). So getsf SED mass/temp
correspond to the grid's 4-band columns:
  getsf TOTL_MASS (col121) <-> grid M_SED4bs ; DUST_TEMP (col123) <-> T_SED4bs
  -> apply ObservableCorrector keyed on (M_SED4bs, SD_emb, FWHM250bs) -> M_BE
     (NOT M_SED3bs; swap to *_4bs / *_4bsl columns).
  i1j2k4 truth: M_SED4bs=0.0709, T_SED4bs=12.12 K.

Weighting (PRE-fit, inside the LM chi^2 minimization -- NOT post-fit):
  fitfun minimizes  sum_b (model_b - flux_b)^2 / sigma_b^2 .
  sigma_b is built in fitfluxes as:
     fxterror = sqrt( fxterro^2 + (fxtbest * adderro)^2 )
  i.e. sigma^2 = sigma_meas^2 + (rel_calib * flux)^2  -- getsf measurement error
  (fxterro) in quadrature with the additional/calibration relative error. This IS
  1/sigma^2 weighting = exactly Konyves+15 ("weighted by 1/sigma^2_err" + elevated
  calibration uncertainties). "Unweighted" was a misnomer; it is fully weighted,
  just no extra weight layered on top. covar from the LM fit -> mass uncertainty
  -> RANDMASS draws (rightmost ~100 cols), so those inherit the same sigma. Keep
  RANDMASS in all runs (no parse cost; gives per-source mass error bars and lets
  noise be separated from environmental scatter).

Calibration (additional relative) errors -- USE THESE:
  command line  1:0.2,1:0.2,1:0.1,1:0.1,1:0.1
  = 070:0.20(ignored)  160:0.20 (PACS)  250/350/500:0.10 (SPIRE)
  Earlier runs used uniform 0.20 on all; change SPIRE 0.20->0.10 to match HGBS.
  Effect on bright injected cores (calib-dominated sigma): SPIRE 4x more
  constraining vs 160 -> fitted T slightly cooler, mass slightly higher, RANDMASS
  spread narrower. Modest on central mass, more visible on error bars; right
  direction for HGBS comparability. Dust-to-gas 0.01+-20%, ref opacity 10+-20%,
  distance 260 pc (exact) also propagate into the mass uncertainty.
  Set this convention for ALL FOUR injection cells + the real-Aquila run so every
  number shares one stated convention.

## 2x2 RESULTS (3 of 4 cells) -- CORRECTS the flattening diagnosis
Truth M_BE=0.0800, model M_SED4bs=0.0709, FWHM250bs=26.9". Match to truth.csv+1,
SED fit 160-500 (col121<->M_SED4bs, col123<->T_SED4bs), surfdens=col53,
FWHM250=geo-mean(col72,col73). Corrected via direct interp on grid
(M_SED4bs,SD_emb=3e21,FWHM250bs)->M_BE.

  cell                         surfdens_bias  SED_bias  T(K)  corrected/M_BE
  hiresSD  (flat + hires)          2.04         1.38    12.6      0.83
  modelSD  (flat + clean col)      1.07         1.38    12.6      0.83
  modelSD.raw (raw + clean col)    1.07         1.44    12.8      0.79
  raw + hires                       PENDING (4th cell)

FINDING 1 -- surfdens deficit is the HIRES-Sigma-FROM-BANDS reconstruction, NOT
flattening. Clean model column recovers 1.07 (94% of the column) on BOTH flat and
raw bg (1.07==1.07); only hires-from-bands loses half (2.04). The earlier
"flattening artifact, getsf bg 22% high" reading was WRONG: that cb=3.65e21 was
getsf's background measured ON THE HIRES IMAGE; over-generalized. On the actual
model column getsf's bg subtraction is fine. Mechanism (to confirm with 4th cell):
hires reconstruction yields a more extended / lower-contrast core than the true
column, so getsf's bg sits high on it. flat-vs-raw is negligible for surfdens.

FINDING 2 -- the correction UNDERCORRECTS by ~17-21% through the real getsf chain.
Hold-out was ~3% (model observables in). But getsf SED mass = 0.058 < model
M_SED4bs 0.0709: getsf's 4-direction linear background over-subtracts ~18% MORE
than the model's idealized interpolated background. Fed 0.058, the corrector
returns 0.066 not 0.080. This is exactly the test target ("does getsf interp-bg ==
model interp-bg?" -> NO). FWHM matches (26.0 vs 26.9"), so it is the background,
not size. Actionable: either make the model's interp-bg match getsf's 4-direction
linear interpolation, or recalibrate the correction on getsf-measured observables.

SED bias is robust (1.38-1.44) across cells; raw bg adds a few % more over-sub vs
flat. Temp recovered 12.6-12.8 K (model 12.12; warm, from including 160).

## getsf background algorithm -- REPLICATED and validated (getsf_bg.py)
getsf footprint bg per interior pixel: linear interpolation along 4 directions
(H, V, 2 diagonals) from the pixels just OUTSIDE the footprint; drop directions
with length >= 2x the shortest (crowding); combine survivors with a POWER-9 mean
(~= the max of the directional values). My replica (getsf_bg.py) reproduces getsf's cb to 0.02% rms (EXACT).

KEY: on a UNIFORM cloud the 4 func are equal and the power-9 mean reduces EXACTLY
to that constant -> getsf bg == the model's rim-constant subtraction (I_bs =
raw - const). So the schemes AGREE on uniform cloud; rebuilding the grid with
getsf's scheme would NOT change M_SED4bs (still 0.0709). The ~18% getsf-vs-model
offset is therefore a REAL-CLOUD FLUCTUATION effect: getsf's power-9 (max-biased)
value read from the structured rim sits ABOVE the local mean -> over-subtracts.
Present in the flattened run too (flattening smooths only INSIDE the footprint;
getsf reads OUTSIDE). => "use getsf's bg in the model" does NOT close the gap;
the injection test itself is the recalibration (getsf/model M_SED4bs ~0.82 here),
and the offset is expected to vary with environment = the bias-vs-Sigma surface.

## getsf bg -- EXACT replica complete (getsf_bg.py)
Final algorithm, validated to 0.02% rms / 0.08% max vs cb using getsf's exact
footprints (foots_injected_hiresSD):
  per pixel, 4 directions; just-outside bracket pixels (iadd=0); each endpoint =
  MEAN of outside pixels in a 3x3 box (nave=1); linear interp between endpoints;
  drop directions with length >= 2x shortest; POWER-9 mean of survivors; then
  circular median filter radius nradmf=beam/pix (~4 for 161), lower-middle median.
The 3x3 endpoint averaging (nave=1) was the missing piece that closed the prior
0.7% (single-pixel endpoints sat high at the core-elevated footprint edge).
Uniform-cloud equivalence and the ~18% real-cloud offset conclusion unchanged.

## 2x2 COMPLETE (all 4 cells)
  cell                     surfdens_bias  SED_bias  T(K)  corrected/M_BE
  flat + hires                 2.04         1.38    12.6      0.83
  flat + clean                 1.07         1.38    12.6      0.83
  raw  + clean                 1.07         1.44    12.8      0.79
  raw  + hires (PAPER)         2.16         1.45    12.8      0.78

- Surfdens deficit = hires-Sigma-from-bands reconstruction, BACKGROUND-INDEPENDENT
  (clean 1.07 both bg; hires 2.04/2.16 flat/raw). Flattening not the cause.
- SED bias 1.38(flat)->1.45(raw); same for clean/hires (SED uses bands). Raw adds
  ~5% more over-sub. Correction undercorrects 0.83(flat)->0.78(raw): getsf
  recovers M_SED4bs=0.055 vs model 0.0709 -> ~22% extra over-sub on real cloud.
- SD-vs-SED: a FAITHFUL column (clean) recovers 0.93*M_BE (7% low) > corrected SED
  0.78*M_BE (22% low). Surfdens-integration would be the better estimator IF the
  Sigma reconstruction preserved the column; standard hires-from-bands loses ~half.
- Fix path for SED undercorrection: extraction transfer function = getsf/model
  M_SED4bs vs Sigma (~0.78 at 3e21), via getsf_bg.py on real-cloud cutouts.

## hires-Sigma deficit FULLY DIAGNOSED -- temperature-driven (for the SD paper)
Decomposition of the surfdens deficit (raw+hires cell), model i1j2k4:
  true M_BE 0.080 (peak excess 4.32e21, model-implied core T = 11.6 K)
   -> hires MAP itself 0.054 = 0.67 M_BE   (peak suppressed to 0.41, mass spread out)
   -> getsf measures      0.039 = 0.49 M_BE (bg over-subtracts the now-extended core)
Radial: hires/model rises 0.42(6")->0.51(18")->0.67(40") = peak suppressed, wings filled.

ROOT CAUSE: color temperatures feeding Eq (mosurfden) are warm-biased at the core.
At the 18 core centers: T4(36.3",160-500)=15.2 K, T3(24.9",160-350)=14.7 K,
T2(18.2",160-250)=15.7 K -- vs model-implied core T = 11.6 K (and T_SED4bs=12.1).
The cold dense center is ~invisible in the 160-250 color (warm cloud + warm core
envelope dominate the LOS), so the fit reads ~15 K. Since Sigma = I/(B(T)k..) and
160um is near the Wien peak, a 3-4 K overestimate suppresses the high-res Sigma
increment ~5x (B(160,11.6)/B(160,15)=0.16-0.19); + smooth 36.3" base -> peak 0.41.

FIX LEVER = the core temperature in Eq (mosurfden). Avenues:
 (a) post-hoc grid transfer function (hires-integrated mass, Sigma, FWHM)->M_BE [easy];
 (b) sharper / less mismatched T for the high-res term [medium];
 (c) RT-model-informed core temperature -> re-derive hires Sigma MAP [the SD paper].
Inputs in hand to reproduce + test the fix: inj_*_raw bands, T2/T3/T4 maps
(thin.+hires.r{250,350,500}.cfg.00.tempers.fits), opacity k0=9.31 lam0=300 beta=2
eta=0.01 mu=2.8. NEXT: reproduce 0.67/0.41 in Python, then swap corrected T.

## hires fix TEST -- temperature is the lever; bg-subtracted core T recovers mass
Reproduced the hires reconstruction (Eqs mosurfden/superdens, Gaussian beams) from
the inj_*_raw bands + T2/T3/T4 maps. Standalone repro matches hires BEHAVIOR (peak
suppression, spreading) and T-sensitivity, but has a residual ~1.24x scale +
baseline offset vs the production pipeline (warm-T fraction 0.45 here vs 0.67
production) -- so use the RELATIVE recovery, which is robust:

  core T used in Sigma conversion        recovered core fraction (my repro)
  hires per-pixel color T (~15 K)              0.45   (cloud-contaminated LOS)
  bg-subtracted SED core T (12.1 K)            0.68   (OBSERVABLE: getsf provides it)
  true core T (11.6 K)                         0.75   (ceiling)

=> Cooling the core conversion temperature from the cloud-contaminated ~15 K toward
the cloud-removed ~12 K lifts recovered mass ~50%, monotonically. THE FIX:
in each core footprint, convert Sigma using the getsf BACKGROUND-SUBTRACTED core
SED temperature (cloud removed) instead of the hires per-pixel color T (cloud-
contaminated). Uses only quantities getsf already produces; RT grid can refine the
core T further (-> 11.6 K ceiling). Keeps hires resolution, fixes its column
accuracy -> the standalone hires-improvement paper.
NEXT: test this substitution in the production hires pipeline (exact normalization),
across the grid (mass, Sigma) to map the recovery.

## Can we get T_avg from temperatures? -- grid-wide test (293 models)
Test on i1j2k4 bands: bg-subtracted CORE color T = 11.97-12.05 K from ANY band pair
(160/250, 250/350, 350/500, 250/500) == Td_avg 11.98 K to <0.1 K, while TOTAL
(cloud+core) color T = 15-16 K. => warm bias is CLOUD CONTAMINATION, removable by
background subtraction; the multi-resolution color Ts (T2/T3/T4 ~15 K) do NOT trend
to T_avg, so resolution/band variation alone cannot give T_avg.

Grid-wide (T_SED4bs = bg-subtracted 4-band SED T vs Td_avg; M_SED4bs/M_BE = residual
mass bias AFTER perfect bg-subtraction):
  whole grid M_SED4bs/M_BE: median 0.906, 16-84% 0.80-1.16
  by Sigma: 3e21->0.825, 6e21->0.900, 1.2e22->1.019, 2.4e22->1.48(outlier),
            4.8e22->0.957, 9.6e22->1.028 ; T_SED4bs-Td_avg in -1.05..+0.19 K
  corr(residual, gradient Td_emb-Td_avg) = -0.13 (weak)

CONCLUSIONS:
- Background subtraction is the dominant, necessary fix (removes ~3-4 K cloud bias).
  bg-subtracted color T tracks T_avg to ~0.2-1 K (vs ~3-4 K contaminated).
- A perfectly mass-conserving T_avg is NOT obtainable from SED/color Ts alone:
  residual ~9% (median) mass bias remains, environment-dependent, NOT cleanly
  gradient-correlated -> needs the RT grid.
- TWO-TIER SCHEME: (1) in hires, bg-subtract bands before deriving T -> color T~T_avg
  -> recover bulk of map-level deficit at full resolution (self-contained hires fix);
  (2) RT grid ObservableCorrector mops up the residual ~10% environment-dependent bias.
- FLAG: Sigma=2.4e22 bin (M_SED4bs/M_BE=1.48, T 1K cold) is a sharp outlier -- likely
  the bias-crossing / min-core-mass region; check for grid artifact before trusting.

## GRID CHECK: Sigma=2.4e22 slice is a localized ARTIFACT (recompute it)
Per-Sigma bad fraction (|M_SED4bs/M_BE - 1| > 0.3):
  3.0e21:0.05  6.0e21:0.06  1.2e22:0.07  2.4e22:0.69(!)  4.8e22:0.00  9.6e22:0.00
Only Sigma=2.4e22 is broken (69% bad); higher/colder columns (4.8e22, 9.6e22) are
clean -> NOT a steep-gradient effect, a computation artifact at that one Sigma.
Signature: T_SED4bs collapses to 5-7 K (vs Td_avg ~8 K), worst at high j; cold T
inflates Sigma -> M_SED4bs/M_BE up to 9. ACTION: recompute/inspect the 45 models at
SD_emb=2.4e22 (RT or SED-fit step) before trusting that slice in the correction
surface. Rest of the grid sound. (Excluding outliers, that bin's median ~1.09.)

## Background subtraction method for hires: USE getimages (Men'shchikov 2017)
getimages = the right tool, unifies the user's option 2 (in-hires median bg) and
option 3 (clean peak clipping): multiscale median filtering, N windows from 2*O_lam
to 4*X_lam, MINIMUM across windows (Eq2) + iterations (Eqs4-5) = a lower-envelope
estimator that clips peaks WITHOUT riding up under sources. One pass, no prior
extraction, single parameter X_lam. Better than option 1 (getsf x2 + power-9
over-sub). Recipe: getimages -> bg-subtracted bands -> derive T -> color T ~ T_M ->
corrected hires Sigma, all before getsf.
Caveat (paper Fig3 / 2016 App B): bg can OVERestimate inside very extended/overlapping
sources (true bg under a source is uncertain) -> mild over-subtraction; but it's ~10%
for compact cores and acts on bands->T (small T error) not directly on the column.
NEXT TEST: implement getimages multiscale-median bg, apply to inj bands, derive core
color T, check ~T_M=11.98 K; compare vs getsf cb and true-field ceiling.

## BG-method test: estimated background recovers T_M (and fast methods work)
Core color T (350/500, bands convolved to common 36.3"), target T_M=11.98 K, 8 cores:
  none (current)            15.12 K   (warm, cloud-contaminated)
  true field (ceiling)      11.97 K
  getimages-style (mscale median) 12.14 K   1x speed
  grey opening (fast)       12.59 K   164x FASTER than getimages
=> Background subtraction with an ESTIMATED bg recovers T_M to <1 K (vs 3-4 K bias),
   confirming the in-hires fix is practical. getimages best (0.2 K) but slow (as A.M.
   found); fast morphological opening within ~0.6 K. Opening under-subtracts slightly
   (lower-envelope -> residual warm cloud -> marginally warm T); residual is RT-grid-
   correctable. RECOMMENDATION for hires (speed matters): fast lower-envelope bg --
   multiscale grey opening (O(N) per scale), OR reuse getsf's multiscale-decomposition
   large-scale component (already in-pipeline, built to replace slow getimages).
Caveat: 8 cores, light getimages (M=2,3 windows); full getimages would tighten;
opening radius (~X_lambda) needs tuning.

## A.M.'s separation scheme VALIDATED (step 2 test) + complete hires algorithm
TEST: hires on cleanly-separated sources-only bands (inj-field), per-pixel 4-band MBB:
  <M_core> = 0.0735 = 0.92 M_BE (= 1.04 M_SED4bs), vs current blended hires 0.054=0.67.
=> separating sources from bg in each band BEFORE hires recovers the cloud-contamination
   loss (0.67 -> 0.92); residual ~8% = single-T-fit gap -> RT grid. (Earlier 0.0007 was
   a bug: dropped eta=0.01 -> x100 mass error; fixed.)

hires algorithm (where T enters):
  (1) per-pixel single-T MBB fits at several resolutions -> T4(36.3",160-500),
      T3(24.9",160-350), T2(18.2",160-250).
  (2) per-resolution column D(O_lam,T_j) = I_lam / [B_lam(T_j) kap_lam eta mu mH].
      <-- ONLY place T enters; at a core T_j is the cloud+core BLEND (~15 K) -> suppressed.
  (3) multiscale combine: D_hires = D(O500,T4) + sum_lam max_j[ dD(O_lam,T_j) ],
      dD = D - G*D (unsharp increment to next-coarser res), positive-clip at 160.

Why separation fixes it: stages (2)-(3) are LINEAR in I at fixed T, so hires(I_src)+
hires(I_bg)=hires(I_total) IF T fixed. Nonlinearity is only in stage (1) T-derivation.
Separating first -> hires(src bands) derives T from cloud-removed emission (T~T_M~12 K)
-> correct source column; hires(bg bands) -> cloud column; SUM = Sigma_src(T_M)+
Sigma_cloud(T_cloud) = true total (columns add, each at its own correct T). No "paint
T-image" step needed -- two-T-in-one-pixel emerges automatically.
Caveats: (a) 0.92 is CEILING (clean separation); real getsf per-band separation leaves
residual cloud -> small residual warm bias -> RT-grid-correctable; (b) per-band
separations must use consistent footprints/scales across bands before adding.

## OUTSKIRT PROBLEM (A.M.'s experience) -- diagnosed & fixed via simulation
Problem: hires on sources-only (per-band bg-subtracted) bands gives bad results at
source outskirts. MECHANISM (2D mass-conserving sim, cloud+Gaussian core, Tc=12 Tcl=18):
independent per-band bg subtraction truncates the source at a band-dependent radius
(sharp 160 vanishes first, broad 500 persists). At outskirts the bg-subtracted SED
RISES toward long wavelengths ([160,250,350,500]=[0.03,0.21,0.47,0.66] at r=40") ->
a free single-T fit reads it as very cold -> spurious huge column.
  M/M_true, per-pixel FREE-T fit: 1.66 (indep footprints), 3.36 (common) -- footprint
    tuning does NOT help (common is worse: keeps more low-S/N outskirt pixels).
  M/M_true, FIXED source T=T_M:   1.01 (indep), 1.02 (common) -- blowup ELIMINATED.
FIX: do NOT float temperature per pixel on sources-only bands. Determine source T from
the high-S/N core (bg-subtracted color T is robust & ~uniform = T_M) and hold it fixed
(or smooth, frozen at low S/N) in the hires column reconstruction
D(O_lam,T_M)=I_src,lam/[B_lam(T_M) kap eta mu mH]. Each increment well-behaved; broad-
band outskirt flux adds small increments at fixed T instead of triggering a cold fit.
Justified: source is ~isothermal at T_M once cloud removed; per-pixel T map is overfit
& harmful here. Consistent w/ scheme: BACKGROUND bands keep per-pixel T (have cloud
floor -> fit OK); SOURCE bands use fixed T_M. Figure: outskirt_problem.png.

## ALGORITHM: smooth source T image for hires (normalized convolution)
Need T in EVERY footprint pixel, well-behaved, from reliable fits near the peak only.
ALGORITHM:
 1. Per-pixel SED fit T_fit(x,y) over footprint; mark pixel INVALID if any band <=0.
 2. Weight w(x,y) = SOURCE BRIGHTNESS (e.g. source 250um flux, or min-band S/N); 0 if
    invalid. *** NOT inverse fit-variance *** -- spurious rising-SED outskirt pixels are
    tightly (confidently) fit by a COLD blackbody -> small formal error -> high weight ->
    they would drag T_smooth cold (tested: inverse-variance weight gave T~9.5K, M=2.9x).
    Brightness weight emphasizes the mass-bearing peak (mass-weighting proxy).
 3. Normalized convolution: T_smooth = G_sig*(w*T_fit) / G_sig*(w), sig ~ footprint scale.
    -> reproduces measured T at peak (high w), smoothly extends it into outskirts (w~0);
    well-behaved in all pixels. Outskirt T has ~no mass leverage (I_src->0 there).
 4. Use T_smooth in hires column reconstruction.
RESULT (2D sim, vs true mass): free-T per-pixel 1.72 (blows up); T_smooth 0.86; const
T_M 0.97. T_smooth recovers correct 12.2 K smoothly across footprint, NO garbage.
Residual ~10-15% is bg-subtraction (footprint) quality, not the T image; RT-grid-correctable.
Gradient core: T_smooth ~ peak-weighted mean (large kernel washes gradient); mass still
fine (0.83). Smaller kernel captures more gradient near peak if wanted.
Figures: Tsmooth_algorithm.png, outskirt_problem.png.

## Does T-smoothing need source positions? NO -- global formulation keeps it in hires
A.M. concern: per-footprint T-smoothing needs source centers -> would force a prior
getsf detection run -> method not self-contained in hires.
RESOLUTION: per-FOOTPRINT smoothing would need positions, but the brightness-weighted
normalized convolution can run GLOBALLY (one fixed kernel over the whole separated
image) with NO source list:
  - weight w = bg-subtracted brightness -> automatically high on sources, ~0 elsewhere;
  - validity mask = "any band <=0" -> per-pixel;
  - empty regions get a smoothed T but carry no mass (I_src~0) -> harmless.
TEST (2 sources, different T & size, positions NEVER used in construction):
  T=12,24": T_smooth@ctr=12.4 K, M/Mtrue=0.85 ; T=15,40": 15.4 K, 0.89.
  One global op + single kernel gave each source its OWN correct T.
=> method stays INSIDE hires, pre-detection. Only getsf input = per-band source/bg
   separation (detection-independent, per A.M.'s scheme step 1). Single full getsf
   extraction happens AFTER, on the corrected surfdens.
Interpretation: G*(w*T_fit)/G*(w) with w=brightness is a locally MASS-WEIGHTED smoothed
T = the mass-averaged T_M that conserves mass (2016). Not ad-hoc.
Caveats: single global kernel ok for moderate size range; multiscale/adaptive sigma
(still position-free) for wide range; strongly overlapping sources of very different T
are the hard case -> fallback = one getsf run for positions, then per-source T (option 1).

## Practical implementation of the T_smooth convolution
Normalized convolution = two ordinary convolutions + a divide, using EXISTING conv code:
  num = conv(w*T_fit, K);  den = conv(w, K);  T_smooth = num/den  (where den>floor).
Define K as a kernel IMAGE and feed it to the same routine used for Herschel kernels.
KERNEL SHAPE IS IRRELEVANT (denominator divides out its normalization): tested Gaussian
s=30", top-hat r=40", heavy-tail(15,60") -> all give same T_smooth (12.5/15.4) & mass
(0.80/0.88). => use a plain GAUSSIAN; no special non-Gaussian kernel needed (unlike the
physical Herschel PSF case).
SINGLE FIXED KERNEL SUFFICES for moderate size range: explicit multiscale (conv at
15/30/60", per-pixel pick finest scale with den>thr) gave 12.8/15.7, 0.77/0.88 --
essentially identical. Multiscale only needed "in principle"; if ever wanted, it's just
a few conv passes + per-pixel finest-reliable-scale selection (still existing conv code,
no spatially-variant machinery).

## Weight definition + position-dependent noise (tested)
w(x,y) = max(S_ref,0) / sigma_ref(x,y)  [local source S/N, ref=250 or min over bands],
validity gate S_lam > k*sigma_lam(x,y) (k~3) using LOCAL noise. sigma_lam(x,y) = the
sd9-type local fluctuation map the flattening step already makes.
TEST (2 sources, src2 region ~6x noisier), T_smooth@center & M/Mtrue:
  brightness-only: 12.5/0.84, 15.5/0.88 ; global-sigma: same (no gain) ;
  LOCAL-sigma:     12.0/1.03, 14.8/1.05  <- correct T & mass in BOTH regions.
=> use LOCAL noise (global constant behaves like ignoring noise). Brightness/S-N weight,
NOT inverse fit-variance (confident cold outskirt fits have small variance -> would
dominate; tested earlier gave T~9.5K, M=2.9x).

## CLEAN ALGORITHM SPEC written: hires_T_algorithm.md  (full method, 6 steps + design table)

## Kernel strategy test: isolated all-fine; blends need per-source (getsf individual images)
Mixed field (compact 16", extended 52", pair 24" at 42" sep), T-fits at common 36.3":
ISOLATED sources (compact & extended): beam(18"), max(60"), per-source ALL recover
mass 1.00 -> kernel choice irrelevant when isolated.
BLENDED pair (C1=11K, C2=16K, 42" apart):
  GLOBAL max kernel: C1 T=14.3(true11) m=0.76 ; C2 T=15.6(true16) m=1.21 -- BLENDS Ts.
  PER-SOURCE (own getsf deblended image, own footprint-sized kernel, own column):
    C1 T=11.0 m=1.00 ; C2 T=16.0 m=1.00 -- exact.
=> Recommended path: run getsf on surfdens to completion (cheap: separation already
done) -> use DEBLENDED individual source images for per-source T AND per-source column.
Kernel size per source = footprint equiv-diameter D~2sqrt(Nvalid*pix^2/pi), floored at
beam; size REACH to footprint not bright core.
*** KEY RULE: combine COLUMNS, not T images. Each source's T converts its own column;
sum the columns. Summing T images re-blends neighbors (my earlier framing was wrong). ***
Bonus: getsf individual images remove the need to cut stamps by hand & solve blending.
Spec updated: hires_T_algorithm.md (per-source variant section).

## Clarification: source = single T_M at all resolutions; background = per-resolution series
hires' per-resolution T series (160-250@18.2", 160-350@24.9", 160-500@36.3") tracks the
scale-dependent effective T of the BLENDED image. After separation:
 - BACKGROUND: keep the full per-resolution series (cloud has real scale/position T
   structure + a stable floor -> per-pixel fits stable). Standard hires on bg bands.
 - SOURCE: isothermal once decloud­ed (bg-sub color T ~ T_M from ANY band interval,
   verified) -> ONE T_smooth ~ T_M used at ALL resolution levels of the source recon.
KEY: a spatially-constant T_M has NO resolution, so using it to convert the sharp 160
increment carries NO resolution-mismatch penalty -- the "sharp band / coarse T
unreliable" issue hires flags does NOT arise for the source. Resolution improvement
comes entirely from the intensity increments; T_smooth just converts them.
Fit source T_fit from bg-sub source bands at common (coarsest) resolution (any interval
agrees; full set = best S/N). Mild real gradients -> small residual -> RT grid.
Spec updated (Steps 2-4).

## REAL-DATA TEST on 18 injected cores (cfg.02 smoothed T = injected raw getsf extraction)
Known truth per core: M_BE=0.080, T_avg=11.98 K, M_SED~0.0716.
Findings (footprint surfdens-weighted, foots_injected_hiresSD):
 - Method runs end-to-end; smoothing FAITHFUL (raw fit == smoothed at cores, +0.0 K).
 - Recovered core T: median 13.8 K vs truth 11.98 -> +1.8 K WARM bias, roughly uniform
   (range 11-15.3 K; cores 168/147 recover ~12 K, cores 154/201/206 ~15 K).
 - Implied mass bias from T error: ~0.53 (250um) / ~0.63 (350um) of M_BE.
 - corr(T_rec, cloud column under core) = -0.16 (weak; cloud col only 2.5-4.4e21, little range).
 - Sanity: injected MODEL column integrated in footprints -> 0.98 M_BE (measurement OK).
DIAGNOSIS: warm bias is in the per-pixel fit (smoothing faithful, cores isolated) =>
the real getsf per-band SEPARATION leaves warm cloud residual under cores. Clean
(true-field) separation gave bg-sub color T=11.97 -> 0.92 M_BE; real getsf -> 13.8 K ->
~0.55. => T-smoothing algorithm is fine; SEPARATION QUALITY is now the dominant error,
LARGER than the ~10% RT-grid residual. Next lever: better background under sources
(lower-envelope/getimages-style bg; per-source deblended images for complex groups),
not the RT grid yet.

## CATALOG VERIFICATION: 2.4e22 block diagnosed (NOT a bug, NOT recompute-fixable)
Files: bes_model_params_catalog (recomputed 2.4e22) vs ..._cSD_04-bad (old).
 - DIFF: only 2.4e22 block differs; all other blocks byte-identical. Recompute changed
   the 45 models by ~1% (photon noise) -> reproduced same values, did NOT fix. Deterministic.
 - CAUSE: 160um zero-crossing in the 4-BAND SED fit only. IC160bs block medians cross zero:
   +23.8,+9.8,+2.76,[+0.40 @2.4e22],-0.94,-1.64. |T4bs-T3bs| spikes 0.81 K at 2.4e22 (~0 elsewhere).
   42 broken models have IC160 near zero (-2.9..+0.8): small residual 160 incompatible with cold
   SED -> fit collapses T->5K, mass up to 9x. 3 survivors have strongly-neg IC160 (-3.5..-7.8) ->
   dropped as non-detection (same as 4.8e22/9.6e22, which are clean).
 - 3-BAND fit (250/350/500) clean in ALL blocks: Chi2~0, 0 flagged; M_SED3bs/M_BE smooth monotonic
   0.85,0.87,0.90,0.91,0.96,1.03. The ~15 4-band-flagged models in other blocks all clean in 3-band.
 - Analytic BE cols (M_BE,R_BE,rho_BE,T_BE) + Td_avg/Td_emb: smooth, finite, positive, monotonic, all blocks.
VERDICT: do NOT recompute other blocks (all clean). Fix is band selection, not recompute:
 S/N-gate 160um (treat insignificant 160 as non-detection) OR use 3-band at/above 2.4e22.
PAPER: this 160 collapse would hit any 4-band fit on real data at Sigma~2.4e22 with near-zero 160 ->
 replicate the target catalog's 160 S/N gating in the grid; 3-band mass is the clean backbone.

## Konyves-based per-band noise -> fixes 2.4e22 160 collapse + defines completeness
Noise model: sigma_lambda = sqrt(sigma_instr^2 + (f*I_lambda,back)^2), f=sigma_N/N_back,
 sigma_N=3.9e20*(N/7e21)^1.6 (Konyves+2015). I_lambda,back from interior-plateau T_back
 (exclude warm rim) via O&H beta=2 (kappa=83.364(100/lambda)^2 /g dust), mu=2.8, dtg=0.01, MJy/sr.
T_back per block (interior): 13.9,11.6,9.74,8.47,7.67,7.02 K for SD_emb 3e21..9.6e22.
f per block: 3.4,5.1,7.7,11.7,17.7,26.8 %.
RESULT (S/N>1 band gate on catalog IC fluxes):
 - 160 KEPT <=1.2e22 (S/N 10->2.1), DROPPED >=2.4e22 (S/N 0.4). At 2.4e22 all 45 drop 160,
   incl all 42 broken -> fall to clean 3-band M_SED3bs/M_BE=0.91. Fix is robust (S/N 0.4, far<1,
   normalization-insensitive). THE 1e-10 placeholder noise had disabled this gate; real noise restores it.
 - High column side-effect: at 4.8e22/9.6e22 cirrus f (18-27%) also sinks 250-500 below S/N 1 for
   LOW-mass/low-peak cores -> median bands kept 1, then 0. CONFIRMED mass/peak-dependent (4.8e22 survivors
   M_BE 0.91-41, dropped 0.01-0.6) = real Konyves completeness limit, NOT uniform artifact.
CAVEATS for completeness BOUNDARY precision (NOT for 160 fix): (1) verify I_back abs normalization vs
 fitfluxes opacity zero-point; (2) band-inclusion sigma for compact cores should be cirrus power at SOURCE
 scale (getsf-filtered), smaller than raw per-pixel rms used -> top-column losses are an upper bound;
 (3) add instrumental floor (matters faint/low-column only).
DECISION pending: let gate define usable grid region vs keep noiseless fits + separate completeness layer.
Deploy 160 gate either way.

## Konyves+2015 Appendix B.2 completeness model -- GRID REPRODUCES IT (validation)
B.2 apparent column-density significance:
 S~ = 1.5 * [1/(1+(HPBW/R_BE)^2)] * [B350(T_core)/B350(T_back)] * 18*(N_back/7e21)^-0.6
 (intrinsic BE contrast 1.5; beam dilution; temperature dilution at 350um fiducial;
  cirrus factor = 1/f, same B.1 law). Detected/~90% complete when S~>5; C=F[S~] (their Fig B.5).
Aligned to catalog: R_BE=R_BE_as, T_core=Td_avg, T_back=interior(13.9..7.02), HPBW=18.2", N_back=SD_emb.
RESULT (grid): completeness mass 0.02 Msun @3e21 -> 6.4 @4.8e22 -> 18 @9.6e22. Survivor counts
 (S~>5): 14/31 @4.8e22, 2/24 @9.6e22 == our independent per-band cirrus gate (14/31, 3/24). MATCH.
 => grid's high-column completeness boundary reproduces published HGBS B.2 for free. Two thresholds,
 one cirrus law: per-band S/N>1 = SED band inclusion (160 fix); S~>5 = source detection/completeness.
 Temperature-dilution term B350(Tcore)/B350(Tback) = same LOS-temperature physics as our mass bias;
 B.2a/b report derived/true mass ~0.8 and SED T high by ~0.8-1 K == our 3-band M_SED/M_BE 0.83-0.91.
 Grid reproduces Konyves bias AND completeness -> validated as correction-surface baseline.
