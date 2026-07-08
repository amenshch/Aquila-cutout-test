# Temperature-corrected high-resolution surface density (hires-T)

A self-contained method that removes the cloud-contamination temperature bias from
hires surface densities, before any getsf detection, while keeping hires resolution.

---

## Why

hires forms each per-resolution column as `D = I / [B(T) kappa eta mu m_H]`, with `T`
fit per pixel. On a source the per-pixel color temperature is the cloud + core blend
(warm, ~15 K), so `B(T)` is too large and the source column is suppressed (peak lost,
mass ~0.67 of true). The cloud-removed core temperature is the mass-averaged `T_M`
(~12 K), which by construction recovers the correct mass. The method derives the source
column with `T_M` instead of the blended temperature.

---

## Inputs

- Observed band images `I_lam` at native resolutions `O_lam` (e.g. 160/250/350/500 um
  at 13.5/18.2/24.9/36.3").
- Per-band LOCAL fluctuation maps `sigma_lam(x,y)`. In practice this IS the per-band
  getsf flattening image `Fflat_lam(x,y)` -- the small-scale local standard deviation
  (`sd9`-type), median-filtered to exclude the sources themselves, so it represents the
  background + noise fluctuation level at every pixel with sources removed. Use each
  band's own `Fflat_lam`. No conversion to physical noise units is needed: only the
  RELATIVE variation of `sigma_lam(x,y)` across pixels matters (the normalized
  convolution and the gate use ratios), so any overall normalization cancels.
  Using `Fflat_lam` also ties the temperature weighting to the SAME significance
  reference getsf uses for detection (its detection image is `S_lam / Fflat_lam`), so
  there is no second, independent noise model to reconcile.

---

## Step 1 -- per-band source / background separation

Run the per-band getsf decomposition (independent per band, no detection needed).
For each band obtain a source image `S_lam` (cloud removed) and a background image
`Bg_lam` (cloud only), with `I_lam = S_lam + Bg_lam`.

## Step 2 -- source temperature image `T_smooth`  (the key new piece)

NOTE on resolution. Standard hires carries a SERIES of temperatures at different
resolutions (160-250 @ 18.2", 160-350 @ 24.9", 160-500 @ 36.3") because the BLENDED
image has a scale-dependent effective temperature. After separation this is needed only
for the BACKGROUND (Step 4). The SOURCE is essentially isothermal once the cloud is
removed -- its background-subtracted color temperature is ~T_M from ANY band interval
(verified: 160/250, 250/350, 350/500, 250/500 all give ~T_M) -- so the source needs a
SINGLE `T_smooth` ~ T_M, used at ALL resolution levels of the source reconstruction.
Because that temperature is spatially ~constant it has no resolution, so using it to
convert the sharp 160 increment carries NO resolution-mismatch penalty (the
"unreliable" sharp-band/coarse-T issue hires flags does not arise for the source). The
resolution improvement still comes entirely from the intensity increments.

a. Per pixel, fit a single-temperature modified blackbody to the source bands `S_lam`
   at the common (coarsest) resolution -> `T_fit(x,y)`. (Any band interval agrees since
   the source is isothermal; use the full set for best S/N.)

b. Validity gate: mark a pixel INVALID if any `S_lam <= 0` or any `S_lam < k*Fflat_lam`
   (k ~ 3, using the per-band flattening image as the LOCAL fluctuation level). These are
   the over-subtracted / faint-outskirt pixels whose rising SED a single-T fit mis-reads
   as spuriously cold.

c. Weight (per valid pixel): `w(x,y) = max(S_ref, 0) / Fflat_ref(x,y)`  -- the local
   source signal-to-fluctuation ratio in a reference band (250 um), or
   `min_lam (S_lam/Fflat_lam)` for a stricter version. Invalid pixels: `w = 0`.
   *** Use brightness / S-N, NOT inverse fit-variance: a confidently-but-wrongly cold
   outskirt fit has a small formal error and would otherwise dominate. ***

d. Smooth by NORMALIZED CONVOLUTION with a kernel image `K`:
       `T_smooth = conv(w * T_fit, K) / conv(w, K)`     (where `conv(w,K)` > small floor)
   - two ordinary convolutions + a divide, using existing convolution code;
   - `K` = plain Gaussian, scale ~ a few beams (or the largest expected source);
   - kernel SHAPE is irrelevant (the denominator divides out its normalization);
   - a single fixed kernel suffices; multiscale is optional (a few passes at growing
     scale, per-pixel pick the finest scale with `conv(w,K)` above threshold).

`T_smooth` is then defined and smooth in every pixel: it equals the locally
mass-weighted (= brightness-weighted) temperature ~ `T_M` on sources, and is a harmless
smooth value in empty regions (where `S_lam ~ 0` carries no mass). It runs GLOBALLY, so
no source positions / prior detection are required.

## Step 3 -- source surface density

Run hires on the SOURCE bands `S_lam`, using the single `T_smooth` at ALL resolution
levels: `D(O_lam) = S_lam / [B_lam(T_smooth) kappa_lam eta mu m_H]`, combined across
resolutions as in standard hires. The multi-resolution improvement comes from the
intensity increments; the temperature is the one flat `T_smooth`. -> `Sigma_source`
(sharp, correct peak and mass).

## Step 4 -- background surface density

Run STANDARD hires on the BACKGROUND bands `Bg_lam` -- including its normal per-
resolution temperature series (160-250, 160-350, 160-500) -- since the cloud has real
scale/position temperature structure and retains its emission floor, so those per-pixel
fits are stable. -> `Sigma_cloud`.

## Step 5 -- combine

`Sigma_total = Sigma_source + Sigma_cloud`.  Columns add along the line of sight, and
each component now carries its own correct temperature -- which the single blended-T map
cannot represent. This is exactly why separation works (hires is linear in `I` at fixed
`T`; the only nonlinearity is the temperature fit).

## Step 6 -- extraction

Run the single full getsf extraction (detection + measurement) on `Sigma_total` as the
surface-density image, as for any other band.

---

## Variant: per-source temperatures from getsf individual images (recommended)

The global formulation (Step 2) needs no detection, but a single kernel BLENDS the
temperatures of close sources with different T (tested: a 11 K and a 16 K core 42" apart
-> global kernel assigns 14.3 K to the cold one, mass 0.76; the hot one 15.6 K, mass
1.21). Since getsf on the surfdens completes cheaply (the costly part, per-band
component separation, is already done), run it and use its DEBLENDED individual source
images. For each source:
  - fit T per pixel and weight on ITS OWN image (cloud + neighbors already removed);
  - size its Gaussian kernel from its OWN footprint extent -- equivalent diameter from
    the valid-pixel area, `D ~ 2*sqrt(N_valid*pix^2/pi)`, floored at the beam; size the
    REACH to the footprint, not just the bright core;
  - build `T_smooth` for that source and convert ITS OWN column `Sigma_k`.
Then ADD THE COLUMNS: `Sigma_source = sum_k Sigma_k`.
*** Combine columns, NOT temperature images: each source's T converts its own column;
summing T images re-blends neighbors. ***
Tested on the same 11/16 K pair: per-source -> T = 11.0 / 16.0 K, mass 1.00 / 1.00.
Isolated sources (compact or extended) recover 1.00 with ANY kernel, so this matters
only in crowded regions -- but it removes blending entirely and needs no stamp-cutting
(getsf supplies the individual images). Cost: it runs AFTER a getsf detection pass
rather than purely inside hires; given the separation is already done, that cost is
negligible.

---

## Two-tier accuracy

- Steps 1-5 remove the DOMINANT cloud-contamination temperature bias
  (recovers ~0.67 -> ~0.9 of true mass at the map level).
- The residual ~10% (single-T-fit vs mass-conserving gap; environment-dependent, not
  cleanly gradient-correlated) is corrected by the RT grid (`ObservableCorrector`).

---

## Key design choices (and why)

| choice | reason |
|---|---|
| separate per band, then hires per component, then add | each component gets its own correct T; the blended single-T map suppresses the source peak |
| source uses smoothed `T_smooth`; background keeps per-pixel T | sources-only bands lack a floor -> per-pixel T blows up at outskirts; background retains the cloud floor |
| weight = local S/N (not fit variance) | spurious cold outskirt fits are tight (small variance) -> would dominate inverse-variance weighting |
| local `Fflat_lam(x,y)` (getsf flattening image), not global | fluctuation level varies by orders of magnitude across the field; tested -- local recovers correct T & mass, global does not. Source-excluded & smoothed; same significance reference getsf uses for detection |
| validity gate on `S_lam > k*Fflat_lam` | removes truncated / sub-fluctuation pixels entirely rather than down-weighting them |
| normalized convolution, global, fixed Gaussian kernel | needs no source positions (brightness localizes); kernel shape/scale not critical due to normalization |

---

## Caveats

- Strongly overlapping sources of very different temperature are the hard case for a
  global kernel (it can blend them). Fallback: one getsf run on the current surfdens for
  positions/footprints, then per-source temperatures.
- The recovered absolute mass also depends on the per-band separation quality
  (background over/under-subtraction) -- a separate, smaller error than the temperature
  bias, and within the RT grid's reach.
