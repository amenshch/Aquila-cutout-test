# Session summary — grid detectability redesign, peak calibration, hole-clip

## Context
Methods paper `latex/rt_mass_correction.tex`. Grid axes are now
**(Σ_emb, M, FWHM)** in physical units, with T_BE, ξ_max, ρ_c **derived** per node.
Repo: github.com/amenshch/Aquila-cutout-test. Distance 260 pc, μ_H2=2.8,
μ_particle=2.33, dust:gas=0.01, opascale=1.0746, hires beam 13.5", ref beam 18.2",
3" common map pixel.

## What was settled this session

### 1. Grid-construction knobs (Fortran side, all resolved)
- **R_inner = r0/100** (i.e. `Rinner = radscale*Rmodel/(100*ximax)`), not R_BE/100 —
  ξ_max cancels, giving a fixed ξ_in=0.01 for all models. Same for the pixel scale below.
- **Cloud zones fixed via config key `nzcloud`** (~12–15), not `fncore`; the cloud log
  span (~1.36 dex = log10(30/ftrans)) is model-independent, so its cell count should be too.
  `fncore` retired (derived, not input). `nzones=50` kept fixed by user choice.
- **Core geometric grid `fgeom=10` kept uniform** — automatically finer cells for small-ξ
  (steep T-gradient, big density drop at R_BE) and coarser for large-ξ. User confirmed the
  finest cells sit near R_BE where the mass and ⟨T_dust⟩ live, so deep models are OK.
- **tanh cloud grid**: with fixed δ the first cloud cell matches the core's last cell at
  one ξ only; user accepted the mild (≤~4×) mismatch since the cloud is radiatively inert.
- **ftrans=1.3 kept** (scale-free fraction of R_BE).

### 2. Pixel-size rule for RADMC imaging (mc3d script, config-param based)
Final rule (user sets pixel outside Fortran, in the mc3d script):
```
pixel = r0 * min(ximax, phi(ximax)) / 7.5     # 7.5 px per min(R_BE, FWHM)
phi(xi) = 1/sqrt(1/(3*xi^2) + 1/6.1^2)        # closed-form projected-FWHM/r0
```
UNITS: r0 comes out in cm from `r0 = c_s/sqrt(4 pi G besden)`, with
`c_s = sqrt(GasConst*bestem/2.33)`. Convert `r0_AU = r0/xAU`, `r0_as = r0_AU/distpc`.
Then `mpixel = r0_as*min(...)/7.5`. Earlier bug: r0 left in cm while comparing to
pixres=3" — must convert to arcsec first. Use `$ximax*$ximax` not `^` in CALCULATE.
The min() rescues compact high-ξ cores from collapsing to ~1 px FWHM.

### 3. run_grid_mc3d.sh — rewritten
Reads params per model from the grid table and writes ONLY four config keys via
`set_param` (awk, preserves trailing `! comments`): **besden(ρ_c), bestem(T_BE),
ximax, embden(Σ_emb)** + nphot/setthreads. Does NOT touch `Rbesph` (a FLAG=1.0, not a
radius), `Router`=0.0 (auto → 30·R_BE), `Rinner`=0.0 (auto). Copies +template/ (incl.
+plots/, carried through by cp -a). Ragged-tree tolerant; per-model lock; skip-completed.
Tree: `cSD_{i:02d}/M_{j:02d}/{k:02d}`.
TODO (user): add one line to **skip detect==0 nodes**.

### 4. add_recoverable_mass.py — updated
- Tree path fixed `tBE_{j}` → `M_{j}` (j is now mass index). THIS was the one real bug.
- `PIX_ARCSEC=3.0` is CORRECT (all stamps resampled to common 3" grid before getsf;
  the variable RT pixel matters only upstream). My per-header-pixel change was reverted.
- **Detection threshold replaced**: the old `5*N_rms` Konyves-cirrus cut was a category error
  (getsf's per-single-scale combination σ, not the fluctuation level of the original image).
  Replaced with an empirical floor calibrated on real Aquila `ok` sources:
  `N_src_min(N_back) = 1.971e21*(N_back/1e22)^1.342`. (Provisional; injection-recovery will
  supersede it.) `recoverable_fraction` physics otherwise unchanged (geometry-agnostic).

### 5. THE KEY PHYSICS RESULT — detectability & the grid-design flaw
Using the user's contrast definition **contrast = PEAK^SBF/PEAK^BGF = CONTRSRC03
= 1 + PEAK^SRC/PEAK^BGF** (dimensionless):
- **Real Aquila detection floor: contrast ≈ 1.16** (min 1.158, 1st pct 1.163) over 136
  reliable `ok`-SED sources. No real core detected below ~16% peak excess over background.
- The old grid built **89 models below this floor**, ALL of which correctly got frac_rec=0.
  So the frac→0 collapse was NEVER a threshold bug — it was the grid building undetectable,
  low-contrast models.
- **Root cause**: the grid samples the SAME FWHM range at every Σ_emb, but the detectable
  FWHM ceiling collapses with column (194"→18") and the min detectable ξ_max rises
  (0.5→7.4). Contrast is driven by **compactness (ξ_max)**, not mass (corr contrast–ξ=0.73);
  peak excess does NOT grow with M (beam dilution spreads massive cores).
- **Physical picture (user, correct)**: cores in the densest environments (filaments) are
  smaller and more concentrated (high ξ_max) — infalling, gravitationally unstable
  configurations. Detectability selection and astrophysics point the same way.

### 6. Peak-excess prediction — now EXACT, CAL=1.0 (long debugging saga, resolved)
The convolved, background-subtracted central column ICSDbs is analytically predictable
because "bs" in the model grid is just a **clean constant-Σ_emb pedestal subtraction**
(NOT getsf's nonlinear interpolated background). The correct recipe:
- Project the truncated BE profile → Σ(θ) (excess above pedestal).
- **2D Hankel convolution** with the 13.5" Gaussian (NOT a 1D radial average — that was
  the first bug, gave 3.6× too high):
  `C(ρ) = ∫ Σ(r) exp(-(r²+ρ²)/2σ²) I0(rρ/σ²) r dr / σ²`, integrated to ≥5σ (beam), not
  just the source radius.
- **Average C over the central 3" pixel** (the resampling step).
- Pass **r0** (BE scale radius) as the angular scale, NOT R_out=ξ_max·r0 (the second bug,
  an argument swap that gave the 3.5× factor).
VALIDATION (user's directly-measured model-1 values):
  intrinsic point peak 5.585e22 (mine) vs 5.56025e22 (user) — 0.4%;
  convolved+resampled 1.0499e21 (mine) vs 1.04961e21 (user) — **0.02%**.
So **CAL ≡ 1.0, no fudge factor.** The user's 2.5085e22 was the *3"-pixel-averaged*
intrinsic peak (lower than the true b=0 point value because a 3" pixel smears a 2.3" core).

### 7. Grid framing — UNIVERSAL, region-independent (final design)
- The grid stays **277 physical BE nodes on the full regular 6×12×10 (i,j,k) lattice**.
- Detectability is a **per-node flag** (`detect`), a function of **background column only**:
  contrast = 1 + peak_exc/Σ_emb ≥ FLOOR. Nothing region-specific enters the grid.
- **FLOOR = 1.1** chosen (slightly below Aquila's 1.16) to include denser regions that may
  detect marginally deeper. A different cloud just changes this one number.
- i/j/k structure preserved; undetectable nodes stay in the catalogue with detect=0 and are
  **skipped at RT time** (empty leaves). User's insight: empty j,k cost nothing.
- Current `bes_grid_final.txt`: 277 nodes, columns incl. peak_exc, contrast, detect, CAL=1.0.
  ~219–236 detectable depending on exact floor. ID 1-based.
- Corrector `fit_mass.py` uses `LinearNDInterpolator` (scattered/Delaunay), so the ragged
  detectable subset triangulates directly — no corrector change needed. Feeding only
  detect=1 nodes improves conditioning.

### 8. OPEN ISSUE — the ρ_edge < ρ_emb "hole" and massive high-column cores (BLOCKING)
Checking whether detectable grid nodes span where real cores live revealed:
- At Σ_emb=4.8e22, real cores have M=0.9–7.6 M☉, but grid detectable nodes top out at
  M≈0.58 M☉. The grid MISSES the massive detectable cores that dense filaments actually
  produce.
- Reason: those massive cores (ξ=7–36) **fail the ρ_edge ≥ ρ_emb hole-clip** — their BE
  edge density falls below the cloud density, so `node()` currently returns None. They are
  FORBIDDEN by construction, not merely unsampled. **This is why "nothing changed" in the
  regenerated grid** — the missing models are excluded by the hole-clip.
- So the grid genuinely under-covers, but the fix is a **model-scope decision**, not sampling.

**User's proposed resolution** (to test first): keep the hole, but build **several zones at
ρ = ρ_edge** just outside R_BE before the density rises to the cloud ρ_emb. Rationale: a real
interpolated background subtraction samples the LOCAL floor (ρ_edge plateau) across the source
footprint, not the distant cloud — so after bg subtraction the convolution+resampling never
see the hole. This would admit the massive high-column cores with clean observables.

**NEXT STEP (user is doing):** test the CURRENT RADMC-3D setup on hole models to see whether
bg subtraction already handles the hole correctly (maybe no plateau needed). Representative
test models provided (all Σ_emb=4.8e22):
| M | ξ_max | edge/emb | besden(ρ_c) | bestem(T_BE) | ximax |
|---|---|---|---|---|---|
| 1.0 | 3.81 | 1.13 (control, no hole) | 7.1957e-20 | 7.6231 | 3.8079 |
| 2.0 | 7.04 | 0.63 | 1.0227e-19 | 7.6005 | 7.0369 |
| 2.0 | 16.52 | 0.39 | 4.8901e-19 | 7.5010 | 16.5238 |
| 4.0 | 16.21 | 0.20 | 1.2256e-19 | 7.5869 | 16.2055 |
| 7.6 | 36.40 | 0.08 | 1.2435e-19 | 7.5843 | 36.3955 |
embden=4.8e22 for all. If bg subtraction works even at edge/emb=0.08, drop the hole-clip
entirely; if it fails progressively, add the ρ_edge plateau (and the edge/emb threshold where
it breaks becomes a construction criterion).

## Decisions pending user input on resume
1. Result of the hole bg-subtraction test → keep hole-clip / add plateau / drop clip.
2. If plateau added: (a) how many zones / how far out the plateau extends; (b) whether Σ_emb
   for the T_BE calibration (N_tot shielding) vs the contrast denominator uses cloud ρ_emb or
   plateau ρ_edge. (Claude's view: N_tot uses cloud column for ISRF shielding; contrast
   denominator uses the local plateau that's actually subtracted.)
3. Then regenerate grid with newly-admitted massive nodes, re-run RT (10^9 first, ~1 day, to
   validate correctness; then 10^10 for deep high-ξ), re-run add_recoverable_mass.py, refit
   the 4-D corrector.

## Draft sections needing update (old grid described)
- sec:grid — still says geometric (Σ_cloud, T_BE, ρ_c) axes; must become (Σ_emb, M, FWHM)
  physical grid with derived ξ_max, detectability flag, universal framing.
- New subsection needed: **detectability & contrast floor** (the contrast=SBF/BGF definition,
  floor 1.1, peak-excess Hankel prediction, why it's region-independent).
- sec:uniformcloud / sec:edge — connect to the hole-clip domain-of-validity discussion.

---

## UPDATE (continuation) — hole-clip dropped, grid → 377, BE-family boundary

### Hole-clip resolved and REMOVED
- User changed the reduction to subtract the constant Σ_emb pedestal **at R_BE
  on the original (unconvolved) image, before convolution** (not the old
  convolve-then-subtract order that created the spurious depression). This
  removes the depression artefact entirely.
- User ran the 5 test hole models (Σ_emb=4.8e22, edge/emb 1.13→0.08). Direct
  surfdens integration:
  | Model | ξ_max | edge/emb | M_true | Crater | Interp (positive px only) |
  |---|---|---|---|---|---|
  | 1 | 3.81 | 1.13 | 1.0 | 0.999 | 0.555 |
  | 2 | 7.04 | 0.63 | 2.0 | 2.000 | 0.914 |
  | 3 | 16.52 | 0.39 | 2.0 | 1.998 | 0.941 |
  | 4 | 16.21 | 0.20 | 4.0 | 4.000 | 1.097 |
  | 5 | 36.40 | 0.08 | 7.6 | 7.598 | 1.147 |
  Crater recovers truth to <0.1% at ALL hole depths → the hole causes no bias.
  Interpolated under-estimates increasingly with concentration/hole depth
  (56%→15%). Interp mass integrated over POSITIVE pixels only (observer would
  not trust negative surfdens). Model 4 < Model 3 because its hole is deeper.
- **RT handles the hole fine** (user confirmed): photons cross the depression
  with few interactions, T neither spikes nor drops. The hole schematically
  resembles the rarefaction of an inside-out collapsing core → plausibly real,
  not an artefact.
- **DECISION: the continuity constraint ρ_edge ≥ ρ_cloud is DROPPED.** Crater
  is the model reference (defines truth); real extractions always interpolate.

### Grid regenerated: 277 → 377 models
- Dropping the constraint admits **100 new "hole" models** (0 at Σ=3e21 rising
  to 39 at 9.6e22), reaching M up to 18 M☉ at high column — filling the
  massive-compact corner the constraint had removed.
- `bes_grid_final.txt` now 377 rows, 23 columns, with **both `detect` and
  `hole` flags** (1/0). ~260 detectable. ID 1-based.
- **Column layout (23):** ID i j k SD_emb T_BE M_BE FWHM_pc xi_max contr_rho
  peak_exc contrast detect hole rho_c rho_edge rho_emb r0_AU R_out_AU
  R_cloud_AU N_tot a_BE stab.
- **run_grid_mc3d.sh FIXED** for the new layout: the `read` line broke because
  peak_exc/contrast/detect/hole were inserted before rho_c (shifted its
  column). Updated read list; added a `detect==0` skip (ALL_NODES=1 overrides).
  The 277 old models' params are UNCHANGED (existing RT stands); only the ~new
  detectable hole models need RT.
- **CAUTION:** add_recoverable_mass.py reads grid columns by position too —
  check its --col-* indices against the new 23-column layout before running.

### The BE-family upper boundary (KEY SCOPE RESULT)
- Checked whether the grid covers observed compact-massive sources: **25/685
  (~4%) lie OUTSIDE the grid hull**, in the compact (FWHM≲30″), massive
  (M≳1 M☉) corner (up to M=19.7 M☉ at 47″).
- These need **ξ_max ~ 70–150** to be BE spheres — and colder T makes it WORSE
  (lower c_s needs even higher concentration). That is central-to-edge contrast
  far outside any quasi-static isothermal / physically sensible single-core
  regime. **They are not BE spheres.**
- Interpretation (user agrees): blended multiple sources at 260 pc, and/or
  interpolated bg sweeping in filament mass → apparent mass inflated with no
  size increase → displaced up-and-left out of the BE locus.
- **DECISION: state the BE domain boundary, do NOT extend the grid to cover
  them.** The grid's most-massive boundary IS the physically plausible upper
  envelope for genuine single, cleanly-measured BE cores. A source above it is
  a diagnostic flag (blended/contaminated), not a core to correct. Method
  applies to ~96% of cores consistent with the BE family.
- Options considered and rejected: (2) add a power-law-envelope model — a
  separate core model / separate paper; (3) push ξ_max→150 — dressing an
  unphysical fit as coverage.

### Draft updated (all compiles clean, 23 pp)
- §Grid construction: continuity constraint retired (subtract-first removes
  depression; hole ≈ collapsing-core rarefaction; RT handles it).
- §Crater vs interpolated: replaced old single-model table with the 5-model
  hole sequence (crater <0.1%, interp 56%→15%).
- §Mass floor and high-mass cap: rewritten — high-mass boundary is the BE
  family limit (not a CMF margin), doubles as a blend/contamination diagnostic,
  25 outliers need ξ~70–150 so are not BE cores, 96% coverage stated as domain
  of validity.
- §uniform-cloud, §low-mass edge, §Summary, §embedding-cloud: all references to
  the dropped constraint reconciled; low-mass edge now attributed to
  detectability.
- Working-ranges \note: 277 → 377 (277 + 100 formerly hole-clipped).

### Still pending on resume
1. Higher-photon (10^10) T_BE recalibration → finalize node count / \note.
2. RT-run the new detectable hole models (run_grid_mc3d.sh picks them up,
   skips the 277 done).
3. Re-run add_recoverable_mass.py (CHECK column indices for 23-col layout) →
   updated _rec catalogue.
4. Refit 4-D corrector on the fuller (377) grid.
5. Figure: M–FWHM for 377 grid vs Aquila showing new nodes + the 25 outliers
   beyond the BE boundary (draft file grid_M_FWHM_new377.pdf exists).
