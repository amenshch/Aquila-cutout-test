# Addendum to session summary — recovery-table method and the C_M≈2 result

Append to `Summary_from_previous_chat_4.md`. This covers the discussion **after**
that summary was written and supersedes it wherever they conflict. Terminology
is stated explicitly throughout because it caused repeated confusion:

- **f_rec** = recovered flux fraction = (flux measured inside the
  detection-limited footprint, after the flat interpolated background is
  subtracted) / (total flux the model emits in that band). Dimensionless, 0–1.
- **C_M** = mass correction factor = M_BE (true Bonnor–Ebert mass) /
  (mass the extraction reports). Dimensionless, normally >1.
- **R_foot** = footprint-to-size ratio = 2r/FWHM in models, FOOA/AFWHM in getsf.
- **C_M(0.1–2)** = median C_M over sources of *measured* mass 0.1–2 M⊙ (a
  matched-range comparison so clouds sampling different masses compare fairly).

---

## HEADLINE RESULT — C_M ≈ 2 for Aquila (this replaces the old open question)

The correction factor for Aquila is **C_M ≈ 1.9–2.2 with f_rec ≈ 0.48**, of
which ≈1.16 is temperature bias and the rest (≈1/0.48≈2.1) is flux loss outside
the footprint. Three independent matching-axis choices agree: (M_rec, Σ,
conc_peakmean)→2.19, (M_rec, Σ, R_foot)→2.03, all four→1.87. This is higher than
the old single-value-frac 1.65 but now rests on a **measured** constraint
(R_foot), not an assumed noise relation. **In the draft** as `sec:footprint` +
`Table~\ref{tab:sevencloud}`.

### How this was reached — the degeneracy and its resolution
- **The degeneracy (restated cleanly):** f_rec is not observable — the total
  flux is exactly the unknown the correction supplies. So a measured source is
  reproduced equally by a well-recovered model (f_rec→1) and a poorly-recovered
  truncated one (f_rec small), and these give C_M differing by up to 10×. During
  the session, letting the parametrisation choose freely gave C_M = 1.15, 1.65,
  4.29, 4.47, 6.43, 10.43, 13.12 — the swing IS the degeneracy, not progress.
- **Alexander's resolution (the key idea):** R_foot cannot approach the beam —
  background fluctuations broaden any measured profile, so getsf never returns a
  footprint much below ~1.5×FWHM. Truncating a model profile at small r
  artificially steepens it and gives an unphysically small R_foot. Imposing the
  observed R_foot removes exactly the low-f_rec rows that produced the large
  factors, and the answer converges.
- **Verified:** measured R_foot = FOOA/AFWHM is **constant at 1.81–2.02 across
  all seven clouds** (Scorpius 2.02, Ophiuchus 1.81, Aquila 1.82, Orion 1.82,
  California 1.93, Cygnus 1.98, W3/W4/W5 1.91), 16–84% ≈ 1.5–2.5. So it is set
  by the extraction, not the sources. This ≈ the η=2 that
  add_recoverable_mass.py already assumed — that assumption was right.
- **Model tables at R_foot≈1.8–2.0 give f_rec≈0.45–0.65** (read directly from
  the v2 recovery tables), consistent with the matched result. A Gaussian at
  R_foot=1.82 keeps 90% of its flux; **BE profiles keep only ~45–50%** because
  of their extended low-surface-brightness envelope — that envelope is the
  physics the correction recovers, and why C_M≈2 not ≈1.1.
- **Adding `peak` as a 5th axis drops it back to 1.65** with inferred f_rec 0.16
  — inconsistent (that axis over-reaches into unphysical small-footprint rows).
  Do NOT use peak as an axis; conc_peakmean and R_foot are the good ones.

### Seven-cloud extension FAILS — distance problem (do NOT put raw C_M in draft)
Applying the 4D corrector (M_rec, Σ, conc, R_foot) to all seven:
- **C_M(0.1–2) is STABLE at 2.09–2.95** across 13× distance — consistent with
  Aquila. This is the trustworthy number.
- **Raw C_M runs backwards with distance** (2.45→0.81) and inferred f_rec falls
  (0.50→0.12). These are mutually contradictory: f_rec=0.12 implies flux-loss
  ~8×, not 0.81. CMF slopes steepen under correction, reversing prior results.
- **Cause:** reported mass is a matching axis but is DISTANCE-DEPENDENT, while
  the grid is fixed at 260 pc. A source at 1700 pc with identical structure has
  a mass 43× larger → matches a structurally wrong node. Only Aquila (grid
  distance) is self-consistent; matched-range repairs the others.
- **Fix needed before quoting multi-cloud raw factors:** a distance-invariant
  mass proxy, OR grids computed at each field's distance. Flagged `\note{TO DO}`
  in `sec:footprint`. The aquila-sim getsf injection (still ~a week out) will
  measure f_rec directly for known-flux sources and settle the normalisation.

---

## THE RECOVERY-TABLE METHOD (new machinery, replaces single-value frac)

**Why:** the old scheme collapsed two broad distributions to single values —
per-source noise spans 3.3–5.7× at fixed background (real variation, not
measurement error), and R_foot spans 1.5–2.5 with 62% of sources outside η=2±10%.
Both collapses inject error. The tables tabulate recovery quantities vs
**truncation radius** so the matching variable is chosen afterwards; everything
(frac, conc_peakmean, R_foot, peak) is monotonic in r.

### `python/tabulate_recovery.py` (NEW, staged to outputs)
Post-processes existing band stamps (no new RT). Per node, per band, per
truncation radius r, writes: r, r/R_conv, peak, I_rim, peak_rim (=peak/I_rim),
frac, fwhm_rec, cfoot (=2r/FWHM), conc_peakmean, conc_slope, rel, absorp.
Key design points, all learned this session:
- **R_conv from the source mask** (`nc.surfdens.bsl.{res}x0.rs3p0as.mask.fits`,
  1|0 pixels), equivalent-area radius √(N_pix·pix²/π). The mask defines the
  footprint the "bs" background was interpolated over, so frac=1 at R_conv is
  exact. Falls back to a profile threshold only if no mask.
- **mc3d builds masks by clipping surfdens <1e16** (Alexander confirmed). Fine —
  captures all flux — but reaches 1.25–1.47× beyond the 99%-flux radius (worst
  for compact/coarse-beam). So DON'T sample radii in r/R_conv; sample by target
  **frac** or target **cfoot** instead (implemented: `--sample {frac,cfoot,both}`).
- **fwhm_rec must be INTERPOLATED to half-max**, not the first bin below it.
  getsf gets AFWHM/BFWHM by interpolating intensity to ½-peak — SAME quantity;
  the bin-centre version quantised to 3″ and was incomparable. FIXED (v1 tables
  had ~6 distinct fwhm_rec values, v2 has 3833). **Claude was wrong twice here**
  — first claimed getsf "fits" FWHM (it interpolates, like the model), then had
  a quantisation bug; both corrected.
- **absorp flag** for negative central pixel (70/100 µm cold cores; also a
  depression-model diagnostic). Validated: absorption in 98% of nodes at 70µm,
  83% at 100, 25% at 160, 0% at 250–500 — clean monotonic wavelength trend.
- `--cfoot-min/max` filter, `--fwhm-floor` (beam units) secondary size floor.
- **7 bands** SD,070,100,160,250,350,500 (matches grid catalogue).

### Band-floor calibration — mostly SUPERSEDED by the R_foot approach
- `python/calibrate_floors.py` (NEW, staged): per-band detection floor from a
  getsf catalogue by the **envelope method** (low percentile of PEAK^SRC in bins
  of PEAK^BGF). Envelope scatter 1.03–1.39× vs 1.5–250× for k·FXP_ERR — the
  envelope is far better, FXP_ERR is NOT a usable noise proxy for this.
- BUT the envelope floor is percentile-dependent (factor 2.3 between 1st and
  20th pct) because it's a property of the detected population, not the map —
  it never converges. This is WHY the floor-based frac was unreliable and gave
  the factor-of-4 discrepancy vs the R_foot answer. The **R_foot constraint
  sidesteps the floor entirely** and is the preferred route now.
- Index p varies per cloud (band02: Aquila 1.44, Ophiuchus 1.05, Orion 1.18,
  Cygnus 1.10) → floors must be calibrated per field, never imported.

### Physical band-floor idea (Alexander's) — VALIDATED as cross-check
`floor_b = I_bg,b · (floor_SD/N_bg) / (θ_b/θ_SD)^0.45`. No dust law, no assumed
T needed — I_bg,b read from the bg band stamp (models) or PEAK^BGF_b (real);
the fractional threshold floor_SD/N_bg is band-independent if the fluctuation is
a column-density fluctuation; the beam term (exponent ≈0.45, forced through
origin) accounts for cirrus smoothing in coarser beams. Reproduces the
independent per-band envelope floors to ~13% (250: 0.86×, 350: 0.88×, 500:
1.13× after beam term). Background temperature per source is obtainable from
PEAK^BGF across bands (falls 16.7→13.9 K with column — warmer than model interior
T_back 13.9→7.0 K, which is column-weighted-SED vs local-at-depth; use the
SED-effective one). Kept as fallback/cross-check, not primary.

### `python/add_recoverable_mass_multiband.py` (NEW, staged)
Per-band frac via floor table or the physical conversion. Superseded
conceptually by tabulate_recovery.py (which tabulates vs r instead of assuming
one floor), but retained. `recoverable_fraction()` now takes the floor as an
argument rather than computing it from the surfdens formula — the one structural
change from the original.

---

## MULTIBAND CHAIN — mechanically validated
Per-band flux correction, then SED refit, then temperature correction:
`F_b/frac_b → refit SED → M_SED → idealised (temperature-bias) correction → M_BE`.
- **The 3-band "bs" SED fit uses 250/350/500** (verified: refit reproduces
  T_SED3bs to 0.006 K; 160/250/350 gives 0.113 K, so it's the SPIRE three).
- Splitting the correction this way isolates the accurate half: the idealised
  temperature correction is ~2% (validated on subdivision models), vs ~10% for
  the entangled recoverable path.
- **Per-band fluxes now in the grid catalogue** (Alexander added them):
  `FXSDbs…FX500bs` at 1-based cols 33–39, plus the bsl set. Catalogue is now
  **67 columns** (was 53). Legend order confirmed: FX{band}bs at 33–39.
- Each band loses a DIFFERENT frac (different beam+noise) → distorts SED colour,
  not just normalisation → biases the fitted temperature. This is real and is
  what per-band correction fixes.

---

## DRAFT CHANGES THIS SESSION (all in `/mnt/user-data/outputs/rt_mass_correction.tex`)
- **NEW `sec:footprint`** — the R_foot degeneracy-breaking finding +
  `Table~\ref{tab:sevencloud}` (7 clouds: d, N_corr, R_foot, f_rec, C_M,
  C_M(0.1–2)), with explicit column definitions in the caption and the plain
  warning that raw C_M is not a distance trend. `\note{TO DO}` on the distance
  proxy / per-distance grids.
- Earlier in the same session (already in prior summary's scope but confirming
  present): `sec:invariant` (distance-invariant corrector), `sec:reliability`
  (f_rec-binned error table), `sec:highsigma` (coverage gap), `app:fields`
  (per-field M(FWHM) + contrast diagnostics, below-grid = beam-diluted
  condensations), Menshchikov2023 bib entry, six→seven cloud fixes.
- Integrity: 25/25 equations, 3/3 tabulars, all refs resolve. Note `\Mrec` is
  NOT a defined macro — draft uses `M^{\mathrm{rec}}_{\mathrm{SED}}`.

---

## VALIDATION STATUS (clean split)
- **Temperature-bias half: VALIDATED ~2.7%** on frac=1 subdivision models (their
  fluxes are measured over the true footprint, so frac=1 by construction —
  validates ONLY the idealised temperature correction, which is the point).
  Fails 50–150% above 0.35 M⊙ but only for the out-of-family critical-BE models.
- **Flux-loss half: converged for Aquila (C_M≈2) via R_foot**, but the
  normalisation still wants the direct injection measurement to be beyond
  dispute. Node-level LOO of the multiband corrector = 28.8% median (17.7% at
  f_rec>0.4) — the realistic per-source uncertainty, worse than the old 14%
  because that assumed the average noise relation held per source (it doesn't).
- **Könyves injection 1.16–1.25 RETRACTED as a benchmark** (Alexander's point,
  re-confirmed): their injected models sat on fainter backgrounds → higher
  contrast → easier measurement → biased low. With it removed, no evidence of a
  systematic model offset remains, and C_M≈2 stands.

---

## CRITICAL-PATH / NEXT STEPS (updated)
1. **aquila-sim getsf extraction (~1 week, Alexander running)** — measures f_rec
   directly for known-flux injected sources; the only assumption-free settle of
   the normalisation. Hold the multi-cloud raw factors until it lands.
2. **Distance problem** — distance-invariant mass proxy OR per-distance grids
   (432, 1700 pc native), before raw multi-cloud C_M can be quoted.
3. **Regenerate v2 tables if the grid bg definition changed again** — Alexander
   updated the bg definition mid-session; current v2 (45200 rows, 518 nodes, 7
   bands, cfoot 5/50/95 = 1.29/1.87/3.22) is on the improved bg. Re-run with
   `--fwhm-floor 1.0` (the measured hard floor is exactly 1.0 beam in all seven
   clouds; the 1.5 I'd guessed is wrong — pile-up at 1.0 is uninformative, a
   secondary cut near 1.3 excludes it) and a cfoot window at analysis time (keep
   full 1.2–3.0 in the file).
4. Injection also gates the frac_rec stochastic floor and the shallow-CMF test.

## KEY LEARNINGS (this session)
- Always say WHAT each factor/fraction is (f_rec, C_M, R_foot defined above).
  Alexander loses the thread otherwise — this is a hard requirement.
- R_foot (footprint/FWHM) is the constraint that breaks the frac degeneracy.
  It's ~1.8–2.0, constant across distance, set by extraction. BE keeps ~45% of
  flux there vs 90% for a Gaussian — that gap is the correction.
- Reported MASS is distance-dependent → a bad matching axis for a fixed-distance
  grid. Use matched mass range (0.1–2) for cross-cloud comparison.
- getsf AFWHM/BFWHM = interpolated half-max = same as model fwhm_rec (once
  interpolated, not binned). Claude was wrong to say getsf "fits" it.
- Source masks clip at surfdens 1e16 → reach 1.25–1.47× past 99%-flux radius →
  sample tables by frac or cfoot, not r/R_conv.
- Envelope detection-floor is percentile-dependent and never converges (it's a
  property of the detected population) → don't rely on absolute floor; R_foot
  sidesteps it.
- Temperature-bias correction (idealised, ~2%) and flux-loss correction (~2×
  for Aquila) are now cleanly SEPARABLE via the multiband chain — this is the
  main methodological gain of the session.

---

## UPDATE 2 — the distance problem is SOLVED (predict C_M, not mass)

This supersedes the "Seven-cloud extension FAILS" and "Distance problem" items
above. Alexander's idea: formulate the corrector to predict the correction
FACTOR C_M = M_BE / M_reported directly, on axes that are all distance-invariant,
instead of predicting an absolute mass from a distance-dependent mass axis.

**Axes (all distance-invariant — ratios or surface brightnesses):**
- Σ = PEAK^BGF (background column)
- conc_peakmean = PEAK^SRC / footprint-mean surface brightness
- R_foot = FOOA / AFWHM (footprint-to-size ratio)

The reported mass is NOT an axis, so nothing in the matching scales with d².
C_M depends only on the source's shape and contrast, which don't depend on
distance. Application: measure the three observables per source → read off C_M →
multiply the reported mass. Same corrector at any distance.

**Result — the backward-running factor is fixed. Seven clouds, invariant C_M:**

| cloud | d(pc) | C_M invariant | C_M(0.1–2) | C_M mass-grid (old, broken) |
|---|---|---|---|---|
| Scorpius | 130 | 3.10 | 3.08 | 2.45 |
| Ophiuchus | 144 | 3.57 | 2.82 | 2.60 |
| Aquila | 260 | 3.41 | 3.69 | 1.87 |
| Orion A | 432 | 3.03 | 3.17 | 1.74 |
| California | 470 | 3.55 | 3.38 | 1.85 |
| Cygnus X | 1150 | 2.93 | 2.85 | 1.26 |
| W3/W4/W5 | 1700 | 3.04 | 2.67 | 0.81 |

Invariant C_M is **flat at 2.9–3.6 with no distance trend** (mass-grid ran
backwards 2.45→0.81). This is the direct demonstration that the correction is a
property of sources+extraction, not distance.

**Validation of the reformulation:**
- Node-level LOO: invariant C_M grid 55% vs mass grid 52% — SAME to within noise.
  Invariance gained at NO cost in accuracy.
- Subdivision models: NOT a valid test of the C_M (flux-recovery) grid — they are
  frac=1 (measured over true footprint), so they have no flux loss to recover;
  passing them through the flux-recovery grid double-corrects them (65% error).
  They validate the IDEALISED temperature-only corrector (2.7%, already done).
  The two validations test the two halves cleanly: subdivision→temperature(2.7%),
  LOO→flux-recovery(same invariant or not).

**IMPORTANT caveats:**
- The invariant level came out ~3, NOT the ~2 the mass version gave for Aquila.
  Reason: dropping the mass axis makes each (Σ,conc,R_foot) cell average over a
  WIDER frac range including the low-frac tail, pulling C_M up. It is a different,
  distance-honest statistic with more scatter — not more wrong than 1.87. Which
  absolute level (~2 vs ~3) is right is STILL the injection question, now cleanly
  ISOLATED from the distance question. Distance-independence won't change with
  the injection normalisation; the level might.
- peak/rim is NOT a usable axis. Claude tried it and got C_M~10–37 (garbage).
  peak/rim in the tables = peak / I_rim(r), a MODEL-INTERNAL truncation coordinate
  (≈ a restatement of frac), with no real-data counterpart. PEAK^SRC/PEAK^BGF is
  a different ratio (~1.2–1.4) and does not match it. Only Σ, conc_peakmean,
  R_foot are genuinely measurable invariant axes. (Claude's earlier "+peak/rim
  improves scatter" test was circular — it smuggled in the truncation.)
- Per-source scatter (~2×) is UNCHANGED by the reformulation. It removes the
  distance contaminant, not the underlying degeneracy floor. So this is the right
  form for POPULATION-level / cross-cloud work; individual masses still uncertain.

**In `fit_mass.py` (staged):** new class `CMInvariantCorrector`.
- `__init__(sigma, conc, rfoot, cm, frac=None)` — build from arrays; cm = M_BE/M_rec.
- `.factor(sigma, conc, rfoot)` → C_M ; `.correct(mass, sigma, conc, rfoot)` →
  mass×C_M ; `.frac(...)` → inferred f_rec if grid carried it. All log10 axes.
- Built from the v2 recovery tables: each (node,band,radius) row gives one
  (Σ, conc, R_foot)→C_M sample with C_M = M_BE/(M_SED·frac).

**Draft (staged rt_mass_correction.tex):** `Table~\ref{tab:sevencloud}` REPLACED
with the invariant C_M version (columns: d, N_corr, R_foot, C_M, C_M(0.1–2),
C_M^mass for contrast). Prose rewritten: flatness = distance-independence
demonstrated; `\note` that absolute level awaits injection. Integrity 25/25 eqs,
tabulars balanced, refs resolve.

**Flowchart (staged, matplotlib PDF — see UPDATE 3):** corrector box now
"(Σ, concentration, R_foot) → C_M, no absolute mass or size"; match box "all
ratios, so distance cancels"; output "same C_M at any distance, median flat at
~3 over 130–1700 pc".

---

## UPDATE 3 — method flowchart created

`fig_method_flowchart.pdf` + `.png` (staged), generated by
`python/make_flowchart_mpl.py` (staged). **matplotlib, true vector PDF — NOT
cairosvg.** History worth keeping: first built as hand-written SVG, but cairosvg
silently ignores baseline-shift and mishandles tspan dy under text-anchor=middle,
so subscripts overlapped/drifted. Rewrote in matplotlib whose mathtext typesets
subscripts correctly. The SVG generator (make_flowchart.py) also exists but the
matplotlib PDF is the one for the paper (matches all other figures' toolchain).

Structure: two converging tracks — CALIBRATION (build corrector once from RADMC-3D
models, M_BE known) and APPLICATION (measure real sources) — meeting at the MATCH
step, then OUTPUT with the C_M = temperature(1.16) × flux-recovery(~2–3)
decomposition, and VALIDATION (temperature 2.7% done, flux injection pending).
To add a subscripted variable to the figure: one entry in the TOK list, e.g.
('X_{sub}', r'$X_{\mathrm{sub}}$').

## REVISED CRITICAL-PATH (supersedes the earlier list)
1. **aquila-sim getsf injection (~1 wk, Alexander running)** — still THE critical
   path, but now settles only the ABSOLUTE LEVEL (~2 vs ~3), no longer the
   distance question. Measures f_rec directly for known-flux sources.
2. Distance problem: **DONE** via CMInvariantCorrector. (Remove old items 1–2.)
3. Regenerate v2 tables only if the grid bg definition changes again.
4. Injection also gates the frac_rec stochastic floor and the shallow-CMF test.
5. Consider building CMInvariantCorrector directly into `collect`/the catalogue
   pipeline so the (Σ,conc,R_foot)→C_M grid ships with the tables.

## KEY LEARNINGS (this update)
- Predict the FACTOR C_M on invariant axes, not the mass — removes distance at no
  accuracy cost. This is now the PRIMARY corrector for multi-cloud work.
- The three real invariant axes are Σ, conc_peakmean, R_foot. peak/rim is model-
  internal and unusable on real data — do not reach for it again.
- Invariant level ~3 vs mass-grid ~2 is a statistic difference (wider averaging),
  not an error. Injection settles the level; distance-independence is already shown.
- Subdivision models validate temperature-only (frac=1); never feed them to the
  flux-recovery corrector.
