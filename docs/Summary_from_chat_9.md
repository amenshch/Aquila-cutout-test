# Handoff summary — session 4 (2026-07-23)

Add to project knowledge before the next chat. Companion to
`Summary_from_previous_chat_1..3.md`. Everything below is **decided but not yet
written into `rt_mass_correction.tex`** unless marked IN DRAFT.

---

## 1. Draft audit: what is in, what is missing

### Already in the draft
- MLE (Clauset+2009) CMF slope method + KS-selected `M_min` — `sec:cmfcompare`
- Seven-cloud universality section + figure — `sec:multicloud`, `fig:sixcloud`
- Shallow-CMF / ALMA-IMF hypothesis — `sec:shallowcmf`, flagged `\note{HYPOTHESIS}`
- Cygnus X mass-shift (2–20 → 6–60 M⊙) and its blending caveat — `sec:multicloud`
- W3/W4/W5 distance justification (1706 ± 44 pc, 12 determinations) — `sec:multicloud`
- Field-star / background-galaxy contamination caveat — `sec:multicloud`
- Per-field diagnostics appendix (M(FWHM), contrast, below-grid population,
  band-by-unit-string) — `app:fields`, `fig:mfwhm`, `fig:contrast`
- `frac_rec` irreducible environmental scatter — `\note{TO DO}` in `sec:accuracy`

### **MISSING — highest priority**
**The distance-invariant corrector is never described as a method.**
`sec:multicloud` line ~1980 asserts "the correction developed above is expressed
entirely in distance-invariant observables" — but the method sections describe
the **FWHM-based** corrector. This is an internal inconsistency a referee will
catch. Needs a proper subsection in `sec:concentration` or a new one covering:
- axes = (`M_SED_rec`, Σ, `conc_peakmean`); all three distance-invariant
- `conc_peakmean` = bg-subtracted peak / bg-subtracted footprint mean → a ratio
  of intensities in the same map, so distance cancels *exactly*
- justification from **Men'shchikov (2023)**, `docs/deconv2023.pdf`
  (= `2023_Men'shchikov_aa46152-23.pdf`): Gaussian deconvolution errs by factors
  up to ~20 (unresolved) and ~6 (power-law), i.e. exactly the regime distant
  clouds occupy — so the FWHM path's real-world error exceeds its LOO value
- LOO comparison: FWHM 10.2% / 92% within-2× / 313 scored;
  `conc_peakmean` 14.0% / 82% / **348 scored** (better coverage)
- statement that the invariant path is now the **primary** method, FWHM a check
- **needs a bib entry for Men'shchikov (2023)** — not currently in the .bib

### Other missing items
- **Photon budget closure.** Stale `\note` at line ~641: "Final node count to be
  confirmed after the higher-photon run", and it still says 377 models (now
  742 grid / 517 catalog). Replace with the measured result (§3 below).
- **`bs` vs `bsl` validity distinction** in the SED-fitting methods text
  (only 1 mention of `bsl` in the whole draft).
- **Negative-pixel clipping treatment** and **median-ring pedestal estimator**
  (§4, §5 below) — both decided this session, neither written up.
- Line ~1588 caption cites a "±2% photon-noise floor" — recheck against the
  measured 0.003–0.027 K shift.

---

## 2. CMF slopes — all eight fits use the identical MLE estimator

Single implementation (`mle_slope` / `ks_mmin` / `fit`, currently `/tmp/mle.py`
— **should be moved into a permanent module**). `α = −n/Σln(M_i/M_min)`,
`σ = √n/S`; `M_min` by minimising the KS distance.

| region | d (pc) | published | corrected | M_min pub/corr |
|---|---|---|---|---|
| Könyves Aquila | 260 | −2.22±0.21 (n111) | −1.85±0.24 (n57) | 0.37 / 1.11 |
| Scorpius | 130 | −2.01±0.61 (n11) | −1.19±0.37 (n10) | 0.13 / 0.21 |
| Ophiuchus | 144 | −1.42±0.15 (n95) | −1.04±0.09 (n126) | 0.06 / 0.13 |
| Aquila (getsf) | 260 | −2.92±0.32 (n85) | −1.23±0.09 (n198) | 0.49 / 0.54 |
| Orion A | 432 | −3.63±0.70 (n27) | −2.08±0.25 (n69) | 1.13 / 2.41 |
| California | 470 | −6.29±1.90 (n11) | −1.05±0.08 (n169) | 2.05 / 0.53 |
| Cygnus X | 1150 | −5.40±0.80 (n46) | −3.93±0.59 (n44) | 1.72 / 6.40 |
| W3/W4/W5 | 1700 | −4.01±0.61 (n43) | −2.43±0.32 (n57) | 2.09 / 3.35 |

**Caveat to state in the paper:** `M_min` is selected independently for each
fit, so published and corrected slopes within a region are not fitted over the
same mass interval (California 2.05 → 0.53; Cygnus 1.72 → 6.40 are extreme).
This is the statistically correct procedure but it is *not* a like-for-like
comparison. Worth adding a fixed-`M_min` cross-check before publication.

---

## 3. Photon budget — SETTLED, 10⁹ is adequate

Two independent lines:

1. **10¹⁰ test** (`cats/bes_model_grid_final2_e10_catalog_rec`, 2 nodes).
   The `bs` columns of node (5,18,9) are **corrupt** — its `bs` background
   over-subtracts to all-negative SPIRE fluxes, the fit is undefined, and the
   collect step carried node (3,18,9)'s `bs` values into the row. Use the
   **`bsl`** columns, which are valid for both:

   | node | ξ_max | T_SED3bsl 10⁹→10¹⁰ | M_SED3bsl |
   |---|---|---|---|
   | (3,18,9) | 28.3 | 9.577 → 9.550 (**−0.027 K**) | 11.00 → 11.00 (0.00%) |
   | (5,18,9) | 36.8 | 7.485 → 7.482 (**−0.003 K**) | 11.88 → 11.88 (0.00%) |

2. **Edge-circle scatter is pixelation, not photons.** Images are ray-traced
   through fixed ρ(r), T(r) (Alexander confirmed), so there is *no* photon noise
   in the image at all; azimuthal variation at r≈R_model is the Cartesian grid
   sampling a sharp radial edge at slightly different true radii. More photons
   would not reduce it.

Node-to-node non-smoothness (2nd differences along the mass axis) for context:
σ(T_SED) ≈ 0.025 K flat in ξ; σ(M_SED) grows 0.37% (ξ<2) → 2.35% (ξ>20);
σ(`frac_rec`) 7.4% median / 39% at 90th pct — an order of magnitude noisier than
anything else, and **not** photon-driven (see §6).

---

## 4. ρ_edge < ρ_cloud models — KEEP them, clip the fluxes

82 detectable models (16% of 520), M_BE 0.25–27.6, ξ 5.2–59.9.

**Do not exclude.** Alexander's reasons, which override the earlier reasoning:
- ρ_edge < ρ_cloud describes real collapse + rarefaction-wave cases
- the high ρ_cloud is a *numerical device* to produce a prescribed low T(edge)
  by ISRF shielding, not a claim about real cloud density — so comparing
  ρ_edge to it is comparing a physical to an instrumental quantity
- the constraint was applied long ago and dropped because it shrank coverage

Dropping them costs almost nothing on real data (−2 corrected sources,
−0.01 in the factor), but that is not the reason to keep them.

**Treatment (decided):**
- **per-pixel** zero-clipping of negatives on the background-subtracted
  **160–500 µm** images, *before* integrating over the footprint and fitting
  `bs`. Pixel-wise, **not** a floor on the integrated flux (a deep negative ring
  must not cancel real peak flux before clamping).
- leave **70/100 µm untouched** — negative flux there is genuine absorption
  against a warm background, a diagnostic Alexander built deliberately.
- **preserve both flux sets** as separate columns (raw with negatives, and
  clipped). The difference per band is a depression-depth diagnostic, and
  keeping both documents the treatment for a referee.
- **images:** do not preserve the large negative-area FITS; keep the RADMC-3D
  output so they can be regenerated, plus **one representative** cropped image
  (node 5,18,9) for a paper figure.
- **validation invariant:** `bs ≤ bsl` always (crater ≥ flat interpolation).
  After clipping, check `bs ≤ bsl` at every node; on the depression nodes the
  clipped `bs` should move *up toward* `bsl` without overshooting. Test node
  (5,18,9) first — its `bsl` mass is 11.88.

---

## 5. `bs` background pedestal — median over an inner ring

- Take the pedestal at **r → R_model** (the BE outer boundary), **not** at the
  field edge / cloud value. For ρ_edge < ρ_cloud models the column at R_model
  lies *below* the cloud level, so using the cloud value over-subtracts.
- Use a **median** (or clipped mean), not a plain mean: the contamination is
  one-sided (pixels straddling the sharp truncation read low and would bias a
  mean downward — the same direction as the depression problem).
- Cleanest definition: a **one-pixel-wide ring just inside R_model**
  (R_model − 1.5 px < r < R_model − 0.5 px), so no pixel crosses the truncation.
- **Validation available:** the image is fully deterministic, so the median
  inner-ring pedestal should equal the analytic edge column of the model at
  R_model; the residual is a direct measure of the pixelation error.

---

## 6. `frac_rec` is intrinsically stochastic (Alexander's point)

Not a noisy estimate of one number — a **distribution**. The extraction
interpolates a background under the footprint, and for a fixed intrinsic model
the interpolated pedestal depends on the **local curvature of the fluctuating
background** at the placement site (convex vs concave environments give
different pedestals → different recovered flux and mass). A single-image grid
node samples one realisation, so the grid delivers the **mean** over
environments; the correction is inherently a population-level statement with an
irreducible scatter set by the sky, not by the method or the photon budget.

This most likely explains both the 7.4% node-to-node scatter and why the
recoverable path (~10% LOO) is so much worse than the idealised one (~2%).

**IN DRAFT** as `\note{TO DO}` in `sec:accuracy`, with the required injection
design: ≥20 placements per (Σ, M, FWHM) cell, each site tagged with local
background curvature (Laplacian of the smoothed map), reporting **both** the
mean recovery (validates the grid) and the **distribution width** (the
irreducible per-source floor, plausibly wider than the 8% LOO figure).

---

## 7. Seven-cloud universality — final numbers

Matched **0.1–2 M⊙** factor (raw medians are contaminated by mass sampling —
nearby clouds detect sub-0.1 M⊙ cores near the floor, inflating the median):

| cloud | d (pc) | matched factor | corrected N | H2 band |
|---|---|---|---|---|
| Scorpius | 130 | 1.51 (n12) | 38 | 02 |
| Ophiuchus | 144 | **2.37** (n42) | 377 | 02 |
| Aquila (getsf) | 260 | 1.55 (n437) | 815 | 02 |
| Orion A | 432 | 1.51 (n468) | 574 | 02 |
| California | 470 | 1.37 (n251) | 282 | **03** |
| Cygnus X | 1150 | 1.60 (n690) | 743 | **03** |
| W3/W4/W5 | 1700 | 1.37 (n165) | 215 | 02 |

Six of seven cluster **1.37–1.60 across a 13× distance range with no distance
information used**. Ophiuchus is an *environmental* outlier (L1688 density +
Sco OB2 illumination), not a distance effect — Scorpius at the same distance is
normal.

Mass-dependence in the high-mass regions: Cygnus X 1.60 (0.1–2) → 2.61 (1–10)
→ 3.09 (2–20) M⊙; W3/W4/W5 1.37 → 1.63 → 1.86.

---

## 8. Tools and files

- **`python/getsf_catalog.py` (NEW, staged).** Robust getsf reader: finds
  columns by **header name**, identifies the 13.5″ H2 band by unit string
  (`H2/cm^2`) + beam + PEAK^BGF magnitude. Tested on all seven catalogs —
  correctly returns band 02 for five, **03 for California and Cygnus X** (both
  have a sub-160 band that shifts the numbering). **Should replace the
  hardcoded-band parsing in the remaining analysis scripts.**
- **`python/fit_mass.py`** — `InvariantCorrector` added (docstring carries the
  Men'shchikov 2023 rationale and the LOO comparison).
- `/tmp/mle.py` — MLE helper; **move somewhere permanent.**
- Figures staged: `fig_sixcloud_universality`, `fig_multicloud_cmf`,
  `fig_mfwhm_allfields`, `fig_contrast_allfields` (all .pdf + .png).
- Cloud catalogs fetched into `<Cloud>-Guoyin/` (Aquila, OrionA, California,
  Ophiuchus, Scorpius, CygnusX, W3W4W5=TriRegion).

**Catalog distance gotcha:** Guoyin's SED fits sometimes used a stale DISTANCE.
Aquila was fit at 260 pc despite a "484" prefix (Alexander: intended, non-issue).
W3/W4/W5 was initially fit at 130 pc and **has since been corrected to 1700 pc**
— the current file is the fixed one. Always check the DISTANCE column against
the header prefix on any new catalog.

---

## 9. Open items / next steps

1. **Alexander is regridding** with the corrected background subtraction (§4, §5)
   and will supply an updated catalog. Validation sequence agreed:
   - `bs ≤ bsl` at every node
   - depression nodes: clipped `bs` moves up toward `bsl`, no overshoot
     (check node 5,18,9 first)
   - re-derive the correction and confirm stability of: LOO accuracy (~8%),
     Aquila factor (1.59–1.65), the seven matched factors above
   - second-difference smoothness test: the median-ring pedestal should reduce
     node-to-node scatter vs the old single-pixel estimator
   - rebuild `InvariantCorrector` on the new catalog, re-run all seven clouds,
     regenerate the universality table and figures
   - flag if the column layout / recovery columns changed in the regrid
2. **Write the distance-invariant corrector into the draft** (§1) — top priority,
   plus the Men'shchikov 2023 bib entry.
3. Close the photon `\note`; update the stale 377-model text.
4. Add the `bs`/`bsl` validity distinction and the §4/§5 treatments to methods.
5. Three j20/k10 nodes still to complete the grid (5 GB images; expected to
   change nothing — removing the analogous i01j20k10 node changed results only
   in the 3rd decimal).
6. **Injection experiment** — now gates three things: the `frac_rec` scatter
   floor (§6), end-to-end accuracy on real backgrounds, and the shallow-CMF
   claim (inject known Salpeter / shallower / steeper slopes and verify the
   recovered slope matches the input; if Salpeter in → Salpeter out, the
   shallow real-data slope is a property of the cores).
7. Fixed-`M_min` cross-check on the CMF slopes (§2 caveat).
8. Consider moving `fig_sixcloud_universality` earlier as a headline figure.
9. Consider native 432 pc / 1700 pc grids to remove hull losses for distant
   clouds (now the *only* remaining resolution penalty, since the invariant
   path removed the deconvolution one).
10. More Guoyin clouds available if wanted.

---

## 10. Paper division (the living document feeds several)

- **Paper I — methods:** RT grid + mass correction + validation.
- **Paper II — universality:** seven clouds, 130–1700 pc, plus the
  **shallow-CMF / ALMA-IMF** hypothesis. Needs the injection test first.
- **Paper III — diagnostics:** blending / BE contrast ceiling.

`\note{}` markers double as the to-do list and as paper-division signposts.
