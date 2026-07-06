# 03: Feasible Growth (Part 3): design notes

Design notes for a notebook-first analysis.

## Current V8 status

This module has been rerun for v8 (2026-06-01) from the current Part-1
pathways and Part-2 screening results. The reproducible pipeline is in
`code/` (`build_part3_input_table.py` → `build_part3_reconstruction.py` →
`build_part3_rescreen.py`); the notebook `notebooks/part3_feasible_growth.ipynb`
loads those outputs and renders the figures.

v8 scope (active): 39 failed `(country, scenario, model)` cases give 78
ordering-specific reconstruction targets across 8 countries (Australia, Brazil,
Canada, China, EU, Thailand, UK, US). All 78 solve; all 78 pass re-screening.

The notes below remain the method-design spec. One sourcing fix relative to
the original notes: the input table is built from `timeseries_central.csv` +
`mcmc_input_*.csv`, not `scalars_central.csv` alone, because scalars omits 2
model-paths whose sibling model's central was infeasible.

## 1. Purpose

Part 3 is the final step after the first two parts. Part 1 defines
deterministic national growth pathways. Part 2 tests whether those pathways can
be delivered by basin performance under ascending and descending allocation
orderings. Part 3 reconstructs a smooth, ordering-specific pathway for the
failed Part-2 cases, then screens that reconstructed curve again to recover the
basin allocation pattern.

The intended role of Part 3 is therefore:

> For each failed ordering-specific case, derive a smooth S-curve
> reconstruction anchored to the same 2030 cumulative amount and the
> same 2050 annual rate used in Part 2, while matching the annual
> peak rate actually delivered by that ordering in Part 2; then
> rescreen that reconstructed curve to see how basin allocation is
> re-organised under the same ordering.

So Part 3 is ordering-specific (ascending and descending are reconstructed
separately) and scenario-aware (the preserved quantity depends on the original
scenario expression, and the 2050 anchor is inherited from the Part-1 / Part-2
deterministic pathway). It is reconstructive rather than purely certifying: the
main output is a smooth screened-pathway reconstruction, not a direct
year-by-year feasible solve. And it is closed-loop, because the reconstructed
curve is screened again with the Part-2 allocation engine so the basin
distribution remains visible.

Part 3 should therefore not be described as pure `fixed-C only`. In this revised
version, `C` becomes a solved quantity for the main reconstruction rather than a
universally fixed one.

## 2. Analytical unit and main question

The analytical unit is:

```text
(country, scenario, model, ordering)
```

That is, a failed ordering-specific Part-2 case. So if the same `(country,
scenario, model)` fails under both orderings, Part 3 reconstructs two separate
smooth pathways, one for ascending and one for descending. This is necessary
because Part 2 already showed that ascending and descending can produce
different screened paths, different severity, and different basin portfolios,
even when the pass/fail verdict is unchanged.

The main Part-3 question:

> Given the Part-2 screened pathway under a specific ordering, what
> smooth same-family S-curve best reconstructs that pathway when
> anchored to the same 2030 cumulative amount, the same 2050 annual
> rate, and the same achieved annual peak rate?

This differs from the earlier fixed-`C` framing, which asked what the slowest or
fastest curve was that remained below the full yearly screened envelope. That
version was mathematically strict, but for several failed Logistic cases it
generated extremely degenerate curves with very small `g`, because fixed `C`
plus the late-horizon envelope forced the curve almost flat.

The revised framing instead aims to recover a smooth interpretable
reconstruction of the screened pathway. It preserves the same 2030 anchor and
the same 2050 target used in Parts 1–2, uses the Part-2 achieved annual peak
rate as the third landmark, then evaluates how well that reconstruction agrees
with the full yearly envelope and screens it again to recover basin allocation.
So the main output is best described as an ordering-specific,
screened-landmark-matched reconstruction rather than a purely analytical
`fixed-C` feasible solve.

## 3. Scenario-aware framing

Part 3 is still scenario-aware, but the meaning has changed relative to the
earlier fixed-`C` notes. For the main reconstruction, all failed
ordering-specific cases preserve `S0` (the 2030 cumulative anchor from Part 1),
`rate_2050_original` (the same 2050 annual rate used in Parts 1 and 2),
`R_peak_order` (the annual peak rate actually delivered by the corresponding
Part-2 screened path under that ordering), and the same curve family (`Logistic`
or `Gompertz`).

Scenario-aware interpretation by family:

| Scenario family | Main interpretation in Part 3 |
|---|---|
| `reference`, `minimum`, `maximum` | Reconstruct a smooth screened counterpart; compare solved `C_recon` back to `C_scenario` |
| `growth10` | Same as above, but additionally flag whether the reconstructed `g_recon` exceeds `g_cap` |
| `Type == A` rows in `policy`, `us1gt`, `ipcc_*` | Interpret the reconstruction as an ordering-specific smooth counterpart to the screened path, not as an externally imposed policy target exercise |
| `Type == B_fixed` rows | The same reconstruction has an additional policy-facing meaning because the 2050 rate is externally anchored in the scenario definition |

So Part 3 remains scenario-aware because the 2050 anchor is inherited from the
original scenario-specific deterministic pathway, the interpretation of the
solved `C_recon` depends on the scenario family, and `growth10` keeps its cap as
a diagnostic consistency check.

If needed later, a true target-preserving / resource-concession extension can
still be added for `Type == B_fixed` cases only. But that should not be the main
notebook path in the first implementation: the clean main workflow is the
landmark-matched reconstruction plus rescreening.

## 4. Reconstruction definition

For each failed ordering-specific case, define three preserved landmarks:

1. `S0`, the 2030 cumulative amount from Part 1
2. `R2050`, the 2050 annual rate from the deterministic Part-1 / Part-2 pathway
3. `R_peak_order`, the annual peak rate achieved by the screened path in Part 2
   under the same ordering

where:

```text
R_peak_order = max_y screened_rate_order(y)
```

and `screened_rate_order(y)` comes from the corresponding ordering-specific
`case_model_screened_paths.csv`.

With the model family fixed and `S0` fixed, the unknowns are `C_recon_order` and
`r_recon_order`, and `g_recon_order` is then recovered through the Part-1
`r ↔ g` mapping. For each ordering-specific failed case solve:

```text
rate_at_2050(model, C_recon, S0, r_recon) = R2050
peak_rate(model, C_recon, r_recon)        = R_peak_order
```

Then recover:

```text
g_recon_order = r_to_g(model, C_recon, S0, r_recon)
```

This is the core Part-3 reconstruction.

### 4.4 Branch selection rule

For Logistic and Gompertz, this two-equation system can in principle admit more
than one mathematical branch. The notebook imposes a deterministic selection
rule. First, reject branches with `peak_year_recon < 2050` unless no post-2050
branch exists. Then, if more than one post-2050 branch remains, choose the
branch with the smaller envelope inconsistency (defined in §5). This keeps the
reconstruction aligned with the study's interpretation of a post-2030 growth
pathway.

Flag a reconstruction as `landmark_infeasible = True` if no real positive
`(C_recon, r_recon)` satisfies both equations, or if every admissible branch has
a pathological shape (for example, invalid timing or non-physical
`C_recon <= S0`). This is different from the earlier `no_solution_in_bracket`
issue in the fixed-`C` notebook: here the problem is the landmark system itself,
not an arbitrary lower bound on `g`.

## 5. Envelope consistency check

The reconstructed curve is not automatically guaranteed to remain under the full
yearly screened envelope, because it is matched only to the 2050 rate and the
peak rate. So a second layer of diagnostics is required.

For each failed ordering-specific case define the ordering-specific yearly
envelope:

```text
E_order(y) = screened_rate_order(y)
```

using the Part-2 screened path for that ordering. For the reconstructed smooth
curve `rate_recon_order(y)`, compute:

```text
diff_order(y) = rate_recon_order(y) - E_order(y)
```

and store `max_envelope_exceedance_mt_yr`, `years_exceeding_envelope`,
`first_exceed_year`, `last_exceed_year`, `rmse_to_envelope_mt_yr`, and
`peak_year_gap = peak_year_recon_order - Y_peak_order`. These answer how
faithfully the smooth landmark-matched reconstruction reproduces the full
screened envelope.

Interpreting the result: if `years_exceeding_envelope = 0`, the reconstruction
is fully consistent with the ordering-specific envelope. If exceedance is small
and local, it remains a useful smooth interpretation of the screened path, but
should not be called fully feasible without qualification. If exceedance is
large, the reconstruction is informative as a smoothed landmark fit, but the
rescreening step will show that it is not fully deliverable under the same
ordering.

## 6. Screening again: basin reallocation on the reconstructed curve

This is the extra layer requested for Part 3. Once the ordering-specific smooth
curve is reconstructed, run the screening allocation again under the same
ordering to recover whether the reconstructed curve now passes or still fails,
how the basin allocation changes under this smoother curve, and how concentrated
or redistributed the basin portfolio becomes.

For each reconstructed ordering-specific case, take the reconstructed annual
curve `rate_recon_order(y)` as the new national demand path, keep the same
country, scenario, model, and ordering, rerun the Part-2 allocation core with
the same Step-1 basin library and the same ordering rule, then rebuild the
screened path and summary metrics exactly as in Part 2.

The reconstruction itself only gives a smooth path; the re-screening step tells
the paper:

> If we replace the failed deterministic path with its screened
> landmark-matched reconstruction, what basin portfolio would actually
> deliver it under the same ordering?

This is the right place to bring basin distribution back in, without turning
Part 3 into another oversized screening section.

Store at least `rescreen_pass`, `rescreen_total_unmet_gt`,
`rescreen_years_of_shortfall`, `rescreen_total_wells_2050`,
`rescreen_total_wells_2100`, `rescreen_active_basins_2050`, and
`rescreen_active_basins_2100`, plus one raw assignment table
`reconstructed_assignment.csv` containing the basin-by-basin allocation of the
reconstructed curve. This raw table can stay as an internal output, while the
summary CSV and country-specific figures carry the main story.

## 7. Inputs

The notebook reads a small fixed set of tables.

From Part 1:

- `01_growth_model/output/v8_2026-06-01/scalars_central.csv`
- `01_growth_model/output/v8_2026-06-01/timeseries_central.csv`

used for `S0`, `rate_2050_original`, `g_original`, `r_original`, `tn_original`,
`tp_original`, `peak_rate_original`, `C_scenario`, `Type`, and `g_cap`.

From Part 2:

- `02_co2block_screening/output/step2_screening/ascending/final/case_model_shortfall_windows.csv`
- `02_co2block_screening/output/step2_screening/descending/final/case_model_shortfall_windows.csv`
- `02_co2block_screening/output/step2_screening/ascending/final/case_model_screened_paths.csv`
- `02_co2block_screening/output/step2_screening/descending/final/case_model_screened_paths.csv`
- the allocation-core inputs needed to rerun the Part-2 screening engine on a
  supplied national time series

For the failed-case universe, keep rows with:

```text
whole_period_pass = False
```

and do not collapse asc and desc. So the main Part-3 table is one row per
`(country, scenario, model, ordering)`. Under the v8 Part-2 outputs this is
exactly 39 failed `(country, scenario, model)` cases and 78 failed
ordering-specific rows (asc/desc fail sets identical; 0 pass-flips).

## 8. Notebook logic

One notebook: `03_feasible_growth/notebooks/part3_feasible_growth.ipynb`.

Cell A, setup and loading: load Part-1 and Part-2 tables, build the failed
ordering-specific universe, attach `Type`, `g_cap`, `C_scenario`,
`rate_2050_original`.

Cell B, curve helpers: reuse or lightly port the Part-1 helpers
`forward(model, C, S0, g, years)` (or an equivalent `r`-based forward
evaluator), `rate_at_year(...)`, `peak_rate_of_curve(...)`, `r_to_g(...)`,
`g_to_r(...)`, `peak_year_of_curve(...)`.

Cell C, ordering-specific screened landmarks: for each failed ordering-specific
case, extract `E_order(y)` from Part 2, compute
`R_peak_order = max_y E_order(y)` and `Y_peak_order = argmax_y E_order(y)`, and
store the full envelope for later comparison.

Cell D, solve the smooth reconstruction: for each failed ordering-specific case,
solve `(C_recon_order, r_recon_order)` from `rate_at_2050(...) = R2050` and
`peak_rate(...) = R_peak_order`, apply the branch-selection rule from §4.4,
recover `g_recon_order`, and generate `rate_recon_order(y)` and the cumulative
series.

Cell E, envelope consistency diagnostics: for each reconstructed case, compare
`rate_recon_order(y)` against `E_order(y)`, compute exceedance counts and
magnitudes, and classify `envelope_consistent` (True/False) plus an optional
softer `small_local_exceedance` label.

Cell F, screening again: for each reconstructed ordering-specific case, pass the
reconstructed annual curve into the Part-2 allocation core, rerun screening
under the same ordering, and capture pass/fail, unmet volume, wells, active
basins, and the assignment table.

Cell G, outputs and figures: write the summary tables and the two figure
families below.

## 9. Outputs

`smooth_reconstruction_summary.csv`, one row per failed ordering-specific case
`(country, scenario, model, ordering)`. Recommended core columns: identifiers
(`country`, `scenario`, `model`, `ordering`, `Type`); original values (`S0_gt`,
`C_scenario_gt`, `g_original`, `r_original`, `rate_2050_original`,
`peak_rate_original`, `tp_original`); screened landmarks (`R_peak_order`,
`Y_peak_order`); reconstructed curve (`C_recon_order_gt`, `g_recon_order`,
`r_recon_order`, `rate_2050_recon_order`, `peak_rate_recon_order`,
`peak_year_recon_order`, `tp_recon_order`); comparisons
(`C_ratio_recon_to_scenario`, `delta_tp_recon_vs_original`,
`delta_peak_year_recon_vs_screened`); envelope consistency
(`max_envelope_exceedance_mt_yr`, `years_exceeding_envelope`,
`first_exceed_year`, `last_exceed_year`, `rmse_to_envelope_mt_yr`,
`envelope_consistent`); and scenario-aware flags (`violates_g_cap`,
`landmark_infeasible`).

`smooth_reconstruction_timeseries.csv`, long format. Columns: `country`,
`scenario`, `model`, `ordering`, `year`, `curve_kind`, `rate_mt_yr`, `cum_gt`.
Recommended `curve_kind` values: `original`, `screened_envelope`,
`reconstructed`, and optionally `rescreened_delivered`.

`reconstructed_rescreen_summary.csv`, one row per reconstructed
ordering-specific case, containing `country`, `scenario`, `model`, `ordering`,
`rescreen_pass`, `rescreen_total_unmet_gt`, `rescreen_years_of_shortfall`,
`rescreen_total_wells_2050`, `rescreen_total_wells_2100`,
`rescreen_active_basins_2050`, `rescreen_active_basins_2100`.

`reconstructed_assignment.csv`, the raw basin allocation table from the
screening-again step. This is mainly for internal traceability and
country-specific figure construction.

## 10. Figures

Three figure families with explicit paper-vs-supplementary roles:

| Figure | Role | Paper / supplementary |
|---|---|---|
| Fig P3-3 | Country-level feasibility summary (forest + sparklines) | Paper main |
| Fig P3-1 | Envelope-fit honesty diagnostic (single panel, all 78 ordering-specific cases) | Supplementary |
| Fig P3-2 | Per-(country, scenario) reconstruction profile (23 figures) | Supplementary |

P3-1 was originally designed as a two-panel scatter, but its first panel
(`C_scenario` vs `C_recon` log-log) duplicates the country-level story of P3-3
panel (a) in a less readable form. The active V8 figure therefore keeps only the
envelope-fit panel (see §10.1), which makes the figure architecture
redundancy-free.

### Figure P3-1: Envelope-fit honesty diagnostic (supplementary)

Single panel. x = `years_exceeding_envelope` (count over 2030–2179), y =
`max_envelope_exceedance_mt_yr`, colour = ordering (`ascending` /
`descending`), marker = model (circle Logistic, triangle Gompertz). It reads as:
how closely does the landmark-matched smooth reconstruction agree with the full
ordering-specific yearly screened envelope?

This is the honesty figure that answers the obvious reviewer question: you said
your reconstruction is anchored to three landmarks, so how much does the rest of
the curve drift from what Part 2 actually delivers? The answer is the
distribution of exceedance years times peak exceedance magnitude across all 78
ordering-specific reconstructed cases.

Subtitle statistics: number of fully envelope-consistent cases (`years_exceeding
== 0`), median `years_exceeding_envelope`, median `rmse_to_envelope_mt_yr`.

The previous `C_scenario` vs `C_recon` panel is omitted because P3-3 panel (a)
presents the same information per country with median markers and min/max
whiskers, which is far more readable than a log-log scatter of ordering-specific
points. Saved to `output/figures/fig_envelope_consistency.{pdf,png}`.

### Figure P3-2: Country-specific reconstruction and basin replay

One figure per failed `(country, scenario)` group, with the left panel Logistic
and the right panel Gompertz.

The upper main axis of each panel shows the smooth reconstruction: the original
Part-1 deterministic curve (grey dashed), the ascending screened envelope (blue
thin dashed), the descending screened envelope (red thin dashed), the
reconstructed ascending smooth curve (blue solid), and the reconstructed
descending smooth curve (red solid). Markers sit at the 2050 point on each
reconstructed curve and at the reconstructed peak point.

The lower inset/strip shows the basin replay after screening again: a compact
basin-composition summary from the re-screening step, for example one stacked
horizontal bar for ascending at 2100 and one for descending at 2100, using top 3
basins plus `other`. This is enough to show how the reconstructed curve
redistributes basin use under each ordering, without reviving the full Part-2
plot family.

A compact text box should report only the highest-signal quantities:

```text
orig : rate2050 / peak / tp / C_scenario
asc  : C_recon / peak_year / exceed_years / rescreen_pass
desc : C_recon / peak_year / exceed_years / rescreen_pass
```

Do not add wells, basin-service timelines, or cumulative panels unless a
concrete paper need appears later.

### Figure P3-3: Country-level feasibility summary (paper main figure)

P3-1 is a technical scatter and P3-2 is a per-(country, scenario) diagnostic
family. Neither is paper-main-text-optimal on its own. Fig P3-3 fills that role:
one composite figure comparing feasible vs original growth at the country level,
suitable for inclusion in the paper main text.

Three design options were documented; only one is implemented, and the default
is Option 1.

Option 1, forest plot plus curve sparklines (default), a two-panel composite
~22 cm × 12 cm. Panel (a), a "feasibility cost" forest plot, with y axis = the 7
failed countries (Australia, Canada, China, Thailand, UK, US, plus any others
that appear in the failed universe), x axis = `C_recon / C_scenario` ratio (0 to
1), two rows per country (Logistic and Gompertz), each row showing all
failed-scenario points plus median marker plus min/max whiskers, colour =
ordering, and a reference line at 1.0 (no reduction). It reads as: which
countries lose the largest share of their long-term storage commitment under
geological reality, decomposed by model and ordering? Panel (b), original vs
feasible curve sparklines in a two rows × N countries grid, with the top row
(Logistic) showing each country's Logistic representative scenario and the
bottom row (Gompertz) showing each country's Gompertz representative scenario;
each sparkline plots original (grey dashed), feasible_asc (blue solid), and
feasible_desc (red solid). The representative scenario is picked per (country,
model), so Logistic and Gompertz may show different scenarios for the same
country (the two families fail in different scenario sets).

Representative-scenario rule: if `reference` is a failed scenario for this
(country, model), use `reference`; else pick the failed scenario with smallest
`C_recon / C_scenario` ratio (worst geological constraint); if no failed
scenarios exist for this (country, model), show a "passed Part 2 in every
scenario" placeholder. This is more honest than picking one scenario per country
(and showing only the model that failed for it), because both model families are
first-class citizens in the screening framework. Panel (b) reads as: what does
the loss of commitment look like as a curve, and how does the reconstruction
reshape the rate trajectory?

The trade-off: Option 1 combines the quantitative panel (a) with the visual
panel (b), is headline-friendly (both "how much" and "what does it look like"),
and fits one paper-press-friendly figure. The one methodological choice is
picking one representative scenario per country in panel (b), which is governed
by the rule above.

Option 2, a country × scenario heatmap, a single dense panel with rows = 7
countries × 2 models = 14 rows, columns = 8 scenarios, cell colour =
`C_recon / C_scenario` ratio (red gradient = heavy reduction, green =
preserved), hatching or split-cell encoding for ordering (`/` asc, `\` desc),
and passed cells in grey. It is comprehensive (every case visible at once) but
dense, better suited to Supplementary Material, and encoding ordering as
hatching can be hard to read at small print sizes.

Option 3, per-country small multiples, a 7-panel grid (one per failed country),
each panel showing original (grey dashed) plus feasible_asc plus feasible_desc
curves for all failed scenarios in that country, curves coloured/styled by
scenario, with a compact per-panel annotation (median ratio, count of failed
scenarios). It is the most visual and story-driven, but panels for high-failure
countries (Australia, Canada, China, each with 6+ failed scenarios) become busy,
and per-scenario colour encoding within each panel may make scenarios hard to
distinguish.

Option 1 was chosen because the forest plot in panel (a) gives the quantitative
paper-claim number at a glance ("country X loses Y% of its commitment"), the
sparklines in panel (b) give the curve-shape intuition that supports the
headline number, and both stories live on one ~22 × 12 cm figure that fits the
paper main text without compression. Save to
`output/figures/fig_country_feasibility_summary.{pdf,png}`. Option 2 and Option
3 are documented as fallbacks; switch only if Option 1 turns out to obscure a
story-relevant comparison.

## 11. What the paper can claim

Part 3 supports these statements. Failed Part-2 pathways can be translated into
ordering-specific smooth reconstructions anchored to the same 2030 cumulative
amount, the same 2050 annual rate, and the same screened annual peak rate. The
solved long-term resource base `C_recon` provides an interpretable bridge
between the original scenario definition and the basin-deliverable screened
path. Agreement with the full yearly envelope is not automatic and must be
reported explicitly. Re-screening the reconstructed curve reveals how the basin
distribution changes under the same ordering once the path is smoothed. And
ascending and descending remain meaningfully distinct in Part 3, because they
provide different screened envelopes and different re-screened basin portfolios.

## 12. What Part 3 should avoid

Don't call the landmark-matched reconstruction automatically "fully feasible"
without checking the full yearly envelope. Don't collapse ascending and
descending into one curve. Don't re-run MC, introduce new basin-ordering
variants, bring back large basin-service figure families, or let the
re-screening outputs dominate the section. Part 3 is smooth reconstruction plus
allocation replay, not a second full screening chapter.

## 13. Folder structure

```text
03_feasible_growth/
  README.md                                             formal results README
  README_NOTES.md                                       this file
  notes/
    REVIEW_NOTES.md                                     future review/audit memos
  docs/
    CSV_DATA_DICTIONARY.md                              column schemas
  notebooks/
    part3_feasible_growth.ipynb                         the one notebook
  output/
    smooth_reconstruction_summary.csv                   main analytical table
    smooth_reconstruction_timeseries.csv                long-format curves
    reconstructed_rescreen_summary.csv                  screening-again summary
    reconstructed_assignment.csv                        internal traceability
    figures/
      fig_country_feasibility_summary.{pdf,png}         Fig P3-3 paper main
      fig_envelope_consistency.{pdf,png}                Fig P3-1 (supplementary)
      country_profiles/
        {country}_{scenario}_reconstruction.{pdf,png}   Fig P3-2 (per-case)
```

## 14. One-sentence workflow summary

> Part 3 takes each failed Part-2 ordering-specific case, solves a
> smooth S-curve reconstruction using the same 2030 anchor, the same
> 2050 annual rate, and the same screened annual peak rate, checks how
> closely that curve matches the full screened envelope, and then runs
> screening again to recover the basin distribution under the same
> ordering.
