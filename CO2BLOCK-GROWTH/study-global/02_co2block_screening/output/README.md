# 02 Output directory

> V8 status (2026-06-01). `step1_precompute/` and
> `step2_screening/` are both regenerated for the active 10-region
> scope: Australia, Brazil, Canada, China, EU, Indonesia, Middle East,
> Thailand, UK, and US. Step 1 contains 123 study basins and Step 2
> uses the V8 fixed-C deterministic central pathways from
> `01_growth_model/output/v8_2026-06-01/timeseries_central.csv`.

Two-step output tree. Both subdirectories are regenerable from the
input data plus notebooks/scripts. The raw per-case allocation folders
under `step2_screening/{order}/{country_slug}/{scenario}_{model}/` are
the slowest to recompute (allocation core, ~10-15 min total); the
post-processing CSVs and figures regenerate in ~2 min if those raw
folders exist.

```
output/
  step1_precompute/                             v8, 10 study regions
    README.md
    precompute/
      basin_period_resource_long.csv            123 basins x 15 periods (+ ScreeningCase)
      basin_metadata.csv                        one row per basin (+ ScreeningCase)
    figures/
      fig1_resource_landscape.{pdf,png}         grouped by study region
      fig2_correlation_heatmap.{pdf,png}

  step2_screening/                              V8 screening results
    README.md                                   Path D-full interpretation
    fig3_dashboard.{pdf,png}                    Q1, 2x2 feasibility scorecards
                                                  (asc/desc rows x Logistic/Gompertz cols)
    fig_ordering_comparison.{pdf,png}           Q3, severity Δ + wells burden Δ
                                                  (separate top-level Q3 figure)
    ascending/                                  smallest-basin-first
      final/
        screening_summary.csv                   one row per case (allocation-core status)
        case_model_screened_paths.csv           year x case (raw + screened rate)
        case_model_portfolio_metrics.csv        per-checkpoint summary (drives Fig 3 Panel c)
        case_model_basin_metrics.csv            per (basin, checkpoint), drives Fig 4
        case_model_shortfall_windows.csv        whole-period verdict + diagnostics (main)
      figures/
        country_profiles/                       Fig 4, one per (country, scenario)
          {country}_{scenario}_profile.{pdf,png}
      {country_slug}/{scenario}_{model}/        raw allocation core outputs:
        growth_curve.csv                        demand path used as input
        Resource_assignment_python.xlsx         basin-by-year assignment table
        Rate_opt_summary_python.xlsx
        Site_no_summary_python.xlsx
    descending/                                 largest-basin-first (same sub-tree)
    extended_data/
      ordering_sensitivity_deltas.csv           Q3 supplement (one row per case,
                                                  Δ wells / Δ n_eff / Δ shortfall at
                                                  each of 2050 / 2100 / 2179)
```

## Scope (V8 paper-ready package)

The screening directory answers three questions, each with its own
figure:

1. Q1: Can the basin portfolio match the growth curve?
   Answered by Fig 3 (2×2 whole-period scorecards, asc/desc rows
   × Logistic/Gompertz cols) plus `case_model_shortfall_windows.csv`.
2. Q2: Why does each case pass or fail?
   Answered by Fig 4 country profiles (per ordering, 76 figures
   each, basin-stack rate curves plus 2050/2100 share bars).
3. Q3: Does ordering change feasibility or only infrastructure?
   Answered by the separate `fig_ordering_comparison.{pdf,png}`
   (2 panels, severity Δ per case plus wells burden Δ at 2050/2100/2179)
   plus `extended_data/ordering_sensitivity_deltas.csv` (full per-case
   audit trail). Headline: feasibility verdict is order-invariant
   (0/150 pass-flips); failure severity moves modestly (19/150 differ
   on shortfall years); operational well fleet differs systematically
   (median |Δ wells| @ 2100 = 202, max = 2 882).

Files no longer kept are listed under "Files removed in earlier cleanup"
below.

## Pass criterion (Path D-full, whole-horizon)

A case passes if the screened path is at-or-above demand for every
year in 2030-2179. The `whole_period_pass` flag from
`case_model_shortfall_windows.csv` is the primary verdict. (Distinct
from `status == "ok"` in `screening_summary.csv`, which only means
the allocation core ran without error.)

## Files removed in earlier cleanup

* `case_model_basin_service.csv` and the `basin_service/` Fig 5
  family: per-basin service-window detail beyond what Fig 4
  shows. Not used in the paper narrative.
* `delay_sensitivity.csv`: legacy checkpoint-based robustness
  diagnostic that Path D-full's whole-period criterion supersedes.
* Old single-ordering Fig 3 dashboards (one each under
  `{ascending,descending}/figures/`): replaced by the 2×2
  scorecard at top level showing both orderings directly.
* The shortfall-timeline panel and the Q3 inset on the old Fig 3.
  Q1 stays in Fig 3 as scorecards only; Q3 moves into the
  dedicated `fig_ordering_comparison.{pdf,png}`.
* Standalone `fig_ED_ordering_sensitivity.{pdf,png}`: retired. Its
  3-checkpoint Δ-wells content is now Panel (b) of
  `fig_ordering_comparison`, and the full per-case table remains at
  `extended_data/ordering_sensitivity_deltas.csv`.

## How to regenerate

1. Stage A (slow, ~10-15 min), only needed if raw per-case
   folders are missing:

   ```bash
   cd 02_co2block_screening/code
   python run_step2_screening.py
   ```

2. Stage B (~2 min), post-processing and all figures:

   ```bash
   cd 02_co2block_screening/notebooks
   # Confirm cell 7 has SKIP_ALLOCATION = True
   jupyter nbconvert --to notebook --execute --inplace \
       step2_allocation_screening.ipynb
   ```
