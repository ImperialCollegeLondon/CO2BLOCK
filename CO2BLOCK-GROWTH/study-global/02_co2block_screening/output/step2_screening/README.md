# Step 2: V8 Whole-Period Screening Outputs

This directory contains the completed V8 basin-screening outputs for the
active 10-region study pool: Australia, Brazil, Canada, China, EU, Indonesia,
Middle East, Thailand, UK, and US. The demand inputs are the fixed-C
deterministic central pathways from
`01_growth_model/output/v8_2026-06-01/timeseries_central.csv`.

The allocation horizon is 2030-2179. Each ordering screens the same 160
scenario-country-model combinations: 150 have an available growth path and
allocation output, and 10 are recorded as `no_growth_path` in
`screening_summary.csv`.

## Main Verdict

The primary feasibility verdict is `whole_period_pass` in
`{ascending,descending}/final/case_model_shortfall_windows.csv`.

| Result | Count |
|---|---:|
| Whole-period pass | 111 / 150 |
| Whole-period fail | 39 / 150 |
| Pass/fail flips between orderings | 0 / 150 |

Ordering does not change the pass/fail verdict. It changes the operational
portfolio, especially the number and identity of active wells and basins.

## Directory Contents

```text
step2_screening/
  fig3_dashboard.{pdf,png}                  # Q1: asc/desc feasibility scorecards
  fig_ordering_comparison.{pdf,png}         # Q3: severity and well-fleet deltas
  ascending/
    final/
      screening_summary.csv
      case_model_screened_paths.csv
      case_model_portfolio_metrics.csv
      case_model_basin_metrics.csv
      case_model_shortfall_windows.csv
    figures/country_profiles/
    {country_slug}/{scenario}_{model}/
  descending/                               # same structure as ascending
  extended_data/
    ordering_sensitivity_deltas.csv
```

## Output Tables

| File | Purpose |
|---|---|
| `screening_summary.csv` | Allocation-core status for every requested run. `status = ok` means the allocation ran; it is not the feasibility verdict. |
| `case_model_shortfall_windows.csv` | Primary whole-horizon verdict and shortfall diagnostics. |
| `case_model_screened_paths.csv` | Annual raw demand, screened delivery, and gap for each successful case. |
| `case_model_portfolio_metrics.csv` | Portfolio diagnostics at 2050, 2100, and 2179. |
| `case_model_basin_metrics.csv` | Per-basin metrics at the checkpoint years, used for country-profile bars. |
| `ordering_sensitivity_deltas.csv` | Descending minus ascending deltas for shortfall and infrastructure metrics. |

## Figures

| Figure | Role |
|---|---|
| `fig3_dashboard.{pdf,png}` | Whole-period feasibility scorecards for ascending and descending orderings. |
| `fig_ordering_comparison.{pdf,png}` | Ordering sensitivity: unmet-volume deltas and well-fleet deltas. |
| `{ascending,descending}/figures/country_profiles/*_profile.{pdf,png}` | Per-country scenario profiles combining Logistic and Gompertz panels. |

## Ordering-Sensitivity Headline

Across the 150 screened cases:

- `pass_flip = False` for all cases.
- 19 / 150 cases differ in `years_of_shortfall`.
- 38 / 150 cases have `|Δ total_unmet_volume| > 0.1 Gt`.
- Median `|Δ total_wells|` is 202 at 2100; the maximum at 2100 is 2 882.

Read ascending and descending as an infrastructure bracket, not as different
feasibility hypotheses.
