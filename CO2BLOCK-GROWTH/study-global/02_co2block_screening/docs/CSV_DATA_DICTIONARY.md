# CSV Data Dictionary

This file explains the main CSV outputs that should be kept for plotting, reporting, and paper figures.

Most active paths below now use the `02_co2block_screening/` layout. Any
remaining `python_version/` references are legacy documentation only.

## 1. Country Capacity CSVs

Produced by:
- Legacy capacity workflow retained only in archived materials; no active script is kept in the screening mainline.

### `country_capacity_static_theoretical.csv`

Purpose:
- Country-level sum of basin static theoretical capacity.

Columns:
- `Country`: country assigned using `Majority Country`
- `Storage_Capacity_Gt`: total static theoretical storage capacity in Gt

Use cases:
- country ranking plots
- comparison to dynamic capacities

### `country_capacity_closed_tank_model.csv`

Purpose:
- Country-level sum of basin closed tank model capacity.

Columns:
- `Country`
- `Storage_Capacity_Gt`

Use cases:
- compare static pressure-limited capacity against theoretical maximum

### `country_capacity_co2block_closed_boundary.csv`

Purpose:
- Country-level sum of dynamic CO₂BLOCK capacities when all basins are forced to closed boundaries.

Columns:
- `Country`
- `Storage_Capacity_Gt`

Use cases:
- dynamic capacity benchmark under conservative boundary assumption

### `country_capacity_co2block_open_boundary.csv`

Purpose:
- Country-level sum of dynamic CO₂BLOCK capacities when all basins are forced to open boundaries.

Columns:
- `Country`
- `Storage_Capacity_Gt`

Use cases:
- dynamic capacity benchmark under optimistic boundary assumption

### `country_capacity_four_types_combined.csv`

Purpose:
- Master country-level comparison table for plotting and reporting.

Columns:
- `Country`
- `Static_Theoretical_Gt`
- `Closed_Tank_Model_Gt`
- `CO2BLOCK_ClosedBoundary_Gt`
- `CO2BLOCK_OpenBoundary_Gt`
- `Open_to_Closed_Ratio`
- `Open_to_StaticTheo_Ratio`

Use cases:
- all country-comparison plots
- supplementary data table
- figure panel construction

Recommended plotting fields:
- x: `Country`
- y: any of the four capacity columns
- compare ratios using `Open_to_Closed_Ratio` and `Open_to_StaticTheo_Ratio`

## 2. Basin Capacity CSVs

### `basin_capacity_static_and_tank.csv`

Purpose:
- Basin-level intermediate table for static capacity calculations.

Columns:
- `SiteNo`
- `Basin`
- `AssignedCountry`
- `Area_km2`
- `Thickness_m`
- `Porosity`
- `CO2Density_kg_m3`
- `TotalCompressibility_1_Pa`
- `PressureLimit_MPa`
- `Static_Theoretical_Gt`
- `Closed_Tank_Model_Gt`

Use cases:
- basin-level supplemental tables
- debugging
- tracing national totals back to basin contributions

### `basin_capacity_co2block_open_boundary.csv`
### `basin_capacity_co2block_closed_boundary.csv`

Purpose:
- Basin-level dynamic CO₂BLOCK outputs for open and closed boundary assumptions.

Columns:
- `SiteNo`
- `Basin`
- `AssignedCountry`
- `BoundaryMode`
- `Max_Q_Mt_per_yr`
- `Optimum_Q_Mt_per_yr`
- `Max_Capacity_Gt`
- `Number_of_Wells`
- `Distance_km`
- `Region_Q_Mt_per_yr`

Use cases:
- identify dominant basins within each country
- case-study figures
- basin-level sensitivity analysis

## 3. Algorithm 2 Precompute CSVs

Produced by:
- `code/run_step1_precompute.py` (or `notebooks/step1_precompute_resource_matrix.ipynb`)

### `basin_period_resource_long.csv`

Purpose:
- Long-format global basin-period resource library for Algorithm 2.

Columns:
- `Region_no`
- `Region_name`
- `AssignedCountry`  (raw majority country, e.g. USA / France / Saudi Arabia)
- `RegionGroup`  (raw region, e.g. Europe / Middle East)
- `ScreeningCase`  (v8 study region the basin rolls up to: US, …, EU, Middle East, Brazil)
- `Period_yr`
- `Region_Q_Mt_yr`
- `Number_of_Wells`
- `Max_Capacity_Gt`
- `Max_Q_Mt_yr`
- `Optimum_Q_Mt_yr`
- `Distance_km`

Use cases:
- flexible plotting by basin, country, or storage period
- country-level aggregation before allocation
- sensitivity analysis over injection duration

This is the most important wide-to-long source table for Algorithm 2.

### `basin_metadata.csv`

Purpose:
- Mapping table for basin-to-country and basin-to-region assignment.

Columns:
- `Region_no`
- `Region_name`
- `AssignedCountry`  (raw majority country)
- `RegionGroup`  (raw region)
- `ScreeningCase`  (v8 study region; step 2 also re-derives this via `resolve_region_nos`)

Use cases:
- country matching
- region matching
- audit trail for country case selection

### `resource_rate_matrix_global.csv`

Purpose:
- Wide matrix of basin sustainable rates across storage periods.

Columns:
- `Region_no`
- `Region_name`
- `Q [Mt/y] for t= 10 y`
- `Q [Mt/y] for t= 20 y`
- ...

Use cases:
- direct Algorithm 2 allocation input
- heatmaps of basin-period rate potential

### `site_number_matrix_global.csv`

Purpose:
- Wide matrix of required site numbers across storage periods.

Columns:
- `Region_no`
- `Region_name`
- `Q [Mt/y] for t= 10 y`
- `Q [Mt/y] for t= 20 y`
- ...

Meaning:
- Each cell gives the number of sites/wells used by the optimum basin-period configuration.

Use cases:
- engineering complexity comparison
- supplementary tables

## 4. Step 2 screening CSVs (current, Path D-full, simplified May 2026)

Produced by `02_co2block_screening/notebooks/step2_allocation_screening.ipynb`
(post-processing cell 9, plotting cells 11 and 13) reading from raw
per-case allocation folders written by `code/run_step2_screening.py`.
Per-ordering files (5 CSVs) live under
`output/step2_screening/{ordering}/final/` where `{ordering}` is
either `ascending` or `descending`. The cross-ordering supplement
(`ordering_sensitivity_deltas.csv`) lives under
`output/step2_screening/extended_data/`.

Two files from earlier versions are no longer regenerated by the current
notebook: `case_model_basin_service.csv` (per-basin service window) and
`delay_sensitivity.csv` (legacy checkpoint diagnostic).

### `screening_summary.csv`

Produced by the allocation core (Stage A). One row per
`(country, scenario, model, allocation_order)`. Does not contain
the whole-period pass verdict; that lives in
`case_model_shortfall_windows.csv`.

| Column | Type | Description |
|---|---|---|
| `allocation_order` | str | `ascend` or `descend` |
| `scenario` | str | Scenario name (`reference`, `maximum`, …) |
| `pathway` | str | Deterministic central pathway identifier for the active V8 growth input |
| `country` | str | Screening case |
| `model` | str | `Logistic` or `Gompertz` |
| `status` | str | `ok` / `no_basin_identified` / `no_growth_path` / `allocation_error: …`; records only whether the allocation core ran without error, not a feasibility pass. |
| `basin_match_basis` | str | The country-basin matching rule applied (e.g. `Majority Country == UK`) |
| `basin_count` | int | Number of basins matched to this country |
| `peak_demand_mt_yr` | float | Maximum annual demand from the input growth curve |
| `allocated_peak_mt_yr` | float | Maximum annual rate the allocation engine actually realised |

### `case_model_shortfall_windows.csv` (primary feasibility table)

Produced by `_compute_whole_period_diagnostics()` in cell 9. One row
per feasible `(country, scenario, model)` case. Carries the
whole-period pass verdict and five severity diagnostics.

| Column | Description |
|---|---|
| `country`, `scenario`, `model`, `allocation_order` | Identifier columns |
| `whole_period_pass` | Primary feasibility flag. `True` iff the screened path is at-or-above demand for every year in 2030-2179. Not equivalent to `status == "ok"` in `screening_summary.csv`. |
| `years_of_shortfall` | Number of years where `gap > 0.01 Mt/yr` |
| `first_shortfall_year` | First year the gap opens (or `None` if pass) |
| `last_shortfall_year` | Last year the gap stays open |
| `shortfall_window_yr` | `last_shortfall_year − first_shortfall_year + 1` |
| `max_gap_rate_mt_yr` | Peak single-year deficit |
| `max_gap_year` | Year the peak deficit happens |
| `total_unmet_volume_gt` | Integrated shortfall (∫ gap(t) dt × 10⁻³) |
| `cumulative_demand_full_gt` | ∫ raw demand over 2030-2179 |
| `cumulative_screened_full_gt` | ∫ screened over 2030-2179 |
| `overall_r_cumulative` | `cumulative_screened / cumulative_demand` over the full horizon |

### `case_model_portfolio_metrics.csv`

Per-checkpoint portfolio summary (legacy checkpoint family, retained
for cross-referencing). One row per
`(country, scenario, model, allocation_order, checkpoint)` where
`checkpoint ∈ {2050, 2100, 2179}`. The third checkpoint at 2179 is
the end of the 150-year allocation horizon and reports the
fully-integrated portfolio state.

| Column | Description |
|---|---|
| `checkpoint` | 2050, 2100, or 2179 |
| `annual_demand_mt_yr` | Demand rate at checkpoint year |
| `annual_screened_mt_yr` | Screened rate at checkpoint year |
| `annual_gap_mt_yr` | Demand − screened at checkpoint |
| `r_annual` | screened / demand at checkpoint |
| `cumulative_demand_gt`, `cumulative_screened_gt`, `cumulative_gap_gt` | Same three quantities integrated 2030 to checkpoint |
| `r_cumulative` | Cumulative-screened / cumulative-demand to checkpoint |
| `annual_pass`, `cumulative_pass`, `overall_pass` | Per-checkpoint pass flags (legacy criterion, not the primary verdict; use `whole_period_pass` instead) |
| `first_shortfall_year`, `years_of_shortfall` | Counted only up to the checkpoint year; for the full horizon see `case_model_shortfall_windows.csv` |
| `max_gap_rate_mt_yr`, `max_gap_year` | Peak deficit up to checkpoint |
| `total_unmet_volume_gt` | Cumulative gap up to checkpoint |
| `active_basins`, `total_basins`, `active_basin_fraction` | Portfolio diversity at checkpoint (active = inside its assigned service window at this year) |
| `total_wells` | Operational well fleet at the checkpoint year: Σ `wells_used` across basins still inside their assigned service window. A snapshot, not cumulative wells ever drilled. The metric tracks the demand curve (rises, peaks ~2100, then declines as service windows close). |
| `n_eff` | Effective basin diversity at the checkpoint: `(Σpᵢ)² / Σpᵢ²` where `pᵢ` is each basin's share of the annual rate at this year |
| `largest_basin_share`, `dominant_basin_name` | Single-basin concentration at the checkpoint |

### `case_model_basin_metrics.csv`

Per (basin, checkpoint) metrics that drive the basin-share proportion
bars in Fig 4. One row per `(case, basin, checkpoint)`.

### `case_model_screened_paths.csv`

Year × case matrix (long format). One row per (case, year) with
`year`, `raw_rate_mt_yr` (demand), `screened_rate_mt_yr` (what the
basin portfolio delivered), and `gap_rate_mt_yr`. Drives the
rate-curve stack in the country-profile figures and supports
whole-horizon shortfall diagnostics.

### `growth_curve.csv` (per-case input mirror)

Country-model-specific demand path as fed into one allocation run.
Stored under `{country_slug}/{scenario}_{model}/growth_curve.csv`.

| Column | Type | Description |
|---|---|---|
| `year` | int | Year (2030 to 2179) |
| `total_rate` | float | Demand rate (Mt CO₂ / yr) |

Use cases:
- traceability: confirm the exact demand curve used for a case
- regenerating the screened path independently of the notebook

### `ordering_sensitivity_deltas.csv` (Q3 supplementary table)

One row per `(country, scenario, model)` case (n = 150) recording
the descending − ascending difference on the screening verdict
plus infrastructure metrics at three checkpoints (2050, 2100,
2179). Stored at
`output/step2_screening/extended_data/ordering_sensitivity_deltas.csv`.

This table is the full per-case audit trail for Q3. It supports:

* `output/step2_screening/fig_ordering_comparison.{pdf,png}`,
  the dedicated Q3 figure (Panel (a) per-case Δ unmet, Panel (b) box-
  and-strip plot of |Δ wells| at the three checkpoints).
* Supplementary tabulation of the asc / desc comparison in the
  paper supplement.

The figure shows the headline statistics; this CSV carries the full
per-case audit trail for any reviewer who wants to reproduce or
re-aggregate them.

Shortfall / verdict columns (independent of checkpoint):

| Column | Type | Description |
|---|---|---|
| `country`, `scenario`, `model` | string | Case identity (join key) |
| `years_of_shortfall_asc`, `years_of_shortfall_desc` | int | Years short under each ordering |
| `total_unmet_volume_gt_asc`, `total_unmet_volume_gt_desc` | float | Integrated shortfall under each ordering (Gt CO₂) |
| `whole_period_pass_asc`, `whole_period_pass_desc` | bool | Pass flag under each ordering |
| `d_years_of_shortfall` | int | `years_of_shortfall_desc − years_of_shortfall_asc` |
| `d_total_unmet_gt` | float | `total_unmet_volume_gt_desc − total_unmet_volume_gt_asc` (Gt CO₂) |
| `pass_flip` | bool | `whole_period_pass_asc ≠ whole_period_pass_desc` (identically False under Path D-full) |

Infrastructure columns repeat at three checkpoints (`CP ∈ {2050, 2100, 2179}`):

| Column template | Type | Description |
|---|---|---|
| `total_wells_asc_{CP}`, `total_wells_desc_{CP}` | int | Operational well fleet at checkpoint year under each ordering (Σ `wells_used` across basins still inside their service window) |
| `n_eff_asc_{CP}`, `n_eff_desc_{CP}` | float | Effective basin diversity (`(Σpᵢ)² / Σpᵢ²`) computed from basin shares of annual rate at the checkpoint |
| `active_basins_asc_{CP}`, `active_basins_desc_{CP}` | int | Number of basins delivering at the checkpoint year |
| `d_total_wells_{CP}` | int | `total_wells_desc_{CP} − total_wells_asc_{CP}`; the main wells-burden difference metric used in Q3 Figure Panel (b) box-and-strip plots |
| `d_n_eff_{CP}` | float | `n_eff_desc_{CP} − n_eff_asc_{CP}` |
| `d_active_basins_{CP}` | int | `active_basins_desc_{CP} − active_basins_asc_{CP}` |

So one row carries checkpoint infrastructure fields plus the
shortfall / identity fields needed to audit the ordering comparison.

Use cases:
- supplementary-material table for the journal paper
- regenerating the per-case scatter without re-loading both orderings
- spotting cases where infrastructure shifts most at each checkpoint
  (`d_total_wells_2050`, `d_total_wells_2100`, `d_total_wells_2179`)
- confirming `pass_flip` is identically False (the headline result)
- tracking how Δ diversity collapses from `median |Δ n_eff_2050|≈1.02`
  to `≈0.18` by 2179 (portfolios converge on basin identity at the
  end of the horizon)

## 5. Allocation Output File

### `Resource_assignment_python.xlsx`

Purpose:
- Step-by-step allocation output for a given country-model run.

Columns:
- `Step`
- `Region no`
- `Region name`
- `Start [y]`
- `End [y]`
- `Duration [y]`
- `Step rate increment [Mt/y]`
- `Step rate cumulative [Mt/y]`
- `No_sites`
- `Curve start [y]`
- `Curve end [y]`

Use cases:
- case-study figures
- stacked-area allocation plots
- narrative explanation of which basins are used first

## Recommended Minimal CSV Set For Plotting

If you want the smallest useful package for later plotting, keep these:

1. `basin_period_resource_long.csv`
2. `basin_metadata.csv`
3. `screening_summary.csv`
4. `case_model_shortfall_windows.csv`
5. `case_model_screened_paths.csv`
6. `case_model_portfolio_metrics.csv`
7. `case_model_basin_metrics.csv`
8. `ordering_sensitivity_deltas.csv`

This set is sufficient for:
- basin resource auditing
- pass/fail scorecards
- country-profile figures
- ordering-sensitivity figures
- supplementary tables
