# MCMC Inputs Derived From Workbook Anchors

This folder contains the scenario-specific input CSVs used by the growth-model
workflow. The CSVs here are the live v8 (2026-06-01) inputs.

Source workbook:
- `01_growth_model/input/2026_global_raw_data_20260601.xlsx`

Country pool (10): UK, US, EU, China, Middle East, Australia, Canada, Indonesia,
Thailand, Brazil. `EU` is sourced from the workbook sheet `EU without UK` but is
labelled `EU` in all CSVs/outputs.

## Maintenance ownership

- Auto-generated from the workbook by
  `01_growth_model/code/build_mcmc_inputs_from_workbook.py`:
  `mcmc_input_reference.csv`, `mcmc_input_minimum.csv`,
  `mcmc_input_maximum.csv`, `mcmc_input_growth10.csv`
- Generated from workbook sheet `IPCCreference` (on top of the reference base
  columns) by `01_growth_model/code/build_anchored_inputs_v8.py`:
  `mcmc_input_us1gt.csv`, `mcmc_input_policy.csv`,
  `mcmc_input_ipcc_low.csv`, `mcmc_input_ipcc_high.csv`.
  These remain reconcilable to the workbook and can also be hand-edited.
- `mcmc_input_ipcc_low.csv` and `mcmc_input_ipcc_high.csv` contain
  identical content by design (dual-file / single-source convention, where
  runtime reads `p2050_lo` for `ipcc_low` and `p2050_hi` for `ipcc_high`).
  Edit both files together when the workbook changes.
- This `README.md` is hand-maintained and is not overwritten by the
  builders.

## Technical scenario rules

- `reference`: `C_scenario = 1 x C_ref`, growth cap `20%` except China `25%`
- `minimum`: `C_scenario = 0.1 x C_ref`, growth cap `10%`
- `maximum`: `C_scenario = 10 x C_ref`, growth cap `20%` except China `25%`
- `growth10`: `C_scenario = 1 x C_ref`, growth cap `10%`

## Anchored scenario rules

All anchored scenarios keep:

- `C_scenario = 1 x C_ref`
- growth cap = reference growth cap (`20%`, China `25%`)

For fixed-target policy scenarios:

- `us1gt` / `policy` use `p2050_lo = p2050_hi = target_fixed`

For countries without workbook-defined targets:

- `p2050_lo` and `p2050_hi` are left blank
- these rows should be treated as unconstrained `Type A`

### `us1gt`

- US fixed at `1.0 Gt/yr`; all other countries (incl. Canada) unconstrained

### `policy`

Published storage targets (workbook `IPCCreference`, column `Published storage target`):

- US fixed at `1.0 Gt/yr`
- EU fixed at `0.33 Gt/yr`
- UK fixed at `0.175 Gt/yr`
- China fixed at `0.216 Gt/yr`
- Canada fixed at `0.06 Gt/yr`  (new in v8; Canada now has a published target)

### `ipcc_low`

- both `p2050_lo` and `p2050_hi` preserve the workbook low/high values
- the `ipcc_low` scenario should read only `p2050_lo`
- low values are:
  - Australia `0.00623012495040894`
  - Canada `0.00178528571128845`
  - China `0.0879427`
  - Indonesia `0.0034121739470087`
  - Brazil `0.0014094`
  - US `0.00363542116016001`
- UK, EU, Middle East, Thailand have no workbook IPCC target, so they are blank (Type A)

### `ipcc_high`

- both `p2050_lo` and `p2050_hi` preserve the workbook low/high values
- the `ipcc_high` scenario should read only `p2050_hi`
- high values are:
  - Australia `0.238633232709833`
  - Canada `0.443912423813518`
  - China `6.57128321180474`
  - Indonesia `3.7238848`
  - Brazil `1.95888264991766`
  - US `3.7238848`
- UK, EU, Middle East, Thailand have no workbook IPCC target, so they are blank (Type A)

## Shared columns

All files use:

- `BASE_YEAR = 2030`
- `END_YEAR = 2430`

The CSV column names are chosen so that the first-stage growth-model notebook can read them directly.
