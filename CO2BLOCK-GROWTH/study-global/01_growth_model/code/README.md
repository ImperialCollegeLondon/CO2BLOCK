# Growth Model Code

This folder contains the active Part-1 growth-model helper scripts.

## Active scripts

- `build_mcmc_inputs_from_workbook.py`
  - reads `../input/2026_global_raw_data_20260601.xlsx`
  - writes the four technical scenario CSVs:
    - `reference`
    - `minimum`
    - `maximum`
    - `growth10`
- `build_anchored_inputs_v8.py`
  - reads workbook sheet `IPCCreference`
  - builds the four anchored scenario CSVs:
    - `us1gt`
    - `policy`
    - `ipcc_low`
    - `ipcc_high`
- `preprocess_growth_model_inputs.py`
  - retained utility for raw country-file preprocessing
- `run_logistic_mcmc_from_scenarios.py`
  - retained standalone logistic-only runner

The active full Logistic/Gompertz workflow is notebook-led:

- `../notebooks/v8_growth_model_2026-06-01.ipynb`

