# 01 Growth Model

National- and global-scale CO2 storage growth modelling under eight
bottom-up scenario families, using Logistic and Gompertz S-curve models.

The active workflow is the v8 2026-06-01 rerun:

- source workbook: `input/2026_global_raw_data_20260601.xlsx`
- active notebook: `notebooks/v8_growth_model_2026-06-01.ipynb`
- active scenario inputs: `output/mcmc_inputs/`
- active outputs: `output/v8_2026-06-01/`

The v8 method uses the validated Logistic/Gompertz equations, sampling logic,
acceptance filters, and global aggregation logic with the final workbook input
set. `EU` is sourced from workbook sheet `EU without UK`, Brazil is included
in the active country pool, and Canada carries a published policy target.

## Active Scope

This module generates country-level growth pathways: annual storage-rate
time-series from 2030 to 2180 for each scenario, country, and growth model.
These pathways are the Part-1 input to the later basin-screening workflow.

Part 1 does not test basin deliverability. That remains the role of
`02_co2block_screening`.

The active 10-country / regional pool is:

| Output label | Workbook sheet | Notes |
|---|---|---|
| UK | `UK` | published policy target available |
| US | `US` | US 1 Gt and policy target available |
| EU | `EU without UK` | output remains `EU`; methods/captions should state EU excludes UK |
| China | `China` | reference and maximum growth cap is 25% |
| Middle East | `Middle East` | no workbook IPCC/policy target |
| Australia | `Australia` | IPCC low/high targets available |
| Canada | `Canada` | IPCC low/high and policy target available |
| Indonesia | `Indonesia` | IPCC low/high targets available |
| Thailand | `Thailand` | no workbook IPCC/policy target |
| Brazil | `Brazil` | active v8 country |

## Method

For each `(country, scenario, growth model)` combination, the workflow:

1. Samples `N = 10,000` candidate parameter sets.
2. Solves or evaluates the implied S-curve trajectory.
3. Keeps only samples satisfying the explicit feasibility criteria.
4. Aggregates accepted country samples into global pathways by Monte Carlo
   recombination.
5. Reports global and country-level rates, cumulative storage, required
   resource, and timing diagnostics.

The workflow is bottom-up throughout: samples are generated at country level,
then combined to global pathways. AR6 categories are used only as post-hoc
benchmarks; they are not sampling constraints.

### Growth Models

Both models describe cumulative storage `S(t)` and annual storage rate
`dS/dt`.

Logistic:

```text
S(t) = C / (1 + k exp[-r(t - t0)]),    k = (C - S0) / S0
dS/dt = r S (1 - S / C)
```

Gompertz:

```text
S(t) = C exp[-b exp(-r tau)],          b = -ln(S0 / C), tau = t - t0
dS/dt = C r b exp(-r tau) exp[-b exp(-r tau)]
```

Where:

- `C` is total storage resource availability (Gt)
- `S0` is cumulative storage at base year 2030 (Gt)
- `t0` is 2030
- `r` is the model-internal growth parameter

### Shared CAGR Definition

Both Logistic and Gompertz use the same growth-rate metric:

```text
g_CAGR = (S(tn) / S0) ** (1 / (tn - t0)) - 1
```

`tn` is the inflection year of the annual storage-rate curve. `tp` is the peak
year of the annual storage-rate curve.

### Sampling and Acceptance

For unconstrained Type-A rows, the sampled pair is:

```text
g_CAGR ~ U(G_MIN, g_cap)
p2050  ~ U(P2050_MIN, p_hi)
```

For anchored Type-B rows (`us1gt`, `policy`, `ipcc_low`, `ipcc_high`), `p2050`
is fixed from the scenario CSV and only `g_CAGR` is sampled.

The active acceptance criteria are:

| # | Criterion |
|---|---|
| 1 | root-finding residual is below tolerance |
| 2 | `C > 0` and `C <= C_upper` |
| 3 | `2050 < tp < 2100` |
| 4 | `C` is finite |
| 5 | `rate_2100` is finite |

`tn` is reported as a diagnostic quantity, but it is not used as an acceptance
filter.

## Scenarios

All eight scenarios use `BASE_YEAR = 2030` and `END_YEAR = 2430`.

| ID | Scenario | Family | `C_upper` | `g_upper` | p2050 constraint |
|---|---|---|---:|---:|---|
| S1 | `reference` | Technical | `1 x C_ref` | 20%, China 25% | none |
| S2 | `minimum` | Technical | `0.1 x C_ref` | 10% | none |
| S3 | `maximum` | Technical | `10 x C_ref` | 20%, China 25% | none |
| S4 | `growth10` | Technical | `1 x C_ref` | 10% | none |
| S5 | `us1gt` | Policy | `1 x C_ref` | 20%, China 25% | US fixed at 1.0 Gt/yr |
| S6 | `policy` | Policy | `1 x C_ref` | 20%, China 25% | published targets |
| S7 | `ipcc_low` | Demand | `1 x C_ref` | 20%, China 25% | workbook low-demand point targets |
| S8 | `ipcc_high` | Demand | `1 x C_ref` | 20%, China 25% | workbook high-demand point targets |

### Reference Inputs

Anchor extraction convention: for each country/region sheet, use the row where
`Year == 2030`; rate is `Storage Rate (Mt/Year)`, cumulative storage is
`Cumulative Storage (Gt)`, and if that column is blank use
`Cumulative Storage (Mt) / 1000`.

| Country | `S0` 2030 cumulative (Gt) | 2030 rate (Mt/yr) | `C_ref` (Gt) |
|---|---:|---:|---:|
| UK | 0.113520 | 46.210000 | 78.000 |
| US | 1.355426 | 190.279749 | 506.000 |
| EU | 0.409187 | 109.514390 | 94.000 |
| China | 0.069866 | 9.062000 | 403.000 |
| Middle East | 0.099609 | 21.401800 | 45.600 |
| Australia | 0.061900 | 7.700000 | 502.400 |
| Canada | 0.109074 | 19.937000 | 404.000 |
| Indonesia | 0.017700 | 5.900000 | 15.900 |
| Thailand | 0.004000 | 1.000000 | 10.500 |
| Brazil | 0.326600 | 14.200000 | 300.000 |

### Anchored Targets

`us1gt` keeps only the US fixed:

| Country | p2050 target (Gt/yr) |
|---|---:|
| US | 1.000 |

`policy` uses published storage targets from workbook sheet `IPCCreference`:

| Country | p2050 target (Gt/yr) |
|---|---:|
| US | 1.000 |
| EU | 0.330 |
| UK | 0.175 |
| China | 0.216 |
| Canada | 0.060 |

`ipcc_low` / `ipcc_high` preserve both workbook demand columns in both CSVs.
The notebook reads `p2050_lo` for `ipcc_low` and `p2050_hi` for `ipcc_high`.

| Country | IPCC low (Gt/yr) | IPCC high (Gt/yr) |
|---|---:|---:|
| Australia | 0.00623012495040894 | 0.238633232709833 |
| Canada | 0.00178528571128845 | 0.443912423813518 |
| China | 0.0879427 | 6.57128321180474 |
| Indonesia | 0.0034121739470087 | 3.7238848 |
| Brazil | 0.0014094 | 1.95888264991766 |
| US | 0.00363542116016001 | 3.7238848 |

UK, EU, Middle East, and Thailand have no workbook IPCC low/high targets and
therefore remain Type-A unconstrained in both IPCC scenarios.

## Active Outputs

The active v8 output folder is:

- `output/v8_2026-06-01/`

Key files:

| File | Description |
|---|---|
| `acceptance.csv` | accepted sample counts by scenario, country, and model |
| `samples_{scenario}.csv` | accepted Monte Carlo samples for each scenario |
| `bands.csv` | per-country time-series bands |
| `scalars_central.csv` | selected central country-level scalar diagnostics |
| `timeseries_central.csv` | selected central country-level time series |
| `feasibility.csv` | empty-pool and unreachable-target diagnostics |
| `global_bands.csv` | global recombined pathway bands |
| `global_central.csv` | central global pathways |
| `sobol_indices.csv` | sensitivity summary |
| `anova_conditional.csv` | conditional ANOVA summary |
| `robustness_sample_paths.csv` | robustness pathway samples |
| `q25_q75_sample_paths.csv` | selected quartile sample paths |

Figures are under:

- `output/v8_2026-06-01/figures/`
- `output/v8_2026-06-01/figures/per_country_dossier/`
- `output/v8_2026-06-01/figures/scenario_dashboards/`

## Partial Pathways

A global pathway is complete when every country has at least one accepted
sample for the given `(scenario, model)` pair. If one or more countries have
zero accepted samples, the global pathway is flagged as partial and the missing
countries are listed in the output.

In v8, the notable zero-accepted combinations are:

| Scenario | Countries / models with zero accepted samples | Main reason |
|---|---|---|
| `ipcc_low` | US, Canada, Brazil; both models | low targets are below feasible model domain under the active filters |
| `ipcc_high` | China and Indonesia; both models | high targets exceed feasible domain at `C_ref` |
| `ipcc_high` | Brazil; Gompertz only | fixed high target cannot satisfy the active Gompertz constraints |

For included countries, country-level values are still valid accepted samples.
Only the global sum is partial.

## Running

To rebuild the four technical scenario inputs from the workbook:

```bash
python 01_growth_model/code/build_mcmc_inputs_from_workbook.py
```

To rebuild the four anchored scenario inputs:

```bash
python 01_growth_model/code/build_anchored_inputs_v8.py
```

To rerun the active notebook:

```bash
jupyter nbconvert --to notebook --execute \
  01_growth_model/notebooks/v8_growth_model_2026-06-01.ipynb \
  --output v8_growth_model_2026-06-01.ipynb
```

## Directory Structure

```text
01_growth_model/
  README.md
  input/
    2026_global_raw_data_20260601.xlsx
  code/
    build_mcmc_inputs_from_workbook.py
    build_anchored_inputs_v8.py
    preprocess_growth_model_inputs.py
    run_logistic_mcmc_from_scenarios.py
  notebooks/
    v8_growth_model_2026-06-01.ipynb
  output/
    mcmc_inputs/
    v8_2026-06-01/
```

## Archived Diagnostics

The historical checkpoint diagnostics are not required for the active v8
Part-1 deliverable because the v8 notebook already exports the acceptance,
feasibility, central, global, sensitivity, and figure outputs needed for the
paper workflow.

If a checkpoint-style QA layer is needed later, create a new v8-specific
diagnostic script from the active v8 outputs.

## Reference

> Zhang, Y., Ringrose, P., Krevor, S. et al.
> *Global assessment of geological CO2 storage resource deployment rates.*
> *Nature Communications* 15, 1071 (2024).
> [https://doi.org/10.1038/s41467-024-44726-0](https://doi.org/10.1038/s41467-024-44726-0)

This implementation keeps the two-model Logistic/Gompertz comparison,
country-level feasible-sample framework, and global recombination workflow,
while refreshing the country pool and workbook anchors for the final v8 paper
inputs.
