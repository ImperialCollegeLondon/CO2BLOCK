# CO₂ Storage Growth Modelling and Basin Screening

This repository develops a two-stage research workflow for evaluating whether rapid geological CO₂ storage deployment is both growth-feasible at national scale and deliverable at basin scale.

If global-scale CO₂ storage expands to climate-relevant levels, which countries are expected to carry that growth, and can their basin portfolios physically sustain those pathways? That is the question the workflow is built to answer.

The workflow is organised in two linked modules:

- [01_growth_model](./01_growth_model): builds national CO₂ storage pathways using Logistic and Gompertz growth models under technical, policy, and demand-constrained scenarios.
- [02_co2block_screening](./02_co2block_screening): screens whether those national pathways remain physically deliverable at basin scale using CO₂BLOCK-style analytical constraints.

The project is informed by Zhang et al. (2024), but it is not intended as a line-by-line reproduction. Instead, it extends that lineage with a stricter feasible-sample framework, a harmonised CAGR definition across growth models, IQR-based primary pathway reporting, and a direct interface to basin screening.

## What This Repository Does

The repository moves in one direction:

1. Generate feasible country-level growth samples. For each country, scenario, and growth model, the code samples `CAGR` and `p2050`, solves for the implied storage resource requirement, and retains only trajectories that satisfy explicit feasibility conditions.

2. Construct global pathways. Accepted country-level sample pools are recombined into global 2050 storage pathways. For each scenario and model, the globally minimum and maximum 2050 pathways are selected, and the corresponding country values are read back.

3. Screen delivery at basin scale. Those country-level trajectories are then passed to basin-screening workflows to test whether national basin portfolios can plausibly sustain the required injection rates over time.

This makes the repository useful for two purposes: interpreting what global low/high deployment futures imply for individual countries, and checking whether those country trajectories remain plausible once geological delivery constraints are applied.

## Study Domain

The current growth model covers ten country or regional aggregates:

- Australia
- Brazil
- Canada
- China
- EU (reported label for EU without UK)
- Indonesia
- Middle East
- Thailand
- UK
- US

## Repository Map

```text
.
  01_growth_model/
    input/                  Workbook inputs and scenario definitions
    code/                   Preprocessing utilities
    notebooks/              Main growth-model notebooks
    output/
      mcmc_inputs/          Scenario-specific input CSVs
      v8_2026-06-01/        Accepted samples, pathway summaries, plots
  02_co2block_screening/
    code/                   Basin screening code and CO₂BLOCK Python port
    input/                  Basin datasets and imported growth paths
    output/                 Screening results and figures
    docs/                   Data dictionaries and notes
  90_reference_matlab/      Reference MATLAB sources and validation notes
  99_shared/environment/    Shared Python requirements
  README.md
```

## Module 1: Growth Modelling

The growth module treats `p2050` as the anchor variable of interest and asks what storage resource base and long-run pathway are implied if a country reaches that rate under Logistic or Gompertz growth.

### Growth equations

Both models describe cumulative storage and the corresponding annual storage rate.

For the Logistic model, cumulative storage is:

$$
S(t) = \frac{C}{1 + k \exp[-r(t-t_0)]},
\qquad
k = \frac{C-S_0}{S_0}
$$

and the annual storage rate is:

$$
\frac{dS}{dt} = rS\left(1-\frac{S}{C}\right)
$$

For the Gompertz model, letting $\tau = t-t_0$, cumulative storage is:

$$
S(\tau) = C \exp[-b \exp(-r\tau)]
$$

and the annual storage rate is:

$$
\frac{dS}{d\tau} = C r b \exp(-r\tau)\exp[-b\exp(-r\tau)]
$$

where:

- $C$ is total storage resource availability
- $S_0$ is cumulative storage at the base year
- $t_0$ is the base year (`2030` in the present framework)
- $r$ is the model-internal growth parameter
- $b = -\ln(S_0/C)$ for the Gompertz model

### Inflection year and peak year of the annual rate curve

The workflow distinguishes two different time markers from the annual storage-rate curve:

- $t_n$: inflection year of the annual storage-rate curve
- $t_p$: peak year of the annual storage-rate curve

For the Logistic model:

$$
S(T_n) = \frac{C}{3+\sqrt{3}}
$$

$$
t_n = t_0 + \frac{A}{r},
\qquad
A = \ln\left(\frac{C-S_0}{S_0}\right) - \ln(2+\sqrt{3})
$$

$$
t_p = t_0 + \frac{1}{r}\ln\left(\frac{C-S_0}{S_0}\right)
$$

For the Gompertz model:

$$
b = -\ln\left(\frac{S_0}{C}\right),
\qquad
\phi = \frac{3+\sqrt{5}}{2}
$$

$$
S(T_n) = C \exp(-\phi)
$$

$$
t_n = t_0 + \frac{B}{r},
\qquad
B = \ln\left(\frac{b}{\phi}\right)
$$

$$
t_p = t_0 + \frac{1}{r}\ln(b)
$$

In this framework, accepted Monte Carlo samples must satisfy:

- $2050 < t_p < 2100$

The inflection year $t_n$ is retained as a diagnostic quantity and for the
shared CAGR definition; it is not used as an additional acceptance filter in
the active V8 workflow.

### Shared CAGR definition

To make Logistic and Gompertz directly comparable, growth is not reported using the model-internal parameter $r$. Instead, both models are expressed through the same CAGR metric:

$$
g_{\mathrm{CAGR}} = \left(\frac{S(T_n)}{S_0}\right)^{1/(t_n-t_0)} - 1
$$

This is the compound annual growth rate of cumulative storage from the base year to the inflection year of the annual storage-rate curve.

### Core modelling choices

- Growth rate is defined consistently across both models as CAGR from 2030 to the inflection year of the annual storage-rate curve.
- Samples are accepted only if they satisfy explicit feasibility filters, including:
  - `2050 < tp < 2100`
  - positive and bounded storage resource requirement
  - small numerical residual in the inverse solve
- Country-level accepted samples are then used to build global-first pathway interpretations.

### Scenario families

The current growth stage includes:

- Technical feasibility scenarios
  - `reference`
  - `minimum`
  - `maximum`
  - `growth10`
- Policy-anchored scenarios
  - `us1gt`
  - `policy`
- Demand-constrained scenarios
  - `ipcc_low`
  - `ipcc_high`

These scenarios are implemented through combinations of:

- upper bounds on storage resource availability
- upper bounds on growth rate
- optional country-specific constraints on `p2050`

### Outputs

The active growth outputs live in [01_growth_model/output/v8_2026-06-01](./01_growth_model/output/v8_2026-06-01).

The primary files currently retained at top level include accepted samples, central pathways, uncertainty bands, feasibility summaries, and global aggregation outputs:

- [samples_reference.csv](./01_growth_model/output/v8_2026-06-01/samples_reference.csv)
- [samples_minimum.csv](./01_growth_model/output/v8_2026-06-01/samples_minimum.csv)
- [samples_maximum.csv](./01_growth_model/output/v8_2026-06-01/samples_maximum.csv)
- [samples_growth10.csv](./01_growth_model/output/v8_2026-06-01/samples_growth10.csv)
- [samples_us1gt.csv](./01_growth_model/output/v8_2026-06-01/samples_us1gt.csv)
- [samples_policy.csv](./01_growth_model/output/v8_2026-06-01/samples_policy.csv)
- [samples_ipcc_low.csv](./01_growth_model/output/v8_2026-06-01/samples_ipcc_low.csv)
- [samples_ipcc_high.csv](./01_growth_model/output/v8_2026-06-01/samples_ipcc_high.csv)
- [timeseries_central.csv](./01_growth_model/output/v8_2026-06-01/timeseries_central.csv)
- [global_central.csv](./01_growth_model/output/v8_2026-06-01/global_central.csv)

## Module 2: Basin Screening

The basin-screening stage asks a harder question: even if a national growth pathway is mathematically feasible, can it actually be delivered by the country's available basin portfolio?

This stage uses CO₂BLOCK-style analytical screening to test:

- injection-rate feasibility
- basin allocation structure
- pressure and storage limits over time

The intention is not to treat national growth curves as an endpoint, but as inputs to a downstream deliverability test.

The active screening outputs live in [02_co2block_screening/output/step2_screening](./02_co2block_screening/output/step2_screening). The V8 screening package evaluates the ten study regions under eight scenarios, two growth models, and two basin-ordering rules. Of the 150 country-scenario-model cases with valid Part-1 central pathways, 111 pass the whole-period basin-screening test and 39 fail under both ordering rules; there are no pass/fail flips between ascending and descending basin allocation.

## How The Two Modules Connect

The logic of the repository runs from global pathway construction, to conditional country trajectories, to basin-level deliverability screening.

That ordering matters.

The growth module identifies what countries would contribute on selected global pathways. The screening module then evaluates whether those country trajectories can be supported by basin-scale geology. This separation keeps the growth stage focused on pathway construction and the basin stage focused on physical delivery constraints.

## Getting Started

### Requirements

- Python 3.10+
- `numpy`
- `pandas`
- `scipy`
- `matplotlib`
- `openpyxl`

Install dependencies with:

```bash
pip install -r 99_shared/environment/requirements.txt
```

### Run the growth notebook

```bash
jupyter nbconvert --to notebook --execute \
  01_growth_model/notebooks/v8_growth_model_2026-06-01.ipynb \
  --output v8_growth_model_2026-06-01.executed.ipynb
```

Outputs will be written to:

- [01_growth_model/output/v8_2026-06-01](./01_growth_model/output/v8_2026-06-01)

## Reference Point

This repository is built in conversation with:

> Zhang, Y., Ringrose, P., Krevor, S. et al.  
> *Global assessment of geological CO₂ storage resource deployment rates.*  
> *Nature Communications* 15, 1071 (2024).  
> [https://doi.org/10.1038/s41467-024-44726-0](https://doi.org/10.1038/s41467-024-44726-0)

The present implementation extends that line of work toward a more explicit feasibility-screening framework and a direct connection to basin-level delivery analysis.

## Status

This is an active research repository. The current V8 paper-ready package is the two-module workflow above: Part 1 growth modelling and Part 2 basin screening. Downstream exploratory workspaces, including `03_feasible_growth` and `04_global_aggregate_2050`, should be treated as design notes until rerun from the current V8 outputs.
