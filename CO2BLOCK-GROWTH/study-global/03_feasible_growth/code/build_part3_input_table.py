"""Part 3 (feasible growth), Step 1: build the input table.

Assembles one robust table covering ALL 150 v8 screened central paths
(country, scenario, model) x 2 allocation orderings = 300 rows, each row
carrying the parameters Part 3 needs (original Part-1 landmarks + scenario
config + the Part-2 ordering-specific screened peak + fail status).

Source-of-truth rule (per the cross-check): never rely on scalars_central.csv
for coverage, since it omits 2 model-paths whose sibling-model central was
infeasible (ipcc_high/Brazil/Logistic, ipcc_low/Canada/Logistic). Instead:
  - per-(c,s,m) landmarks  <- timeseries_central.csv      (complete, 150 paths)
  - per-(c,s)  config      <- mcmc_input_<scenario>.csv    (C_scenario, g_cap, Type)
  - per-(c,s,m) analytic   <- scalars_central.csv          (g/r/tn/tp; may be NaN
                                                            for the 2 omitted paths,
                                                            both of which PASS so are
                                                            not reconstruction targets)
  - per-(c,s,m,order)      <- case_model_screened_paths.csv (R_peak_order, Y_peak_order)
  - fail status            <- case_model_shortfall_windows.csv

Reconstruction targets = rows with whole_period_pass == False (78 = 39 x 2).
Asserts zero param gaps among those targets.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
GM = BASE / "01_growth_model" / "output" / "v8_2026-06-01"
MI = BASE / "01_growth_model" / "output" / "mcmc_inputs"
S2 = BASE / "02_co2block_screening" / "output" / "step2_screening"
OUT = BASE / "03_feasible_growth" / "output" / "part3_input_table.csv"

SCENARIOS = ["reference", "minimum", "maximum", "growth10", "us1gt", "policy", "ipcc_low", "ipcc_high"]
ORDERINGS = ["ascending", "descending"]


def per_path_landmarks(ts: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (s, c, m), g in ts.groupby(["Scenario", "Country", "Model"]):
        g = g.sort_values("Year")
        rows.append(dict(
            country=c, scenario=s, model=m,
            S0_gt=float(g.loc[g.Year == 2030, "cumulative_central"].iloc[0]),
            rate_2050_original_mt=float(g.loc[g.Year == 2050, "rate_central"].iloc[0]) * 1000.0,
            peak_rate_original_mt=float(g["rate_central"].max()) * 1000.0,
            Y_peak_original=int(g.loc[g["rate_central"].idxmax(), "Year"]),
        ))
    return pd.DataFrame(rows)


def scenario_config() -> pd.DataFrame:
    rows = []
    for scen in SCENARIOS:
        d = pd.read_csv(MI / f"mcmc_input_{scen}.csv")
        for _, r in d.iterrows():
            typ = "B_fixed" if ("p2050_lo" in d.columns and pd.notna(r.get("p2050_lo"))) else "A"
            rows.append(dict(country=str(r["Country"]).strip(), scenario=scen,
                             C_scenario_gt=float(r["Storage capacity (Gt)"]),
                             g_cap=float(r["Growth_to_Tn_%"]) / 100.0, Type=typ))
    return pd.DataFrame(rows)


def analytic_params(sc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in sc.iterrows():
        for m, gk, rk, tnk, tpk in [
            ("Logistic", "g_star_L", "r_L", "tn_L", "tp_L"),
            ("Gompertz", "g_star_G", "r_G", "tn_G", "tp_G")]:
            rows.append(dict(country=r.Country, scenario=r.Scenario, model=m,
                             g_original=r[gk], r_original=r[rk],
                             tn_original=r[tnk], tp_original=r[tpk]))
    return pd.DataFrame(rows)


def main() -> None:
    ts = pd.read_csv(GM / "timeseries_central.csv")
    sc = pd.read_csv(GM / "scalars_central.csv")

    base = (per_path_landmarks(ts)
            .merge(scenario_config(), on=["country", "scenario"], how="left")
            .merge(analytic_params(sc), on=["country", "scenario", "model"], how="left"))
    assert len(base) == 150, f"expected 150 base paths, got {len(base)}"

    out = []
    for o in ORDERINGS:
        sp = pd.read_csv(S2 / o / "final" / "case_model_screened_paths.csv")
        sw = pd.read_csv(S2 / o / "final" / "case_model_shortfall_windows.csv")
        pk = (sp.groupby(["country", "scenario", "model"])
                .apply(lambda g: pd.Series({
                    "R_peak_order_mt": g.screened_rate_mt_yr.max(),
                    "Y_peak_order": int(g.loc[g.screened_rate_mt_yr.idxmax(), "year"])}),
                       include_groups=False)
                .reset_index())
        swk = sw[["country", "scenario", "model", "whole_period_pass",
                  "years_of_shortfall", "total_unmet_volume_gt"]]
        b = base.assign(ordering=o).merge(pk, on=["country", "scenario", "model"], how="left") \
                                    .merge(swk, on=["country", "scenario", "model"], how="left")
        out.append(b)
    full = pd.concat(out, ignore_index=True)
    full["is_reconstruction_target"] = ~full["whole_period_pass"].astype(bool)

    # --- coverage guarantees ---
    crit = ["S0_gt", "rate_2050_original_mt", "peak_rate_original_mt", "C_scenario_gt",
            "g_cap", "Type", "R_peak_order_mt", "whole_period_pass"]
    assert len(full) == 300, f"expected 300 rows, got {len(full)}"
    assert full[crit].isna().sum().sum() == 0, "null in critical columns"
    tgt = full[full.is_reconstruction_target]
    assert len(tgt) == 78, f"expected 78 reconstruction targets, got {len(tgt)}"
    # reconstruction needs the analytic landmarks too -> must be complete among targets
    assert tgt[["g_original", "r_original", "tp_original"]].isna().sum().sum() == 0, \
        "reconstruction target missing analytic params"

    cols = ["country", "scenario", "model", "ordering", "Type", "is_reconstruction_target",
            "S0_gt", "C_scenario_gt", "g_cap", "g_original", "r_original", "tn_original", "tp_original",
            "rate_2050_original_mt", "peak_rate_original_mt", "Y_peak_original",
            "R_peak_order_mt", "Y_peak_order",
            "whole_period_pass", "years_of_shortfall", "total_unmet_volume_gt"]
    full = full[cols].sort_values(["ordering", "country", "scenario", "model"]).reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    full.to_csv(OUT, index=False)
    print(f"wrote {OUT}  ({len(full)} rows; {len(tgt)} reconstruction targets)")
    print("reconstruction targets by country:",
          dict(tgt[tgt.ordering == "ascending"].country.value_counts()))


if __name__ == "__main__":
    main()
