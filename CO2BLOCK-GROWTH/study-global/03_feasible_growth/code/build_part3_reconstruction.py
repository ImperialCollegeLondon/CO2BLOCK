"""Part 3 v8 smooth-reconstruction core (solve + envelope diagnostics).

Faithful port of the validated part3 notebook cells A–F, repointed to the v8
inputs. Reads the robust Step-1 input table (part3_input_table.csv), loads the
v8 ordering-specific screened envelopes, solves the 3-landmark reconstruction
(S0, R2050, R_peak_order) for the 78 failed cases, scores envelope consistency,
and writes smooth_reconstruction_summary.csv + smooth_reconstruction_timeseries.csv.

Math is identical to Part 1 (CAGR g ↔ ODE r) and to the prior notebook; only the
data source (v8) and the input-table sourcing (no scalars-only join) changed.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

T_BASE, T_END = 2030, 2179
YEARS = np.arange(T_BASE, T_END + 1)
PHI = (3.0 + math.sqrt(5.0)) / 2.0
SQRT3 = math.sqrt(3.0)
EXCEED_EPS_GT = 1e-4
KEY = ["country", "scenario", "model", "ordering"]

BASE = Path(__file__).resolve().parents[2]
P2_DIR = BASE / "02_co2block_screening" / "output" / "step2_screening"
OUT_DIR = BASE / "03_feasible_growth" / "output"
INPUT_TABLE = OUT_DIR / "part3_input_table.csv"

# curve helpers (Part-1 identical)
def logistic_inv(C, S0, t0, g):
    k = (C - S0) / S0; S_Tn = C / (3.0 + SQRT3); A = math.log(k) - math.log(2.0 + SQRT3)
    r = A * math.log(1.0 + g) / math.log(S_Tn / S0); return r, t0 + A / r, t0 + math.log(k) / r

def gompertz_inv(C, S0, t0, g):
    b = -math.log(S0 / C); P_Tn = C * math.exp(-PHI); B = math.log(b / PHI)
    r = B * math.log(1.0 + g) / math.log(P_Tn / S0); return r, b, t0 + B / r, t0 + math.log(b) / r

def logistic_rate_at(y, C, r, S0, t0=T_BASE):
    k = (C - S0) / S0; S = C / (1.0 + k * np.exp(-r * (y - t0))); return r * S * (1.0 - S / C)

def gompertz_rate_at(y, C, b, r, t0=T_BASE):
    s = b * np.exp(-r * (y - t0)); return C * r * s * np.exp(-s)

def logistic_peak_year(C, r, S0, t0=T_BASE):
    return t0 + math.log((C - S0) / S0) / r

def gompertz_peak_year(C, S0, r, t0=T_BASE):
    return t0 + math.log(-math.log(S0 / C)) / r

def forward_g(model, C, S0, g, years=YEARS):
    if model == "Logistic":
        r, tn, tp = logistic_inv(C, S0, T_BASE, g); rate = logistic_rate_at(years, C, r, S0)
    else:
        r, b, tn, tp = gompertz_inv(C, S0, T_BASE, g); rate = gompertz_rate_at(years, C, b, r)
    cum = S0 + np.concatenate([[0], np.cumsum(rate[:-1])]); return {"r": r, "tp": tp, "rate": rate, "cum": cum}

def forward_r(model, C, S0, r, years=YEARS):
    if model == "Logistic":
        rate = logistic_rate_at(years, C, r, S0); tp = logistic_peak_year(C, r, S0)
    else:
        b = -math.log(S0 / C); rate = gompertz_rate_at(years, C, b, r); tp = gompertz_peak_year(C, S0, r)
    cum = S0 + np.concatenate([[0], np.cumsum(rate[:-1])]); return {"rate": rate, "cum": cum, "tp": tp}

def r_to_g(model, C, S0, r):
    if model == "Logistic":
        k = (C - S0) / S0; S_Tn = C / (3.0 + SQRT3); A = math.log(k) - math.log(2.0 + SQRT3)
        if A <= 0 or math.log(S_Tn / S0) <= 0: return float("nan")
        return math.exp(r * math.log(S_Tn / S0) / A) - 1.0
    b = -math.log(S0 / C); P_Tn = C * math.exp(-PHI); B = math.log(b / PHI)
    if B <= 0 or math.log(P_Tn / S0) <= 0: return float("nan")
    return math.exp(r * math.log(P_Tn / S0) / B) - 1.0

# 3-landmark solver (peak_rate = r*C/4 logistic, C*r/e gompertz, gives C(r))
def _C_of_r(model, R_peak, r):
    return (4.0 * R_peak / r) if model == "Logistic" else (math.e * R_peak / r)

def _peak_year_of(model, C, S0, r):
    return logistic_peak_year(C, r, S0) if model == "Logistic" else gompertz_peak_year(C, S0, r)

def _rate_at_2050(model, C, S0, r):
    if model == "Logistic":
        return float(logistic_rate_at(2050.0, C, r, S0))
    return float(gompertz_rate_at(2050.0, C, -math.log(S0 / C), r))

def _bound_constraint(model, S0, R_peak):
    ratio = (3.0 + SQRT3) if model == "Logistic" else math.exp(PHI)
    coef = 4.0 if model == "Logistic" else math.e
    return coef * R_peak / (S0 * ratio * 1.01)

_NAN_COLS = ("C_recon_order_gt", "r_recon_order", "g_recon_order", "rate_2050_recon_order",
             "peak_rate_recon_order", "peak_year_recon_order", "tp_recon_order", "tn_recon_order")

def solve_reconstruction(row, env):
    model, S0 = row["model"], row["S_0"]
    R2050, R_peak = float(row["rate_2050_original"]), float(row["R_peak_order"])
    if R_peak < R2050 * (1.0 - 1e-6):
        return {"landmark_infeasible": True, "landmark_reason": "R_peak < R2050", **{c: np.nan for c in _NAN_COLS}}
    r_upper = _bound_constraint(model, S0, R_peak); r_grid = np.geomspace(1e-5, r_upper * 0.999, 600)
    f = []
    for r in r_grid:
        C = _C_of_r(model, R_peak, r)
        if C <= S0: f.append(np.nan); continue
        try: f.append(_rate_at_2050(model, C, S0, r) - R2050)
        except (ValueError, OverflowError): f.append(np.nan)
    f = np.asarray(f); valid = ~np.isnan(f)
    if not valid.any():
        return {"landmark_infeasible": True, "landmark_reason": "no valid r grid", **{c: np.nan for c in _NAN_COLS}}
    roots = []
    for i in range(len(r_grid) - 1):
        if valid[i] and valid[i + 1] and f[i] * f[i + 1] < 0:
            try:
                roots.append(brentq(lambda r: _rate_at_2050(model, _C_of_r(model, R_peak, r), S0, r) - R2050,
                                    r_grid[i], r_grid[i + 1], rtol=1e-7))
            except (ValueError, RuntimeError): pass
    if not roots:
        return {"landmark_infeasible": True, "landmark_reason": "no root in scan", **{c: np.nan for c in _NAN_COLS}}
    post = [(rr, _C_of_r(model, R_peak, rr), _peak_year_of(model, _C_of_r(model, R_peak, rr), S0, rr)) for rr in roots]
    post = [(rr, C, tp) for rr, C, tp in post if tp >= 2050.0]
    if not post:
        return {"landmark_infeasible": True, "landmark_reason": "all roots have peak < 2050", **{c: np.nan for c in _NAN_COLS}}
    if len(post) > 1:
        post.sort(key=lambda t: float(np.sqrt(np.mean((forward_r(model, t[1], S0, t[0])["rate"] - env) ** 2))))
    r_sel, C_sel, tp_sel = post[0]; fwd = forward_r(model, C_sel, S0, r_sel)
    if model == "Logistic":
        A = math.log((C_sel - S0) / S0) - math.log(2.0 + SQRT3); tn_sel = T_BASE + A / r_sel
    else:
        b = -math.log(S0 / C_sel); tn_sel = T_BASE + math.log(b / PHI) / r_sel
    return {"landmark_infeasible": False, "landmark_reason": "",
            "C_recon_order_gt": C_sel, "r_recon_order": r_sel, "g_recon_order": r_to_g(model, C_sel, S0, r_sel),
            "rate_2050_recon_order": _rate_at_2050(model, C_sel, S0, r_sel), "peak_rate_recon_order": R_peak,
            "peak_year_recon_order": int(YEARS[int(np.argmax(fwd["rate"]))]), "tp_recon_order": tp_sel, "tn_recon_order": tn_sel}


def main():
    it = pd.read_csv(INPUT_TABLE)
    failed = it[it.is_reconstruction_target].copy().rename(
        columns={"S0_gt": "S_0", "C_scenario_gt": "C_scenario"})
    failed["rate_2050_original"] = failed["rate_2050_original_mt"] / 1000.0  # -> Gt/yr
    assert len(failed) == 78, f"expected 78 targets, got {len(failed)}"

    # envelopes (Gt/yr) from v8 screened paths
    env_gt = {}
    for o in ["ascending", "descending"]:
        sp = pd.read_csv(P2_DIR / o / "final" / "case_model_screened_paths.csv")
        sp["ordering"] = o
        idx = sp.set_index(["country", "scenario", "model", "ordering", "year"])["screened_rate_mt_yr"].sort_index()
        for _, row in failed[failed.ordering == o].iterrows():
            ck = tuple(row[k] for k in KEY)
            arr = idx.loc[ck].reindex(YEARS).interpolate(limit_direction="both").to_numpy() / 1000.0
            env_gt[ck] = arr
    failed["R_peak_order"] = [float(np.max(env_gt[tuple(r[k] for k in KEY)])) for _, r in failed.iterrows()]
    failed["Y_peak_order"] = [int(YEARS[int(np.argmax(env_gt[tuple(r[k] for k in KEY)]))]) for _, r in failed.iterrows()]

    # solve + envelope diagnostics
    rec, diag = [], []
    for _, row in failed.iterrows():
        ck = tuple(row[k] for k in KEY); env = env_gt[ck]
        s = solve_reconstruction(row, env); rec.append({**{k: row[k] for k in KEY}, **s})
        if s["landmark_infeasible"]:
            diag.append({**{k: row[k] for k in KEY}, "max_envelope_exceedance_mt_yr": np.nan,
                         "years_exceeding_envelope": np.nan, "first_exceed_year": np.nan,
                         "last_exceed_year": np.nan, "rmse_to_envelope_mt_yr": np.nan, "envelope_consistent": False})
            continue
        fwd = forward_r(row["model"], s["C_recon_order_gt"], row["S_0"], s["r_recon_order"])
        d = (fwd["rate"] - env) * 1000.0; m = (fwd["rate"] - env) > EXCEED_EPS_GT; n = int(m.sum())
        diag.append({**{k: row[k] for k in KEY},
                     "max_envelope_exceedance_mt_yr": float(np.max(d)), "years_exceeding_envelope": n,
                     "first_exceed_year": int(YEARS[m.argmax()]) if n else None,
                     "last_exceed_year": int(YEARS[len(YEARS) - 1 - m[::-1].argmax()]) if n else None,
                     "rmse_to_envelope_mt_yr": float(np.sqrt(np.mean(d ** 2))), "envelope_consistent": bool(n == 0)})
    failed = failed.merge(pd.DataFrame(rec), on=KEY, how="left").merge(pd.DataFrame(diag), on=KEY, how="left")

    failed["C_ratio_recon_to_scenario"] = failed["C_recon_order_gt"] / failed["C_scenario"]
    failed["delta_tp_recon_vs_original"] = failed["tp_recon_order"] - failed["tp_original"]
    failed["delta_peak_year_recon_vs_screened"] = failed["peak_year_recon_order"] - failed["Y_peak_order"]
    failed["violates_g_cap"] = failed.apply(
        lambda r: bool(r["scenario"] == "growth10" and pd.notna(r["g_recon_order"]) and r["g_recon_order"] > r["g_cap"]), axis=1)

    # summary
    cols = ["country", "scenario", "model", "ordering", "Type", "S_0", "C_scenario", "g_cap",
            "g_original", "r_original", "rate_2050_original", "tn_original", "tp_original",
            "R_peak_order", "Y_peak_order", "C_recon_order_gt", "g_recon_order", "r_recon_order",
            "rate_2050_recon_order", "peak_rate_recon_order", "peak_year_recon_order", "tp_recon_order",
            "tn_recon_order", "C_ratio_recon_to_scenario", "delta_tp_recon_vs_original",
            "delta_peak_year_recon_vs_screened", "max_envelope_exceedance_mt_yr", "years_exceeding_envelope",
            "first_exceed_year", "last_exceed_year", "rmse_to_envelope_mt_yr", "envelope_consistent",
            "landmark_infeasible", "landmark_reason", "violates_g_cap"]
    out = failed[cols].rename(columns={"S_0": "S0_gt", "C_scenario": "C_scenario_gt"})
    out.to_csv(OUT_DIR / "smooth_reconstruction_summary.csv", index=False)

    # timeseries (original / screened_envelope / reconstructed), Mt/yr
    ts = []
    for _, row in failed.iterrows():
        ck = tuple(row[k] for k in KEY); base = {k: row[k] for k in KEY}; env = env_gt[ck]
        fo = forward_g(row["model"], row["C_scenario"], row["S_0"], row["g_original"])
        for y, rt, cm in zip(YEARS, fo["rate"], fo["cum"]):
            ts.append({**base, "year": int(y), "curve_kind": "original", "rate_mt_yr": float(rt) * 1000.0, "cum_gt": float(cm)})
        screened_cum = row["S_0"] + np.concatenate([[0.0], np.cumsum(env[:-1])])
        for y, rt, cm in zip(YEARS, env, screened_cum):
            ts.append({**base, "year": int(y), "curve_kind": "screened_envelope", "rate_mt_yr": float(rt) * 1000.0, "cum_gt": float(cm)})
        if not row["landmark_infeasible"]:
            fr = forward_r(row["model"], row["C_recon_order_gt"], row["S_0"], row["r_recon_order"])
            for y, rt, cm in zip(YEARS, fr["rate"], fr["cum"]):
                ts.append({**base, "year": int(y), "curve_kind": "reconstructed", "rate_mt_yr": float(rt) * 1000.0, "cum_gt": float(cm)})
    pd.DataFrame(ts).to_csv(OUT_DIR / "smooth_reconstruction_timeseries.csv", index=False)

    solved = failed[~failed.landmark_infeasible]
    print(f"reconstructed {len(solved)}/{len(failed)} (landmark_infeasible: {int(failed.landmark_infeasible.sum())})")
    print(f"  median C_recon/C_scenario : {solved.C_ratio_recon_to_scenario.median():.3f}")
    print(f"  envelope_consistent       : {int(solved.envelope_consistent.sum())}/{len(solved)}")
    print(f"  median years_exceeding    : {solved.years_exceeding_envelope.median():.0f}")
    print(f"  g_cap violations          : {int(failed.violates_g_cap.sum())}/{int((failed.scenario=='growth10').sum())} growth10")
    if int(failed.landmark_infeasible.sum()):
        print("  infeasible reasons:", dict(failed[failed.landmark_infeasible].landmark_reason.value_counts()))
    print(f"wrote smooth_reconstruction_summary.csv ({len(out)}), smooth_reconstruction_timeseries.csv")


if __name__ == "__main__":
    main()
