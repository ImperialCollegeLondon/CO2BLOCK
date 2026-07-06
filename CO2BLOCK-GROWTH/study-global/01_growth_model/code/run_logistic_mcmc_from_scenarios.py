import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[2]
INPUT_DIR = BASE_DIR / "01_growth_model" / "output" / "mcmc_inputs"
OUTPUT_DIR = BASE_DIR / "01_growth_model" / "output" / "mcmc_output"

SCENARIO_FILES = {
    "reference": INPUT_DIR / "mcmc_input_reference.csv",
    "minimum": INPUT_DIR / "mcmc_input_minimum.csv",
    "maximum": INPUT_DIR / "mcmc_input_maximum.csv",
    "growth10": INPUT_DIR / "mcmc_input_growth10.csv",
}

OUTPUT_FILES = {
    "reference": OUTPUT_DIR / "reference_results.csv",
    "minimum": OUTPUT_DIR / "minimum_results.csv",
    "maximum": OUTPUT_DIR / "maximum_results.csv",
    "growth10": OUTPUT_DIR / "growth10_results.csv",
}

TARGET_ACCEPTED_SAMPLES = 1000
MAX_ATTEMPTS_PER_COUNTRY = 100000
RATE_TOLERANCE = 0.01
MIN_GROWTH_RATE = 1e-4
COARSE_GRID_SIZE = 250
REFINE_GRID_SIZE = 120
REFINE_STEPS = 3


def stable_rng(*parts):
    payload = "::".join(str(part) for part in parts).encode("utf-8")
    seed = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big") % (2**32)
    return np.random.default_rng(seed)


def parse_numeric(x):
    if x is None:
        return np.nan
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x)
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "null"}:
        return np.nan
    s = s.replace(",", "")
    if s.endswith("%"):
        try:
            return float(s[:-1].strip()) / 100.0
        except Exception:
            return np.nan
    try:
        return float(s)
    except Exception:
        return np.nan


def parse_growth_value(x, default=np.nan):
    v = parse_numeric(x)
    if not np.isfinite(v):
        return default
    return v / 100.0 if v > 1.0 else v


def calculate_cagr(p_t0, p_t1, t0, t1):
    if p_t0 <= 0 or p_t1 <= 0 or t1 <= t0:
        return np.nan
    return (p_t1 / p_t0) ** (1.0 / (t1 - t0)) - 1.0


def logistic_cumulative(t, c, r, t0, s0):
    t = np.asarray(t, dtype=float)
    k = (c - s0) / s0
    return c / (1.0 + k * np.exp(-r * (t - t0)))


def logistic_rate(t, c, r, t0, s0):
    s = logistic_cumulative(t, c, r, t0, s0)
    return r * s * (1.0 - s / c)


def logistic_inverse_r_from_growth_to_tn(c, s0, t0, g):
    if not (0 < s0 < c):
        raise ValueError("Logistic: require 0 < S0 < C.")
    s_tn = c / (3.0 + math.sqrt(3.0))
    if s0 >= s_tn:
        raise ValueError(f"Logistic: need S0 < C/(3+sqrt3)≈{s_tn:.6f} so Tn is after t0.")
    a = math.log((c - s0) / s0) - math.log(2.0 + math.sqrt(3.0))
    r = a * math.log(1.0 + g) / math.log(s_tn / s0)
    tn = t0 + a / r
    k = (c - s0) / s0
    tp = t0 + math.log(k) / r
    return float(r), float(tn), float(tp), float(s_tn)


def best_capacity_for_growth(s0_gt, p0_gt_yr, c_upper, g, base_year):
    lower = max(s0_gt * 1.0001, 1e-10)
    if lower >= c_upper:
        return None

    left = lower
    right = c_upper
    best = None

    for step in range(REFINE_STEPS + 1):
        grid_size = COARSE_GRID_SIZE if step == 0 else REFINE_GRID_SIZE
        capacities = np.geomspace(left, right, grid_size)
        local_best = None
        local_idx = None

        for idx, c_gt in enumerate(capacities):
            try:
                r, tn, tp, s_tn = logistic_inverse_r_from_growth_to_tn(c_gt, s0_gt, base_year, g)
            except Exception:
                continue
            if not np.isfinite(tp) or tp <= 2050.0 or tp >= 2100.0:
                continue

            p0_model = float(logistic_rate(base_year, c_gt, r, base_year, s0_gt))
            rate_error = abs(p0_model - p0_gt_yr) / p0_gt_yr if p0_gt_yr > 0 else np.inf
            row = {
                "c_gt": float(c_gt),
                "tn": float(tn),
                "tp": float(tp),
                "r": float(r),
                "s_tn": float(s_tn),
                "rate_error": float(rate_error),
                "rate_2050": float(logistic_rate(2050, c_gt, r, base_year, s0_gt)),
            }
            if local_best is None or row["rate_error"] < local_best["rate_error"]:
                local_best = row
                local_idx = idx

        if local_best is None:
            return None

        best = local_best
        if step == REFINE_STEPS:
            break

        lo_idx = max(0, local_idx - 1)
        hi_idx = min(len(capacities) - 1, local_idx + 1)
        left = float(capacities[lo_idx])
        right = float(capacities[hi_idx])
        if left == right:
            break

    if best is None or best["rate_error"] > RATE_TOLERANCE:
        return None
    return best


def run_country_mcmc(country, c_upper, s0_gt, p0_gt_yr, g_cap, base_year):
    rng = stable_rng(country, c_upper, s0_gt, p0_gt_yr, g_cap, base_year)
    accepted = []
    attempts = 0

    while len(accepted) < TARGET_ACCEPTED_SAMPLES and attempts < MAX_ATTEMPTS_PER_COUNTRY:
        attempts += 1
        g = float(rng.uniform(MIN_GROWTH_RATE, g_cap))
        best = best_capacity_for_growth(
            s0_gt=s0_gt,
            p0_gt_yr=p0_gt_yr,
            c_upper=c_upper,
            g=g,
            base_year=base_year,
        )
        if best is None:
            continue

        cagr = calculate_cagr(s0_gt, best["s_tn"], base_year, best["tn"])
        if not np.isfinite(cagr):
            continue

        accepted.append(
            {
                "Country": country,
                "Growth Model": "Logistic",
                "Growth Rate [%]": round(float(100.0 * cagr), 6),
                "Storage Resource Base Required [Gt]": round(best["c_gt"], 6),
                "Modelled Storage Rate 2050 [Gt/yr]": round(best["rate_2050"], 9),
            }
        )

    return accepted


def build_results(input_csv):
    df = pd.read_csv(input_csv)
    rows = []
    for _, row in df.iterrows():
        rows.extend(
            run_country_mcmc(
                country=str(row["Country"]).strip(),
                c_upper=float(row["Storage capacity (Gt)"]),
                s0_gt=float(row["2030 Cumulative Storage (Gt)"]),
                p0_gt_yr=float(row["2030 Storage Rate (Mt/Year)"]) / 1000.0,
                g_cap=parse_growth_value(row["Growth_to_Tn_%"]),
                base_year=int(row["BASE_YEAR"]),
            )
        )

    out = pd.DataFrame(
        rows,
        columns=[
            "Country",
            "Growth Model",
            "Growth Rate [%]",
            "Storage Resource Base Required [Gt]",
            "Modelled Storage Rate 2050 [Gt/yr]",
        ],
    )
    return out.sort_values(["Country", "Growth Rate [%]"], kind="stable").reset_index(drop=True)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for scenario, input_csv in SCENARIO_FILES.items():
        result = build_results(input_csv)
        result.to_csv(OUTPUT_FILES[scenario], index=False)


if __name__ == "__main__":
    main()
