"""Step 1 V8 incremental precompute helper.

This helper extends an existing V8-compatible basin-period matrix by computing
the active regional additions: EU (Europe without UK), Middle East, and Brazil.
It keeps a cache-trust gate for reused rows and writes a complete matrix to
temporary files for manual installation.
"""
from __future__ import annotations

import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

CODE_ROOT = Path(__file__).resolve().parent
SCREENING_ROOT = CODE_ROOT.parent
sys.path.insert(0, str(CODE_ROOT))
from co2block_py.core import CalculationConfig, build_regional_summary, calculate_site
from co2block_py.screening import assign_screening_case

BASIN_FILE = SCREENING_ROOT / "input" / "basin_data" / "Global.xlsx"
LIVE_LONG = SCREENING_ROOT / "output" / "step1_precompute" / "precompute" / "basin_period_resource_long.csv"
OUT_LONG = Path("/tmp/step1_v8_long.csv")
OUT_META = Path("/tmp/step1_v8_metadata.csv")

INJ_DURATION_YR = 150
STORAGE_PERIOD_STEP_YR = 10
NR_DIST = 100
MAX_Q_MT_YR = 20.0
MIN_Q_MT_YR = 1.0
NEW_CASES = {"EU", "Middle East", "Brazil"}
PERIODS = list(range(STORAGE_PERIOD_STEP_YR, INJ_DURATION_YR + 1, STORAGE_PERIOD_STEP_YR))
WORKERS = 6


def safe_float(v):
    f = float(v)
    return f if math.isfinite(f) else 0.0


def evaluate_task(args):
    region_no, period = args
    calc = CalculationConfig(
        data_path=str(BASIN_FILE), site_no=region_no, correction="off",
        dist_min_km=2.0, dist_max_km="auto", nr_dist=NR_DIST, nr_well_max="auto",
        rw_m=0.2, time_yr=period, max_q_mt_per_year=MAX_Q_MT_YR, min_q_mt_per_year=MIN_Q_MT_YR,
    )
    result = calculate_site(calc)
    summary, _, _ = build_regional_summary(result, MIN_Q_MT_YR)
    return {
        "Region_no": region_no, "Region_name": summary.region, "Period_yr": period,
        "Region_Q_Mt_yr": safe_float(summary.region_q_mt_per_year),
        "Number_of_Wells": int(summary.number_of_wells),
        "Max_Capacity_Gt": safe_float(summary.max_capacity_gt),
        "Max_Q_Mt_yr": safe_float(summary.max_q_mt_per_year),
        "Optimum_Q_Mt_yr": safe_float(summary.optimum_q_mt_per_year),
        "Distance_km": safe_float(summary.distance_km),
    }


def main():
    basin_df = pd.read_excel(BASIN_FILE).copy()
    basin_df.insert(0, "region_no", range(1, len(basin_df) + 1))
    basin_df["AssignedCountry"] = basin_df["Majority Country"]
    basin_df["RegionGroup"] = basin_df["Region"]
    basin_df["screening_case"] = basin_df.apply(
        lambda r: assign_screening_case(r.get("Majority Country"), r.get("Region")), axis=1)

    existing = pd.read_csv(LIVE_LONG)
    existing_basins = set(existing["Region_no"].unique())
    print(f"Existing cached matrix: {len(existing)} rows, {len(existing_basins)} basins")

    new_basins = basin_df[basin_df["screening_case"].isin(NEW_CASES)]["region_no"].astype(int).tolist()
    print(f"Regional-addition basins (EU/Middle East/Brazil): {len(new_basins)}")
    for case in sorted(NEW_CASES):
        n = (basin_df["screening_case"] == case).sum()
        print(f"  {case:12s} {n:3d} basins")

    # no overlap with existing
    overlap = existing_basins & set(new_basins)
    assert not overlap, f"region_no overlap between existing and new: {overlap}"

    # Cache-trust gate: recompute 2 existing basins, must match cache exactly
    sample = sorted(existing_basins)[:2] + sorted(existing_basins)[-1:]
    print(f"\nCache-trust gate: recomputing existing basins {sample} ...")
    for rno in sample:
        for p in (10, 150):
            got = evaluate_task((int(rno), p))["Region_Q_Mt_yr"]
            ref = float(existing[(existing.Region_no == rno) & (existing.Period_yr == p)]["Region_Q_Mt_yr"].iloc[0])
            assert abs(got - ref) <= 1e-9 + 1e-9 * abs(ref), f"MISMATCH basin {rno} @ {p}yr: {got} vs cached {ref}"
    print("  gate PASSED: cached baseline rows are reproducible; reusing them.\n")

    # Parallel compute of new basins
    tasks = [(int(r), p) for r in new_basins for p in PERIODS]
    print(f"Computing {len(tasks)} tasks ({len(new_basins)} basins x {len(PERIODS)} periods) on {WORKERS} workers...")
    t0 = time.time()
    rows = []
    done = 0
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(evaluate_task, t): t for t in tasks}
        for fut in as_completed(futs):
            rows.append(fut.result())
            done += 1
            if done % 30 == 0 or done == len(tasks):
                el = time.time() - t0
                print(f"  {done}/{len(tasks)} done | {el/60:.1f} min | eta {el/done*(len(tasks)-done)/60:.1f} min", flush=True)

    new_long = pd.DataFrame(rows)
    meta = basin_df.set_index("region_no")[["AssignedCountry", "RegionGroup"]].to_dict("index")
    new_long["AssignedCountry"] = new_long["Region_no"].map(lambda x: meta.get(x, {}).get("AssignedCountry", "Unknown"))
    new_long["RegionGroup"] = new_long["Region_no"].map(lambda x: meta.get(x, {}).get("RegionGroup", "Unknown"))
    new_long = new_long[existing.columns.tolist()]

    full = pd.concat([existing, new_long], ignore_index=True).sort_values(
        ["Region_no", "Period_yr"]).reset_index(drop=True)
    full.to_csv(OUT_LONG, index=False)
    metadata = full[["Region_no", "Region_name", "AssignedCountry", "RegionGroup"]].drop_duplicates().sort_values("Region_no")
    metadata.to_csv(OUT_META, index=False)

    print(f"\nDONE in {(time.time()-t0)/60:.1f} min")
    print(f"  complete V8 matrix: {len(full)} rows, {full['Region_no'].nunique()} basins -> {OUT_LONG}")
    print(f"  metadata: {len(metadata)} basins -> {OUT_META}")
    # quick region resolution sanity
    from co2block_py.screening import resolve_region_nos
    for case in ["US", "China", "Indonesia", "Australia", "UK", "Canada", "Thailand", "EU", "Middle East", "Brazil"]:
        nos, why = resolve_region_nos(metadata, case)
        print(f"    {case:12s} -> {len(nos):2d} basins ({why})")


if __name__ == "__main__":
    main()
