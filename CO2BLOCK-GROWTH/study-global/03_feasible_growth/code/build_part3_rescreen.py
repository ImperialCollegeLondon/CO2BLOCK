"""Part 3 re-screening (README_NOTES section 6).

Runs the Part-2 CO2BLOCK allocation core on each reconstructed smooth curve
(under its own ordering, reusing Part-2's cached per-country resource matrices
and the MATLAB-faithful allocation), reconstructs the screened path with Part-2's
exact method, and records whether the smoothed curve now passes, how much is
unmet, and how the basin portfolio redistributes.

Outputs: reconstructed_rescreen_summary.csv + reconstructed_assignment.csv.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
S2_CODE = BASE / "02_co2block_screening" / "code"
sys.path.insert(0, str(S2_CODE))
from co2block_py.allocation import AllocationConfig, run_allocation_workflow
from co2block_py.screening import slugify

P2_DIR = BASE / "02_co2block_screening" / "output" / "step2_screening"
BASIN_FILE = BASE / "02_co2block_screening" / "input" / "basin_data" / "Global.xlsx"
OUT_DIR = BASE / "03_feasible_growth" / "output"
CHECKPOINTS = [2050, 2100, 2179]
ORDER_MAP = {"ascending": "ascend", "descending": "descend"}


def _active_mask(years, start, end):
    return (years >= start) & (years <= end)


def reconstruct_screened(years_w, raw_w, assignment_path):
    """Part-2 Stage-B method: assignment table -> per-year delivered (clipped to demand)."""
    rows = []
    if assignment_path.exists():
        a = pd.read_excel(assignment_path)
        for _, r in a.iterrows():
            for yr in years_w[_active_mask(years_w.astype(float), float(r["Start [y]"]), float(r["End [y]"]))]:
                rows.append({"year": int(yr), "basin": str(r["Region name"]), "region_no": int(r["Region no"]),
                             "allocated_rate_mt_yr": float(r["Step rate increment [Mt/y]"]),
                             "delivered_rate_mt_yr": float(r["Step rate increment [Mt/y]"]),
                             "wells_used": float(r["No_sites"]), "step": int(r["Step"]),
                             "duration": float(r["Duration [y]"])})
    ba = pd.DataFrame(rows)
    if not ba.empty:
        ba["delivered_rate_mt_yr"] = 0.0
        y2raw = dict(zip(years_w.tolist(), raw_w.tolist()))
        for yr in years_w:
            ym = ba["year"] == yr
            if not ym.any():
                continue
            remaining = float(y2raw[int(yr)])
            for idx in ba.loc[ym].sort_values(["step", "basin"]).index:
                d = min(float(ba.at[idx, "allocated_rate_mt_yr"]), max(remaining, 0.0))
                ba.at[idx, "delivered_rate_mt_yr"] = d
                remaining -= d
        annual = ba.groupby("year", as_index=False)["delivered_rate_mt_yr"].sum()
    else:
        annual = pd.DataFrame({"year": years_w, "delivered_rate_mt_yr": 0.0})
    annual = pd.DataFrame({"year": years_w}).merge(annual, on="year", how="left").fillna(0.0)
    sc = pd.DataFrame({"year": years_w, "raw_rate_mt_yr": raw_w,
                       "screened_rate_mt_yr": annual["delivered_rate_mt_yr"].values})
    sc["gap_rate_mt_yr"] = sc["raw_rate_mt_yr"] - sc["screened_rate_mt_yr"]
    return ba, sc


def cp_metrics(sc, ba, cp):
    s = sc[sc.year <= cp]
    if s.empty:
        return {}
    short = s[s.gap_rate_mt_yr > 1e-9]
    out = {"checkpoint": cp,
           "years_of_shortfall": int(len(short)),
           "total_unmet_volume_gt": float(s.gap_rate_mt_yr.clip(lower=0).sum() / 1000.0),
           "total_wells": 0, "active_basins": 0}
    if not ba.empty:
        at = ba[ba.year == cp]
        if at.empty:
            vy = ba[ba.year <= cp]["year"]
            if not vy.empty:
                at = ba[ba.year == vy.max()]
        if not at.empty:
            sh = at.groupby("basin")["delivered_rate_mt_yr"].sum()
            out["active_basins"] = int((sh > 0).sum())
            out["total_wells"] = int(at["wells_used"].sum())
    return out


def main():
    summ = pd.read_csv(OUT_DIR / "smooth_reconstruction_summary.csv")
    ts = pd.read_csv(OUT_DIR / "smooth_reconstruction_timeseries.csv")
    recon = ts[ts.curve_kind == "reconstructed"]
    solved = summ[~summ.landmark_infeasible]
    print(f"re-screening {len(solved)} reconstructed cases ...")

    rows, assign_rows = [], []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for _, r in solved.iterrows():
            c, s, m, o = r.country, r.scenario, r.model, r.ordering
            curve = recon[(recon.country == c) & (recon.scenario == s) & (recon.model == m) & (recon.ordering == o)].sort_values("year")
            years = curve.year.to_numpy(int); raw = curve.rate_mt_yr.to_numpy(float)
            slug = slugify(c)
            cdir = P2_DIR / o / slug
            rc, sc_cache = cdir / "resource_rate_matrix_cache.csv", cdir / "site_number_matrix_cache.csv"
            if not rc.exists():
                rows.append({"country": c, "scenario": s, "model": m, "ordering": o, "rescreen_status": "no_cache"}); continue
            nr = len(pd.read_csv(rc))
            mdir = td / f"{slug}_{s}_{m.lower()}_{o}"; mdir.mkdir(parents=True, exist_ok=True)
            gc = mdir / "growth_curve.csv"
            pd.DataFrame({"year": years, "total_rate": raw}).to_csv(gc, index=False)
            try:
                run_allocation_workflow(AllocationConfig(
                    data_path=BASIN_FILE, output_dir=mdir, nr_region=nr, correction="off",
                    dist_min_km=2.0, dist_max_km="auto", nr_dist=100, nr_well_max="auto", rw_m=0.2,
                    max_q_mt_per_year=20.0, min_q_mt_per_year=1.0, inj_duration_yr=150,
                    allocation_duration_yr=150, allocation_order=ORDER_MAP[o],
                    storage_resource_calculation="savedfile", storage_period_step_yr=10,
                    growth_curve_path=gc, resource_rate_cache_path=rc, site_no_cache_path=sc_cache))
            except Exception as e:
                rows.append({"country": c, "scenario": s, "model": m, "ordering": o, "rescreen_status": f"error: {e}"}); continue
            ba, scr = reconstruct_screened(years, raw, mdir / "Resource_assignment_python.xlsx")
            cpm = {cp: cp_metrics(scr, ba, cp) for cp in CHECKPOINTS}
            full = cpm[2179]
            rows.append({"country": c, "scenario": s, "model": m, "ordering": o, "rescreen_status": "ok",
                         "rescreen_pass": bool(full.get("years_of_shortfall", 1) == 0),
                         "rescreen_total_unmet_gt": round(full.get("total_unmet_volume_gt", float("nan")), 4),
                         "rescreen_years_of_shortfall": full.get("years_of_shortfall"),
                         "rescreen_total_wells_2050": cpm[2050].get("total_wells"),
                         "rescreen_total_wells_2100": cpm[2100].get("total_wells"),
                         "rescreen_active_basins_2050": cpm[2050].get("active_basins"),
                         "rescreen_active_basins_2100": cpm[2100].get("active_basins")})
            if not ba.empty:
                for cp in CHECKPOINTS:
                    at = ba[ba.year == cp]
                    for b, g in at.groupby("basin"):
                        assign_rows.append({"country": c, "scenario": s, "model": m, "ordering": o,
                                            "checkpoint": cp, "basin": b,
                                            "delivered_rate_mt_yr": round(float(g.delivered_rate_mt_yr.sum()), 4),
                                            "wells_used": int(g.wells_used.iloc[0])})

    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "reconstructed_rescreen_summary.csv", index=False)
    pd.DataFrame(assign_rows).to_csv(OUT_DIR / "reconstructed_assignment.csv", index=False)
    ok = out[out.rescreen_status == "ok"]
    print(f"  ok: {len(ok)}/{len(out)}")
    if len(ok):
        print(f"  rescreen_pass: {int(ok.rescreen_pass.sum())}/{len(ok)}")
        print(f"  median rescreen unmet (Gt): {ok.rescreen_total_unmet_gt.median():.2f}")
    bad = out[out.rescreen_status != "ok"]
    if len(bad):
        print("  non-ok:", dict(bad.rescreen_status.value_counts()))
    print("wrote reconstructed_rescreen_summary.csv, reconstructed_assignment.csv")


if __name__ == "__main__":
    main()
