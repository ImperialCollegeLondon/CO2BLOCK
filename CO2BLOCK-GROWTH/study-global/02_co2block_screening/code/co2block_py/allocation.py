from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .core import CalculationConfig, build_regional_summary, calculate_site, calculate_site_cached, rows_to_frame

# MATLAB fidelity (Algorithm2_with modified last allocation step / CO2BLOCK.m):
# once cumulative allocated rate reaches the growth peak, MATLAB recomputes the
# last region with minimum sites and immediately breaks the inner period loop
# (`if break_max_rate == 1; break`). The original Python port omitted that break,
# letting the period loop continue and overwrite the last step at a shorter
# period. Set True to match MATLAB exactly (recommended); False reproduces the
# legacy port behaviour for A/B comparison.
MATLAB_LAST_STEP_BREAK = True


def integrate_trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


@dataclass(frozen=True)
class AllocationConfig:
    data_path: Path
    output_dir: Path
    nr_region: int
    correction: str = "off"
    dist_min_km: float = 2.0
    dist_max_km: float | str = "auto"
    nr_dist: int = 100
    nr_well_max: int | str = "auto"
    rw_m: float = 0.2
    max_q_mt_per_year: float = 20.0
    min_q_mt_per_year: float = 1.0
    inj_duration_yr: int = 150
    allocation_duration_yr: int = 150
    allocation_order: str = "descend"
    storage_resource_calculation: str = "calculate"
    storage_period_step_yr: int = 10
    growth_curve_path: Path | None = None
    resource_rate_cache_path: Path | None = None
    site_no_cache_path: Path | None = None


def _interpolate_growth_value(years: np.ndarray, rates: np.ndarray, t: float) -> float:
    return float(np.interp(t, years, rates))


def _load_growth_curve(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        frame = pd.read_csv(path)
    elif path.suffix.lower() in {".xlsx", ".xls"}:
        frame = pd.read_excel(path)
    else:
        raise ValueError(f"Unsupported growth curve file format: {path}")
    if frame.shape[1] < 2:
        raise ValueError("Growth curve file must contain at least two columns: year and total rate.")
    return frame.iloc[:, :2].rename(columns={frame.columns[0]: "year", frame.columns[1]: "total_rate"})


def _load_cached_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    raise ValueError(f"Unsupported cache file format: {path}")


def _storage_periods(config: AllocationConfig) -> list[int]:
    return list(range(config.storage_period_step_yr, config.inj_duration_yr + 1, config.storage_period_step_yr))


def _calculate_resource_tables(config: AllocationConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    rate_rows: list[list[object]] = []
    site_rows: list[list[object]] = []
    periods = _storage_periods(config)

    for region_no in range(1, config.nr_region + 1):
        rate_row: list[object] = [region_no]
        site_row: list[object] = [region_no]
        region_name: str | None = None

        for time_yr in periods:
            calc = CalculationConfig(
                data_path=config.data_path,
                site_no=region_no,
                correction=config.correction,
                dist_min_km=config.dist_min_km,
                dist_max_km=config.dist_max_km,
                nr_dist=config.nr_dist,
                nr_well_max=config.nr_well_max,
                rw_m=config.rw_m,
                time_yr=time_yr,
                max_q_mt_per_year=config.max_q_mt_per_year,
                min_q_mt_per_year=config.min_q_mt_per_year,
                m0_guess_mt_per_year=None,
            )
            result = calculate_site_cached(calc)
            summary, _, _ = build_regional_summary(result, config.min_q_mt_per_year)
            if region_name is None:
                region_name = summary.region
                rate_row.append(region_name)
                site_row.append(region_name)
            rate_row.append(summary.region_q_mt_per_year)
            site_row.append(summary.number_of_wells)

        rate_rows.append(rate_row)
        site_rows.append(site_row)

    columns = ["Region_no", "Region_name"] + [f"Q [Mt/y] for t= {period} y" for period in periods]
    return pd.DataFrame(rate_rows, columns=columns), pd.DataFrame(site_rows, columns=columns)


def _sort_indices(matrix: np.ndarray, allocation_order: str) -> np.ndarray:
    if allocation_order == "descend":
        return np.argsort(np.sum(matrix, axis=1))[::-1]
    if allocation_order == "ascend":
        return np.argsort(np.sum(matrix, axis=1))
    if allocation_order == "random":
        rng = np.random.default_rng(0)
        return rng.permutation(matrix.shape[0])
    raise ValueError("allocation_order must be 'descend', 'ascend', or 'random'.")


def run_allocation_workflow(config: AllocationConfig) -> dict[str, pd.DataFrame]:
    if config.growth_curve_path is None:
        raise FileNotFoundError("A growth curve file is required for the allocation workflow.")

    config.output_dir.mkdir(parents=True, exist_ok=True)
    growth_curve = _load_growth_curve(config.growth_curve_path)
    years = growth_curve["year"].to_numpy(dtype=float)
    total_rate = growth_curve["total_rate"].to_numpy(dtype=float)
    limit_year = years[0] + config.allocation_duration_yr - 1
    mask = years <= limit_year
    years = years[mask]
    total_rate = total_rate[mask]

    if config.storage_resource_calculation == "calculate":
        rate_table, site_no_table = _calculate_resource_tables(config)
    elif config.storage_resource_calculation == "savedfile":
        if config.resource_rate_cache_path is None or config.site_no_cache_path is None:
            raise FileNotFoundError(
                "savedfile mode requires resource_rate_cache_path and site_no_cache_path."
            )
        rate_table = _load_cached_table(config.resource_rate_cache_path)
        site_no_table = _load_cached_table(config.site_no_cache_path)
    else:
        raise ValueError("storage_resource_calculation must be 'calculate' or 'savedfile'.")
    periods = _storage_periods(config)
    resource_rate_matrix = rate_table.iloc[:, 2:].to_numpy(dtype=float)
    site_no_matrix = site_no_table.iloc[:, 2:].to_numpy(dtype=float)
    region_name_array = rate_table.iloc[:, 1].to_numpy(dtype=object)
    region_no_array = rate_table.iloc[:, 0].to_numpy(dtype=int)

    sort_idx = _sort_indices(resource_rate_matrix, config.allocation_order)
    resource_rate_matrix_sorted = resource_rate_matrix[sort_idx]
    site_no_matrix_sorted = site_no_matrix[sort_idx]
    region_name_array_sorted = region_name_array[sort_idx]
    region_no_array_sorted = region_no_array[sort_idx]

    t_remaining = config.allocation_duration_yr
    total_rate_remain = total_rate.copy()

    assignments: list[dict[str, object]] = []
    storage_resource_total: list[float] = []
    cum_vol: list[float] = []
    step_curve_start: list[float] = []
    step_curve_end: list[float] = []

    for allocation_step in range(resource_rate_matrix_sorted.shape[0]):
        idx_candidates = np.where(np.array(periods) >= t_remaining)[0]
        if len(idx_candidates) == 0:
            break
        idx = int(idx_candidates[0])
        vol_possible = np.array(periods[: idx + 1], dtype=float) * resource_rate_matrix_sorted[allocation_step, : idx + 1]
        storage_resource_step = 0.0
        storage_resource_total_step = storage_resource_total[-1] if storage_resource_total else 0.0
        selected_period = float(periods[idx])
        selected_site_no = 0.0
        selected_curve_start = math.nan
        selected_curve_end = math.nan
        breaking = False
        break_max_rate = False

        for i in range(idx, -1, -1):
            total_rate_remain_test = total_rate_remain - resource_rate_matrix_sorted[allocation_step, i]
            non_negative = np.where(total_rate_remain_test >= 0)[0]
            left_idx = int(non_negative[0]) if len(non_negative) else None
            right_idx = int(non_negative[-1]) if len(non_negative) else None

            provisional_total = (
                resource_rate_matrix_sorted[allocation_step, i]
                if not storage_resource_total
                else storage_resource_total[-1] + resource_rate_matrix_sorted[allocation_step, i]
            )

            if left_idx is not None:
                if years[left_idx] > 2030:
                    provisional_curve_start = years[left_idx - 1] + (
                        (provisional_total - total_rate[left_idx - 1])
                        * (years[left_idx] - years[left_idx - 1])
                        / (total_rate[left_idx] - total_rate[left_idx - 1])
                    )
                else:
                    provisional_curve_start = 2030.0
                if years[right_idx] < years[-1]:
                    provisional_curve_end = years[right_idx] + (
                        (provisional_total - total_rate[right_idx])
                        * (years[right_idx + 1] - years[right_idx])
                        / (total_rate[right_idx + 1] - total_rate[right_idx])
                    )
                else:
                    provisional_curve_end = years[-1]

                area_time = np.concatenate(
                    ([provisional_curve_start], years[left_idx : right_idx + 1], [provisional_curve_end])
                )
                area_rate = np.concatenate(([0.0], total_rate_remain_test[non_negative], [0.0]))
                if not cum_vol:
                    area_growthrate = integrate_trapezoid(total_rate, years) - integrate_trapezoid(area_rate, area_time)
                else:
                    area_growthrate = integrate_trapezoid(total_rate, years) - cum_vol[-1] - integrate_trapezoid(area_rate, area_time)
            else:
                area_growthrate = integrate_trapezoid(total_rate, years) if not cum_vol else integrate_trapezoid(total_rate, years) - cum_vol[-1]

            if vol_possible[i] >= area_growthrate - 0.01 and resource_rate_matrix_sorted[allocation_step, i] >= storage_resource_step:
                storage_resource_step = float(resource_rate_matrix_sorted[allocation_step, i])
                storage_resource_total_step = float(provisional_total)
                selected_period = float(periods[i])
                selected_site_no = float(site_no_matrix_sorted[allocation_step, i])
                if left_idx is not None:
                    selected_curve_start = float(provisional_curve_start)
                    selected_curve_end = float(provisional_curve_end)
                if not cum_vol:
                    cum_volume = float(area_growthrate)
                else:
                    cum_volume = float(cum_vol[-1] + area_growthrate)

                if provisional_total >= float(np.max(total_rate)):
                    break_max_rate = True
                    remaining_rate_final = float(np.max(total_rate) - (storage_resource_total[-1] if storage_resource_total else 0.0))
                    calc = CalculationConfig(
                        data_path=config.data_path,
                        site_no=int(region_no_array_sorted[allocation_step]),
                        correction=config.correction,
                        dist_min_km=config.dist_min_km,
                        dist_max_km=config.dist_max_km,
                        nr_dist=config.nr_dist,
                        nr_well_max=config.nr_well_max,
                        rw_m=config.rw_m,
                        time_yr=periods[i],
                        max_q_mt_per_year=config.max_q_mt_per_year,
                        min_q_mt_per_year=config.min_q_mt_per_year,
                    )
                    result = calculate_site_cached(calc)
                    data_q = result.table_q.to_numpy(copy=True)
                    data_q_total = data_q * result.well_list.reshape(-1, 1)
                    data_v = result.table_v.to_numpy(copy=True)
                    valid_mask = (data_q_total > remaining_rate_final) & (data_v * 1000.0 > area_growthrate)
                    row_idx, col_idx = np.where(valid_mask)
                    if len(row_idx):
                        min_row = int(np.min(row_idx))
                        first_col = int(col_idx[np.where(row_idx == min_row)[0][0]])
                        storage_resource_step = float(data_q[min_row, first_col] * result.well_list[min_row])
                        storage_resource_total_step = float((storage_resource_total[-1] if storage_resource_total else 0.0) + storage_resource_step)
                        selected_site_no = float(result.well_list[min_row])

                if left_idx is None:
                    breaking = True

                final_cum_volume = cum_volume

            # MATLAB: `if break_max_rate == 1; break`. Lock the recomputed last
            # region and exit the period loop instead of scanning shorter periods.
            if MATLAB_LAST_STEP_BREAK and break_max_rate:
                break

        if storage_resource_step == 0.0:
            break

        total_rate_remain = total_rate_remain - storage_resource_step
        storage_resource_total.append(storage_resource_total_step)
        cum_vol.append(final_cum_volume)
        step_start = 2030.0 if allocation_step == 0 else step_curve_start[-1]
        step_end = 2030.0 + config.allocation_duration_yr - 1 if allocation_step == 0 else step_curve_end[-1]

        assignments.append(
            {
                "Step": allocation_step + 1,
                "Region no": int(region_no_array_sorted[allocation_step]),
                "Region name": str(region_name_array_sorted[allocation_step]),
                "Start [y]": step_start,
                "End [y]": step_end,
                "Duration [y]": selected_period,
                "Step rate increment [Mt/y]": storage_resource_step,
                "Step rate cumulative [Mt/y]": storage_resource_total_step,
                "No_sites": selected_site_no,
                "Curve start [y]": selected_curve_start,
                "Curve end [y]": selected_curve_end,
            }
        )

        if breaking or math.isnan(selected_curve_start) or math.isnan(selected_curve_end):
            break

        step_curve_start.append(selected_curve_start)
        step_curve_end.append(selected_curve_end)
        t_remaining = selected_curve_end - selected_curve_start + 1.0

    rate_output = config.output_dir / "Rate_opt_summary_python.xlsx"
    site_output = config.output_dir / "Site_no_summary_python.xlsx"
    assignment_output = config.output_dir / "Resource_assignment_python.xlsx"
    rate_cache_output = config.output_dir / "resource_rate_matrix_cache.xlsx"
    site_cache_output = config.output_dir / "site_number_matrix_cache.xlsx"
    if config.storage_resource_calculation == "calculate":
        rate_table.to_excel(rate_cache_output, index=False)
        site_no_table.to_excel(site_cache_output, index=False)
    rate_table.to_excel(rate_output, index=False)
    site_no_table.to_excel(site_output, index=False)
    assignment_table = pd.DataFrame(assignments)
    assignment_table.to_excel(assignment_output, index=False)
    return {
        "rate_table": rate_table,
        "site_no_table": site_no_table,
        "assignment_table": assignment_table,
    }
