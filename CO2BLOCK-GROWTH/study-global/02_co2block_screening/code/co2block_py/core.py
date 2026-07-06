from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    from scipy.special import lambertw as scipy_lambertw
except Exception:  # pragma: no cover
    scipy_lambertw = None


SECONDS_PER_YEAR = 86400 * 365


@dataclass(frozen=True)
class SiteData:
    site_name: str
    thickness_m: float
    area_km2: float
    permeability_m2: float
    porosity: float
    co2_density_kg_m3: float
    co2_viscosity_pa_s: float
    brine_viscosity_pa_s: float
    total_compressibility_pa_inv: float
    pressure_limit_mpa: float
    closed_boundary_radius_m: float
    viscosity_ratio: float
    delta: float
    omega: float


@dataclass(frozen=True)
class CalculationConfig:
    data_path: Path
    site_no: int
    correction: str = "off"
    dist_min_km: float = 2.0
    dist_max_km: float | str = "auto"
    nr_dist: int = 100
    nr_well_max: int | str = "auto"
    rw_m: float = 0.2
    time_yr: float = 100.0
    max_q_mt_per_year: float = 20.0
    min_q_mt_per_year: float = 1.0
    m0_guess_mt_per_year: float | None = None


@dataclass(frozen=True)
class CalculationResult:
    site_name: str
    d_list_km: np.ndarray
    well_list: np.ndarray
    d_max_km: np.ndarray
    q_m_each_mt_per_year: np.ndarray
    v_m_gt: np.ndarray
    table_q: pd.DataFrame
    table_v: pd.DataFrame
    p_sup_vec_mpa: np.ndarray


@dataclass(frozen=True)
class RegionalSummaryRow:
    region: str
    max_q_mt_per_year: float
    optimum_q_mt_per_year: float
    max_capacity_gt: float
    number_of_wells: int
    distance_km: float
    region_q_mt_per_year: float


def eos(temperature_c: float, pressure_mpa: float, salinity_fraction: float, co2_density_kg_m3: float) -> tuple[float, float, float]:
    """Return brine viscosity, CO2 density, and CO2 viscosity.

    This follows the MATLAB formulas directly. The CO2 density calculation
    solves the Redlich-Kwong cubic and selects the smallest positive real
    molar volume, which corresponds to the dense CO2 phase relevant here.
    """

    brine_viscosity = (
        0.1
        + 0.333 * salinity_fraction
        + (1.65 + 91.9 * salinity_fraction**3)
        * math.exp(-(0.42 * (salinity_fraction**0.8 - 0.17) ** 2 + 0.045) * temperature_c**0.8)
    ) / 1e3

    temperature_k = temperature_c + 273.15
    pressure_pa = pressure_mpa * 1e6

    a0 = 7.54
    a1 = -4.13e-3
    b = 2.78e-5
    a = a0 + a1 * temperature_k
    gas_constant = 8.314472

    coeffs = [
        1.0,
        -(gas_constant * temperature_k / pressure_pa),
        -(gas_constant * temperature_k * b / pressure_pa - a / pressure_pa / math.sqrt(temperature_k) + b**2),
        -(a * b / pressure_pa / math.sqrt(temperature_k)),
    ]
    roots = np.roots(coeffs)
    real_positive_roots = sorted(root.real for root in roots if abs(root.imag) < 1e-9 and root.real > 0)
    if not real_positive_roots:
        raise ValueError("No positive real root found for the CO2 density cubic equation.")
    molar_volume = real_positive_roots[0]
    solved_co2_density = 0.044 / molar_volume

    a10 = 0.248566120
    a11 = 0.004894942
    a20 = -0.373300660
    a21 = 1.22753488
    a30 = 0.363854523
    a31 = -0.774229021
    a40 = -0.0639070755
    a41 = 0.142507049
    reduced_temperature = temperature_k / 304.0
    reduced_density = co2_density_kg_m3 / 468.0
    mu_0 = math.sqrt(reduced_temperature) * (
        27.2246461 - 16.6346068 / reduced_temperature + 4.66920556 / (reduced_temperature**2)
    ) * 1e-6
    co2_viscosity = mu_0 * math.exp(
        a10 * reduced_density
        + a11 * reduced_density / reduced_temperature
        + a20 * reduced_density**2
        + a21 * reduced_density**2 / reduced_temperature
        + a30 * reduced_density**3
        + a31 * reduced_density**3 / reduced_temperature
        + a40 * reduced_density**4
        + a41 * reduced_density**4 / reduced_temperature
    )

    return brine_viscosity, solved_co2_density, co2_viscosity


def _matlab_scalar(value: object) -> float:
    if pd.isna(value):
        return float("nan")
    return float(value)


def read_site_data(data_path: Path | str, site_no: int) -> SiteData:
    data = pd.read_excel(data_path)
    row = data.iloc[site_no - 1]

    litho_grad = 23.0
    hydro_grad = 10.0
    temp_grad = 33.0
    def_k0 = 0.7
    def_friction_angle = 30.0
    def_cohesion = 0.0
    def_cr = 5e-4
    def_cw = 3e-4
    def_salinity = 180000.0

    site_name = str(row.iloc[0])
    domain_type = str(row.iloc[4]).strip().lower()
    depth = _matlab_scalar(row.iloc[5])
    depth_mean = _matlab_scalar(row.iloc[6])
    thickness_m = _matlab_scalar(row.iloc[7])
    area_km2 = _matlab_scalar(row.iloc[8])
    permeability_m2 = _matlab_scalar(row.iloc[9]) * 1e-15
    porosity = _matlab_scalar(row.iloc[10])
    rock_compressibility = _matlab_scalar(row.iloc[11]) / 1e6
    water_compressibility = _matlab_scalar(row.iloc[12]) / 1e6
    co2_density_kg_m3 = _matlab_scalar(row.iloc[13]) * 1e3
    co2_viscosity_pa_s = _matlab_scalar(row.iloc[14]) / 1e3
    brine_viscosity_pa_s = _matlab_scalar(row.iloc[15]) / 1e3
    pressure_top_mpa = _matlab_scalar(row.iloc[19])
    pressure_mean_mpa = _matlab_scalar(row.iloc[16])
    temperature_mean_c = _matlab_scalar(row.iloc[17])
    salinity_fraction = _matlab_scalar(row.iloc[18]) / 1e6
    principal_stress_total_mpa = _matlab_scalar(row.iloc[20])
    stress_ratio = _matlab_scalar(row.iloc[21])
    friction_deg = _matlab_scalar(row.iloc[22])
    cohesion_mpa = _matlab_scalar(row.iloc[23])
    tensile_strength_mpa = _matlab_scalar(row.iloc[24])

    if domain_type == "open":
        closed_boundary_radius_m = math.inf
    elif domain_type == "closed":
        closed_boundary_radius_m = math.sqrt(area_km2 * 1e6 / math.pi)
    else:
        raise ValueError(f"Unsupported domain type {domain_type!r} in {data_path}.")

    if pressure_top_mpa == 0 or math.isnan(pressure_top_mpa):
        pressure_top_mpa = hydro_grad * depth / 1000.0
    if pressure_mean_mpa == 0 or math.isnan(pressure_mean_mpa):
        pressure_mean_mpa = hydro_grad * depth_mean / 1000.0
    if temperature_mean_c == 0 or math.isnan(temperature_mean_c):
        temperature_mean_c = temp_grad * depth_mean / 1000.0 + 15.0
    if principal_stress_total_mpa == 0 or math.isnan(principal_stress_total_mpa):
        principal_stress_total_mpa = litho_grad * depth / 1000.0
    if stress_ratio == 0 or math.isnan(stress_ratio):
        stress_ratio = def_k0
    if friction_deg == 0 or math.isnan(friction_deg):
        friction_deg = def_friction_angle
    if cohesion_mpa == 0 or math.isnan(cohesion_mpa):
        cohesion_mpa = def_cohesion
    if tensile_strength_mpa == 0 or math.isnan(tensile_strength_mpa):
        tensile_strength_mpa = cohesion_mpa / 2.0
    if rock_compressibility == 0 or math.isnan(rock_compressibility):
        rock_compressibility = def_cr / 1e6
    if water_compressibility == 0 or math.isnan(water_compressibility):
        water_compressibility = def_cw / 1e6
    if salinity_fraction == 0 or math.isnan(salinity_fraction):
        salinity_fraction = def_salinity / 1e6
    if co2_density_kg_m3 == 0 or math.isnan(co2_density_kg_m3):
        _, co2_density_kg_m3, _ = eos(temperature_mean_c, pressure_mean_mpa, salinity_fraction, 0.0)
    if co2_viscosity_pa_s == 0 or math.isnan(co2_viscosity_pa_s):
        _, _, co2_viscosity_pa_s = eos(temperature_mean_c, pressure_mean_mpa, salinity_fraction, co2_density_kg_m3)
    if brine_viscosity_pa_s == 0 or math.isnan(brine_viscosity_pa_s):
        brine_viscosity_pa_s, _, _ = eos(temperature_mean_c, pressure_mean_mpa, salinity_fraction, 0.0)

    effective_max_stress_mpa = principal_stress_total_mpa - pressure_top_mpa
    effective_min_stress_mpa = stress_ratio * effective_max_stress_mpa
    friction_rad = math.radians(friction_deg)
    theta = (1.0 - math.sin(friction_rad)) / (1.0 + math.sin(friction_rad))
    pressure_limit_shear_mpa = (
        (stress_ratio - theta) / (1.0 - theta) * effective_max_stress_mpa
        + cohesion_mpa * math.cos(friction_rad) / math.sin(friction_rad)
    )
    pressure_limit_tensile_mpa = effective_min_stress_mpa + tensile_strength_mpa
    pressure_limit_mpa = min(pressure_limit_shear_mpa, pressure_limit_tensile_mpa)

    viscosity_ratio = co2_viscosity_pa_s / brine_viscosity_pa_s
    delta = (brine_viscosity_pa_s - co2_viscosity_pa_s) / brine_viscosity_pa_s
    omega = (
        (co2_viscosity_pa_s + brine_viscosity_pa_s)
        / (co2_viscosity_pa_s - brine_viscosity_pa_s)
        * math.log(math.sqrt(co2_viscosity_pa_s / brine_viscosity_pa_s))
        - 1.0
    )
    total_compressibility_pa_inv = rock_compressibility + porosity * water_compressibility

    return SiteData(
        site_name=site_name,
        thickness_m=thickness_m,
        area_km2=area_km2,
        permeability_m2=permeability_m2,
        porosity=porosity,
        co2_density_kg_m3=co2_density_kg_m3,
        co2_viscosity_pa_s=co2_viscosity_pa_s,
        brine_viscosity_pa_s=brine_viscosity_pa_s,
        total_compressibility_pa_inv=total_compressibility_pa_inv,
        pressure_limit_mpa=pressure_limit_mpa,
        closed_boundary_radius_m=closed_boundary_radius_m,
        viscosity_ratio=viscosity_ratio,
        delta=delta,
        omega=omega,
    )


def fd_nor(x_m: float, radius_m: float, external_radius_m: float) -> float:
    if x_m < external_radius_m:
        if radius_m <= external_radius_m:
            return math.log(radius_m / x_m)
        return math.log(external_radius_m / x_m) + 2.0 / 2.25 * (radius_m / external_radius_m) ** 2 - 0.75
    return 0.0


def nordbotten_solution(r_m: float, radius_m: float, psi_m: float, external_radius_m: float, gamma: float) -> float:
    if r_m <= psi_m:
        return gamma * math.log(psi_m / r_m) + fd_nor(psi_m, radius_m, external_radius_m)
    if psi_m < r_m <= radius_m:
        return fd_nor(r_m, radius_m, external_radius_m)
    return 0.0


# Vectorized versions for numpy arrays (same math, ~50-100x faster)

def _fd_nor_vec(x: np.ndarray, radius: float, ext_radius: float) -> np.ndarray:
    result = np.zeros_like(x)
    in_ext = x < ext_radius
    if not np.any(in_ext):
        return result
    if radius <= ext_radius:
        result[in_ext] = np.log(radius / x[in_ext])
    else:
        result[in_ext] = np.log(ext_radius / x[in_ext]) + 2.0 / 2.25 * (radius / ext_radius) ** 2 - 0.75
    return result


def _nordbotten_vec(r: np.ndarray, radius: float, psi: float, ext_radius: float, gamma: float) -> np.ndarray:
    result = np.zeros_like(r)
    fd_psi = _fd_nor_vec(np.full_like(r, psi), radius, ext_radius)
    mask1 = r <= psi
    if np.any(mask1):
        result[mask1] = gamma * np.log(psi / r[mask1]) + fd_psi[mask1]
    mask2 = (r > psi) & (r <= radius)
    if np.any(mask2):
        result[mask2] = _fd_nor_vec(r[mask2], radius, ext_radius)
    return result


def _pressure_superposition_vec(
    dist_vec: np.ndarray, influence_radius: float, psi: float,
    ext_radius: float, gamma: float, char_pressure: float,
) -> float:
    r_flat = dist_vec.ravel()
    return float(np.sum(_nordbotten_vec(r_flat, influence_radius, psi, ext_radius, gamma)) * char_pressure)


def _lambertw_branch_minus_one_scalar(x: float, max_iter: int = 100, tol: float = 1e-12) -> float:
    if x == 0.0:
        return float("-inf")

    if not (-1.0 / math.e <= x < 0.0):
        raise ValueError(f"Lambert W branch -1 requires x in [-1/e, 0), got {x}.")

    if abs(x + 1.0 / math.e) < 1e-15:
        return -1.0

    if x < -0.05:
        q = math.sqrt(2.0 * (math.e * x + 1.0))
        w = -1.0 - q - q**2 / 3.0
    else:
        log_term = math.log(-x)
        log_log_term = math.log(-log_term)
        w = log_term - log_log_term + log_log_term / log_term

    for _ in range(max_iter):
        e_w = math.exp(w)
        f = w * e_w - x
        denom = e_w * (w + 1.0) - (w + 2.0) * f / (2.0 * (w + 1.0))
        w_next = w - f / denom
        if abs(w_next - w) <= tol * (1.0 + abs(w_next)):
            return w_next
        w = w_next

    if abs(w * math.exp(w) - x) <= 1e-10 * max(1.0, abs(x)):
        return w
    raise RuntimeError(f"Lambert W branch -1 did not converge for x={x}.")


def lambertw_branch_minus_one(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    clipped = np.clip(arr, -1.0 / math.e, 0.0)
    if scipy_lambertw is not None:
        solved = scipy_lambertw(clipped, k=-1)
        return np.real(solved).astype(float, copy=False)
    flat = clipped.ravel()
    solved = np.array([_lambertw_branch_minus_one_scalar(value) for value in flat], dtype=float)
    return solved.reshape(arr.shape)


def calculate_site(config: CalculationConfig) -> CalculationResult:
    site = read_site_data(config.data_path, config.site_no)
    time_s = config.time_yr * SECONDS_PER_YEAR
    influence_radius_m = math.sqrt(
        2.246 * site.permeability_m2 * time_s / (site.brine_viscosity_pa_s * site.total_compressibility_pa_inv)
    )

    nr_well_max = int(math.floor(site.area_km2 / (config.dist_min_km**2))) if config.nr_well_max == "auto" else int(config.nr_well_max)
    dist_max_km = math.sqrt(2.0 * site.area_km2) / 2.0 if config.dist_max_km == "auto" else float(config.dist_max_km)
    d_list_km = np.linspace(config.dist_min_km, dist_max_km, int(config.nr_dist))
    # Match Iman's MATLAB calculate.m exactly: M0 = perm / 1e-13 [Mt/y]
    m0_guess = (
        config.m0_guess_mt_per_year
        if config.m0_guess_mt_per_year is not None
        else site.permeability_m2 / 1e-13
    )

    well_list: list[int] = []
    d_max_per_scenario_km: list[float] = []
    b_terms: list[float] = []
    q0_list: list[float] = []
    p_sup_rows: list[np.ndarray] = []

    max_side = int(math.sqrt(nr_well_max))
    for x_grid_num in range(1, max_side + 1):
        plus = 1 if x_grid_num * (x_grid_num + 1) < nr_well_max else 0
        for y_grid_num in range(x_grid_num, x_grid_num + plus + 1):
            well_count = x_grid_num * y_grid_num
            well_list.append(well_count)
            d_max_per_scenario_km.append(math.sqrt(site.area_km2 / well_count))

            p_sup_for_distance = np.zeros_like(d_list_km, dtype=float)
            q0 = m0_guess * 1e9 / site.co2_density_kg_m3 / SECONDS_PER_YEAR / well_count
            q0_list.append(q0)
            csi_m = math.sqrt(q0 * time_s / math.pi / site.porosity / site.thickness_m)
            psi_m = math.exp(site.omega) * csi_m
            characteristic_pressure_mpa = (q0 * site.brine_viscosity_pa_s) / (
                2.0 * math.pi * site.thickness_m * site.permeability_m2
            ) / 1e6

            if config.correction == "off":
                b_term = (site.brine_viscosity_pa_s - site.co2_viscosity_pa_s) / (
                    4.0 * math.pi * site.permeability_m2 * site.thickness_m
                )
            elif config.correction == "on":
                b_term = None
            else:
                raise ValueError("correction must be 'off' or 'on'.")

            # Precompute normalized grid distances (scale-independent)
            central_x = math.ceil(x_grid_num / 2) - 1
            central_y = math.ceil(y_grid_num / 2) - 1
            ix = np.tile(np.arange(x_grid_num, dtype=float), (y_grid_num, 1))
            iy = np.tile(np.arange(y_grid_num, dtype=float).reshape(-1, 1), (1, x_grid_num))
            norm_dist = np.sqrt((ix - central_x) ** 2 + (iy - central_y) ** 2).ravel()

            for d_idx, d_km in enumerate(d_list_km):
                distance_m = d_km * 1000.0
                dist_vec = norm_dist * distance_m
                dist_vec[central_y * x_grid_num + central_x] = config.rw_m

                p_sup = _pressure_superposition_vec(
                    dist_vec, influence_radius_m, psi_m,
                    site.closed_boundary_radius_m, site.viscosity_ratio,
                    characteristic_pressure_mpa,
                )

                if config.correction == "off":
                    sup_error = 0.0
                else:
                    if well_count < 9 or influence_radius_m * csi_m / (distance_m**2) < 1.0:
                        sup_error = 0.0
                        b_term = (site.brine_viscosity_pa_s - site.co2_viscosity_pa_s) / (
                            4.0 * math.pi * site.permeability_m2 * site.thickness_m
                        )
                    else:
                        sup_error = well_count * site.delta / 4.0 * math.log(influence_radius_m * csi_m / (distance_m**2))
                        b_term = (site.brine_viscosity_pa_s - site.co2_viscosity_pa_s) / (
                            4.0 * math.pi * site.permeability_m2 * site.thickness_m
                        ) * (1.0 + well_count / 4.0)

                p_sup_for_distance[d_idx] = p_sup - sup_error * characteristic_pressure_mpa

            b_terms.append(float(b_term))
            p_sup_rows.append(p_sup_for_distance)

    well_array = np.array(well_list, dtype=int)
    d_max_km = np.array(d_max_per_scenario_km, dtype=float)
    b = np.repeat(np.array(b_terms, dtype=float).reshape(-1, 1), len(d_list_km), axis=1)
    q1 = np.repeat(np.array(q0_list, dtype=float).reshape(-1, 1), len(d_list_km), axis=1)
    well_matrix = np.repeat(well_array.reshape(-1, 1), len(d_list_km), axis=1)
    p1_pa = np.vstack(p_sup_rows) * 1e6
    p2_pa = site.pressure_limit_mpa * 1e6

    lambert_argument = -p2_pa / q1 / b * np.exp(-p1_pa / q1 / b)
    lambert_value = lambertw_branch_minus_one(lambert_argument)
    q2_m3_s = -p2_pa / b / lambert_value
    q_m_each_mt_per_year = q2_m3_s * SECONDS_PER_YEAR * site.co2_density_kg_m3 / 1e9

    for d_idx, d_km in enumerate(d_list_km):
        distance_m = d_km * 1000.0
        if q_m_each_mt_per_year[0, d_idx] > 0.9999 * config.max_q_mt_per_year:
            q_m_each_mt_per_year[0, d_idx] = 0.9999 * config.max_q_mt_per_year
        for scenario_idx in range(1, len(well_array)):
            plume_cap = (distance_m**2) * math.pi * site.porosity * site.thickness_m / 4.001 / time_s
            if q_m_each_mt_per_year[scenario_idx, d_idx] > 0.9999 * config.max_q_mt_per_year or q2_m3_s[scenario_idx, d_idx] > plume_cap:
                q_m_each_mt_per_year[scenario_idx, d_idx] = min(
                    0.9999 * config.max_q_mt_per_year,
                    plume_cap * SECONDS_PER_YEAR * site.co2_density_kg_m3 / 1e9,
                )

    q_m_total_mt_per_year = q_m_each_mt_per_year * well_matrix
    v_m_gt = q_m_total_mt_per_year * config.time_yr / 1000.0

    d_mat = np.repeat(d_list_km.reshape(1, -1), len(well_array), axis=0)
    d_max_check = d_mat / d_max_km.reshape(-1, 1)
    possible = d_max_check < 1.0
    q_possible = q_m_each_mt_per_year * possible
    v_possible = v_m_gt * possible

    q_columns = [f"Q_M_for_d_{distance_km * 1000:.0f}_m" for distance_km in d_list_km]
    v_columns = [f"V_M_for_d_{distance_km * 1000:.0f}_m" for distance_km in d_list_km]
    index = pd.Index([str(well_count) for well_count in well_array], name="number_of_wells")

    table_q = pd.DataFrame(q_possible, index=index, columns=q_columns)
    table_v = pd.DataFrame(v_possible, index=index, columns=v_columns)

    return CalculationResult(
        site_name=site.site_name,
        d_list_km=d_list_km,
        well_list=well_array,
        d_max_km=d_max_km,
        q_m_each_mt_per_year=q_m_each_mt_per_year,
        v_m_gt=v_m_gt,
        table_q=table_q,
        table_v=table_v,
        p_sup_vec_mpa=np.vstack(p_sup_rows),
    )


# In-process memo for calculate_site. The function is deterministic in its
# config, so caching just skips recomputing identical (site, period) tables
# that the allocation inner loop would otherwise repeat.
_SITE_MEMO: dict = {}


def calculate_site_cached(config: CalculationConfig) -> CalculationResult:
    key = (
        str(config.data_path), int(config.site_no), config.correction,
        config.dist_min_km, str(config.dist_max_km), int(config.nr_dist),
        str(config.nr_well_max), config.rw_m, int(config.time_yr),
        config.max_q_mt_per_year, config.min_q_mt_per_year,
        config.m0_guess_mt_per_year,
    )
    result = _SITE_MEMO.get(key)
    if result is None:
        result = calculate_site(config)
        _SITE_MEMO[key] = result
    return result


def build_regional_summary(result: CalculationResult, min_q_mt_per_year: float) -> tuple[RegionalSummaryRow, pd.DataFrame, pd.DataFrame]:
    data_q = np.nan_to_num(result.table_q.to_numpy(copy=True), nan=0.0, posinf=0.0, neginf=0.0)
    data_v = np.nan_to_num(result.table_v.to_numpy(copy=True), nan=0.0, posinf=0.0, neginf=0.0)

    mask = data_q < min_q_mt_per_year
    mask[0, 0] = False
    data_q[mask] = 0.0
    data_v[mask] = 0.0

    max_q = float(np.max(data_q))
    linear_idx = int(np.argmax(data_v))
    max_vol = float(np.max(data_v))
    row_idx, col_idx = np.unravel_index(linear_idx, data_v.shape)
    optimum_q = float(data_q[row_idx, col_idx])
    well_count = int(result.well_list[row_idx])
    distance_km = round(float(result.d_list_km[col_idx]), 1)

    row = RegionalSummaryRow(
        region=result.site_name[:31],
        max_q_mt_per_year=max_q,
        optimum_q_mt_per_year=optimum_q,
        max_capacity_gt=max_vol,
        number_of_wells=well_count,
        distance_km=distance_km,
        region_q_mt_per_year=well_count * optimum_q,
    )

    q_columns = list(result.table_q.columns)
    v_columns = list(result.table_v.columns)
    q_masked = pd.DataFrame(data_q, index=result.table_q.index, columns=q_columns)
    v_masked = pd.DataFrame(data_v, index=result.table_v.index, columns=v_columns)
    return row, q_masked, v_masked


def rows_to_frame(rows: Iterable[RegionalSummaryRow]) -> pd.DataFrame:
    frame = pd.DataFrame([row.__dict__ for row in rows])
    if frame.empty:
        return frame
    return frame.rename(
        columns={
            "region": "Region",
            "max_q_mt_per_year": "Max Q [Mt/y]",
            "optimum_q_mt_per_year": "Optimum Q [Mt/y]",
            "max_capacity_gt": "Max capacity [Gt]",
            "number_of_wells": "Number of wells",
            "distance_km": "Distance",
            "region_q_mt_per_year": "Region Q [Mt/y]",
        }
    )
