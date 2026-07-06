from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from co2block_py.core import CalculationConfig, build_regional_summary, calculate_site, read_site_data, rows_to_frame


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    data_path = base_dir / "Global.xlsx"

    nr_region = 203
    correction = "off"
    dist_min = 2.0
    dist_max = "auto"
    nr_dist = 200
    nr_well_max = "auto"
    rw = 0.2
    time_yr = 100
    max_q = 10.0
    min_q = 1.0

    tank_rows: list[dict[str, float | str]] = []
    summary_rows = []

    for region_no in range(1, nr_region + 1):
        site = read_site_data(data_path, region_no)
        max_vol_theory = site.area_km2 * 1e6 * site.thickness_m * site.porosity * site.co2_density_kg_m3 / 1e12
        max_vol_tank = site.area_km2 * 1e6 * site.thickness_m * site.co2_density_kg_m3 * site.total_compressibility_pa_inv * site.pressure_limit_mpa * 1e6 / 1e12
        tank_rows.append(
            {
                "Region": site.site_name[:31],
                "Max theoretical capacity [Gt]": max_vol_theory,
                "Max tank capacity [Gt]": max_vol_tank,
            }
        )

        calc = CalculationConfig(
            data_path=data_path,
            site_no=region_no,
            correction=correction,
            dist_min_km=dist_min,
            dist_max_km=dist_max,
            nr_dist=nr_dist,
            nr_well_max=nr_well_max,
            rw_m=rw,
            time_yr=time_yr,
            max_q_mt_per_year=max_q,
            min_q_mt_per_year=min_q,
            m0_guess_mt_per_year=20.0,
        )
        result = calculate_site(calc)
        summary, _, _ = build_regional_summary(result, min_q)
        summary_rows.append(summary)
        print(f"Calculations of region {region_no} completed")

    tank_table = pd.DataFrame(tank_rows)
    regional_table = rows_to_frame(summary_rows)
    tank_table.to_excel(base_dir / "Tank_model_summary_python.xlsx", index=False)
    regional_table.to_excel(base_dir / f"Regional_storage_summary_minQ={min_q}_{time_yr}y_python.xlsx", index=False)

    total_capacity = regional_table["Max capacity [Gt]"].sum()
    print(f"Total storage capacity = {total_capacity} [Gt]")


if __name__ == "__main__":
    main()
