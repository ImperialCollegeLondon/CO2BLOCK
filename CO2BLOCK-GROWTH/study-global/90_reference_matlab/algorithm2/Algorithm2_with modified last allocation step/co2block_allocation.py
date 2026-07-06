from __future__ import annotations

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from co2block_py.allocation import AllocationConfig, run_allocation_workflow


def main() -> None:
    base_dir = Path(__file__).resolve().parent

    config = AllocationConfig(
        data_path=base_dir / "Australia_sites.xlsx",
        output_dir=base_dir,
        nr_region=11,
        correction="off",
        dist_min_km=2.0,
        dist_max_km="auto",
        nr_dist=100,
        nr_well_max="auto",
        rw_m=0.2,
        max_q_mt_per_year=20.0,
        min_q_mt_per_year=1.0,
        inj_duration_yr=150,
        allocation_duration_yr=150,
        allocation_order="descend",
        storage_resource_calculation="calculate",
        storage_period_step_yr=10,
        growth_curve_path=base_dir / "growth_curve.csv",
    )

    outputs = run_allocation_workflow(config)
    print(f"Wrote {len(outputs['assignment_table'])} allocation steps to {base_dir}")


if __name__ == "__main__":
    main()
