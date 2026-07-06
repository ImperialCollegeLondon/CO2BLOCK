#!/usr/bin/env python3
"""Step 1: Precompute basin-period resource matrix (parallel, with tqdm).

Run from terminal:
    cd 02_co2block_screening/code
    python run_step1_precompute.py              # 4 workers (default)
    python run_step1_precompute.py --workers 8  # 8 workers
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from tqdm import tqdm

# Project paths
CODE_ROOT = Path(__file__).resolve().parent
SCREENING_ROOT = CODE_ROOT.parent
sys.path.insert(0, str(CODE_ROOT))

from co2block_py.core import CalculationConfig, build_regional_summary, calculate_site
from co2block_py.screening import assign_screening_case

# Paths
BASIN_FILE = SCREENING_ROOT / 'input' / 'basin_data' / 'Global.xlsx'
OUTPUT_ROOT = SCREENING_ROOT / 'output' / 'step1_precompute'
PRECOMPUTE_DIR = OUTPUT_ROOT / 'precompute'

# Parameters (aligned with Iman MATLAB Algorithm2)
INJ_DURATION_YR = 150
STORAGE_PERIOD_STEP_YR = 10
NR_DIST = 100
MAX_Q_MT_YR = 20.0
MIN_Q_MT_YR = 1.0

# V8 (2026-06-01) active study pool. EU is Europe without UK.
SCREENING_COUNTRIES = ['US', 'China', 'Indonesia', 'Australia', 'UK', 'Canada', 'Thailand', 'EU', 'Middle East', 'Brazil']


def safe_float(v):
    f = float(v)
    return f if math.isfinite(f) else 0.0


def evaluate_task(args):
    """Compute CO2BLOCK for one (basin, period) pair. Must be top-level for pickle."""
    region_no, period, basin_file = args
    calc = CalculationConfig(
        data_path=basin_file,
        site_no=region_no,
        correction='off',
        dist_min_km=2.0,
        dist_max_km='auto',
        nr_dist=NR_DIST,
        nr_well_max='auto',
        rw_m=0.2,
        time_yr=period,
        max_q_mt_per_year=MAX_Q_MT_YR,
        min_q_mt_per_year=MIN_Q_MT_YR,
    )
    result = calculate_site(calc)
    summary, _, _ = build_regional_summary(result, MIN_Q_MT_YR)
    return {
        'Region_no': region_no,
        'Region_name': summary.region,
        'Period_yr': period,
        'Region_Q_Mt_yr': safe_float(summary.region_q_mt_per_year),
        'Number_of_Wells': int(summary.number_of_wells),
        'Max_Capacity_Gt': safe_float(summary.max_capacity_gt),
        'Max_Q_Mt_yr': safe_float(summary.max_q_mt_per_year),
        'Optimum_Q_Mt_yr': safe_float(summary.optimum_q_mt_per_year),
        'Distance_km': safe_float(summary.distance_km),
    }


def main():
    parser = argparse.ArgumentParser(description='Step 1: Precompute basin resource matrix')
    parser.add_argument('--workers', type=int, default=4, help='Parallel workers (default: 4)')
    args = parser.parse_args()

    PRECOMPUTE_DIR.mkdir(parents=True, exist_ok=True)

    # Load basins
    basin_df = pd.read_excel(BASIN_FILE).copy()
    basin_df.insert(0, 'region_no', range(1, len(basin_df) + 1))
    basin_df['screening_case'] = basin_df.apply(
        lambda row: assign_screening_case(
            row.get('Majority Country'), row.get('Region'), scope='country'),
        axis=1,
    )
    assigned = basin_df[basin_df['screening_case'].isin(SCREENING_COUNTRIES)].copy()

    print(f'Basin database:  {BASIN_FILE}')
    print(f'Output dir:      {PRECOMPUTE_DIR}')
    print(f'Parameters:      nr_dist={NR_DIST}, maxQ={MAX_Q_MT_YR}, minQ={MIN_Q_MT_YR}')
    print(f'Workers:         {args.workers}')
    print()
    for c in SCREENING_COUNTRIES:
        n = len(assigned[assigned['screening_case'] == c])
        print(f'  {c:20s} {n:3d} basins')

    periods = list(range(STORAGE_PERIOD_STEP_YR, INJ_DURATION_YR + 1, STORAGE_PERIOD_STEP_YR))
    basins_to_compute = assigned['region_no'].unique()
    tasks = [(int(rno), period, str(BASIN_FILE)) for rno in basins_to_compute for period in periods]
    total = len(tasks)

    print(f'\nTotal tasks: {total}  ({len(basins_to_compute)} basins × {len(periods)} periods)')
    print()

    # Parallel compute with tqdm
    rows = []
    basin_names = basin_df.set_index('region_no').iloc[:, 0].to_dict()

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(evaluate_task, t): t for t in tasks}
        pbar = tqdm(as_completed(futures), total=total,
                    desc='CO₂BLOCK precompute',
                    bar_format='{l_bar}{bar:35}{r_bar}',
                    unit='task')
        basins_done = set()
        for future in pbar:
            result = future.result()
            rows.append(result)
            basins_done.add(result['Region_no'])
            name = basin_names.get(result['Region_no'], '?')
            pbar.set_postfix_str(
                f'basins: {len(basins_done)}/{len(basins_to_compute)} | '
                f'{name} @ {result["Period_yr"]}yr')

    # Save
    resource_long = pd.DataFrame(rows).sort_values(['Region_no', 'Period_yr']).reset_index(drop=True)

    meta_map = basin_df.set_index('region_no')[['Majority Country', 'Region']].to_dict('index')
    resource_long['AssignedCountry'] = resource_long['Region_no'].map(
        lambda x: meta_map.get(x, {}).get('Majority Country', 'Unknown'))
    resource_long['RegionGroup'] = resource_long['Region_no'].map(
        lambda x: meta_map.get(x, {}).get('Region', 'Unknown'))

    long_csv = PRECOMPUTE_DIR / 'basin_period_resource_long.csv'
    resource_long.to_csv(long_csv, index=False)

    metadata_df = resource_long[['Region_no', 'Region_name', 'AssignedCountry', 'RegionGroup']
        ].drop_duplicates().sort_values('Region_no')
    metadata_df.to_csv(PRECOMPUTE_DIR / 'basin_metadata.csv', index=False)

    print(f'\nDone: {len(resource_long):,} rows saved to {long_csv}')
    print(f'      {resource_long["Region_no"].nunique()} basins × {resource_long["Period_yr"].nunique()} periods')


if __name__ == '__main__':
    main()
