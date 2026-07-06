#!/usr/bin/env python3
"""Step 2: Multi-scenario allocation core (standalone script).

Runs the CO₂BLOCK greedy allocation for all scenario × country × model
combinations and writes screening_summary.csv + Resource_assignment xlsx
files. Post-processing (checkpoint metrics, delay sensitivity, figures)
is done in step2_allocation_screening.ipynb with SKIP_ALLOCATION = True.

Run from code/ directory:
    python run_step2_screening.py
    python run_step2_screening.py --scenarios reference maximum
    python run_step2_screening.py --orders descend
"""
from __future__ import annotations

import argparse
import math
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings('ignore', category=FutureWarning)

CODE_ROOT = Path(__file__).resolve().parent
SCREENING_ROOT = CODE_ROOT.parent
sys.path.insert(0, str(CODE_ROOT))

from co2block_py.allocation import AllocationConfig, run_allocation_workflow
from co2block_py.screening import resolve_region_nos, slugify

# Paths
GROWTH_MODEL_ROOT = SCREENING_ROOT.parent / '01_growth_model'
GROWTH_TIMESERIES = GROWTH_MODEL_ROOT / 'output' / 'v8_2026-06-01' / 'timeseries_central.csv'
BASIN_FILE = SCREENING_ROOT / 'input' / 'basin_data' / 'Global.xlsx'
PRECOMPUTE_DIR = SCREENING_ROOT / 'output' / 'step1_precompute' / 'precompute'
OUTPUT_ROOT = SCREENING_ROOT / 'output' / 'step2_screening'

# Parameters (150yr, matches Step 1 & MATLAB Algorithm2)
INJ_DURATION_YR = 150
ALLOCATION_DURATION_YR = 150
STORAGE_PERIOD_STEP_YR = 10
NR_DIST = 100
MAX_Q_MT_YR = 20.0
MIN_Q_MT_YR = 1.0
SCREENING_SCOPE = 'country'

LAST_ALLOC_YEAR = 2030 + ALLOCATION_DURATION_YR - 1
CHECKPOINTS = [2050, 2100, LAST_ALLOC_YEAR]

GROWTH_PATHWAY_LABEL = 'v8_fixedC_deterministic'

ALL_SCENARIOS = [
    'reference', 'minimum', 'maximum', 'growth10',
    'us1gt', 'policy', 'ipcc_low', 'ipcc_high',
]

# V8 (2026-06-01) active study pool. EU is Europe without UK.
SELECTED_CASES = [
    'UK', 'US', 'EU', 'China', 'Middle East',
    'Australia', 'Canada', 'Indonesia', 'Thailand', 'Brazil',
]

MODELS = ['Logistic', 'Gompertz']
ORDER_DIR_MAP = {'descend': 'descending', 'ascend': 'ascending'}


def year_active_mask(years, start, end):
    return (years >= start) & (years <= end)


def compute_checkpoint_metrics(screened_df, basin_alloc_df, checkpoint_year, resource_long_df=None):
    mask = screened_df['year'] <= checkpoint_year
    sc = screened_df[mask]
    if sc.empty:
        return {}

    cp_row = sc[sc['year'] == checkpoint_year]
    if cp_row.empty:
        cp_row = sc.iloc[[-1]]

    annual_demand = float(cp_row['raw_rate_mt_yr'].iloc[0])
    annual_screened = float(cp_row['screened_rate_mt_yr'].iloc[0])
    annual_gap = float(cp_row['gap_rate_mt_yr'].iloc[0])
    r_annual = annual_screened / annual_demand if annual_demand > 0 else 0.0

    cum_demand = float(sc['raw_rate_mt_yr'].sum() / 1000.0)
    cum_screened = float(sc['screened_rate_mt_yr'].sum() / 1000.0)
    cum_gap = cum_demand - cum_screened
    r_cumulative = cum_screened / cum_demand if cum_demand > 0 else 0.0

    annual_pass = r_annual >= 1.0 - 1e-6
    cumulative_pass = r_cumulative >= 1.0 - 1e-6
    overall_pass = annual_pass and cumulative_pass

    shortfalls = sc[sc['gap_rate_mt_yr'] > 1e-9]
    first_shortfall_year = int(shortfalls['year'].iloc[0]) if not shortfalls.empty else None
    years_of_shortfall = len(shortfalls)
    max_gap_rate = float(sc['gap_rate_mt_yr'].max())
    max_gap_year = int(sc.loc[sc['gap_rate_mt_yr'].idxmax(), 'year']) if max_gap_rate > 1e-9 else None
    total_unmet_gt = float(sc['gap_rate_mt_yr'].clip(lower=0).sum() / 1000.0)

    active_basins = 0
    total_wells = 0
    n_eff = 0.0
    largest_share = 0.0
    dominant_basin = None
    total_basins = 0

    if not basin_alloc_df.empty:
        total_basins = basin_alloc_df['basin'].nunique()
        at_cp = basin_alloc_df[basin_alloc_df['year'] == checkpoint_year]
        if at_cp.empty:
            valid_years = basin_alloc_df[basin_alloc_df['year'] <= checkpoint_year]['year']
            if not valid_years.empty:
                at_cp = basin_alloc_df[basin_alloc_df['year'] == valid_years.max()]
        if not at_cp.empty:
            shares = at_cp.groupby('basin')['delivered_rate_mt_yr'].sum()
            active_basins = int((shares > 0).sum())
            total_wells = int(at_cp['wells_used'].sum())
            total = shares.sum()
            if total > 0:
                s = shares / total
                n_eff = float(1.0 / (s ** 2).sum())
                largest_share = float(s.max())
                dominant_basin = s.idxmax()

    return {
        'checkpoint': checkpoint_year,
        'annual_demand_mt_yr': round(annual_demand, 4),
        'annual_screened_mt_yr': round(annual_screened, 4),
        'annual_gap_mt_yr': round(annual_gap, 4),
        'r_annual': round(r_annual, 4),
        'cumulative_demand_gt': round(cum_demand, 4),
        'cumulative_screened_gt': round(cum_screened, 4),
        'cumulative_gap_gt': round(cum_gap, 4),
        'r_cumulative': round(r_cumulative, 4),
        'annual_pass': annual_pass,
        'cumulative_pass': cumulative_pass,
        'overall_pass': overall_pass,
        'first_shortfall_year': first_shortfall_year,
        'years_of_shortfall': years_of_shortfall,
        'max_gap_rate_mt_yr': round(max_gap_rate, 4),
        'max_gap_year': max_gap_year,
        'total_unmet_volume_gt': round(total_unmet_gt, 4),
        'active_basins': active_basins,
        'total_basins': total_basins,
        'active_basin_fraction': round(active_basins / total_basins, 4) if total_basins > 0 else 0.0,
        'total_wells': total_wells,
        'n_eff': round(n_eff, 2),
        'largest_basin_share': round(largest_share, 4),
        'dominant_basin_name': dominant_basin,
    }


def main():
    parser = argparse.ArgumentParser(description='Step 2: Multi-scenario screening')
    parser.add_argument('--scenarios', nargs='+', default=ALL_SCENARIOS)
    parser.add_argument('--orders', nargs='+', default=['descend', 'ascend'])
    args = parser.parse_args()

    scenarios = args.scenarios
    alloc_orders = args.orders

    # Load inputs
    resource_long = pd.read_csv(PRECOMPUTE_DIR / 'basin_period_resource_long.csv')
    metadata_df = pd.read_csv(PRECOMPUTE_DIR / 'basin_metadata.csv')
    growth_raw = pd.read_csv(GROWTH_TIMESERIES)
    required_growth_cols = {'Country', 'Scenario', 'Year', 'Model', 'rate_central'}
    missing_growth_cols = required_growth_cols.difference(growth_raw.columns)
    if missing_growth_cols:
        raise ValueError(
            f'Growth timeseries is missing required V8 columns: {sorted(missing_growth_cols)}'
        )

    screening_input_all = pd.DataFrame({
        'Country': growth_raw['Country'],
        'Scenario': growth_raw['Scenario'],
        'Year': growth_raw['Year'].astype(int),
        'Model': growth_raw['Model'],
        'Rate_Mt_yr': growth_raw['rate_central'] * 1000.0,
    })
    screening_input_all = screening_input_all[
        screening_input_all['Country'].isin(SELECTED_CASES)
    ].sort_values(['Scenario', 'Country', 'Model', 'Year']).reset_index(drop=True)

    # Build pivots
    periods = list(range(STORAGE_PERIOD_STEP_YR, INJ_DURATION_YR + 1, STORAGE_PERIOD_STEP_YR))
    rate_pivot = resource_long.pivot_table(
        index=['Region_no', 'Region_name'], columns='Period_yr',
        values='Region_Q_Mt_yr', aggfunc='first').reset_index()
    site_pivot = resource_long.pivot_table(
        index=['Region_no', 'Region_name'], columns='Period_yr',
        values='Number_of_Wells', aggfunc='first').reset_index()
    rate_pivot.columns = ['Region_no', 'Region_name'] + [
        f'Q [Mt/y] for t= {int(c)} y' for c in rate_pivot.columns[2:]]
    site_pivot.columns = ['Region_no', 'Region_name'] + [
        f'Q [Mt/y] for t= {int(c)} y' for c in site_pivot.columns[2:]]

    print(f'Growth:    {GROWTH_PATHWAY_LABEL}')
    print(f'Source:    {GROWTH_TIMESERIES}')
    print(f'Scenarios: {scenarios}')
    print(f'Orders:    {alloc_orders}')
    print(f'Countries: {SELECTED_CASES}')
    print(f'Models:    {MODELS}')
    print()

    # Allocation loop
    all_summaries = []
    total = len(alloc_orders) * len(scenarios) * len(SELECTED_CASES) * len(MODELS)

    for alloc_order in alloc_orders:
        order_label = ORDER_DIR_MAP[alloc_order]
        order_root = OUTPUT_ROOT / order_label
        for d in ['intermediate', 'final', 'figures']:
            (order_root / d).mkdir(parents=True, exist_ok=True)

        pbar = tqdm(total=len(scenarios) * len(SELECTED_CASES) * len(MODELS),
                    desc=f'{order_label}', unit='run')

        for scenario in scenarios:
            scenario_input = screening_input_all[screening_input_all['Scenario'] == scenario]
            for case in SELECTED_CASES:
                region_nos, basis = resolve_region_nos(metadata_df, case, scope=SCREENING_SCOPE)
                if region_nos:
                    cr = rate_pivot[rate_pivot['Region_no'].isin(region_nos)].copy()
                    cs = site_pivot[site_pivot['Region_no'].isin(region_nos)].copy()
                    cdir = order_root / slugify(case)
                    cdir.mkdir(parents=True, exist_ok=True)
                    rc = cdir / 'resource_rate_matrix_cache.csv'
                    sc = cdir / 'site_number_matrix_cache.csv'
                    cr.to_csv(rc, index=False)
                    cs.to_csv(sc, index=False)

                for model in MODELS:
                    pbar.update(1)
                    if not region_nos:
                        all_summaries.append({
                            'allocation_order': alloc_order, 'scenario': scenario,
                            'pathway': GROWTH_PATHWAY_LABEL, 'country': case, 'model': model,
                            'status': 'no_basin_identified', 'basin_match_basis': basis,
                            'basin_count': 0, 'peak_demand_mt_yr': None, 'allocated_peak_mt_yr': 0.0,
                        })
                        continue

                    cm = scenario_input[(scenario_input['Country'] == case) &
                                        (scenario_input['Model'] == model)].sort_values('Year')
                    if cm.empty:
                        all_summaries.append({
                            'allocation_order': alloc_order, 'scenario': scenario,
                            'pathway': GROWTH_PATHWAY_LABEL, 'country': case, 'model': model,
                            'status': 'no_growth_path', 'basin_match_basis': basis,
                            'basin_count': len(cr), 'peak_demand_mt_yr': None, 'allocated_peak_mt_yr': 0.0,
                        })
                        continue

                    mdir = order_root / slugify(case) / f'{scenario}_{model.lower()}'
                    mdir.mkdir(parents=True, exist_ok=True)
                    gc = mdir / 'growth_curve.csv'
                    cm[['Year', 'Rate_Mt_yr']].rename(
                        columns={'Year': 'year', 'Rate_Mt_yr': 'total_rate'}).to_csv(gc, index=False)

                    try:
                        out = run_allocation_workflow(AllocationConfig(
                            data_path=BASIN_FILE, output_dir=mdir, nr_region=len(cr),
                            correction='off', dist_min_km=2.0, dist_max_km='auto',
                            nr_dist=NR_DIST, nr_well_max='auto', rw_m=0.2,
                            max_q_mt_per_year=MAX_Q_MT_YR, min_q_mt_per_year=MIN_Q_MT_YR,
                            inj_duration_yr=INJ_DURATION_YR,
                            allocation_duration_yr=ALLOCATION_DURATION_YR,
                            allocation_order=alloc_order,
                            storage_resource_calculation='savedfile',
                            storage_period_step_yr=STORAGE_PERIOD_STEP_YR,
                            growth_curve_path=gc,
                            resource_rate_cache_path=rc, site_no_cache_path=sc,
                        ))
                        at = out['assignment_table']
                        pa = float(at['Step rate cumulative [Mt/y]'].max()) if not at.empty else 0.0
                        pd_val = float(cm['Rate_Mt_yr'].max())
                        status = 'ok'
                    except Exception as e:
                        pa, pd_val, status = 0.0, float(cm['Rate_Mt_yr'].max()), f'allocation_error: {e}'

                    all_summaries.append({
                        'allocation_order': alloc_order, 'scenario': scenario,
                        'pathway': GROWTH_PATHWAY_LABEL, 'country': case, 'model': model,
                        'status': status, 'basin_match_basis': basis,
                        'basin_count': len(cr), 'peak_demand_mt_yr': pd_val, 'allocated_peak_mt_yr': pa,
                    })

        pbar.close()

    # Save summaries
    sdf = pd.DataFrame(all_summaries)
    for ao in alloc_orders:
        ol = ORDER_DIR_MAP[ao]
        sdf[sdf['allocation_order'] == ao].to_csv(
            OUTPUT_ROOT / ol / 'final' / 'screening_summary.csv', index=False)

    ok = (sdf['status'] == 'ok').sum()
    print(f'\nDone: {ok} successful, {len(sdf) - ok} skipped')
    print('Run the notebook for post-processing, metrics, and figures.')


if __name__ == '__main__':
    main()
