#!/usr/bin/env python3
"""Step 2: parallel-by-country allocation driver (v8 speedup).

Identical allocation and outputs to ``run_step2_screening.py``, but runs each
country in its own worker process. Combined with the in-process
``calculate_site`` memo (co2block_py.core.calculate_site_cached), each worker
computes its country's (basin, period) CO2BLOCK site tables once and reuses
them across that country's 32 runs (2 orders x 8 scenarios x 2 models), which
collapses the pathological 'descend' cost. Results are byte-for-byte identical
to the serial script (the cache is transparent and the allocation is
deterministic).

Writes the same layout: output/step2_screening/{order}/{country}/... plus a
per-order final/screening_summary.csv.
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

CODE_ROOT = Path(__file__).resolve().parent
SCREENING_ROOT = CODE_ROOT.parent
sys.path.insert(0, str(CODE_ROOT))
from co2block_py.allocation import AllocationConfig, run_allocation_workflow
from co2block_py.screening import resolve_region_nos, slugify

GROWTH_MODEL_ROOT = SCREENING_ROOT.parent / '01_growth_model'
GROWTH_TIMESERIES = GROWTH_MODEL_ROOT / 'output' / 'v8_2026-06-01' / 'timeseries_central.csv'
BASIN_FILE = SCREENING_ROOT / 'input' / 'basin_data' / 'Global.xlsx'
PRECOMPUTE_DIR = SCREENING_ROOT / 'output' / 'step1_precompute' / 'precompute'
OUTPUT_ROOT = SCREENING_ROOT / 'output' / 'step2_screening'

INJ_DURATION_YR = 150
ALLOCATION_DURATION_YR = 150
STORAGE_PERIOD_STEP_YR = 10
NR_DIST = 100
MAX_Q_MT_YR = 20.0
MIN_Q_MT_YR = 1.0
SCREENING_SCOPE = 'country'
GROWTH_PATHWAY_LABEL = 'v8_fixedC_deterministic'
ALL_SCENARIOS = ['reference', 'minimum', 'maximum', 'growth10',
                 'us1gt', 'policy', 'ipcc_low', 'ipcc_high']
SELECTED_CASES = ['UK', 'US', 'EU', 'China', 'Middle East',
                  'Australia', 'Canada', 'Indonesia', 'Thailand', 'Brazil']
MODELS = ['Logistic', 'Gompertz']
ALLOC_ORDERS = ['descend', 'ascend']
ORDER_DIR_MAP = {'descend': 'descending', 'ascend': 'ascending'}


def _load_pivots():
    resource_long = pd.read_csv(PRECOMPUTE_DIR / 'basin_period_resource_long.csv')
    metadata_df = pd.read_csv(PRECOMPUTE_DIR / 'basin_metadata.csv')
    rate_pivot = resource_long.pivot_table(index=['Region_no', 'Region_name'], columns='Period_yr',
                                           values='Region_Q_Mt_yr', aggfunc='first').reset_index()
    site_pivot = resource_long.pivot_table(index=['Region_no', 'Region_name'], columns='Period_yr',
                                           values='Number_of_Wells', aggfunc='first').reset_index()
    rate_pivot.columns = ['Region_no', 'Region_name'] + [f'Q [Mt/y] for t= {int(c)} y' for c in rate_pivot.columns[2:]]
    site_pivot.columns = ['Region_no', 'Region_name'] + [f'Q [Mt/y] for t= {int(c)} y' for c in site_pivot.columns[2:]]
    return metadata_df, rate_pivot, site_pivot


def _load_growth():
    g = pd.read_csv(GROWTH_TIMESERIES)
    si = pd.DataFrame({'Country': g['Country'], 'Scenario': g['Scenario'],
                       'Year': g['Year'].astype(int), 'Model': g['Model'],
                       'Rate_Mt_yr': g['rate_central'] * 1000.0})
    return si[si['Country'].isin(SELECTED_CASES)].sort_values(
        ['Scenario', 'Country', 'Model', 'Year']).reset_index(drop=True)


def process_country(country: str):
    """Run all order x scenario x model allocations for one country. Top-level for pickling."""
    metadata_df, rate_pivot, site_pivot = _load_pivots()
    growth = _load_growth()
    region_nos, basis = resolve_region_nos(metadata_df, country, scope=SCREENING_SCOPE)
    summaries = []
    for alloc_order in ALLOC_ORDERS:
        order_root = OUTPUT_ROOT / ORDER_DIR_MAP[alloc_order]
        for d in ('intermediate', 'final', 'figures'):
            (order_root / d).mkdir(parents=True, exist_ok=True)
        rc = sc = None
        cr = None
        if region_nos:
            cr = rate_pivot[rate_pivot['Region_no'].isin(region_nos)].copy()
            cs = site_pivot[site_pivot['Region_no'].isin(region_nos)].copy()
            cdir = order_root / slugify(country)
            cdir.mkdir(parents=True, exist_ok=True)
            rc = cdir / 'resource_rate_matrix_cache.csv'
            sc = cdir / 'site_number_matrix_cache.csv'
            cr.to_csv(rc, index=False)
            cs.to_csv(sc, index=False)
        for scenario in ALL_SCENARIOS:
            si = growth[growth['Scenario'] == scenario]
            for model in MODELS:
                base = {'allocation_order': alloc_order, 'scenario': scenario,
                        'pathway': GROWTH_PATHWAY_LABEL, 'country': country, 'model': model,
                        'basin_match_basis': basis}
                if not region_nos:
                    summaries.append({**base, 'status': 'no_basin_identified', 'basin_count': 0,
                                      'peak_demand_mt_yr': None, 'allocated_peak_mt_yr': 0.0})
                    continue
                cm = si[(si['Country'] == country) & (si['Model'] == model)].sort_values('Year')
                if cm.empty:
                    summaries.append({**base, 'status': 'no_growth_path', 'basin_count': len(cr),
                                      'peak_demand_mt_yr': None, 'allocated_peak_mt_yr': 0.0})
                    continue
                mdir = order_root / slugify(country) / f'{scenario}_{model.lower()}'
                mdir.mkdir(parents=True, exist_ok=True)
                gc = mdir / 'growth_curve.csv'
                cm[['Year', 'Rate_Mt_yr']].rename(columns={'Year': 'year', 'Rate_Mt_yr': 'total_rate'}).to_csv(gc, index=False)
                try:
                    out = run_allocation_workflow(AllocationConfig(
                        data_path=BASIN_FILE, output_dir=mdir, nr_region=len(cr),
                        correction='off', dist_min_km=2.0, dist_max_km='auto',
                        nr_dist=NR_DIST, nr_well_max='auto', rw_m=0.2,
                        max_q_mt_per_year=MAX_Q_MT_YR, min_q_mt_per_year=MIN_Q_MT_YR,
                        inj_duration_yr=INJ_DURATION_YR, allocation_duration_yr=ALLOCATION_DURATION_YR,
                        allocation_order=alloc_order, storage_resource_calculation='savedfile',
                        storage_period_step_yr=STORAGE_PERIOD_STEP_YR,
                        growth_curve_path=gc, resource_rate_cache_path=rc, site_no_cache_path=sc))
                    at = out['assignment_table']
                    pa = float(at['Step rate cumulative [Mt/y]'].max()) if not at.empty else 0.0
                    status = 'ok'
                except Exception as e:
                    pa, status = 0.0, f'allocation_error: {e}'
                summaries.append({**base, 'status': status, 'basin_count': len(cr),
                                  'peak_demand_mt_yr': float(cm['Rate_Mt_yr'].max()), 'allocated_peak_mt_yr': pa})
    return summaries


def main():
    import os
    workers = min(8, max(1, (os.cpu_count() or 4) - 1))
    print(f'Parallel step-2 v8: {len(SELECTED_CASES)} countries on {workers} workers')
    all_summaries = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(process_country, c): c for c in SELECTED_CASES}
        for fut in as_completed(futs):
            c = futs[fut]
            res = fut.result()
            all_summaries.extend(res)
            print(f'  done: {c}  ({len(res)} runs)', flush=True)
    sdf = pd.DataFrame(all_summaries)
    cols = ['allocation_order', 'scenario', 'pathway', 'country', 'model', 'status',
            'basin_match_basis', 'basin_count', 'peak_demand_mt_yr', 'allocated_peak_mt_yr']
    sdf = sdf[cols]
    for ao in ALLOC_ORDERS:
        sub = sdf[sdf['allocation_order'] == ao].sort_values(['scenario', 'country', 'model']).reset_index(drop=True)
        sub.to_csv(OUTPUT_ROOT / ORDER_DIR_MAP[ao] / 'final' / 'screening_summary.csv', index=False)
    ok = int((sdf['status'] == 'ok').sum())
    print(f'\nDone: {ok} ok, {len(sdf) - ok} skipped/error. Run the notebook (SKIP_ALLOCATION=True) for figures.')


if __name__ == '__main__':
    main()
