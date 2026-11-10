# Geological carbon storage growth and basin allocation

Code for the paper:

> Gao, X., Rahimzadeh Kivi, I., Elizarev, M. & Krevor, S. *The geological limits to subsurface carbon dioxide storage emerge after 2050.* Manuscript ready for submission.

Contact: Xiaowei Gao (xg320@ic.ac.uk)

These notebooks reproduce the results and main figures of the paper from the input files in `input/`. Open the notebooks in `code/` and run them in order, with the working directory set to this `Global Study` folder or its `code/` folder. Each notebook writes to `Global Study/output/`, and the later notebooks read those files.

| Notebook | Content | Runtime |
|---|---|---|
| `code/01_growth_model.ipynb` | Logistic and Gompertz growth, Monte Carlo sampling and design pathways | ~1.5 min |
| `code/02_basin_allocation.ipynb` | Basin allocation, reconstruction and geological lifetime caps (closed and open boundaries) | ~10 s |
| `code/03_monte_carlo.ipynb` | Global distributions of injection rate and cumulative storage, and AR6 comparison data | ~15 s |
| `code/04_main_figures.ipynb` | Main-text Figures 1 to 5 (`output/figures/`) | ~35 s |

## Requirements

All four notebooks were tested in one environment:

Python 3.12, numpy 2.2, pandas 3.0, scipy 1.16, matplotlib 3.10, geopandas 1.0, shapely 2.0, pyproj 3.7, openpyxl 3.1, jupyter

Other package versions may change the last digits of the results or the figure rendering.

## Inputs

- `scenario_inputs.xlsx`: initial storage, storage resources and prescribed targets by jurisdiction
- `co2block_basin_resources_{closed,open}.csv`: CO2BLOCK output summaries of basin injection rates and resources for injection periods of 10 to 150 years
- `co2block_well_configurations_{closed,open}.csv`: alternative CO2BLOCK well arrays used to refine the last allocated basin
- `basin_properties_{closed,open}.xlsx`: geological and physical inputs to CO2BLOCK; also used to check basin order and boundary conditions
- `ar6_*`: IPCC AR6 scenario CCS rates
- `basin_map/global_sedimentary_basins.shp`: basin polygons, with `Basin Name` retained in the `.dbf` for matching to the model inputs; keep the `.dbf`, `.shx` and `.prj` files together
- `basin_map/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp`: Natural Earth country outlines, with the accompanying shapefile components and source documentation

## Scope

The notebooks start from precomputed CO2BLOCK outputs and do not run CO2BLOCK. With the inputs unchanged, they reproduce the published results. The well-configuration files cover only the basin and period combinations that the published pathways reach, so changed scenario inputs may require regenerating them with CO2BLOCK. See the last section of `code/02_basin_allocation.ipynb`.

The basin-level injection rates and storage resources can be reproduced using the [MATLAB implementation of CO2BLOCK](../../CO2BLOCK/README.md) with the supplied basin properties and the same model settings. These notebooks use the precomputed results provided in `input/`.
