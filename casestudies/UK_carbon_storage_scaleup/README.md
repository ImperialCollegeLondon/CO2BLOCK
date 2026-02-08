# UK Carbon Storage Scale-up

This case study demonstrates the use of **CO2BLOCK** for estimating subsurface CO<sub>2</sub> storage resources in UK offshore storage units and for allocating these resources to national-scale carbon storage scale-up trajectories. The scale-up trajectories are derived from logistic growth and Gompertz models for the UK. The original CO2BLOCK code has been adapted to address this case study.

## Resource Estimates

In the folder `Storage_capacity_UK`, CO2BLOCK is used to estimate storage
capacity for 27 selected offshore storage units in the UK. Input data for all units are provided in the Excel spreadsheet `UK_sites`. Additional input parameters can be set directly in the code, as explained in the main CO2BLOCK README file.

The code produces
- Maximum storage capacities evolving with time for all units
- Contour plots of injectivity for two example units (Argyll 038 14 and Mey 5)
- A summary table of storage resources provided as an Excel spreadsheet

## Resource Allocation

In the folder `Resource_allocation_UK`, CO2BLOCK is used to screen the allocation of the estimated storage resources to UK CCS scale-up trajectories. Storage resources must be calculated for different injection durations. This can be done directly by running the code with the variable `storage_resource_calculation` set to `calculate`. The resulting storage resource estimates can be saved and reused. To run the allocation using previously calculated storage resources, the variable needs to be set to `savedfile`.

Storage resource estimates are already provided for five scenarios:
- `Storage_resources_400y_initialparameters`: 400 years of injection using most-likely input parameter values as described in Kivi et al. (2026),
- `Storage_resources_130y_initialparameters_maxsiteno100`: same input parameters as the reference scenario, but with a maximum of 100 injection sites per unit and an injection duration of 130 years,
- `Storage_resources_400y_reducedperm2`,
- `Storage_resources_400y_reducedperm5`, and
- `Storage_resources_400y_reducedperm10`.  
These scenarios use reduced permeability values (by factors of 2, 5, and 10) relative to the reference case to account for subsurface uncertainty, assuming 400 years of injection.

Scale-up trajectories for geological carbon storage in the UK are provided in the folder `Growth_Gompertz_models`. Each trajectory describes the temporal evolution of injection rate and cumulative stored CO<sub>2</sub>. Trajectories differ in their growth model (logistic or Gompertz), target cumulative storage volume, 2050 injection rate, and peak injection rate. Details are provided in `Format of file names.txt`.

The code allocates storage resources to the scale-up trajectories using the algorithm described in Kivi et al. (2026). By default, storage units with the largest estimated capacities are prioritised. Resource allocation figures are plotted by the code.


