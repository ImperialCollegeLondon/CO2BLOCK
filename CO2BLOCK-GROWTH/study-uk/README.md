# CO2BLOCK-GROWTH

This folder contains **CO2BLOCK-GROWTH**, an extension of CO2BLOCK that evaluates feasible geological CO<sub>2</sub> storage scaleup pathways by coupling technoeconomic growth models with geophysical constraints imposed by reservoir pressurisation.

## Code structure

The tool comprises a number of scripts and functions as detailed below:

- **CO2BLOCKGROWTH.m**: this is the main MATLAB script, where all input data and calculation options are specified. This script calls the required functions to perform the resource estimation and allocation calculations.
- **calculate.m**: this function estimates pressure buildup and pressure-limited storage resources for different arrangements of injection sites by varying the number of sites and the spacing between them for a specified injection duration.
- **calculate_timeseries.m**: this function performs the same calculations as **calculate.m**, but over a series of injection durations to evaluate the evolution of storage resources through time. This generates the matrix of storage resources for all geological storage units across different injection durations, which is the key input file for resource allocation. 
- **resource_allocation.m**: this function allocates pressure-constrained storage resources estimated above to logistic or Gompertz growth trajectories while satisfying the specified allocation criteria.
- **get_allocation_polygon.m**: this function extracts the allocation parameters required to generate resource allocation figures.
- **read_data.m**: this function reads the hydraulic and geological input parameters from the input Excel file.
- **Nordbotten_solution.m** and **FD_Nor.m**: these functions calculate pressure changes in a saline aquifer in response to two-phase flow of CO<sub>2</sub> and water around a single injection site based on the analytical solution by Nordbotten et al. (2005).
- **eos.m**: this function returns brine viscosity and CO<sub>2</sub> density and viscosity using appropriate equations of state.
- **Python scripts**: these scripts generate logistic and Gompertz growth trajectories, which are subsequently imported into MATLAB and combined with the geophysical calculations.

## Input

Input parameters describing the geological storage units are introduced either directly within **CO2BLOCKGROWTH.m** or through the input `.xlsx` files. The structure of these input data is generally the same as in CO2BLOCK and CO2BLOCK-SEISM, although data for all storage units are read simultaneously to enable regional calculations and resource allocation.

In addition to the geological input data, the tool requires growth trajectories generated using the accompanying Python scripts. These trajectories are imported into MATLAB and combined with the pressure-constrained storage resource calculations to evaluate feasible scaleup pathways.

Storage resources can be directly calculated or loaded from previous calculations. Previously saved storage resources list as well as growth trajectories are called from the folder "output" inside the main path. 

The example case study provided with the code contains input data for offshore saline aquifers in the United Kingdom.

## Output

The main outputs of the tool are pressure-constrained storage resources together with the corresponding resource allocation parameters, including the start and end times, injection duration, and contribution of each storage unit to the injection-rate and cumulative-storage growth trajectories. A summary report of all resource allocation parameters is provided in a `.xlsx` file.

The tool also produces the main figures used to analyse storage deployment, including:

- resource allocation along individual logistic or Gompertz growth trajectories, and
- comparison of 100 random resource allocation sequences with the two reference allocation strategies based on descending and ascending storage capacity.

The example input data included with this repository reproduce the UK offshore CO<sub>2</sub> storage scaleup case study presented in the accompanying publication.