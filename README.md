# CO2BLOCK: A tool for estimating storage resource capacity in saline aquifers

CO2BLOCK provides the automatic estimate of the storage capacity of an aquifer, based on the geomechanical limits to pressure build-up. 

Users only need to open and modify the script "CO2BLOCK.m". There is no need to open or modify the other scripts, and modifying the other scripts ‎not recommended.

Please refer to the user guide for more detailed instructions. 
If you want to learn more about the theory behind CO2BLOCK, please refer to the following papers:
- De Simone, S., & Krevor, S., 2021. A tool for first order estimates and optimisation of dynamic storage resource capacity in saline aquifers. _International Journal of Greenhouse Gas Control_, 106, 103258. <https://doi.org/10.1016/j.ijggc.2021.103258> 
- De Simone S., Jackson S. J. & Krevor S., 2019. The error in using superposition to estimate pressure during multi-site subsurface CO<sub>2</sub> storage. _Geophysical Research Letters_, 46, 6525-6533.  <https://doi.org/10.1029/2019GL082738>

# CO2BLOCK-SEISM
CO2BLOCK-SEISM extends CO2BLOCK in screening regional storage resources constrained by injection-induced seismicity. This extension primarily involves:
- calculation of the spatial and temporal evolution of fault slip potential,
- estimating potential magnitude and frequency of seismic events as key parameters controlling seismicity risk and hazard, and
- using a probabilistic framework to address uncertainties in subsurface geomechanical attributes.

To learn more about the theory behind CO2BLOCKSEISM, please refer to the following paper:
- Kivi, I.R., De Simone, S. and Krevor, S., 2025. A simplified physics model for estimating subsurface CO<sub>2</sub> storage resources constrained ‎by fault slip potential. _Journal of Rock Mechanics and Geotechnical Engineering_, 18(4), 2546-2560. <https://doi.org/10.1016/j.jrmge.2025.06.031>

# CO2BLOCK-GROWTH
CO2BLOCK-GROWTH extends CO2BLOCK by coupling techno-economic growth models with geophysical constraints to evaluate feasible CO<sub>2</sub> storage scaleup pathways. The framework combines Python scripts for generating logistic and Gompertz growth trajectories with MATLAB scripts for assessing reservoir pressure constraints and allocating storage resources. This extension primarily involves:
- generation of CO<sub>2</sub> storage deployment trajectories using logistic and Gompertz growth models (Python),
- assessment of reservoir pressurisation and injectivity constraints along the growth trajectories (MATLAB), and
- allocation of storage resources and identification of feasible storage pathways under coupled techno-economic growth and geophysical constraints (MATLAB).
To learn more about the theory behind CO2BLOCK-GROWTH, please refer to the following paper:
Kivi, I.R., Gao, X. and Krevor, S., 2026. Coupled geophysical and technoeconomic growth constraints on geological carbon storage scaleup with an application to the UK. (submitted).

Should you have any question/comment/suggestion, please contact Silvia De Simone at silviadesi@gmail.com, Iman R Kivi at iman.rahimzadeh@uib.es, or Sam Krevor at s.krevor@imperial.ac.uk
These tools are free licensed software.If you use it for academic purposes, please cite the reference paper De Simone and Krevor (2021) for CO2BLOCK, Kivi et al. (2025) for CO2BLOCK-SEISM and Kivi et al. (2026) for CO2BLOCK-GROWTH.

## Releases

Major versions of the software modules contained in this repository are archived using GitHub Releases.

For detailed descriptions of the capabilities, updates, and version history of each module, see the corresponding documentation:

- [CO2BLOCK](CO2BLOCK/README.md)
- [CO2BLOCKSEISM](CO2BLOCKSEISM/README.md)

Users requiring a specific version for reproducibility should download the corresponding GitHub Release.
