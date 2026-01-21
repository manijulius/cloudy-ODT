# Overview
A Cloudy One-Dimensional Turbulence model with warm-cloud microphysics was developed to simulate  Rayleigh-Bénard convection (RBC) in the Pi cloud chamber at Michigan Technological University (MTU). The clouds in the chamber form when warm, moist air from the bottom plate mixes with cold, moist air from the top plate. Such clouds are referred to as mixing clouds, as opposed to clouds forming from the adiabatic cooling of a rising air parcel. The Pi cloud chamber, producing clouds in a laboratory setting, has been instrumental in understanding the effects of turbulence and aerosols on the cloud's microphysical properties. 

The moist-ODT model simulated temperature and water vapor fields, but has no aerosol or cloud droplets. Therefore, cloudy-ODT was developed over moist-ODT, implementing various processes, such as particle injection, condensational growth, gravitational settling, and particle displacement in eddies, to simulate a cloudy RBC in the chamber. 
* **Particle injection:** The model offers options to inject dry aerosol, deliquent aerosol, or small cloud droplets. The particle size distribution can be empirical, gamma, or monodisperse and supports three aerosol compositions (sodium chloride, ammonium sulphate, and ammonium bisulfate). These particles are randomly distributed within the specified domain volume, which can be the entire domain or a subvolume.
* **Condensation/Evaporation:** A particle in a grid cell either grows or reduces in size based on its local supersaturation.
* **Gravitational settling:** Gravitational settling based on Stokes' formulation and the density of air at the grid cell.
* **Eddy displacement:** triplet map used to move air parcel/grid cell scalar is also used to move the particles and assume no eddy hopping

In addition to the implementation of the above processes, various statistical outputs are available at regular intervals from simulations: frequency distributions of scalar and particles, mean and variance profiles of scalar and particle concentrations, evolution of domain-mean/domain-variance of scalar and particle concentrations, and particle trajectories with their ambient conditions.

This cloudy ODT implementation lacks side walls and therefore doesn't account for the effects of sidewall boundary conditions (temperature and moisture flux). However, it is useful to simulate a cloud chamber with an adiabatic sidewall.  The microphysical processes are based on the cloud Linear Eddy Model (LEM; Su et al  1998) and Explicit Mixing Parcel Model (EMPM; Krueger et al 1997)

# Contributions
* **Mani Rajagopal and Steve Krueger** implemented the aerosol injection and warm-cloud microphysical processes over moist-ODT code.

# References
* Su, C. W., Krueger, S. K., McMurtry, P. A., & Austin, P. H. (1998). Linear eddy modeling of droplet spectral evolution during entrainment and mixing in cumulus clouds. Atmospheric research, 47, 41-58.
* Krueger, S. K., Su, C. W., & McMurtry, P. A. (1997). Modeling entrainment and finescale mixing in cumulus clouds. Journal of the atmospheric sciences, 54(23), 2697-2712.
* Rajagopal 2024. cloudy ODT documentation v1.0 [here](https://github.com/manijulius/cloudy-ODT/blob/main/cloudy_ODT_documentation.pdf)

