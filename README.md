# Overview
A Cloudy One-Dimensional Turbulence model with warm-cloud microphysics was developed to simulate  Rayleigh-Bénard convection (RBC) in the Pi cloud chamber at Michigan Technological University (MTU). The clouds in the chamber form when warm, moist air from the bottom plate mixes with cold, moist air from the top plate. Such clouds are referred to as mixing clouds, as opposed to clouds forming from the adiabatic cooling of a rising air parcel. The Pi cloud chamber, producing clouds in a laboratory setting, has been instrumental in understanding the effects of turbulence and aerosols on the cloud's microphysical properties. 

The moist-ODT model simulated temperature and water vapor fields, but has no aerosol or cloud droplets. Therefore, cloudy-ODT was developed over moist-ODT, implementing various processes, such as particle injection, condensational growth, gravitational settling, and particle displacement in eddies, to simulate a cloudy RBC in the chamber. 
* **Particle injection:** The model offers options to inject dry aerosol, deliquent aerosol, or small cloud droplets. The particle size distribution can be empirical, gamma, or monodisperse and supports three aerosol compositions (sodium chloride, ammonium sulphate, and ammonium bisulfate). These particles are randomly distributed within the specified domain volume, which can be the entire domain or a subvolume.
* **Condensation/Evaporation:** A particle in a grid cell either grows or reduces in size based on its local supersaturation.
* **Gravitational settling:** Gravitational settling based on Stokes' formulation and the density of air at the grid cell.
* **Eddy displacement:** triplet map used to move air parcel/grid cell scalar is also used to move the particles and assume no eddy hopping

In addition to the implementation of the above processes, various statistical outputs are available at regular intervals from simulations: frequency distributions of scalar and particles, mean and variance profiles of scalar and particle concentrations, evolution of domain-mean/domain-variance of scalar and particle concentrations, and particle trajectories with their ambient conditions.

This implementation of cloudy ODT lacks side walls and therefore doesn't account for the effects of sidewall boundary conditions (temperature and moisture flux). However, it is useful to simulate a cloud chamber with an adiabatic sidewall. 

# Contributions

# References
