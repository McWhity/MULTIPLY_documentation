Introduction
============
Priors are an essential component in the MULTIPLY inference engine as they provide a priori information
on different components of the unknown state vector of the system, helping to constrain the ill-posed problem given that the information content from the observations alone is insufficient. A series of prior models with different levels of complexity is therefore required and will be developed and implemented as part of the MULTIPLY platform.

The priors to be implemented are:

- Differential characterisation of the traits of vegetation types or (crop) species
- Vegetation phenology
- Surface soil moisture dynamics
- Surface disturbances

Background
-----------
A seamless and gap free integration of SENTINEL data streams requires the transfer of information across temporal and spatial scales. Typically data gaps are filled using low pass filters and different interpolation techniques (e.g. Savitzky-Golay filter; Savitzky & Golay, 1964) directly on parameter space (e.g. Yuan et al., 2011; Kandasamy et al. 2013). However, this approach is inconsistent, as the ill-posed nature of the inversion problem results in strong correlations between parameters: smoothing one parameter breaks that relationship with other retrieved parameters. Additionally, the role of uncertainty is usually ignored in filtering. Given that filtering methods originate from a prior belief in the smoothness of the processes that control the evolution of the parameters, it makes sense to implement these smoothness constraints consistently as priors within the retrieval process. These so-called regularisation constraints encompass our prior belief in the spatial and temporal correlation of the parameter fields. These constraints are implemented within the MULTIPLY platform as a weak constraint. The added benefit of having these constraints is that they not only result in smoother and more consistent series (an added benefit is an important reduction in parameter uncertainty), but also in spatially and temporally gap free estimates of biophysical parameters.

However, other prior information should be used to better constrain the inversion, and make sure that the inferences on the parameters are consistent with our understanding of biogeochemical processes and their effect on the state of the land surface.


Goal
-----
The major objectives of this software are i) to implement the required technical infrastructures to
provide the prior information at appropriate temporal and spatial scales in relation to the SENTINEL
observations, and ii) implement a flexible user interface which allows user to integrate own prior
models as a MULTIPLY plugin.

..
