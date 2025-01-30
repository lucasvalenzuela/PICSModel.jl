# PICSModel.jl

This Julia package is the implementation of PICS (Planetary nebulae In Cosmological Simulations),
a framework for modeling planetary nebulae (PNe) in single stellar populations (SSPs). This can be
applied to the stellar particles in cosmological simulations as well as measured star formation
histories obtained from observations. PICS can also be used to test the parameter space of different
models and their effect on the resulting PN populations.


## Documentation

Currently the code is only presented in its raw form. For more information on how to use the code at the moment,
please contact [Lucas Valenzuela](https://lucasvalenzuela.de/).
The documentation for PICS is still a work in progress and will be updated on Github as progress is made.


## Installation

PICSModel.jl is still not part of the Julia General registry. For the moment it can be
installed with the Julia package manager, but also requires three additional non-registered
packages. The following code will add PICSModel.jl and its dependencies to your environment
(to activate a new environment in the current directory, type `] activate .`, followed by
a backspace to return to the regular Julia mode):

```julia
using Pkg
Pkg.add(url="https://github.com/lucasvalenzuela/FastUnitfulOperations.jl")
Pkg.add(url="https://github.com/lucasvalenzuela/CosmoParticles.jl")
Pkg.add(url="https://github.com/lucasvalenzuela/AstroIMFs.jl")
Pkg.add(url="https://github.com/lucasvalenzuela/PICSModel.jl")
```


## References

The model is presented in the proceedings of the [IAU symposium 384](https://iaus384-pne.ncac.torun.pl/)
and will be introduced in detail in an upcoming paper.

If you use PICS in your work, please cite the following work:

[Valenzuela et al. 2024](https://ui.adsabs.harvard.edu/abs/2024arXiv241208702V/abstract):
```
@ARTICLE{2024arXiv241208702V,
       author = {{Valenzuela}, Lucas M. and {Remus}, Rhea-Silvia and {Miller Bertolami}, Marcelo M. and {M{\'e}ndez}, Roberto H.},
        title = "{PICS: Planetary Nebulae in Cosmological Simulations -- Revelations of the Planetary Nebula Luminosity Function from Realistic Stellar Populations}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Astrophysics of Galaxies, Astrophysics - Solar and Stellar Astrophysics},
         year = 2024,
        month = dec,
          eid = {arXiv:2412.08702},
        pages = {arXiv:2412.08702},
          doi = {10.48550/arXiv.2412.08702},
archivePrefix = {arXiv},
       eprint = {2412.08702},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024arXiv241208702V},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
