module PICSModel

using AstroIMFs
using CosmoParticles
using DataDeps
using Dierckx
using HDF5
using Interpolations
using NPZ
using QuadGK
using Roots
# using Surrogates
using Unitful
using UnitfulAstro
using UnzipLoops

import DataInterpolations

export PNMethod,
    pn_population,
    pn_initial_mass,
    pn_number,
    pn_inds,
    pn_luminosity_temperature,
    pn_luminosity,
    pn_temperature,
    final_mass,
    pn_absorbing_factor,
    pn_I5007,
    absolute_magnitude_5007,
    InitialFinalMass_Cummings18,
    InitialFinalMass_Cunningham24,
    InitialFinalMass_Hollands24,
    PostAGBTracks_MB16_Z01,
    AbsorbingFactor_V19,
    I5007_V19,
    LifetimeFunction_P93_R81,
    LifetimeFunction_MB16_extrapolated,
    LifetimeFunction_MB16_He_linear,
    # LifetimeFunction_MB16_Z01,
    InitialFinalMass_MB16,
    # InitialFinalMass_MB16_Z01,
    PostAGBTracks_MB16


# load data
const datadepname = "PICS_data"

function __init__()
    # set up data dependency and download
    url = "https://www.usm.lmu.de/~lval/data/datadeps/PICS_data.zip"
    register(
        DataDep(
            datadepname,
            "The lookup table data used by the PICS model from Valenzuela et al. (2025), consisting of data from Miller Bertolami (2016) for the stellar evolution tracks and data from Valenzuela et al. (2019) for the planetary nebula model.",
            url,
            "ac9600e1c7c55693ee9c8e7d37027673a5809807cf3101390739a618b0b6c26d";
            post_fetch_method=unpack,
        ),
    )
    @datadep_str(datadepname)
end

_datadir(path...) = joinpath(@datadep_str(datadepname), path...)


include("api.jl")
include("particles.jl")
include("imf.jl")
include("cunningham24.jl")
include("ifmr.jl")
include("valenzuela2019.jl")
include("padovani1993_renzini1981.jl")
include("millerbertolami2016.jl")

end # module
