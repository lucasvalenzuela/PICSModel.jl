# TODO: obtain data from their paper
struct InitialFinalMass_MarigoGirardi07 <: AbstractInitialFinalMassRelation end

# TODO: obtain data from their paper
struct InitialFinalMass_Choi16 <: AbstractInitialFinalMassRelation end


struct InitialFinalMass_Cummings18 <: AbstractInitialFinalMassRelation
    mist_parsec::Symbol
    extrapolated::Bool

    function InitialFinalMass_Cummings18(mist_parsec::Symbol=:mist; extrapolated=true)
        @assert mist_parsec === :mist || mist_parsec === :parsec
        return new(mist_parsec, extrapolated)
    end
end

function final_mass(ifmr::InitialFinalMass_Cummings18, minit::Real, met::Real)
    if ifmr.mist_parsec === :mist
        _final_mass_cummings18_mist(ifmr, minit)
    else
        _final_mass_cummings18_parsec(ifmr, minit)
    end
end

function _final_mass_cummings18_mist(ifmr, minit)
    ifmr.extrapolated || 0.83 ≤ minit ≤ 7.20 || return NaN

    if minit < 2.85
        0.080minit + 0.489
    elseif minit < 3.60
        0.187minit + 0.184
    else
        0.107minit + 0.471
    end
end

function _final_mass_cummings18_parsec(ifmr, minit)
    ifmr.extrapolated || 0.87 ≤ minit ≤ 8.20 || return NaN

    if minit < 2.8
        0.0873minit + 0.476
    elseif minit < 3.65
        0.181minit + 0.210
    else
        0.0835minit + 0.565
    end
end


# from https://warwick.ac.uk/fac/sci/physics/research/astro/people/cunningham/simulation-data/
# const _minitC24 = [1.173, 1.205, 1.233, 1.261, 1.289, 1.319, 1.349, 1.377, 1.407, 1.436, 1.468, 1.499, 1.530, 1.564, 1.597, 1.633, 1.669, 1.707, 1.747, 1.788, 1.829, 1.871, 1.918, 1.966, 2.017, 2.071, 2.127, 2.187, 2.250, 2.315, 2.389, 2.463, 2.541, 2.625, 2.725, 2.826, 2.939, 3.050, 3.183, 3.327, 3.484, 3.673, 3.880, 4.114, 4.385, 4.698, 5.103, 5.551, 6.111, 6.817]
# const _mfinalC24 = [0.560, 0.564, 0.568, 0.572, 0.575, 0.577, 0.579, 0.583, 0.586, 0.589, 0.592, 0.595, 0.598, 0.601, 0.603, 0.605, 0.608, 0.611, 0.614, 0.616, 0.620, 0.624, 0.627, 0.634, 0.637, 0.641, 0.646, 0.651, 0.658, 0.663, 0.669, 0.677, 0.688, 0.694, 0.700, 0.712, 0.725, 0.738, 0.762, 0.782, 0.801, 0.813, 0.821, 0.845, 0.859, 0.877, 0.917, 0.976, 1.039, 1.330]


struct InitialFinalMass_Cunningham24{F} <: AbstractInitialFinalMassRelation
    extrapolated::Bool
    func::F

    function InitialFinalMass_Cunningham24(; extrapolated=true)
        func = linear_interpolation(_minitC24, _mfinalC24; extrapolation_bc=Line())
        return new{typeof(func)}(extrapolated, func)
    end
end

function final_mass(ifmr::InitialFinalMass_Cunningham24, minit::Real, met::Real)
    ifmr.extrapolated || 1.0 ≤ minit ≤ 7.67 || return NaN

    # from paper, but has discontinuities
    # if minit < 2.5
        # 0.086minit + 0.469
    # elseif minit < 3.4
        # 0.10minit + 0.40
    # elseif minit < 5.03
        # 0.06minit + 0.57
    # else
        # 0.17minit + 0.04
    # end

    return ifmr.func(minit)
end


const _minitH24 = [0.75, 1.00, 1.25, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00, 6.00, 8.00]
const _mfinalH24 = [
    [0.297, 0.552, 0.595, 0.614, 0.632, 0.666, 0.727, 0.803, 0.861, 0.909, 1.236, 4.679],
    [0.295, 0.551, 0.594, 0.613, 0.631, 0.657, 0.711, 0.760, 0.835, 0.875, 0.912, 1.101],
    [0.485, 0.552, 0.594, 0.613, 0.630, 0.658, 0.711, 0.763, 0.828, 0.872, 0.909, 1.101]
]

struct InitialFinalMass_Hollands24{F} <: AbstractInitialFinalMassRelation
    fit::Int
    extrapolated::Bool
    func::F

    function InitialFinalMass_Hollands24(fit=1; extrapolated=true)
        @assert fit == 1 || fit == 2 || fit == 3

        mfinal = _mfinalH24[fit]
        func = linear_interpolation(_minitH24, mfinal; extrapolation_bc=Line())
        return new{typeof(func)}(fit, extrapolated, func)
    end
end

function final_mass(ifmr::InitialFinalMass_Hollands24, minit::Real, met::Real)
    ifmr.extrapolated || 0.75 ≤ minit ≤ 8.00 || return NaN

    return ifmr.func(minit)
end
