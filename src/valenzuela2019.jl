# interpolation is one of :interpolations, :dierckx
struct PostAGBTracks_MB16_Z01{T} <: AbstractPostAGBTracks
    interpolation::Symbol
    lookup::T
end

function PostAGBTracks_MB16_Z01(interpolation::Symbol=:interpolations)
    PostAGBTracks_MB16_Z01(interpolation, _get_lookup_2019(interpolation))
end

@kwdef struct AbsorbingFactor_V19{T} <: AbstractAbsorbingFactorMethod
    mumax::T=1
end

struct I5007_V19 <: AbstractI5007Method end

function _get_lookup_2019(interpolation::Symbol=:interpolations)
    if interpolation === :interpolations
        return _get_lookup_2019(LinearInterpolation)
    elseif interpolation === :dierckx
        return _get_lookup_2019(Spline2D)
    # elseif interpolation === :interpolations_mfinal
        # return get_lookup_2019_mfinal(LinearInterpolation)
    else
        error("The value `interpolation=:$(interpolation)` is invalid. Please pass `:interpolations` or `:dierckx`.")
    end
end

_getnpy(name) = npzread(_datadir("$name.npy"))
_savenpy(name, a) = npzwrite(_datadir("$name.npy"), a)


function _get_lookup_2019(::typeof(LinearInterpolation))
    filenames = ["lookup", "massinit", "age"]
    table, masses, ages = _getnpy.(filenames)
    itpLogT = LinearInterpolation((masses, ages), @view table[:, :, 1])
    itpLogL = LinearInterpolation((masses, ages), @view table[:, :, 2])
    return (; itpLogT, itpLogL)
end

function _get_lookup_2019_mfinal(::typeof(LinearInterpolation))
    filenames = ["lookup", "massx", "age"]
    table, masses, ages = _getnpy.(filenames)
    itpLogT = LinearInterpolation((masses, ages), @view table[:, :, 1])
    itpLogL = LinearInterpolation((masses, ages), @view table[:, :, 2])
    return (; itpLogT, itpLogL)
end

function _get_lookup_2019(::typeof(Spline2D))
    filenames = ["lookup", "massx", "age"]
    table, masses, ages = _getnpy.(filenames)
    itpLogT = Spline2D(masses, ages, @view table[:, :, 1])
    itpLogL = Spline2D(masses, ages, @view table[:, :, 2])
    return (; itpLogT, itpLogL)
end


_get_logT(lookup, mass, age) = lookup.itpLogT(mass, age * 1e9)
_get_logL(lookup, mass, age) = lookup.itpLogL(mass, age * 1e9)



function pn_luminosity(tracks::PostAGBTracks_MB16_Z01, minit::Real, mfinal::Real, postagb_age::Real, met::Real)
    (1 ≤ minit ≤ 3 && 0 ≤ postagb_age ≤ 3e-5) || return NaN
    return _get_logL(tracks.lookup, minit, postagb_age)
end
function pn_temperature(tracks::PostAGBTracks_MB16_Z01, minit::Real, mfinal::Real, postagb_age::Real, met::Real) 
    (1 ≤ minit ≤ 3 && 0 ≤ postagb_age ≤ 3e-5) || return NaN
    return _get_logT(tracks.lookup, minit, postagb_age)
end

function pn_luminosity_temperature(tracks::PostAGBTracks_MB16_Z01, minit::Real, mfinal::Real, postagb_age::Real, met::Real)
    logL = pn_luminosity(tracks, minit, mfinal, postagb_age, met)
    logT = pn_temperature(tracks, minit, mfinal, postagb_age, met)
    return logL, logT
end


function pn_absorbing_factor(mumethod::AbsorbingFactor_V19, minit::Real, mfinal::Real, postagb_age::Real, met::Real, logL::Real, logT::Real; rng=nothing)
    μ_generator = isnothing(rng) ? rand() : rand(rng)
    mumax = mumethod.mumax

    postagb_age_yr = postagb_age * 1e9

    # time to tip of the knee, interpolated over estimates for the five used masses
    time_to_knee = 0.02392279 * mfinal^-21.82379067

    # ranges to compute overlap fractions
    mfinal_range = (0.535, 0.557)
    # minit_range = (1.0215, 1.1607)
    temp_range = (4.51, 4.68)
    temp_range2 = (4.85, 4.95)
    age_range = (time_to_knee, 1.1time_to_knee)

    mass_overlap = _overlap(mfinal, mfinal_range)
    temp_overlap = _overlap(logT, temp_range)
    temp_overlap2 = _overlap(logT, temp_range2)
    age_overlap = _overlap(postagb_age_yr, age_range)

    μ_cooling = (0.1 + 0.9μ_generator^0.9) * (1 - postagb_age_yr / 30_000) # cooling tracks go transparent

    μ_heating_start = min(0.4 + μ_generator, 1) # beginning of heating tracks have more opaque PNe
    μ_heating_mid = # lighter ones almost transparent
        (0.2 + (mumax - 0.2) * (1 - μ_generator^2)) * mass_overlap +
        (0.03 + 0.17μ_generator) * (1 - mass_overlap)
    μ_heating_end = # lighter ones almost transparent
        (0.1 + (mumax - 0.1) * μ_generator) * mass_overlap + (0.03 + 0.17μ_generator) * (1 - mass_overlap)

    μ =
        μ_cooling * age_overlap +
        (
            (1 - temp_overlap) * μ_heating_start +
            temp_overlap * (1 - temp_overlap2) * μ_heating_mid +
            temp_overlap2 * μ_heating_end
        ) * (1 - age_overlap)
    μ = clamp(μ, 1e-8, 1)

    return μ
end

# relative position of x in xrange, returns a value in 0..1
function _overlap(x, xrange)
    x ≤ xrange[1] && return 0
    x ≥ xrange[2] && return 1
    return (x - xrange[1]) / (xrange[2] - xrange[1])
end

function pn_Ibeta(::I5007_V19, L::Real, Teff::Real)
    a = 1e-7
    b = 157_805.52 / Teff
    Sar = 2.404114 - quadgk(_fu, a, b)[1] # integration of fu from a to b

    Iβ = 1.70127e-4 * L * Sar / Teff # for a distance of 10pc
    # Mβ = magnitude(max(Iβ, 1e-15)) - 13.74 + magnitude(μ)

    return Iβ
end


function pn_I5007(method::I5007_V19, minit::Real, mfinal::Real, postagb_age::Real, met::Real, logL::Real, logT::Real; rng=nothing)
    Teff = exp10(logT)

    Iβ = pn_Ibeta(method, exp10(logL), Teff)

    if 30_152 < Teff ≤ 60_000 # 4.479 < logT ≤ 4.778
        I5007 = 0.0396Teff - 1194
    elseif 60_000 < Teff # logT > 4.778
        uniform_generator = isnothing(rng) ? rand() : rand(rng)
        gauss_generator = isnothing(rng) ? randn() : randn(rng)

        #### Random generation of variety in 5007, smaller for hot luminous stars with masses
        #### less than 0.575 (Stasinska 1989, AA 213, 274)
        I5gh = max(375gauss_generator + 950, 0)
        if I5gh > 1200
            I5gh -= 100(max(I5gh - 1200, 0) / 800)^1.5
        end
        if I5gh > 1500
            I5gh -= 100(max(I5gh - 1500, 0) / 500)
        end
        if I5gh > 1800
            I5gh = 1900 - 100uniform_generator
        end

        I5gh = min(I5gh, 1850)

        if mfinal ≤ 0.55 && logL > 3 && Teff > 75_000 # logT > 4.875
            I5gh *= 0.4
        end

        I5007 = I5gh
    else
        I5007 = 0.0
    end

    # intensity_ratio = I5007
    I5007 *= 0.01Iβ

    return I5007
end

_fu(x) = x^2 / (exp(x) - 1)
