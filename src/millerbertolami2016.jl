struct LifetimeFunction_MB16_extrapolated <: AbstractLifetimeFunction end
struct LifetimeFunction_MB16_He_linear <: AbstractLifetimeFunction end
struct LifetimeFunction_MB16_Z01 <: AbstractLifetimeFunction end
struct InitialFinalMass_MB16{F} <: AbstractInitialFinalMassRelation
    final_mass::F
    extrapolate::Bool
end
struct InitialFinalMass_MB16_Z01 <: AbstractInitialFinalMassRelation end



struct TrackInterpolator{T}
    met::Float64
    mfinals::Vector{Float64}
    tables::Vector{T}
end

struct MetallicityTrackInterpolator{T}
    mets::Vector{Float64}
    tracks::Vector{TrackInterpolator{T}}
end

struct PostAGBTracks_MB16{T} <: AbstractPostAGBTracks
    lookup::MetallicityTrackInterpolator{T}
end

function PostAGBTracks_MB16()
    filename = _datadir("tracks.h5")
    return PostAGBTracks_MB16(_loadinterpolator_MB16(filename))
end


const _minitMB16_Z01 = [1, 1.25, 1.5, 2, 2.5, 3]
const _ageMB16_Z01 = [
    7908.5 + 1697.2 + 113.53 + 11.965 + 0.73881 + 0
    3259.5 + 988.43 + 110.20 + 11.177 + 1.0972 + 0.14698e-1
    1889.2 + 400.18 + 101.54 + 11.439 + 1.2904 + 0.75163e-1
    921.34 + 27.072 + 282.61 + 12.747 + 1.6179 + 0.37651
    509.27 + 8.6930 + 148.01 + 7.6400 + 0.71424 + 0.41896
    318.67 + 4.0683 + 73.307 + 4.1321 + 0.31492 + 0.14332
] ./ 1000
const _pn_initial_mass_MB16_Z01 = extrapolate(interpolate(reverse(_ageMB16_Z01), reverse(_minitMB16_Z01), SteffenMonotonicInterpolation()), 0.0)

function pn_initial_mass(l::LifetimeFunction_MB16_Z01, ssp_age::Real, ssp_met=nothing)
    ssp_age < 0.015 && return NaN
    minit = _pn_initial_mass_MB16_Z01(ssp_age)
    1 ≤ minit ≤ 3 || return NaN

    return minit
end

# function ms_pms_lifetime(::LifetimeFunction_MB16_Z01, mass::Real, met=nothing)
    # mass < 0 && return Inf
# 
    # # values for met=0.01
    # a1 = 9.612654251129054
    # a2 = 6.44852715356046
    # m1 = -3.5260675850534087
    # m2 = -2.475677813865734
    # if mass < exp(log(a1 / a2) / (m2 - m1))
        # return a1 * mass^m1
    # else
        # return a2 * mass^m2
    # end
# end

# function pn_initial_mass(l::LifetimeFunction_MB16_Z01, ssp_age::Real, ssp_met=nothing)
    # func = let l=l
        # (mass, ssp_age) -> begin
            # val = ms_pms_lifetime(l, mass) - ssp_age
            # return val
        # end
    # end
# 
    # # use analytic MS dying mass approximation as starting point
    # p = ZeroProblem(func, ms_dyingmass(LifetimeFunction_P93_R81(), ssp_age))
    # return solve(p, Order0(); p=ssp_age) # cannot use Order1 here (results in very wrong values)
# end

function ms_pms_lifetime(::LifetimeFunction_MB16_He_linear, mass, met::Real)
    ms_pms_lifetime(LifetimeFunction_MB16_extrapolated(), mass, min(met, 0.02))
end


function ms_pms_lifetime(::LifetimeFunction_MB16_extrapolated, mass::Real, met::Real)
    mass < 0 && return Inf
    a1 = _get_lower_norm(met)
    a2 = _get_upper_norm(met)
    k1 = _get_lower_slope(met)
    k2 = _get_upper_slope(met)

    return max(a1 * mass^k1, a2 * mass^k2)
end
ms_pms_lifetime(l::LifetimeFunction_MB16_extrapolated, mass::Unitful.Mass, met::Real) = ms_pms_lifetime(l, ustrip(u"Msun", mass), met) * u"Gyr"

_get_lower_slope(met) = -0.9429375326887256met^0.3996632975377153 - 3.376390151587331
function _get_upper_slope(met)
    logmet = log10(met)
    return -3.2902508670185093 - 0.6543566409433365logmet - 0.1235350571834744logmet^2
end
_get_lower_norm(met) = 78.26687651498105met^0.648740756160864 + 5.6672047421884715
_get_upper_norm(met) = 107.5008350340392met^0.8803138562789512 + 4.5830749538317805

function pn_initial_mass(l::Union{LifetimeFunction_MB16_extrapolated,LifetimeFunction_MB16_He_linear}, ssp_age::Real, ssp_met::Real)
    ssp_age < 0.015 && return NaN

    func = let l=l
        (mass, (ssp_age, ssp_met)) -> begin
            val = ms_pms_lifetime(l, mass, ssp_met) - ssp_age
            return val
        end
    end

    # use analytic MS dying mass approximation as starting point
    p = ZeroProblem(func, ms_dyingmass(LifetimeFunction_P93_R81(), ssp_age))
    minit = solve(p, Order0(); p=(ssp_age, ssp_met)) # cannot use Order1 here (results in very wrong values)
    0.8 ≤ minit ≤ 4 || return NaN

    return minit
end





# create monotonic function minit -> mfinal for Z=0.01 for MB16
const _minitMB16_Z01_mono = [1, 1.25, 1.5, 2.5, 3]
const _mfinalMB16_Z01_mono = [0.5319, 0.566, 0.5832, 0.616, 0.7061]
const _final_mass_MB16_Z01 = extrapolate(interpolate(_minitMB16_Z01_mono, _mfinalMB16_Z01_mono, SteffenMonotonicInterpolation()), Flat())

function final_mass(::InitialFinalMass_MB16_Z01, minit::Real, met::Real) 
    1 ≤ minit ≤ 3 || return NaN
    return _final_mass_MB16_Z01(minit)
end




# full Miller Bertolami 2016 tables
const _minitMB16 = [ # in Msun
    # Z = 0.02
    [0.800 # fake point for better interpolation
        1.00
        1.25
        1.50
        2.00
        3.00
        4.00],
    # Z = 0.01
    [0.800 # fake point for better interpolation
        1.00
        1.25
        1.50
        2.00
        2.50
        3.00
        4.00], # fake point for better interpolation
    # Z = 0.001
    [0.800
        0.900
        1.00
        1.25
        1.50
        1.75
        2.00
        2.50
        3.00
        4.00], # fake point for better interpolation
    # Z = 0.0001
    [0.800
        0.850
        1.00
        1.50
        2.20
        2.50
        4.00], # fake point for better interpolation
]

const _mfinalMB16 = [ # in Msun
    # Z = 0.02
    [0.48 # fake point for better interpolation
        0.5281
        0.5615
        0.5760
        0.5804
        0.6573
        0.8328],
    # Z = 0.01
    [0.49 # fake point for better interpolation
        0.5319
        0.5660
        0.5832
        0.5826
        0.6160
        0.7061
        0.95], # fake point for better interpolation
    # Z = 0.001
    [0.4971
        0.5340
        0.5517
        0.5849
        0.5995
        0.5867
        0.6182
        0.7101
        0.8314
        1.20], # fake point for better interpolation
    # Z = 0.0001
    [0.5183
        0.5328
        0.5631
        0.6024
        0.7130
        0.7543
        1.26], # fake point for better interpolation
]

const _metMB16 = [
    0.02
    0.01
    0.001
    0.0001
]

const _ageMB16 = [ # in Gyr
    # Z = 0.02
    [25000
        9626.5 + 2043.9 + 131.60 + 11.161 + 0.65512 + 0.0000
        3857.6 + 1262.8 + 124.96 + 10.302 + 0.5268 + 0.89729 + 0.0000
        2212.4 + 480.67 + 114.64 + 11.508 + 0.5267 + 1.2959 + 0.16583e-1
        1055.0 + 33.928 + 307.42 + 19.420 + 0.4873 + 2.9196 + 0.97759e-1
        347.06 + 4.6762 + 83.942 + 6.2381 + 0.6103 + 0.79992 + 0.10622
        163.85 + 1.5595 + 28.196 + 1.9150 + 0.8086 + 0.32029 + 0.0000] ./ 1000,
    # Z = 0.01
    [20000
        7908.5 + 1697.2 + 113.53 + 11.965 + 0.73881 + 0.0000
        3259.5 + 988.43 + 110.20 + 11.177 + 1.0972 + 0.14698e-1
        1889.2 + 400.18 + 101.54 + 11.439 + 1.2904 + 0.75163e-1
        921.34 + 27.072 + 282.61 + 12.747 + 1.6179 + 0.37651
        509.27 + 8.6930 + 148.01 + 7.6400 + 0.71424 + 0.41896
        318.67 + 4.0683 + 73.307 + 4.1321 + 0.31492 + 0.14332
        170] ./ 1000,
    # Z = 0.001
    [12430 + 1313.5 + 116.80 + 13.089 + 0.4918 + 0.35321 + 0.0000
        8020.9 + 1024.7 + 102.38 + 10.345 + 0.5167 + 0.77626 + 0.0000
        5413.9 + 808.66 + 95.343 + 10.531 + 0.5282 + 0.85648 + 0.75481e-1
        2368.4 + 533.69 + 103.72 + 8.2357 + 0.5439 + 1.0575 + 0.15402
        1300.9 + 288.15 + 99.209 + 7.7074 + 0.5573 + 0.74257 + 0.46337
        951.65 + 46.973 + 210.27 + 8.2381 + 0.5298 + 1.1446 + 1.1904
        671.98 + 17.363 + 151.06 + 7.5069 + 0.5717 + 0.60336 + 1.0910
        389.98 + 5.9308 + 70.605 + 3.6796 + 0.6976 + 0.11696 + 0.49306
        255.51 + 2.9146 + 41.077 + 2.0364 + 0.8244 + 0.22804 + 0.82625e-1
        130] ./ 1000,
    # Z = 0.0001
    [11752 + 891.50 + 98.484 + 10.089 + 0.5049 + 0.77058 + 0.0000
        9425.1 + 768.48 + 93.850 + 9.5246 + 0.5128 + 1.1515 + 0.28265e-1
        5255.3 + 525.95 + 94.476 + 8.3594 + 0.5280 + 1.3575 + 0.92408e-1
        1261.4 + 241.51 + 112.01 + 8.1327 + 0.5526 + 0.76920 + 0.89965
        470.74 + 10.813 + 98.665 + 3.2623 + 0.7022 + 0.12220 + 0.48387
        349.85 + 5.8835 + 66.012 + 2.6371 + 0.7568 + 0.20994 + 0.20654
        100] ./ 1000,
]

function InitialFinalMass_MB16(; extrapolate=true)
    minits = _minitMB16
    mets = _metMB16
    mfinals = _mfinalMB16

    # individual 1D fits per metallicity
    itps = [DataInterpolations.QuadraticInterpolation(mfinal, minit) for (mfinal, minit) in zip(mfinals, minits)]

    # create additional data points on each metallicity curve
    m = range(0.8, 4, 50)
    mf = [itp(m) for itp in itps]

    # Interpolations method
    mfinal = reduce(vcat, mf)
    func = Interpolations.LinearInterpolation((m, reverse(log10.(mets))), reverse(reshape(mfinal, length(m), length(mets)); dims=2); extrapolation_bc=Line())

    #= Surrogates Kriging method (very slow)
    # prepare arrays for 2D grid interpolation
    minit = repeat(m; outer=length(mets))
    met = repeat(mets; inner=length(m))
    mfinal = reduce(vcat, mf)

    # 2D interpolation
    points = zip(minit, log10.(met)) |> collect
    ks = Kriging(points, mfinal, [minimum(minit), log10(minimum(mets))], [maximum(minit), log10(maximum(mets))]; p=[1.0, 1.0]) # p is a kind of smoothing parameter
    func = let ks = ks
        (minit, met) -> begin
            ks((minit, log10(met)))
        end
    end
    =#

    return InitialFinalMass_MB16(func, extrapolate)
end

function final_mass(ifm::InitialFinalMass_MB16, minit::Real, met::Real)
    0.8 ≤ minit ≤ 4 || return NaN
    # met_clamp = ifm.extrapolate ? met : clamp(met, 0.0001, 0.02)
    met_clamp = ifm.extrapolate ? max(met, 0.0001) : clamp(met, 0.0001, 0.02) # do not extrapolate towards low end
    return ifm.final_mass(minit, log10(met_clamp))
end






function Base.show(io::IO, mime::MIME"text/plain", itp::TrackInterpolator)
    printstyled(io, "TrackInterpolator"; bold=true)
    print(io, " at Z = $(itp.met)\n")
    mlow, mhigh = first(itp.mfinals), last(itp.mfinals)
    print(io, " for Mfinal = $mlow-$mhigh M⊙")
end

function Base.show(io::IO, itp::TrackInterpolator)
    print(io, "TrackInterpolator at Z = $(itp.met)")
end

function (itp::TrackInterpolator)(mfinal, t)
    first(itp.mfinals) ≤ mfinal ≤ last(itp.mfinals) || return NaN, NaN

    # find the tracks for the two relevant masses
    mind1 = min(searchsortedlast(itp.mfinals, mfinal), length(itp.mfinals) - 1)
    m1 = @inbounds itp.mfinals[mind1]
    m2 = @inbounds itp.mfinals[mind1 + 1]
    df1 = @inbounds itp.tables[mind1]
    df2 = @inbounds itp.tables[mind1 + 1]

    mfrac2 = _overlap(mfinal, (m1, m2))
    mfrac1 = 1 - mfrac2

    # search for index in tables at which the time fits
    t_lower = t_upper = zero(eltype(df1.t))
    tind1 = 0
    found = false
    for (t1, t2) in zip(df1.t, df2.t)
        t_upper = t1 * mfrac1 + t2 * mfrac2
        if t_upper > t
            found = true
            break
        end
        t_lower = t_upper
        tind1 += 1
    end
    found || return NaN, NaN

    tfrac2 = _overlap(t, (t_lower, t_upper))
    tfrac1 = 1 - tfrac2

    logL = @inbounds mfrac1 * (df1.logL[tind1] * tfrac1 + df1.logL[tind1 + 1] * tfrac2) + mfrac2 * (df2.logL[tind1] * tfrac1 + df2.logL[tind1 + 1] * tfrac2)
    logT = @inbounds mfrac1 * (df1.logTeff[tind1] * tfrac1 + df1.logTeff[tind1 + 1] * tfrac2) + mfrac2 * (df2.logTeff[tind1] * tfrac1 + df2.logTeff[tind1 + 1] * tfrac2)

    return logL, logT
end


function Base.show(io::IO, mime::MIME"text/plain", itp::MetallicityTrackInterpolator)
    printstyled(io, "MetallicityTrackInterpolator\n"; bold=true)
    metlow, methigh = first(itp.mets), last(itp.mets)
    print(io, " for Z = $metlow-$methigh")
end

function Base.show(io::IO, itp::MetallicityTrackInterpolator)
    print(io, "MetallicityTrackInterpolator")
end


# interpolate metallicity tracks
function _get_mets_MB16()
    return [0.0001, 0.001, 0.01, 0.02]
end

function _get_mfinals_MB16(met)
    if met == 0.0001
        return range(0.525, 0.8, 12)
    else
        return range(0.525, 1.0, 20)
    end
end

function _loadinterpolator_MB16(filename)
    mets = _get_mets_MB16()

    f = h5open(filename, "r")
    tracks = [_loadtrack_MB16(f["Z$met"], met) for met in mets]
    close(f)

    return MetallicityTrackInterpolator(collect(mets), tracks)
end

function _loadtrack_MB16(gmet, met)
    mfinals = _get_mfinals_MB16(met)
    tables = [_loadtable_MB16(gmet["M$mfinal"]) for mfinal in mfinals]
    return TrackInterpolator(met, collect(mfinals), tables)
end

function _loadtable_MB16(g) 
    t = read(g["t"])
    logL = read(g["logL"])
    logTeff = read(g["logTeff"])
    return (; t, logL, logTeff)
end


_is_mfinal_in_track_MB16(mfinal, track) = first(track.mfinals) ≤ mfinal ≤ last(track.mfinals)

function _get_mind_mfracs_MB16(mfinal, track)
    mind1 = min(searchsortedlast(track.mfinals, mfinal), length(track.mfinals) - 1)
    m1 = @inbounds track.mfinals[mind1]
    m2 = @inbounds track.mfinals[mind1 + 1]
    df1 = @inbounds track.tables[mind1]
    df2 = @inbounds track.tables[mind1 + 1]

    mfrac2 = _overlap(mfinal, (m1, m2))
    mfrac1 = 1 - mfrac2

    return mind1, mfrac1, mfrac2
end

function (itp::MetallicityTrackInterpolator)(mfinal, t, met)
    met = clamp(met, first(itp.mets), last(itp.mets)) # TODO: perhaps change to throw away low mets?

    # find the tracks for the two relevant metallicities
    zind1 = min(searchsortedlast(itp.mets, met), length(itp.mets) - 1)
    z1 = @inbounds itp.mets[zind1]
    z2 = @inbounds itp.mets[zind1 + 1]
    track1 = @inbounds itp.tracks[zind1]
    track2 = @inbounds itp.tracks[zind1 + 1]

    zfrac2 = _overlap(log10(met), (log10(z1), log10(z2)))
    zfrac1 = 1 - zfrac2

    if !_is_mfinal_in_track_MB16(mfinal, track1) || !_is_mfinal_in_track_MB16(mfinal, track2)
        return NaN, NaN
    end

    # find the tracks for the two relevant masses for each of the two metallicity tracks
    mind_z1_m1, mfrac_z1_m1, mfrac_z1_m2 = _get_mind_mfracs_MB16(mfinal, track1)
    mind_z2_m1, mfrac_z2_m1, mfrac_z2_m2 = _get_mind_mfracs_MB16(mfinal, track2)
    df_z1_m1, df_z1_m2 = @inbounds itp.tracks[zind1].tables[mind_z1_m1], itp.tracks[zind1].tables[mind_z1_m1 + 1]
    df_z2_m1, df_z2_m2 = @inbounds itp.tracks[zind1 + 1].tables[mind_z2_m1], itp.tracks[zind1 + 1].tables[mind_z2_m1 + 1]

    # search for index in tables at which the time fits
    t_lower = t_upper = zero(eltype(df_z1_m1.t))
    tind1 = 0
    found = false
    for (t11, t12, t21, t22) in zip(df_z1_m1.t, df_z1_m2.t, df_z2_m1.t, df_z2_m2.t)
        t_upper = zfrac1 * (t11 * mfrac_z1_m1 + t12 * mfrac_z1_m2) + zfrac2 * (t21 * mfrac_z2_m1 + t22 * mfrac_z2_m2)
        if t_upper > t
            found = true
            break
        end
        t_lower = t_upper
        tind1 += 1
    end
    found || return NaN, NaN

    tfrac2 = _overlap(t, (t_lower, t_upper))
    tfrac1 = 1 - tfrac2


    logL = @inbounds zfrac1 * (mfrac_z1_m1 * (df_z1_m1.logL[tind1] * tfrac1 + df_z1_m1.logL[tind1 + 1] * tfrac2) +
                               mfrac_z1_m2 * (df_z1_m2.logL[tind1] * tfrac1 + df_z1_m2.logL[tind1 + 1] * tfrac2)) +
                     zfrac2 * (mfrac_z2_m1 * (df_z2_m1.logL[tind1] * tfrac1 + df_z2_m1.logL[tind1 + 1] * tfrac2) +
                               mfrac_z2_m2 * (df_z2_m2.logL[tind1] * tfrac1 + df_z2_m2.logL[tind1 + 1] * tfrac2))
    logT = @inbounds zfrac1 * (mfrac_z1_m1 * (df_z1_m1.logTeff[tind1] * tfrac1 + df_z1_m1.logTeff[tind1 + 1] * tfrac2) +
                               mfrac_z1_m2 * (df_z1_m2.logTeff[tind1] * tfrac1 + df_z1_m2.logTeff[tind1 + 1] * tfrac2)) +
                     zfrac2 * (mfrac_z2_m1 * (df_z2_m1.logTeff[tind1] * tfrac1 + df_z2_m1.logTeff[tind1 + 1] * tfrac2) +
                               mfrac_z2_m2 * (df_z2_m2.logTeff[tind1] * tfrac1 + df_z2_m2.logTeff[tind1 + 1] * tfrac2))

    return logL, logT
end

function pn_luminosity(tracks::PostAGBTracks_MB16, minit::Real, mfinal::Real, postagb_age::Real, met::Real) 
    pn_luminosity_temperature(tracks, minit, mfinal, postagb_age, met)[1]
end

function pn_temperature(tracks::PostAGBTracks_MB16, minit::Real, mfinal::Real, postagb_age::Real, met::Real) 
    pn_luminosity_temperature(tracks, minit, mfinal, postagb_age, met)[2]
end

function pn_luminosity_temperature(tracks::PostAGBTracks_MB16, minit::Real, mfinal::Real, postagb_age::Real, met::Real)
    (0.525 ≤ mfinal ≤ 1 && -1e-6 ≤ postagb_age ≤ 1e-4) || return NaN, NaN
    lookup = tracks.lookup
    logL, logT = lookup(mfinal, postagb_age * 1e9, met)
    return logL, logT
end
