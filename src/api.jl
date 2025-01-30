abstract type AbstractPNMethodType end
Broadcast.broadcastable(m::AbstractPNMethodType) = Ref(m)

abstract type AbstractLifetimeFunction <: AbstractPNMethodType end
abstract type AbstractPostAGBTracks <: AbstractPNMethodType end
abstract type AbstractInitialFinalMassRelation <: AbstractPNMethodType end
abstract type AbstractAbsorbingFactorMethod <: AbstractPNMethodType end
abstract type AbstractI5007Method <: AbstractPNMethodType end

struct PNMethod{IMF<:AbstractIMF,LF<:AbstractLifetimeFunction,PT<:AbstractPostAGBTracks,IF<:AbstractInitialFinalMassRelation,AF<:AbstractAbsorbingFactorMethod,I5<:AbstractI5007Method}
    imf::IMF
    lifetime::LF
    postagb_tracks::PT
    initial_final_mass::IF
    absorbing_factor::AF
    i5007::I5
end

function PNMethod(imf::AbstractIMF, type::Symbol)
    postagb_tracks = PostAGBTracks_MB16_Z01(:interpolations)
    initial_final_mass = InitialFinalMass_MB16_Z01()
    absorbing_factor = AbsorbingFactor_V19(1)
    i5007 = I5007_V19()
    if type === :padovani
        lifetime = LifetimeFunction_P93_R81()
    elseif type === :millerbertolami
        lifetime = LifetimeFunction_MB16_Z01()
    elseif type === :millerbertolamiz
        postagb_tracks = PostAGBTracks_MB16()
        initial_final_mass = InitialFinalMass_MB16()
        lifetime = LifetimeFunction_MB16_extrapolated()
    else
        error("The type `:$type` is not supported. Available types are `:padovani`, `:millerbertolami`, and `:millerbertolamiz`.")
    end
    return PNMethod(imf, lifetime, postagb_tracks, initial_final_mass, absorbing_factor, i5007)
end

function Base.show(io::IO, mime::MIME"text/plain", itp::PNMethod)
    printstyled(io, "PNMethod"; bold=true)
end



function ms_lifetime(::AbstractLifetimeFunction, minit, met) end
ms_lifetime(l::AbstractLifetimeFunction, minit::Unitful.Mass, met) = ms_lifetime(l, ustrip(u"Msun", minit), met)

function pms_lifetime(::AbstractLifetimeFunction, minit, met) end
pms_lifetime(l::AbstractLifetimeFunction, minit::Unitful.Mass, met) = ms_lifetime(l, ustrip(u"Msun", minit), met)

function ms_pms_lifetime(l::AbstractLifetimeFunction, minit, met)
    ms_lifetime(l, minit, met) + pms_lifetime(l, minit, met)
end
ms_pms_lifetime(l::AbstractLifetimeFunction, minit::Unitful.Mass, met) = ms_lifetime(l, ustrip(u"Msun", minit), met)

function pn_initial_mass(::AbstractLifetimeFunction, ssp_age, ssp_met) end
function pn_initial_mass(l::AbstractLifetimeFunction, ssp_age::Unitful.Time, ssp_met)
    pn_initial_mass(l, ustrip(u"Gyr", ssp_age), ssp_met) * u"Msun"
end

function final_mass(::AbstractInitialFinalMassRelation, minit, met) end
function final_mass(fmr::AbstractInitialFinalMassRelation, minit::Unitful.Mass, met)
    final_mass(fmr, ustrip(u"Msun", minit), met) * u"Msun"
end

function pn_luminosity(::AbstractPostAGBTracks, minit, mfinal, postagb_age, met) end
function pn_luminosity(pt::AbstractPostAGBTracks, minit::Unitful.Mass, mfinal::Unitful.Mass, postagb_age::Unitful.Time, met)
    pn_luminosity(pt, ustrip(u"Msun", minit), ustrip(u"Msun", mfinal), ustrip(u"Gyr", postagb_age), met)
end

function pn_temperature(::AbstractPostAGBTracks, minit, mfinal, postagb_age, met) end
function pn_temperature(pt::AbstractPostAGBTracks, minit::Unitful.Mass, mfinal::Unitful.Mass, postagb_age::Unitful.Time, met)
    pn_temperature(pt, ustrip(u"Msun", minit), ustrip(u"Msun", mfinal), ustrip(u"Gyr", postagb_age), met)
end

function pn_luminosity_temperature(pt::AbstractPostAGBTracks, minit, mfinal, postagb_age, met)
    pn_luminosity(pt, minit, mfinal, postagb_age, met), pn_temperature(pt, minit, mfinal, postagb_age, met)
end
function pn_luminosity_temperature(pt::AbstractPostAGBTracks, minit::Unitful.Mass, mfinal::Unitful.Mass, postagb_age::Unitful.Time, met)
    pn_luminosity_temperature(pt, ustrip(u"Msun", minit), ustrip(u"Msun", mfinal), ustrip(u"Gyr", postagb_age), met)
end

function pn_absorbing_factor(::AbstractAbsorbingFactorMethod, minit, mfinal, postagb_age, met, logL, logT) end
function pn_absorbing_factor(pt::AbstractAbsorbingFactorMethod, minit::Unitful.Mass, mfinal::Unitful.Mass, postagb_age::Unitful.Time, met, logL, logT; kwargs...)
    pn_absorbing_factor(pt, ustrip(u"Msun", minit), ustrip(u"Msun", mfinal), ustrip(u"Gyr", postagb_age), met, logL, logT; kwargs...)
end

function pn_I5007(pt::AbstractI5007Method, minit, mfinal, postagb_age, met, logL, logT) end
function pn_I5007(pt::AbstractI5007Method, minit::Unitful.Mass, mfinal::Unitful.Mass, postagb_age::Unitful.Time, met, logL, logT; kwargs...)
    pn_I5007(pt, ustrip(u"Msun", minit), ustrip(u"Msun", mfinal), ustrip(u"Gyr", postagb_age), met, logL, logT; kwargs...)
end

absolute_magnitude_5007(I5007) = magnitude(max(I5007, 1e-15)) - 13.74

magnitude(I) = -2.5log10(I)
