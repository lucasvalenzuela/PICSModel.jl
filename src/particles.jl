function pn_inds(pn_numbers::AbstractVector{T}) where {T<:Integer}
    n = sum(pn_numbers)
    inds = Vector{T}(undef, n)

    iinds = 1
    @inbounds for (i, number) in enumerate(pn_numbers)
        for _ in 1:number
            inds[iinds] = i
            iinds += 1
        end
    end

    return inds
end


_pn_lifetime(::T) where {T<:Real} = float(T)(3e-5)
_pn_lifetime(::Unitful.Time) = 30_000u"yr" |> u"Gyr"
_pn_initial_mass_upper(l, ssp_age, ssp_met) = pn_initial_mass(l, ssp_age - _pn_lifetime(ssp_age), ssp_met)

function pn_population(method::PNMethod, ssp_minit::AbstractVector, ssp_age::AbstractVector, ssp_met::AbstractVector; rng=nothing)
    @assert length(ssp_minit) == length(ssp_age) == length(ssp_met)

    minit = pn_initial_mass.(method.lifetime, ssp_age, ssp_met)
    minit_upper = _pn_initial_mass_upper.(method.lifetime, ssp_age, ssp_met)
    pn_numbers = pn_number.(method.imf, ssp_minit, minit, minit_upper; rng)
    ind = pn_inds(pn_numbers)

    n_pne = length(ind)
    minit = minit[ind]
    met = ssp_met[ind]

    postagb_age = _pn_lifetime(first(ssp_age)) .* (isnothing(rng) ? rand(n_pne) : rand(rng, n_pne))

    mfinal = final_mass.(method.initial_final_mass, minit, met)
    if eltype(minit) <: Quantity
        logL = pn_luminosity.(method.postagb_tracks, minit, mfinal, postagb_age, met)
        logT = pn_temperature.(method.postagb_tracks, minit, mfinal, postagb_age, met)
    else
        logL, logT = broadcast_unzip(pn_luminosity_temperature, method.postagb_tracks, minit, mfinal, postagb_age, met)
    end

    # throw away any particles that are clearly not PNe
    mask = isfinite.(logL)
    ind = ind[mask]
    minit = minit[mask]
    mfinal = mfinal[mask]
    met = met[mask]
    postagb_age = postagb_age[mask]
    logL = logL[mask]
    logT = logT[mask]

    I5007 = @. pn_I5007(method.i5007, minit, mfinal, postagb_age, met, logL, logT; rng) * pn_absorbing_factor(method.absorbing_factor, minit, mfinal, postagb_age, met, logL, logT; rng)
    M5007 = absolute_magnitude_5007.(I5007)

    return (; ind, met, postagb_age, minit, mfinal, logL, logT, M5007)
end

function pn_population(stars::Particles, method::PNMethod; pntype=:pn, initialmassprop=:iM, ageprop=:age, metprop=:met, m5007prop=:m5007, pn_initialmassprop=:pnimass, pn_finalmassprop=:pnfmass, pn_postagb_ageprop=:postagb_age, logLprop=:logL, logTprop=:logT, rng=nothing)
    pn = Particles(pntype)

    ssp_minit = stars[initialmassprop]
    ssp_age = stars[ageprop]
    ssp_met = stars[metprop]

    vals = pn_population(method, ssp_minit, ssp_age, ssp_met; rng)

    # save PNe properties to PN Particle object
    pn[pn_postagb_ageprop] = vals.postagb_age
    pn[pn_initialmassprop] = vals.minit
    pn[pn_finalmassprop] = vals.mfinal
    pn[logLprop] = vals.logL
    pn[logTprop] = vals.logT
    pn[m5007prop] = vals.M5007

    # copy over stellar props to pn particles (pos, vel, etc.)
    ind = vals.ind # indices to index into stars
    pn[metprop] = vals.met # make use of metallicities already having been indexed into
    for prop in keys(stars)
        prop === metprop && continue
        pn[prop] = CosmoParticles._applyind(stars[prop], ind)
    end

    return pn
end
