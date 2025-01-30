struct LifetimeFunction_P93_R81 <: AbstractLifetimeFunction end


@doc raw"""
    ms_lifetime(::LifetimeFunction_P93_R81, m::Number, met=nothing)

Returns the main sequence lifetime of a star of mass `m` (in M⊙) in Gyr.

Lifetime function from [Padovani & Matteucci (1993)](https://ui.adsabs.harvard.edu/abs/1993ApJ...416...26P), Eq. 3 and 4. Note that the -9 at the end of Eq. 4 is supposed to be in the exponent.

```math
\tau_\mathrm{MS}(m) =
\begin{cases}
    1.2m^{-1.85} + 0.003 \quad\text{if}\ \ m > 6.6\;\mathrm{M}_\odot, \\
    10^{ (1.338 - \sqrt{ 1.79 - 0.2232 (7.764 - \log{m}) }) / 0.1116 - 9 } \quad \text{otherwise}. \\
\end{cases}
```
"""
function ms_lifetime(::LifetimeFunction_P93_R81, m::Real, met=nothing)
    if m ≤ 6.6
        return 10^((1.338 - sqrt(1.79 - 0.2232(7.764 - log10(m)))) / 0.1116 - 9)
    else
        return 1.2m^-1.85 + 0.003
    end
end
ms_lifetime(l::LifetimeFunction_P93_R81, m::Unitful.Mass, met=nothing) = ms_lifetime(l, ustrip(u"Msun", m)) * u"Gyr"


"""
    ms_dyingmass(::LifetimeFunction_P93_R81, τ::Number, met=nothing)

Returns the mass of stars reaching the end of the main sequence after the time `τ` (in Gyr) in M⊙.

Inverse of [`ms_lifetime`](@ref).
"""
function ms_dyingmass(::LifetimeFunction_P93_R81, τ::Real, met=nothing)
    if τ > 0.03976531865906481 # = lifetime(6.6)
        return exp10(7.764 - (1.79 - (1.338 - 0.1116 * (9 + log10(τ)))^2) / 0.2232)
    elseif τ > 0.003
        return ((τ - 0.003) / 1.2)^(-1 / 1.85)
    else
        return Inf
    end
end
ms_dyingmass(l::LifetimeFunction_P93_R81, τ::Unitful.Time, met=nothing) = ms_dyingmass(l, ustrip(u"Gyr", τ)) * u"Msun"


@doc raw"""
    ms_dm_dτ(::LifetimeFunction_P93_R81, τ::Number[, m::Number])

Returns the time derivative, ``dm/dτ``, of [`ms_dyingmass`](@ref) in M⊙/Gyr.

If the dyingmass at time `τ` is already known, it can be passed as a second argument for better performance.
"""
function ms_dm_dτ(l::LifetimeFunction_P93_R81, τ::Real, m::Real=ms_dyingmass(l, τ))
    if τ > 0.03976531865906481 # = lifetime(6.6)
        return -m / τ * (1.338 - 0.1116 * (9 + log10(τ)))
    elseif τ > 0.003
        return -1 / 1.85 * m / (τ - 0.003)
    else
        return -Inf
    end
end
ms_dm_dτ(l::LifetimeFunction_P93_R81, τ::Unitful.Time) = ms_dm_dτ(l, ustrip(u"Gyr", τ)) * u"Msun/Gyr"
function ms_dm_dτ(l::LifetimeFunction_P93_R81, τ::Unitful.Time, m::Unitful.Mass)
    ms_dm_dτ(l, ustrip(u"Gyr", τ), ustrip(u"Msun", m)) * u"Msun/Gyr"
end




@doc raw"""
    pms_lifetime(::LifetimeFunction_P93_R81, mass::Number, met=nothing)

Returns the post-main sequence lifetime of a star of mass `m` (in M⊙) in Gyr.

Lifetime function from [Renzini (1981)](https://ui.adsabs.harvard.edu/abs/1981AnPh....6...87R) and
[Renzini & Buzzoni (1986)](https://ui.adsabs.harvard.edu/abs/1986ASSL..122..195R).

```math
\tau_\mathrm{PMS}(m) = 1.66 m^{-2.72}
```
"""
pms_lifetime(::LifetimeFunction_P93_R81, mass::Real, met=nothing) = 1.66mass^-2.72
pms_lifetime(l::LifetimeFunction_P93_R81, mass::Unitful.Mass, met=nothing) = pms_lifetime(l, ustrip(u"Msun", mass)) * u"Gyr"


function ms_pms_lifetime(l::LifetimeFunction_P93_R81, mass::Number, met=nothing)
    ms_lifetime(l, mass, met) + pms_lifetime(l, mass, met)
end

# function used for root solving of the problem t(M) = t_MS(M) + t_PMS(M)
# with the parameter p representing the age
_problem_ms_pms_lifetime_p93_r81(mass, p) = ms_pms_lifetime(LifetimeFunction_P93_R81(), mass) - p


function pn_initial_mass(l::LifetimeFunction_P93_R81, ssp_age::Real, ssp_met::Real)
    p = ZeroProblem(_problem_ms_pms_lifetime_p93_r81, ms_dyingmass(l, ssp_age))
    return solve(p, Order1(); p=ssp_age)
end
