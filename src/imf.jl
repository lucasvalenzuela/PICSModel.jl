function pn_number(imf::AbstractIMF, ssp_minit::Number, minit_lower::Number, minit_upper::Number; rng=nothing)
    if !isfinite(minit_lower) || !isfinite(minit_upper) || iszero(minit_lower) || iszero(minit_upper)
        return 0
    end

    n = nstars_in_massinterval(imf, ssp_minit, (minit_lower, minit_upper))

    r = isnothing(rng) ? rand() : rand(rng)
    Δn = r > n - floor(n) ? 0 : 1
    n = floor(Int, n) + Δn

    return n
end
