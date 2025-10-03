using Supergrassi
using Distributions
using Enzyme
using Roots
using Peaks
using Random
using LinearAlgebra

"""
Wrapper around Delta that exposes logOmega, price_uk and zOC as the arguments for computation
of derivatives H.3 - H.5 (on the RHS of H.1 and H.2)
"""
function Delta_wrapper(logOmega, price_uk, zOC, mu, muBar, sigmaBar, delta, tau, gammaK, xi, lambda, R)

    Bval = B(price_uk, mu, zOC, delta, tau, gammaK, xi)
    bval = b(price_uk, mu, zOC, tau, gammaK, xi)

    return Delta(logOmega, Bval, bval, mu, muBar, sigmaBar, lambda, R)

end

"""
Function Delta as a function of logOmega, B and b as written in the paper (C. 47).
"""
function Delta(logOmega::T, B::T, b::T, μ::T, muBar::T, sigmaBar::T, λ::T, R::T) where {T <: Real}

    ζ = (logOmega - muBar) / sigmaBar
    F = cdf(Normal(), ζ)

    ω0 = -B / b
    if (ω0 > 0)
        ζ0 = (log(ω0) - muBar) / sigmaBar
        F0 = cdf(Normal(), (log(ω0) - muBar) / sigmaBar)
        ω̂ = μ * (cdf(Normal(), ζ - sigmaBar) - cdf(Normal(), ζ0 - sigmaBar))
    else
        F0 = 0.0
        ω̂ = μ * cdf(Normal(), ζ - sigmaBar)
    end

    x = B * (1.0 - F0 - λ * (F - F0)) / R
    y = b * (exp(logOmega) * (1.0 - F) + (1.0 - λ) * ω̂) / R

    return x + y

end

"""
Second part of C.47

B(iNZ) = 1 - parms.delta(iNZ) + pd(iNZ).*(1-parms.tau(iNZ)).*parms.chi0(iNZ)/parms.q0 ...
      - rk(iNZ).*ROCTheta(iNZ)./(1-ROCTheta(iNZ));
"""
function B(price_uk, μ, zOC, δ, τ, γK, χ0, ξ, q0)

    return 1.0 - δ + (1.0 - τ) * price_uk * χ0 / q0 - rk(price_uk, μ, zOC, τ, γK, ξ, q0) * TzOC(zOC) / (1.0 - TzOC(zOC))

end

"""
Second part of C.47

b(iNZ) = rk(iNZ)./(parms.mu(iNZ).*(1-ROCTheta(iNZ)));
"""
function b(price_uk, μ, zOC, τ, γK, ξ, q0)

    return rk(price_uk, μ, zOC, τ, γK, ξ, q0) / (μ * (1 - TzOC(zOC)))

end

"""
As far as I can tell, in the Matlab code chi0 is always 0, q0 is always 1.
Therefore this method assumes chi0 = 0, q0 = 1 and simplifies the calculation.
"""
function B(price_uk, μ, zOC, δ, τ, γK, ξ)

    return 1.0 - δ - rk(price_uk, μ, zOC, τ, γK, ξ, 1) * TzOC(zOC) / (1.0 - TzOC(zOC))

end

"""
This method assumes q0 = 1 and simplifies the calculation
"""
function b(price_uk, μ, zOC, τ, γK, ξ)

    return rk(price_uk, μ, zOC, τ, γK, ξ, 1) / (μ * (1 - TzOC(zOC)))

end


"""
Above C.1
"""
function TzOC(zOC)

    return exp(zOC) / (1 + exp(zOC))

end

"""
C. 48
"""
function rk(price_uk, μ, zOC, τ, γK, ξ, q0)

    return (1 - τ) * price_uk * μ * γK ^ (1 / (ξ - 1)) * (1 - TzOC(zOC)) ^ (1 / (1 - ξ)) / q0

end

function residual(logOmega::T, L::T, fun::Function) where {T<:Real}

    return L - fun(logOmega)

end

"""
C. 44
"""
function G(price_uk::T, zOC::T, mu::T, gammaK::T, delta::T,
           tau::T, logOmega::T, chi0::T, xi::T, q0::T, k0::T, muBar::T, sigmaBar::T) where {T <: Real}

    Bval = B(price_uk, mu, zOC, delta, tau, gammaK, chi0, xi , q0)
    bval = b(price_uk, mu, zOC, tau, gammaK, xi, q0)

    dist = Normal()

    if any(Bval .< 0)

        ζ0 = logOmega - muBar / sigmaBar

        G = q0 * k0 * (
            B * cdf(dist, ζ0)
            + b * mu * cdf(dist, ζ0 - sigmaBar)
        )

    else

        G = zeros(length(price_uk))

    end

    return G

end

function omegaV(zOC::T, μ::T, χ1::T, γK::T, ξ::T) where {T <: Real}

    ωV = μ - χ1 * ((1 - TzOC(zOC)) / γK ^ (1/ξ)) ^ (ξ/(ξ-1))
    return max(ωV, 0.0)

end

# using Plots

function compute_capital_market(price_uk::Vector{T}, zOC::Vector{T}, data::IndustryData, params::Parameters) where {T <: Real}

    tau = compute_advalorem_tax(data)
    muBar = compute_muBar(params)
    industry_names = data.tax.industry

    capital_liquidated = Vector{T}(undef, params.constants.number_of_industries)
    capital_demand = Vector{T}(undef, params.constants.number_of_industries)
    free_cash_flow = Vector{T}(undef, params.constants.number_of_industries)
    logOmegaBar = Vector{T}(undef, params.constants.number_of_industries)

    # plots = []

    for i in 1:params.constants.number_of_industries

        @debug "### industry ", industry_names[i]

        df = filter(row -> row.SIC16 == industry_names[i], data.assets_liabilities.current_year)

        grid = compute_grid(100, muBar[i], data.shock_stdev.val[i])

        Bval = B(price_uk[i], params.production.shock_mean[i], zOC[i], data.depreciation.val[i], tau[i],
                 params.production.capital[i], params.constants.elasticities.production.substitution)
        bval = b(price_uk[i], params.production.shock_mean[i], zOC[i], tau[i],
                 params.production.capital[i], params.constants.elasticities.production.substitution)

        fun(logOmega) = Delta(logOmega, Bval, bval, params.production.shock_mean[i], muBar[i],
                              data.shock_stdev.val[i], params.constants.loss_given_default,
                              params.constants.interest_rate)

        #fun(logOmega) = cos(2 * logOmega) + sin(4 * logOmega) + 0.3 * abs(logOmega)

        ind, logOmegaBar, P = compute_logOmegaBar(bval, Bval, grid, df.Ratio, fun)

        grad = logOmegaBar_gradients(logOmegaBar, price_uk[i], zOC[i],
                                     params.production.shock_mean[i], # mu
                                     data.depreciation.val[i], # delta
                                     tau[i],
                                     params.production.capital[i], # gammaK
                                     params.constants.elasticities.production.substitution)

        # push!(plots, P)

        KL, KD, FCF = capital_market_terms(price_uk[i],
                                           zOC[i],
                                           params.production.shock_mean[i], # mu
                                           muBar[i], # muBar
                                           data.shock_stdev.val[i], # sigmaBar
                                           data.depreciation.val[i], # delta
                                           tau[i], # tau
                                           params.production.capital[i], # gammaK
                                           params.constants.elasticities.production.substitution, # xi
                                           params.constants.loss_given_default, # lambda
                                           params.constants.interest_rate, # R
                                           data.capital.current_year[i], # k0
                                           data.capital.next_year[i], # k1
                                           abs(data.regional.delta_v.agg[i]) / data.capital.next_year[i], #chi1
                                           df.Assets, df.Ratio, ind, logOmegaBar, fun, grid)

        capital_liquidated[i] = KL
        capital_demand[i] = KD
        free_cash_flow[i] = FCF

    end

    # P = plot((plots[i] for i in 1:16)...; layout=16)
    # display(P)

    return capital_liquidated, capital_demand, free_cash_flow

end

function logOmegaBar_gradients(logOmegaBar, price_uk, zOC, mu, delta, tau, gammaK, xi)

    grad = gradient(ForwardWithPrimal, Delta_wrapper, logOmegaBar, price_uk, zOC, Const(mu),
                    Const(delta), Const(tau), Const(gammaK), Const(xi))

    # eqns H.1 and H.2
    ∂logω_∂pdj = - grad.derivs[2] / grad.derivs[1]
    ∂logω_∂zOC = - grad.derivs[3] / grad.derivs[1]

    return ∂logω_∂pdj, ∂logω_∂zOC

end

function compute_grid(grid_size::Int, muBar::T, sigmaBar::T) where {T <: Real}

    return range(muBar - 4 * sigmaBar, muBar + 4 * sigmaBar, grid_size)

end


function compute_logOmegaBar(bval::T, Bval::T, grid, liabilities::Vector{T}, fun::Function) where {T <: Real}

    # Find local maxima and minima of fun on grid
    iMin, DeltaMin = findminima(fun.(grid))
    iMax, DeltaMax = findmaxima(fun.(grid))

    # Add first grid point to minima if it is a minima (we don't care if it's a maxima)
    if fun(grid[1]) < fun(grid[2])
        pushfirst!(iMin, 1)
        pushfirst!(DeltaMin, fun(grid[1]))
    end

    # Add last grid point to maxima if it is a maxima (we don't care if it's a minima)
    if fun(grid[end]) > fun(grid[end-1])
        push!(iMax, length(grid))
        push!(DeltaMax, fun(grid[end]))
    end

    @debug iMin, DeltaMin
    @debug iMax, DeltaMax
    # @debug fun.(grid)

    # Keep track of the highest max found so far
    global_max = -Inf

    # P = plot(grid, fun.(grid))
    # plot!(P, grid[iMin], DeltaMin, seriestype=:scatter)
    # plot!(P, grid[iMax], DeltaMax, seriestype=:scatter)

    logOmegaBar = Vector{T}(undef, 0)
    nonzero_indices = Vector{Int}(undef, 0)

    # Loop through each interval from local minima to local maxima where we know
    # fun is increasing. Starting from the first local min, find the interval
    # above the global max to the next local max where fun is increasing.
    for i = 1:length(iMin)

        if (DeltaMin[i] > global_max)
            interval = (DeltaMin[i], DeltaMax[i])
            global_max = DeltaMax[i]
        elseif (DeltaMax[i] > global_max)
            interval = (global_max, DeltaMax[i])
            global_max = DeltaMax[i]
        else
            interval = (0.0,0.0)
        end

        @info "Interval", i, ": ", interval

        # Find liabilities data points and indices that are within values of fun in the found interval.
        iL = findall(x -> x >= interval[1] && x < interval[2], liabilities)
        L = liabilities[iL]

        @info "Length(L): ",length(L)

        if (length(L) > 0)
            for il = 1:length(L)
                # Define residual function that returns the difference between fun(Ω) and L
                Δres(Ω) = residual(Ω, L[il], fun)
                # Find Ω such that L - fun(Ω) == 0 and store the solution and the corresponding index of liabilities
                push!(logOmegaBar, find_zero(Δres, (grid[iMin[i]], grid[iMax[i]]), Bisection(), verbose=false))
                push!(nonzero_indices, iL[il])
            end
        else
            @info "no roots found in interval", i
        end
    end

    # plot!(P, logOmegaBar, liabilities[nonzero_indices], seriestype=:scatter)

    return nonzero_indices, logOmegaBar #, P

end

function capital_market_terms(price_uk::T, zOC::T, mu::T, muBar::T, sigmaBar::T, delta::T,
                              tau::T, gammaK::T, xi::T, lambda::T, R::T, k0::T, k1::T, chi1::T,
                              assets::Vector{T}, liabilities::Vector{T},
                              nonzero_indices::Vector{Int}, logOmegaBar::Vector{T}, fun::Function, grid) where {T <: Real}

    bval = b(price_uk, mu, zOC, tau, gammaK, xi)

    ζ = (logOmegaBar .- muBar) ./ sigmaBar
    ζV = (log(omegaV(zOC, mu, chi1, gammaK, xi)) - muBar) / sigmaBar

    F = cdf.(Normal(), ζ)
    FV = cdf(Normal(), ζV)
    omegaHat2 = mu .* cdf(Normal(), ζ .- sigmaBar) ./ F

    # Find all liabilities above the maximum of Delta
    iGlobalMax = findall(x -> x > maximum(fun.(grid)), liabilities)

    @debug mean(assets[nonzero_indices]), mean(F)
    @debug dot(assets[nonzero_indices], F)

    iωV = omegaV(zOC, mu, chi1, gammaK, xi) .<= exp.(logOmegaBar)
    F2 = copy(F)
    F2[iωV] .= FV
    omegaV(zOC, mu, chi1, gammaK, xi) > 0.0 ? omegaHat2[iωV] .= mu * cdf(Normal(), ζV - sigmaBar / FV) : omegaHat2[iωV] .= 0.0

    append!(nonzero_indices, iGlobalMax)
    append!(F, ones(length(iGlobalMax)))
    append!(F2, ones(length(iGlobalMax)))
    append!(omegaHat2, ones(length(iGlobalMax)) * mu)
    append!(logOmegaBar, ones(length(iGlobalMax)) * Inf)

    FCFTerm2 = mu * (k0 - dot(assets[nonzero_indices], F - F2))
    FCFTerm3 = dot(assets[nonzero_indices], logOmegaBar .* (1 .- F) + omegaHat2 .* F2)
    FCFTerm4 = (1 - tau) * price_uk * chi1 * dot(assets[nonzero_indices], 1 .- (F .- FV))

    # @show FCFTerm2, FCFTerm3, FCFTerm4

    capital_liquidated = (1 - lambda) * (1 - delta) * dot(assets[nonzero_indices], F)
    capital_demand = k1 - (1 - delta) * dot(assets[nonzero_indices], (1 .- F))
    free_cash_flow = -dot(assets[nonzero_indices] .* (1 .- liabilities[nonzero_indices]), 1 .- F) + bval * (FCFTerm2 - FCFTerm3) - FCFTerm4

    return capital_liquidated, capital_demand, free_cash_flow

end

# TODO: Fix this function with correct parameters
# function compute_dividends(FCF::Vector{T}, params::Parameters) where {T <: Real}

#     return FCF + params.D1 - params.D0 - q1 * (params.k1 - params.k0);

# end
