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
function Delta_wrapper(logOmega, price_uk, zOC, mu, delta, tau, gammaK, xi)

    muBar = 1.0
    sigmaBar = 1.0
    lambda = 1.0
    R = 1.0

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

using Plots

function compute_capital_market(price_uk::Vector{T}, zOC::Vector{T}, data::IndustryData, params::Parameters) where {T <: Real}

    tau = compute_advalorem_tax(data)
    muBar = compute_muBar(params)
    industry_names = data.tax.industry

    capital_liquidated = Vector{T}(undef, params.constants.number_of_industries)
    capital_demand = Vector{T}(undef, params.constants.number_of_industries)
    free_cash_flow = Vector{T}(undef, params.constants.number_of_industries)

    for i in 1:params.constants.number_of_industries
    #for i in 10:10
        # @show price_uk[i]
        # @show zOC[i]
        # @show TzOC(zOC[i])
        # @show params.production.shock_mean[i] # mu
        # @show muBar[i] # muBar
        # @show params.production.shock_stdev[i] #sigma
        # @show data.shock_stdev.val[i] # sigmaBar
        # @show data.depreciation.val[i] # delta
        # @show tau[i] # tau
        # @show params.production.capital[i] # gammaK
        # @show params.constants.elasticities.production.substitution # xi
        # @show params.constants.loss_given_default # lambda
        # @show params.constants.interest_rate # R
        # @show data.capital.next_year[i] # k1

        @show "### industry ", industry_names[i]

        df = filter(row -> row.SIC16 == industry_names[i], data.assets_liabilities.current_year)
        out = compute_capital_market(price_uk[i],
                                     zOC[i],
                                     params.production.shock_mean[i], # mu
                                     muBar[i], # muBar
                                     params.production.shock_stdev[i], #sigma
                                     data.shock_stdev.val[i], # sigmaBar
                                     data.depreciation.val[i], # delta
                                     tau[i], # tau
                                     params.production.capital[i], # gammaK
                                     params.constants.elasticities.production.substitution, # xi
                                     params.constants.loss_given_default, # lambda
                                     params.constants.interest_rate, # R
                                     data.capital.next_year[i], # k1
                                     df.Assets,
                                     df.Ratio,
                                     Delta)

        capital_liquidated[i] = out[1]
        capital_demand[i] = out[2]
        free_cash_flow[i] = out[3]
    end

    return capital_liquidated, capital_demand, free_cash_flow

end

function compute_capital_market(price_uk::T, zOC::T, mu::T, muBar::T, sigma::T, sigmaBar::T,
                                delta::T, tau::T, gammaK::T, xi::T, lambda::T, R::T, k1::T,
                                assets::Vector{T}, liabilities::Vector{T}, DeltaFun::Function) where {T <: Real}

    grid_size = 100
    grid = range(muBar - 4 * sigmaBar, muBar + 4 * sigmaBar, grid_size)

    Bval = B(price_uk, mu, zOC, delta, tau, gammaK, xi)
    bval = b(price_uk, mu, zOC, tau, gammaK, xi)

    # @show Bval
    # @show bval

    fun(logOmega) = DeltaFun(logOmega, Bval, bval, mu, muBar, sigmaBar, lambda, R)
    #fun(logOmega) = cos(2 * logOmega) + sin(4 * logOmega) + 0.3 * abs(logOmega)

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
        push!(iMax, grid_size)
        push!(DeltaMax, fun(grid[end]))
    end

    @show iMin, DeltaMin
    @show iMax, DeltaMax
    #@show fun.(grid)

    # Keep track of the highest max found so far
    global_max = -Inf

    P = plot(grid, fun.(grid))
    plot!(P, grid[iMin], DeltaMin, seriestype=:scatter)
    plot!(P, grid[iMax], DeltaMax, seriestype=:scatter)

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

        @show i, interval

        # Find liabilities data points and indices that are within values of fun in the found interval.
        iL = findall(x -> x >= interval[1] && x < interval[2], liabilities)
        L = liabilities[iL]

        @show length(L)

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

    zeta = (logOmegaBar .- muBar) ./ sigmaBar

    F = cdf.(Normal(), zeta)

    # Find all liabilities above the maximum of Delta
    iGlobalMax = findall(x -> x > maximum(fun.(grid)), liabilities)

    append!(nonzero_indices, iGlobalMax)
    append!(F, ones(length(iGlobalMax)))

    plot!(P, logOmegaBar, liabilities[nonzero_indices], seriestype=:scatter)
    display(P)

    @show mean(assets[nonzero_indices]), mean(F)
    @show dot(assets[nonzero_indices], F)

    capital_liquidated = (1 - lambda) * (1 - delta) * dot(assets[nonzero_indices], F)
    capital_demand = k1 - (1 - delta) * dot(assets[nonzero_indices], (1 .- F))

    # TODO: Implement FCFTerms from
    # https://github.com/UCL/Supergrassi/blob/b807e57cd49da2c5f907accf3863dfa312bd39f5/code/matlab/macro_v2/finance/CapitalMarket.m#L322-L324
    FCFTerm2 = 0.0
    FCFTerm3 = 0.0
    FCFTerm4 = 0.0
    free_cash_flow = -dot(assets[nonzero_indices] .* (1 .- liabilities[nonzero_indices]), 1 .- F) + bval * (FCFTerm2 - FCFTerm3) - FCFTerm4

    return capital_liquidated, capital_demand, free_cash_flow

end

# function compute_dividends(FCF::Vector{T}, params::Parameters) where {T <: Real}

#     return FCF + params.D1 - params.D0 - q1 * (params.k1 - params.k0);

# end

# Random.seed!(1238)

# price_uk = rand()
# mu = rand()
# zOC = rand()
# delta = rand()
# tau = rand()
# gammaK = rand()
# chi0 = rand()
# xi = rand()
# q0 = rand()
# k0 = rand()
# k1 = rand()

# # L = 1.0
# σ = 1.4

# muBar = 1.0
# sigmaBar = 1.0
# lambda = 1.0
# R = 1.0

# assets = randn(10000) .* 10
# liabilities = randn(10000) .* 10

# KL, KD, FCF = compute_capital_market(price_uk, mu, muBar, σ, assets, liabilities, Delta)

# Δ(Ω) = Delta_wrapper(Ω, price_uk, zOC, mu, delta, tau, gammaK, chi0, xi, q0)
# Δres(Ω) = residual(Ω, L, Δ)
# logOmega = find_zero(Δres, (mu - 4 * σ, mu + 4 * σ), Bisection(), verbose=true)

# grad = gradient(ForwardWithPrimal, Delta_wrapper, logOmega, price_uk, zOC, Const(mu),
#                 Const(delta), Const(tau), Const(gammaK), Const(chi0), Const(xi), Const(q0))

# # eqns H.1 and H.2
# ∂logω_∂pdj = - grad.derivs[1] / grad.derivs[2]
# ∂logω_∂zOC = - grad.derivs[1] / grad.derivs[3]
