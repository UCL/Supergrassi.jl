using Distributions
using Enzyme
using Roots
using Peaks
using Random

"""
Wrapper around Delta that exposes logOmega, price_uk and zOC as the arguments for computation
of derivatives H.3 - H.5 (on the RHS of H.1 and H.2)
"""
function Delta_wrapper(logOmega, price_uk, zOC, mu, delta, tau, gammaK, chi0, xi, q0)

    muBar = 1.0
    sigmaBar = 1.0
    lambda = 1.0
    R = 1.0

    Bval = B(price_uk, mu, zOC, delta, tau, gammaK, chi0, xi , q0)
    bval = b(price_uk, mu, zOC, tau, gammaK, xi, q0)

    return Delta(logOmega, Bval, bval, mu, muBar, sigmaBar, lambda, R)

end

"""
Function Delta as a function of logOmega, B and b as written in the paper (C. 47).
"""
function Delta(logOmega::T, B::T = 0.0, b::T = 1.0, μ::T = 1.0, muBar::T = 1.0, sigmaBar::T = 1.0, λ::T = 1.0, R::T = 1.0) where {T <: Real}

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

    return rk(price_uk, μ, zOC, τ, γK, ξ, q0) / μ * (1 - TzOC(zOC))

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

function compute_capital_market(price_uk::T, mu::T, muBar::T, sigma::T, liabilities::Vector{T},
                                fun::Function) where {T <: Real}

    grid_size = 100
    grid = range(muBar - 4 * sigma, muBar + 4 * sigma, grid_size)

    Bval = B(price_uk, mu, zOC, delta, tau, gammaK, chi0, xi , q0)
    bval = b(price_uk, mu, zOC, tau, gammaK, xi, q0)

    #DeltaFun(logOmega) = fun(logOmega, Bval, bval, mu, muBar, sigmaBar, lambda, R)
    DeltaFun(logOmega) = cos(2 * logOmega) + sin(4 * logOmega) + 0.3 * abs(logOmega)

    # Find local maxima and minima of DeltaFun on grid
    iMin, DeltaMin = findminima(DeltaFun.(grid))
    iMax, DeltaMax = findmaxima(DeltaFun.(grid))

    # Add first grid point to minima (we don't care if it's a maxima)
    if DeltaFun(grid[1]) < DeltaFun(grid[2])
        pushfirst!(iMin, 1)
        pushfirst!(DeltaMin, DeltaFun(grid[1]))
    end

    # Add last grid point to maxima if it is a maxima (we don't care if it's a minima)
    if DeltaFun(grid[end]) > DeltaFun(grid[end-1])
        push!(iMax, grid_size)
        push!(DeltaMax, DeltaFun(grid[end]))
    end

    @show iMin, DeltaMin
    @show iMax, DeltaMax
    #@show DeltaFun.(grid)

    # Keep track of the highest max found so far
    global_max = -Inf

    P = plot(grid, DeltaFun.(grid))
    plot!(P, grid[iMin], DeltaMin, seriestype=:scatter)
    plot!(P, grid[iMax], DeltaMax, seriestype=:scatter)

    # Loop through each interval from local minima to local maxima where we know
    # DeltaFun is increasing. Starting from the first local min, find the interval
    # above the global max to the next local max where DeltaFun is increasing.
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

        # Find liabilities data points that are within values of DeltaFun in the found interval.
        L = liabilities[liabilities .>= interval[1] .&& liabilities .<= interval[2]]
        logOmegaBar = zeros(length(L))
        @show length(L)

        if (length(L) > 0)
            for il = 1:length(L)
                # Define residual function that returns the difference between DeltaFun(Ω) and L
                Δres(Ω) = residual(Ω, L[il], DeltaFun)
                # Find Ω such that the difference is 0
                logOmegaBar[il] = find_zero(Δres, (grid[iMin[i]], grid[iMax[i]]), Bisection(), verbose=false)
            end
            plot!(P, logOmegaBar, L, seriestype=:scatter)
            # @show logOmegaBar
        else
            @info "no roots found in interval", i
        end

        # TODO: Accumulate logOmegaBar into one vector.
        # Do we need to keep track of the mapping to liabilities?
        
    end

    display(P)

    # return KL, KD, FCF
    # return logOmegaBar

end

# function compute_dividends(FCF::Vector{T}, params::Parameters) where {T <: Real}

#     return FCF + params.D1 - params.D0 - q1 * (params.k1 - params.k0);

# end

Random.seed!(1238)

price_uk = rand()
mu = rand()
zOC = rand()
delta = rand()
tau = rand()
gammaK = rand()
chi0 = rand()
xi = rand()
q0 = rand()
k0 = rand()

# L = 1.0
σ = 1.4

muBar = 1.0
sigmaBar = 1.0
lambda = 1.0
R = 1.0

liabilities = randn(10000) .* 10

compute_capital_market(price_uk, mu, muBar, σ, liabilities, Delta)

# Δ(Ω) = Delta_wrapper(Ω, price_uk, zOC, mu, delta, tau, gammaK, chi0, xi, q0)
# Δres(Ω) = residual(Ω, L, Δ)
# logOmega = find_zero(Δres, (mu - 4 * σ, mu + 4 * σ), Bisection(), verbose=true)

# grad = gradient(ForwardWithPrimal, Delta_wrapper, logOmega, price_uk, zOC, Const(mu),
#                 Const(delta), Const(tau), Const(gammaK), Const(chi0), Const(xi), Const(q0))

# # eqns H.1 and H.2
# ∂logω_∂pdj = - grad.derivs[1] / grad.derivs[2]
# ∂logω_∂zOC = - grad.derivs[1] / grad.derivs[3]
