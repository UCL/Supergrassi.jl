using Distributions
using Enzyme

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
function Delta(logOmega::T, B::T = 0.0, b::T = 0.0, μ::T = 1.0, muBar::T = 1.0, sigmaBar::T = 1.0, λ::T = 1.0, R::T = 1.0) where {T <: Real}

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


logOmega = 1.0
price_uk = rand()
mu = rand()
zOC = rand()
delta = rand()
tau = rand()
gammaK = rand()
chi0 = rand()
xi = rand()
q0 = rand()

grad = gradient(ForwardWithPrimal, Delta_wrapper, logOmega, price_uk, zOC, Const(mu),
                Const(delta), Const(tau), Const(gammaK), Const(chi0), Const(xi), Const(q0))

# eqns H.1 and H.2
∂logω_∂pdj = - grad.derivs[1] / grad.derivs[2]
∂logω_∂zOC = - grad.derivs[1] / grad.derivs[3]
