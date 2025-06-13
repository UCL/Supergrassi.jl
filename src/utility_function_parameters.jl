"""
  Compute an utility function parameter for a single region
"""
function parameter_by_region(elasticity::T, quantity_region::T, logP_region::T, logP::T) where {T <: Real}

    return quantity_region == 0.0 ? 0.0 : quantity_region * exp((elasticity - 1) * (logP_region - logP))

end

"""
  Compute log of an utility function parameter for a single region
"""
function log_parameter_by_region(elasticity::T, quantity_region::T, logP_region::T, logP::T) where {T <: Real}

    return quantity_region == 0.0 ? 0.0 : log(quantity_region) + (elasticity - 1) * (logP_region - logP)

end

"""
  Compute utility function parameters by region considered in the model (uk, eu, rest of world).
  Parameters of this type appear in multiple utility functions in the paper, and are annotated
  (at least) α, β1, β2, γ and ρ.

  This is refactored from the Matlab code in e.g. ComputeTheta lines 62-64

  # Arguments
   - elasticity [ϵ, χ1, χ2, η, ξ]
   - log_price_{uk, eu, w} [log(p)]
   - quantity_{uk, eu, w} [f, x1, x2, I, m]

  # Outputs
  - weight_{uk, eu, w}

"""
function parameters_by_region(elasticity::T,
                              log_price_uk::T,log_price_eu::T,log_price_world::T,
                              quantity_uk::T,quantity_eu::T,quantity_world::T) where {T <: Real}

    logP = log_price_index(elasticity, log_price_uk, log_price_eu, log_price_world, quantity_uk, quantity_eu, quantity_world)

    weight_uk = parameter_by_region(elasticity, quantity_uk, log_price_uk, logP)
    weight_eu = parameter_by_region(elasticity, quantity_eu, log_price_eu, logP)
    weight_world = parameter_by_region(elasticity, quantity_world, log_price_world, logP)

    return weight_uk, weight_eu, weight_world

end

"""
  Compute log of utility function parameters. See parameters_by_region.
"""
function log_parameters_by_region(elasticity::T,
                               log_price_uk::T,log_price_eu::T,log_price_world::T,
                               quantity_uk::T,quantity_eu::T,quantity_world::T) where {T <: Real}

    logP = log_price_index(elasticity, log_price_uk, log_price_eu, log_price_world, quantity_uk, quantity_eu, quantity_world)

    log_weight_uk = log_parameter_by_region(elasticity, quantity_uk, log_price_uk, logP)
    log_weight_eu = log_parameter_by_region(elasticity, quantity_eu, log_price_eu, logP)
    log_weight_world = log_parameter_by_region(elasticity, quantity_world, log_price_world, logP)

    return log_weight_uk, log_weight_eu, log_weight_world

end

"""
  Compute the total utility function parameter. Parameters of this type appear in multiple utility function
  in the paper, and are annotated (at least) α, β1, β2, γ and ρ.

  This is refactored from the Matlab code in e.g. ComputeTheta.m line 59

  # Arguments
  - log_price_index: price index computed by log_price_index()
  - quantity: total quantity [f, x1, x2, I, m]
  - elasticity: [ϵ, χ1, χ2, η, ξ]

  # Outputs
  - parameters
"""
function total_parameters(log_price_index::Vector{T}, quantity::Vector{T}, elasticity::T ) where {T <: Real}

    length(log_price_index) == length(quantity) || error()

    logPBar = log_total_price_index(elasticity, log_price_index, quantity)

    parameters = Vector{T}(undef, length(log_price_index))
    for i in axes(log_price_index,1)
        parameters[i] = weight_kernel(quantity[i], exp(log_price_index[i] - logPBar), elasticity)
    end
    return parameters

    #return [parameters[i] = weight_kernel(quantity[i], exp(log_price_index[i] - logPBar), elasticity) for i in eachindex(quantity, log_price_index)]

end

"""
  Compute the share of EU expenditure on UK exports, defined in section 2.2.1 of the main paper.
"""
function log_eu_expenditure_on_uk_exports(log_price_index::Vector{T}, quantity::Vector{T},
                                          Ex::T, ETilde::T, ePx::T, PTilde::T, elasticity::T,
                                          elasticity_tilde::T) where {T <: Real}

    length(log_price_index) == length(quantity) || error()
    logPBar = log_total_price_index(elasticity, log_price_index, quantity)

    return log_weight_kernel(Ex/ETilde, exp(logPBar) * ePx / PTilde, elasticity_tilde)

end

"""
  Compute the the consumer price index defined in equation 2.3 of the main paper as \bar{p}.
  Matlab code reference e.g. ComputeTheta.m line 49.
"""
function price_index(elasticity::T,
                     log_price_uk::T, log_price_eu::T, log_price_world::T,
                     demand_uk::T, demand_eu::T, demand_world::T) where {T <: Real}

    return (
        demand_uk ^ (1 / elasticity) * exp((elasticity - 1) * log_price_uk / elasticity) +
        demand_eu ^ (1 / elasticity) * exp((elasticity - 1) * log_price_eu / elasticity) +
        demand_world  ^ (1 / elasticity) * exp((elasticity - 1) * log_price_world  / elasticity)
    ) ^ ( elasticity / (elasticity - 1) )

end

function log_price_index(elasticity::T,
                         log_price_uk::T, log_price_eu::T, log_price_world::T,
                         demand_uk::T, demand_eu::T, demand_world::T) where {T <: Real}

    if(demand_uk == demand_eu == demand_world == 0.0) return 0.0 end

    return (elasticity / (elasticity - 1)) * log(
        demand_uk ^ (1 / elasticity) * exp((elasticity - 1) * log_price_uk / elasticity) +
        demand_eu ^ (1 /elasticity) * exp((elasticity - 1) * log_price_eu / elasticity) +
        demand_world  ^ (1 / elasticity) * exp((elasticity - 1) * log_price_world  / elasticity)
    )

end

"""
  Compute the consumer price index defined in equation 2.7 of the main paper as \bar{P}.
  Matlab code reference e.g. ComputeTheta.m line 58.
"""
function log_total_price_index(elasticity::T, log_price_index::Vector{T}, quantity::Vector{T}) where {T <: Real}

    s = sum_kernel(quantity, log_price_index, elasticity)
    weight = elasticity / (elasticity - 1) * log(s)

    return weight

end

"""
  Compute utility function parameters (γ) by region considered in the model for the firms production function.
  This is a wrapper around parameters_by_region, but here we have to consider inputs to each firm from each firm
  which increases the dimension of the parameter array.

  Matlab code reference ComputeTheta.m lines 257-259
"""
function firms_parameters_by_region(elasticity::T,
                                 log_price_uk::Vector{T},log_price_eu::Vector{T},log_price_world::Vector{T},
                                 input_uk::Matrix{T},input_eu::Matrix{T},input_world::Matrix{T}) where {T <: Real}

    parameters = Array{Float64}(undef, length(log_price_uk), length(log_price_uk), 3)

    for i in axes(input_uk, 1)
        for j in axes(input_uk, 2)
            parameters[i,j,:] .= parameters_by_region(elasticity, log_price_uk[j], log_price_eu[j], log_price_world[j],
                                                input_uk[i,j], input_eu[i,j], input_world[i,j])
        end
    end

    return parameters

end

"""
  Compute the total parameter (γM) for the firms input utility function.

  Matlab code reference ComputeTheta.m line 251
"""
function total_input_parameters(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, log_wages, tau)

    parameters = Vector{T}(undef, length(log_price_index))

    for i in axes(log_price_index,1)
        parameters[i] = weight_kernel(input[i], exp(log_price_index[i])/tauP, elasticity)
    end

    return parameters

    #return [weight_kernel(input[i], exp(log_price_index[i]) / tauP, elasticity) for i in eachindex(input, log_price_index)]

end

"""
  Compute the total parameter (γH) for the firms labor utility function.

  Matlab code reference ComputeTheta.m line 249
"""
function total_labor_parameters(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, log_wages, tau)
    return weight_kernel(labor, exp(log_wages) / tauP, elasticity)

end

"""
  Compute the total parameter (γK) for the firms capital utility function.

  Matlab code reference ComputeTheta.m line 254
"""
function total_capital_parameters(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, log_wages, tau)
    tauY = (1 - tau) * output / demand0
    return capital ^ elasticity * (tauY / tauP) ^ (elasticity - 1)

end

function weight_kernel(a::T, b::T, elasticity::T) where {T <: Real}

    return a / b ^ (1 - elasticity)

end

function log_weight_kernel(a::T, b::T, elasticity::T) where {T <: Real}

    return log(a) + (elasticity - 1) * log(b)

end

function sum_kernel(var::Vector{T}, logP::Vector{T}, elasticity::T) where {T <: Real}

    s = 0.0
    for i in axes(logP,1)
        s += var[i] ^ (1.0 / elasticity) * exp((elasticity - 1.0) * logP[i] / elasticity)
    end
    return s

    #return sum(x -> x[1] ^ (1 / elasticity) * exp((elasticity - 1) * x[2] / elasticity), zip(var, logP))

end

"""
  Compute the productivity shock mean μ
"""
function productivity_shock_mean(elasticity::T, log_price_uk::T, log_price_index::Vector{T}, input::Vector{T}, capital::T, demand0::T, output::T, labor::T, log_wages::T, tau::T) where {T <: Real}

    tauP = tauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, log_wages, tau)
    return tauP / ((1-tau) * exp(log_price_uk))

    # logTauP = logTauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, log_wages, tau)
    # logMu = logTauP - log(1.0 - tau) - log_price_uk
    # return exp(logMu)

end

function logTauPdMu(elasticity::T, log_price_index::Vector{T}, input::Vector{T}, capital::T, demand0::T, output::T, labor::T, log_wages::T, tau::T) where {T <: Real}

    s = Supergrassi.sum_kernel(input, log_price_index, elasticity)
    k = capital_fun(capital, tau, output, demand0, elasticity)
    h = labor_fun(labor, log_wages, elasticity)

    return elasticity / ( elasticity - 1 ) * log(s + k + h)

end

function tauPdMu(elasticity::T, log_price_index::Vector{T}, input::Vector{T}, capital::T, demand0::T, output::T, labor::T, log_wages::T, tau::T) where {T <: Real}

    s = Supergrassi.sum_kernel(input, log_price_index, elasticity)
    k = capital_fun(capital, tau, output, demand0, elasticity)
    h = labor_fun(labor, log_wages, elasticity)

    return (s + k + h) ^ (elasticity / ( elasticity - 1 ) )

end

function capital_fun(capital::T, tau::T, output::T, demand0::T, elasticity::T) where T

    return capital * exp((elasticity - 1) / elasticity * log((1 - tau) * output / demand0))

end

function labor_fun(labor::T, log_wages::T, elasticity::T) where T

    return labor ^ ( 1 / elasticity ) * exp(( elasticity - 1 ) / elasticity * log_wages)

end
