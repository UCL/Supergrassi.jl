"""
    function parameter_by_region(elasticity::T, quantity_region::T, logP_region::T, logP::T) where {T <: Real}

Compute an utility function parameter for a single region

# Arguments
- `elasticity::Real` : substitution elasticity
- `quantitity_region::Real` : quantity for the region. One of [f, x1, x2, I, m]
- `logP_region::Real` : log of price index for region
- `logP::Real` :  log of aggregate price index
"""
function parameter_by_region(elasticity::T, quantity_region::T, logP_region::T, logP::T) where {T <: Real}

    return quantity_region == 0.0 ? 0.0 : quantity_region * exp((elasticity - 1) * (logP_region - logP))

end

"""
    function log_parameter_by_region(elasticity::T, quantity_region::T, logP_region::T, logP::T) where {T <: Real}

Compute log of an utility function parameter for a single region. See [parameter_by_region](@ref)
"""
function log_parameter_by_region(elasticity::T, quantity_region::T, logP_region::T, logP::T) where {T <: Real}

    return quantity_region == 0.0 ? 0.0 : log(quantity_region) + (elasticity - 1) * (logP_region - logP)

end

"""
    function parameters_by_region(elasticity::T,log_price_uk::T,log_price_eu::T,log_price_world::T,quantity_uk::T,quantity_eu::T,quantity_world::T) where {T <: Real}

Compute utility function parameters by region considered in the model (uk, eu, rest of world).
Parameters of this type appear in multiple utility functions in the paper, and are annotated
(at least) α, β1, β2, γ and ρ.

This is refactored from the Matlab code in e.g. ComputeTheta lines 62-64

# Arguments
 - `elasticity::Real` : substitution elasticity. One of [ϵ, χ1, χ2, η, ξ]
 - `log_price_{uk, eu, w}::Real` : log price indices for each region [log(p)]
 - `quantity_{uk, eu, w}::Real` :  quantity for the region. One of [f, x1, x2, I, m]
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
    function log_parameters_by_region(elasticity::T,log_price_uk::T,log_price_eu::T,log_price_world::T,quantity_uk::T,quantity_eu::T,quantity_world::T) where {T <: Real}


Compute log of utility function parameters. See [parameters_by_region](@ref).
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
    function total_parameters(log_price_index::Vector{T}, quantity::Vector{T}, elasticity::T ) where {T <: Real}

Compute the aggregate utility function parameter. Parameters of this type appear in multiple utility function
in the paper, and are annotated (at least) α, β1, β2, γ and ρ.

This is refactored from the Matlab code in e.g. ComputeTheta.m line 59

# Arguments
- `log_price_index::Vector{Real}` : price index computed by [log_price_index](@ref)
- `quantity:Vector{Real}` : aggregate quantity [f, x1, x2, I, m]
- `elasticity::Real` : substitution elasticity. One of [ϵ, χ1, χ2, η, ξ]
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
    function log_eu_expenditure_on_uk_exports(log_price_index::Vector{T}, quantity::Vector{T},Ex::T, ETilde::T, ePx::T, PTilde::T, elasticity::T,elasticity_tilde::T) where {T <: Real}

Compute the share of EU expenditure on UK exports, defined in section 2.2.1 of the main paper.

for the export parameters β, this is the share of foreign expenditure on UK exports.

# Arguments
- `log_price_index::Vector{Real}` : price index computed by [log_price_index](@ref)
- `quantity::Vector{Real}` : Quantity of domestically-produced good i exported to the EU
- `Ex::Real` : sum of imports
- `Etilde::Real`: EU expenditure on UK exports
- `ePx::Real` : exchange rate to foreign currency
- `Ptilde::Real` : UK exportprice index
- `elasticity::Real` : substitution elasticity
- `elasticity_tilde::Real` : substitution from uk to other elasticity

# See also [ParamsStruct](@ref).
"""
function log_eu_expenditure_on_uk_exports(log_price_index::Vector{T}, quantity::Vector{T},
                                          Ex::T, ETilde::T, ePx::T, PTilde::T, elasticity::T,
                                          elasticity_tilde::T) where {T <: Real}

    length(log_price_index) == length(quantity) || error()
    logPBar = log_total_price_index(elasticity, log_price_index, quantity)

    return log_weight_kernel(Ex/ETilde, exp(logPBar) * ePx / PTilde, elasticity_tilde)

end

"""
    function price_index(elasticity::T,log_price_uk::T, log_price_eu::T, log_price_world::T,demand_uk::T, demand_eu::T, demand_world::T) where {T <: Real}

Compute the the consumer price index defined in equation 2.3 of the main paper as \bar{p}.
Matlab code reference e.g. ComputeTheta.m line 49.

# Arguments
- `elasticity::Real` : substitution elasticity
- `log_price_{uk, eu, w}::Real` : log price indices for each region [log(p)]
- `demand_{uk, eu, w}::Real` :  demand for the region. One of [f, x1, x2, I, m]
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

"""
    function log_price_index(elasticity::T,log_price_uk::T, log_price_eu::T, log_price_world::T,demand_uk::T, demand_eu::T, demand_world::T) where {T <: Real}

Compute the log of the consumer price index. See [price_index](@ref)
"""
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
    function log_total_price_index(elasticity::T, log_price_index::Vector{T}, quantity::Vector{T}) where {T <: Real}

Compute the total consumer price index defined in equation 2.7 of the main paper as \bar{P}.
Matlab code reference e.g. ComputeTheta.m line 58.

# Arguments
- `elasticity::Real` : substitution elasticity
- `log_price_index::Vector{Real}` : price index computed by [log_price_index](@ref)
- `quantity:Vector{Real}` : aggregate quantity [f, x1, x2, I, m]
"""
function log_total_price_index(elasticity::T, log_price_index::Vector{T}, quantity::Vector{T}) where {T <: Real}

    s = sum_kernel(quantity, log_price_index, elasticity)
    weight = elasticity / (elasticity - 1) * log(s)

    return weight

end

"""
    function total_input_parameters(log_price_index::Vector{T}, input::Vector{T}, surplus::T, capital::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

Compute the total parameter (γM) for the firms input utility function.

Matlab code reference ComputeTheta.m line 251

# Arguments
- `log_price_index::Vector{Real}` : price index computed by [log_price_index](@ref)
- `input::Vector{Real}` : intermediate input from firm i to firms j. mValue in matlab code
- `suprlus::Real` : firms suprlus, data.kValue
- `capital::Real` : firms capital, data.k0
- `output::Real` : total household use, data.yValue
- `labor::Real` : payments for labor, data.hValue
- `log_wages::Real` : aggregate wages of households
- `elasticity::Real` : substitution elasticity
- `tau::Real` : ratio of taxes to household use

# See also
[compute_agg_wages](@ref)
"""
function total_input_parameters(log_price_index::Vector{T}, input::Vector{T},
                            surplus::T, capital::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, surplus, capital, output, labor, log_wages, tau)

    parameters = Vector{T}(undef, length(log_price_index))

    for i in axes(log_price_index,1)
        parameters[i] = weight_kernel(input[i], exp(log_price_index[i])/tauP, elasticity)
    end

    return parameters

    #return [weight_kernel(input[i], exp(log_price_index[i]) / tauP, elasticity) for i in eachindex(input, log_price_index)]

end

"""
    function total_labor_parameters(log_price_index::Vector{T}, input::Vector{T}, surplus::T, capital::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

Compute the total parameter (γH) for the firms labor utility function.

Matlab code reference ComputeTheta.m line 249


# Arguments
- `log_price_index::Vector{Real}` : price index computed by [log_price_index](@ref)
- `input::Vector{Real}` : intermediate input from firm i to firms j. mValue in matlab code
- `suprlus::Real` : firms suprlus, data.kValue
- `capital::Real` : firms capital, data.k0
- `output::Real` : total household use, data.yValue
- `labor::Real` : payments for labor, data.hValue
- `log_wages::Real` : aggregate wages of households
- `elasticity::Real` : substitution elasticity
- `tau::Real` : ratio of taxes to household use

# See also
[compute_agg_wages](@ref)
"""
function total_labor_parameters(log_price_index::Vector{T}, input::Vector{T},
                            surplus::T, capital::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, surplus, capital, output, labor, log_wages, tau)
    return weight_kernel(labor, exp(log_wages) / tauP, elasticity)

end

"""
    function total_surplus_parameters(log_price_index::Vector{T}, input::Vector{T}, surplus::T, capital::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

Compute the total parameter (γK) for the firms capital utility function.

Matlab code reference ComputeTheta.m line 254

# Arguments
- `log_price_index::Vector{Real}` : price index computed by [log_price_index](@ref)
- `input::Vector{Real}` : intermediate input from firm i to firms j. mValue in matlab code
- `suprlus::Real` : firms suprlus, data.kValue
- `capital::Real` : firms capital, data.k0
- `output::Real` : total household use, data.yValue
- `labor::Real` : payments for labor, data.hValue
- `log_wages::Real` : aggregate wages of households
- `elasticity::Real` : substitution elasticity
- `tau::Real` : ratio of taxes to household use

# See also
[compute_agg_wages](@ref)
"""
function total_surplus_parameters(log_price_index::Vector{T}, input::Vector{T},
                            surplus::T, capital::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, surplus, capital, output, labor, log_wages, tau)
    tauY = (1 - tau) * output / capital
    return surplus ^ elasticity * (tauY / tauP) ^ (elasticity - 1)

end

"""
    weight_kernel(a::T, b::T, elasticity::T) where {T <: Real}

Computes a/b^(1-elasticity)

# Argument:
- `a`: Numerator value.
- `b`: Denominator value.
- `elasticity`: substitution elasticity parameter
"""
function weight_kernel(a::T, b::T, elasticity::T) where {T <: Real}

    return a / b ^ (1 - elasticity)

end

"""
    log_weight_kernel(a::T, b::T, elasticity::T) where {T <: Real}

computes log(a/b^(1-elasticity))

# Argument:
- `a`: Numerator value (in log form).
- `b`: Denominator base value (in log form).
- `elasticity`: substitution elasticity parameter
"""
function log_weight_kernel(a::T, b::T, elasticity::T) where {T <: Real}

    return log(a) + (elasticity - 1) * log(b)

end

"""
    sum_kernel(var::Vector{T}, logP::Vector{T}, elasticity::T) where {T <: Real}

# Argument:
- `var`: Vector of variable values.
- `logP`: Vector of log prices.
- `elasticity`: substitution elasticity parameter
"""
function sum_kernel(var::Vector{T}, logP::Vector{T}, elasticity::T) where {T <: Real}

    s = 0.0
    for i in axes(logP,1)
        s += var[i] ^ (1.0 / elasticity) * exp((elasticity - 1.0) * logP[i] / elasticity)
    end
    return s

    #return sum(x -> x[1] ^ (1 / elasticity) * exp((elasticity - 1) * x[2] / elasticity), zip(var, logP))

end


function productivity_shock_mean(elasticity::T, log_price_uk::T, log_price_index::Vector{T}, input::Vector{T}, surplus::T, capital::T, output::T, labor::T, log_wages::T, tau::T) where {T <: Real}    

    tauP = tauPdMu(elasticity, log_price_index, input, surplus, capital, output, labor, log_wages, tau)
    return tauP / ((1-tau) * exp(log_price_uk))

    # logTauP = logTauPdMu(elasticity, log_price_index, input, surplus, capital, output, labor, log_wages, tau)
    # logMu = logTauP - log(1.0 - tau) - log_price_uk
    # return exp(logMu)

end

function logTauPdMu(elasticity::T, log_price_index::Vector{T}, input::Vector{T}, surplus::T, capital::T, output::T, labor::T, log_wages::T, tau::T) where {T <: Real}

    s = Supergrassi.sum_kernel(input, log_price_index, elasticity)
    k = capital_fun(surplus, tau, output, capital, elasticity)
    h = labor_fun(labor, log_wages, elasticity)

    return elasticity / ( elasticity - 1 ) * log(s + k + h)

end

function tauPdMu(elasticity::T, log_price_index::Vector{T}, input::Vector{T}, surplus::T, capital::T, output::T, labor::T, log_wages::T, tau::T) where {T <: Real}

    s = Supergrassi.sum_kernel(input, log_price_index, elasticity)
    k = capital_fun(surplus, tau, output, capital, elasticity)
    h = labor_fun(labor, log_wages, elasticity)

    return (s + k + h) ^ (elasticity / ( elasticity - 1 ) )

end

"""
    capital_fun(surplus::T, tau::T, output::T, capital::T, elasticity::T) where T

# Argument:
- `surplus`: Surplus value.
- `tau`: Tax parameter.
- `output`: total household use, data.yValue.
- `capital`: firms capital.
- `elasticity`: substitution elasticity parameter.
"""
function capital_fun(surplus::T, tau::T, output::T, capital::T, elasticity::T) where T

    return surplus * exp((elasticity - 1) / elasticity * log((1 - tau) * output / capital))

end

"""
    labor_fun(labor::T, log_wages::T, elasticity::T) where T

# Argument:
- `labor`: Labor input.
- `log_wages`: Log of wages.
- `elasticity`: Elasticity parameter.

# See also
[compute_agg_wages](@ref)
"""
function labor_fun(labor::T, log_wages::T, elasticity::T) where T

    return labor ^ ( 1 / elasticity ) * exp(( elasticity - 1 ) / elasticity * log_wages)

end
