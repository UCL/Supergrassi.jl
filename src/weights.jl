function weights_by_region(elasticity::T,
                           log_price_uk::T,log_price_eu::T,log_price_world::T,
                           consumption_uk::T,consumption_eu::T,consumption_world::T) where {T <: Real}

    # Inputs
    # elasticity # parms.epsilon_a
    # log_price_{uk, eu, w} # logPd, parms.logP{eu,w}
    # consumption_{uk, eu, w} #data.fValue{UK, EU, W}

    # Outputs
    # weight_{uk, eu, w}

    logP = log_price_index(elasticity, log_price_uk, log_price_eu, log_price_world, consumption_uk, consumption_eu, consumption_world)

    weight_uk = consumption_uk * exp((elasticity - 1.0) * (log_price_uk - logP))
    weight_eu = consumption_eu * exp((elasticity - 1.0) * (log_price_eu - logP))
    weight_world = consumption_world * exp((elasticity - 1.0) * (log_price_world - logP))

    return weight_uk, weight_eu, weight_world

end

function log_weights_by_region(elasticity::T,
                               log_price_uk::T,log_price_eu::T,log_price_world::T,
                               consumption_uk::T,consumption_eu::T,consumption_world::T) where {T <: Real}

    # Inputs
    # elasticity # parms.epsilon_a
    # log_price_{uk, eu, w} # logPd, parms.logP{eu,w}
    # consumption_{uk, eu, w} #data.fValue{UK, EU, W}

    # Outputs
    # log(weight_{uk, eu, w})

    logP = log_price_index(elasticity, log_price_uk, log_price_eu, log_price_world, consumption_uk, consumption_eu, consumption_world)

    log_weight_uk = log(consumption_uk) + (elasticity - 1.0) * (log_price_uk - logP)
    log_weight_eu = log(consumption_eu)  + (elasticity - 1.0) * (log_price_eu - logP)
    log_weight_world = log(consumption_world) + (elasticity - 1.0) * (log_price_world - logP)

    return log_weight_uk, log_weight_eu, log_weight_world

end

function total_weights(log_price_index::Vector{T}, consumption::Vector{T}, elasticity::T ) where {T <: Real}

    # Both e.g. parms.beta1 and parms.beta1Tilde are of this form.
    # in beta1, consumption = data.x1Value and log_total_price_index = LogPX1Bar
    # in beta1Tilde, consumption = EX1/E1Tilde and log_total_price_index = LogPX1Bar + parms.fx_EUR.
    # TODO: Think about if this should be separated into different functions
    # TODO: Check if this is still true

    length(log_price_index) == length(consumption) || error()

    logPBar = log_total_price_index(elasticity, log_price_index, consumption)
    weights = Vector{T}(undef, length(log_price_index))

    for i in axes(log_price_index,1)
        weights[i] = weight_kernel(consumption[i], exp(log_price_index[i] - logPBar), elasticity)
    end

    return weights

end

function log_eu_expenditure_on_uk_exports(log_price_index::Vector{T}, consumption::Vector{T},
                                          Ex::T, ETilde::T, ePx::T, PTilde::T, elasticity::T,
                                          elasticity_tilde::T) where {T <: Real}

    length(log_price_index) == length(consumption) || error()
    logPBar = log_total_price_index(elasticity, log_price_index, consumption)

    return log_weight_kernel(Ex/ETilde, exp(logPBar) * ePx / PTilde, elasticity_tilde)

end

function log_price_index(elasticity::T,
                         log_price_uk::T, log_price_eu::T, log_price_world::T,
                         demand_uk::T, demand_eu::T, demand_world::T) where {T <: Real}

    # Abbreviated logP in papers and Matlab code

    return elasticity / (elasticity - 1.0) * log(
        demand_uk ^ (1.0/elasticity) * exp((elasticity - 1.0) * log_price_uk / elasticity) +
        demand_eu ^ (1.0/elasticity) * exp((elasticity - 1.0) * log_price_eu / elasticity) +
        demand_world  ^ (1.0/elasticity) * exp((elasticity - 1.0) * log_price_world  / elasticity)
    )

end

function price_index(elasticity::T,
                     log_price_uk::T, log_price_eu::T, log_price_world::T,
                     demand_uk::T, demand_eu::T, demand_world::T) where {T <: Real}

    # Abbreviated logP in papers and Matlab code

    return (
        demand_uk ^ (1.0/elasticity) * exp((elasticity - 1.0) * log_price_uk / elasticity) +
        demand_eu ^ (1.0/elasticity) * exp((elasticity - 1.0) * log_price_eu / elasticity) +
        demand_world  ^ (1.0/elasticity) * exp((elasticity - 1.0) * log_price_world  / elasticity)
    ) ^ ( elasticity / (elasticity - 1.0) )

end


function log_total_price_index(elasticity::T, log_price_index::Vector{T}, consumption::Vector{T}) where {T <: Real}

    # Abbreviated LogPBar in papers and Matlab code

    s = sum_kernel(consumption, log_price_index, elasticity)
    weight = elasticity/(elasticity - 1.0) * log(s)

    return weight

end

function firms_weights_by_region(elasticity::T,
                                 log_price_uk::Vector{T},log_price_eu::Vector{T},log_price_world::Vector{T},
                                 input_uk::Matrix{T},input_eu::Matrix{T},input_world::Matrix{T}) where {T <: Real}

    # Inputs
    # elasticity # parms.xi_a
    # log_price_{uk, eu, w} # logPd, parms.logP{eu,w}
    # input_{uk, eu, w} #data.mValue{UK, EU, W}

    # Outputs
    # weights{uk, eu, w}

    weights = Array{Float64}(undef, length(log_price_uk), length(log_price_uk), 3)

    for i in axes(input_uk, 1)
        for j in axes(input_uk, 2)
            weights[i,j,:] .= weights_by_region(elasticity, log_price_uk[j], log_price_eu[j], log_price_world[j],
                                                input_uk[i,j], input_eu[i,j], input_world[i,j])
        end
    end

    return weights

end

function total_input_weights(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, log_wages, tau)

    weights = Vector{T}(undef, length(log_price_index))

    for i in axes(log_price_index,1)
        weights[i] = weight_kernel(input[i], exp(log_price_index[i])/tauP, elasticity)
    end

    return weights

end

function total_labor_weights(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, log_wages, tau)
    return weight_kernel(labor, exp(log_wages) / tauP, elasticity)

end

function total_capital_weights(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, log_wages, tau)
    tauY = (1 - tau) * output / demand0
    return capital ^ elasticity * (tauY / tauP) ^ (elasticity - 1.0)

end

function weight_kernel(a::T, b::T, elasticity::T) where {T <: Real}

    return a / b ^ (1.0 - elasticity)

end

function log_weight_kernel(a::T, b::T, elasticity::T) where {T <: Real}

    return log(a) + (elasticity - 1.0) * log(b)

end

function sum_kernel(var::Vector{T}, logP::Vector{T}, elasticity::T) where {T <: Real}

    s = 0.0
    for i in axes(logP,1)
        s += var[i] ^ (1.0 / elasticity) * exp((elasticity - 1.0) * logP[i] / elasticity)
    end

    return s
    
end

function logTauPdMu(elasticity::T, log_price_index::Vector{T}, input::Vector{T}, capital::T, demand0::T, output::T, labor::T, wages::T, tau::T) where {T <: Real}

    s = Supergrassi.sum_kernel(input, log_price_index, elasticity)
    k = capital_fun(capital, tau, output, demand0, elasticity)
    h = labor_fun(labor, wages, elasticity)

    return elasticity / ( elasticity - 1.0 ) * log(s + k + h)

end

function tauPdMu(elasticity::T, log_price_index::Vector{T}, input::Vector{T}, capital::T, demand0::T, output::T, labor::T, log_wages::T, tau::T) where {T <: Real}

    s = Supergrassi.sum_kernel(input, log_price_index, elasticity)
    k = capital_fun(capital, tau, output, demand0, elasticity)
    h = labor_fun(labor, log_wages, elasticity)

    return (s + k + h) ^ (elasticity / ( elasticity - 1.0 ) )

end

function capital_fun(capital::T, tau::T, output::T, demand0::T, elasticity::T) where T

    return capital * exp((elasticity - 1.0) / elasticity * log((1 - tau) * output / demand0))

end

function labor_fun(labor::T, log_wages::T, elasticity::T) where T

    return labor ^ ( 1.0 / elasticity ) * exp(( elasticity - 1.0 ) / elasticity * log_wages)

end
