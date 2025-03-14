function log_price_index(elasticity::T,
                         log_price_uk::T, log_price_eu::T, log_price_world::T,
                         demand_uk::T, demand_eu::T, demand_world::T) where {T <: Real}

    # Abbreviated logPf in papers and Matlab code

    return elasticity / (elasticity - 1.0) * log(
        demand_uk ^ (1.0/elasticity) * exp((elasticity - 1.0) * log_price_uk / elasticity) +
        demand_eu ^ (1.0/elasticity) * exp((elasticity - 1.0) * log_price_eu / elasticity) +
        demand_world  ^ (1.0/elasticity) * exp((elasticity - 1.0) * log_price_world  / elasticity)
    )

end

function consumption_weights(elasticity::T,
                             log_price_uk::T,log_price_eu::T,log_price_world::T,
                             consumption_uk::T,consumption_eu::T,consumption_world::T) where {T <: Real}

    # Inputs
    # elasticity # parms.epsilon_a
    # log_price_{uk, eu, w} # logPd, parms.logP{eu,w}
    # consumption_{uk, eu, w} #data.fValue{UK, EU, W}

    # Outputs
    # weight_{uk, eu, w}

    logPf = log_price_index(elasticity, log_price_uk, log_price_eu, log_price_world, consumption_uk, consumption_eu, consumption_world)

    weight_uk = consumption_uk * exp((elasticity - 1.0) * (log_price_uk - logPf))
    weight_eu = consumption_eu * exp((elasticity - 1.0) * (log_price_eu - logPf))
    weight_world = consumption_world * exp((elasticity - 1.0) * (log_price_world - logPf))

    return weight_uk, weight_eu, weight_world, log(weight_uk), log(weight_eu), log(weight_world)

end

function log_total_price_index(elasticity::T, log_price_index::Vector{T}, consumption::Vector{T}) where {T <: Real}

    # Abbreviated LogPBar in papers and Matlab code

    s = 0.0
    for i in axes(log_price_index,1)
        s += consumption[i] ^ (1.0 / elasticity) * exp((elasticity - 1.0) * log_price_index[i] / elasticity)
    end

    weight = elasticity/(elasticity - 1.0) * log(s)

    return weight

end

function total_consumption_weight(log_price_index::Vector{T}, consumption::Vector{T}, elasticity::T ) where {T <: Real}

    # Both e.g. parms.beta1 and parms.beta1Tilde are of this form.
    # in beta1, consumption = data.x1Value and log_total_price_index = LogPX1Bar
    # in beta1Tilde, consumption = EX1/E1Tilde and log_total_price_index = LogPX1Bar + parms.fx_EUR.
    # TODO: Think about if this should be separated into different functions
    # TODO: Check if this is still true

    length(log_price_index) == length(consumption) || error()

    logPBar = log_total_price_index(elasticity, log_price_index, consumption)
    weight = Vector{T}(undef, length(log_price_index))

    for i in axes(log_price_index,1)
        #weight[i] = consumption[i] * exp((elasticity - 1.0)  * (log_price_index[i] - logPBar))
        weight[i] = weight_kernel(consumption[i], exp(log_price_index[i] - logPBar), elasticity)
    end

    return weight

end

function log_beta_tilde(log_price_index::Vector{T}, consumption::Vector{T},
                        Ex::T, ETilde::T, ePx::T, PTilde::T, elasticity::T, elasticity_tilde::T) where {T <: Real}

    #TODO: Rename this function
    #TODO: Should be possible to refactor this and total_consumption_weight into one

    length(log_price_index) == length(consumption) || error()
    logPBar = log_total_price_index(elasticity, log_price_index, consumption)

    return log_weight_kernel(Ex/ETilde, exp(logPBar) * ePx / PTilde, elasticity_tilde)

end

function weight_kernel(a::T, b::T, elasticity::T) where T

    return a / b ^ (1.0 - elasticity)

end

function log_weight_kernel(a::T, b::T, elasticity::T) where T

    return log(a) + (elasticity - 1.0) * log(b)

end
