function log_price_index(elasticity::T, log_price_uk::T, log_price_eu::T, log_price_world::T,
                         input_uk::Vector{T}, input_eu::Vector{T}, input_world::Vector{T}) where T

    logPm = elasticity/(elasticity - 1.0) * log(
        input_uk ^ (1.0/elasticity) * exp((elasticity - 1.0) / elasticity * log_price_uk)
        + input_eu ^ (1.0/elasticity) * exp((elasticity - 1.0) / elasticity * log_price_eu)
        + input_world ^ (1.0/elasticity) * exp((elasticity - 1.0) / elasticity * log_price_world)
    )

    return logPm

end

function logTauPdMu(elasiticty::T, log_price_index::Vector{T}, input::Vector{T}, capital::T, demand0::T, output::T, labor::T, wages::T, tau::T)

    logTauPdMu = elasticity / ( elasticity - 1.0 ) * log(
        sum(input ^( 1.0 / elasticity ) * exp(( elasticity - 1.0 ) / elasticity * log_price_index))
        + capital * exp((elasticity - 1.0) / elasticity * log((1 - tau) * output / demand0))
        + labor ^ ( 1.0 / elasticity ) * exp(( elasticity - 1.0 ) / elasticity * wages)
    )

    return logTauPdMu

end

function tauPdMu(elasiticty::T, log_price_index::Vector{T}, input::Vector{T}, capital::T, demand0::T, output::T, labor::T, wages::T, tau::T)

    tauPdMu =  (
        sum(input ^( 1.0 / elasticity ) * exp(( elasticity - 1.0 ) / elasticity * log_price_index))
        + capital * exp((elasticity - 1.0) / elasticity * log((1 - tau) * output / demand0))
        + labor ^ ( 1.0 / elasticity ) * exp(( elasticity - 1.0 ) / elasticity * wages)
    ) ^ ( elasticity / ( elasticity - 1.0 ) )

    return tauPdMu

end

function input_weights(elasticity::T,
                       log_price_uk::T,log_price_eu::T,log_price_world::T,
                       input_uk::Vector{T},input_eu::Vector{T},input_world::Vector{T}) where {T <: Real}

    # Inputs
    # elasticity # parms.epsilon_a
    # log_price_{uk, eu, w} # logPd, parms.logP{eu,w}
    # input_{uk, eu, w} #data.mValue{UK, EU, W}

    # Outputs
    # weight_{uk, eu, w}

    logP = log_price_index(elasticity, log_price_uk, log_price_eu, log_price_world, input_uk, input_eu, input_world)

    weight_uk = input_uk * exp((elasticity - 1.0) * (log_price_uk - logP))
    weight_eu = input_eu * exp((elasticity - 1.0) * (log_price_eu - logP))
    weight_world = input_world * exp((elasticity - 1.0) * (log_price_world - logP))

    return weight_uk, weight_eu, weight_world, log(weight_uk), log(weight_eu), log(weight_world)

end

function total_input_weight(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    logTauPdMu = logTauPdMu(elasticity, log_price_index, input, capital, demand0, output, lavor, wages)

    weight = Vector{T}(undef, length(log_price_index))
    for i in axes(log_price_index,1)
        weight[i] = weight_kernel(input[i], exp(log_price_index[i] - logTauPdMu), elasticity)
    end

end

function total_labor_weight(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    TauPdMu = logTauPdMu(elasticity, log_price_index, input, capital, demand0, output, lavor, wages)

    weight = weight_kernel(labor, wages - TauPdMu, elasticity)

end

function total_capital_weight(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    TauPdMu = TauPdMu(elasticity, log_price_index, input, capital, demand0, output, lavor, wages)

    weight = weight_kernel(capital, (1 - tau) * output / demand0 - TauPdMu, elasticity)

end
