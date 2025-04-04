function logTauPdMu(elasticity::T, log_price_index::Vector{T}, input::Vector{T}, capital::T, demand0::T, output::T, labor::T, wages::T, tau::T) where {T <: Real}

    s = Supergrassi.sum_kernel(input, log_price_index, elasticity)
    k = capital_fun(capital, tau, output, demand0, elasticity)
    h = labor_fun(labor, wages, elasticity)

    logTauPdMu = elasticity / ( elasticity - 1.0 ) * log(s + k + h)

    return logTauPdMu

end

function tauPdMu(elasticity::T, log_price_index::Vector{T}, input::Vector{T}, capital::T, demand0::T, output::T, labor::T, log_wages::T, tau::T) where {T <: Real}

    s = Supergrassi.sum_kernel(input, log_price_index, elasticity)
    k = capital_fun(capital, tau, output, demand0, elasticity)
    h = labor_fun(labor, log_wages, elasticity)

    tauPdMu = (s + k + h) ^ (elasticity / ( elasticity - 1.0 ) )

    return tauPdMu

end

function capital_fun(capital::T, tau::T, output::T, demand0::T, elasticity::T) where T

    return capital * exp((elasticity - 1.0) / elasticity * log((1 - tau) * output / demand0))

end

function labor_fun(labor::T, log_wages::T, elasticity::T) where T

    return labor ^ ( 1.0 / elasticity ) * exp(( elasticity - 1.0 ) / elasticity * log_wages)

end

function input_weights(elasticity::T,
                       log_price_uk::Vector{T},log_price_eu::Vector{T},log_price_world::Vector{T},
                       input_uk::Matrix{T},input_eu::Matrix{T},input_world::Matrix{T}) where {T <: Real}

    # Inputs
    # elasticity # parms.xi_a
    # log_price_{uk, eu, w} # logPd, parms.logP{eu,w}
    # input_{uk, eu, w} #data.mValue{UK, EU, W}

    # Outputs
    # weights{uk, eu, w}

    weights = Array{Float64}(undef, length(log_price_uk), length(log_price_uk), 6)

    for i in axes(input_uk, 1)
        for j in axes(input_uk, 2)
            weights[i,j,:] .= Supergrassi.consumption_weights(elasticity, log_price_uk[j], log_price_eu[j], log_price_world[j],
                                                              input_uk[i,j], input_eu[i,j], input_world[i,j])
        end
    end

    return weights

end

function total_input_weight(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    logTau = logTauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, wages, tau)

    weight = Vector{T}(undef, length(log_price_index))

    for i in axes(log_price_index,1)
        weight[i] = weight_kernel(input[i], exp(log_price_index[i] - logTau), elasticity)
    end

    return weight

end

function total_labor_weight(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, log_wages, tau)
    weight = weight_kernel(labor, exp(log_wages) / tauP, elasticity)

    return weight

end

function total_capital_weight(log_price_index::Vector{T}, input::Vector{T},
                            capital::T, demand0::T, output::T, labor::T, log_wages::T, elasticity::T, tau::T) where T

    length(log_price_index) == length(input) || error()
    tauP = tauPdMu(elasticity, log_price_index, input, capital, demand0, output, labor, log_wages, tau)
    tauY = (1 - tau) * output / demand0
    weight = capital ^ elasticity * (tauY / tauP) ^ (elasticity - 1.0)

    return weight

end
