function intermediate_goods_price_index(log_price_uk::Vector{T}, zOC::Vector{T},
                                        tau::Vector{T}, mu::Vector{T}, gammaK::Vector{T},
                                        K0::Vector{T}, xi::T) where {T <: Real}

    # Computes the intermediate goods price index, Step 1 of ExcessDemand.m

    # pdYBar = Vector{T}(undef, length(log_price_uk))
    # for i in axes(log_price_uk, 1)
    #     pdYBar[i] = intermediate_goods_price_index(log_price_uk[i], zOC[i], tau[i], mu[i], gammaK[i], K0[i], xi)
    # end
    # return pdYBar

    return [intermediate_goods_price_index(log_price_uk[i], zOC[i], tau[i], mu[i], gammaK[i], K0[i], xi) for i in eachindex(log_price_uk, zOC, tau, mu, gammaK, K0)]

end

function intermediate_goods_price_index(log_price_uk::T, zOC::T, tau::T, mu::T, gammaK::T, K0::T, xi::T) where {T <: Real}

    # Computes the intermediate goods price index, Step 1 of ExcessDemand.m

    # log_price_uk: part of x
    # zOC: part of x
    # tau: parameter_interface L173
    # mu: shock_mean
    # gammaK: in params
    # K0: in industry_data_capital
    # xi: production elasticities


    return ( (1 - tau) * exp(log_price_uk) * mu * K0 * gammaK ^ (1 / (xi - 1) )
             * (1 - (exp(zOC) ) / ( 1 + exp(zOC) ) ) ^ (xi / (1 - xi) ) / (1 - tau) )

end
