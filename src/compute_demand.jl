# function compute_demand(n,elasticity,log_price_uk,log_price_eu,log_price_world,goods_consumption)
function compute_demand(elasticity::T,
                        log_price_uk::T,
                        log_price_eu::T,
                        log_price_world::T,
                        goods_consumption_uk::T,
                        goods_consumption_eu::T,
                        goods_consumption_world::T) where {T <: Real}

    # Inputs
    # elasticity # parms.epsilon_a
    # log_price_{uk, eu, w} # logPd, parms.logP{eu,w}
    # goods_consumption_{uk, eu, w} #data.fValue{UK, EU, W}

    # Outputs
    # alpha.{uk, eu, w}
    # # log_alpha.{uk, eu, w}
    # # ∂_log_alpha.{uk, eu, w}
    
    log_consumer_price_index = elasticity / (elasticity - 1.0) * log(
        goods_consumption_uk ^ -elasticity * exp((elasticity - 1.0) * log_price_uk / elasticity) +
        goods_consumption_eu ^ -elasticity * exp((elasticity - 1.0) * log_price_eu / elasticity) +
        goods_consumption_world  ^ -elasticity * exp((elasticity - 1.0) * log_price_world  / elasticity)
    ) #logPf

    alpha_uk = goods_consumption_uk * exp((elasticity - 1.0) * (log_price_uk - log_consumer_price_index))
    alpha_eu = goods_consumption_eu * exp((elasticity - 1.0) * (log_price_eu - log_consumer_price_index))
    alpha_world = goods_consumption_world * exp((elasticity - 1.0) * (log_price_world - log_consumer_price_index))

    return alpha_uk, alpha_eu, alpha_world
    
end
