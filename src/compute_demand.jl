function compute_demand(n::Int,
                        elasticity::@NamedTuple{_::T, a::T},
                        log_price_uk::Vector{T},
                        log_price_eu::Vector{T},
                        log_price_world::Vector{T},
                        goods_consumption::@NamedTuple{composite::Vector{T},
                                                       uk::Vector{T},
                                                       eu::Vector{T},
                                                       world::Vector{T}}) where {T <: Real}


    # Inputs
    # elasticity.{a} # parms.epsilon
    # log_price_{uk, eu, w} # logPd, parms.logP{eu,w}
    # goods_consumption.{uk, eu, w} #data.fValue{UK, EU, W}

    # Outputs
    # alpha.{uk, eu, w}
    # log_alpha.{uk, eu, w}
    # ∂_log_alpha.{uk, eu, w}
    
    alpha = Dict("composite" => zeros(n), "uk" => zeros(n), "eu" => zeros(n), "world" => zeros(n))
    #log_alpha = Dict(composite => ones(n,1), uk => ones(n,1), eu => ones(n,1), world => ones(n,1))
    #∂_log_alpha = Dict(composite => zeros(n,1), uk => zeros(n,1), eu => zeros(n,1), world => zeros(n,1))

    log_consumer_price_index = zeros(n)
    composite_price_index = zeros(n)
    
    @. log_consumer_price_index = elasticity.a / (elasticity.a - 1.0) * log(
        goods_consumption.uk ^ -elasticity.a * exp((elasticity.a - 1.0) * log_price_uk / elasticity.a) +
        goods_consumption.eu ^ -elasticity.a * exp((elasticity.a - 1.0) * log_price_eu / elasticity.a) +
        goods_consumption.world  ^ -elasticity.a * exp((elasticity.a - 1.0) * log_price_world  / elasticity.a)
    ) #logPf

    @. composite_price_index = goods_consumption.composite ^ -elasticity._ * exp((elasticity._ - 1.0) * log_consumer_price_index / elasticity._)
    
    @. alpha["composite"] = goods_consumption.composite * exp( (elasticity.a - 1.0) *
        (log_consumer_price_index - elasticity._/(elasticity._ - 1.0) * log( sum( composite_price_index) ) ) )

    @. alpha["uk"] = goods_consumption.uk * exp((elasticity.a - 1.0) * (log_price_uk - log_consumer_price_index))
    @. alpha["eu"] = goods_consumption.uk * exp((elasticity.a - 1.0) * (log_price_eu - log_consumer_price_index))
    @. alpha["world"] = goods_consumption.uk * exp((elasticity.a - 1.0) * (log_price_world - log_consumer_price_index))

    return alpha
    
end
