function compute_objective_function(log_price_uk::Vector{<:Number}, zOC::Vector{<:Number}, data::CleanData, params::Parameters)

    tau = (data.industry.tax.products .+ data.industry.tax.production) ./ data.industry.regional.total_use.agg

    mu = params.production.shock_mean
    gammaK = params.production.input_capital
    k0 = data.industry.capital.current_year
    xi = params.constants.elasticities.production.substitution

    excess_demand = intermediate_goods_price_index(log_price_uk, zOC, tau, mu, gammaK, k0, xi)

    w = sqrt.(data.industry.regional.total_use.agg)
    e = w.*(excess_demand .- data.industry.regional.total_use.agg)
    objective_value = 0.5 * sum(e.^2)

    return objective_value

end

function compute_objective_function(x::Vector{<:Number}, data::CleanData, params::Parameters)

    n = length(data.industry.regional.total_use.agg)

    log_price_uk = x[1:n]
    zOC = x[(n+1):end]

    return compute_objective_function(log_price_uk, zOC, data, params)

end