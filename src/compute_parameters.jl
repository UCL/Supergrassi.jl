using Supergrassi
using DataFrames
using Enzyme

function compute_all_parameters(data::RegionalData, elasticities::Elasticities, prices::DataFrame)

    α, ∂α, Jacα = compute_parameter(data.consumption, elasticities.consumption, prices)
    β1, ∂β1, Jacβ1 = compute_parameter(data.export_eu, elasticities.export_eu, prices)
    β2, ∂β2, Jacβ2 = compute_parameter(data.export_world, elasticities.export_world, prices)
    ρ, ∂ρ, Jacρ = compute_parameter(data.investment, elasticities.investment, prices)

    γM, ∂γM = compute_parameter(data.input_matrices, elasticities.production, prices)

end

function compute_parameter(demand::DataFrame, elasticity::Elasticity, prices)

    n = nrow(demand)
    empty_col = missings(Float64, n)
    val = DataFrame([demand.industry, empty_col, empty_col, empty_col], names(demand)[1:4])
    grad = DataFrame([demand.industry, empty_col, empty_col, empty_col], names(demand)[1:4])

    for row in 1:n

        param_regional = gradient(ForwardWithPrimal,
                                  Supergrassi.parameters_by_region,
                                  Const(elasticity.armington),
                                  prices.uk[row],
                                  Const(prices.eu[row]),
                                  Const(prices.world[row]),
                                  Const(demand.uk[row]),
                                  Const(demand.eu[row]),
                                  Const(demand.world[row]))

        val.uk[row] = param_regional.val[1]
        val.eu[row] = param_regional.val[2]
        val.world[row] = param_regional.val[3]
        grad.uk[row] = param_regional.derivs[2][1]
        grad.eu[row] = param_regional.derivs[2][2]
        grad.world[row] = param_regional.derivs[2][3]

    end

    logPf = Supergrassi.log_price_index.(elasticity.armington,
                                         prices.uk, prices.eu, prices.world,
                                         demand.uk, demand.eu, demand.world)
    param = jacobian(ForwardWithPrimal, Supergrassi.total_parameters, logPf,
                     Const(demand.agg), Const(elasticity.substitution))

    val.agg .= param.val

    jac = DataFrame(first(param.derivs), demand.industry)

    return val, grad, jac

end

function compute_parameter(demand::InputMatrices, elasticity::Elasticity, prices)

    n = nrow(demand.agg)

    val = zeros(n,n,3)
    grad = zeros(n,n,3)

    for row in 1:n
        for col in 1:n

            param_regional = gradient(ForwardWithPrimal,
                                      Supergrassi.parameters_by_region,
                                      Const(elasticity.armington),
                                      prices.uk[col],
                                      Const(prices.eu[col]),
                                      Const(prices.world[col]),
                                      Const(demand.uk[row, col]),
                                      Const(demand.eu[row, col]),
                                      Const(demand.world[row, col]))

            val[row, col, :] .= param_regional.val
            grad[row, col, :] .= param_regional.derivs[2]

        end
    end

    return val, grad

end
