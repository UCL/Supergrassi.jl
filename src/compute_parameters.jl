using Supergrassi
using DataFrames
using Enzyme

struct Values{T<:Real}

    uk::Vector{T}
    eu::Vector{T}
    world::Vector{T}
    tilde::Vector{T}
    agg::Vector{T}

end

struct Derivatives{T<:Real}

    uk::Vector{T}
    eu::Vector{T}
    world::Vector{T}
    tilde::Vector{T}
    agg::Matrix{T}

end


"""
  Compute all utility function parameters from regional data, elasticities and prices.
  Currently missing most of the γ parameters and a data structure for the output
"""
function compute_all_parameters(data::RegionalData, elasticities::Elasticities, prices::DataFrame)

    α, ∂α = compute_parameter(data.consumption, elasticities.consumption, prices)

    β1, ∂β1 = compute_parameter(data.export_eu, elasticities.export_eu, prices)
    β2, ∂β2 = compute_parameter(data.export_world, elasticities.export_world, prices)
    ρ, ∂ρ = compute_parameter(data.investment, elasticities.investment, prices)

    γM, ∂γM = compute_parameter(data.input_matrices, elasticities.production, prices)

end

"""
  Compute 1d utility function parameters from a regional demand data frame and the corresponding elasticity.
  Currently missing the tilde parameters for exports.
"""
function compute_parameter(demand::DataFrame, elasticity::Elasticity, prices, fun = Supergrassi.parameters_by_region)

    n = nrow(demand)
    v0 = Vector{Float64}(undef, n)
    m0 = Matrix{Float64}(undef, n, n)

    val = Values(v0, similar(v0), similar(v0), similar(v0), similar(v0))
    grad = Derivatives(similar(v0), similar(v0), similar(v0), similar(v0), m0)

    for row in 1:n

        param_regional = gradient(ForwardWithPrimal,
                                  fun,
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
    grad.agg .= first(param.derivs)

    return val, grad

end

function compute_EU_share(demand::DataFrame, elasticity::Elasticity, prices, Ex::T, Etilde::T, Ptilde::T) where {T<:Real}

    tilde = gradient(ForwardWithPrimal,
                     Supergrassi.log_eu_expenditure_on_uk_exports,
                     prices.uk,
                     Const(demand.agg),
                     Const(Ex),
                     Const(ETilde),
                     Const(exchange_rates.usd),
                     Const(PTilde),
                     Const(elasticity.substitution),
                     Const(elasticity.substitution_uk_other))

    return tilde.val, first(tilde.derivs)

end

"""
  Compute the 2d utility function parameters from regional InputMatrices and the corresponding elasticity.
  Currently missing the jacobian, and the 1d parameters gammaL and gammaH
"""
function compute_production_parameter(clean::CleanData, prices)

    n = nrow(clean.industry.regional.input_matrices.agg)

    val = Array{Float64}(undef, n, n, 4)
    grad = Array{Float64}(undef, n, n, 3)
    jac = Array{Float64}(undef, n, n ,n)

    for row in 1:n
        for col in 1:n

            param_regional = gradient(ForwardWithPrimal,
                                      Supergrassi.parameters_by_region,
                                      Const(clean.constants.elasticities.production.armington),
                                      prices.uk[col],
                                      Const(prices.eu[col]),
                                      Const(prices.world[col]),
                                      Const(clean.industry.regional.input_matrices.uk[row, col]),
                                      Const(clean.industry.regional.input_matrices.eu[row, col]),
                                      Const(clean.industry.regional.input_matrices.world[row, col]))

            val[row, col, 1:3] .= param_regional.val
            grad[row, col, 1:3] .= param_regional.derivs[2]

        end

        logPm = Supergrassi.log_price_index.(clean.constants.elasticities.production.armington,
                                             prices.uk, prices.eu, prices.world,
                                             clean.industry.regional.input_matrices.uk[row,:],
                                             clean.industry.regional.input_matrices.eu[row,:],
                                             clean.industry.regional.input_matrices.world[row,:])
        logPm[isinf.(logPm)] .= 0.0

        # We are missing this from our data struct. Temporarily compute it on the fly.
        tau = (clean.industry.tax.products .+ clean.industry.tax.production) ./ clean.industry.regional.total_use.agg

        jacM = jacobian(ForwardWithPrimal,
                        Supergrassi.total_input_parameters,
                        logPm,
                        Const(clean.industry.regional.input_matrices.agg[row,:]),
                        Const(clean.industry.capital.next_year[row]),
                        Const(clean.industry.capital.current_year[row]),
                        Const(clean.industry.total_use.agg[row]),
                        Const(clean.household.income[row]),
                        Const(clean.household.wages[row]),
                        Const(clean.constants.elasticities.production.substitution),
                        Const(tau))

        val[row,:,4] .= jacM.val
        jac[row,:,:] .= first(jacM.derivs)

    end

    return val, grad, jac

end
