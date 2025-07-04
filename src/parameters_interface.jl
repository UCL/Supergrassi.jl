using Supergrassi
using DataFrames
using Enzyme

"""
  Compute all utility function parameters from regional data, elasticities and prices.
  Currently missing some of the γ parameters

  fun parameter should be either parameters_by_region or log_parameters_by_region
"""
function compute_all_parameters(data::CleanData, prices::DataFrame, log_scale::Bool = false)

    #@assert fun in [parameters_by_region, log_parameters_by_region] "fun argument should be either parameters_by_region or log_parameters_by_region"

    reg = data.industry.regional
    constants = data.constants

    imports_uk_share_eu = constants.total_imports_from_uk.eu / (constants.total_imports_from_all_sources.eu
                                                                / constants.exchange_rates.eur )
    imports_uk_share_world = constants.total_imports_from_uk.world / (constants.total_imports_from_all_sources.world
                                                                      / constants.exchange_rates.usd )

    α, ∂α = compute_parameter(reg.consumption, constants.elasticities.consumption, prices, log_scale)
    β1, ∂β1 = compute_parameter(reg.export_eu, constants.elasticities.eu_export_demand, prices, log_scale)
    β1, ∂β1 = compute_foreign_share(β1, ∂β1, reg.export_eu, constants.elasticities.eu_export_demand, prices,
                           imports_uk_share_eu, reg.totals.imports.eu, 1.0, constants.exchange_rates.eur)
    β2, ∂β2 = compute_parameter(reg.export_world, constants.elasticities.world_export_demand, prices, log_scale)
    β2, ∂β2 = compute_foreign_share(β2, ∂β2, reg.export_world, constants.elasticities.world_export_demand, prices,
                           imports_uk_share_world, reg.totals.imports.world, 1.0, constants.exchange_rates.usd)
    ρ, ∂ρ = compute_parameter(reg.investment, constants.elasticities.investment, prices, log_scale)

    γ, ∂γ = compute_production_parameter(data, prices, log_scale)

    loss_given_default = 0.12 # TODO: This should be in constants

    consts = ParameterConstants(constants.elasticities, loss_given_default, constants.interest_rate)

    vals = Parameters(consts, α, β1, β2, γ, ρ)
    derivs = Parameters(consts, ∂α, ∂β1, ∂β2, ∂γ, ∂ρ)

    return vals, derivs

end

"""
  Compute 1d utility function parameters from a regional demand data frame and the corresponding elasticity.
  Currently missing the tilde parameters for exports.
"""
function compute_parameter(demand::DataFrame, elasticity::Elasticity, prices::DataFrame, log_scale::Bool)

    n = nrow(demand)
    v0 = Vector{Float64}(undef, n)
    m0 = Matrix{Float64}(undef, n, n)

    val = ParamsStruct(similar(v0), similar(v0), similar(v0), similar(v0), nothing)
    grad = ParamsStruct(similar(v0), similar(v0), similar(v0), similar(m0), nothing)

    fun = log_scale ? log_parameters_by_region : parameters_by_region
    
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

    logPf = log_price_index.(elasticity.armington,
                             prices.uk, prices.eu, prices.world,
                             demand.uk, demand.eu, demand.world)
    param = jacobian(ForwardWithPrimal, total_parameters, logPf,
                     Const(demand.agg), Const(elasticity.substitution))

    val.agg .= param.val
    grad.agg .= first(param.derivs)

    return val, grad

end

using Accessors

# Ptilde = 1
# Ex = totals.imports.{eu, world}
# Etilde = Ex / E = totals.imports.{eu, world} / totals.savings
function compute_foreign_share(param::ParamsStruct, dparam::ParamsStruct, demand::DataFrame, elasticity::Elasticity, prices::DataFrame, E::T, Ex::T, PTilde::T, exchange_rate::T) where {T<:Real}

    logPf = Supergrassi.log_price_index.(elasticity.armington, prices.uk, prices.eu, prices.world,
                                         demand.uk, demand.eu, demand.world)

    tilde = gradient(ForwardWithPrimal,
                     log_eu_expenditure_on_uk_exports,
                     logPf,
                     Const(demand.agg),
                     Const(Ex),
                     Const(Ex / E),
                     Const(exchange_rate),
                     Const(PTilde),
                     Const(elasticity.substitution),
                     Const(elasticity.substitution_uk_other))

    @reset param.tilde = [tilde.val]
    @reset dparam.tilde = Vector(first(tilde.derivs))

    return param, dparam

end

"""
  Compute the 2d utility function parameters from regional InputMatrices and the corresponding elasticity.
  Currently missing the jacobian, and the 1d parameters gammaL and gammaH
"""
function compute_production_parameter(data::CleanData, prices::DataFrame, log_scale::Bool)

    n = size(data.industry.regional.input_matrices.agg, 1)
    v0 = Vector{Float64}(undef, n)
    m0 = Matrix{Float64}(undef, n, n)
    t0 = Array{Float64}(undef, n, n, n)

    val = ParamsProduction(similar(v0), similar(v0), similar(v0), similar(v0),
                           similar(v0), data.industry.shock_stdev.val,
                           similar(m0), similar(m0), similar(m0), similar(m0))
    grad = ParamsProduction(similar(m0), similar(m0), similar(v0), similar(v0),
                            similar(m0), similar(v0),
                            similar(m0), similar(m0), similar(m0), similar(t0))

    compute_agg_wages!(data.household, data.constants.elasticities.production)

    fun = log_scale ? log_parameters_by_region : parameters_by_region
    
    for row in 1:n
        for col in 1:n

            param_regional = gradient(ForwardWithPrimal,
                                      fun,
                                      Const(data.constants.elasticities.production.armington),
                                      prices.uk[col],
                                      Const(prices.eu[col]),
                                      Const(prices.world[col]),
                                      Const(data.industry.regional.input_matrices.uk[row, col]),
                                      Const(data.industry.regional.input_matrices.eu[row, col]),
                                      Const(data.industry.regional.input_matrices.world[row, col]))

            val.input_uk[row, col] = param_regional.val[1]
            val.input_eu[row, col] = param_regional.val[2]
            val.input_world[row, col] = param_regional.val[3]

            grad.input_uk[row, col] = param_regional.derivs[2][1]
            grad.input_eu[row, col] = param_regional.derivs[2][2]
            grad.input_world[row, col] = param_regional.derivs[2][3]

        end

        # logPm = log_price_index.(data.constants.elasticities.production.armington,
        #                                      prices.uk, prices.eu, prices.world,
        #                                      m_uk[row,:],
        #                                      m_eu[row,:],
        #                                      m_world[row,:])
        # logPm[isinf.(logPm)] .= 0.0

        # We are missing this from our data struct. Temporarily compute it on the fly.
        tau = (data.industry.tax.products[row] .+ data.industry.tax.production[row]) ./ data.industry.regional.total_use.agg[row]

        jacM = jacobian(ForwardWithPrimal,
                        total_input_parameters,
                        prices.uk,
                        Const(prices.eu),
                        Const(prices.world),
                        Const(data.industry.regional.input_matrices.uk[row,:]),
                        Const(data.industry.regional.input_matrices.eu[row,:]),
                        Const(data.industry.regional.input_matrices.world[row,:]),
                        Const(data.industry.regional.input_matrices.agg[row,:]),
                        Const(data.industry.surplus.val[row]),
                        Const(data.industry.capital.current_year[row]),
                        Const(data.industry.regional.total_use.agg[row]),
                        Const(data.household.payments.agg[row]),
                        Const(data.household.wages.logW[row]),
                        Const(data.constants.elasticities.production),
                        Const(tau),
                        Const(log_scale))

        val.input_agg[row, :] = jacM.val
        grad.input_agg[row,:,:] .= first(jacM.derivs)

        jacH = jacobian(ForwardWithPrimal,
                        total_labor_parameters,
                        prices.uk,
                        Const(prices.eu),
                        Const(prices.world),
                        Const(data.industry.regional.input_matrices.uk[row,:]),
                        Const(data.industry.regional.input_matrices.eu[row,:]),
                        Const(data.industry.regional.input_matrices.world[row,:]),
                        Const(data.industry.regional.input_matrices.agg[row,:]),
                        Const(data.industry.surplus.val[row]),
                        Const(data.industry.capital.current_year[row]),
                        Const(data.industry.regional.total_use.agg[row]),
                        Const(data.household.payments.agg[row]),
                        Const(data.household.wages.logW[row]),
                        Const(data.constants.elasticities.production),
                        Const(tau),
                        Const(log_scale))

        val.input_human[row] = jacH.val
        grad.input_human[row,:] .= first(jacH.derivs)

        jacK = jacobian(ForwardWithPrimal,
                        total_capital_parameters,
                        prices.uk,
                        Const(prices.eu),
                        Const(prices.world),
                        Const(data.industry.regional.input_matrices.uk[row,:]),
                        Const(data.industry.regional.input_matrices.eu[row,:]),
                        Const(data.industry.regional.input_matrices.world[row,:]),
                        Const(data.industry.regional.input_matrices.agg[row,:]),
                        Const(data.industry.surplus.val[row]),
                        Const(data.industry.capital.current_year[row]),
                        Const(data.industry.regional.total_use.agg[row]),
                        Const(data.household.payments.agg[row]),
                        Const(data.household.wages.logW[row]),
                        Const(data.constants.elasticities.production),
                        Const(tau),
                        Const(log_scale))

        val.input_capital[row] = jacK.val
        grad.input_capital[row,:] .= first(jacK.derivs)

        mu = jacobian(ForwardWithPrimal,
                      productivity_shock_mean,
                      Const(data.constants.elasticities.production),
                      prices.uk,
                      Const(prices.eu),
                      Const(prices.world),
                      Const(data.industry.regional.input_matrices.uk[row,:]),
                      Const(data.industry.regional.input_matrices.eu[row,:]),
                      Const(data.industry.regional.input_matrices.world[row,:]),
                      Const(data.industry.regional.input_matrices.agg[row,:]),
                      Const(data.industry.surplus.val[row]),
                      Const(data.industry.capital.current_year[row]),
                      Const(data.industry.regional.total_use.agg[row]),
                      Const(data.household.payments.agg[row]),
                      Const(data.household.wages.logW[row]),
                      Const(tau),
                      Const(row),
                      Const(log_scale))
        
        val.shock_mean[row] = mu.val
        grad.shock_mean[row,:] .= mu.derivs[2]

    end

    return val, grad

end

function compute_agg_wages!(data::HouseholdData, elasticity::Elasticity)

    # ref. https://github.com/UCL/Supergrassi/blob/29510a8c9f50427068a475be01583b544975bd5c/code/matlab/macro_v2/B1_SetupParameters.m#L328-L332

    ξh = elasticity.skill_substitution
    data.wages.logW = ξh/(ξh-1) * log.(
        (data.payments.low.^(1/ξh) .* data.wages.low.^((ξh-1)/ξh))
        + (data.payments.high.^(1/ξh) .* data.wages.high.^((ξh-1)/ξh))
    )
    replace!(data.wages.logW, NaN => 0.0)

end
