using Supergrassi
using DataFrames
using Enzyme

"""
    function compute_all_parameters(data::CleanData, prices::DataFrame, fun::Function = parameters_by_region)

Compute all utility function parameters and their derivativesfrom regional data, elasticities and prices.

For description of the parameters, see Table 2  of Nesheim et al (in prep). Return two [`Parameters`](@ref) structs,
the first one contains values of parameters α, β1, β2, ρ and γ, the second one contains their derivatives with
respect to `prices.uk`.

See the tests in [test_parameters_with_data.jl](https://github.com/UCL/Supergrassi.jl/blob/main/test/test_parameters_with_data.jl) for an example of use.

# Arguments
- `data::CleanData`: Data structure created by [`clean_data`](@ref)
- `prices::DataFrame`: Logarithm of price index equilibrium variable. Has three components (uk, eu, world)
- `fun::Function=parameters_by_region`: Function that computes three parameters by region. Either [`parameters_by_region`](@ref) or [`log_parameters_by_region`](@ref)

See also [`compute_parameter`](@ref), [`compute_foreign_share`](@ref), [`compute_production_parameter`](@ref), [`Parameters`](@ref)
"""
function compute_all_parameters(data::CleanData, prices::DataFrame, log_scale::Bool = false)

    reg = data.industry.regional
    constants = data.constants

    imports_uk_share_eu, imports_uk_share_world = compute_imports_shares(constants)

    α, ∂α = compute_parameter(reg.consumption, constants.elasticities.consumption, prices, log_scale)
    β1, ∂β1 = compute_parameter(reg.export_eu, constants.elasticities.export_eu, prices, log_scale)
    β1, ∂β1 = compute_foreign_share(β1, ∂β1, reg.export_eu, constants.elasticities.export_eu, prices,
                           imports_uk_share_eu, reg.totals.imports.eu, 1.0, constants.exchange_rates.eur)
    β2, ∂β2 = compute_parameter(reg.export_world, constants.elasticities.export_world, prices, log_scale)
    β2, ∂β2 = compute_foreign_share(β2, ∂β2, reg.export_world, constants.elasticities.export_world, prices,
                           imports_uk_share_world, reg.totals.imports.world, 1.0, constants.exchange_rates.usd)
    ρ, ∂ρ = compute_parameter(reg.investment, constants.elasticities.investment, prices, log_scale)

    γ, ∂γ = compute_production_parameter(data, prices, log_scale)

    loss_given_default = 0.12 # TODO: This should be in constants

    consts = ParameterConstants(constants.elasticities, loss_given_default, constants.interest_rate)

    vals = Parameters(consts, α, β1, β2, γ, ρ, log_scale, false)
    derivs = Parameters(consts, ∂α, ∂β1, ∂β2, ∂γ, ∂ρ, log_scale, true)


    return vals, derivs

end

"""
    function compute_parameter(demand::DataFrame, elasticity::Elasticity, prices::DataFrame, fun::Function = parameters_by_region)

Compute 1d utility function parameters and their derivatives from a regional demand DataFrame and the corresponding elasticity.

# Arguments
- `demand::DataFrame`: Demand data disaggregated by region
- `elasticity::Elasticity`: Values for elasticity of substitution
- `prices::DataFrame`: Logarithm of price index equilibrium variable. Disagregated by region.
- `fun::Function = parameters_by_region`: Function that computes three parameters by region. Either [`parameters_by_region`](@ref) or [`log_parameters_by_region`](@ref)
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



"""
    compute_foreign_share(param::ParamsStruct, dparam::ParamsStruct, demand::DataFrame,
                          elasticity::Elasticity, prices::DataFrame, E::T, Ex::T,
                          PTilde::T, exchange_rate::T) where {T<:Real}

Compute the share of foreign expenditure on UK exports, β̃  in the paper by Nesheim et al.

Computes the value and derivative of β̃  and adds them to the structures `param` and `dparam`, respectively. Because
the strucs are immutable, a copy is created. Requires input outside of the regional demand data which is why this has
been split into a different function from [`compute_parameter`](@ref).

# Arguments
- `param::ParamsStruct`: Parameters.
- `dparam::ParamsStruct`: Derivatives of parameters.
- `demand::DataFrame`: Demand data disaggregated by region.
- `elasticity::Elasticity`: Values for elasticity of substitution
- `prices::DataFrame`: Logarithm of price index equilibrium variable. Disagregated by region.
- `E::Real`: Household expenditure
- `Ex::Real`: Foreign expenditure on UK exports
- `PTilde::Real`: The foreign price index.
- `exchange_rate::Real`: The exchange rate between domestic and foreign currency.
"""
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
    function compute_production_parameter(data::CleanData, prices::DataFrame, fun::Function = parameters_by_region)

Compute the 2d utility function parameters γM, γH, γK from regional InputMatrices and the corresponding elasticity.

# Arguments
- `data::CleanData`: Data structure created by [`clean_data`](@ref)
- `prices::DataFrame`: Logarithm of price index equilibrium variable. Has three components (uk, eu, world)
- `fun::Function=parameters_by_region`: Function that computes three parameters by region. Either [`parameters_by_region`](@ref) or [`log_parameters_by_region`](@ref)

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


    logW = compute_wage_index(data.household, data.constants.elasticities.production)
    tau = compute_advalorem_tax(data.industry)
    ih,il = compute_input_by_skill(data.household, data.constants.elasticities.production)
    val.low_skill .= il
    val.high_skill .= ih

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

            val.uk[row, col] = param_regional.val[1]
            val.eu[row, col] = param_regional.val[2]
            val.world[row, col] = param_regional.val[3]

            grad.uk[row, col] = param_regional.derivs[2][1]
            grad.eu[row, col] = param_regional.derivs[2][2]
            grad.world[row, col] = param_regional.derivs[2][3]

        end

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
                        Const(logW[row]),
                        Const(data.constants.elasticities.production),
                        Const(tau[row]),
                        Const(log_scale))

        val.agg[row, :] = jacM.val
        grad.agg[row,:,:] .= first(jacM.derivs)

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
                        Const(logW[row]),
                        Const(data.constants.elasticities.production),
                        Const(tau[row]),
                        Const(log_scale))

        val.human[row] = jacH.val
        grad.human[row,:] .= first(jacH.derivs)

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
                        Const(logW[row]),
                        Const(data.constants.elasticities.production),
                        Const(tau[row]),
                        Const(log_scale))

        val.capital[row] = jacK.val
        grad.capital[row,:] .= first(jacK.derivs)

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
                      Const(logW[row]),
                      Const(tau[row]),
                      Const(row),
                      Const(log_scale))

        val.shock_mean[row] = mu.val
        grad.shock_mean[row,:] .= mu.derivs[2]

    end

    return val, grad

end

"""
    compute_wage_index(data::HouseholdData, elasticity::Elasticity)

Computes aggregate wage index for households based on household data and elasticity of production.

[ref](https://github.com/UCL/Supergrassi/blob/29510a8c9f50427068a475be01583b544975bd5c/code/matlab/macro_v2/B1_SetupParameters.m#L328-L332)

# Arguments
- `data::HouseholdData`: Data structure containing `wages` and `payments` DataFrames [`HouseholdData`](@ref)
- `elasticity::Elasticity`: Elasticity of substitution for production parameters (ξ)
"""
function compute_wage_index(data::HouseholdData, elasticity::Elasticity)

    ξh = elasticity.skill_substitution
    logW = ξh/(ξh-1) * log.(
        (data.payments.low.^(1/ξh) .* data.wages.low.^((ξh-1)/ξh))
        + (data.payments.high.^(1/ξh) .* data.wages.high.^((ξh-1)/ξh))
    )
    replace!(logW, NaN => 0.0)

    return logW

end

"""
    compute_advalorem_tax(data::IndustryData)

Compute the ad valorem tax rate by combining product and production taxes, normalized by total use.

# Arguments
- `data::IndustryData`: Industry data containing tax and regional total use information
"""
function compute_advalorem_tax(data::IndustryData)

    n = nrow(data.tax)
    tau = Vector{Float64}(undef, n)
    
    for i = 1:n
        tau[i] = (data.tax.products[i] + data.tax.production[i]) / data.regional.total_use.agg[i]
    end

    if any(x -> x < 0 && x >= 1, tau)
        throw(ArgumentError("Expected 0 <= τ < 1, got: $tau"))
    end

    return tau

end

"""
    compute_input_by_skill(data::HouseholdData, elasticity::Elasticity)

Compute input shares by skill level (low and high skill) based on household payments and wages.

# Arguments
- `data::HouseholdData`: Household data containing payments and wages information
- `elasticity::Elasticity`: Elasticity of substitution for production parameters
"""
function compute_input_by_skill(data::HouseholdData, elasticity::Elasticity)

    logW = compute_wage_index(data, elasticity)

    input_low_skill = data.payments.low .* (data.wages.low ./ (exp.(logW))) .^ (elasticity.skill_substitution - 1)
    input_high_skill = data.payments.high .* (data.wages.high ./ (exp.(logW))) .^ (elasticity.skill_substitution - 1)

    replace!(input_low_skill, NaN => 0.0)
    replace!(input_high_skill, NaN => 0.0)

    return input_low_skill, input_high_skill

end

"""
    compute_imports_shares(constants::Constants)

Compute the share of UK imports in total imports for EU and world regions.

# Arguments
- `constants::Constants`: Constants containing import and exchange rate information
"""
function compute_imports_shares(constants::Constants)

    imports_uk_share_eu = constants.total_imports_from_uk.eu / (constants.total_imports_from_all_sources.eu
                                                                / constants.exchange_rates.eur )
    imports_uk_share_world = constants.total_imports_from_uk.world / (constants.total_imports_from_all_sources.world
                                                                      / constants.exchange_rates.usd )

    return imports_uk_share_eu, imports_uk_share_world

end
