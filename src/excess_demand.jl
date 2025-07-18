using Supergrassi

"""
    intermediate_goods_price_index(log_price_uk::Vector{T}, zOC::Vector{T}, tau::Vector{T}, mu::Vector{T}, gammaK::Vector{T}, K0::Vector{T}, xi::T) where {T <: Real}

Computes the intermediate goods price index for multiple industries, Step 1 of ExcessDemand.m.

# Arguments
- `log_price_uk::Vector{T}`: Logarithm of UK prices
- `zOC::Vector{T}`: Operating cost parameters
- `tau::Vector{T}`: Ad Valorem tax rates
- `mu::Vector{T}`: Productivity shock means
- `gammaK::Vector{T}`: Capital input parameters
- `K0::Vector{T}`: Current year capital values
- `xi::T`: Production substitution elasticity
"""
function intermediate_goods_price_index(log_price_uk::Vector{T}, zOC::Vector{T}, tau::Vector{T}, mu::Vector{T},
                                        gammaK::Vector{T}, K0::Vector{T}, xi::T) where {T <: Real}

    pdYBar = Vector{T}(undef, length(log_price_uk))
    for i in 1:length(log_price_uk)
        pdYBar[i] = intermediate_goods_price_index(log_price_uk[i], zOC[i], tau[i], mu[i], gammaK[i], K0[i], xi)
    end
    return pdYBar

    # return [intermediate_goods_price_index(log_price_uk[i], zOC[i], tau[i], mu[i], gammaK[i], K0[i], xi) for i in eachindex(log_price_uk, zOC, tau, mu, gammaK, K0)]

end

"""
    intermediate_goods_price_index(log_price_uk::T, zOC::T, tau::T, mu::T, gammaK::T, K0::T, xi::T) where {T <: Real}

Computes the intermediate goods price index for a single industry, Step 1 of ExcessDemand.m.

# Arguments
- `log_price_uk::T`: Logarithm of UK price
- `zOC::T`: Operating cost parameter
- `tau::T`: Ad Valorem tax rate
- `mu::T`: Productivity shock mean
- `gammaK::T`: Capital input parameter
- `K0::T`: Current year capital value
- `xi::T`: Production substitution elasticity
"""
function intermediate_goods_price_index(log_price_uk::T, zOC::T, tau::T, mu::T, gammaK::T, K0::T, xi::T) where {T <: Real}

    return ( (1 - tau) * exp(log_price_uk) * mu * K0 * gammaK ^ (1 / (xi - 1) )
             * (1 - (exp(zOC) ) / ( 1 + exp(zOC) ) ) ^ (xi / (1 - xi) ) / (1 - tau) )

end

# Below three functions calculate equations 4.1, 4.3 and 4.4 of main paper that describe the equilibrium
# Equation 4.2 is the log_total_price_index in utility_function_paramerers.jl


"""
    function market_clearing_price(prices::DataFrame, operating_cost::Vector{T}, household_expenditure::T,
                               elasticities::Elasticities, params::Parameters, data::IndustryData) where {T <: Real}


Calculate the market clearing price, equation 4.1.

# Arguments

- `prices::DataFrame` : prices for uk, eu, rest of world (pd, peu, pw)
- `operating_cost::Vector{Real}` : (zOC)
- `household_expenditure::Real` : (E)
- `elasticities::Elasticities`
- `params::Parameters`
- `data::IndustryData`
"""
function market_clearing_price(price_uk::Vector{T}, operating_cost::Vector{T}, household_expenditure::T,
                               price_eu::Vector{T}, price_world::Vector{T},
                               params::Parameters, data::IndustryData, constants::Constants) where {T <: Real}

    # Needs:

    # PdYBar: intermediate_goods_price_index()

    # expenditure (by region)
    ## EFd:  α_UK  * exp(logEF  + (1 - ϵ_a)  * (logPd - logPf) )
    ## EX1d: β1_UK * exp(logEX1 + (1 - ζ1_a) * (logPd - logPX1) )
    ## EX2d: β2_UK * exp(logEX2 + (1 - ζ2_a) * (logPd - logPX2) )
    ## EId:  ρ_UK  * exp(logEI  + (1 - η_a)  * (logPd - logPI) )
    ## EMd:  γM_UK * exp(logEM  + (1 - ξ_a)  * (logPd - logPm) ) Note, this is a 2d Matrix
    ## DeltaV: From input data

    # expenditure (by commodity)
    ## logEF:  log(α)  + log(E)                + (1 - ϵ)  * (logPf - logPBar)
    ## logEX1: log(β1) + log(EX1~)             + (1 - ζ1) * (logPX1 - logPX1Bar)
    ## logEX2: log(β1) + log(EX2~)             + (1 - ζ2) * (logPX2 - logPX1Bar)
    ## logEI:  log(ρ)  + logPIBar + log(KS/μI) + (1 - η)  * (logPI - logPIBar)
    ## logEM:  log(γM) + logTauPdYBar          + (1 - ξ)  * (logPm - logTauPdMu)

    # log price index (by commodity)
    ## logPf   =              1/(1 - ϵ_a)  * log(sum(N,2))
    ## logPX1  = log(1+τx1) + 1/(1 - ζ1_a) * log(sum(N,2))
    ## logPX2  = log(1+τx2) + 1/(1 - ζ2_a) * log(sum(N,2))
    ## logPI   =              1/(1 - η_a)  * log(sum(N,2))
    ## logPm   =              1/(1 - ξ_a)  * log(sum(N,3))

    # log price index (aggregate)
    ## logPBar   = 1/(1 - ϵ_a)  * log(sum(α  * exp(1 - ϵ)   * logPf))
    ## logPX1Bar = 1/(1 - ζ1_a) * log(sum(β1 * exp(1 - ζ1) * logPX1))
    ## logPX2Bar = 1/(1 - ζ2_a) * log(sum(β2 * exp(1 - ζ2) * logPX2))
    ## logPIBar  = 1/(1 - η_a)  * log(sum(ρ  * exp(1 - η)  * logPI))

    # logTauPdMu = log(1 - τ) + logPd + log(μ)
    # logTauPdYBar = logTauPdMu + (log(γK) + ξ * log(1 - TOCθ)) / (ξ - 1) + logK0

    elasticity = constants.elasticities

    tau = compute_advalorem_tax(data)

    # PdYBar = intermediate_goods_price_index(price_uk,
    #                                         data.surplus.val,
    #                                         tau,
    #                                         params.production.shock_mean,
    #                                         params.production.capital,
    #                                         data.capital.current_year,
    #                                         elasticity.production.substitution)

    imports_uk_share_eu, imports_uk_share_world = compute_imports_shares(constants)
    E1tilde = data.regional.totals.imports.eu / imports_uk_share_eu
    E2tilde = data.regional.totals.imports.world / imports_uk_share_world
    Ptilde = 1.0

    n = length(price_uk)
    keys = (:consumption, :export_eu, :export_world, :investment)
    logP = NamedTuple{keys}((zeros(n), zeros(n), zeros(n), zeros(n)))
    taux = NamedTuple{keys}((0.0, constants.export_costs.eu, constants.export_costs.world, 0.0))
    logPBar = Vector{T}(undef, length(keys))

    # Compute logP and logPBar for paramters in keys. Store values in NamedTuples.

    for (i,key) in enumerate(keys)
        param = getfield(params, key)
        el = getfield(elasticity, key)
        getfield(logP, key) .= log_price_index(param.uk, param.eu, param.world, price_uk, price_eu, price_world, el.armington, getfield(taux, key))
        logPBar[i] = log_agg_price_index(param.agg, getfield(logP, key), el.substitution)
    end
    logPBar = NamedTuple{keys}(logPBar)

    # Household

    logEf = log_expenditure(params.consumption.agg, log(household_expenditure), elasticity.consumption.substitution,
                            logP.consumption, logPBar.consumption)
    EF_uk = expenditure_by_region(params.consumption.uk, params.consumption.eu, params.consumption.world,
                                  price_uk, price_eu, price_world, logEf, elasticity.consumption)

    # Exports from eu

    eu_spending = export_spending(elasticity.export_eu.substitution_uk_other, params.export_eu.tilde, logPBar.export_eu,
                                  constants.exchange_rates.eur, E1tilde, Ptilde)
    logEX1 = log_expenditure(params.export_eu.agg, eu_spending, elasticity.export_eu.substitution,
                             logP.export_eu, logPBar.export_eu)
    EX1_uk = expenditure_by_region(params.export_eu.uk, params.export_eu.eu, params.export_eu.world,
                                   price_uk, price_eu, price_world, logEX1, elasticity.export_eu)
    # Exports from rest of world

    world_spending = export_spending(elasticity.export_world.substitution_uk_other, params.export_world.tilde,
                                     logPBar.export_world, constants.exchange_rates.usd, E2tilde, Ptilde)
    logEX2 = log_expenditure(params.export_world.agg, world_spending, elasticity.export_world.substitution,
                             logP.export_world, logPBar.export_world)
    EX2_uk = expenditure_by_region(params.export_world.uk, params.export_world.eu, params.export_world.world,
                                   price_uk, price_eu, price_world, logEX2, elasticity.export_world)

    # Investment

    new_capital_supply = capital_market()
    muI = compute_muI(data, elasticity.investment)
    logEI = log_expenditure(params.investment.agg, logPBar.investment .+ log(new_capital_supply ./ muI),
                            elasticity.investment.substitution, logP.investment, logPBar.investment)
    EI_uk = expenditure_by_region(params.investment.uk, params.investment.eu, params.investment.world,
                                  price_uk, price_eu, price_world, logEI, elasticity.investment)

    # Production intermediates

    EM_uk = zeros(n)
    pdYBar = Vector{T}(undef, n)

    for i = 1:n

        TOCTheta = exp(operating_cost[i]) / (1 + exp(operating_cost[i]))
        logTauPdMu = log(1 - tau[i]) + log(price_uk[i]) + log(params.production.shock_mean[i])

        # TODO: This is already computed as part of intermediate_goods_price_index. Refactor and reuse.
        logTauPdYBar = (logTauPdMu
                        + log(params.production.capital[i]) / (elasticity.production.substitution - 1)
                        + elasticity.production.substitution / (1 - elasticity.production.substitution) * log(1 - TOCTheta)
                        + log(data.capital.current_year[i])
                        )
        pdYBar[i] = exp(logTauPdYBar - log(1 - tau[i]))
        
        logPM = log_price_index(params.production.uk[i,:], params.production.eu[i,:],
                                params.production.world[i,:],
                                price_uk, price_eu, price_world, elasticity.production.armington)

        logEM = log_expenditure(params.production.agg[i,:], logTauPdYBar, elasticity.production.substitution,
                                logPM, logTauPdMu)
        EM_uk += expenditure_by_region(params.production.uk[i,:], params.production.eu[i,:],
                                       params.production.world[i,:],
                                       price_uk, price_eu, price_world, logEM, elasticity.production)
    end

    return pdYBar, EF_uk, EX1_uk, EX2_uk, EI_uk, EM_uk

end

"""
    function log_price_index(param::ParamsStruct, prices::DataFrame, elasticity::T, export_cost::T) where {T<:Real}

Compute the log of final price index (logP) for export commodities, taking into account the export cost

# Arguments

- `param::ParamsStruct` : parameters by region, one of [α, β1, β2, ρ, γ], see [`ParamsStruct`](@ref)
- `prices::DataFrame` : prices for uk, eu, rest of world
- `elasticity::Real : armington elasticity
- `export_cost::Real` : export cost parameter, see [`Constants`](@ref)

# See also

[`ParamsStruct`](@ref), [`Constants`](@ref)
"""
function log_price_index(param_uk::Vector{T}, param_eu::Vector{T}, param_world::Vector{T}, price_uk::Vector{T}, price_eu::Vector{T}, price_world::Vector{T}, elasticity::T, export_cost::T) where {T<:Real}

    logP = log_price_index(param_uk, param_eu, param_world, price_uk, price_eu, price_world, elasticity)
    logP .+= log(1 + export_cost)

    return logP

end

"""
    function log_price_index(param::ParamsStruct, prices::DataFrame, elasticity::T) where {T<:Real}

Compute the log of final price index (logP) for each commodity

# Arguments

- `param::ParamsStruct` : parameters by region, one of [α, β1, β2, ρ, γ], see [`ParamsStruct`](@ref)
- `prices::DataFrame` : prices for uk, eu, rest of world
- `elasticity::Real : armington elasticity

# See also
[`ParamsStruct`](@ref)
"""
function log_price_index(param_uk::Vector{T}, param_eu::Vector{T}, param_world::Vector{T}, price_uk::Vector{T}, price_eu::Vector{T}, price_world::Vector{T}, elasticity::T) where {T<:Real}

    n = length(price_uk)
    logP = Vector{T}(undef, n)

    for i in 1:n
        logP[i] = 1 / (1 - elasticity) * log(
            param_uk[i] * price_uk[i] ^ (1 - elasticity)
            + param_eu[i] * price_eu[i] ^ (1 - elasticity)
            + param_world[i] * price_world[i] ^ (1 - elasticity)
        )
    end
    replace!(logP, Inf => 0.0)

    return logP

end

"""
    function log_agg_price_index(param_agg::Vector{T}, log_price_index::Vector{T}, elasticity::T) where {T<:Real}

Compute the log of aggregate price index (logPBar)

# Arguments

- `param_agg::Vector{Real}` : aggregate parameter, one of [α, β1, β2, ρ, γ]
- `log_price_index::Vector{Real}` : log of price index (by commodity, logP), see [`log_price_index`](@ref)
- `elasticity::Real : substitution elasticity

# Output

- `logPBar::Real`

# See also
[`ParamsStruct`](@ref), [`log_price_index`](@ref)
"""
function log_agg_price_index(param_agg::Vector{T}, log_price_index::Vector{T}, elasticity::T) where {T<:Real}

    n = length(param_agg)
    N = Vector{T}(undef, n)

    for i in 1:n
        N[i] = param_agg[i] * exp((1.0 - elasticity) * log_price_index[i])
    end

    logPBar = 1.0 / (1.0 - elasticity) * log(sum(N))
    return logPBar

end

# EF (logEF, E)
"""
    function expenditure_by_region(param::ParamsStruct, elasticity::Elasticity, prices::DataFrame, expenditure::Vector{T}, region::Symbol=:uk) where {T<:Real}

Compute expenditure for region `region`. Called [EF, EX1, EX2, EI, EM] in Matlab code.

# Arguments
- `param::ParamsStruct` : parameter with components [uk, eu, world]. One of [α, β1, β2, ρ], see [`ParamsStruct`](@ref)
- `elasticity::Elasticity` : elasticity struct corresponding to parameter. Must contain substitution and armington.
- `prices::DataFrame` : prices for uk, eu, rest of world
- `expenditure::Vector{T}` : expenditure by commodity, see [`log_expenditure`](@ref)
"""
#function expenditure_by_region(param::ParamsStruct, elasticity::Elasticity, prices::DataFrame, log_expenditure::Vector{T}, region::Symbol) where {T<:Real}
function expenditure_by_region(param_uk::Vector{T}, param_eu::Vector{T}, param_world::Vector{T}, price_uk::Vector{T}, price_eu::Vector{T}, price_world::Vector{T}, log_expenditure::Vector{T}, elasticity::Elasticity) where {T<:Real}

    n = length(price_uk)
    EF = Vector{T}(undef, n)

    logP = log_price_index(param_uk, param_eu, param_world, price_uk, price_eu, price_world, elasticity.armington)

    for i = 1:n
        EF[i] = param_uk[i] * exp(log_expenditure[i] + (1.0 - elasticity.armington) * (log(price_uk[i]) - logP[i]))
    end

    return EF

end

# EF (logEF, E)
function expenditure_by_region(param::ParamsProduction, row::Int, elasticity::Elasticity, prices::DataFrame, log_expenditure::Vector{T}, region::Symbol=:uk) where {T<:Real}

    n = nrow(prices)
    EF = Vector{T}(undef, n)

    logP = log_price_index(param, row, prices, elasticity.substitution)

    fieldname = Symbol("input_", region)
    for i = 1:n
        EF[i] = getfield(param, fieldname)[row,i] * exp(log_expenditure[i] .+ (1.0 .- elasticity.armington) .* (prices[i,region] - logP[i]))
    end

    return EF

end


# LogEF (logPf, logPBar)
function log_expenditure(param_agg::Vector{T}, expenditure::T, elasticity::T, logPf::Vector{T}, logPBar::T) where {T <: Real}

    n = length(param_agg)
    logE = Vector{T}(undef, n)

    for i in 1:n
        logE[i] = log(param_agg[i]) + expenditure + (1.0 - elasticity) * (logPf[i] - logPBar)
    end

    return logE

end

"""
    function export_spending(elasticity, tilde_param, logPXBar, exchange_rate, ETilde, PTilde)

Computes foreign spending on UK exports (logEXTilde).

# Arguments
- `elasticity::Real` : substitution elasticity
- `tilde_param::Vector{Real}` : share of foreign expenditure on UK exports, see [`ParamsStruct`](@ref). Vector of size 1.
- `logPXBar::Real` : log of aggregate price index, see [`log_agg_price_index`](@ref)
- `exchange_rate::Real` : exhange rate to foreign currency, see [`Constants`](@ref)
- `Etilde::Real`: EU expenditure on UK exports, see [`Constants`](@ref)
- `Ptilde::Real` : UK export price index, see [`Constants`](@ref)

"""
function export_spending(elasticity::T, tilde_param::Vector{T}, logPXBar::T, exchange_rate::T, ETilde::T, PTilde::T) where {T <: Real}

    logPXTilde = log(exchange_rate) .+ logPXBar
    logEXTilde = first(tilde_param) .+ log(ETilde) .+ (1 - elasticity) .* (logPXTilde .- log(PTilde))
    return logEXTilde

end

"""
    function capital_market()

Dummy function for capital market. Returns a number.
"""
function capital_market()

    KS = 2.0712
    return KS

end


"""
    function compute_muI

Compute μ_I parameter from industry data. Ref [L429 of B1_SetupParameters.m](https://github.com/UCL/Supergrassi/blob/2f147384f02f8eef85e3cb59c73fd64ebfc19f82/code/matlab/macro_v2/B1_SetupParameters.m#L425)
"""
function compute_muI(data::IndustryData, elasticity::Elasticity)

    DeltaK = sum(data.capital.next_year) - sum((1 .- data.depreciation.val) .* data.capital.current_year)
    muI = sum((data.regional.totals.investments .* data.regional.investment.agg ./ DeltaK) .^ (1 / elasticity.substitution))

    return muI

end
