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
function intermediate_goods_price_index(log_price_uk::Vector{T}, zOC::Vector{T},
                                        tau::Vector{T}, mu::Vector{T}, gammaK::Vector{T},
                                        K0::Vector{T}, xi::T) where {T <: Real}

    pdYBar = Vector{T}(undef, length(log_price_uk))
    for i in axes(log_price_uk, 1)
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

    prices = DataFrame([price_uk, price_eu, price_world], ["uk", "eu", "world"])

    elasticity = constants.elasticities

    tau = compute_advalorem_tax(data)

    PdYBar = intermediate_goods_price_index(prices.uk,
                                            data.surplus.val,
                                            tau,
                                            params.production.shock_mean,
                                            params.production.input_capital,
                                            data.capital.current_year,
                                            elasticity.production.substitution)

    imports_uk_share_eu, imports_uk_share_world = compute_imports_shares(constants)
    E1tilde = data.regional.totals.imports.eu / imports_uk_share_eu
    E2tilde = data.regional.totals.imports.world / imports_uk_share_world
    Ptilde = 1.0

    # Household
    # TODO: Could probably wrap each of these blocks in a single function (or two), requires some thought.

    logPf = log_price_index(params.consumption, prices, elasticity.consumption.armington)
    logPBar = log_agg_price_index(params.consumption.agg, logPf, elasticity.consumption.substitution)
    logEf = log_expenditure(params.consumption.agg, household_expenditure,
                            elasticity.consumption.substitution, logPf, logPBar)
    EF_uk = expenditure_by_region(params.consumption, elasticity.consumption, prices, logEf)

    # Exports from eu

    logPX1 = log_price_index(params.export_eu, prices, elasticity.export_eu.armington, constants.export_costs.eu)
    logPX1Bar = log_agg_price_index(params.export_eu.agg, logPX1, elasticity.export_eu.substitution)
    eu_spending = export_spending(elasticity.export_eu.substitution_uk_other, params.export_eu.tilde, logPX1Bar, constants.exchange_rates.eur, E1tilde, Ptilde)
    logEX1 = log_expenditure(params.export_eu.agg, eu_spending, elasticity.export_eu.substitution, logPX1, logPX1Bar)
    EX1_uk = expenditure_by_region(params.export_eu, elasticity.export_eu, prices, logEX1)

    # Exports from rest of world

    logPX2 = log_price_index(params.export_world, prices, elasticity.export_world.armington, constants.export_costs.world)
    logPX2Bar = log_agg_price_index(params.export_world.agg, logPX2, elasticity.export_world.substitution)
    world_spending = export_spending(elasticity.export_world.substitution_uk_other, params.export_world.tilde, logPX2Bar, constants.exchange_rates.usd, E2tilde, Ptilde)
    logEX2 = log_expenditure(params.export_world.agg, world_spending, elasticity.export_world.substitution, logPX2, logPX2Bar)
    EX2_uk = expenditure_by_region(params.export_world, elasticity.export_world, prices, logEX2)

    # Investment

    logPI = log_price_index(params.investment, prices, elasticity.investment.armington)
    logPIBar = log_agg_price_index(params.investment.agg, logPI, elasticity.investment.substitution)
    new_capital_supply = capital_market()
    logEI = log_expenditure(params.investment.agg, logPIBar .+ new_capital_supply,
                            elasticity.investment.substitution, logPI, logPIBar)
    EI_uk = expenditure_by_region(params.investment, elasticity.investment, prices, logEI)

    # Production intermediates

    n = nrow(prices)
    EM_uk = zeros(n)

    for i = 1:n

        TOCTheta = exp(operating_cost[i]) / (1 + exp(operating_cost[i]))
        logTauPdMu = log(1 - tau[i]) + prices.uk[i] + log(params.production.shock_mean[i])
        logTauPdYBar = (logTauPdMu
                        + log(params.production.input_capital[i]) / (elasticity.production.substitution - 1)
                        + elasticity.production.substitution / (1 - elasticity.production.substitution) * log(1 - TOCTheta)
                        + data.capital.current_year[i]
                        )

        logPM = log_price_index(params.production, i, prices, elasticity.production.armington)
        logEM = log_expenditure(params.production.input_agg[i,:], logTauPdYBar, elasticity.production.substitution, logPM, logTauPdMu)
        EM_uk += expenditure_by_region(params.production, i, elasticity.production, prices, logEM)
    end

    # @show PdYBar
    # @show EF_uk
    # @show EX1_uk
    # @show EX2_uk
    # @show EI_uk
    # @show EM_uk
    # @show data.regional.delta_v.agg

    F = PdYBar + EF_uk + EX1_uk + EX2_uk + EI_uk + EM_uk + data.regional.delta_v.agg

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
function log_price_index(param::ParamsStruct, prices::DataFrame, elasticity::T) where {T<:Real}

    return log_price_index(param.uk, param.eu, param.world, prices, elasticity)

end

function log_price_index(param::ParamsProduction, row::Int, prices::DataFrame, elasticity::T) where {T<:Real}

    return log_price_index(param.input_uk[row,:], param.input_eu[row,:], param.input_world[row,:], prices, elasticity)

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
function log_price_index(param::ParamsStruct, prices::DataFrame, elasticity::T, export_cost::T) where {T<:Real}

    logP = log_price_index(param, prices, elasticity)
    logP .+= log.(1 + export_cost)

    return logP

end

function log_price_index(uk::Vector{T}, eu::Vector{T}, world::Vector{T}, prices::DataFrame, elasticity::T) where {T<:Real}

    n = length(uk)
    logP = Vector{T}(undef, n)

    for i in 1:n
        logP[i] = 1 / (1 - elasticity) * log(
            uk[i] * prices.uk[i] ^ (1 - elasticity)
            + eu[i] * prices.eu[i] ^ (1 - elasticity)
            + world[i] * prices.world[i] ^ (1 - elasticity)
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
function expenditure_by_region(param::ParamsStruct, elasticity::Elasticity, prices::DataFrame, log_expenditure::Vector{T}, region::Symbol=:uk) where {T<:Real}

    n = nrow(prices)
    EF = Vector{T}(undef, n)

    logP = log_price_index(param, prices, elasticity.substitution)

    for i = 1:n
        EF[i] = getfield(param, region)[i] * exp(log_expenditure[i] .+ (1.0 .- elasticity.armington) .* (prices[i,region] - logP[i]))
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
        logE[i] = log(param_agg[i]) + log(expenditure) + (1.0 - elasticity) * (logPf[i] - logPBar)
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

    KS = 2.0172
    return KS

end
