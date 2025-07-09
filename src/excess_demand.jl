using Supergrassi

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

# Below three functions calculate equations 4.1, 4.3 and 4.4 of main paper that describe the equilibrium
# Equation 4.2 is the log_total_price_index in utility_function_paramerers.jl


"""
function market_clearing_price(prices::DataFrame, elasticities::Elasticities, params::Parameters)

    Calculate the market clearing price, equation 4.1.
"""
function market_clearing_price(prices::DataFrame, operating_cost::DataFrame, elasticities::Elasticities, params::Parameters, data::IndustryData)

    # Needs:
    # PdYBar: intermediate_goods_price_index()
    # EFd:  α_UK  * exp(logEF  + (1 - ϵ_a)  * (logPd - logPf) )
    # EX1d: β1_UK * exp(logEX1 + (1 - ζ1_a) * (logPd - logPX1) )
    # EX2d: β2_UK * exp(logEX2 + (1 - ζ2_a) * (logPd - logPX2) )
    # EId:  ρ_UK  * exp(logEI  + (1 - η_a)  * (logPd - logPI) )
    # EMd:  γM_UK * exp(logEM  + (1 - ξ_a)  * (logPd - logPm) ) Note, this is a 2d Matrix
    # DeltaV: From input data

    # logEF:  log(α)  + log(E)                + (1 - ϵ)  * (logPf - logPBar)
    # logEX1: log(β1) + log(EX1~)             + (1 - ζ1) * (logPX1 - logPX1Bar)
    # logEX2: log(β1) + log(EX2~)             + (1 - ζ2) * (logPX2 - logPX1Bar)
    # logEI:  log(ρ)  + logPIBar + log(KS/μI) + (1 - η)  * (logPI - logPIBar)
    # logEM:  log(γM) + logTauPdYBar          + (1 - ξ)  * (logPm - logTauPdMu)

    # logPf   =              1/(1 - ϵ_a)  * log(sum(N,2))
    # logPX1  = log(1+τx1) + 1/(1 - ζ1_a) * log(sum(N,2))
    # logPX2  = log(1+τx2) + 1/(1 - ζ2_a) * log(sum(N,2))
    # logPI   =              1/(1 - η_a)  * log(sum(N,2))
    # logPm   =              1/(1 - ξ_a)  * log(sum(N,3))

    # logPBar   = 1/(1 - ϵ_a)  * log(sum(α  * exp(1 - ϵ)  * logPf))
    # logPX1Bar = 1/(1 - ζ1_a) * log(sum(β1 * exp(1 - ζ1) * logPX1))
    # logPX2Bar = 1/(1 - ζ2_a) * log(sum(β2 * exp(1 - ζ2) * logPX2))
    # logPIBar  = 1/(1 - η_a)  * log(sum(ρ  * exp(1 - η)  * logPI))

    # logTauPdMu = log(1 - τ) + logPd + log(μ)
    # logTauPdYBar = logTauPdMu + (log(γK) + ξ * log(1 - TOCθ)) / (ξ - 1) + logK0

    
    
    tau = (data.tax.products .+ data.tax.production) ./ data.regional.total_use.agg

    PdYBar = intermediate_goods_price_index(prices.uk,
                                            data.surplus.val,
                                            tau,
                                            params.production.shock_mean,
                                            params.production.input_capital,
                                            data.capital.current_year,
                                            elasticities.production.substitution)

    # Household expenditure

    logPf = log_price_index(params.consumption, prices, elasticity.consumption.armington)
    
    # Expenditure on exports from eu

    logPx1 = log_price_index(params.export_eu, prices, elasticity.eu_export_demand.armington,
                             clean.constants.export_costs.eu)

    # Expenditure on exports from rest of world

    logPx2 = log_price_index(params.export_world, prices, elasticity.world_export_demand.armington,
                             clean.constants.export_costs.world)
    
    # Investment expenditure

    logPI = log_price_index(params.investment, prices, elasticity.investment.armington)

    # Production intermediates expenditure

    for i = 1:n
        logPm = log_price_index(params.production, prices, elasticity.production,armington)

    end

    return PdYBar + EFd + EX1d + EX2d + EId + EMd + DeltaV

end

# EF (logEF, E)
function expenditure_by_region(param::ParamsStruct, elasticity::Elasticity, prices::DataFrame, expenditure::T, region::Symbol=:uk) where {T<:Real}

    n = nrow(prices)
    EF = Vector{T}(undef, n)

    logP = log_price_index(param, prices, elasticity.substitution)
    
    for i = 1:n
        EF[n] = getfield(param, region)[n] * exp(log(expenditure) .+ (1.0 .- elasticity.armington) .* (prices[n,region] - logPf[n]))
    end

    return EF

end

#logP
function log_price_index(param::ParamsStruct, prices::DataFrame, elasticity::T) where {T<:Real}

    logP = 1.0 ./ (1.0 - elasticity) .* (
        param.uk .* prices.uk .^ (1 - elasticity)
        + param.eu .* prices.eu .^ (1 - elasticity)
        + param.world .* prices.world .^ (1 - elasticity)
    )

    return logP

end

function log_price_index(param::ParamsStruct, prices::DataFrame, elasticity::T, export_cost::T) where {T<:Real}

    logP = log_price_index(param, prices, elasticity)
    logP =+ log(1 + export_cost)

    return logP
    
end

#logPBar (logPf)
function log_agg_price(param_agg::Vector{T}, log_final_price_by_commodity::Vector{T}, elasticity::T) where {T<:Real}

    return 1.0 / (1.0 - elasticity) * log(sum(param_agg .* exp.((1.0 - elasticity) .* log_final_price_by_commodity)))

end

# LogEF (logPf, logPBar)
function log_expenditure_kernel(param_agg, expenditure, elasticity, logPf, logPBar)

    return log.(param_agg) .+ log.(expenditure) .+ (1.0 .- elasticity) .* (logPf .- logPBar)

end

function operating_cost()

end

function household_budget_constraint()

end
