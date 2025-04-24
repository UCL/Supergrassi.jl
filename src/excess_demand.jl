function intermediate_goods_price_index(log_price_uk::Vector{T}, zOC::Vector{T},
                                        tau::Vector{T}, mu::Vector{T}, gammaK::Vector{T},
                                        K0::Vector{T}, xi::T) where {T <: Real}

    # Computes the intermediate goods price index, Step 1 of ExcessDemand.m

    pdYBar = Vector{T}(undef, length(log_price_uk))

    for i in axes(log_price_uk, 1)

        pdYBar[i] = intermediate_goods_price_index(log_price_uk[i], zOC[i], tau[i], mu[i], gammaK[i], K0[i], xi)

    end

    return pdYBar

end

function intermediate_goods_price_index(log_price_uk::T, zOC::T, tau::T, mu::T, gammaK::T, K0::T, xi::T) where {T <: Real}

    return ( (1 - tau) * exp(log_price_uk) * mu * K0 * gammaK ^ (1 / (xi - 1) )
             * (1 - (exp(zOC) ) / ( 1 + exp(zOC) ) ) ^ (xi / (1 - xi) ) / (1 - tau) )

end

# Below three functions calculate equations 4.1, 4.3 and 4.4 of main paper that describe the equilibrium
# Equation 4.2 is the log_total_price_index in utility_function_paramerers.jl

function market_clearing_price(log_price_uk::Vector{T}, zOC::Vector{T}, household_expenditure::T,
                               log_price_eu::Vector{T}, log_price_world::Vector{T},
                               consumption::Vector{T}, exports1::Vector{T}, exports2::Vector{T},
                               imports::Vector{T}, inputs::Vector{T},
                               ϵ::T, ζ1::T, ζ2::T, η::T, ξ::T) where T
    
    # Needs:
    # PdYBar: intermediate_goods_price_index()
    # EFd: α_UK * exp(logEF + (1-ϵ_a) * (logPd - logPf) )
    # EX1d: β1_UK * exp(logEX1 + (1-ζ1_a) * (logPd - logPX1) )
    # EX2d: β2_UK * exp(logEX2 + (1-ζ2_a) * (logPd - logPX2) )
    # EId: ρ_UK * exp(logEI + (1 - η_a) * (logPd - logPI) )
    # EMd: γM_UK * exp(logEM + (1 - ξ_a) * (logPd - logPm) ) Note, this is a 2d Matrix
    # DeltaV: From input data

    # logEF: log(α) + log(E) + (1 - ϵ) * (logPf - logPBar)
    # logEX1: log(β1) + log(EX1~) + (1 - ζ1) * (logPX1 - logPX1Bar)
    # logEX2: log(β1) + log(EX2~) + (1 - ζ2) * (logPX2 - logPX1Bar)
    # logEI: log(ρ) + log(PIBar) + log(KS/μI) + (1 - η) * (logPI - logPIBar)
    # logEM: log(γM) + (1 - ξ) * logPm - (1 - ξ) * logTauPdMu + logTauPdYBar

    α_reg = parameters_by_region(ϵ, log_price_uk, log_price_eu, log_price_world,
                                consumption.uk, consumption.eu, consumption.world)
    β1_reg = parameters_by_region(ζ1, log_price_uk, log_price_eu, log_price_world,
                                exports1.uk, exports1.eu, exports1.world)
    β2_reg = parameters_by_region(ζ2, log_price_uk, log_price_eu, log_price_world,
                                exports2.uk, exports2.eu, exports2.world)
    ρ_reg = parameters_by_region(η, log_price_uk, log_price_eu, log_price_world,
                                imports.uk, imports.eu, imports.world)
    
    PdYBar = intermediate_goods_price_index()

    # Household expenditure
    logPf = log_price_by_commodity(α_reg[1], α_reg[2], α_reg[3], log_price_uk, log_price_eu, log_price_world, ϵ)
    logPBar = log_agg_price(α, logPf, ϵ)
    EFd = expenditure_by_region(α_reg[1], α_agg, ϵ, ϵ_a, log_price_uk, household_expenditure, logPf, logPBar)

    # Expenditure on exports from eu
    logPX1 = log(1.0 + τx1) + log_price_by_commodity(β1_reg[1], β1_reg[2], β1_reg[3],
                                                   log_price_uk, log_price_eu, log_price_world, ζ1_a)
    logPX1Bar = log_agg_price(β1, logPX1, ζ1)
    logEX1Tilde = log_expenditure_kernel(β1Tilde, E1Tilde, ζ1Tilde, log(fx_EUR) + logPX1Bar, logP1Tilde)
    EX1d = expenditure_by_region(β1_reg[1], β1_agg, ζ1, ζ1_a, log_price_uk, logEX1Tilde, logPX1, logPX1Bar)

    # Expenditure on exports from rest of world
    logPX2 = log(1.0 + τx2) + log_price_by_commodity(β2_reg[1], β2_reg[2], β2_reg[3],
                                                   log_price_uk, log_price_eu, log_price_world, ζ2_a)
    logPX2Bar = log_agg_price(β2, logPX2, ζ2)
    logEX1Tilde = log_expenditure_kernel(β2Tilde, E2Tilde, ζ2Tilde, log(fx_EUR) + logPX2Bar, logP2Tilde)
    EX2d = expenditure_by_region(β2_reg[1], β2_agg, ζ2, ζ2_a, log_price_uk, E2Tilde, logPX2, logPX2Bar)

    # Investment expenditure
    logPI = log_price_by_commodity(ρ_reg[1], ρ_reg[2], ρ_reg[3], log_price_uk, log_price_eu, log_price_world, η_a)
    logPIBar = log_agg_price(ρ, logPI, η)
    logPITilde = logPIBar + log(KS/μI) # TODO next, KS requires CapitalMarket.m
    EId = expenditure_by_region(ρ_reg[1], ρ_agg, η, η_a, log_price_uk, logPITilde, logPI, logPIBar)

    # # Production intermediates expenditure
    # EMd = zeros(n,n)
    # for j in axes(log_price_uk,1)
    #     γM_reg = parameters_by_region(ξ_a, log_price_uk, log_price_eu, log_price_world,
    #                                   inputs.uk[:,j], inputs.eu[:,j], inputs.world[:,j])
    #     logPM = log_price_by_commodity(γM_reg[1,j], γM_reg[2,j], γM_reg[3,j],
    #                                    log_price_uk[:,j], log_price_eu[:,j], log_price_world[:,j], ξ_a)
    #     logTauPdYBar = ...
    #     EMd[j,:] = expenditure_by_region(γM_reg[1], γM_agg, ξ, ξ_a, log_price_uk, ..., logPM)
    # end
    
    return PdYBar + EFd + EX1d + EX2d + EId + EMd + DeltaV
    
end

# EFd (logEF, E)
function expenditure_by_region(param_region, param_agg, elasticity, elasticity_a, log_price_region, expenditure, logPf, logPBar)

    EF = log_expenditure_kernel(param_agg, expenditure, elasticity, logPf, logPBar)
    EFd = param_region * exp(log(EF) + (1.0 - elasticity_a) * (log_price_region - logPf))
    return EFd
    
end

# LogEF (logPf, logPBar)
function log_expenditure_kernel(param_agg, expenditure, elasticity, logPf, logPBar)

    return log(param_agg) + log(expenditure) + (1.0 - elasticity) * (logPf - logPBar)
    
end

#logPf
function log_price_by_commodity(param_uk, param_eu, param_w, price_uk, price_eu, price_w, elasticity)
    
    return 1.0 / (1.0 - elasticity) * (
        param_uk * price_uk ^ (1 - elasticity)
        + param_eu * price_eu ^ (1 - elasticity)
        + param_w * price_w ^ (1 - elasticity)
    )
    
end

#logPBar (logPf)
function log_agg_price(param_agg, log_final_price_by_commodity, elasticity)

    return 1.0 / (1.0 - elasticity) * log(sum(param_agg * exp((1.0 - elasticity) * log_final_price_by_commodity)))
    
end

function operating_cost()

end

function household_budget_constraint()

end
