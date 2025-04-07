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

function market_clearing_price_index(log_price_uk::Vector{T}, zOC::Vector{T}, household_expenditure::T,
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

    α_uk = parameters_by_region(ϵ, log_price_uk, log_price_eu, log_price_world,
                                consumption.uk, consumption.eu, consumption.world)[1]
    β1_uk = parameters_by_region(ζ1, log_price_uk, log_price_eu, log_price_world,
                                exports1.uk, exports1.eu, exports1.world)[1]
    β2_uk = parameters_by_region(ζ2, log_price_uk, log_price_eu, log_price_world,
                                exports2.uk, exports2.eu, exports2.world)[1]
    ρ_uk = parameters_by_region(η, log_price_uk, log_price_eu, log_price_world,
                                imports.uk, imports.eu, imports.world)[1]
    
    PdYBar = intermediate_goods_price_index()
    EF = α + log(household_expenditure) + (1 - ϵ) * (log_total_price_index - 1)
    EFd = foo(α_uk, ϵ, EF, log_price_uk, log_total_price_index(ϵ_a, log_price_index(), consumption.total))
    EX1d = foo(β1_uk, ζ1, EX1, log_price_uk, log_total_price_index(ζ1_a, log_price_index(), exports1.total))
    EX2d = foo(β2_uk, ζ2, EX2, log_price_uk, log_total_price_index(ζ2_a, log_price_index(), exports2.total))
    EId = foo(ρ_uk, η, EI, log_price_uk, log_total_price_index(η_a, log_price_index(), imports.total))

    EMd = zeros(n,n)
    for j in axes(log_price_uk,1)
        γM_uk = parameters_by_region(ξ_a, log_price_uk, log_price_eu, log_price_world,
                                     inputs.uk[:,j], inputs.eu[:,j], inputs.world[:,j])
        EMd[j,:] = foo(γM_uk, ξ, EM, log_price_uk, log_total_price_index(ξ_a, log_price_index(), input))
    end
    
    return PdYBar + EFd + EX1d + EX2d + EId + EMd + DeltaV
    
end

#
function foo(parameter, elasticity, pBarF, log_price_uk, log_agg_price_index)

    return parameter * exp(log(pBarF) + (1.0 - elasticity) * (log_price_uk - log_agg_price_index))
    
end

#logPf
function log_final_price_by_commodity(param_uk, param_eu, param_w, price_uk, price_eu, price_w, elasticity)
    
    return 1.0 / (1.0 - elasticity) * (param_uk * price_uk ^ (1 - elasticity)
        + param_eu * price_eu ^ (1 - elasticity)
        + param_w * price_w ^ (1 - elasticity)
    )
    
end

#logPBar
function log_final_agg_price(param, log_final_price_by_commodity, elasticity)

    return 1.0 / (1.0 - elasticity) * log(sum(param * exp((1.0 - elasticity) * log_final_price_by_commodity)))
    
end

function operating_cost()

end

function household_budget_constraint()

end
