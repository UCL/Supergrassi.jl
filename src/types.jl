using DataFrames

struct ExchangeRates
    usd::Float64
    eur::Float64
end

struct TotalImports
    eu::Float64
    world::Float64
end

struct Elasticity
    substitution::Float64
    armington::Float64 # _a
    substitution_uk_other::Union{Float64, Nothing} # ~
    skill_substitution::Union{Float64, Nothing} # _h
end

struct Elasticities
    production::Elasticity # ξ
    world_export_demand::Elasticity # β2
    eu_export_demand::Elasticity # β1
    consumption::Elasticity # α
    investment::Elasticity # ρ
end

struct Constants
    data_year::Int64
    exchange_rates::ExchangeRates
    interest_rate::Float64

    total_imports_from_uk::TotalImports
    total_imports_from_all_sources::TotalImports

    elasticities::Elasticities
end

struct Totals

    savings::Float64 # E
    investments::Float64 # ISum
    imports::TotalImports # EX1, EX2

end


struct InputMatrices

    # These correspond to mValue matrices in the Matlab code

    uk::Matrix{Float64}
    eu::Matrix{Float64}
    world::Matrix{Float64}
    agg::Matrix{Float64}

end

struct AssetsLiabilities

    current_year::DataFrame
    next_year::DataFrame

end

struct RegionalData

    total_use::DataFrame # data.yValue
    consumption::DataFrame # data.fValue
    delta_v::DataFrame # data.deltaVValue
    export_eu::DataFrame # data.x1Value
    export_world::DataFrame # data.x2Value
    investment::DataFrame  # data.IValue Also called "payments to capital" in the code
    input_matrices::InputMatrices # data.mValue
    totals::Totals # data.{E, ISum, EX1, EX2}

end

struct HouseholdData

    income::DataFrame
    payments::DataFrame # data.{hValueLO, hValueHI, hValue}
    hours::DataFrame
    wages::DataFrame

end

struct IndustryData

    depreciation::DataFrame
    tax::DataFrame # data.taxValue{1,2}
    capital::DataFrame # data.k{0,1}
    surplus::DataFrame # data.kValue
    shock_stdev::DataFrame # sigmaBar
    assets_liabilities::AssetsLiabilities
    regional::RegionalData

end

struct CleanData

    household::HouseholdData
    industry::IndustryData
    constants::Constants

end

struct ParamsStruct

    uk::Vector{Float64}
    eu::Vector{Float64}
    world::Vector{Float64}
    agg::Array{Float64}
    tilde::Union{Vector{Float64}, Nothing} # Optional, used for some parameters that have a tilde version

end

struct ParamsProduction

    input_human::Array{Float64}          # gamma_hi
    input_capital::Array{Float64}        # gamma_ki
    input_low_skill::Vector{Float64}      # gamma_Li
    input_high_skill::Vector{Float64}     # gamma_Hi
    shock_mean::Array{Float64}           # mu
    shock_stddev::Array{Float64}         # sigmaBar

    input_uk::Matrix{Float64}             # gamma_mdij
    input_eu::Matrix{Float64}             # gamma_meuij
    input_world::Matrix{Float64}          # gamma_mwij
    input_agg::Array{Float64}             # gamma_mij

end

struct ParameterConstants

    elasticities::Elasticities
    loss_given_default::Float64
    risk_free_interest_rate::Float64

end

struct Parameters

    constants::ParameterConstants

    consumption::ParamsStruct # alpha
    export_eu::ParamsStruct # beta1
    export_world::ParamsStruct # beta2
    production::ParamsProduction # gamma
    investment::ParamsStruct # rho

end
