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
    armington::Float64
    substitution_uk_other::Union{Float64, Nothing}
    skill_substitution::Union{Float64, Nothing}
end

struct Elasiticities
    production::Elasticity
    world_export_demand::Elasticity
    eu_export_demand::Elasticity
    consumption::Elasticity
    investment::Elasticity
end

struct Constants
    data_year::Int64
    exchange_rates::ExchangeRates
    interest_rate::Float64

    total_imports_from_uk::TotalImports
    total_imports_from_all_sources::TotalImports

    elasticities::Elasiticities
end

struct Totals

    savings::Float64 # E
    investments::Float64 # ISum
    imports::TotalImports # EX1, EX2

end


struct InputMatrices

    # These correspond to mValue matrices in the Matlab code

    uk::DataFrame
    eu::DataFrame
    world::DataFrame
    agg::DataFrame

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
    income_share::DataFrame
    payments::DataFrame
    hours::DataFrame
    wages::DataFrame

end

struct IndustryData

    depreciation::DataFrame
    tax::DataFrame
    capital::DataFrame
    surplus::DataFrame
    shock_stdev::DataFrame
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
    agg::Vector{Float64}
    tilde::Union{Vector{Float64}, Nothing} # Optional, used for some parameters that have a tilde version

end

struct ParamsProduction

    var1::Vector{Float64}
    var2::Vector{Float64}
    var3::Vector{Float64}
    var4::Vector{Float64}
    var5::Vector{Float64}

    mat1::Matrix{Float64}
    mat2::Matrix{Float64}
    mat3::Matrix{Float64}
    mat4::Matrix{Float64}

end

struct ParameterConstants

    elasiticities::Elasiticities
    loss_given_default::Float64
    risk_free_interest_rate::Float64
   
end

struct Parameters

    constants::ParameterConstants

    consumption::ParamsStruct
    export_eu::ParamsStruct
    export_world::ParamsStruct

    intermediate_input::ParamsProduction

    investment::ParamsStruct

end
