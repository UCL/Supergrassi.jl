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
    imports::DataFrame
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
