using DataFrames

"""
    struct ExchangeRates

Stores exchange rates from GBP to other currencies
"""
struct ExchangeRates
    usd::Float64
    eur::Float64
end

"""
    struct ForeignRegionalValues

Stores eu and rest of world components of a scalar quantity.
"""
struct ForeignRegionalValues
    eu::Float64
    world::Float64
end

"""
    struct Elasticity

Stores components of an elasticity constant.

Names for reference in paper/matlab code for e.g. production elasticity ξ:

- substitution: ξ
- armington: ξ_a
- substitution_uk_other: ~ξ
- skill_substitution: ξ_h
"""
struct Elasticity
    substitution::Float64
    armington::Float64
    substitution_uk_other::Union{Float64, Nothing}
    skill_substitution::Union{Float64, Nothing}
end

"""
    struct Elasticities

Stores sets of elasticity constants.

Names for reference in paper/matlab code:

- consumption: α
- eu_export_demand: β1
- world_export_demand: β2
- investment: ρ
- production: ξ
"""
struct Elasticities
    production::Elasticity
    world_export_demand::Elasticity
    eu_export_demand::Elasticity
    consumption::Elasticity
    investment::Elasticity
end

"""
    struct Constants

Stores constants read from settings file.
"""
struct Constants
    data_year::Int64
    exchange_rates::ExchangeRates
    interest_rate::Float64

    total_imports_from_uk::ForeignRegionalValues
    total_imports_from_all_sources::ForeignRegionalValues

    import_tariffs::ForeignRegionalValues
    export_costs::ForeignRegionalValues

    elasticities::Elasticities
end

"""
    struct Totals

Stores sums of quantities over industries before renormalisation.

Names for reference with the paper/matlab code:

- savings: E
- investments: ISum
- imports: EX1, EX2
"""
struct Totals

    savings::Float64
    investments::Float64
    imports::ForeignRegionalValues

end

"""
    struct InputMatrices

Stores intermediate input production matrices .

In the matlab code these matrices are called `mValue_*`
"""
struct InputMatrices

    uk::Matrix{Float64}
    eu::Matrix{Float64}
    world::Matrix{Float64}
    agg::Matrix{Float64}

end

"""
    struct AssetsLiabilities

In the matlab code these are called `assets0` and `assets1`
"""
struct AssetsLiabilities

    current_year::DataFrame
    next_year::DataFrame

end

"""
    struct RegionalData

Stores industries data that is split between uk/eu/rest of the world.

The data is vectors with one value per industry stored in a DataFrame. Names for reference with
the matlab code/paper:

- total_use: data.yValue / y
- consumption:data.fValue / f
- delta_v: data.deltaVValue / Δv
- export_eu: data.x1Value / x1
- export_world: data.x2Value / x2
- investment:: data.IValue / I
- input_matrices:: data.mValue / m
- totals: data.{E, ISum, EX1, EX2}

"""
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

"""
    struct HouseholdData

Stores data on households

Names reference for matlab code/paper

- income: data.income_{lo, hi}
- payments: data.{hValueLO, hValueHI, hValue}
- hours: data.h{LO, HI}
- wages: data.w{LO, HI}

"""
struct HouseholdData

    income::DataFrame
    payments::DataFrame
    hours::DataFrame
    wages::DataFrame

end

"""
    struct IndustryData

Stores data on industries

Names reference for matlab code/paper

- depreciation: data.depreciation
- tax: data.taxValue{1,2}
- capital: data.k{0,1}
- surplus: data.kValue
- shock_stdev: data.sigmaBar
- assets_liabilities: [AssetsLiabilities](@ref)
- regional: [RegionalData](@ref)

"""
struct IndustryData

    depreciation::DataFrame
    tax::DataFrame
    capital::DataFrame
    surplus::DataFrame
    shock_stdev::DataFrame
    assets_liabilities::AssetsLiabilities
    regional::RegionalData

end

"""
    struct CleanData

Top level structure for cleaned up input data.

Contains

- [HouseholdData](@ref)
- [IndustryData](@ref)
- [Constants](@ref)
"""
struct CleanData

    household::HouseholdData
    industry::IndustryData
    constants::Constants

end

"""
    struct ParamsStruct

Stores parameters that are split between uk/eu/rest of the world.

Contains an optional field "tilde" which, for the export parameters β,
is the share of foreign expenditure on UK exports.
"""
struct ParamsStruct

    uk::Vector{Float64}
    eu::Vector{Float64}
    world::Vector{Float64}
    agg::Array{Float64}
    tilde::Union{Vector{Float64}, Nothing} # Optional, used for some parameters that have a tilde version

end

"""
    struct ParamsProduction

Stores firms intermediate production parameters

These require a different structure from [ParamsStruct](@ref) because firms trade with all other firms
and the parameter is a matrix that contains a value for each combination of firms.

Names reference
- `input_human::Vector` : γ_h
- `input_capital::Vector` :  γ_k
- `input_low_skill::Vector` : γ_L
- `input_high_skill::Vector` : γ_H
- `shock_mean::Vector` : μ
- `shock_stddev::Vector` : ̄σ
- `input_uk::Matrix` : γ_Md
- `input_eu::Matrix` : γ_Meu
- `input_world::Matrix` : γ_Mw
- `input_agg::Matrix` : γ_M
"""
struct ParamsProduction

    input_human::Array{Float64}
    input_capital::Array{Float64}
    input_low_skill::Vector{Float64}
    input_high_skill::Vector{Float64}
    shock_mean::Vector{Float64}
    shock_stddev::Vector{Float64}

    input_uk::Matrix{Float64}
    input_eu::Matrix{Float64}
    input_world::Matrix{Float64}
    input_agg::Array{Float64}

end

"""
    struct ParameterConstants

Stores constant values that are stored in parameters but not updated
"""
struct ParameterConstants

    elasticities::Elasticities
    loss_given_default::Float64
    risk_free_interest_rate::Float64

end

"""
    struct Parameters

Main structure for parameters

Contains

- `constants::ParameterConstants` : constant values
- `consumption::ParamsStruct` : α
- `export_eu::ParamsStruct` : β1
- `export_world::ParamsStruct` : β2
- `production::ParamsProduction` : γ
- `investment::ParamsStruct` : ρ
"""
struct Parameters

    constants::ParameterConstants

    consumption::ParamsStruct
    export_eu::ParamsStruct
    export_world::ParamsStruct
    production::ParamsProduction
    investment::ParamsStruct

end
