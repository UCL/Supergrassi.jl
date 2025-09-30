# Supergrassi

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/UCL/Supergrassi.jl/actions/workflows/UnitTests.yml/badge.svg?branch=main)](https://github.com/UCL/Supergrassi.jl/actions/workflows/UnitTests.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/UCL/Supergrassi.jl/graph/badge.svg?token=LU852FO3FP)](https://codecov.io/gh/UCL/Supergrassi.jl)

Multi-sector dynamic macroeconomics model with debt and default

## Read data into Julia

```julia
using Supergrassi

settings_path = create_filepath("config/settings.yml")
settings = read_settings(settings_path)
filepaths = check_file_availability(settings)
data = read_data(filepaths, settings)

```

```mermaid
classDiagram
    %% ========= Core data =========
    class CleanData {
        household : HouseholdData
        industry : IndustryData
        constants : Constants
    }

    class HouseholdData {
        income : DataFrame
        payments : DataFrame
        hours : DataFrame
        wages : DataFrame
    }

    class IndustryData {
        depreciation : DataFrame
        tax : DataFrame
        capital : DataFrame
        surplus : DataFrame
        shock_stdev : DataFrame
        assets_liabilities : AssetsLiabilities
        regional : RegionalData
    }

    class AssetsLiabilities {
        current_year : DataFrame
        next_year : DataFrame
    }

    class RegionalData {
        total_use : DataFrame        %% y
        consumption : DataFrame      %% f
        delta_v : DataFrame          %% Δv
        export_eu : DataFrame        %% x1
        export_world : DataFrame     %% x2
        investment : DataFrame       %% I
        input_matrices : InputMatrices
        totals : Totals
    }

    class InputMatrices {
        uk : Matrix{Float64}
        eu : Matrix{Float64}
        world : Matrix{Float64}
        agg : Matrix{Float64}
    }

    class Totals {
        expenditure : Float64        %% E
        investments : Float64        %% ISum
        imports : ForeignRegionalValues  %% (EX1, EX2)
    }

    %% ========= Constants & helpers =========
    class Constants {
        data_year : Int64
        exchange_rates : ExchangeRates
        interest_rate : Float64
        total_imports_from_uk : ForeignRegionalValues
        total_imports_from_all_sources : ForeignRegionalValues
        import_tariffs : ForeignRegionalValues
        export_costs : ForeignRegionalValues
        elasticities : Elasticities
        loss_given_default : Float64
        number_of_industries : Int64
    }

    class ExchangeRates {
        usd : Float64
        eur : Float64
    }

    class ForeignRegionalValues {
        eu : Float64
        world : Float64
    }

    class Elasticity {
        substitution : Float64            %% ξ
        armington : Float64               %% ξ_a
        substitution_uk_other : Float64?  %% ~ξ
        skill_substitution : Float64?     %% ξ_h
    }

    class Elasticities {
        production : Elasticity    %% ξ
        export_world : Elasticity  %% β2
        export_eu : Elasticity     %% β1
        consumption : Elasticity   %% α
        investment : Elasticity    %% ρ
    }

    %% ========= Parameters =========
    class ParamsStruct {
        uk : Vector{Float64}
        eu : Vector{Float64}
        world : Vector{Float64}
        agg : Vector{Float64}
        tilde : Vector{Float64}?    %% optional
    }

    class ParamsProduction {
        human : Vector{Float64}         %% γ_h
        capital : Vector{Float64}       %% γ_k
        low_skill : Vector{Float64}     %% γ_L
        high_skill : Vector{Float64}    %% γ_H
        shock_mean : Vector{Float64}    %% μ
        shock_stddev : Vector{Float64}  %% σ̄
        uk : Matrix{Float64}            %% γ_Md
        eu : Matrix{Float64}            %% γ_Meu
        world : Matrix{Float64}         %% γ_Mw
        agg : Matrix{Float64}           %% γ_M
    }

    class Parameters {
        constants : Constants
        consumption : ParamsStruct      %% α
        export_eu : ParamsStruct        %% β1
        export_world : ParamsStruct     %% β2
        production : ParamsProduction   %% γ
        investment : ParamsStruct       %% ρ
        log : Bool
    }

    class ParameterSubset {
        constants : Constants
        production : ParamsProduction
    }

    %% ========= Relationships =========
    CleanData --> HouseholdData
    CleanData --> IndustryData
    CleanData --> Constants

    IndustryData --> AssetsLiabilities
    IndustryData --> RegionalData

    RegionalData --> InputMatrices
    RegionalData --> Totals

    Totals --> ForeignRegionalValues

    Constants --> ExchangeRates
    Constants --> ForeignRegionalValues
    Constants --> Elasticities

    Elasticities --> Elasticity

    Parameters --> Constants
    Parameters --> ParamsStruct : consumption
    Parameters --> ParamsStruct : export_eu
    Parameters --> ParamsStruct : export_world
    Parameters --> ParamsProduction : production
    Parameters --> ParamsStruct : investment

    ParameterSubset --> Constants
    ParameterSubset --> ParamsProduction

```
