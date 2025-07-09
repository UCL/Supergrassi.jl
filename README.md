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
    class CleanData {
        **household** : HouseholdData
        **industry** : IndustryData
        **constants** : Constants
    }

    class HouseholdData {
        **income** : DataFrame
        **income_share** : DataFrame
        **payments** : DataFrame
        **hours** : DataFrame
        **wages** : DataFrame
    }

    class IndustryData {
        **depreciation** : DataFrame
        **tax** : DataFrame
        **capital** : DataFrame
        **surplus** : DataFrame
        **shock_stdev** : DataFrame
        **assets_liabilities** : AssetsLiabilities
        **regional** : RegionalData
    }

    class RegionalData {
        **total_use** : DataFrame
        **consumption** : DataFrame
        **delta_v** : DataFrame
        **export_eu** : DataFrame
        **export_world** : DataFrame
        **investment** : DataFrame
        **input_matrices** : InputMatrices
        **totals** : Totals
    }

    class InputMatrices {
        **uk** : DataFrame
        **eu** : DataFrame
        **world** : DataFrame
        **imports** : DataFrame
        **agg** : DataFrame
    }

    class Totals {
        **savings** : Float64
        **investments** : Float64
        **imports** : TotalImports
    }

    class TotalImports {
        **eu** : Float64
        **world** : Float64
    }

    class AssetsLiabilities {
        **current_year** : DataFrame
        **next_year** : DataFrame
    }

    class Constants {
        **data_year** : Int64
        **exchange_rates** : ExchangeRates
        **interest_rate** : Float64
        **total_imports_from_uk** : TotalImports
        **total_imports_from_all_sources** : TotalImports
    }

    class ExchangeRates {
        **usd** : Float64
        **eur** : Float64
    }

    CleanData --> HouseholdData
    CleanData --> IndustryData
    CleanData --> Constants

    IndustryData --> AssetsLiabilities
    IndustryData --> RegionalData

    RegionalData --> InputMatrices
    RegionalData --> Totals

    Totals --> TotalImports
    Constants --> ExchangeRates
    Constants --> TotalImports
```
