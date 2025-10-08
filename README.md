# Supergrassi

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/UCL/Supergrassi.jl/actions/workflows/UnitTests.yml/badge.svg?branch=main)](https://github.com/UCL/Supergrassi.jl/actions/workflows/UnitTests.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/UCL/Supergrassi.jl/graph/badge.svg?token=LU852FO3FP)](https://codecov.io/gh/UCL/Supergrassi.jl)

A comprehensive Julia package for multi-sector dynamic macroeconomic modeling with debt and default, featuring constrained optimization of prices and other parameters.

## Overview

Supergrassi.jl provides a complete framework for macroeconomic analysis, including:

- **Multi-regional economic modeling** (UK, EU, World regions)
- **Industry-level analysis** with input-output matrices
- **Household and firm behavior modeling**
- **Capital market equilibrium**
- **Constrained optimization** using Ipopt solver
- **Automatic differentiation** with Enzyme.jl
- **Batch estimation** with error handling and logging

## Installation

```bash
git clone https://github.com/UCL/Supergrassi.jl.git
cd Supergrassi.jl
julia --project=.
```

## Quick Start

### Basic Data Pipeline

```julia
using Supergrassi

# Load configuration and data
settings_path = create_filepath("config/settings.yml")
settings = read_settings(settings_path)
filepaths = check_file_availability(settings)
data = read_data(filepaths, settings)

# Clean and process data
clean_data = clean_data(data, settings)
postprocess_clean_data!(clean_data)
```

### Run Economic Estimation

All the data processing steps are encapsulated in the `estimate` function. For batch processing, use `batch_estimation`.

```julia
# Single estimation
status, prob = estimate()

# Batch estimation with logging
batch_estimation(
    batch_size=100,
    log_errors=true,
    log_results=true,
    log_results_filepath="results.csv"
)
```

## Configuration

All configurations are managed through a YAML file (`config/settings.yml`).

The data files should be organized in the `input/` directory as specified in the configuration file.

## Data Structure

The main data structures used in Supergrassi.jl are illustrated below:

```mermaid
classDiagram
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
        total_use : DataFrame
        consumption : DataFrame
        delta_v : DataFrame
        export_eu : DataFrame
        export_world : DataFrame
        investment : DataFrame
        input_matrices : InputMatrices
        totals : Totals
    }

    class InputMatrices {
        uk : Matrix(Float64)
        eu : Matrix(Float64)
        world : Matrix(Float64)
        agg : Matrix(Float64)
    }

    class Totals {
        expenditure : Float64
        investments : Float64
        imports : ForeignRegionalValues
    }

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
        substitution : Float64
        armington : Float64
        substitution_uk_other : Float64?
        skill_substitution : Float64?
    }

    class Elasticities {
        production : Elasticity
        export_world : Elasticity
        export_eu : Elasticity
        consumption : Elasticity
        investment : Elasticity
    }

    class ParamsStruct {
        uk : Vector(Float64)
        eu : Vector(Float64)
        world : Vector(Float64)
        agg : Vector(Float64)
        tilde : Vector(Float64)?
    }

    class ParamsProduction {
        human : Vector(Float64)
        capital : Vector(Float64)
        low_skill : Vector(Float64)
        high_skill : Vector(Float64)
        shock_mean : Vector(Float64)
        shock_stddev : Vector(Float64)
        uk : Matrix(Float64)
        eu : Matrix(Float64)
        world : Matrix(Float64)
        agg : Matrix(Float64)
    }

    class Parameters {
        constants : Constants
        consumption : ParamsStruct
        export_eu : ParamsStruct
        export_world : ParamsStruct
        production : ParamsProduction
        investment : ParamsStruct
        log : Bool
    }

    class ParameterSubset {
        constants : Constants
        production : ParamsProduction
    }

    %% === Arrow labels ===
    CleanData --> HouseholdData : household
    CleanData --> IndustryData : industry
    CleanData --> Constants : constants

    IndustryData --> AssetsLiabilities : assets_liabilities
    IndustryData --> RegionalData : regional

    RegionalData --> InputMatrices : input_matrices
    RegionalData --> Totals : totals

    Totals --> ForeignRegionalValues : imports

    Constants --> ExchangeRates : exchange_rates
    Constants --> ForeignRegionalValues : total_imports_from_uk
    Constants --> ForeignRegionalValues : total_imports_from_all_sources
    Constants --> ForeignRegionalValues : import_tariffs
    Constants --> ForeignRegionalValues : export_costs
    Constants --> Elasticities : elasticities

    Elasticities --> Elasticity : production
    Elasticities --> Elasticity : export_world
    Elasticities --> Elasticity : export_eu
    Elasticities --> Elasticity : consumption
    Elasticities --> Elasticity : investment

    Parameters --> Constants : constants
    Parameters --> ParamsStruct : consumption
    Parameters --> ParamsStruct : export_eu
    Parameters --> ParamsStruct : export_world
    Parameters --> ParamsProduction : production
    Parameters --> ParamsStruct : investment

    ParameterSubset --> Constants : constants
    ParameterSubset --> ParamsProduction : production

```
