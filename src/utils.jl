using Dates

"""
    function select_year(data::DataFrame, year::Int64)

Select rows from a DataFrame based on a specific year.

# Arguments
- `data::DataFrame`: The DataFrame containing the data.
- `year::Int64`: The year to filter the DataFrame by.

# Returns
- `DataFrame`: A new DataFrame containing only the rows for the specified year, with columns converted to Float64 if they are strings.
"""
function select_year(data::DataFrame, year::Int64)
    rr = data[data.year .== year, 4:end]

    for col in names(rr)
        if eltype(rr[!, col]) <: AbstractString
            rr[!, col] = parse.(Float64, rr[!, col])
        end
    end

    return rr
end

"""
    function parse_string_dataframe!(df::DataFrame, T::Type, default_val=nothing)

Parse string columns in a DataFrame to a specified type, replacing missing values with a default value.

# Arguments
- `df::DataFrame`: The DataFrame containing the string columns to be parsed.
- `T::Type`: The type to which the string columns should be parsed (e.g., `Float64`).
- `default_val`: The value to replace missing values with (default is `nothing`).
"""
function parse_string_dataframe!(df::DataFrame, T::Type, default_val=nothing)

    for v in names(df)
        df[!,v] = tryparse.(T, df[!,v])
        df[!,v] = replace(df[!,v], nothing => default_val)
    end


end

"""
    function merge_quarterly_data(df::DataFrame, year::Int64, industry_names::Vector, fun::Function)

Merge quarterly data for a specific year by applying a function to each row.

# Arguments
- `df::DataFrame`: The DataFrame containing quarterly data.
- `year::Int64`: The year for which the data is to be merged.
- `industry_names::Vector`: A vector of industry names to be used as row names.
- `fun::Function`: A function to combine the quarterly data (e.g., `sum`).

# Returns
- `DataFrame`: A new DataFrame with the merged data for the specified year, with industry names as row names.
"""
function merge_quarterly_data(df::DataFrame, year::Int64, industry_names::Vector, fun::Function)

    dd = select_year(df, year)
    dd = combine_dataframe_row_wise(dd, fun)
    dd = DataFrame(dd, industry_names)
    return dd

end

"""
    function add_aggregate!(df::DataFrame)

Add an aggregate column to a DataFrame that sums the UK, EU, and World columns.

# Arguments
- `df::DataFrame`: The DataFrame to which the aggregate column will be added.
"""
function add_aggregate!(df::DataFrame)

    df[!,:agg] = df.uk + df.eu + df.world

end

"""
    generate_constants(data::Data, settings::Dict{String, Any})

Generate a Constants struct from the provided data and settings.

# Arguments
- `data::Data`: The Data struct containing various economic data.
- `settings::Dict{String, Any}`: A dictionary containing constants and settings for the model.

# Returns
- `Constants`: A Constants struct containing the year, exchange rates, interest rate, total imports, import tariffs, export costs, and elasticities.
"""
function generate_constants(data::Data, settings::Dict{String, Any})

    year::Int64 = settings["constants"]["data_year"]

    exchange_rates = ExchangeRates(settings["constants"]["exchange_rates"]["usd"],
                                   settings["constants"]["exchange_rates"]["eur"])

    total_imports_from_uk = ForeignRegionalValues(settings["constants"]["total_imports"]["from_uk"]["eu"],
                                         settings["constants"]["total_imports"]["from_uk"]["world"])
    total_imports_from_all_sources = ForeignRegionalValues(settings["constants"]["total_imports"]["from_all_sources"]["eu"],
                                                  settings["constants"]["total_imports"]["from_all_sources"]["world"])

    import_tariffs = ForeignRegionalValues(settings["constants"]["import_tariff"]["eu"], 
                                  settings["constants"]["import_tariff"]["world"])
    export_costs = ForeignRegionalValues(settings["constants"]["export_costs"]["eu"],
                                settings["constants"]["export_costs"]["world"])

    elasticities = settings["constants"]["elasticities"]

    loss_given_default = settings["constants"]["loss_given_default"]

    production_elasticity = Elasticity(elasticities["production"][1],
                                       elasticities["production"][3],
                                       nothing,
                                       elasticities["production"][2])

    world_export_demand_elasticity = Elasticity(elasticities["rest_of_world_export_demand"][1],
                                                elasticities["rest_of_world_export_demand"][2],
                                                elasticities["rest_of_world_export_demand"][3],
                                                nothing)

    eu_export_demand_elasticity = Elasticity(elasticities["eu_export_demand"][1],
                                             elasticities["eu_export_demand"][2],
                                             elasticities["eu_export_demand"][3],
                                             nothing)

    consumption_elasticity = Elasticity(elasticities["consumption"][1],
                                        elasticities["consumption"][2],
                                        nothing,
                                        nothing)

    investment_elasticity = Elasticity(elasticities["investment"][1],
                                       elasticities["investment"][2],
                                       nothing,
                                       nothing)

    elasticities_struct = Elasticities(production_elasticity,
                                       world_export_demand_elasticity,
                                       eu_export_demand_elasticity,
                                       consumption_elasticity,
                                       investment_elasticity)

    interest_rates = data.risk_free_rate[Dates.year.(data.risk_free_rate.date) .== year, 2:end]
    parse_string_dataframe!(interest_rates, Float64)
    interest_rate = 1 + geomean(interest_rates[!, 1] / 100)

    return Constants(year, exchange_rates, interest_rate, total_imports_from_uk, total_imports_from_all_sources, 
                     import_tariffs, export_costs, elasticities_struct, loss_given_default)

end

