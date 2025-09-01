using StatsBase

"""
    function clean_rows(data::DataFrame, row_names::String, col_names::Array{String, 1}, mapping::Dict{String, String})

Clean selected subset of rows in a DataFrame based on row names and rename columns.

# Arguments
- `data::DataFrame`: The DataFrame containing the data to be cleaned.
- `row_names::String`: The name of the row to select.
- `col_names::Array{String, 1}`: An array of new column names to rename the selected columns.
- `mapping::Dict{String, String}`: A dictionary mapping old column names to new column names for reduction.

# Returns
- `DataFrame`: A new DataFrame with selected rows, renamed columns, and reduced columns based on the mapping.
"""
function clean_rows(data::DataFrame, row_names::String, col_names::Array{String, 1}, mapping::Dict{String, String})
    rr = data[data.x2 .== row_names, 3:end]

    rename!(rr, col_names)

    rr = reduce_columns_by_group_sum(rr, mapping)

    return rr
end

"""
    function clean_vector(data::Vector{<:Number}, industry_names::Array{String, 1}, mapping::Dict{String, String})

Clean a vector of numbers by creating a DataFrame, renaming columns, and reducing columns based on a mapping.

# Arguments
- `data::Vector{<:Number}`: A vector of numbers to be cleaned.
- `industry_names::Array{String, 1}`: An array of industry names to be used as column names.
- `mapping::Dict{String, String}`: A dictionary mapping old column names to new column names for reduction.

# Returns
- `DataFrame`: A new DataFrame with the vector data, renamed columns, and reduced columns based on the mapping.
"""
function clean_vector(data::Vector{<:Number}, industry_names::Array{String, 1}, mapping::Dict{String, String})
    # Create a DataFrame from the vector
    df = DataFrame([data], :auto)
    df = DataFrame(permutedims(df), industry_names)
    rr = reduce_columns_by_group_sum(df, mapping)

    return rr
end

"""
    function clean_matrix(data::DataFrame, industry_names::Array{String, 1}, mapping::Dict{String, String})

Clean a matrix of data by reducing based on a mapping and renaming them with industry names.

# Arguments
- `data::DataFrame`: A DataFrame containing the matrix data to be cleaned.
- `industry_names::Array{String, 1}`: An array of industry names to be used as column names.
- `mapping::Dict{String, String}`: A dictionary mapping old column names to new column names for reduction.

# Returns
- `DataFrame`: A new DataFrame with the matrix data, renamed columns, and reduced columns based on the mapping.
"""
function clean_matrix(data::DataFrame, industry_names::Array{String, 1}, mapping::Dict{String, String})


    rr = reduce_columns_by_group_sum(data, mapping)

    rr = DataFrame(permutedims(rr), industry_names)

    rr = reduce_columns_by_group_sum(rr, mapping)

    return rr
end

"""
    funvtion clean_assets_liabilities(assets::DataFrame, year::Int64, map_to_16::Dict{String, String}, n_samples::Int64 = nrow(assets))

Wrapper function to clean assets and liabilities data.

# Arguments
- `assets::DataFrame`: The DataFrame containing assets and liabilities data.
- `year::Int64`: The year for which the data is to be cleaned.
- `map_to_16::Dict{String, String}`: A dictionary mapping SIC 64 industry names to SIC 16 industry names.
- `n_samples::Int64`: The number of samples to limit per SIC 64 industry (default is the number of rows in `assets`.).

# Returns
- `DataFrame`: A cleaned DataFrame containing assets, liabilities, and their ratio, limited to `n_samples` per SIC 64 industry.
"""
function clean_assets_liabilities(assets::DataFrame, year::Int64, map_to_16::Dict{String, String}, n_samples::Int64 = nrow(assets))

    # Step 1: Extract the relevant columns as strings
    year_str = string(year)
    assets_col = Symbol("total_assets_" * year_str)
    liabilities_col = Symbol("total_liabilities_" * year_str)

    # Step 2: Create a DataFrame from those columns
    df_raw = DataFrame(
        SIC64 = assets[!, :sic64],
        Assets = assets[!, assets_col],
        Liabilities = assets[!, liabilities_col],
    )

    # Step 3: Filter out rows with "NA"
    df_clean = filter(row -> row.Assets != "NA" && row.Liabilities != "NA" && row.SIC64 != "NA", df_raw)

    # Step 4: Convert strings to Float64
    df_clean.Assets = parse.(Float64, df_clean.Assets)
    df_clean.Liabilities = parse.(Float64, df_clean.Liabilities)

    # Step 5: Compute Liabilities/Assets, handling Assets == 0
    ratio = [a == 0.0 ? 0.0 : abs(l) / a for (l, a) in zip(df_clean.Liabilities, df_clean.Assets)]

    # Step 6: Build final DataFrame and filter on ratio <= 1
    df_final = DataFrame(
        SIC64 = df_clean.SIC64,
        Assets = df_clean.Assets,
        Ratio = ratio
    )

    # Step 7: Drop rows where Ratio > 1
    df_final = filter(row -> row.Ratio ≤ 1.0, df_final)

    # display(df_final)

    # Step 8: Limit to n_samples rows per SIC64
    grouped = groupby(df_final, :SIC64)
    limited = combine(grouped) do sdf
        n = nrow(sdf)
        n > n_samples ? sdf[rand(1:n, n_samples), :] : sdf
    end

    limited.SIC16 = missings(String, nrow(limited))

    for i in 1:nrow(limited)

        limited[i,"SIC16"] = map_to_16["SIC_64_" * @sprintf("%03i", limited.SIC64[i])]

    end

    return limited

end


"""
    function clean_2d_values(data::Data, split_factor::Float64)

Function to process the 2D values that are split between uk, eu, world and stored in matrices.

# Arguments
- `data::Data`: The Data struct containing input-output matrices and imports.
- `split_factor::Float64`: The factor by which to split the imports between EU and World.

# Returns
- `InputMatrices`: A struct containing the cleaned input-output matrices for UK, EU, World, and aggregate values.
"""
function clean_2d_values(data::Data, split_factor::Float64)

    industry_names = data.input_output.industry_names
    mapping_105_to_64 = create_map_105_to_64(data.merge_codes_105)
    mapping_64_to_16 = create_map_64_to_16(data.merge_codes_64)

    import_export = clean_matrix(data.input_output.input_output_matrix, industry_names, mapping_105_to_64)
    import_export = Matrix(clean_matrix(import_export, names(import_export), mapping_64_to_16))

    imports_import_export = clean_matrix(data.imports.input_output_matrix, industry_names, mapping_105_to_64)
    imports_import_export = Matrix(clean_matrix(imports_import_export, names(imports_import_export), mapping_64_to_16))

    eu_import_export = imports_import_export .* split_factor
    world_import_export = imports_import_export .* (1 - split_factor)

    return InputMatrices(import_export, eu_import_export, world_import_export, import_export .+ eu_import_export .+ world_import_export)


end

"""
    function clean_1d_values(val::Vector{<:Number}, val_imp::Vector{<:Number}, map_64::Dict{String, String}, map_16::Dict{String, String}, names_105::Vector, names_16::Vector, split_factor::Float64, industries_in_cols::Bool)

Function to process the values that are split between uk, eu, world and stored in vectors.

# Arguments
- `val::Vector{<:Number}`: The vector containing the values for UK.
- `val_imp::Vector{<:Number}`: The vector containing the values for imports.
- `map_64::Dict{String, String}`: A dictionary mapping SIC 64 industry names to SIC 16 industry names.
- `map_16::Dict{String, String}`: A dictionary mapping SIC 16 industry names to their final names.
- `names_105::Vector`: A vector of industry names for SIC 105.
- `names_16::Vector`: A vector of industry names for SIC 16.
- `split_factor::Float64`: The factor by which to split the imports between EU and World.
- `industries_in_cols::Bool`: If true, industries will be on columns; otherwise, they will be on rows.

# Returns
- `DataFrame`: A DataFrame containing the cleaned values for UK, EU, World, and aggregate values.
"""
function clean_1d_values(val::Vector{<:Number}, val_imp::Vector{<:Number}, map_64::Dict{String, String}, map_16::Dict{String, String}, names_105::Vector, names_16::Vector, split_factor::Float64, industries_in_cols::Bool)

    val_uk = clean_vector(val, names_105, map_64)
    val_imports = clean_vector(val_imp, names_105, map_64)
    val_eu = val_imports .* split_factor
    val_world = val_imports .* (1 - split_factor)

    df =  group_dataframes([val_uk, val_eu, val_world],
                           ["uk", "eu", "world"], names_16,
                           industries_in_cols, reduce_columns_by_group_sum, map_16)

    add_aggregate!(df)
    return df

end

"""
    function correct_exports_with_services!(export_to_eu::DataFrame, export_to_world::DataFrame, services_export::DataFrame)

Function to correct exports by accounting for NaN values and scaling appropriately.

# Arguments
- `export_to_eu::DataFrame`: DataFrame containing exports to the EU.
- `export_to_world::DataFrame`: DataFrame containing exports to the World.
- `services_export::DataFrame`: DataFrame containing service exports.
"""
function correct_exports_with_services!(export_to_eu::DataFrame, export_to_world::DataFrame, services_export::DataFrame)

    export_ratio_eu_vs_eu_and_world = export_to_eu ./ (export_to_eu .+ export_to_world)
    export_ratio_eu_vs_eu_and_world .= ifelse.(isnan.(export_ratio_eu_vs_eu_and_world), 0.5, export_ratio_eu_vs_eu_and_world)

    export_to_eu .= export_to_eu .+ export_ratio_eu_vs_eu_and_world .* services_export
    export_to_world .= export_to_world .+ (1 .- export_ratio_eu_vs_eu_and_world) .* services_export

end

"""
    function clean_exports(input_output::InputOutput, imports::InputOutput, split::Float64, names_16::Vector, industries_in_cols::Bool, map_64::Dict{String, String}, map_16::Dict{String, String})

Wrapper function to clean exports data.

Exports need a special treatment because they are a sum of export and services_export data frames.
Service export is scaled by the sum of eu and world exports. Further refactoring definitely possible here.

# Arguments
- `input_output::InputOutput`: The InputOutput struct containing export data.
- `imports::InputOutput`: The InputOutput struct containing import data.
- `split::Float64`: The factor by which to split the imports between EU and World.
- `names_16::Vector`: A vector of industry names for SIC 16.
- `industries_in_cols::Bool`: If true, industries will be on columns; otherwise, they will be on rows.
- `map_64::Dict{String, String}`: A dictionary mapping SIC 64 industry names to SIC 16 industry names.
- `map_16::Dict{String, String}`: A dictionary mapping SIC 16 industry names to their final names.

# Returns
- `DataFrame, DataFrame`: Two DataFrames containing exports to the EU and World, respectively.
"""
function clean_exports(input_output::InputOutput, imports::InputOutput, split::Float64, names_16::Vector, industries_in_cols::Bool, map_64::Dict{String, String}, map_16::Dict{String, String})

    eu = clean_vector(input_output.exports_eu_to_uk, input_output.industry_names, map_64)
    world = clean_vector(input_output.export_world_to_uk, input_output.industry_names, map_64)
    services = clean_vector(input_output.services_export, input_output.industry_names, map_64)

    correct_exports_with_services!(eu, world, services)

    eu_imp = clean_vector(imports.exports_eu_to_uk, imports.industry_names, map_64)
    world_imp = clean_vector(imports.export_world_to_uk, imports.industry_names, map_64)
    services_imp = clean_vector(imports.services_export, imports.industry_names, map_64)

    correct_exports_with_services!(eu_imp, world_imp, services_imp)

    exports_to_eu = group_dataframes([eu, eu_imp .* split, eu_imp .* (1 - split)],
                                     ["uk", "eu", "world"], names_16,
                                     industries_in_cols, reduce_columns_by_group_sum, map_16)

    exports_to_world = group_dataframes([world, world_imp .* split, world_imp .* (1 - split)],
                                        ["uk", "eu", "world"], names_16,
                                        industries_in_cols, reduce_columns_by_group_sum, map_16)

    add_aggregate!(exports_to_eu)
    add_aggregate!(exports_to_world)

    return exports_to_eu, exports_to_world

end


"""
    function clean_household(data::Data, year::Int64, map_64::Dict{String, String}, map_16::Dict{String, String}, names_105::Vector, names_64::Vector, names_16::Vector, industries_in_cols::Bool)

Function to process the household incomes and their derived data.

# Arguments
- `data::Data`: The Data struct containing household income and hours data.
- `year::Int64`: The year for which the data is to be cleaned.
- `map_64::Dict{String, String}`: A dictionary mapping SIC 105 industry names to SIC 64 industry names.
- `map_16::Dict{String, String}`: A dictionary mapping SIC 64 industry names to SIC 16 industry names.
- `names_105::Vector`: A vector of industry names for SIC 105.
- `names_64::Vector`: A vector of industry names for SIC 64.
- `names_16::Vector`: A vector of industry names for SIC 16.
- `industries_in_cols::Bool`: If true, industries will be on columns; otherwise, they will be on rows.

# Returns
- `HouseholdData`: A struct containing cleaned household income, payments, hours, and wages data.
"""
function clean_household(data::Data, year::Int64, map_64::Dict{String, String}, map_16::Dict{String, String}, names_105::Vector, names_64::Vector, names_16::Vector, industries_in_cols::Bool)

    compensation_employees = clean_rows(data.others, "Compensation of employees", names_105, map_64)

    low_income = merge_quarterly_data(data.household.income.low, year, names_64, sum)
    high_income = merge_quarterly_data(data.household.income.high, year, names_64, sum)

    high_income_share = high_income ./ (high_income .+ low_income)

    payments_to_low_skilled = (1 .- high_income_share) .* compensation_employees
    payments_to_high_skilled = high_income_share .* compensation_employees

    low_hours = merge_quarterly_data(data.household.hours.low, year, names_64, sum)
    high_hours = merge_quarterly_data(data.household.hours.high, year, names_64, sum)

    income = group_dataframes([low_income, high_income], ["low", "high"], names_16,
                              industries_in_cols, reduce_columns_by_group_sum, map_16)
    payments = group_dataframes([payments_to_low_skilled, payments_to_high_skilled, compensation_employees],
                                ["low", "high", "agg"], names_16, industries_in_cols, reduce_columns_by_group_sum,
                                map_16)

    hours = group_dataframes([low_hours, high_hours], ["low", "high"], names_16,
                             industries_in_cols, reduce_columns_by_group_sum, map_16)

    wages = DataFrame([names_16, payments.low ./ hours.low, payments.high ./ hours.high], ["industries", "low", "high"])
    mapcols(col -> replace!(col, NaN=>0.0), wages)

    return HouseholdData(income, payments, hours, wages)

end


"""
    clean_sigma_bar(sigma_data::Vector, zero_list::Vector{Int64}, sic64::Vector{String})
Clean the sigma bar data by converting it to a DataFrame, parsing strings to Float64, and setting specified indices to zero.

# Arguments
- `sigma_data::Vector`: A vector containing sigma bar data.
- `zero_list::Vector{Int64}`: A vector of indices where the sigma bar values should be set to zero.
- `sic64::Vector{String}`: A vector of SIC 64 industry names to be used as column names.

# Returns
- `DataFrame`: A DataFrame containing the cleaned sigma bar data, with specified indices set to zero.
"""
function  clean_sigma_bar(sigma_data::Vector, zero_list::Vector{Int64}, sic64::Vector{String})

    sigma_bar = DataFrame(permutedims(sigma_data), sic64)

    parse_string_dataframe!(sigma_bar, Float64, 0.0)

    for idx in zero_list
        if idx ≥ 1 && idx ≤ size(sigma_bar, 2)
            sigma_bar[1, idx] = 0.0
        end
    end

    return sigma_bar

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

"""
    clean_data(data::Data, settings::Dict{String, Any})

Main function for data cleaning. Should take in a Data struct and return a CleanData struct.

# Arguments
- `data::Data`: The Data struct containing various economic data.
- `settings::Dict{String, Any}`: A dictionary containing constants and settings for the model.
"""
function clean_data(data::Data, settings::Dict{String, Any})

    year::Int64 = settings["constants"]["data_year"]

    industry_names = data.input_output.industry_names
    aggregated_names = unique(data.merge_codes_64.x7[2:end])

    mapping_105_to_64 = create_map_105_to_64(data.merge_codes_105)
    mapping_64_to_16 = create_map_64_to_16(data.merge_codes_64)

    industries_in_cols = false
    split = settings["constants"]["imports_split"]

    # Clean dataframes from the data struct.
    # These contain data from 105 industries, aggregate to 64 industries

    tax_products = clean_rows(data.others, "Taxes less subsidies on products", industry_names, mapping_105_to_64)
    sic64 = names(tax_products)
    tax_production = clean_rows(data.others, "Taxes less subsidies on production", industry_names, mapping_105_to_64)

    total_use = clean_1d_values(data.input_output.total_use, data.imports.total_use, mapping_105_to_64, mapping_64_to_16, industry_names, aggregated_names, split, industries_in_cols)

    # ROW is not added to the aggregate value of total_use for some reason
    total_use.agg .-= (total_use.eu .+ total_use.world)

    consumption = clean_1d_values(data.input_output.final_consumption, data.imports.final_consumption, mapping_105_to_64, mapping_64_to_16, industry_names, aggregated_names, split, industries_in_cols)
    delta_v = clean_1d_values(data.input_output.delta_v_value_uk, data.imports.delta_v_value_uk, mapping_105_to_64, mapping_64_to_16, industry_names, aggregated_names, split, industries_in_cols)
    investment = clean_1d_values(data.input_output.gross_fixed_capital_formation, data.imports.gross_fixed_capital_formation, mapping_105_to_64, mapping_64_to_16, industry_names, aggregated_names, split, industries_in_cols)

    export_to_eu, export_to_world = clean_exports(data.input_output, data.imports, split, aggregated_names, industries_in_cols, mapping_105_to_64, mapping_64_to_16)

    input_matrices = clean_2d_values(data, split)

    # Store sums of some variables that are needed later.
    # These variables get re-scaled to sum to 1 in the postprocessing step and this information gets otherwise lost.
    total_vals = Totals(sum(consumption.agg) / mean(total_use.uk),
                        sum(investment.agg) / mean(total_use.uk),
                        ForeignRegionalValues(sum(export_to_eu.agg) / mean(total_use.uk),
                                     sum(export_to_world.agg) / mean(total_use.uk)))

    household = clean_household(data, year, mapping_105_to_64, mapping_64_to_16, industry_names, sic64, aggregated_names, industries_in_cols)


    ###################################################################################
    ###################################################################################

    gross_operating_surplus_and_mixed_income = clean_rows(data.others, "Gross operating surplus and mixed income", industry_names, mapping_105_to_64)


    regional = RegionalData(total_use,
        consumption,
        delta_v,
        export_to_eu,
        export_to_world,
        investment,
        input_matrices,
        total_vals
        )

    asset_liability_current_year = clean_assets_liabilities(data.assets, year, mapping_64_to_16, 1000)
    asset_liability_next_year = clean_assets_liabilities(data.assets, year + 1, mapping_64_to_16)
    assets_liabilities = AssetsLiabilities(asset_liability_current_year,
                                           asset_liability_next_year)
    # Process dataframes from the data struct.
    # These have quarterly data and are merged to yearly.
    # Already contain data from 64 industries

    mean_capital_current_year = merge_quarterly_data(data.industry.capital, year, sic64, mean)
    mean_capital_next_year = merge_quarterly_data(data.industry.capital, year + 1, sic64, mean)
    mean_capital = group_dataframes([mean_capital_current_year, mean_capital_next_year],
                                    ["current_year", "next_year"], aggregated_names, industries_in_cols,
                                    reduce_columns_by_group_sum, mapping_64_to_16)

    # Process data frames from the data struct.
    # These require other custom processing.

    depreciation = DataFrame(data.depreciation)[!, string(year)]
    depreciation = DataFrame(permutedims(depreciation), sic64)
    depreciation = group_dataframes([depreciation], ["val"], aggregated_names, industries_in_cols, reduce_columns_by_group_weighted_mean, mapping_64_to_16, weights = mean_capital_current_year)

    total_use_uk = clean_vector(data.input_output.total_use, industry_names, mapping_105_to_64)

    sigma_bar = clean_sigma_bar(data.model_results.sigma, settings["constants"]["zero_sigma_bar_industry_codes"], sic64)
    sigma_bar = group_dataframes([sigma_bar], ["val"], aggregated_names, industries_in_cols, reduce_columns_by_group_weighted_mean, mapping_64_to_16, weights = total_use_uk)

    # Join data it is split either by "low"/"high", or "uk"/"eu"/"world"/"imports"
    # Aggregate from 64 to 16 industries.
    # I did these in the same function to reduce the number of function calls

    tax = group_dataframes([tax_products, tax_production], ["products", "production"], aggregated_names,
                           industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    surplus = group_dataframes([gross_operating_surplus_and_mixed_income], ["val"], aggregated_names,
                               industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)

    # total_use_uk before aggregation is needed as weights for sigma_bar

    industry = IndustryData(depreciation, tax, mean_capital, surplus, sigma_bar, assets_liabilities, regional)

    constants = generate_constants(data, settings)

    return CleanData(
        household,
        industry,
        constants,
    )

end
