using DataFrames
using Printf
using StatsBase
using Dates


function create_map_105_to_64(merge_codes::DataFrame)

    initial_industry_names = merge_codes[!, :sic105]
    final_industry_names = merge_codes[!, :sic64]


    println("Initial industry names: ", initial_industry_names)
    println("Final industry names: ", final_industry_names)

    map_105_to_64 = Dict{String, String}()

    for i in eachindex(initial_industry_names)
        initial_name = initial_industry_names[i]
        final_name = final_industry_names[i]
        map_105_to_64[initial_name] = "SIC_64_" * @sprintf("%03i", final_name)
    end

    return map_105_to_64

end

function create_map_64_to_16(merge_codes::DataFrame)

    initial_industry_names = merge_codes[2:end, :x1]
    final_industry_names = merge_codes[2:end, :x7]

    println("Initial industry names: ", initial_industry_names)
    println("Final industry names: ", final_industry_names)

    map_64_to_16 = Dict{String, String}()

    for i in eachindex(initial_industry_names)
        initial_name = "SIC_64_" * @sprintf("%03i", parse(Int64, initial_industry_names[i]))
        final_name = final_industry_names[i]
        map_64_to_16[initial_name] = final_name
    end

    return map_64_to_16

end

function reduce_columns_by_group_sum(df::DataFrame, mapping::Dict{String, String})

    grouped = group_columns_by_new_name(mapping)

    # Sum columns per group
    new_cols = Dict{Symbol, Vector{eltype(df[!, 1])}}()
    for (new_name, old_syms) in grouped
        new_cols[Symbol(new_name)] = sum(eachcol(df[!, old_syms]))
    end

    return DataFrame(new_cols)
end

function reduce_columns_by_group_weighted_mean(df::DataFrame, mapping::Dict{String, String}, weights::DataFrame)

    grouped = group_columns_by_new_name(mapping)

    dd = df .* weights

    # Sum columns per group
    new_cols = Dict{Symbol, Vector{eltype(df[!, 1])}}()

    for (new_name, old_syms) in grouped
        new_cols[Symbol(new_name)] = sum(eachcol(dd[!, old_syms])) ./ sum(eachcol(weights[!, old_syms]))
    end

    return DataFrame(new_cols)

end

function group_columns_by_new_name(mapping::Dict{String, String})

    grouped = Dict{String, Vector{Symbol}}()
    for (old, new) in mapping
        push!(get!(grouped, new, Vector{Symbol}()), Symbol(old))
    end

    return grouped

end

function group_dataframes(dfs, col_names, industry_names, industries_on_cols = true, reduction_fun = nothing, reduction_map = nothing)

    dd = DataFrame([industry_names], ["industry"])
    for (df, col) in zip(dfs, col_names)

        # Do data aggregation if reduction function and map are passed in
        if (!isnothing(reduction_fun) & !isnothing(reduction_map))
            df = reduction_fun(df, reduction_map)
        end

        dd[!,col] = Float64.(permutedims(df).x1)

    end

    if industries_on_cols
        return permutedims(dd, :industry, "key")
    else
        return dd
    end

end

function select_year(data::DataFrame, year::Int64)
    rr = data[data.year .== year, 4:end]

    for col in names(rr)
        if eltype(rr[!, col]) <: AbstractString
            rr[!, col] = parse.(Float64, rr[!, col])
        end
    end

    return rr
end

function combine_dataframe_row_wise(data::DataFrame, func::Function)
    return combine(data, names(data, Real) .=> func, renamecols=false)
end

function clean_rows(data::DataFrame, row_names::String, col_names::Array{String, 1}, mapping::Dict{String, String})
    rr = data[data.x2 .== row_names, 3:end]

    rename!(rr, col_names)

    rr = reduce_columns_by_group_sum(rr, mapping)

    return rr
end

function clean_vector(data::Vector{<:Number}, industry_names::Array{String, 1}, mapping::Dict{String, String})
    # Create a DataFrame from the vector
    df = DataFrame([data], :auto)
    df = DataFrame(permutedims(df), industry_names)
    rr = reduce_columns_by_group_sum(df, mapping)

    return rr
end

function clean_matrix(data::DataFrame, industry_names::Array{String, 1}, mapping::Dict{String, String})


    rr = reduce_columns_by_group_sum(data, mapping)
    final_names = names(rr)

    rr = DataFrame(permutedims(rr), industry_names)

    rr = reduce_columns_by_group_sum(rr, mapping)
    rr = DataFrame(permutedims(rr), final_names)

    return rr
end

function clean_assets_liabilities(assets::DataFrame, year::Int64, n_samples::Int64 = nrow(assets))

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

    return limited

end

struct CleanData

    average_interest::Float64

    income::DataFrame # data.income_*
    income_share::DataFrame # data.*_income_share
    payments::DataFrame # data.hValue*
    depreciation::DataFrame # data.depreciation
    tax::DataFrame # data.taxValue*
    mean_capital::DataFrame # data.k*

    total_use::DataFrame # data.yValue
    consumption::DataFrame # data.fValue
    delta_v::DataFrame # data.deltaVValue
    export_eu::DataFrame # data.x1Value
    export_world::DataFrame # data.x2Value
    operating_surplus::DataFrame # data.kValue
    capital_formation::DataFrame  # data.IValueAlso called "payments to capital" in the code

    import_export_matrix::DataFrame # data.mValue
    import_export_matrix_eu::DataFrame # data.mValueEU
    import_export_matrix_world::DataFrame # data.mValueW
    import_export_matrix_imports::DataFrame # data.mValueIMP

    asset_liability_current_year::DataFrame
    asset_liability_next_year::DataFrame

end

function clean_data(data::Data, year::Int64)

    ################################################################

    industry_names = data.input_output.industry_names

    mapping_105_to_64 = create_map_105_to_64(data.merge_codes_105)
    mapping_64_to_16 = create_map_64_to_16(data.merge_codes_64)

    low_income = select_year(data.household.income.low, year)
    low_income = combine_dataframe_row_wise(low_income, sum)

    high_income = select_year(data.household.income.high, year)
    high_income = combine_dataframe_row_wise(high_income, sum)

    capital_current_year = select_year(data.industry.capital, year)
    mean_capital_current_year = combine_dataframe_row_wise(capital_current_year, mean)

    capital_next_year = select_year(data.industry.capital, year + 1)
    mean_capital_next_year = combine_dataframe_row_wise(capital_next_year, mean)

    ######################################
    tax_products = clean_rows(data.others, "Taxes less subsidies on products", industry_names, mapping_105_to_64)
    tax_production = clean_rows(data.others, "Taxes less subsidies on production", industry_names, mapping_105_to_64)
    compensation_employees = clean_rows(data.others, "Compensation of employees", industry_names, mapping_105_to_64)
    gross_operating_surplus_and_mixed_income = clean_rows(data.others, "Gross operating surplus and mixed income", industry_names, mapping_105_to_64)
    #######################################

    final_consumption = clean_vector(data.input_output.final_consumption, industry_names, mapping_105_to_64)
    gross_fixed_capital_formation = clean_vector(data.input_output.gross_fixed_capital_formation, industry_names, mapping_105_to_64)
    delta_v_value_uk = clean_vector(data.input_output.delta_v_value_uk, industry_names, mapping_105_to_64)

    export_eu = clean_vector(data.input_output.exports_eu_to_uk, industry_names, mapping_105_to_64)
    export_world = clean_vector(data.input_output.export_world_to_uk, industry_names, mapping_105_to_64)

    total_use_uk = clean_vector(data.input_output.total_use, industry_names, mapping_105_to_64)
    services_export = clean_vector(data.input_output.services_export, industry_names, mapping_105_to_64)

    export_ratio_eu_vs_eu_and_world = DataFrame(export_eu ./ (export_eu .+ export_world))
    export_ratio_eu_vs_eu_and_world .= ifelse.(isnan.(export_ratio_eu_vs_eu_and_world), 0.5, export_ratio_eu_vs_eu_and_world)

    export_eu = export_eu .+ export_ratio_eu_vs_eu_and_world .* services_export
    export_world = export_world .+ (1 .- export_ratio_eu_vs_eu_and_world) .* services_export

    ##########################################

    nms = names(tax_products)

    low_income = DataFrame(low_income, nms)
    high_income = DataFrame(high_income, nms)
    mean_capital_current_year = DataFrame(mean_capital_current_year, nms)
    mean_capital_next_year = DataFrame(mean_capital_current_year, nms)

    high_income_share = high_income ./ (high_income .+ low_income)
    low_income_share = low_income ./ (high_income .+ low_income)

    import_export_matrix = clean_matrix(data.input_output.input_output_matrix, industry_names, mapping_105_to_64)

#    depreciation = DataFrame([nms data.depreciation[!, string(year)]], ["industry", "val"])
    depreciation = DataFrame(data.depreciation)[!, string(year)]
    depreciation = DataFrame(permutedims(depreciation), nms)

    #################################################

    payments_to_low_skilled = (1 .- high_income_share) .* compensation_employees
    payments_to_high_skilled = high_income_share .* compensation_employees

    #################################################

    imports_import_export_matrix = clean_matrix(data.imports.input_output_matrix, industry_names, mapping_105_to_64)
    imports_final_consumption = clean_vector(data.imports.final_consumption, industry_names, mapping_105_to_64)
    imports_gross_fixed_capital_formation = clean_vector(data.imports.gross_fixed_capital_formation, industry_names, mapping_105_to_64)
    imports_delta_v_value_uk = clean_vector(data.imports.delta_v_value_uk, industry_names, mapping_105_to_64)
    imports_export_eu = clean_vector(data.imports.exports_eu_to_uk, industry_names, mapping_105_to_64)
    imports_export_world = clean_vector(data.imports.export_world_to_uk, industry_names, mapping_105_to_64)
    imports_total_use = clean_vector(data.imports.total_use, industry_names, mapping_105_to_64)
    imports_services_export = clean_vector(data.imports.services_export, industry_names, mapping_105_to_64)
    imports_export_ratio_eu_vs_eu_and_world = DataFrame(imports_export_eu ./ (imports_export_eu .+ imports_export_world))
    imports_export_ratio_eu_vs_eu_and_world .= ifelse.(isnan.(imports_export_ratio_eu_vs_eu_and_world), 0.5, imports_export_ratio_eu_vs_eu_and_world)
    imports_export_eu = imports_export_eu .+ imports_export_ratio_eu_vs_eu_and_world .* imports_services_export
    imports_export_world = imports_export_world .+ (1 .- imports_export_ratio_eu_vs_eu_and_world) .* imports_services_export

    ######################################################

    split_factor = 0.5

    eu_import_export_matrix = imports_import_export_matrix .* split_factor
    eu_final_consumption = imports_final_consumption .* split_factor
    eu_gross_fixed_capital_formation = imports_gross_fixed_capital_formation .* split_factor
    eu_delta_v_value_uk = imports_delta_v_value_uk .* split_factor
    eu_export_eu = imports_export_eu .* split_factor
    eu_export_world = imports_export_world .* split_factor
    eu_total_use = imports_total_use .* split_factor
    eu_services_export = imports_services_export .* split_factor

    world_import_export_matrix = imports_import_export_matrix .* (1 - split_factor)
    world_final_consumption = imports_final_consumption .* (1 - split_factor)
    world_gross_fixed_capital_formation = imports_gross_fixed_capital_formation .* (1 - split_factor)
    world_delta_v_value_uk = imports_delta_v_value_uk .* (1 - split_factor)
    world_export_eu = imports_export_eu .* (1 - split_factor)
    world_export_world = imports_export_world .* (1 - split_factor)
    world_total_use = imports_total_use .* (1 - split_factor)
    world_services_export = imports_services_export .* (1 - split_factor)

    #############################################################

    interest_rates = data.risk_free_rate[Dates.year.(data.risk_free_rate.date) .== year, 2:end]

    for interest_rate in names(interest_rates)
        interest_rates[!, interest_rate] = parse.(Float64, interest_rates[!, interest_rate])
    end

    R = 1 + geomean(interest_rates[!, 1] / 100)

    #############################################################

    asset_liability_current_year = clean_assets_liabilities(data.assets, year, 1000)
    asset_liability_next_year = clean_assets_liabilities(data.assets, year + 1)

    #############################################################

    # Join data it is split either by "low"/"high", or "uk"/"eu"/"world"/"imports"
    # Aggregate from 64 to 16 industries.
    # I did these in the same function to reduce the number of function calls

    sic16 = unique(data.merge_codes_64.x7[2:end]) # TODO: Check the order stays correct
    industries_in_cols = false

    income = group_dataframes([low_income, high_income],
                              ["low", "high"], sic16,
                              industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    income_share = group_dataframes([low_income_share, high_income_share],
                                    ["low", "high"], sic16,
                                    industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    payments = group_dataframes([payments_to_low_skilled, payments_to_high_skilled, compensation_employees],
                                ["low", "high", "agg"], sic16,
                                industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    mean_capital = group_dataframes([mean_capital_current_year, mean_capital_next_year],
                                    ["current", "next"], sic16,
                                    industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    tax = group_dataframes([tax_products, tax_production],
                           ["products", "production"], sic16,
                           industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    consumption = group_dataframes([final_consumption, eu_final_consumption,
                                    world_final_consumption, imports_final_consumption],
                                   ["uk", "eu", "world", "imports"], sic16,
                                   industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    delta_v = group_dataframes([delta_v_value_uk, eu_delta_v_value_uk, world_delta_v_value_uk, imports_delta_v_value_uk],
                               ["uk", "eu", "world", "imports"], sic16,
                               industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    export_to_eu = group_dataframes([export_eu, eu_export_eu, world_export_eu, imports_export_eu],
                                    ["uk", "eu", "world", "imports"], sic16,
                                    industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    export_to_world = group_dataframes([export_world, eu_export_world, world_export_world, imports_export_world],
                                       ["uk", "eu", "world", "imports"], sic16,
                                       industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    total_use = group_dataframes([total_use_uk, eu_total_use, world_total_use, imports_total_use],
                                 ["uk", "eu", "world", "imports"], sic16,
                                 industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    capital_formation = group_dataframes([gross_fixed_capital_formation, eu_gross_fixed_capital_formation,
                                          world_gross_fixed_capital_formation, imports_gross_fixed_capital_formation],
                                         ["uk", "eu", "world", "imports"], sic16,
                                         industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    services_export = group_dataframes([services_export, eu_services_export, world_services_export, imports_services_export],
                                       ["uk", "eu", "world", "imports"], sic16,
                                       industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    surplus = group_dataframes([gross_operating_surplus_and_mixed_income],
                               ["val"], sic16,
                               industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)
    depreciation = group_dataframes([depreciation],
                                    ["val"], sic16,
                                    industries_in_cols, reduce_columns_by_group_sum, mapping_64_to_16)

    return CleanData(
        R,
        income,
        income_share,
        payments,
        depreciation,
        tax,
        mean_capital,
        total_use,
        consumption,
        delta_v,
        export_to_eu,
    export_to_world,
    surplus,
    capital_formation,
    import_export_matrix,
    eu_import_export_matrix,
    world_import_export_matrix,
    imports_import_export_matrix,
    asset_liability_current_year,
    asset_liability_next_year)

end
