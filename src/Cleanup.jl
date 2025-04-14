using DataFrames
using Printf
using Statistics


function create_map_105_to_64(data::Data)

    initial_industry_names = data.merge_codes_105[!, :sic105]
    final_industry_names = data.merge_codes_105[!, :sic64]


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

function create_map_64_to_16(data::Data)

    initial_industry_names = data.merge_codes_64[2:end, :x1]
    final_industry_names = data.merge_codes_64[2:end, :x7]

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
    # Group old columns by new name
    grouped = Dict{String, Vector{Symbol}}()
    for (old, new) in mapping
        push!(get!(grouped, new, Vector{Symbol}()), Symbol(old))
    end

    # Sum columns per group
    new_cols = Dict{Symbol, Vector{eltype(df[!, 1])}}()
    for (new_name, old_syms) in grouped
        new_cols[Symbol(new_name)] = sum(eachcol(df[!, old_syms]))
    end

    return DataFrame(new_cols)
end

function reduce_columns_by_group_weighted_mean(df::DataFrame, mapping::Dict{String, String}, weights::DataFrame)
    # Group old columns by new name
    grouped = Dict{String, Vector{Symbol}}()
    for (old, new) in mapping
        push!(get!(grouped, new, Vector{Symbol}()), Symbol(old))
    end

    dd = df .* weights

    # Sum columns per group
    new_cols = Dict{Symbol, Vector{eltype(df[!, 1])}}()
    for (new_name, old_syms) in grouped
        new_cols[Symbol(new_name)] = sum(eachcol(dd[!, old_syms])) / sum(eachcol(weights[!, old_syms]))
    end

    return DataFrame(new_cols)

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


struct CleanData

    low_income::DataFrame
    high_income::DataFrame

    low_income_share::DataFrame
    high_income_share::DataFrame

    mean_capital_current_year::DataFrame
    mean_capital_next_year::DataFrame

    tax_products::DataFrame
    tax_production::DataFrame

    compensation_employees::DataFrame
    gross_operating_surplus_and_mixed_income::DataFrame


    final_consumption::DataFrame
    gross_fixed_capital_formation::DataFrame

    delta_v_value_uk::DataFrame

    exports_eu_to_uk::DataFrame
    export_world_to_uk::DataFrame

    total_use::DataFrame

    services_export::DataFrame

    import_export_matrix::DataFrame

    function CleanData(data::Data, year::Int64)



        low_income = select_year(data.household.income.low, year)
        low_income = combine_dataframe_row_wise(low_income, sum)

        high_income = select_year(data.household.income.high, year)
        high_income = combine_dataframe_row_wise(high_income, sum)

        capital_current_year = select_year(data.industry.capital, year)
        mean_capital_current_year = combine_dataframe_row_wise(capital_current_year, mean)

        capital_next_year = select_year(data.industry.capital, year + 1)
        mean_capital_next_year = combine_dataframe_row_wise(capital_next_year, mean)

        ################################################################

        industry_names = data.input_output.industry_names

        mapping_105_to_64 = create_map_105_to_64(data)

        ######################################
        tax_products = clean_rows(data.others, "Taxes less subsidies on products", industry_names, mapping_105_to_64)
        tax_production = clean_rows(data.others, "Taxes less subsidies on production", industry_names, mapping_105_to_64)
        compensation_employees = clean_rows(data.others, "Compensation of employees", industry_names, mapping_105_to_64)
        gross_operating_surplus_and_mixed_income = clean_rows(data.others, "Gross operating surplus and mixed income", industry_names, mapping_105_to_64)
        #######################################

        final_consumption = clean_vector(data.input_output.final_consumption, industry_names, mapping_105_to_64)
        gross_fixed_captital_formation = clean_vector(data.input_output.gross_fixed_capital_formation, industry_names, mapping_105_to_64)
        delta_v_value_uk = clean_vector(data.input_output.delta_v_value_uk, industry_names, mapping_105_to_64)
        export_eu = clean_vector(data.input_output.exports_eu_to_uk, industry_names, mapping_105_to_64)
        export_world = clean_vector(data.input_output.export_world_to_uk, industry_names, mapping_105_to_64)
        total_use = clean_vector(data.input_output.total_use, industry_names, mapping_105_to_64)
        services_export = clean_vector(data.input_output.services_export, industry_names, mapping_105_to_64)


        nms = names(tax_products)

        low_income = DataFrame(low_income, nms)
        high_income = DataFrame(high_income, nms)
        mean_capital_current_year = DataFrame(mean_capital_current_year, nms)
        mean_capital_next_year = DataFrame(mean_capital_current_year, nms)

        high_income_share = high_income ./ (high_income .+ low_income)
        low_income_share = low_income ./ (high_income .+ low_income)


        import_export_matrix = clean_matrix(data.input_output.input_output_matrix, industry_names, mapping_105_to_64)

        return new(low_income, high_income, low_income_share, high_income_share, mean_capital_current_year, mean_capital_next_year, tax_products, tax_production, compensation_employees, gross_operating_surplus_and_mixed_income, final_consumption, gross_fixed_captital_formation, delta_v_value_uk, export_eu, export_world, total_use, services_export, import_export_matrix)
        
    end

end
