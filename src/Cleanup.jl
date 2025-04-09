using DataFrames
using Statistics


safe_parse_float(x) = x isa AbstractString ? parse(Float64, x) : x


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

function get_rows(data::DataFrame, row_names::String, col_names::Array{String, 1}, mapping::Dict{String, Int64})
    rr = data[data.x2 .== row_names, 3:end]

    # # Remove trailing spaces from column names
    # col_names = strip.(col_names)
    rename!(rr, col_names)

    @show rr[!, :B05]
    rr = reduce_columns_by_group(rr, mapping)
    @show rr[!, "4"]

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

    function CleanData(low_income::DataFrame, high_income::DataFrame, mean_capital_current_year::DataFrame, mean_capital_next_year::DataFrame, others::DataFrame, industry_names::Array{String, 1}, mapping_105_to_64::Dict{String, Int64})
        high_income_share = high_income ./ (high_income .+ low_income)
        low_income_share = low_income ./ (high_income .+ low_income)

        tax_products = get_rows(others, "Taxes less subsidies on products", industry_names, mapping_105_to_64)
        tax_production = get_rows(others, "Taxes less subsidies on production", industry_names, mapping_105_to_64)
        compensation_employees = get_rows(others, "Compensation of employees", industry_names, mapping_105_to_64)
        gross_operating_surplus_and_mixed_income = get_rows(others, "Gross operating surplus and mixed income", industry_names, mapping_105_to_64)

        return new(low_income, high_income, low_income_share, high_income_share, mean_capital_current_year, mean_capital_next_year, tax_products, tax_production, compensation_employees, gross_operating_surplus_and_mixed_income)
        
    end

end

function cleanup(data::Data, year::Int64)
    low_income = select_year(data.household.income.low, year)
    low_income = combine_dataframe_row_wise(low_income, sum)

    high_income = select_year(data.household.income.high, year)
    high_income = combine_dataframe_row_wise(high_income, sum)

    capital_current_year = select_year(data.industry.capital, year)
    capital_current_year = combine_dataframe_row_wise(capital_current_year, mean)

    capital_next_year = select_year(data.industry.capital, year + 1)
    capital_next_year = combine_dataframe_row_wise(capital_next_year, mean)

    industry_names = data.input_output.industry_names

    mapping_105_to_64 = create_map_105_to_64(data)

    return CleanData(low_income, high_income, capital_current_year, capital_next_year, data.others, industry_names, mapping_105_to_64)
end

function create_map_105_to_64(data::Data)

    initial_industry_names = data.merge_codes_105[!, :sic105]
    final_industry_names = data.merge_codes_105[!, :sic64]


    println("Initial industry names: ", initial_industry_names)
    println("Final industry names: ", final_industry_names)

    map_105_to_64 = Dict{String, Int64}()

    for i in eachindex(initial_industry_names)
        initial_name = initial_industry_names[i]
        final_name = final_industry_names[i]
        map_105_to_64[initial_name] = final_name
    end
    
    return map_105_to_64

end

using DataFrames

function reduce_columns_by_group(df::DataFrame, mapping::Dict{String, Int64})
    # Group old columns by new name
    grouped = Dict{Int64, Vector{Symbol}}()
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
