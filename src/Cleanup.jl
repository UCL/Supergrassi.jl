using DataFrames
using Statistics

function select_year(data::DataFrame, year::Int64)
    return data[data.year .== year, 4:end]
end

function combine_dataframe_row_wise(data::DataFrame, func::Function)
    return combine(data, names(data, Real) .=> func, renamecols=false)
end

function get_rows(data::DataFrame, row_names::String, col_names::Array{String, 1})
    rr = data[data.x2 .== row_names, 4:end]
    rename!(rr, col_names)
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

    function CleanData(low_income::DataFrame, high_income::DataFrame, mean_capital_current_year::DataFrame, mean_capital_next_year::DataFrame, others::DataFrame, industry_names::Array{String, 1})
        high_income_share = high_income ./ (high_income .+ low_income)
        low_income_share = low_income ./ (high_income .+ low_income)

        tax_products = get_rows(others, "Taxes less subsidies on products", industry_names)
        tax_production = get_rows(others, "Taxes less subsidies on production", industry_names)

        compensation_employees = get_rows(others, "Compensation of employees", industry_names)
        gross_operating_surplus_and_mixed_income = get_rows(others, "Gross operating surplus and mixed income", industry_names)

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

    return CleanData(low_income, high_income, capital_current_year, capital_next_year, data.others, industry_names)
end