using DataFrames
using Statistics

function select_year(data::DataFrame, year::Int64)
    return data[data.year .== year, 4:end]
end

function combine_dataframe_row_wise(data::DataFrame, func::Function)
    return combine(data, names(data, Real) .=> func, renamecols=false)
end


struct CleanData

    low_income::DataFrame
    high_income::DataFrame

    low_income_share::DataFrame
    high_income_share::DataFrame

    mean_capital_current_year::DataFrame
    mean_capital_next_year::DataFrame

    function CleanData(low_income::DataFrame, high_income::DataFrame, mean_capital_current_year::DataFrame, mean_capital_next_year::DataFrame)
        high_income_share = high_income ./ (high_income .+ low_income)
        low_income_share = low_income ./ (high_income .+ low_income)

        return new(low_income, high_income, low_income_share, high_income_share, mean_capital_current_year, mean_capital_next_year)
        
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

    return CleanData(low_income, high_income, capital_current_year, capital_next_year)
end