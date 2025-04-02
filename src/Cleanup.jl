using DataFrames

function select_year(data::DataFrame, year::Int64)
    return data[data.year .== year, 4:end]
end

function sum_dataframe_row_wise(data::DataFrame)
    return combine(data, names(data, Real) .=> sum)
end


struct CleanData

    low_income::DataFrame
    high_income::DataFrame

    low_income_share::DataFrame
    high_income_share::DataFrame

    function CleanData(low_income::DataFrame, high_income::DataFrame)
        high_income_share = high_income ./ (high_income .+ low_income)
        low_income_share = low_income ./ (high_income .+ low_income)

        return new(low_income, high_income, low_income_share, high_income_share)
        
    end

end

function cleanup(data::Data, year::Int64)
    low_income = select_year(data.household.income.low, year)
    low_income = sum_dataframe_row_wise(low_income)

    high_income = select_year(data.household.income.high, year)
    high_income = sum_dataframe_row_wise(high_income)

    return CleanData(low_income, high_income)
end