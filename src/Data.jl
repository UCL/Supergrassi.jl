
struct IncomeData
    high::DataFrame
    low::DataFrame
end

struct HouseHoldData
    income::IncomeData
    hours::IncomeData
end



struct Data

    household::HouseHoldData

    capital::DataFrame
    turnover::DataFrame

    inventory::DataFrame
    depreciation::DataFrame

    risk_free_rate::DataFrame
    assets::DataFrame
    model_results::DataFrame

    merge_codes_105::DataFrame

    # gdp::DataFrame

end