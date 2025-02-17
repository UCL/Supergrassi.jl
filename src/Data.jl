
struct IncomeData
    high::DataFrame
    low::DataFrame
end

struct HouseHoldData
    income::IncomeData
    hours::IncomeData
end

struct IndustryData
    capital::DataFrame
    turnover::DataFrame
    inventory::DataFrame
end

struct Data

    household::HouseHoldData
    industry::IndustryData

    depreciation::DataFrame

    risk_free_rate::DataFrame
    assets::DataFrame
    model_results::DataFrame

    merge_codes_105::DataFrame

    # gdp::DataFrame

end