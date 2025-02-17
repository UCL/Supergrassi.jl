
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

struct InputOutputData
    # todo
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

    function Data(data_struct::Dict{String, DataFrame})

        household = HouseHoldData(
            IncomeData(data_struct["hi_income"], data_struct["lo_income"]),
            IncomeData(data_struct["hi_hours"], data_struct["lo_hours"])
        )     

        industry = IndustryData(
            data_struct["capital"],
            data_struct["turnover"],
            data_struct["inventory"]
        )

        depreciation = data_struct["depreciation"]
        risk_free_rate = data_struct["risk_free_rate"]
        assets = data_struct["assets"]
        model_results = data_struct["model_results"]
        merge_codes_105 = data_struct["merge_codes_105"]

        return new(household, industry, depreciation, risk_free_rate, assets, model_results, merge_codes_105)
    end

end

function organise_data(data_struct::Dict{String, DataFrame})

    data = Data(data_struct)

    return data
    
end