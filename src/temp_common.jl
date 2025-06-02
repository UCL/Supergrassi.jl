struct CommonVariables
    year::Int64
    industry_names::Vector{String}
    aggregated_names::Vector{String}
    mapping_105_to_64::Dict{String, String}
    mapping_64_to_16::Dict{String, String}

    split::Float64
    industries_in_cols::Bool

    tax_products::DataFrame
    tax_production::DataFrame
    sic64::Vector{String}


    exchange_rates::ExchangeRates
    interest_rate::Float64

    total_imports_from_uk::TotalImports
    total_imports_from_all_sources::TotalImports
end

