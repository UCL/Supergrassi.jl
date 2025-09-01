using Printf

"""
    function create_map_105_to_64(merge_codes::DataFrame, verbose::Bool = false)

Create a mapping from SIC 105 industry names to SIC 64 industry names.

# Arguments
- `merge_codes::DataFrame`: A DataFrame containing the mapping between SIC 105 and SIC 64 industry names. It should have two columns: `sic105` containing SIC 105 industry names and `sic64` containing SIC 64 industry names.
- `verbose::Bool`: If true, prints the initial and final industry names.

# Returns
- `Dict{String, String}`: A dictionary mapping SIC 105 industry names to SIC 64 industry names.
"""
function create_map_105_to_64(merge_codes::DataFrame, verbose::Bool = false)

    initial_industry_names = merge_codes[!, :sic105]
    final_industry_names = merge_codes[!, :sic64]

    if verbose
        println("Initial industry names: ", initial_industry_names)
        println("Final industry names: ", final_industry_names)
    end

    map_105_to_64 = Dict{String, String}()

    for i in eachindex(initial_industry_names)
        initial_name = initial_industry_names[i]
        final_name = final_industry_names[i]
        map_105_to_64[initial_name] = "SIC_64_" * @sprintf("%03i", final_name)
    end

    return map_105_to_64

end

"""
    function create_map_105_to_64(merge_codes::DataFrame, verbose::Bool = false)

Create a mapping from SIC 64 industry names to SIC 16 industry names.
# Arguments
- `merge_codes::DataFrame`: A DataFrame containing the mapping between SIC 64 and SIC 16 industry names. It should have two columns: `x1` containing SIC 64 industry names and `x7` containing SIC 16 industry names.
- `verbose::Bool`: If true, prints the initial and final industry names.

# Returns
- `Dict{String, String}`: A dictionary mapping SIC 64 industry names to SIC 16 industry names.
"""
function create_map_64_to_16(merge_codes::DataFrame, verbose::Bool = false)

    initial_industry_names = merge_codes[2:end, :x1]
    final_industry_names = merge_codes[2:end, :x7]

    if verbose
        println("Initial industry names: ", initial_industry_names)
        println("Final industry names: ", final_industry_names)
    end

    map_64_to_16 = Dict{String, String}()

    for i in eachindex(initial_industry_names)
        initial_name = "SIC_64_" * @sprintf("%03i", parse(Int64, initial_industry_names[i]))
        final_name = final_industry_names[i]
        map_64_to_16[initial_name] = final_name
    end

    return map_64_to_16

end

"""
    function reduce_columns_by_group_sum(df::DataFrame, mapping::Dict{String, String})

Reduce columns in a DataFrame by summing them based on a mapping.

# Arguments
- `df::DataFrame`: The DataFrame containing the columns to be reduced.
- `mapping::Dict{String, String}`: A dictionary mapping old column names to new column names.

# Returns
- `DataFrame`: A new DataFrame with reduced columns, where each new column is the sum of the old columns that map to it.
"""
function reduce_columns_by_group_sum(df::DataFrame, mapping::Dict{String, String})

    grouped = group_columns_by_new_name(mapping)

    # Sum columns per group
    new_cols = Dict{Symbol, Vector{eltype(df[!, 1])}}()
    for (new_name, old_syms) in grouped
        new_cols[Symbol(new_name)] = sum(eachcol(df[!, old_syms]))
    end

    return DataFrame(new_cols)
end

"""
    fubction reduce_columns_by_group_weighted_mean(df::DataFrame, mapping::Dict{String, String}; weights::DataFrame = DataFrame(ones(1,ncol(df)), names(df)))

Reduce columns in a DataFrame by calculating the weighted mean based on a mapping.
# Arguments
- `df::DataFrame`: The DataFrame containing the columns to be reduced.
- `mapping::Dict{String, String}`: A dictionary mapping old column names to new column names.
- `weights::DataFrame`: A DataFrame containing weights for each column. Defaults to a DataFrame of ones.

# Returns
- `DataFrame`: A new DataFrame with reduced columns, where each new column is the weighted mean of the old columns that map to it.
"""
function reduce_columns_by_group_weighted_mean(df::DataFrame, mapping::Dict{String, String}; weights::DataFrame = DataFrame(ones(1,ncol(df)), names(df)))

    grouped = group_columns_by_new_name(mapping)

    dd = df .* weights

    # Sum columns per group
    new_cols = Dict{Symbol, Vector{eltype(df[!, 1])}}()

    for (new_name, old_syms) in grouped
        new_cols[Symbol(new_name)] = sum(eachcol(dd[!, old_syms])) ./ sum(eachcol(weights[!, old_syms]))
    end

    return DataFrame(new_cols)

end

"""
    function group_columns_by_new_name(mapping::Dict{String, String})

Group columns by their new names based on a mapping.

# Arguments
- `mapping::Dict{String, String}`: A dictionary mapping old column names to new column names.

# Returns
- `Dict{String, Vector{Symbol}}`: A dictionary where keys are new column names and values are vectors of old column names (as Symbols) that map to them.
"""
function group_columns_by_new_name(mapping::Dict{String, String})

    grouped = Dict{String, Vector{Symbol}}()
    for (old, new) in mapping
        push!(get!(grouped, new, Vector{Symbol}()), Symbol(old))
    end

    return grouped

end

"""
    function group_dataframes(dfs::AbstractArray, col_names::AbstractArray, industry_names::AbstractArray, industries_on_cols::Bool = true, reduction_fun::Function = nothing, mapping::Dict{String, String} = Dict(); kwargs...)

Group multiple DataFrames by industry names and aggregate their values.

# Arguments
- `dfs::AbstractArray`: An array of DataFrames to be grouped.
- `col_names::AbstractArray`: An array of column names for the resulting DataFrame.
- `industry_names::AbstractArray`: An array of industry names to be used as row names.
- `industries_on_cols::Bool`: If true, industries will be on columns; otherwise, they will be on rows.
- `reduction_fun::Function`: A function to apply for data aggregation (default is `nothing`).
- `mapping::Dict{String, String}`: A dictionary mapping old column names to new column names.
- `kwargs...`: Additional keyword arguments to be passed to the `reduction_fun`.

# Returns
- `DataFrame`: A DataFrame with industries as rows or columns, depending on `industries_on_cols`, and aggregated values.
"""
function group_dataframes(dfs::AbstractArray, col_names::AbstractArray, industry_names::AbstractArray, industries_on_cols::Bool = true, reduction_fun::Function = nothing, mapping::Dict{String, String} = Dict(); kwargs...)

    dd = DataFrame([industry_names], ["industry"])
    for (df, col) in zip(dfs, col_names)

        # Do data aggregation if reduction function and maping are passed in
        if (!isnothing(reduction_fun) & !isempty(mapping))
            df = reduction_fun(df, mapping; kwargs...)
        end

        dd[!,col] = Float64.(permutedims(df).x1)

    end

    if industries_on_cols
        return permutedims(dd, :industry, "key")
    else
        return dd
    end

end


"""
    function combine_dataframe_row_wise(data::DataFrame, func::Function)

Combine rows of a DataFrame by applying a function to each row.

# Arguments
- `data::DataFrame`: The DataFrame containing the data to be combined.
- `func::Function`: A function to apply to each row of the DataFrame.

# Returns
- `DataFrame`: A new DataFrame with the combined results, where each column is named after the original column names.
"""
function combine_dataframe_row_wise(data::DataFrame, func::Function)
    return combine(data, names(data, Real) .=> func, renamecols=false)
end

