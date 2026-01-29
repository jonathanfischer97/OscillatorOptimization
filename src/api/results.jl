"""
    Results

A struct holding the results of an optimization run.
"""
struct Results 
    "The hitrate of oscillatory individuals in the final population."
    hitrate::Float64
    "A DataFrame holding the results of the optimization."
    df::DataFrame 
    "A DimArray holding the final population."
    dimarray::DimArray
end

"""
    extract_trace_data(results)

Extract population data, objective values, and lineage information from optimization trace.

Returns a NamedTuple containing:
- populations: Combined DimArray of all populations across generations
- objectives: Combined matrix of objective values across generations
- gen_numbers: Vector mapping each individual to its generation number
- indices: Vector mapping each individual to its position within its generation

Note: The returned populations and objectives matrices have individuals arranged columnwise,
with generations concatenated horizontally.
"""
function extract_trace_data(results)
    # Calculate total number of individuals across all generations
    total_individuals = sum(size(tr.metadata["population"], 2) for tr in results.trace)
    n_generations = length(results.trace)
    
    # Pre-allocate arrays
    # Store each generation's population/objectives before combining
    # populations = Vector{typeof(first(results.trace).metadata["population"])}(undef, n_generations)
    @debug "Type of first(results.trace).metadata[\"population\"]: $(typeof(first(results.trace).metadata["population"]))"
    # full_population = similar(first(results.trace).metadata["population"], total_individuals, size(first(results.trace).metadata["population"], 1))
    full_population = Matrix{Float64}(undef, size(first(results.trace).metadata["population"], 1), total_individuals)
    @debug "Type of full_population: $(typeof(full_population))"
    # @debug "Size of full_population: $(size(full_population))"
    # objective_values = Vector{typeof(first(results.trace).metadata["objective_values"])}(undef, n_generations)
    objective_values = Matrix{Float64}(undef, size(first(results.trace).metadata["objective_values"], 1), total_individuals)
    # Track generation number and position for each individual
    gen_numbers = Vector{Int}(undef, total_individuals)
    indices = Vector{Int}(undef, total_individuals)
     
    # Single pass through trace to extract all metadata
    current_idx = 1
    for (gen, tr) in enumerate(results.trace)
        pop_size = size(tr.metadata["population"], 2)
        range_idx = current_idx:(current_idx + pop_size - 1)
        
        # Store population and objective values for this generation
        # populations[gen] = tr.metadata["population"]  # Direct assignment for array elements
        # @debug "Type of tr.metadata[\"population\"]: $(typeof(tr.metadata["population"]))"
        @debug "Size of tr.metadata[\"population\"]: $(size(tr.metadata["population"]))"
        full_population[:, range_idx] .= tr.metadata["population"]
        # objective_values[gen] = tr.metadata["objective_values"]  
        objective_values[:, range_idx] .= tr.metadata["objective_values"]
        
        # Record generation number and position for each individual
        gen_numbers[range_idx] .= gen  # Broadcast needed for scalar-to-range assignment
        indices[range_idx] .= 1:pop_size  # Broadcast needed for range-to-range assignment
        
        current_idx += pop_size
    end
    
    # Combine data across generations using reduce
    # combined_populations = reduce(cat, populations; dims = Individuals)
    # combined_objectives = reduce(hcat, objective_values)

    full_population_dimarray = DimArray(full_population, (dims(first(results.trace).metadata["population"], 1), Individuals); metadata = metadata(first(results.trace).metadata["population"]))

    @debug "Size of full_population_dimarray: $(size(full_population_dimarray))"
    
    return (
        populations = full_population_dimarray,
        objectives = objective_values,
        gen_numbers = gen_numbers,
        indices = indices
    )
end

"""
    Results(results::EvolutionaryOptimizationResults{<:QD, <:DimArray, <:Any})

Construct a Results object from optimization results by:
1. Extracting and combining data from all generations
2. Identifying unique individuals across all generations
3. Creating a DataFrame with fitness metrics and lineage information
4. Verifying the validity of critical values (fitness, period, amplitude)

Returns a Results object containing:
- hitrate: Proportion of successful oscillatory solutions found
- df: DataFrame with unique solutions and their properties
- dimarray: Combined population data across all generations
"""
function Results(results::EvolutionaryOptimizationResults{<:QD, <:DimArray, <:Any})
    # Extract all trace data in a single pass
    trace_data = extract_trace_data(results)
    isempty(trace_data.populations) && return Results(0.0, DataFrame(), trace_data.populations)
    
    # Find unique individuals and create views of their data
    unique_data = extract_unique_data(trace_data)
    
    # Create DataFrame and verify critical values
    results_df = create_verified_dataframe(unique_data)
    
    # Calculate success rate
    hitrate = length(unique_data.indices) / results.f_calls
    println("Hitrate: ", hitrate)
    
    Results(hitrate, results_df, trace_data.populations)
end

"""
Extract unique individuals from the population and create views of their associated data.
"""
function extract_unique_data(trace_data)
    # sliced_dimarray = eachslice(trace_data.populations; dims = Individuals)
    sliced_dimarray = eachslice(trace_data.populations, dims = 2)
    unique_idxs = unique(i -> sliced_dimarray[i], eachindex(sliced_dimarray))

    @debug "Type of trace_data.populations: $(typeof(trace_data.populations))"
    
    (
        population = @view(trace_data.populations[:, unique_idxs]),
        objectives = @view(trace_data.objectives[:, unique_idxs]),
        gen_numbers = trace_data.gen_numbers[unique_idxs],
        indices = trace_data.indices[unique_idxs]
    )
end

"""
Create a DataFrame from unique individual data and verify critical values.
Throws AssertionError if any critical values are invalid.
"""
function create_verified_dataframe(unique_data)
    # Create initial DataFrame
    results_df = dimarray_to_df(unique_data.population, unique_data.objectives)
    
    # Add generation and index information
    insertcols!(results_df, 1,
        :Generation => unique_data.gen_numbers,
        :Index => unique_data.indices
    )
    
    # Verify critical values
    @assert all(unique_data.objectives .> 0.0) "Oscillatory objective values are zero"
    @assert all(results_df.Fitness .> 0.0) "Fitness values are zero"
    @assert all(results_df.Period .> 0.0) "Period values are zero"
    @assert all(results_df.Amplitude .> 0.0) "Amplitude values are zero"
    
    results_df
end






