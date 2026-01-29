# This file contains functions for generating the initial population of individuals for a genetic algorithm.
# It includes functions for creating the initial population, ensuring that the individuals are valid (i.e. they satisfy any constraints), and for converting the population between different data structures. 
# The main data structure used here is the `DimArray` type from the DimensionalData package, which is a multi-dimensional array with labeled dimensions and zero-cost labeled indexing. The two dimensions are `Genes`, for each parameter and species that can be tuned during optimzation, and `Individuals`, for each individual in the population.
# SlicedDimArray is a convenience type that is just the equivalent of ColumnSlices but with a DimArray as the parent, as the genetic algorithm expects a single dimension to iterate over.

function Base.similar(dimslices::SlicedDimArray)
    return eachslice(similar(parent(dimslices)), dims = Individuals)
end

function Base.copy(dimslices::SlicedDimArray)
    return eachslice(similar(parent(dimslices)), dims = Individuals)
end


"""
    make_dims(tunable_symbols::Vector{Symbol})

Make the `Genes` dimension for a `DimArray` that defines the symbolic indexing for the tunable parameters and initial conditions of a system.
"""
function make_dims(tunable_symbols::Vector{Symbol})
    cateogorical_lookup = DimensionalData.Categorical(tunable_symbols; order=DimensionalData.ForwardOrdered())

    # Here genes are the all the tunable input parameters and initial conditions for the system
    return Genes(cateogorical_lookup)
end

function make_dims(sys::ReactionSystem)
    return make_dims(get_tunable_symbols(sys))
end


"""
    initialize_default_population(default_values::Vector{Float64}, n::Int)

Initialize a population of `n` individuals with default values.
"""
function initialize_default_population(default_values::Vector{Float64}, n::Int)
    population = ones(length(default_values), n) .* default_values
    return population
end


function get_u0(dimvector::DimVector)
    return dimvector[DimensionalData.val(dimvector.metadata).u0]
end

function get_p(dimvector::DimVector)
    return dimvector[DimensionalData.val(dimvector.metadata).p]
end

"""
    get_u0_p(dimvector::DimVector)

Get the initial conditions and parameters from a DimVector individual to be fed into `remake` during ODE solving.
"""
function get_u0_p(dimvector::DimVector)
    u0 = get_u0(dimvector)
    p = get_p(dimvector)
    return (u0 = u0, p = p)
end


"""
    make_tunable_distributions(tunable_symbols, sys)

Make a Distributions.jl Product distribution from the bounds of the tunable symbols in the system.
"""
function make_tunable_distributions(tunable_symbols, sys)
    tunable_bounds_dictionary = get_tunable_bounds_dictionary(tunable_symbols, sys)
    return Product([Uniform(log10.(bound)...) for bound in tunable_bounds_dictionary])
end


"""
    generate_population(rx_sys::ReactionSystem, n::Int, fixed_inputs_dict::Dict, constraint=nothing)

Generate a population of `n` individuals. Each individual is sampled from a log-uniform distribution within the valid range for each parameter or initial condition.
If a constraint function is provided, it is used to validate the generated individuals.
Fixed inputs are excluded from the population generation process and resulting DimArray.
"""
function generate_population(rx_sys::ReactionSystem, n::Int, fixed_inputs_dict::Dict{Symbol, Float64} = Dict{Symbol, Float64}(), constraint_test = nothing; rng = default_rng())

    default_values = get_default_values(rx_sys, keys(fixed_inputs_dict))
    # @info "Default values: $default_values"
    # Initialize default population matrix
    population = initialize_default_population(default_values, n)

    # Get the tunable symbols excluding the fixed inputs
    unfixed_tunable_symbols = get_unfixed_tunable_symbols(fixed_inputs_dict, rx_sys)

    # Convert the initial population to a DimArray with the tunable symbols as the Genes dimension
    population_dimarray = convert_population_to_dimarray(population, unfixed_tunable_symbols, rx_sys, fixed_inputs_dict)

    # Create a multivariate distribution from the array of distributions
    multi_distribution = make_tunable_distributions(unfixed_tunable_symbols, rx_sys)

    # Make a vector of column slices of the population DimArray, because the genetic algorithm expects a single dimension to iterate over
    sliced_dimarray = eachslice(population_dimarray, dims = Individuals)

    if constraint_test !== nothing
        for ind in sliced_dimarray
            generate_valid_individual!(ind, multi_distribution, constraint_test; rng = rng)
        end
    else
        for ind in sliced_dimarray
            generate_valid_individual!(ind, multi_distribution; rng = rng)
        end
    end

    return sliced_dimarray
end

function generate_valid_individual!(individual::DimArray, tunable_input_distribution::Product; rng = default_rng())
    individual .= exp10.(rand(rng, tunable_input_distribution))
end

function generate_valid_individual!(individual::DimArray, tunable_input_distribution::Product, constraint_test; rng = default_rng())
    valid = false
    while !valid
        generate_valid_individual!(individual, tunable_input_distribution; rng = rng)
        if constraint_test(individual)
            valid = true
        end
    end
end


"""
    convert_population_to_dimarray(population, unfixed_tunable_symbols::Vector{Symbol}, rx_sys::ReactionSystem, fixed_inputs_dict::Dict{Symbol, Float64})

Convert a population of individuals to a DimArray with the tunable symbols as the Genes dimension.
"""
function convert_population_to_dimarray(population, unfixed_tunable_symbols::Vector{Symbol}, rx_sys::ReactionSystem, fixed_inputs_dict::Dict{Symbol, Float64})

    genes_dim = make_dims(unfixed_tunable_symbols)
    @info "Genes dimension: $genes_dim"

    tunable_species_syms = intersect(get_species_symbols(rx_sys), unfixed_tunable_symbols)
    @info "Tunable species symbols: $tunable_species_syms"
    tunable_parameter_syms = intersect(get_parameter_symbols(rx_sys), unfixed_tunable_symbols)
    @info "Tunable parameter symbols: $tunable_parameter_syms"

    # Need to make a new dictionary for the fixed inputs that replaces the value with a (value, column_idx) named tuple
    fixed_column_info_dict = Dict(key => (value, findfirst(==(key), get_tunable_symbols(rx_sys))) for (key, value) in fixed_inputs_dict)

    tunable_dim_selector_metadata = DimensionalData.Metadata(u0 = At(tunable_species_syms), p = At(tunable_parameter_syms), fixed_inputs_column_info = fixed_column_info_dict)
    dimarray = DimArray(population, (genes_dim, Individuals); metadata = tunable_dim_selector_metadata)
    
    return dimarray
end


# Helper function to convert population state to DataFrame
function dimarray_to_df(population::DimArray, objective_values)
    # Convert population DimArray to DataFrame 
    pop_df = DataFrame(population)  # Creates long format
    unstacked_df = unstack(pop_df, :Individuals, :Genes, :value)  # Convert to wide format

    # Select the results columns
    results_cols = dims(population, :Genes) |> collect
    selected_results_df = select!(unstacked_df, Cols(results_cols))

    # Add the fixed inputs columns to the results DataFrame
    fixed_inputs_column_info = population.metadata.val.fixed_inputs_column_info
    for (column_name, value_and_idx) in fixed_inputs_column_info
        fixed_value, column_idx = value_and_idx
        insertcols!(selected_results_df, column_idx, column_name => fixed_value)
    end
    
    # Add objective values
    insertcols!(selected_results_df, 1,
        :Fitness => objective_values[1,:],
        :Period => objective_values[2,:],
        :Amplitude => objective_values[3,:])
    
    return selected_results_df
end



# # Converts a SlicedDimArray to a raw matrix
# convert_dimarray_to_raw(population::SlicedDimArray) = collect.(eachcol(parent(parent(population))))

# # Returns the unique columns of a DimArray
# function unique_dimarray(dimarray::DimArray)
#     sliced_dimarray = eachslice(dimarray, dims = Individuals)
#     unique_idxs = unique(i -> sliced_dimarray[i], eachindex(sliced_dimarray))
#     return dimarray[:, unique_idxs]
# end


