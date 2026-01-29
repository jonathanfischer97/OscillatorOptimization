"""
Welcome to OscillatorOptimization!

This package implements computational modeling of a biological oscillator that can be
tuned by adjusting the surface-area to volume ratio. The main functionality allows
you to optimize reaction systems using genetic algorithms to find parameter sets
that produce desired oscillatory behavior.

To use this file:
1. Run from command line: julia --project examples/main.jl
2. Or from Julia REPL: include("examples/main.jl")
"""

#= Step 1: Setup Environment =#
using OscillatorOptimization

begin 
    using LinearAlgebra
    BLAS.set_num_threads(1)

    using DataFrames
    using CSV
    using Dates
end

#= Step 2: Choose your model and parameters =#
# Available models (from simple to complex):
#- fullrn:          Unrestricted 16 variable model that extends base_rn with peripheral_rn
#- base_rn:         Base 12 variable model with only the core oscillator components, no cytsolic binding between enzymes and adaptor
#- trimer_rn:       fullrn extended with heterotrimer system via allowing adaptor binding for each of the monomers

# Set up the model with fixed parameters that won't be optimized
fixed_params = Dict{Symbol, Float64}(
    :DF => 20.0,  # Dimensionality Factor: dimensionless surface-area to volume ratio scaling factor
    # Add other fixed parameters as needed
)

# Create an optimization system
opt_sys = OptimizationReactionSystem(fullrn; fixed_params)

#= Step 3: Run Genetic Algorithm Optimization =#
results = run_optimization(10000, opt_sys;
    # Population and Generation Parameters
    mutationRate = 0.95,     # Probability of applying mutation to an individual (0-1)
    crossoverRate = 0.75,    # Probability of applying crossover between parents (0-1)
    
    # Polynomial Mutation (PLM) Parameters
    mutation_δ = 1.0,        # Mutation step size - larger values = bigger mutations
    pm = 0.25,              # Probability of mutating each gene
    η = 1,                  # Distribution index - lower values = more diverse mutations
    
    # Simulated Binary Crossover (SBX) Parameters
    sbx_pm = 0.3,           # Probability of gene-wise mutation during crossover
    sbx_η = 1,              # Distribution index for crossover - lower = more diverse offspring
    
    # Selection Parameters
    num_tournament_groups = 20,  # Number of tournament groups - affects selection pressure
    
    # Termination Conditions
    n_points = Inf,         # Stop when this many unique oscillatory solutions found
    
    # Optional Parameters
    seed = 1234,           # Random seed for reproducibility, uses MersenneTwister
    show_trace = false,    # Disable trace display to avoid Term.jl display issues
    show_every = 1,        # Update frequency for progress display
    parallelization = :thread  # Use multi-threading for evaluations and show progress bar in stdout
)

#= Step 4: Save Results =#
date_string = Dates.format(today(), "mm-dd-yy")
# git_string = gitdescribe(; warn=false)  # Requires DrWatson
CSV.write("results_$(date_string).csv", results.df)

#= Step 5: Basic Analysis =#
println("\nOptimization Results Summary:")
println("----------------------------")
println("Best fitness: ", maximum(results.df.Fitness))
println("Number of unique solutions: ", nrow(results.df))

#= 
Advanced Usage Examples:

# 1. Run with different surface-area to volume ratio
low_df_params = Dict(:DF => 200.0)
low_df_sys = OptimizationReactionSystem(fullrn; low_df_params)
low_df_results = run_optimization(1000, low_df_sys)

# 2. More explorative search (higher mutation, lower selection pressure)
explorative_results = run_optimization(1000, opt_sys;
    mutationRate = 0.98,
    mutation_δ = 2.0,
    num_tournament_groups = 40  # More groups = less selection pressure
)

# 3. More exploitative search (lower mutation, higher selection pressure)
exploitative_results = run_optimization(1000, opt_sys;
    mutationRate = 0.8,
    mutation_δ = 0.5,
    num_tournament_groups = 10  # Fewer groups = higher selection pressure
)
=#

"""
For more information:
- Documentation: https://jonathanfischer97.github.io/OscillatorOptimization/dev/
- Source: https://github.com/jonathanfischer97/OscillatorOptimization
- Authors: Jonathan Fischer, Margaret Johnson, Ezra Greenberg
"""


