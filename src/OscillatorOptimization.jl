module OscillatorOptimization 

    #- Reaction system modeling and ODE solving stuff
    using Catalyst: ReactionSystem, @reaction_network, @parameters, @species, @variables, @named, setdefaults!, ODESystem, @unpack, conservationlaw_constants, complete, species
    using OrdinaryDiffEq: ODEProblem, ODESolution, Rosenbrock23, Rodas4P, Rodas5P, ReturnCode, SciMLBase, solve, remake
    using SciMLBase: successful_retcode
    using DiffEqCallbacks
    using ModelingToolkit: setmetadata, VariableTunable, VariableBounds, istunable, parameters, unknowns, getbounds, tunable_parameters, AbstractSystem, extend, structural_simplify, independent_variable, hasbounds, get_defaults, varmap_to_vars, get_observed, observed, defaults, getname
    using SymbolicIndexingInterface: getu, getp, setp, setu, parameter_symbols, parameter_values, all_variable_symbols
    using Symbolics
    using ADTypes

    #- Population generation
    using DimensionalData: DimensionalData, YDim, XDim, Metadata, DimArray, DimVector, At, dims, metadata
    DimensionalData.@dim Genes XDim "Genes"
    DimensionalData.@dim Individuals YDim "Individuals"
    const SlicedDimArray = ColumnSlices{D, AX, S} where {D <: DimArray, AX, S}
    using Random: rand, default_rng, seed!, AbstractRNG, MersenneTwister, Xoshiro
    using Distributions: Uniform, Product

    #- Fitness function stuff
    using FFTW
    using Statistics: mean, std
    using StatsBase
    using Peaks

    #- Evolutionary optimization stuff
    using Evolutionary: Evolutionary, GA, AbstractOptimizerState, AbstractOptimizer, OptimizationTraceRecord, OptimizationResults, EvolutionaryOptimizationResults, ConstraintBounds, AbstractConstraints, WorstFitnessConstraints, optimize, value, isfeasible, PLM, SBX, tournament, ConvergenceMetric, AbsDiff, apply!
    import Evolutionary: show, value!, trace!, initial_state, update_state!, recombine!, mutate!, evaluate!, ismultiobjective, summary, minimizer, random, default_values, penalty!, EvolutionaryObjective, default_options

    #- Trace and logging stuff
    using Term
    using UnicodePlots: lineplot, lineplot!, scatterplot!, histogram, heatmap, annotate!, panel
    # using CairoMakie
    using PrettyTables
    using ProgressMeter
    using Dates

    #- Results and saving stuff
    using DataFrames: DataFrame, AbstractDataFrame, DataFrameRow, Not, Cols, Between, unstack, select!, insertcols!, unique!, select, nrow
    using CSV: File, read, write
    using CategoricalArrays: categorical
    using Dictionaries

    # Custom logrange function for 10^x
    function logrange(start, stop, steps)
        10 .^ range(log10(start), log10(stop), length=steps)
    end
    export logrange 

    # Re-export commonly used functions from dependencies
    export ODEProblem, ODESolution  # from OrdinaryDiffEq
    export read, write  # from CSV
    export BLAS  # from LinearAlgebra
    export DataFrame  # from DataFrames


    #* Source code loading
    # Catalyst ReactionSystem model definitions
    include("models/full_model.jl")
    export base_rn, peripheral_rn, fullrn, make_odeprob

    include("models/trimer_model.jl")
    export trimer_rn, make_trimer_odeprob

    # Accessor utility functions for names and properties of the model
    include("models/accessor_functions.jl")
    export get_parameter_symbols, get_species_symbols, get_bounded_symbols, get_default_values, get_tunable_symbols, get_tunable_bounds_dictionary, get_unfixed_tunable_symbols, get_amem_symbols

    # ------------------------------------------------------------------
    # Model-specific traits and helper metadata
    # ------------------------------------------------------------------
    include("api/traits.jl")
    export fitness_observables

    # ODE solving and objective functions
    include("api/evaluate_individual.jl")
    export evaluate_individual!, solve_odes, remake_odeprob

    # Population generation of initial population prior to optimization
    include("api/population_generation.jl")
    export generate_population, Genes, Individuals, SlicedDimArray

    # Fitness function helper functions
    include("api/fitness_functions/helpers/fitness_function_helpers.jl")
    export compute_period, compute_amplitude, getFrequencies, getDif, getSTD

    # Fitness function definition
    include("api/fitness_functions/FitnessFunction.jl")
    export calculate_fitness, calculate_fitness!, find_amem_peaks, find_fft_peaks, find_amem_peaks_no_simd, check_oscillatory, check_oscillation_regularity

    # Evolutionary optimization overloads for quality-diversity and trace
    include("evolutionary_overloads/quality-diversity.jl")
    include("evolutionary_overloads/trace.jl")

    # Results and saving of optimization results
    include("api/results.jl")
    export Results

    # Constraints for optimization
    include("api/constraints.jl")
    export ConstraintSet, ConstraintFunction, null_constraint, km_constraint, kp_constraint, apply_constraint, make_WorstFitnessConstraints

    # Optimization loop
    include("api/optimization.jl")
    export run_optimization, OptimizationReactionSystem, ObjectiveFunction

    # Data handling
    include("utils/datahandling.jl")
    export read_all_csvs, solve_row
end