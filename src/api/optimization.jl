"""
OptimizationReactionSystem: Wrapper for components used in reaction system optimization.
Encapsulates the reaction system, ODE problem, and various utility functions.
"""
struct OptimizationReactionSystem{RS, OP, P, U, G, S, F, A, O} 
    rx_sys::RS                # Original reaction system
    oprob::OP                 # ODE problem derived from the reaction system
    alg::A                    # ODE algorithm, ex. Rodas5P(chunk_size = 16)
    parameter_setter::P       # Function to set parameters in the ODE problem
    species_setter::U         # Function to set species in the ODE problem
    amem_getter::G           # getter for the SymbolicIndex for A/membrane ratio observable
    amem_saver::S             # SavingCallback function to save A/membrane ratio observable to SavedValues
    rfft_plan::F #FFTW.rFFTWPlan{Float64, -1, false, 1, Tuple{Int64}} # Plan for real-valued FFT
    t::Vector{Float64}        # Time vector
    fixed_params::Dict{Symbol, Float64}       # Dictionary of fixed parameters
    constraint_set::ConstraintSet  # Set of constraints for the optimization
    obs_symbols::O                 # Tuple of observables (fft, time-domain)
end

# Utility function to extract the ODE system from an OptimizationReactionSystem
get_odesystem(opt_sys::OptimizationReactionSystem) = opt_sys.oprob.f.sys

"""
Constructor for OptimizationReactionSystem.
Converts a ReactionSystem to an ODESystem and sets up components for optimization.

Args:
- rx_sys: The reaction system to optimize
- constraint_set: Set of constraints for the optimization (default: null constraint)
- remove_conserved: Flag to remove conserved quantities (default: false)

Returns:
An OptimizationReactionSystem object
"""
function OptimizationReactionSystem(rx_sys::ReactionSystem, alg::T = Rodas5P(autodiff = AutoForwardDiff(chunksize=length(species(rx_sys)))), constraint_set::ConstraintSet = ConstraintSet([null_constraint()]);
        remove_conserved = false, fixed_params = Dict{Symbol, Float64}(), tspan = (0.0, 1968.2), dt = 0.1, abstol = 1e-7, reltol = 1e-6) where T
    
    # Convert ReactionSystem to ODESystem
    osys = convert(ODESystem, rx_sys; remove_conserved = remove_conserved, remove_conserved_warn = false) |> complete;

    # Create ODE problem
    oprob = ODEProblem{true, SciMLBase.FullSpecialize}(osys, Float64[], tspan; 
                jac = true, save_on = false, saveat = dt, verbose=false, abstol = abstol, reltol = reltol, remove_conserved = remove_conserved, remove_conserved_warn = false);

    time_vec = collect(range(oprob.tspan[1], oprob.tspan[2], step = dt))


    oprob = remake(oprob, p = fixed_params, u0 = copy(oprob.u0))

    # Create real-valued FFT plan
    # Needs to be a 2D array because I use Amem_old for the actual FFT fitness evaluation, but need the new Amem as well for the time-domain checks in the same fitness function.
    plan_array = @view reinterpret(reshape, Float64, fill((0.0, 0.0), length(time_vec)))[1, :]
    rfft_plan = plan_rfft(plan_array; flags = FFTW.PATIENT)

    # Create utility functions for parameter and species setting
    settable_parameter_symbols = setdiff(get_parameter_symbols(rx_sys), keys(fixed_params))
    # @info "Settable parameter symbols: $settable_parameter_symbols"
    settable_species_symbols = setdiff(get_species_symbols(rx_sys), keys(fixed_params))
    # @info "Settable species symbols: $settable_species_symbols"
    set_parameters! = setp(oprob, settable_parameter_symbols)
    set_species! = setu(oprob, settable_species_symbols)

    # Determine fitness observables via trait
    fft_sym, td_sym = fitness_observables(rx_sys)

    # Validation: ensure both symbols are actual observables in the model BEFORE building getu
    observed_eqns = get_observed(rx_sys)
    lhs_vars = getproperty.(observed_eqns, :lhs)
    observed_syms = unique(Symbol.(getname.(lhs_vars)))
    @assert all(s -> s in observed_syms, (fft_sym, td_sym)) "Trait returned observables ($(fft_sym), $(td_sym)) not found in model observables $(observed_syms)."

    @debug "OptimizationReactionSystem observables" (;fft_sym, td_sym, model_observables = observed_syms)

    # Create the get_observables function using getu
    get_obs = getu(osys, (fft_sym, td_sym))

    # Create the SavingCallback function to save the observables to SavedValues
    amem_saver = (u, t, integrator) -> get_obs(u)

    return OptimizationReactionSystem(rx_sys, oprob, alg, set_parameters!, set_species!, get_obs, amem_saver, rfft_plan, time_vec, fixed_params, constraint_set, (fft_sym, td_sym))
end

# Show method for OptimizationReactionSystem
function Base.show(io::IO, opt_sys::OptimizationReactionSystem)
    println(io, "> OptimizationReactionSystem:")
    println(io, "  Reaction System: ", opt_sys.rx_sys.name)
    println(io, "  Fixed Parameters: ", opt_sys.fixed_params)
    if !is_null_constraint(opt_sys.constraint_set)
        println(io, "  Constraint Set: ", opt_sys.constraint_set)
    end
end

"""
ObjectiveFunction: Wrapper for the objective function used in the genetic algorithm optimization.
"""
struct ObjectiveFunction{T<:OptimizationReactionSystem}
    opt_sys::T
end

(of::ObjectiveFunction)(phenotype, individual) = evaluate_individual!(phenotype, individual, of.opt_sys)

"""
Main function to run the genetic algorithm optimization.

Args:
- popsize: Size of the population
- opt_sys: OptimizationReactionSystem object
- mutation_δ, pm, η: Mutation parameters
- num_tournament_groups: Number of groups for tournament selection
- n_points: Threshold for unique oscillatory points
- rng: Random number generator
- fixed_params: Fixed parameters for the optimization

Returns:
A Results object with the hit rate, dataframe, and underlying DimArray
"""
function run_optimization(popsize, opt_sys::OptimizationReactionSystem; 
                          mutationRate = 0.95, mutation_δ = 1.0, pm = 0.75, η = 2, 
                          crossoverRate = 0.75, sbx_pm = 0.3, sbx_η = 2,
                          num_tournament_groups = 20, 
                          n_points = Inf, callback = (tr) -> n_points_callback(tr, n_points), seed = 1234, rng = MersenneTwister(seed), kwargs...)
    
    # @info "SEED NUMBER: $seed"
    @debug "RNG: $rng"

    # Record the start time so that we can track the duration of the optimization
    start_time = now()
    @info "STARTING OPTIMIZATION $(Dates.format(start_time, "mm-dd-yy I:MM:SS p"))"

    fixed_params_dict = opt_sys.fixed_params
    rx_sys = opt_sys.rx_sys
    cons = make_WorstFitnessConstraints(rx_sys, fixed_params_dict, opt_sys.constraint_set)

    # Generate initial population
    feasibility_test = (x) -> isfeasible(cons.bounds, x, value(cons, x))
    population = generate_population(rx_sys, popsize, fixed_params_dict, feasibility_test; rng = rng)

    # Display a sample of the initial population
    population_dimarray = parent(population)
    @debug begin
        io = IOBuffer()
        show(io, MIME("text/plain"), population_dimarray[Individuals = 1:5])
        String(take!(io))
    end
    
    # Verify that fixed parameters are indeed fixed in the population
    # @assert all([all(population_dimarray[Genes = At(param)] .== fixed_params_dict[param]) for param in keys(fixed_params_dict)])

    # Set up genetic algorithm components
    # PLM (Polynomial Mutation):
    # - mutation_δ: Controls the magnitude of mutations. Larger values lead to bigger mutations.
    # - pm: Probability of mutation for each gene. Higher values increase mutation frequency.
    # - η: Distribution index. Lower values produce offspring further from parents, higher values closer.
    # How it works: For each gene in an individual's vector:
    # 1. With probability pm, the gene is selected for mutation.
    # 2. If selected, a polynomial distribution is used to generate a perturbation.
    # 3. The perturbation is scaled by mutation_δ and added to the gene value.
    # 4. η controls the shape of the distribution: lower η allows more diverse mutations.
    mutator = PLM(mutation_δ; pm = pm, η = η)

    # SBX (Simulated Binary Crossover):
    # - crossoverRate (0.75): Probability of crossover occurring. Higher values increase genetic mixing.
    # - pm (0.3): Probability of mutation for each gene. Higher values increase mutation frequency.
    # - η (2 in SBX constructor): Distribution index. Lower values produce offspring further from parents.
    # How it works: For each pair of parent vectors:
    # 1. With probability crossoverRate, crossover is performed.
    # 2. If performed, for each gene position:
    #    a. A random number u is generated.
    #    b. A spread factor β is calculated based on u and η.
    #    c. Two child values are created by moving away from the parents' midpoint.
    # 3. η controls how similar children are to parents: lower η allows more diverse offspring.
    crossover = SBX(sbx_pm, sbx_η)
    
    # Tournament Selection:
    # - num_tournament_groups: Number of groups to split the population into.
    group_size = cld(popsize, num_tournament_groups)
    # How it works:
    # 1. For each selection, group_size individuals are randomly chosen, with replacement, into num_tournament_groups groups.
    # 2. The fittest individual from each group is selected.
    # 3. This process is repeated until the desired number of parents is selected.
    # 4. Smaller group_size decreases selection pressure, allowing more diversity.
    #    Larger groups increase selection pressure, favoring more representation of the fittest genotypes in parenting the next generation.
    selection = tournament(group_size, select=argmax)

    qd_method = QD(;populationSize=popsize, crossover = crossover, crossoverRate = crossoverRate, 
                   mutation = mutator, mutationRate = mutationRate, selection = selection, optsys = opt_sys)

    #- Configure optimization options
    # Filter kwargs to only include valid Evolutionary.Options parameters
    valid_option_kwargs = filter(kwargs) do (k,_)
        k in fieldnames(Evolutionary.Options)
    end

    # Configure optimization options with only valid kwargs
    options = Evolutionary.Options(;
        default_options(qd_method)..., 
        rng = rng, 
        callback = callback, 
        valid_option_kwargs...
    )

    @debug options

    # Define the objective function
    # objective_function! = (phenotype, individual) -> evaluate_individual!(phenotype, individual, opt_sys)
    objective_function! = ObjectiveFunction(opt_sys)

    # Run the optimization
    # The optimization process:
    # 1. Evaluate the fitness of each individual in the population.
    # 2. Select parents using tournament selection.
    # 3. Create offspring through SBX crossover.
    # 4. Mutate offspring using PLM.
    # 5. Evaluate new individuals and form the next generation.
    # 6. Repeat until termination criteria are met.
    evolutionary_results = Evolutionary.optimize(objective_function!, zeros(3), cons, qd_method, population, options)

    end_time = now()
    @info "OPTIMIZATION COMPLETE" Dates.format(end_time, "mm-dd-yy I:MM:SS p")

    # Calculate duration
    duration = end_time - start_time

    @info "Duration: $(Dates.canonicalize(duration))"

    # Return the Results object
    return Results(evolutionary_results)
end







