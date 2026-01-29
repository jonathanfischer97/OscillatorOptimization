"""
Implementation of Quality Diversity algorithm (uses GA)

The constructor takes following keyword arguments:

- `populationSize`: The size of the population
- `crossoverRate`: The fraction of the population at the next generation, not including elite children, that is created by the crossover function.
- `mutationRate`: Probability of chromosome to be mutated
- `ɛ`/`epsilon`: Positive integer specifies how many individuals in the current generation are guaranteed to survive to the next generation.
Floating number specifies fraction of population.
- `selection`: [Selection](@ref) function (default: [`tournament`](@ref))
- `crossover`: [Crossover](@ref) function (default: [`genop`](@ref))
- `mutation`: [Mutation](@ref) function (default: [`genop`](@ref))
- `after_op`: a function that is executed on each individual after mutation operations (default: `identity`)
- `metrics` is a collection of convergence metrics.
"""
struct QD{T1,T2,T3,T4, T5} <: AbstractOptimizer
    populationSize::Int
    crossoverRate::Float64
    mutationRate::Float64
    ɛ::Real
    selection::T1
    crossover::T2
    mutation::T3
    metrics::T4
    optsys::T5

    QD(; populationSize::Int=10000, crossoverRate::Float64=0.75, mutationRate::Float64=0.75,
        ɛ::Real=0, epsilon::Real=ɛ,
        # num_tournament_groups = 20,
        selection::T1=tournament(cld(populationSize, 20), select=argmax),
        crossover::T2=TPX,
        mutation::T3=PLM(1.0; pm = 0.75, η = 2),
        metrics = ConvergenceMetric[AbsDiff(-Inf)],
        optsys::T5) where {T1, T2, T3, T5} =
        new{T1,T2,T3,typeof(metrics),T5}(populationSize, crossoverRate, mutationRate, epsilon, selection, crossover, mutation, metrics, optsys)
end
Evolutionary.population_size(method::QD) = method.populationSize

"""
Callback function to check if the number of unique oscillatory points exceeds a threshold.
Used to potentially terminate the optimization early if sufficient diversity is achieved.

Args:
- trace: The optimization trace containing population history
- n_points: The threshold number of unique points

Returns:
Boolean indicating if the number of unique points exceeds the threshold
"""
function n_points_callback(trace::Evolutionary.OptimizationTrace, n_points)
    unique_oscillatory_points = Set{Vector{Float64}}()
    for tr in trace
        for col in eachcol(parent(tr.metadata["population"]))
            push!(unique_oscillatory_points, col) 
        end
    end
    return length(unique_oscillatory_points) > n_points
end


Evolutionary.default_options(method::QD) = (abstol=-Inf, reltol=-Inf, successive_f_tol = 0, iterations=5, parallelization = :thread, show_trace=true, show_every=1, store_trace=true, callback = (tr) -> n_points_callback(tr, Inf))
Evolutionary.summary(m::QD) = "QD[P=$(m.populationSize),x=$(m.crossoverRate),μ=$(m.mutationRate)]"
Evolutionary.show(io::IO,m::QD) = print(io, summary(m))
Evolutionary.ismultiobjective(obj::EvolutionaryObjective{<:Any,<:Any,<:DimArray,<:Val}) = false


"""QD state type that captures additional data from the objective function in `valarray`"""
mutable struct QDState{T} <: AbstractOptimizerState 
    fittestValue::Float64  #* fitness of the fittest individual
    fittestChromosome::T  #* fittest chromosome (vector of gene values)
    objective_values::Matrix{Float64} #* array to store fitness, period, and amplitude of the population
end  
Evolutionary.value(s::QDState) = s.fittestValue #return the fitness of the fittest individual
Evolutionary.minimizer(s::QDState) = s.fittestChromosome #return the fittest individual
get_objective_values(s::QDState) = s.objective_values #return the objective values of the population


"""Get the fitness values from the objective_values matrix. View of the first row."""
function get_fitness(objective_values::AbstractMatrix)
    return @view objective_values[1, :]
end

"""Get the periods from the objective_values matrix. View of the second row."""
function get_periods(objective_values::AbstractMatrix)
    return @view objective_values[2, :]
end

"""Get the amplitudes from the objective_values matrix. View of the third row."""
function get_amplitudes(objective_values::AbstractMatrix)
    return @view objective_values[3, :]
end



"""Initialization of my custom QD algorithm state that captures additional data from the objective function\n
    - `method` is the QD method\n
    - `options` is the options dictionary\n
    - `objfun` is the objective function\n
    - `population` is the initial population
"""
function Evolutionary.initial_state(method::QD, options, objfun, population) 

    #- Initialize the main output array
    objective_values = zeros((3, method.populationSize))
    fitvals = get_fitness(objective_values)

    #- Evaluate population fitness, period and amplitude
    # population_matrix = parent(parent(population))
    # @info "Type of objective function: ", typeof(objfun)
    # @info "Type of objective values: ", typeof(objective_values)
    # @info "Type of population: ", typeof(population)
    value!(objfun, objective_values, population)

    #- Get the maximum fitness and index of the fittest individual
    maxfit, maxfitidx = findmax(fitvals)

    #- Initialize the state object
    return QDState(maxfit, copy(population[maxfitidx]), objective_values)
end

"""Update state function that captures additional data from the objective function"""
function Evolutionary.update_state!(objfun, constraints, state::QDState, parents, method::QD, options, itr)
    populationSize = method.populationSize
    rng = options.rng
    offspring = deepcopy(parents) 

    fitvals = get_fitness(state.objective_values)


    #* select offspring via tournament selection
    selected = method.selection(fitvals, populationSize, rng=rng)

    #* perform mating with TPX
    recombine!(offspring, parents, selected, method, rng=rng)

    #* perform mutation with BGA
    mutate!(offspring, method, constraints, rng=rng) #* only mutate descendants of the selected

    #* calculate fitness, period, and amplitude of the population
    evaluate!(objfun, state.objective_values, offspring, constraints)

    #* select the best individual
    maxfit, maxfitidx = findmax(fitvals)
    state.fittestChromosome .= offspring[maxfitidx]
    state.fittestValue = maxfit

    #* replace population
    parents .= offspring

    return false
end


function Evolutionary.recombine!(offspring, parents, selected, method::QD;
                    rng::AbstractRNG=default_rng())
    n = length(selected)
    mates = ((i,i == n ? i-1 : i+1) for i in 1:2:n)
    for (i,j) in mates
        p1, p2 = parents[selected[i]], parents[selected[j]]
        if rand(rng) < method.crossoverRate
            offspring[i], offspring[j] = method.crossover(p1, p2, rng=rng)
        else
            offspring[i], offspring[j] = p1, p2
        end
    end
end

function Evolutionary.mutate!(population, method::QD, constraints;
                 rng::AbstractRNG=default_rng())
    n = length(population)
    for i in 1:n
        if rand(rng) < method.mutationRate
            method.mutation(population[i], rng=rng)
        end
        population[i] .= abs.(population[i])
        apply!(constraints, population[i])
    end
end



function Evolutionary.evaluate!(objfun, objective_values, population, constraints::WorstFitnessConstraints)
    # calculate fitness of the population
    value!(objfun, objective_values, population)
    # apply penalty to fitness
    penalty!(get_fitness(objective_values), constraints, population)
end



function Evolutionary.value!(obj::Evolutionary.EvolutionaryObjective{TC,TF,TX,Val{:threadprogress}},
    F::AbstractVector, xs::AbstractVector{TX}) where {TC,TF<:Real,TX}

    n = length(xs)
    @showprogress Threads.@threads for i in 1:n
        F[i] = value(obj, xs[i])
        end
    F
end

function Evolutionary.value!(obj::Evolutionary.EvolutionaryObjective{TC,TF,TX,Val{:threadprogress}},
            F::AbstractMatrix, xs::AbstractVector{TX}) where {TC,TF,TX}

    n = length(xs)
    @showprogress Threads.@threads for i in 1:n
        fv = view(F, :, i)
        value(obj, fv, xs[i])
    end
    F
end




function Evolutionary.recombine!(offspring::SlicedDimArray, parents::SlicedDimArray, selected, method::QD;
                    rng::AbstractRNG=default_rng())
    n = length(selected)
    @assert typeof(offspring) == typeof(parents)
    mates = [(i, i == n ? 1 : i+1) for i in 1:2:n]  # Adjusted for circular pairing
    # parents_tunable = view_tunable(parents)
    # offspring_tunable = view_tunable(offspring)

    # Pre-determine crossover events
    crossover_indices = findall(rand(rng, length(mates)) .< method.crossoverRate)

    # Perform crossover only on selected pairs
    for idx in crossover_indices
        i, j = mates[idx]
        p1, p2 = parents[selected[i]], parents[selected[j]]
        o1, o2 = method.crossover(p1, p2, rng=rng)
        offspring[i] .= o1
        offspring[j] .= o2
    end

    # Directly copy non-crossover pairs
    non_crossover_indices = setdiff(1:length(mates), crossover_indices)
    for idx in non_crossover_indices
        i, j = mates[idx]
        offspring[i] .= parents[selected[i]]
        offspring[j] .= parents[selected[j]]
    end
end



function Evolutionary.mutate!(population::SlicedDimArray, method::QD, constraints;
                 rng::AbstractRNG=default_rng())
    # tunable_population = view_tunable(population)
    mutation_indices = findall(rand(rng, length(population)) .< method.mutationRate)

    # Apply mutations only to selected individuals
    for idx in mutation_indices
        method.mutation(population[idx], rng=rng)
        population[idx] .= abs.(population[idx])
        apply!(constraints, population[idx])
    end

    # Apply constraints and abs to all individuals
    for ind in population
        ind .= abs.(ind)
        apply!(constraints, ind)
    end
end



function Evolutionary.EvolutionaryObjective(f::TC, x::TX,
                               F::TF,
                               ;eval::Symbol = :thread) where {TC, TX <: DimArray, TF <: AbstractArray{<:Real}}
    # defval = default_values(x)
    # @info "Creating EvolutionaryObjective with eval = $eval"
    EvolutionaryObjective{TC,TF,TX,Val{eval}}(f, F, x, 0)
end



function compute_constants!(population) end

function Evolutionary.evaluate!(objfun, objective_values, population::SlicedDimArray, constraints::WorstFitnessConstraints)
    # Need to do in-place computation of constants for entire population prior to threaded evaluation, to avoid race conditions
    # compute_constants!(parent(population))

    # Reset the objective_values to 0.0
    objective_values .= 0.0

    #- calculate fitness of the population
    Evolutionary.value!(objfun, objective_values, population)

    # @assert !any((col[2] != 0.0 || col[3] != 0.0) && col[1] == 0.0 for col in eachcol(objective_values))

    # apply penalty to fitness from constraints
    # penalty!(get_fitness(objective_values), constraints, population)
    # @assert !any((col[2] != 0.0 || col[3] != 0.0) && col[1] == 0.0 for col in eachcol(objective_values))
end


Evolutionary.default_values(x::DimArray) = x

function penalty!(fitness::AbstractVector{T}, c::WorstFitnessConstraints{T,F}, population::SlicedDimArray) where {T,F}
    worst = minimum(fitness)
    p = zeros(size(value(c, first(population))))
    for (i,x) in enumerate(population)
        cv = value(c, x)
        # p = zeros(size(cv))
        # x_tunable = view_tunable(x)
        if !isfeasible(c.bounds, x, cv)
            p .= zero(p)
            for (i,j) in enumerate(c.bounds.eqc)
                p[j] = cv[j] - c.bounds.valc[i] - eps()
            end
            for (i,j) in enumerate(c.bounds.ineqc)
                p[j] = c.bounds.σc[i]*(c.bounds.bc[i]-cv[j])
            end
            fitness[i] = worst - sum(abs, p)
        end
    end
    return fitness
end










