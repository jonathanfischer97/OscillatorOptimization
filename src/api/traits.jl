"""
    fitness_observables(::Type)

Trait that returns a tuple `(fft_observable, time_domain_observable)` identifying the
observables that the optimisation pipeline should use for a given `ReactionSystem`
sub-type.

* `fft_observable` - used for FFT-based frequency scoring (row 1 of the saved array)
* `time_domain_observable` - used for time-domain period/amplitude scoring (row 2)

Add methods for each model you want to optimise.  A reasonable default is provided
that throws an informative error so missing cases fail fast.
"""
fitness_observables(::Type) = error("No fitness_observables defined for this ReactionSystem type - extend OscTools.fitness_observables(::Type{typeof(your_model)})")

# ---------------------------------------------------------------------------
# Built-in mappings for the models distributed with OscTools
# ---------------------------------------------------------------------------

#= Replace the type-based dispatch (which collided because all ReactionSystems
   share the same concrete type) with value-based dispatch keyed on the
   ReactionSystem's `name` field. =#

# Convenience: dispatch via Val{symbol(name)} to avoid repeated conditionals
fitness_observables(rs::ReactionSystem) = fitness_observables(Val(rs.name))

# Model-specific specialisations
fitness_observables(::Val{:fullrn})      = (:Amem_old, :Amem)
fitness_observables(::Val{:trimer_rn})   = (:Amem_old, :TrimerYield) 