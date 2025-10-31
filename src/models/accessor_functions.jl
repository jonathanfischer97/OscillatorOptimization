using Symbolics
# function get_parameter_symbols(rs::ReactionSystem; filter_tunable = true)
#     if filter_tunable
#         return Symbolics.tosymbol.(tunable_parameters(rs, parameters(rs)); escape = false)
#     else
#         return Symbolics.tosymbol.(parameters(rs); escape = false)
#     end
# end

# function get_species_symbols(rs::ReactionSystem; filter_tunable = true)
#     if filter_tunable
#         return Symbolics.tosymbol.(tunable_parameters(rs, unknowns(rs)); escape = false)
#     else
#         return Symbolics.tosymbol.(unknowns(rs); escape = false)
#     end
# end


# function get_bounded_symbols(rs::ReactionSystem)
#     all_syms = [parameters(rs); unknowns(rs)]
#     filter!(s -> hasbounds(s), all_syms)
#     return Symbolics.tosymbol.(all_syms; escape = false)
# end

# get_tunable_symbolics(rx_sys) = tunable_parameters(rx_sys, [parameters(rx_sys); unknowns(rx_sys)])

# get_tunable_symbols(rx_sys) = Symbolics.tosymbol.(get_tunable_symbolics(rx_sys); escape = false)

# function get_default_values(rx_sys::ReactionSystem, fixed_symbols)
#     default_dict = get_defaults(rx_sys)
#     default_dict = Dict(Symbolics.tosymbol(k; escape = false) => v for (k,v) in default_dict)
#     tunable_symbols = setdiff(get_tunable_symbols(rx_sys), fixed_symbols)
#     return [default_dict[sym] for sym in tunable_symbols]
# end


"""
    get_parameter_symbols(rs::ReactionSystem; filter_tunable = true)

Returns a vector of Symbols corresponding to the parameters in the ReactionSystem, in the order they appear in the system.

# Arguments
- `rs::ReactionSystem`: The ReactionSystem to extract parameters from.
- `filter_tunable::Bool`: If true, only tunable parameters are returned (as defined by the `Tunable` Symbolics.jl metadata).

# Returns
- A vector of Symbols corresponding to the parameters in the ReactionSystem.
"""
function get_parameter_symbols(rs::ReactionSystem; filter_tunable = true)
    if filter_tunable
        return Symbolics.tosymbol.(tunable_parameters(rs, parameters(rs)); escape = false)
    else
        return Symbolics.tosymbol.(parameters(rs); escape = false)
    end
end

"""
    get_species_symbols(rs::ReactionSystem; filter_tunable = true)

Returns a vector of Symbols corresponding to the species in the ReactionSystem, in the order they appear in the system.

# Arguments
- `rs::ReactionSystem`: The ReactionSystem to extract species from.
- `filter_tunable::Bool`: If true, only tunable species are returned (as defined by the `Tunable` Symbolics.jl metadata).

# Returns
- A vector of Symbols corresponding to the species in the ReactionSystem.
"""
function get_species_symbols(rs::ReactionSystem; filter_tunable = true)
    if filter_tunable
        return Symbolics.tosymbol.(tunable_parameters(rs, unknowns(rs)); escape = false)
    else
        return Symbolics.tosymbol.(unknowns(rs); escape = false)
    end
end

"""
    get_tunable_symbolics(rs::ReactionSystem)

Returns a vector of the tunable Symbolics.jl objects (not native Symbols!) for the parameters and species in the ReactionSystem.
"""
get_tunable_symbolics(rs::ReactionSystem) = tunable_parameters(rs, [parameters(rs); unknowns(rs)])

get_tunable_symbols(rs::ReactionSystem) = Symbolics.tosymbol.(get_tunable_symbolics(rs); escape = false)

get_unfixed_tunable_symbols(fixed_inputs, rs::ReactionSystem) = setdiff(get_tunable_symbols(rs), fixed_inputs)

function get_unfixed_tunable_symbols(fixed_inputs_dict::Dict{Symbol, Float64}, rs::ReactionSystem)
    fixed_inputs = keys(fixed_inputs_dict)
    return get_unfixed_tunable_symbols(fixed_inputs, rs)
end

"""
    get_default_values(rs::ReactionSystem, fixed_symbols)

Returns a vector of the default values for the tunable symbols in the ReactionSystem, excluding the fixed symbols.
"""
function get_default_values(rs::ReactionSystem, fixed_symbols)
    default_dict = get_defaults(rs)
    default_dict = Dict(Symbolics.tosymbol(k; escape = false) => v for (k,v) in default_dict)
    tunable_symbols = get_unfixed_tunable_symbols(fixed_symbols, rs)
    return [default_dict[sym] for sym in tunable_symbols]
end

"""
    get_tunable_bounds_dictionary(tunable_symbols::Vector{Symbol}, rs::ReactionSystem)

Returns a Dictionaries.jl dictionary of the bounds for the tunable symbols in the ReactionSystem.
"""
function get_tunable_bounds_dictionary(tunable_symbols::Vector{Symbol}, rs::ReactionSystem)
    symbol_to_symbolic_dict = rs.var_to_name
    return dictionary([sym => getbounds(symbol_to_symbolic_dict[sym]) for sym in tunable_symbols])
end


# Returns a tuple of all the symbols with "LpA" in them (numerator) and all the symbols with "A" in them (denominator)
function get_amem_symbols(rs::AbstractSystem)
    # Get vector of symbol strings before filtering into numerator or denominator
    symbol_strings = string.(Symbolics.tosymbol.(all_variable_symbols(rs), escape = false)) 
    numerator_symbols = Symbol.(filter(x -> occursin("LpA", x), symbol_strings))
    denominator_symbols = Symbol.(filter(x -> occursin("A", x), symbol_strings))
    return (numerator_symbols, denominator_symbols)
end









    






    
