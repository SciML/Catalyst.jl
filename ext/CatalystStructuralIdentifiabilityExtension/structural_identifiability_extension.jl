### Structural Identifiability ODE Creation ###

# For a reaction system, list of measured quantities and known parameters, generate a StructuralIdentifiability compatible ODE. 
"""
    si_ode(rs::ReactionSystem; measured_quantities=observed(rs), known_p = Num[], ignore_no_measured_warn=false)

Creates a ODE system of the form used within the StructuralIdentifiability.jl package. The output system is compatible with all StructuralIdentifiability functions. 

Arguments:
- `rs::ReactionSystem`; The reaction system we wish to convert to an ODE.
- `measured_quantities=observed(rs)`: The quantities of the system we can measure. May either be equations (e.g. `x1 + x2`), or single species (e.g. `:x`). Defaults to the systems observables.
- `known_p = Num[]`: List of parameters which values are known. 
- `ignore_no_measured_warn=false`: If set to `true`, no warning is provided when the `measured_quantities` vector is empty. 

Example:
```julia
rs = @reaction_network begin
    (p,d), 0 <--> X
end
si_ode(rs; measured_quantities = [:X], known_p = [:p])

Notes:
This function is part of the StructuralIdentifiability.jl extension. StructuralIdentifiability.jl must be imported to access it.
```
"""
function Catalyst.make_si_ode(rs::ReactionSystem; measured_quantities = [], known_p = [], 
                              ignore_no_measured_warn=false, remove_conserved = true)
    # Creates a MTK ODESystem, and a list of measured quantities (there are equations).
    # Gives these to SI to create an SI ode model of its preferred form.
    osys, conseqs, _ = make_osys(rs; remove_conserved)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p, conseqs; ignore_no_measured_warn)
    return SI.preprocess_ode(osys, measured_quantities)[1]
end

### Structural Identifiability Wrappers ###

# Creates dispatch for SI's local identifiability analysis function.
function SI.assess_local_identifiability(rs::ReactionSystem, args...; measured_quantities = Num[], 
                                         known_p = Num[], funcs_to_check = Vector(), remove_conserved = true, 
                                         ignore_no_measured_warn=false, kwargs...)
    # Creates a ODESystem, list of measured quantities, and functions to check, of SI's preferred form.
    osys, conseqs, vars = make_osys(rs; remove_conserved)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p, conseqs; ignore_no_measured_warn)
    funcs_to_check = make_ftc(funcs_to_check, conseqs, vars)

    # Computes identifiability and converts it to a easy to read form.
    out = SI.assess_local_identifiability(osys, args...; measured_quantities, funcs_to_check, kwargs...)
    return make_output(out, funcs_to_check, reverse.(conseqs))
end

# Creates dispatch for SI's global identifiability analysis function.
function SI.assess_identifiability(rs::ReactionSystem, args...; measured_quantities = Num[], known_p = Num[], 
                                   funcs_to_check = Vector(), remove_conserved = true, ignore_no_measured_warn=false, 
                                   kwargs...)
    # Creates a ODESystem, list of measured quantities, and functions to check, of SI's preferred form.
    osys, conseqs, vars = make_osys(rs; remove_conserved)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p, conseqs; ignore_no_measured_warn)
    funcs_to_check = make_ftc(funcs_to_check, conseqs, vars)

    # Computes identifiability and converts it to a easy to read form.
    out = SI.assess_identifiability(osys, args...; measured_quantities, funcs_to_check, kwargs...)
    return make_output(out, funcs_to_check, reverse.(conseqs))
end

# Creates dispatch for SI's function to find all identifiable functions.
function SI.find_identifiable_functions(rs::ReactionSystem, args...; measured_quantities = Num[], 
                                        known_p = Num[], remove_conserved = true, ignore_no_measured_warn=false, 
                                        kwargs...)
    # Creates a ODESystem, and list of measured quantities, of SI's preferred form.
    osys, conseqs = make_osys(rs; remove_conserved)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p, conseqs; ignore_no_measured_warn)

    # Computes identifiable functions and converts it to a easy to read form.
    out = SI.find_identifiable_functions(osys, args...; measured_quantities, kwargs...)
    return vector_subs(out, reverse.(conseqs))
end

### Helper Functions ###

# From a reaction system, creates the corresponding MTK-style ODESystem for SI application 
# Also compute the, later needed, conservation law equations and list of system symbols (states and parameters).
function make_osys(rs::ReactionSystem; remove_conserved=true)
    # Creates the ODESystem corresponding to the ReactionSystem (expanding functions and flattening it).
    # Creates a list of the systems all symbols (states and parameters).
    rs = Catalyst.expand_registered_functions(flatten(rs))
    osys = convert(ODESystem, rs; remove_conserved)
    vars = [states(rs); parameters(rs)]

    # Computes equations for system conservation laws.
    # These cannot be computed for hierarchical systems (and hence this is skipped). 
    # If there are no conserved equations, the `conseqs` variable must still have the `Vector{Pair{Any, Any}}` type.
    if !remove_conserved
        conseqs = Vector{Pair{Any, Any}}[]
    else
        conseqs = [ceq.lhs => ceq.rhs for ceq in conservedequations(rs)]
        isempty(conseqs) && (conseqs = Vector{Pair{Any, Any}}[])
    end

    return osys, conseqs, vars
end

# Creates a list of measured quantities of a form that SI can read.
# Each measured quantity must have a form like:
# `obs_var ~ X` # (Here, `obs_var` is a variable, and X is whatever we can measure).
function make_measured_quantities(rs::ReactionSystem, measured_quantities::Vector{T}, known_p::Vector{S},
                                  conseqs; ignore_no_measured_warn=false) where {T,S}
    # Warning if the user didn't give any measured quantities.
    if ignore_no_measured_warn || isempty(measured_quantities) 
        @warn "No measured quantity provided to the `measured_quantities` argument, any further identifiability analysis will likely fail. You can disable this warning by setting `ignore_no_measured_warn=true`."
    end

    # Appends the known parameters to the measured_quantities vector. Converts any Symbols to symbolics.
    measured_quantities = [measured_quantities; known_p]
    measured_quantities = [(q isa Symbol) ? Catalyst._symbol_to_var(rs, q) : q for q in measured_quantities]
    measured_quantities = vector_subs(measured_quantities, conseqs)

    # Creates one internal observation variable for each measured quantity (`___internal_observables`).
    # Creates a vector of equations, setting each measured quantity equal to one observation variable.
    @variables t (___internal_observables(t))[1:length(measured_quantities)]
    return Equation[(q isa Equation) ? q : (___internal_observables[i] ~ q) for (i,q) in enumerate(measured_quantities)] 
end

# Creates the functions that we wish to check for identifiability.
# If no `funcs_to_check` are given, defaults to checking identifiability for all states and parameters.
# Also, for conserved equations, replaces these in (creating a system without conserved quantities).
# E.g. for `X1 <--> X2`, replaces `X2` with `Γ[1] - X2`.
# Removing conserved quantities makes SI's algorithms much more performant.
function make_ftc(funcs_to_check, conseqs, vars)
    isempty(funcs_to_check) && (funcs_to_check = vars)
    return vector_subs(funcs_to_check, conseqs)
end

# Processes the outputs to a better form.
# Replaces conservation law equations back in the output (so that e.g. Γ are not displayed).
# Sorts the output according to their input order (defaults to the `[states; parameters]` order).
function make_output(out, funcs_to_check, conseqs)
    funcs_to_check = vector_subs(funcs_to_check, conseqs)
    out = Dict(zip(vector_subs(keys(out), conseqs), values(out)))
    out = sort(out; by = x -> findfirst(isequal(x, ftc) for ftc in funcs_to_check))
    return out    
end

# For a vector of expressions and a conservation law, substitutes the law into every equation.
vector_subs(eqs, subs) = [substitute(eq, subs) for eq in eqs]