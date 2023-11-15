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
```
"""
function Catalyst.make_si_ode(rs::ReactionSystem; measured_quantities = [], known_p = [], ignore_no_measured_warn=false, remove_conserved = true)
    osys, conseqs, _ = make_osys(rs; remove_conserved)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p, conseqs; ignore_no_measured_warn)
    return StructuralIdentifiability.preprocess_ode(osys, measured_quantities)[1]
end

### Structural Identifiability Wrappers ###

# Local identifiability.
function StructuralIdentifiability.assess_local_identifiability(rs::ReactionSystem, args...; measured_quantities = Num[], known_p = Num[], funcs_to_check = Vector(), remove_conserved = true, ignore_no_measured_warn=false, kwargs...)
    osys, conseqs, vars = make_osys(rs; remove_conserved)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p, conseqs; ignore_no_measured_warn)
    funcs_to_check = make_ftc(funcs_to_check, conseqs, vars)
    out = StructuralIdentifiability.assess_local_identifiability(osys, args...; measured_quantities, funcs_to_check, kwargs...)
    return make_output(out, funcs_to_check, reverse.(conseqs))
end

# Global identifiability.
function StructuralIdentifiability.assess_identifiability(rs::ReactionSystem, args...; measured_quantities = Num[], known_p = Num[], funcs_to_check = Vector(), remove_conserved = true, ignore_no_measured_warn=false, kwargs...)
    osys, conseqs, vars = make_osys(rs; remove_conserved)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p, conseqs; ignore_no_measured_warn)
    funcs_to_check = make_ftc(funcs_to_check, conseqs, vars)
    out = StructuralIdentifiability.assess_identifiability(osys, args...; measured_quantities, funcs_to_check, kwargs...)
    return make_output(out, funcs_to_check, reverse.(conseqs))
end

# Identifiable functions.
function StructuralIdentifiability.find_identifiable_functions(rs::ReactionSystem, args...; measured_quantities = Num[], known_p = Num[], remove_conserved = true, ignore_no_measured_warn=false, kwargs...)
    osys, conseqs, vars = make_osys(rs; remove_conserved)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p, conseqs; ignore_no_measured_warn)
    out = StructuralIdentifiability.find_identifiable_functions(osys, args...; measured_quantities, kwargs...)
    return vector_subs(out, reverse.(conseqs))
end

### Helper Functions ###

# From a reaction system, creates the corresponding ODESystem for SI application (and also compute the, later needed, conservation law equations and list of system symbols).
function make_osys(rs::ReactionSystem; remove_conserved=true)
    rs = Catalyst.expand_registered_functions(rs)
    osys = convert(ODESystem, rs; remove_conserved)
    vars = [states(rs); parameters(rs)]

    # Fixes conservation law equations. These cannot be computed for hierarchical systems (and hence this is skipped). If none is found, still have to put on the right form.
    if !isempty(Catalyst.get_systems(rs)) || !remove_conserved
        conseqs = Vector{Pair{Any, Any}}[]
    else
        conseqs = [ceq.lhs => ceq.rhs for ceq in conservedequations(rs)]
        isempty(conseqs) && (conseqs = Vector{Pair{Any, Any}}[])
    end
    return osys, conseqs, vars
end

# For input measured quantities, if this is not a vector of equations, convert it to a proper form.
function make_measured_quantities(rs::ReactionSystem, measured_quantities::Vector{T}, known_p::Vector{S}, conseqs; ignore_no_measured_warn=false) where {T,S}
    ignore_no_measured_warn || isempty(measured_quantities) && @warn "No measured quantity provided to the `measured_quantities` argument, any further identifiability analysis will likely fail. You can disable this warning by setting `ignore_no_measured_warn=true`."
    all_quantities = [measured_quantities; known_p]
    all_quantities = [(quant isa Symbol) ? Catalyst._symbol_to_var(rs, quant) : quant for quant in all_quantities]
    all_quantities = vector_subs(all_quantities, conseqs)
    @variables t (___internal_observables(t))[1:length(all_quantities)]
    return Equation[(all_quantities[i] isa Equation) ? all_quantities[i] : (___internal_observables[i] ~ all_quantities[i]) for i in 1:length(all_quantities)] 
end

# Creates the functions that we wish to check for identifiability (if none give, by default, a list of parameters and species). Also replaces conservation law equations in.
function make_ftc(funcs_to_check, conseqs, vars)
    isempty(funcs_to_check) && (funcs_to_check = vars)
    return vector_subs(funcs_to_check, conseqs)
end

# Replaces conservation law equations back in the output, and also sorts it according to their input order (defaults to [states; parameters] order).
function make_output(out, funcs_to_check, conseqs)
    funcs_to_check = vector_subs(funcs_to_check, conseqs)
    out = Dict(zip(vector_subs(keys(out), conseqs), values(out)))
    out = sort(out; by = x -> findfirst(isequal(x, ftc) for ftc in funcs_to_check))
    return out    
end

# For a vector of expressions and a conservation law, replaces the law in.
vector_subs(eqs, subs) = [substitute(eq, subs) for eq in eqs]