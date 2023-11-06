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
function Catalyst.make_si_ode(rs::ReactionSystem; measured_quantities=observed(rs), known_p = [], ignore_no_measured_warn=false)
    ignore_no_measured_warn || isempty(measured_quantities) && @warn "No measured quantity provided to the `measured_quantities` argument, any further identifiability analysis will likely fail."
    known_quantities = make_measured_quantities(rs, measured_quantities, known_p)
    return StructuralIdentifiability.preprocess_ode(convert(ODESystem, rs), known_quantities)[1]
end

# For input measured quantities, if this is not a vector of equations, convert it to a proper form.
function make_measured_quantities(rs::ReactionSystem, measured_quantities::Vector{T}, known_p::Vector{S}) where {T,S}
    (measured_quantities isa Vector{Symbol}) && (measured_quantities = [Catalyst._symbol_to_var(rs, sym) for sym in measured_quantities])
    (known_p isa Vector{Symbol}) && (known_p = [Catalyst._symbol_to_var(rs, sym) for sym in known_p])
    all_quantities = [measured_quantities; known_p]
    @variables t (___internal_observables(t))[1:length(all_quantities)]
    return Equation[(all_quantities[i] isa Equation) ? all_quantities[i] : (___internal_observables[i] ~ all_quantities[i]) for i in 1:length(all_quantities)] 
end

### Structural Identifiability Wrappers ###

# Local identifiability.
function StructuralIdentifiability.assess_local_identifiability(rs::ReactionSystem, args...; measured_quantities=observed(rs), known_p = Num[], kwargs...)
    return StructuralIdentifiability.assess_local_identifiability(Catalyst.make_si_ode(rs; measured_quantities, known_p), args...; kwargs...)
end

# Global identifiability.
function StructuralIdentifiability.assess_identifiability(rs::ReactionSystem, args...; measured_quantities=observed(rs), known_p = Num[], kwargs...)
    return StructuralIdentifiability.assess_identifiability(Catalyst.make_si_ode(rs; measured_quantities, known_p), args...; kwargs...)
end

# Identifiable functions.
function StructuralIdentifiability.find_identifiable_functions(rs::ReactionSystem, args...; measured_quantities=observed(rs), known_p = Num[], kwargs...)
    return StructuralIdentifiability.find_identifiable_functions(Catalyst.make_si_ode(rs; measured_quantities, known_p), args...; kwargs...)
end

