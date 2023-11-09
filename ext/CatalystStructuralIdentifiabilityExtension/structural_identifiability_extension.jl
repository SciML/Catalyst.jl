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
function Catalyst.make_si_ode(rs::ReactionSystem; measured_quantities = [], known_p = [], ignore_no_measured_warn=false)
    known_quantities = make_measured_quantities(rs, measured_quantities, known_p; ignore_no_measured_warn)
    return StructuralIdentifiability.preprocess_ode(convert(ODESystem, rs; expand_functions = true), known_quantities)[1]
end

# For input measured quantities, if this is not a vector of equations, convert it to a proper form.
function make_measured_quantities(rs::ReactionSystem, measured_quantities::Vector{T}, known_p::Vector{S}; ignore_no_measured_warn=false) where {T,S}
    ignore_no_measured_warn || isempty(measured_quantities) && @warn "No measured quantity provided to the `measured_quantities` argument, any further identifiability analysis will likely fail. You can disable this warning by setting `ignore_no_measured_warn=true`."
    all_quantities = [measured_quantities; known_p]
    all_quantities = [(quant isa Symbol) ? Catalyst._symbol_to_var(rs, quant) : quant for quant in all_quantities]
    @variables t (___internal_observables(t))[1:length(all_quantities)]
    return Equation[(all_quantities[i] isa Equation) ? all_quantities[i] : (___internal_observables[i] ~ all_quantities[i]) for i in 1:length(all_quantities)] 
end

### Structural Identifiability Wrappers ###

# Local identifiability.
function StructuralIdentifiability.assess_local_identifiability(rs::ReactionSystem, args...; measured_quantities = Num[], known_p = Num[], ignore_no_measured_warn=false, kwargs...)
    known_quantities = make_measured_quantities(rs, measured_quantities, known_p; ignore_no_measured_warn)
    osys = convert(ODESystem, rs; expand_functions = true)
    return StructuralIdentifiability.assess_local_identifiability(osys, args...; measured_quantities=known_quantities, kwargs...)
end

# Global identifiability.
function StructuralIdentifiability.assess_identifiability(rs::ReactionSystem, args...; measured_quantities = Num[], known_p = Num[], ignore_no_measured_warn=false, kwargs...)
    known_quantities = make_measured_quantities(rs, measured_quantities, known_p; ignore_no_measured_warn)
    osys = convert(ODESystem, rs; expand_functions = true)
    return StructuralIdentifiability.assess_identifiability(osys, args...; measured_quantities=known_quantities, kwargs...)
end

# Identifiable functions.
function StructuralIdentifiability.find_identifiable_functions(rs::ReactionSystem, args...; measured_quantities = Num[], known_p = Num[], ignore_no_measured_warn=false, kwargs...)
    known_quantities = make_measured_quantities(rs, measured_quantities, known_p; ignore_no_measured_warn)
    osys = convert(ODESystem, rs; expand_functions = true)
    return StructuralIdentifiability.find_identifiable_functions(osys, args...; measured_quantities=known_quantities, kwargs...)
end

