### Structural Identifiability Analysis ###

# Local identifiability.
function StructuralIdentifiability.assess_local_identifiability(rs::ReactionSystem, args...; measured_quantities=observed(rs), known_p = [], kwargs...)
    return StructuralIdentifiability.assess_local_identifiability(convert(ODESystem, rs), args...; measured_quantities=get_measured_quantities(rs, measured_quantities, known_p), kwargs...)
end

# Global identifiability.
function StructuralIdentifiability.assess_identifiability(rs::ReactionSystem, args...; measured_quantities=observed(rs), known_p = [], kwargs...)
    return StructuralIdentifiability.assess_identifiability(convert(ODESystem, rs), args...; measured_quantities=get_measured_quantities(rs, measured_quantities, known_p), kwargs...)
end

# For a reaction system, list of measured quantities and known parameters, generate a StructuralIdentifiability compatible ODE. 

# For input measured quantities, if this is not a vector of equations, convert it to a proper form.
get_measured_quantities(rs::ReactionSystem, measured_quantities::Vector{Equation}) = measured_quantities # If form is already a vector of equations.
function get_measured_quantities(rs::ReactionSystem, measured_quantities::Vector{T}, known_p::Vector{S}) where {T,S}
    (measured_quantities isa Vector{Symbol}) && (measured_quantities = [Catalyst._symbol_to_var(rs, sym) for sym in measured_quantities])
    (known_p isa Vector{Symbol}) && (known_p = [Catalyst._symbol_to_var(rs, sym) for sym in known_p])
    @variables t ___internal_observables[1:length(measured_quantities)](t)
    measured_quantities = [eq[1] ~ eq[2] for eq in zip(___internal_observables, measured_quantities)] 
    return [measured_quantities; known_p]
end

