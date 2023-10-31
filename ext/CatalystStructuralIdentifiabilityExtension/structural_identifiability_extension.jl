### Structural Identifiability Analysis ###

# Local identifiability.
function StructuralIdentifiability.assess_local_identifiability(rs::ReactionSystem, args...; measured_quantities=observed(rs), kwargs...)
    return StructuralIdentifiability.assess_local_identifiability(convert(ODESystem, rs), args...; measured_quantities=get_measured_quantities(rs, measured_quantities), kwargs...)
end

# Global identifiability.
function StructuralIdentifiability.assess_identifiability(rs::ReactionSystem, args...; measured_quantities=observed(rs), kwargs...)
    return StructuralIdentifiability.assess_identifiability(convert(ODESystem, rs), args...; measured_quantities=get_measured_quantities(rs, measured_quantities), kwargs...)
end

# For input measured quantities, if this is not a vector of equations, convert it to a proper form.
get_measured_quantities(rs::ReactionSystem, measured_quantities::Vector{Equation}) = measured_quantities # If form is already a vector of equations.
function get_measured_quantities(rs::ReactionSystem, measured_quantities::Vector{Symbol})
    (measured_quantities isa Vector{Symbol}) && (measured_quantities = [Catalyst._symbol_to_var(rs, sym) for sym in measured_quantities])
    @variables t ___internal_observables[1:length(measured_quantities)](t)
    return [eq[1] ~ eq[2] for eq in zip(___internal_observables, measured_quantities)]           
end

