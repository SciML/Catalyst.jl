### Structural Identifiability Analysis ###

# Local identifiability.
function StructuralIdentifiability.assess_local_identifiability(rs::ReactionSystem, args...; measured_quantities=observed(rs), kwargs...)

    return StructuralIdentifiability.assess_local_identifiability(convert(ODESystem, rs), args...; measured_quantities=measured_quantities, kwargs...)
end

# Global identifiability.
function StructuralIdentifiability.assess_identifiability(rs::ReactionSystem, args...; measured_quantities=observed(rs), kwargs...)

    return StructuralIdentifiability.assess_identifiability(convert(ODESystem, rs), args...; measured_quantities=measured_quantities, kwargs...)
end