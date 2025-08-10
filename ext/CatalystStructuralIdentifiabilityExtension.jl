module CatalystStructuralIdentifiabilityExtension

# Fetch packages.
using Catalyst: Catalyst, ReactionSystem, ODESystem, ModelingToolkit, Equation,
                complete, conservationlaw_constants, conservedequations,
                species, unknowns, parameters, substitute, equations,
                observed, flatten, defaults, make_si_ode
import DataStructures: OrderedDict
import StructuralIdentifiability as SI

# Creates and exports make_si_ode function.
include("CatalystStructuralIdentifiabilityExtension/structural_identifiability_extension.jl")

end
