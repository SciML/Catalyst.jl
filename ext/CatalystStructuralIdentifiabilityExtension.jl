module CatalystStructuralIdentifiabilityExtension

# Fetch packages.
using Catalyst: Catalyst, ReactionSystem, conservationlaw_constants, conservedequations,
    ode_model
using ModelingToolkitBase: ModelingToolkitBase, System, complete, flatten, parameters,
    unknowns
using SymbolicUtils: substitute
using Symbolics: Equation, @variables
import DataStructures.OrderedDict
import StructuralIdentifiability as SI

# Creates and exports make_si_ode function.
include("CatalystStructuralIdentifiabilityExtension/structural_identifiability_extension.jl")

end
