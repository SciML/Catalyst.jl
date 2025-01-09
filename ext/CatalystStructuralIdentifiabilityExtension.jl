module CatalystStructuralIdentifiabilityExtension

# Fetch packages.
using Catalyst
import DataStructures.OrderedDict
import StructuralIdentifiability as SI

# Creates and exports make_si_ode function.
if VERSION >= v"1.11.2"
    include("CatalystStructuralIdentifiabilityExtension/structural_identifiability_extension.jl")
else
    @warn "Due to https://github.com/SciML/StructuralIdentifiability.jl/issues/360, CatalystStructuralIdentifiabilityExtension.jl is only compatible with Julia v1.11.2 or later."
end
