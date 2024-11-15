module CatalystStructuralIdentifiabilityExtension

# Fetch packages.
using Catalyst
import DataStructures.OrderedDict
import StructuralIdentifiability as SI

# Creates and exports make_si_ode function.
include("CatalystStructuralIdentifiabilityExtension/structural_identifiability_extension.jl")

end
