module CatalystStructuralIdentifiabilityExtension

# Fetch packages.
using Catalyst
import DataStructures.OrderedDict
import StructuralIdentifiability as SI

# Creates and exports hc_steady_states function.
include("CatalystStructuralIdentifiabilityExtension/structural_identifiability_extension.jl")

end
