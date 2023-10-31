module CatalystHomotopyContinuationExtension

# Fetch packages.
using Catalyst
import StructuralIdentifiability: assess_identifiability, assess_local_identifiability


# Creates and exports hc_steady_states function.
include("CatalystStructuralIdentifiabilityExtension/structural_identifiability_extension.jl")

end