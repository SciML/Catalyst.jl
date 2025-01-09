module CatalystStructuralIdentifiabilityExtension

# Fetch packages.
using Catalyst
import DataStructures.OrderedDict
import StructuralIdentifiability as SI

@static if VERSION < v"1.11.2"
    @warn "Due to https://github.com/SciML/StructuralIdentifiability.jl/issues/360, CatalystStructuralIdentifiabilityExtension.jl is only tested with Julia v1.11.2 or later. Its use is not recommended on earlier Julia versions."
end

# Creates and exports make_si_ode function.
include("CatalystStructuralIdentifiabilityExtension/structural_identifiability_extension.jl")