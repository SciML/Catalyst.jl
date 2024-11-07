module CatalystBifurcationKitExtension

# Fetch packages.
using Catalyst
import BifurcationKit as BK

# Creates and exports hc_steady_states function.
include("CatalystBifurcationKitExtension/bifurcation_kit_extension.jl")

end
