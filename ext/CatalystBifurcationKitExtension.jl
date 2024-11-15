module CatalystBifurcationKitExtension

# Fetch packages.
using Catalyst
import BifurcationKit as BK

# Extends BifurcationProblem to work for ReactionSystem.
include("CatalystBifurcationKitExtension/bifurcation_kit_extension.jl")

end
