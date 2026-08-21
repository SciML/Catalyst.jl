module CatalystBifurcationKitExtension

# Fetch packages.
using Catalyst: Catalyst, ReactionSystem, conservationlaw_constants, isautonomous,
    ss_ode_model
using ModelingToolkitBase: ModelingToolkitBase, complete, get_iv
using Symbolics: Symbolics
import BifurcationKit as BK

# Extends BifurcationProblem to work for ReactionSystem.
include("CatalystBifurcationKitExtension/bifurcation_kit_extension.jl")

end
