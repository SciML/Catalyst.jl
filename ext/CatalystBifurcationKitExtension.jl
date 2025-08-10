module CatalystBifurcationKitExtension

# Fetch packages.
using Catalyst: Catalyst, ReactionSystem, NonlinearSystem, ModelingToolkit,
                isautonomous, get_iv, symmap_to_varmap, conservationlaw_constants,
                complete, conservationlaw_errorcheck, get_networkproperties
import BifurcationKit as BK

# Extends BifurcationProblem to work for ReactionSystem.
include("CatalystBifurcationKitExtension/bifurcation_kit_extension.jl")

end
