#__precompile__()

#module DiffEqBiological
module tmpMod
using DiffEqJump

using DifferentialEquations #Added, need access to constant rate jump because it is used in the jump part of the ReactionNetwork. Probably a better way to get this though.

using Compat
abstract type AbstractReaction end

import DataStructures: OrderedDict

include("reactions.jl")
include("ReactionNetwork.jl")   # New stuff
include("problem.jl")

export GillespieProblem

export VariableRateReaction, Reaction, ReactionSet, build_jumps_from_reaction

export @reaction_network
export newODEProblem, newSDEProblem, newJumpProblem

#New exports
export @reaction_network_new
export ReactionNetwork
export hill, mm

end # module
