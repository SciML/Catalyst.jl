#__precompile__()

#module DiffEqBiological
module tmpMod
using DiffEqJump
using DiffEqBase
using SymEngine
using DataStructures

using Compat
abstract type AbstractReaction end

import DataStructures: OrderedDict

include("reactions.jl")
include("ReactionNetwork.jl")   # New stuff
include("maketype.jl")
include("problem.jl")

export GillespieProblem

export VariableRateReaction, Reaction, ReactionSet, build_jumps_from_reaction

export @reaction_network
export newODEProblem, newSDEProblem, newJumpProblem

#New exports
export @reaction_network_new

end # module
