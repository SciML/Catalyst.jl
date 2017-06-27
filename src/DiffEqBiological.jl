__precompile__()

module DiffEqBiological

using DiffEqJump

using Compat
@compat abstract type AbstractReaction end

import DataStructures: OrderedDict

include("reactions.jl")
include("problem.jl")

export GillespieProblem

export VariableRateReaction, Reaction, ReactionSet, build_jumps_from_reaction

export @reaction_network

end # module
