__precompile__()

module DiffEqBiological

using Reexport
using DiffEqBase                # New stuff
using SymEngine
using DataStructures            # New stuff
@reexport using DiffEqJump

using Compat
abstract type AbstractReaction end

include("reactions.jl")
include("ReactionNetwork.jl")   # New stuff
include("maketype.jl")          # New stuff
include("problem.jl")

export GillespieProblem

export VariableRateReaction, Reaction, ReactionSet, build_jumps_from_reaction, AbstractReaction

export @reaction_network

#New exports
export @reaction_network_new
export ODEProblem, SDEProblem, JumpProblem

end # module
