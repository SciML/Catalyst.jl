__precompile__()

module DiffEqBiological

using DiffEqJump

abstract AbstractReaction

include("reactions.jl")
include("problem.jl")

export GillespieProblem

export VariableRateReaction, Reaction, ReactionSet, build_jumps_from_reaction

end # module
