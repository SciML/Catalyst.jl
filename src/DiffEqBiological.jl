__precompile__()

module DiffEqBiological

using DiffEqJump

using Compat
@compat abstact type AbstractReaction end

include("reactions.jl")
include("problem.jl")

export GillespieProblem

export VariableRateReaction, Reaction, ReactionSet, build_jumps_from_reaction

end # module
