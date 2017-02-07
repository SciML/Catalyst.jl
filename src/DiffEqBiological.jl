module DiffEqBiological

using DiffEqJump

include("reactions.jl")
include("problem.jl")

export GillespieProblem

export Reaction, ReactionSet, build_jumps_from_reaction

end # module
