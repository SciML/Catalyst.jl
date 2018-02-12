__precompile__()

module DiffEqBiological

using Reexport
using DiffEqBase
using DiffEqJump
using SymEngine
using MacroTools
using DataStructures
@reexport using DiffEqJump

using Compat
abstract type AbstractReactionNetwork end

include("ReactionNetwork.jl")
include("maketype.jl")
include("problem.jl")

export AbstractReactionNetwork

export @reaction_network, @reaction_func
export ODEProblem, SDEProblem, JumpProblem

Reaction(args...) = error("""
 The old Reaction DSL is deprecated for a new
 macro-based DSL which supports parameters, regulation,
 etc. Please see the documentation for details
"""

end # module
