module DiffEqBiological

using DataStructures, MacroTools, Parameters, Reexport, SparseArrays, SymEngine
using DiffEqBase, DiffEqJump, Compat
@reexport using DiffEqBase, DiffEqJump
using DynamicPolynomials, LinearAlgebra, RecipesBase, Latexify

@reexport using ModelingToolkit

import HomotopyContinuation

import Base: (==)

const ExprValues = Union{Expr,Symbol,Float64,Int}

include("expression_utils.jl")
include("reaction_network.jl")

include("MT_tmp/reaction_network.jl")
include("MT_tmp/maketype.jl")

# reaction network macro
export @reaction_network
export @reaction_func

end # module
