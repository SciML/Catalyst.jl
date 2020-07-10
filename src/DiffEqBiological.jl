module DiffEqBiological

using Reexport
using MacroTools, LinearAlgebra, DataStructures
using RecipesBase, Latexify

using DiffEqBase, DiffEqJump, ModelingToolkit
@reexport using DiffEqBase, DiffEqJump, ModelingToolkit

import Base: (==), merge!

const ExprValues = Union{Expr,Symbol,Float64,Int}

include("expression_utils.jl")
include("reaction_network.jl")

# reaction network macro
export @reaction_network, @add_reactions
export @reaction_func


# functions to query network properties
include("networkapi.jl")
export species, params, speciesmap, paramsmap, numspecies, numreactions, numparams
export make_empty_network, addspecies!, addparam!, addreaction!

# depreciated currently
#export oderhsfun, jacfun, paramjacfun, odefun, noisefun, sdefun, jumps, regularjumps
#export odeexprs, jacobianexprs, noiseexprs, jumpexprs, rateexpr, oderatelawexpr, ssaratelawexpr
#export substratestoich, productstoich, netstoich, ismassaction, dependants, dependents, substrates, products, substratesymstoich, productsymstoich
#export rxtospecies_depgraph, speciestorx_depgraph, rxtorx_depgraph
#export @add_reactions, @add_reactions!

end # module
