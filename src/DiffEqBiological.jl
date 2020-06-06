module DiffEqBiological

using Reexport
using MacroTools, LinearAlgebra, DataStructures
using RecipesBase, Latexify

using DiffEqBase, DiffEqJump, ModelingToolkit
@reexport using DiffEqBase, DiffEqJump, ModelingToolkit

import Base: (==)

const ExprValues = Union{Expr,Symbol,Float64,Int}

include("expression_utils.jl")
include("reaction_network.jl")

# reaction network macro
export @reaction_network
export @reaction_func

end # module

### Things possibly added back later. ###

#include("modify_reaction_network.jl")
#export @add_reactions, @add_reactions!

#include("network_properties.jl")

# functions to query network properties
#export species, params, speciesmap, paramsmap, numspecies, numreactions, numparams
#export oderhsfun, jacfun, paramjacfun, odefun, noisefun, sdefun, jumps, regularjumps
#export odeexprs, jacobianexprs, noiseexprs, jumpexprs, rateexpr, oderatelawexpr, ssaratelawexpr
#export substratestoich, productstoich, netstoich, ismassaction, dependants, dependents, substrates, products, substratesymstoich, productsymstoich
#export rxtospecies_depgraph, speciestorx_depgraph, rxtorx_depgraph
