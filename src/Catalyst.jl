module Catalyst

using Reexport, ModelingToolkit
@reexport using ModelingToolkit
import MacroTools
import Base: (==), merge!, merge
using Latexify

const ExprValues = Union{Expr,Symbol,Float64,Int}

include("expression_utils.jl")
include("reaction_network.jl")

# reaction network macro
export @reaction_network, @add_reactions
export @reaction_func


# functions to query network properties
include("networkapi.jl")
export species, params, reactions, speciesmap, paramsmap, numspecies, numreactions, numparams
export make_empty_network, addspecies!, addparam!, addreaction!
export dependants, dependents

# for Latex printing of ReactionSystems
include("latexify_recipes.jl")

# for making and saving graphs
import Base.Iterators: flatten
import Catlab.Graphics.Graphviz: Graph, Edge, Attributes, Node, Digraph, run_graphviz
include("graphs.jl")
export Graph, savegraph

end # module
