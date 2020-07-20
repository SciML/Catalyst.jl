import Base.Iterators: flatten
using Catlab.Graphics.Graphviz
import Catlab.Graphics.Graphviz: Graph, Edge

graph_attrs = Attributes(:rankdir=>"LR")
node_attrs  = Attributes(:shape=>"plain", :style=>"filled", :color=>"white")
edge_attrs  = Attributes(:splines=>"splines")

function edgify(δ, i, reverse::Bool)
    attr = Attributes()
    return map(δ) do p
        val = String(p[1].op.name)
      weight = "$(p[2])"
      attr = Attributes(:label=>weight, :labelfontsize=>"6")
      return Edge(reverse ? ["rx_$i", "$val"] :
                            ["$val", "rx_$i"], attr)
    end
end

"""
    Graph(model::Model)

convert a Model into a GraphViz Graph. Transition are green boxes and states are blue circles. Arrows go from the input states to the output states for each transition.
"""
function Graph(model::ReactionSystem)
    rxs = reactions(model)
    statenodes = [Node(string(s.name), Attributes(:shape=>"circle", :color=>"#6C9AC3")) for s in species(model)]
    transnodes = [Node(string("rx_$i"), Attributes(:shape=>"square", :color=>"#E28F41")) for (i,r) in enumerate(rxs)]

    stmts = vcat(statenodes, transnodes)
    edges = map(enumerate(rxs)) do (i,r)
      vcat(edgify(zip(r.substrates,r.substoich), i, false),
           edgify(zip(r.products,r.prodstoich), i, true))
    end |> flatten |> collect
    stmts = vcat(stmts, edges)
    g = Graphviz.Graph("G", true, stmts, graph_attrs, node_attrs,edge_attrs)
    return g
end
