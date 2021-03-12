#######################################################################
# taken from Catlab.jl:
# https://raw.githubusercontent.com/AlgebraicJulia/Catlab.jl/master/src/graphics/Graphviz.jl
#######################################################################
""" AST and pretty printer for Graphviz's DOT language.

References:

- DOT grammar: http://www.graphviz.org/doc/info/lang.html
- DOT language guide: http://www.graphviz.org/pdf/dotguide.pdf
"""

# AST
#####

abstract type Expression end
abstract type Statement <: Expression end

""" AST type for Graphviz's "HTML-like" node labels.

For now, the HTML content is just a string.
"""
struct Html
  content::String
end
Base.print(io::IO, html::Html) = print(io, html.content)

const AttributeValue = Union{String,Html}
const Attributes = OrderedDict{Symbol,AttributeValue}

as_attributes(attrs::Attributes) = attrs
as_attributes(d::OrderedDict) = Attributes(Symbol(k) => d[k] for k in keys(d))
as_attributes(d::AbstractDict) =
  Attributes(Symbol(k) => d[k] for k in sort!(collect(keys(d))))

@with_kw_noshow struct Graph <: Expression
  name::String
  directed::Bool
  prog::String = "dot"
  stmts::Vector{Statement} = Statement[]
  graph_attrs::Attributes = Attributes()
  node_attrs::Attributes = Attributes()
  edge_attrs::Attributes = Attributes()
end

Graph(name::String, stmts::Vector{Statement}; kw...) =
  Graph(; name=name, directed=false, stmts=stmts, kw...)
Graph(name::String, stmts::Vararg{Statement}; kw...) =
  Graph(; name=name, directed=false, stmts=collect(stmts), kw...)
Digraph(name::String, stmts::Vector{Statement}; kw...) =
  Graph(; name=name, directed=true, stmts=stmts, kw...)
Digraph(name::String, stmts::Vararg{Statement}; kw...) =
  Graph(; name=name, directed=true, stmts=collect(stmts), kw...)

@with_kw_noshow struct Subgraph <: Statement
  name::String = "" # Subgraphs can be anonymous
  stmts::Vector{Statement} = Statement[]
  graph_attrs::Attributes = Attributes()
  node_attrs::Attributes = Attributes()
  edge_attrs::Attributes = Attributes()
end

Subgraph(stmts::Vector{Statement}; kw...) = Subgraph(; stmts=stmts, kw...)
Subgraph(stmts::Vararg{Statement}; kw...) = Subgraph(; stmts=collect(stmts), kw...)
Subgraph(name::String, stmts::Vector{Statement}; kw...) =
  Subgraph(; name=name, stmts=stmts, kw...)
Subgraph(name::String, stmts::Vararg{Statement}; kw...) =
  Subgraph(; name=name, stmts=collect(stmts), kw...)

struct Node <: Statement
  name::String
  attrs::Attributes
end
Node(name::String, attrs::AbstractDict) = Node(name, as_attributes(attrs))
Node(name::String; attrs...) = Node(name, attrs)

struct NodeID <: Expression
  name::String
  port::String
  anchor::String
  NodeID(name::String, port::String="", anchor::String="") = new(name, port, anchor)
end

struct Edge <: Statement
  path::Vector{NodeID}
  attrs::Attributes
end
Edge(path::Vector{NodeID}, attrs::AbstractDict) = Edge(path, as_attributes(attrs))
Edge(path::Vector{NodeID}; attrs...) = Edge(path, attrs)
Edge(path::Vararg{NodeID}; attrs...) = Edge(collect(path), attrs)
Edge(path::Vector{String}, attrs::AbstractDict) = Edge(map(NodeID, path), attrs)
Edge(path::Vector{String}; attrs...) = Edge(map(NodeID, path), attrs)
Edge(path::Vararg{String}; attrs...) = Edge(map(NodeID, collect(path)), attrs)

# Bindings
##########

""" Run a Graphviz program.

Invokes Graphviz through its command-line interface. If the `Graphviz_jll`
package is installed and loaded, it is used; otherwise, Graphviz must be
installed on the local system.

For bindings to the Graphviz C API, see the the package
[GraphViz.jl](https://github.com/Keno/GraphViz.jl). At the time of this writing,
GraphViz.jl is unmaintained.
"""
function run_graphviz(io::IO, graph::Graph; prog::Union{String,Nothing}=nothing,
  format::String="json0")
  if isnothing(prog)
    prog = graph.prog
  end
  @assert prog in ("dot", "neato", "fdp", "sfdp", "twopi", "circo")
  if USE_GV_JLL[]
    print("here!!!")
    fun = getfield(Graphviz_jll, Symbol(prog))
    prog = fun(identity)
  end
  open(`$prog -T$format`, io, write=true) do gv
    pprint(gv, graph)
  end
end
function run_graphviz(graph::Graph; kw...)
  io = IOBuffer()
  run_graphviz(io, graph; kw...)
  seekstart(io)
end

function Base.show(io::IO, ::MIME"image/svg+xml", graph::Graph)
  run_graphviz(io, graph, format="svg")
end

# Pretty-print
##############

""" Pretty-print the Graphviz expression.
"""
pprint(expr::Expression) = pprint(stdout, expr)
pprint(io::IO, expr::Expression) = pprint(io, expr, 0)

function pprint(io::IO, graph::Graph, n::Int)
  indent(io, n)
  print(io, graph.directed ? "digraph " : "graph ")
  print(io, graph.name)
  println(io, " {")
  pprint_attrs(io, graph.graph_attrs, n + 2; pre="graph", post=";\n")
  pprint_attrs(io, graph.node_attrs, n + 2; pre="node", post=";\n")
  pprint_attrs(io, graph.edge_attrs, n + 2; pre="edge", post=";\n")
  for stmt in graph.stmts
    pprint(io, stmt, n + 2, directed=graph.directed)
    println(io)
  end
  indent(io, n)
  println(io, "}")
end

function pprint(io::IO, subgraph::Subgraph, n::Int; directed::Bool=false)
  indent(io, n)
  if isempty(subgraph.name)
    println(io, "{")
        else
    print(io, "subgraph ")
    print(io, subgraph.name)
    println(io, " {")
  end
  pprint_attrs(io, subgraph.graph_attrs, n + 2; pre="graph", post=";\n")
  pprint_attrs(io, subgraph.node_attrs, n + 2; pre="node", post=";\n")
  pprint_attrs(io, subgraph.edge_attrs, n + 2; pre="edge", post=";\n")
  for stmt in subgraph.stmts
    pprint(io, stmt, n + 2, directed=directed)
    println(io)
  end
  indent(io, n)
  print(io, "}")
end

function pprint(io::IO, node::Node, n::Int; directed::Bool=false)
  indent(io, n)
  print(io, node.name)
  pprint_attrs(io, node.attrs)
  print(io, ";")
end

function pprint(io::IO, node::NodeID, n::Int)
  print(io, node.name)
  if !isempty(node.port)
    print(io, ":")
    print(io, node.port)
  end
  if !isempty(node.anchor)
    print(io, ":")
    print(io, node.anchor)
end
end

function pprint(io::IO, edge::Edge, n::Int; directed::Bool=false)
  indent(io, n)
  for (i, node) in enumerate(edge.path)
    if i > 1
      print(io, directed ? " -> " : " -- ")
    end
    pprint(io, node, n)
  end
  pprint_attrs(io, edge.attrs)
  print(io, ";")
end

function pprint_attrs(io::IO, attrs::Attributes, n::Int=0;
                      pre::String="", post::String="")
  if !isempty(attrs)
    indent(io, n)
        print(io, pre)
    print(io, " [")
    for (i, (key, value)) in enumerate(attrs)
      if (i > 1) print(io, ",") end
      print(io, key)
      print(io, "=")
      print(io, value isa Html ? "<" : "\"")
      print(io, value)
      print(io, value isa Html ? ">" : "\"")
    end
    print(io, "]")
    print(io, post)
  end
end

indent(io::IO, n::Int) = print(io, " "^n)


#######################################################################
# following is adapted from Petri.jl
# https://github.com/mehalter/Petri.jl
#######################################################################

graph_attrs = Attributes(:rankdir => "LR")
node_attrs  = Attributes(:shape => "plain", :style => "filled", :color => "white")
edge_attrs  = Attributes(:splines => "splines")

function edgify(δ, i, reverse::Bool)
    attr = Attributes()
    return map(δ) do p        
        val = String(p[1].f.name)
      weight = "$(p[2])"
      attr = Attributes(:label => weight, :labelfontsize => "6")
      return Edge(reverse ? ["rx_$i", "$val"] :
                            ["$val", "rx_$i"], attr)
    end
end

# make distinguished edge based on rate constant
function edgifyrates(rxs, specs)
    es = Edge[]
    for (i, rx) in enumerate(rxs)
        deps = rx.rate isa Number ? Any[] : get_variables(rx.rate, specs) 
        for dep in deps
            val = String(dep.f.name)
            attr = Attributes(:color => "#d91111", :style => "dashed")
            e = Edge(["$val", "rx_$i"], attr)
            push!(es, e)
        end
    end
    es
end

"""
    Graph(rn::ReactionSystem)

Converts a [`ReactionSystem`](@ref) into a Graphviz graph.
Reactions correspond to small green circles, and species to blue circles.

Notes:
- Black arrows from species to reactions indicate reactants, and are labelled
  with their input stoichiometry.
- Black arrows from reactions to species indicate products, and are labelled
  with their output stoichiometry.
- Red arrows from species to reactions indicate that species is used within the
  rate expression. For example, in the reaction `k*A, B --> C`, there would be a
  red arrow from `A` to the reaction node. In `k*A, A+B --> C`, there would be
  red and black arrows from `A` to the reaction node.
- Requires Graphviz to be installed and commandline accessible.
"""
function Graph(rn::ReactionSystem)
    rxs   = reactions(rn)
    specs = species(rn)
    statenodes = [Node(string(s.f.name), Attributes(:shape => "circle", :color => "#6C9AC3")) for s in specs]
    transnodes = [Node(string("rx_$i"), Attributes(:shape => "point", :color => "#E28F41", :width => ".1")) for (i, r) in enumerate(rxs)]

    stmts = vcat(statenodes, transnodes)
    edges = map(enumerate(rxs)) do (i, r)
      vcat(edgify(zip(r.substrates, r.substoich), i, false),
           edgify(zip(r.products, r.prodstoich), i, true))
    end
    es = edgifyrates(rxs, specs)
    (!isempty(es)) && push!(edges, es)

    stmts = vcat(stmts, collect(flatten(edges)))
    g = Digraph("G", stmts; graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)
    return g
end


"""
    savegraph(g::Graph, fname, fmt="png")

Given a `Graph` generated by [`Graph`](@ref), save the graph to the file with
name `fname` and extension `fmt`.

Notes:
- `fmt="png"` is the default output format.
- Requires Graphviz to be installed and commandline accessible.
"""
function savegraph(g::Graph, fname, fmt="png")
    open(fname, "w") do io
        run_graphviz(io, g, format=fmt)
    end
    nothing
end
