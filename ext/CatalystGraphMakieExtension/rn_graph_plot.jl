#################################
# Adapted from https://github.com/MakieOrg/GraphMakie.jl/issues/52#issuecomment-1018527479
#################################

"""
    MultiGraphWrap{T}

Wrapper intended to allow plotting of multiple edges. This is needed in the following cases: 
- For the species-reaction graph, multiple edges can exist when a reaction depends on some species for its rate, and if that species is produced by the reaction.
- For the complex graph, multiple edges can exist between a pair of nodes if there are multiple reactions between the same complexes. This might include a reversible pair of reactions - we might have three total edges if one reaction is reversible, and we have a separate reaction going from one complex to the other.

`gen_distances` sets the distances between the edges that allows multiple to be visible on the plot at the same time. 
"""
struct MultiGraphWrap{T} <: Graphs.AbstractGraph{T}
   g::SimpleDiGraph{T}
   multiedges::Vector{Graphs.SimpleEdge{T}}
   """Sets the drawing order of the edges. Needed because multiedges need to be consecutive to be drawn properly."""
   edgeorder::Vector{Int64} 
end

# Create the SimpleDiGraph corresponding to the species and reactions, the species-reaction graph
function SRGraphWrap(rn::ReactionSystem) 
    srg = species_reaction_graph(rn)
    multiedges = Vector{Graphs.SimpleEdge{Int}}()
    sm = speciesmap(rn)
    specs = species(rn)

    deps = Set()
    for (i, rx) in enumerate(reactions(rn)) 
        empty!(deps)
        get_variables!(deps, rx.rate, specs)
        if !isempty(deps)
            for spec in deps 
                specidx = sm[spec]
                has_edge(srg, specidx, i + length(specs)) ? 
                    push!(multiedges, Graphs.SimpleEdge(specidx, i + length(specs))) : 
                    add_edge!(srg, Graphs.SimpleEdge(specidx, i + length(specs))) 
            end
        end
    end
    edgelist = vcat(collect(Graphs.edges(srg)), multiedges)
    edgeorder = sortperm(edgelist)
    MultiGraphWrap(srg, multiedges, edgeorder)
end

# Automatically set edge drawing order if not supplied
function MultiGraphWrap(g::SimpleDiGraph{T}, multiedges::Vector{Graphs.SimpleEdge{T}}) where T
    edgelist = vcat(collect(Graphs.edges(g)), multiedges)
    edgeorder = sortperm(edgelist)
    MultiGraphWrap(g, multiedges, edgeorder)
end

# Return the multigraph and reaction order corresponding to the complex graph. The reaction order is the order of reactions(rn) that would match the edge order given by g.edgeorder.
function ComplexGraphWrap(rn::ReactionSystem) 
    img = incidencematgraph(rn)
    D = incidencemat(rn; sparse=true)
    specs = species(rn)
    rxs = reactions(rn)

    deps = Set()
    edgelist = Vector{Graphs.SimpleEdge{Int}}()
    rows = rowvals(D)
    vals = nonzeros(D)

    # Construct the edge order for reactions.
    for (i, rx) in enumerate(rxs)
        inds = nzrange(D, i)
        val = vals[inds]
        row = rows[inds]
        (sub, prod) = val[1] == -1 ? (row[1], row[2]) : (row[2], row[1])
        push!(edgelist, Graphs.SimpleEdge(sub, prod))

        empty!(deps)
        get_variables!(deps, rx.rate, specs)
    end

    rxorder = sortperm(edgelist)
    edgelist = edgelist[rxorder]
    multiedges = Vector{Graphs.SimpleEdge{Int}}()
    for i in 2:length(edgelist)
        isequal(edgelist[i], edgelist[i-1]) && push!(multiedges, edgelist[i])
    end
    MultiGraphWrap(img, multiedges), rxorder
end

Base.eltype(g::MultiGraphWrap) = eltype(g.g)
Graphs.edgetype(g::MultiGraphWrap) = edgetype(g.g)
Graphs.has_edge(g::MultiGraphWrap, s, d) = has_edge(g.g, s, d)
Graphs.has_vertex(g::MultiGraphWrap, i) = has_vertex(g.g, i)
Graphs.inneighbors(g::MultiGraphWrap{T}, i) where T = inneighbors(g.g, i)
Graphs.outneighbors(g::MultiGraphWrap{T}, i) where T = outneighbors(g.g, i)
Graphs.ne(g::MultiGraphWrap) = length(g.multiedges) + length(Graphs.edges(g.g))
Graphs.nv(g::MultiGraphWrap) = nv(g.g)
Graphs.vertices(g::MultiGraphWrap) = vertices(g.g)
Graphs.is_directed(::Type{<:MultiGraphWrap}) = true
Graphs.is_directed(g::MultiGraphWrap) = is_directed(g.g)
Graphs.is_connected(g::MultiGraphWrap) = is_connected(g.g)

function Graphs.adjacency_matrix(g::MultiGraphWrap) 
    adj = Graphs.adjacency_matrix(g.g)
    for e in g.multiedges
        adj[src(e), dst(e)] = 1
    end
    adj
end

function Graphs.edges(g::MultiGraphWrap)
    edgelist = vcat(collect(Graphs.edges(g.g)), g.multiedges)[g.edgeorder]
end

function gen_distances(g::MultiGraphWrap; inc = 0.2) 
    edgelist = edges(g)
    distances = zeros(length(edgelist))
    edgedict = Dict(edgelist[1] => [1])
    for (i, e) in enumerate(@view edgelist[2:end])
        if edgelist[i] != edgelist[i+1]
            edgedict[e] = [i+1]
        else
            push!(edgedict[e], i+1)
        end
    end

    for (edge, inds) in edgedict
        if haskey(edgedict, Edge(dst(edge), src(edge)))
            distances[inds[1]] != 0. && continue
            inds_ = edgedict[Edge(dst(edge), src(edge))]

            len = length(inds) + length(inds_)
            sp = -inc/2*(len-1)
            ep = sp + inc*(len-1)
            dists = collect(sp:inc:ep)
            distances[inds] = dists[1:length(inds)]
            distances[inds_] = -dists[length(inds)+1:end]
        else
            sp = -inc/2*(length(inds)-1)
            ep = sp + inc*(length(inds)-1)
            distances[inds] = collect(sp:inc:ep)
        end
    end
    distances
end

"""
    plot_network(rn::ReactionSystem; kwargs...)

Converts a [`ReactionSystem`](@ref) into a GraphMakie plot of the species reaction graph
(or Petri net representation). Reactions correspond to small green circles, and 
species to blue circles.

Notes:
- Black arrows from species to reactions indicate reactants, and are labelled
  with their input stoichiometry.
- Black arrows from reactions to species indicate products, and are labelled
  with their output stoichiometry.
- Red arrows from species to reactions indicate that species is used within the
  rate expression. For example, in the reaction `k*A, B --> C`, there would be a
  red arrow from `A` to the reaction node. In `k*A, A+B --> C`, there would be
  red and black arrows from `A` to the reaction node.

For a list of accepted keyword arguments to the graph plot, please see the [GraphMakie documentation](https://graph.makie.org/stable/#The-graphplot-Recipe).
"""  
function Catalyst.plot_network(rn::ReactionSystem; kwargs...)
    srg = SRGraphWrap(rn)
    ns = length(species(rn))
    nodecolors = vcat([:skyblue3 for i in 1:ns], 
                      [:green for i in ns+1:nv(srg)])
    ilabels = vcat(map(s -> String(tosymbol(s, escape=false)), species(rn)),
                   ["R$i" for i in 1:nv(srg)-ns])

    ssm = substoichmat(rn)
    psm = prodstoichmat(rn)
    # Get stoichiometry of reaction
    edgelabels = map(Graphs.edges(srg.g)) do e
        string(src(e) > ns ? 
            psm[dst(e), src(e)-ns] :
            ssm[src(e), dst(e)-ns]) 
    end 
    edgecolors = [:black for i in 1:ne(srg)]

    num_e = ne(srg.g)
    # Handle the rate edges
    for i in 1:length(srg.edgeorder)
        # If there are stoichiometry and rate edges from the same species to reaction
        if srg.edgeorder[i] > num_e 
            edgecolors[i] = :red
            insert!(edgelabels, i, "")
        elseif edgelabels[i] == "0"
            edgecolors[i] = :red
            edgelabels[i] = ""
        end
    end

    layout = if !haskey(kwargs, :layout) 
        Stress()
    end
    f = graphplot(srg; 
              layout,
              edge_color = edgecolors,
              elabels = edgelabels,
              elabels_rotation = 0,
              ilabels = ilabels,
              node_color = nodecolors,
              arrow_shift = :end,
              arrow_size = 20,
              curve_distance_usage = true,
              curve_distance = gen_distances(srg),
              kwargs...
            )

    f.axis.xautolimitmargin = (0.15, 0.15)
    f.axis.yautolimitmargin = (0.15, 0.15)
    hidedecorations!(f.axis)
    hidespines!(f.axis)
    f.axis.aspect = DataAspect()
    
    f
end

"""
    plot_complexes(rn::ReactionSystem; show_rate_labels = false, kwargs...)

Creates a GraphMakie plot of the [`Catalyst.ReactionComplex`](@ref)s in `rn`. Reactions
correspond to arrows and reaction complexes to blue circles.

Notes:
- Black arrows from complexes to complexes indicate reactions whose rate is a
  parameter or a `Number`. i.e. `k, A --> B`.
- Red arrows from complexes to complexes indicate reactions whose rate constants
depends on species. i.e. `k*C, A --> B` for `C` a species.
- The `show_rate_labels` keyword, if set to `true`, will annotate each edge
with the rate constant for the reaction.

For a list of accepted keyword arguments to the graph plot, please see the [GraphMakie documentation](https://graph.makie.org/stable/#The-graphplot-Recipe).
"""
function Catalyst.plot_complexes(rn::ReactionSystem; show_rate_labels = false, kwargs...)
    rxs = reactions(rn)
    specs = species(rn)
    edgecolors = [:black for i in 1:length(rxs)]
    edgelabels = [repr(rx.rate) for rx in rxs]

    deps = Set()
    for (i, rx) in enumerate(rxs)
        empty!(deps)
        get_variables!(deps, rx.rate, specs)
        (!isempty(deps)) && (edgecolors[i] = :red)
    end

    # Get complex graph and reaction order for edgecolors and edgelabels. rxorder gives the order of reactions(rn) that would match the edge order in edges(cg).
    cg, rxorder = ComplexGraphWrap(rn)

    layout = if !haskey(kwargs, :layout)  
        Stress()
    end
    f = graphplot(cg;
              layout,
              edge_color = edgecolors[rxorder],
              elabels = show_rate_labels ? edgelabels[rxorder] : [], 
              ilabels = complexlabels(rn), 
              node_color = :skyblue3,
              elabels_rotation = 0,
              arrow_shift = :end, 
              curve_distance_usage = true,
              curve_distance = gen_distances(cg),
              kwargs...
            )

    f.axis.xautolimitmargin = (0.15, 0.15)
    f.axis.yautolimitmargin = (0.15, 0.15)
    hidedecorations!(f.axis)
    hidespines!(f.axis)
    f.axis.aspect = DataAspect()

    f
end

function complexelem_tostr(e::Catalyst.ReactionComplexElement, specstrs) 
    if e.speciesstoich == 1
        return "$(specstrs[e.speciesid])"  
    else
        return "$(e.speciesstoich)$(specstrs[e.speciesid])"
    end
end

# Get the strings corresponding to the reaction complexes
function complexlabels(rn::ReactionSystem)
    labels = String[]

    specstrs = map(s -> String(tosymbol(s, escape=false)), species(rn))
    complexes, B = reactioncomplexes(rn)

    for complex in complexes
        if isempty(complex) 
            push!(labels, "∅")
        elseif length(complex) == 1
            push!(labels, complexelem_tostr(complex[1], specstrs))
        else
            elems = map(c -> complexelem_tostr(c, specstrs), complex)
            str = reduce((e1, e2) -> *(e1, " + ", e2), @view elems[2:end]; init = elems[1])
            push!(labels, str)
        end
    end
    labels
end
