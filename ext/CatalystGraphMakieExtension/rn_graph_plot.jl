#############
# Adapted from https://github.com/MakieOrg/GraphMakie.jl/issues/52#issuecomment-1018527479
#############

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
   edgeorder::Vector{Int64}
end

# Create the SimpleDiGraph corresponding to the species and reactions, the species-reaction graph
function MultiGraphWrap(rn::ReactionSystem) 
    srg = species_reaction_graph(rn)
    rateedges = Vector{Graphs.SimpleEdge{Int}}()
    sm = speciesmap(rn); specs = species(rn)

    deps = Set()
    for (i, rx) in enumerate(reactions(rn)) 
        empty!(deps)
        get_variables!(deps, rx.rate, specs)
        if !isempty(deps)
            for spec in deps 
                specidx = sm[spec]
                push!(rateedges, Graphs.SimpleEdge(specidx, i + length(specs)))
            end
        end
    end
    edgelist = vcat(collect(Graphs.edges(srg)), rateedges)
    edgeorder = sortperm(edgelist)
    MultiGraphWrap(srg, rateedges, edgeorder)
end

# Automatically set edge order if not supplied
function MultiGraphWrap(g::SimpleDiGraph{T}, multiedges::Vector{Graphs.SimpleEdge{T}}) where T
    edgelist = vcat(collect(Graphs.edges(g)), multiedges)
    edgeorder = sortperm(edgelist)
    MultiGraphWrap(g, multiedges, edgeorder)
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

function Graphs.edges(g::MultiGraphWrap)
    edgelist = vcat(collect(Graphs.edges(g.g)), g.multiedges)[g.edgeorder]
end

function gen_distances(g::MultiGraphWrap; inc = 0.4) 
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
    plot_network(rn::ReactionSystem)

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
"""  
function Catalyst.plot_network(rn::ReactionSystem)
    srg = MultiGraphWrap(rn)
    ns = length(species(rn))
    nodecolors = vcat([:skyblue3 for i in 1:ns], 
                      [:green for i in ns+1:nv(srg)])
    ilabels = vcat(map(s -> String(tosymbol(s, escape=false)), species(rn)),
                   fill("", nv(srg.g) - ns))
    nodesizes = vcat([30 for i in 1:ns],
                     [10 for i in ns+1:nv(srg)])

    ssm = substoichmat(rn); psm = prodstoichmat(rn)
    # Get stoichiometry of reaction
    edgelabels = map(Graphs.edges(srg.g)) do e
        string(src(e) > ns ? 
            psm[dst(e), src(e)-ns] :
            ssm[src(e), dst(e)-ns]) 
    end 
    edgecolors = [:black for i in 1:ne(srg)]

    num_e = ne(srg.g)
    for i in 1:length(srg.edgeorder)
        if srg.edgeorder[i] > num_e
            edgecolors[i] = :red
            insert!(edgelabels, i, "")
        end
    end

    graphplot(srg; 
              edge_color = edgecolors,
              elabels = edgelabels,
              elabels_rotation = 0,
              ilabels = ilabels,
              node_color = nodecolors,
              node_size = nodesizes,
              arrow_shift = :end,
              arrow_size = 20,
              curve_distance_usage = true,
              curve_distance = gen_distances(srg),
            )
end

"""
    plot_complexes(rn::ReactionSystem)

    Creates a GraphMakie plot of the [`ReactionComplex`](@ref)s in `rn`. Reactions
    correspond to arrows and reaction complexes to blue circles.

    Notes:
    - Black arrows from complexes to complexes indicate reactions whose rate is a
      parameter or a `Number`. i.e. `k, A --> B`.
    - Red arrows from complexes to complexes indicate reactions whose rate
    depends on species. i.e. `k*C, A --> B` for `C` a species.
"""
function Catalyst.plot_complexes(rn::ReactionSystem)
    img = incidencematgraph(rn); D = incidencemat(rn; sparse=true)
    specs = species(rn); rxs = reactions(rn)
    edgecolors = [:black for i in 1:length(rxs)]
    edgelabels = [repr(rx.rate) for rx in rxs]

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
        (!isempty(deps)) && (edgecolors[i] = :red)
    end

    # Resolve differences between reaction order and edge order in the incidencematgraph.
    rxorder = sortperm(edgelist); edgelist = edgelist[rxorder]
    multiedges = Vector{Graphs.SimpleEdge{Int}}()
    for i in 2:length(edgelist)
        isequal(edgelist[i], edgelist[i-1]) && push!(multiedges, edgelist[i])
    end
    img_ = MultiGraphWrap(img, multiedges)

    graphplot(img_;
              edge_color = edgecolors[rxorder],
              elabels = edgelabels[rxorder], 
              ilabels = complexlabels(rn), 
              node_color = :skyblue3,
              elabels_rotation = 0,
              arrow_shift = :end, 
              curve_distance_usage = true,
              curve_distance = gen_distances(img_)
            )
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
