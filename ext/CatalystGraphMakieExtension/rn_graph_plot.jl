#############
# Adapted from https://github.com/MakieOrg/GraphMakie.jl/issues/52#issuecomment-1018527479
#############

"""
    SRGraphWrap{T}

Wrapper for the species-reaction graph containing edges for rate-dependence on species. Intended to allow plotting of multiple edges. 
"""
struct SRGraphWrap{T} <: Graphs.AbstractGraph{T}
   g::SimpleDiGraph{T}
   rateedges::Vector{Graphs.SimpleEdge{T}}
   edgeorder::Vector{Int64}
end

# Create the SimpleDiGraph corresponding to the species and reactions
function SRGraphWrap(rn::ReactionSystem) 
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
    SRGraphWrap(srg, rateedges, edgeorder)
end

Base.eltype(g::SRGraphWrap) = eltype(g.g)
Graphs.edgetype(g::SRGraphWrap) = edgetype(g.g)
Graphs.has_edge(g::SRGraphWrap, s, d) = has_edge(g.g, s, d)
Graphs.has_vertex(g::SRGraphWrap, i) = has_vertex(g.g, i)
Graphs.inneighbors(g::SRGraphWrap{T}, i) where T = inneighbors(g.g, i)
Graphs.outneighbors(g::SRGraphWrap{T}, i) where T = outneighbors(g.g, i)
Graphs.ne(g::SRGraphWrap) = length(g.rateedges) + length(Graphs.edges(g.g))
Graphs.nv(g::SRGraphWrap) = nv(g.g)
Graphs.vertices(g::SRGraphWrap) = vertices(g.g)
Graphs.is_directed(g::SRGraphWrap) = is_directed(g.g)

function Graphs.edges(g::SRGraphWrap)
    edgelist = vcat(collect(Graphs.edges(g.g)), g.rateedges)[g.edgeorder]
end

function gen_distances(g::SRGraphWrap; inc = 0.2) 
    edgelist = edges(g)
    distances = zeros(length(edgelist))
    for i in 2:Base.length(edgelist)
        edgelist[i] == edgelist[i-1] && (distances[i] = inc)
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
    srg = SRGraphWrap(rn)
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
             curve_distance = gen_distances(srg)
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
    img = incidencematgraph(rn)
    specs = species(rn); rxs = reactions(rn)
    edgecolors = [:black for i in 1:ne(img)]
    nodelabels = complexlabels(rn)
    edgelabels = [repr(rx.rate) for rx in rxs]

    deps = Set()
    for (i, rx) in enumerate(rxs)
        empty!(deps)
        get_variables!(deps, rx.rate, specs)
        (!isempty(deps)) && (edgecolors[i] = :red)
    end

    graphplot(img; 
             edge_color = edgecolors,
             elabels = edgelabels, 
             ilabels = complexlabels(rn), 
             node_color = :skyblue3,
             elabels_rotation = 0,
             arrow_shift = :end, 
             curve_distance = 0.2
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
