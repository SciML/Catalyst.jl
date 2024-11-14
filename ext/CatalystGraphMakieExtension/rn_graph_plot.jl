#############
# Adapted from https://github.com/MakieOrg/GraphMakie.jl/issues/52#issuecomment-1018527479
#############

"""
    SRGraphWrap{T}

Wrapper for the species-reaction graph containing edges for rate-dependence on species. Intended to allow plotting of multiple edges. 
"""
struct SRGraphWrap{T} <: AbstractGraph{T}
   g::SimpleDiGraph{T}
   rateedges::Vector{SimpleEdge}
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
    edgelist = vcat(collect(Graphs.edges(g.g)), g.rateedges)
    edgeorder = sortperm(edgelist)
    edgelist = edgelist[edgeorder]
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
    PetriNet(rn::ReactionSystem)

    See the documentation for [`SRGraph`](@ref).
"""
function PetriNet(rn::ReactionSystem) 
    SRGraph(rn)
end

"""
    SRGraph(rn::ReactionSystem; interactive=false)

Converts a [`ReactionSystem`](@ref) into a GraphMakie plot of the species reaction graph.
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
- The `interactive` flag sets the ability to interactively drag nodes and edges in the generated plot. 
    Only allowed if `GLMakie` is the loaded Makie backend.
"""  
function SRGraph(rn::ReactionSystem; interactive = false) 
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

    elist = Graphs.edges(srg)
    for i in 2:length(elist)
        elist[i] == elist[i-1] && begin
            edgecolors[i] = :red
            insert!(edgelabels, i, "")
        end
    end

    f, ax, p = graphplot(srg; 
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

    interactive && begin
        deregister_interaction!(ax, :rectanglezoom)
        register_interaction!(ax, :ndrag, NodeDrag(p))
        register_interaction!(ax, :edrag, EdgeDrag(p))
    end
    display(f)
    f
end

# Create the SimpleDiGraph corresponding to the species and reactions
function SRGraphWrap(rn::ReactionSystem) 
    srg = speciesreactiongraph(rn)
    rateedges = Vector{SimpleEdge}()
    sm = speciesmap(rn); specs = species(rn)

    for (i, rx) in enumerate(reactions(rn)) 
        deps = get_variables(rx.rate, specs)
        if deps != Any[]
            for spec in deps 
                specidx = sm[spec]
                push!(rateedges, SimpleEdge(specidx, i + length(specs)))
            end
        end
    end
    SRGraphWrap(srg, rateedges)
end

"""
    ComplexGraph(rn::ReactionSystem; interactive=false)

    Creates a GraphMakie plot of the [`ReactionComplex`](@ref)s in `rn`. Reactions
    correspond to arrows and reaction complexes to blue circles.

    Notes:
    - Black arrows from complexes to complexes indicate reactions whose rate is a
      parameter or a `Number`. i.e. `k, A --> B`.
    - Red arrows from complexes to complexes indicate reactions whose rate
    depends on species. i.e. `k*C, A --> B` for `C` a species.
    - The `interactive` flag sets the ability to interactively drag nodes and edges in the generated plot.
    Only allowed if `GLMakie` is the loaded Makie backend.
"""
function ComplexGraph(rn::ReactionSystem; interactive = false) 
    img = incidencematgraph(rn)
    specs = species(rn); rxs = reactions(rn)
    edgecolors = [:black for i in 1:ne(img)]
    nodelabels = complexlabels(rn)
    edgelabels = [repr(rx.rate) for rx in rxs]

    for (i, rx) in enumerate(rxs)
        deps = get_variables(rx.rate, specs)
        if deps != Any[]
            edgecolors[i] = :red
        end
    end

    f, ax, p = graphplot(img; 
             edge_color = edgecolors,
             elabels = edgelabels, 
             ilabels = complexlabels(rn), 
             node_color = :skyblue3,
             elabels_rotation = 0,
             arrow_shift = :end, 
             curve_distance = 0.2
            )

    interactive && begin
        deregister_interaction!(ax, :rectanglezoom)
        register_interaction!(ax, :ndrag, NodeDrag(p))
        register_interaction!(ax, :edrag, EdgeDrag(p))
    end
    display(f)
    f
end

function complexelemtostr(e::Catalyst.ReactionComplexElement, specstrs) 
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
            push!(labels, complexelemtostr(complex[1], specstrs))
        else
            elems = map(c -> complexelemtostr(c, specstrs), complex)
            str = reduce((e1, e2) -> *(e1, " + ", e2), elems[2:end]; init = elems[1])
            push!(labels, str)
        end
    end
    labels
end
