### Lattice Reaction Network Structure ###
# Describes a spatial reaction network over a graph.
struct LatticeReactionSystem{S,T} # <: MT.AbstractTimeDependentSystem # Adding this part messes up show, disabling me from creating LRSs
    # Input values.
    """The reaction system within each compartment."""
    rs::ReactionSystem{S}
    """The spatial reactions defined between individual nodes."""
    spatial_reactions::Vector{T}
    """The graph on which the lattice is defined."""
    lattice::SimpleDiGraph{Int64}

    # Derived values.
    """The number of compartments."""
    num_verts::Int64
    """The number of edges."""
    num_edges::Int64
    """The number of species."""
    num_species::Int64
    """Whenever the initial input was a digraph."""
    init_digraph::Bool
    """Species that may move spatially."""
    spat_species::Vector{BasicSymbolic{Real}}
    """Parameters which values are tied to edges (adjacencies), e.g. (possibly) have a unique value at each edge of the system."""
    edge_parameters::Vector{BasicSymbolic{Real}}

    function LatticeReactionSystem(rs::ReactionSystem{S},
                                   spatial_reactions::Vector{T},
                                   lattice::DiGraph; init_digraph = true) where {S, T}
        (T <: AbstractSpatialReaction) || error("The second argument must be a vector of AbstractSpatialReaction subtypes.") # There probably some better way to ascertain that T has that type. Not sure how.
        spat_species = unique(vcat(spatial_species.(spatial_reactions)...))
        rs_edge_parameters = filter(isedgeparameter, parameters(rs))
        srs_edge_parameters = setdiff(vcat(parameters.(spatial_reactions)...), parameters(rs))
        edge_parameters = unique([rs_edge_parameters; srs_edge_parameters])

        foreach(sr -> check_spatial_reaction_validity(rs, sr; edge_parameters=edge_parameters), spatial_reactions)   
        return new{S,T}(rs, spatial_reactions, lattice, nv(lattice), ne(lattice), length(unique([species(rs); spat_species])), init_digraph, spat_species, edge_parameters)
    end
end
function LatticeReactionSystem(rs, srs, lat::SimpleGraph)
    return LatticeReactionSystem(rs, srs, DiGraph(lat); init_digraph = false)
end

### Lattice ReactionSystem Getters ###

# Get all species.
species(lrs::LatticeReactionSystem) = unique([species(lrs.rs); lrs.spat_species])
# Get all species that may be transported.
spatial_species(lrs::LatticeReactionSystem) = lrs.spat_species

# Get all parameters.
ModelingToolkit.parameters(lrs::LatticeReactionSystem) = unique([parameters(lrs.rs); lrs.edge_parameters])
# Get all parameters which values are tied to vertexes (compartments).
vertex_parameters(lrs::LatticeReactionSystem) = setdiff(parameters(lrs), edge_parameters(lrs))
# Get all parameters which values are tied to edges (adjacencies).
edge_parameters(lrs::LatticeReactionSystem) = lrs.edge_parameters

# Gets the lrs name (same as rs name).
ModelingToolkit.nameof(lrs::LatticeReactionSystem) = nameof(lrs.rs)

# Checks if a lattice reaction system is a pure (linear) transport reaction system.
is_transport_system(lrs::LatticeReactionSystem) = all(sr -> sr isa TransportReaction, lrs.spatial_reactions)
