### Lattice Reaction Network Structure ###
# Describes a spatial reaction network over a graph.
# Adding the "<: MT.AbstractTimeDependentSystem" part messes up show, disabling me from creating LRSs.
struct LatticeReactionSystem{S, T} # <: MT.AbstractTimeDependentSystem 
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
    """
    All parameters related to the lattice reaction system
    (both with spatial and non-spatial effects).
    """
    parameters::Vector{BasicSymbolic{Real}}
    """
    Parameters which values are tied to vertexes (adjacencies), 
    e.g. (possibly) have an unique value at each vertex of the system.
    """
    vertex_parameters::Vector{BasicSymbolic{Real}}
    """
    Parameters which values are tied to edges (adjacencies), 
    e.g. (possibly) have an unique value at each edge of the system.
    """
    edge_parameters::Vector{BasicSymbolic{Real}}

    function LatticeReactionSystem(rs::ReactionSystem{S}, spatial_reactions::Vector{T},
            lattice::DiGraph; init_digraph = true) where {S, T}
        # There probably some better way to ascertain that T has that type. Not sure how.
        if !(T <: AbstractSpatialReaction)
            error("The second argument must be a vector of AbstractSpatialReaction subtypes.")
        end

        if isempty(spatial_reactions)
            spat_species = Vector{BasicSymbolic{Real}}[]
        else
            spat_species = unique(reduce(
                vcat, [spatial_species(sr) for sr in spatial_reactions]))
        end
        num_species = length(unique([species(rs); spat_species]))
        rs_edge_parameters = filter(isedgeparameter, parameters(rs))
        if isempty(spatial_reactions)
            srs_edge_parameters = Vector{BasicSymbolic{Real}}[]
        else
            srs_edge_parameters = setdiff(
                reduce(vcat, [parameters(sr) for sr in spatial_reactions]), parameters(rs))
        end
        edge_parameters = unique([rs_edge_parameters; srs_edge_parameters])
        vertex_parameters = filter(!isedgeparameter, parameters(rs))
        # Ensures the parameter order begins similarly to in the non-spatial ReactionSystem.
        ps = [parameters(rs); setdiff([edge_parameters; vertex_parameters], parameters(rs))]

        foreach(
            sr -> check_spatial_reaction_validity(rs, sr; edge_parameters = edge_parameters),
            spatial_reactions)
        return new{S, T}(
            rs, spatial_reactions, lattice, nv(lattice), ne(lattice), num_species,
            init_digraph, spat_species, ps, vertex_parameters, edge_parameters)
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
ModelingToolkit.parameters(lrs::LatticeReactionSystem) = lrs.parameters
# Get all parameters which values are tied to vertexes (compartments).
vertex_parameters(lrs::LatticeReactionSystem) = lrs.vertex_parameters
# Get all parameters which values are tied to edges (adjacencies).
edge_parameters(lrs::LatticeReactionSystem) = lrs.edge_parameters

# Gets the lrs name (same as rs name).
ModelingToolkit.nameof(lrs::LatticeReactionSystem) = nameof(lrs.rs)

# Checks if a lattice reaction system is a pure (linear) transport reaction system.
function is_transport_system(lrs::LatticeReactionSystem)
    all(sr -> sr isa TransportReaction, lrs.spatial_reactions)
end
