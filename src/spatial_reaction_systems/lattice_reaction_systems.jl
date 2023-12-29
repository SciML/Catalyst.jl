### Lattice Reaction Network Structure ###
# Describes a spatial reaction network over a graph.
# Adding the "<: MT.AbstractTimeDependentSystem" part messes up show, disabling me from creating LRSs.
struct LatticeReactionSystem{Q,R,S,T} # <: MT.AbstractTimeDependentSystem 
    # Input values.
    """The reaction system within each compartment."""
    rs::ReactionSystem{Q}
    """The spatial reactions defined between individual nodes."""
    spatial_reactions::Vector{R}
    """The graph on which the lattice is defined."""
    lattice::S

    # Derived values.
    """The number of compartments."""
    num_verts::Int64
    """The number of edges."""
    num_edges::Int64
    """The number of species."""
    num_species::Int64

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

    """Whenever the initial input was a digraph."""
    init_digraph::Bool
    """
    A list of all the edges on the lattice. 
    The format depends on the type of lattice (Cartesian grid, grid, or graph).
    """
    edge_list::T

    function LatticeReactionSystem(rs::ReactionSystem{Q}, spatial_reactions::Vector{R},
                                   lattice::S, num_verts, edge_list::T; init_digraph = true) where {Q, R, S, T}
        # Error checks.
        if !(R <: AbstractSpatialReaction) |
            error("The second argument must be a vector of AbstractSpatialReaction subtypes.") 
        end

        # Computes derived values spatial species. Counts the total number of species.
        if isempty(spatial_reactions)
            spat_species = Vector{BasicSymbolic{Real}}[]
        else
            spat_species = unique(reduce(vcat, [spatial_species(sr) for sr in spatial_reactions]))
        end

        # Computes the sets of vertex, edge, and all, parameters.
        rs_edge_parameters = filter(isedgeparameter, parameters(rs))
        if isempty(spatial_reactions)
            srs_edge_parameters = Vector{BasicSymbolic{Real}}[]
        else
            srs_edge_parameters = setdiff(reduce(vcat, [parameters(sr) for sr in spatial_reactions]), parameters(rs))
        end
        edge_parameters = unique([rs_edge_parameters; srs_edge_parameters])
        vertex_parameters = filter(!isedgeparameter, parameters(rs))

        # Counts the number of edges and species types (number of vertexes already counted).
        num_edges = false
        num_species = length(unique([species(rs); spat_species]))

        # Ensures the parameter order begins similarly to in the non-spatial ReactionSystem.
        ps = [parameters(rs); setdiff([edge_parameters; vertex_parameters], parameters(rs))]    

        # Checks that all spatial reactions are valid for this reactions system.
        foreach(sr -> check_spatial_reaction_validity(rs, sr; edge_parameters=edge_parameters), spatial_reactions)   

        return new{Q,R,S,T}(rs, spatial_reactions, lattice, num_verts, ne(lattice), num_species, 
                            spat_species, ps, vertex_parameters, edge_parameters, init_digraph, edge_list)
    end
end

# Creates a LatticeReactionSystem from a CartesianGrid lattice.
function LatticeReactionSystem(rs, srs, lattice::CartesianGridRej{S,T}; diagonal_connections=false) where {S,T}
    num_verts = prod(lattice.dims)
    edge_list = false
    return LatticeReactionSystem(rs, srs, lattice, num_verts, edge_list; init_digraph = false, edge_list)
end

# Creates a LatticeReactionSystem from a Boolean Array lattice.
function LatticeReactionSystem(rs, srs, lattice::Array{Bool, T}; diagonal_connections=false) where {T}
    num_verts = count(lattice)
    edge_list = false
    return LatticeReactionSystem(rs, srs, lattice, num_verts, edge_list; init_digraph = false, edge_list)
end

# Creates a LatticeReactionSystem from a DiGraph lattice.
function LatticeReactionSystem(rs, srs, lattice::SimpleGraph)
    num_verts = nv(lattice)
    edge_list = false
    return LatticeReactionSystem(rs, srs, lattice, num_verts, edge_list; init_digraph = false, edge_list)
end

# Creates a LatticeReactionSystem from a Graph lattice.
function LatticeReactionSystem(rs, srs, lattice::SimpleGraph)
    return LatticeReactionSystem(rs, srs, DiGraph(lattice); init_digraph = false)
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
is_transport_system(lrs::LatticeReactionSystem) = all(sr -> sr isa TransportReaction, lrs.spatial_reactions)

# Checks if a LatticeReactionSystem have a Cartesian grid lattice.
has_cartesian_lattice(lrs::LatticeReactionsSystem) = lrs.lattice isa CartesianGridRej{S,T} where {S,T}
# Checks if a LatticeReactionSystem have a regular grid lattice.
has_cartesian_lattice(lrs::LatticeReactionsSystem) = lrs.lattice isa Array{Bool, T} where T
# Checks if a LatticeReactionSystem have a graph lattice.
has_cartesian_lattice(lrs::LatticeReactionsSystem) = lrs.lattice isa SimpleDiGraph