### Lattice Reaction Network Structure ###
# Describes a spatial reaction network over a graph.
# Adding the "<: MT.AbstractTimeDependentSystem" part messes up show, disabling me from creating LRSs.
struct LatticeReactionSystem{R,S,T} # <: MT.AbstractTimeDependentSystem 
    # Input values.
    """The reaction system within each compartment."""
    rs::ReactionSystem{R}
    """The spatial reactions defined between individual nodes."""
    spatial_reactions::Vector{S}
    """The graph on which the lattice is defined."""
    lattice::T

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

    """
    Whenever the initial input was directed. True for Digraph input, false for other graphs and all grids.
    Used later to know whenever duplication of edge parameter should be duplicated
    (i.e. allow parameter for edge i => j to be used for j => i).
    Even if false, giving separate parameters for both directions is still permitted. 
    """
    directed_edges::Bool
    """
    A list of all the edges on the lattice. 
    The format depends on the type of lattice (Cartesian grid, grid, or graph).
    """
    edge_list::Vector{Vector{Int64}}

    function LatticeReactionSystem(rs::ReactionSystem{R}, spatial_reactions::Vector{S},
                                   lattice::S, num_verts, edge_list::T; directed_edges = false) where {R, S, T}
        # Error checks.
        if !(S <: AbstractSpatialReaction) |
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

        return new{R,S,T}(rs, spatial_reactions, lattice, num_verts, ne(lattice), num_species, 
                            spat_species, ps, vertex_parameters, edge_parameters, init_digraph, edge_list)
    end
end

# Creates a LatticeReactionSystem from a CartesianGrid lattice.
function LatticeReactionSystem(rs, srs, lattice::CartesianGridRej{S,T}; diagonal_connections=false) where {S,T}
    (length(lattice.dims) > 3) && error("Grids of higher dimension than 3 is currently not supported.")
    diagonal_connections && error("Diagonal connections is currently unsupported for Cartesian grid lattices.")
    num_verts = prod(lattice.dims)

    # Prepares the list of edges.
    edge_list = [Vector{Int64}() for i in 1:num_verts]

    # Ensures that the matrix is on a 3d form
    lattice = reshape(vals,dims..., fill(1,3- length(dims))...)


    edge_list = false
    return LatticeReactionSystem(rs, srs, lattice, num_verts, edge_list; edge_list)
end

# Creates a LatticeReactionSystem from a Boolean Array lattice.
function LatticeReactionSystem(rs, srs, lattice::Array{Bool, T}; diagonal_connections=false) where {T}  
    dims = size(vals)
    (length(dims) > 3) && error("Grids of higher dimension than 3 is currently not supported.")
    diagonal_connections && error("Diagonal connections is currently unsupported for grid lattices.")
    num_verts = count(lattice)

    # Prepares the list of edges.
    edge_list = [Vector{Int64}() for i in 1:num_verts]

    # Ensures that the matrix is on a 3d form
    lattice = CartesianGrid((lattice.dims..., fill(1,3- length(lattice.dims))...))

    # Loops through, simultaneously, the coordinates of each position in the grid, as well as that 
    # coordinate's (scalar) flat index.
    indices = [(i, j, k) for i in 1:lattice.dims[1], j in 1:lattice.dims[2], k in 1:lattice.dims[3]]
    flat_indices = 1:prod(lattice.dims)
    for ((i,j,k), idx) in zip(indices, flat_indices)
        for ii in i-1:i+1, jj in j-1:j+1, kk in k-1:k+1
            # Finds the (scalar) flat index of this neighbour. If it is a valid neighbour, add it to edge_list.
            n_idx = (k - 1) * (l * m) + (j - 1) * l + i
            (1 <= n_idx <= flat_indices[end]) || continue
            (n_idx == idx) && continue
            push!(edge_list[cur_vert], n_idx)
        end
    end

    return LatticeReactionSystem(rs, srs, lattice, num_verts, edge_list; edge_list)
end

# Creates a LatticeReactionSystem from a DiGraph lattice.
function LatticeReactionSystem(rs, srs, lattice::SimpleGraph; directed_edges = true)
    num_verts = nv(lattice)
    edge_list = false
    return LatticeReactionSystem(rs, srs, lattice, num_verts, edge_list; directed_edges, edge_list)
end

# Creates a LatticeReactionSystem from a Graph lattice.
function LatticeReactionSystem(rs, srs, lattice::SimpleGraph)
    return LatticeReactionSystem(rs, srs, DiGraph(lattice); directed_edges)
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
has_cartesian_lattice(lrs::LatticeReactionSystem) = lrs.lattice isa CartesianGridRej{S,T} where {S,T}
# Checks if a LatticeReactionSystem have a regular grid lattice.
has_regular_lattice(lrs::LatticeReactionSystem) = lrs.lattice isa Array{Bool, T} where T
# Checks if a LatticeReactionSystem have a grid lattice (cartesian or regular).
has_grid_lattice(lrs::LatticeReactionSystem) = (has_cartesian_lattice(lrs) || has_regular_lattice(lrs))
# Checks if a LatticeReactionSystem have a graph lattice.
has_graph_lattice(lrs::LatticeReactionSystem) = lrs.lattice isa SimpleDiGraph

# Returns the dimensions of a LatticeReactionNetwork with a grid lattice.
function grid_dims(lrs::LatticeReactionSystem)
    has_cartesian_lattice(lrs) && (return length(lrs.lattice.dims))
    has_regular_lattice(lrs) && (return length(size(lrs.lattice)))
    error("Dimensions are only defined for LatticeReactionSystems with grid-based lattices (not graph-based).")
end