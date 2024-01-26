### Lattice Reaction Network Structure ###

# Describes a spatial reaction network over a lattice.
# Adding the "<: MT.AbstractTimeDependentSystem" part messes up show, disabling me from creating LRSs. Should be fixed some time.
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

    """
    Whenever the initial input was directed. True for Digraph input, false for other graphs and all grids.
    Used later to know whenever duplication of edge parameter should be duplicated
    (i.e. allow parameter for edge i => j to be used for j => i).
    Even if false, giving separate parameters for both directions is still permitted. 
    """
    directed_edges::Bool
    """
    An iterator over all the edges on the lattice. 
    The format depends on the type of lattice (Cartesian grid, grid, or graph).
    """
    edge_iterator::T

    function LatticeReactionSystem(rs::ReactionSystem{Q}, spatial_reactions::Vector{R}, lattice::S, 
                                   num_verts, num_edges, edge_iterator::T; directed_edges = false) where {Q,R, S, T}
        # Error checks.
        if !(R <: AbstractSpatialReaction)
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

        # Counts the number of species types. The number of vertexes and compartments is given in input.
        # `num_edges` cannot be computed here, because how it is computed depends on the lattice and `typeof(edge_iterator)`.
        num_species = length(unique([species(rs); spat_species]))

        # Ensures the parameter order begins similarly to in the non-spatial ReactionSystem.
        ps = [parameters(rs); setdiff([edge_parameters; vertex_parameters], parameters(rs))]    

        # Checks that all spatial reactions are valid for this reactions system.
        foreach(sr -> check_spatial_reaction_validity(rs, sr; edge_parameters=edge_parameters), spatial_reactions)   

        return new{Q,R,S,T}(rs, spatial_reactions, lattice, num_verts, num_edges, num_species, 
                            spat_species, ps, vertex_parameters, edge_parameters, directed_edges, edge_iterator)
    end
end

# Creates a LatticeReactionSystem from a CartesianGrid lattice (cartesian grid).
function LatticeReactionSystem(rs, srs, lattice::CartesianGridRej{S,T}; diagonal_connections=false) where {S,T}
    # Error checks.
    (length(lattice.dims) > 3) && error("Grids of higher dimension than 3 is currently not supported.")

    # Ensures that the matrix is on a 3d form
    lattice = CartesianGrid((lattice.dims..., fill(1, 3-length(lattice.dims))...))

    # Counts vertexes and edges. The `num_edges` count formula counts the number of internal, side, 
    # edge, and corner vertexes (on the grid). The number of edges from each depend on whether diagonal
    # connections are allowed. The formula holds even if l, m, and/or n are 1.
    l,m,n = 2,3,4
    num_verts = l * m * n
    (ni, ns, ne, nc) = diagonal_connections ? (26,17,11,7) : (6,5,4,3)
    num_edges = ni*(l-2)*(m-2)*(n-2) +                            # Internal vertexes.
                ns*(2(l-2)*(m-2) + 2(l-2)*(n-2) + 2(m-2)*(n-2)) + # Side vertexes.
                ne*(4(l-2) + 4(m-2) + 4(n-2)) +                   # Edge vertexes.
                nc*8                                             # Corner vertexes.
    
    # Creates an iterator over all edges.
    # Currently creates a full vector. Future version might be for efficient.
    edge_iterator = Vector{Pair{Int64,Int64}}(undef, num_edges)
    # Loops through, simultaneously, the coordinates of each position in the grid, as well as that 
    # coordinate's (scalar) flat index.  For each grid point, loops through all potential neighbours.  
    indices = [(L, M, N) for L in 1:l, M in 1:m, N in 1:n]
    flat_indices = 1:num_verts
    next_vert = 0
    for ((L, M, N), idx) in zip(indices, flat_indices)
        for LL in max(L - 1, 1):min(L + 1, l), 
            MM in max(M - 1, 1):min(M + 1, m), 
            NN in max(N - 1, 1):min(N + 1, n)
    
            # Which (LL,MM,NN) indexes are valid neighbours depends on whether diagonal connects are permitted.
            !diagonal_connections && (count([L==LL, M==MM, N==NN]) == 2) || continue
            diagonal_connections && (L==LL) && (M==MM) && (N==NN) && continue
            
            # Computes the neighbours scalar index. Add that connection to `edge_iterator`.
            neighbour_idx = LL + (MM - 1) * l + (NN - 1) * m * l
            edge_iterator[next_vert += 1] = (idx => neighbour_idx)
        end
    end

    return LatticeReactionSystem(rs, srs, lattice, num_verts, num_edges, edge_iterator)
end

# Creates a LatticeReactionSystem from a Boolean Array lattice (masked grid).
function LatticeReactionSystem(rs, srs, lattice::Array{Bool, T}; diagonal_connections=false) where {T}  
    # Error checks.
    dims = size(lattice)
    (length(dims) > 3) && error("Grids of higher dimension than 3 is currently not supported.")

    # Ensures that the matrix is on a 3d form
    lattice = reshape(lattice, [dims...; fill(1, 3-length(dims))]...)
    
    # Counts vertexes (edges have to be counted after the iterator have been created).
    num_verts = count(lattice)

    # Creates an iterator over all edges.Currently a full vector of all edges (as pairs).
    edge_iterator = Vector{Pair{Int64,Int64}}()
    # Loops through, simultaneously, the coordinates of each position in the grid, as well as that 
    # coordinate's (scalar) flat index.  For each grid point, loops through all potential neighbours.  
    l,m,n = size(lattice)
    indices = [(L, M, N) for L in 1:l, M in 1:m, N in 1:n]
    flat_indices = 1:num_verts
    for ((L, M, N), idx) in zip(indices, flat_indices)
        for LL in max(L - 1, 1):min(L + 1, l), 
            MM in max(M - 1, 1):min(M + 1, m), 
            NN in max(N - 1, 1):min(N + 1, n)

            # Ensures that the neighbour is a valid lattice point.
            lattice[LL,MM,NN] || continue

            # Which (LL,MM,NN) indexes are valid neighbours depends on whether diagonal connects are permitted.
            !diagonal_connections && (count([L==LL, M==MM, N==NN]) == 2) || continue
            diagonal_connections && (L==LL) && (M==MM) && (N==NN) && continue
            
            # Computes the neighbours scalar index. Add that connection to `edge_iterator`.
            neighbour_idx = LL + (MM - 1) * l + (NN - 1) * m * l
            push!(edge_iterator, idx => neighbour_idx)
        end
    end
    num_edges = length(edge_iterator)

    return LatticeReactionSystem(rs, srs, lattice, num_verts, num_edges, edge_iterator)
end

# Creates a LatticeReactionSystem from a DiGraph lattice (graph grid).
function LatticeReactionSystem(rs, srs, lattice::DiGraph; directed_edges = true)
    num_verts = nv(lattice)
    num_edges = ne(lattice)
    edge_iterator = edges(lattice)
    return LatticeReactionSystem(rs, srs, lattice, num_verts, num_edges, edge_iterator; directed_edges)
end

# Creates a LatticeReactionSystem from a Graph lattice (graph grid).
function LatticeReactionSystem(rs, srs, lattice::SimpleGraph)
    return LatticeReactionSystem(rs, srs, DiGraph(lattice); directed_edges = false)
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
# Checks if a LatticeReactionSystem have a masked grid lattice.
has_masked_lattice(lrs::LatticeReactionSystem) = lrs.lattice isa Array{Bool, T} where T
# Checks if a LatticeReactionSystem have a grid lattice (cartesian or masked).
has_grid_lattice(lrs::LatticeReactionSystem) = (has_cartesian_lattice(lrs) || has_masked_lattice(lrs))
# Checks if a LatticeReactionSystem have a graph lattice.
has_graph_lattice(lrs::LatticeReactionSystem) = lrs.lattice isa SimpleDiGraph

# Returns the dimensions of a LatticeReactionNetwork with a grid lattice.
function grid_dims(lrs::LatticeReactionSystem)
    has_cartesian_lattice(lrs) && (return length(lrs.lattice.dims))
    has_masked_lattice(lrs) && (return length(size(lrs.lattice)))
    error("Dimensions are only defined for LatticeReactionSystems with grid-based lattices (not graph-based).")
end