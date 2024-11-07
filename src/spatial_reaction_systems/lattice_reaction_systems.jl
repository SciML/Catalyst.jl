### New Type Unions ###

# Cartesian and masked grids share several traits, hence we declare a common (union) type for them.
const GridLattice{N, T} = Union{Array{Bool, N}, CartesianGridRej{N, T}}

### Lattice Reaction Network Structure ###

"""
$(TYPEDEF)

A representation of a spatial system of chemical reactions on a discrete (lattice) space. 

# Fields
$(FIELDS)

Arguments:
- `rs`: The non-spatial [`ReactionSystem`](@ref) model that is expanded to a spatial model.
- `srs`: A vector of spatial reactions. These provide the rules for how species may move spatially.
- `lattice`: Either a Cartesian grid, a masked grid, or a graph. This describes the discrete space
to which the non-spatial model is expanded.

Keyword Arguments:
- `diagonal_connections = false`: Only relevant for Cartesian and masked lattices. If `true`, 
diagonally adjacent compartments are considered adjacent, and spatial reactions in between these
are possible.

Example:
```julia
# Fetch packages.
using Catalyst, OrdinaryDiffEq
import CairoMakie

# Creates the `LatticeReactionSystem` model.
rs = @reaction_network begin
    (p,d), 0 <--> X
end
diffusion_rx = @transport_reaction D X
lattice = CartesianGrid((5,5))
lrs = LatticeReactionSystem(rs, [diffusion_rx], lattice)

# Simulates the model (using ODE and jumps).
u0 = [:X => rand(5,5)]
tspan = (0.0, 1.0)
ps = [:p => 1.0, :d => 0.5, :D => 0.1]
oprob = ODEProblem(lrs, u0, tspan, ps)
osol = solve(oprob)

# Saves an animation of the solution to the file "lattice_animation.mp4".
lattice_animation(osol, :X, lrs, "lattice_animation.mp4")
```

Notes:
- Spatial modelling in Catalyst is still a work in progress, any feedback (or contributions) to this
is highly welcome.
- `LatticeReactionSystem`s are primarily intended to model systems in discrete space. Modelling
continuous space systems with them is possible, but requires the user to determine the discretisation
(the lattice). Better support for continuous space models is a work in progress.
- Catalyst contains extensive documentation on spatial modelling, which can be found [here](https://docs.sciml.ai/Catalyst/stable/spatial_modelling/lattice_reaction_systems/).
"""
struct LatticeReactionSystem{Q, R, S, T} <: MT.AbstractTimeDependentSystem
    # Input values.
    """The (non-spatial) reaction system within each vertex."""
    reactionsystem::ReactionSystem{Q}
    """The spatial reactions defined between individual vertices."""
    spatial_reactions::Vector{R}
    """The lattice on which the (discrete) spatial system is defined."""
    lattice::S

    # Derived values.
    """The number of vertices (compartments)."""
    num_verts::Int64
    """The number of edges."""
    num_edges::Int64
    """The number of species."""
    num_species::Int64

    """List of species that may move spatially."""
    spatial_species::Vector{BasicSymbolic{Real}}
    """
    All parameters related to the lattice reaction system
    (both those whose values are tied to vertices and edges).
    """
    parameters::Vector{Any}
    """
    Parameters which values are tied to vertices, 
    e.g. that possibly could have unique values at each vertex of the system.
    """
    vertex_parameters::Vector{Any}
    """
    Parameters whose values are tied to edges (adjacencies), 
    e.g. that possibly could have unique values at each edge of the system.
    """
    edge_parameters::Vector{Any}
    """
    An iterator over all the lattice's edges. Currently, the format is always a Vector{Pair{Int64,Int64}}.
    However, in the future, different types could potentially be used for different types of lattice
    (E.g. for a Cartesian grid, we do not technically need to enumerate each edge)
    """
    edge_iterator::T

    function LatticeReactionSystem(rs::ReactionSystem{Q}, spatial_reactions::Vector{R},
            lattice::S, num_verts::Int64, num_edges::Int64,
            edge_iterator::T) where {Q, R, S, T}
        # Error checks.
        if !(R <: AbstractSpatialReaction)
            throw(ArgumentError("The second argument must be a vector of AbstractSpatialReaction subtypes."))
        end
        if !iscomplete(rs)
            throw(ArgumentError("A non-complete `ReactionSystem` was used as input, this is not permitted."))
        end
        if !isempty(MT.get_systems(rs))
            throw(ArgumentError("A non-flattened (hierarchical) `ReactionSystem` was used as input. `LatticeReactionSystem`s can only be based on non-hierarchical `ReactionSystem`s."))
        end
        if length(reactions(rs)) != length(equations(rs))
            throw(ArgumentError("The `ReactionSystem` used as input contain equations (in addition to reactions). This is not permitted."))
        end
        if length(species(rs)) != length(unknowns(rs))
            throw(ArgumentError("The `ReactionSystem` used as input contain variable unknowns (in addition to species unknowns). This is not permitted (the input `ReactionSystem` must contain species unknowns only)."))
        end
        if !isempty(MT.continuous_events(rs)) || !isempty(MT.discrete_events(rs))
            throw(ArgumentError("The `ReactionSystem` used as input to `LatticeReactionSystem contain events. These will be ignored in any simulations based on the created `LatticeReactionSystem`."))
        end
        if !isempty(observed(rs))
            @warn "The `ReactionSystem` used as input to `LatticeReactionSystem contain observables. It will not be possible to access these from the created `LatticeReactionSystem`."
        end

        # Computes the species which are parts of spatial reactions. Also counts the total number of 
        # species types.
        if isempty(spatial_reactions)
            spat_species = Vector{BasicSymbolic{Real}}[]
        else
            spat_species = unique(reduce(vcat,
                [spatial_species(sr) for sr in spatial_reactions]))
        end
        num_species = length(unique([species(rs); spat_species]))

        # Computes the sets of vertex, edge, and all, parameters.
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

        # Checks that all spatial reactions are valid for this reaction system.
        foreach(
            sr -> check_spatial_reaction_validity(rs, sr; edge_parameters = edge_parameters),
            spatial_reactions)

        return new{Q, R, S, T}(
            rs, spatial_reactions, lattice, num_verts, num_edges, num_species,
            spat_species, ps, vertex_parameters, edge_parameters, edge_iterator)
    end
end

# Creates a LatticeReactionSystem from a (directed) Graph lattice (graph grid).
function LatticeReactionSystem(rs, srs, lattice::DiGraph)
    num_verts = nv(lattice)
    num_edges = ne(lattice)
    edge_iterator = [e.src => e.dst for e in edges(lattice)]
    return LatticeReactionSystem(rs, srs, lattice, num_verts, num_edges, edge_iterator)
end
# Creates a LatticeReactionSystem from a (undirected) Graph lattice (graph grid).
function LatticeReactionSystem(rs, srs, lattice::SimpleGraph)
    LatticeReactionSystem(rs, srs, DiGraph(lattice))
end

# Creates a LatticeReactionSystem from a CartesianGrid lattice (cartesian grid) or a Boolean Array 
# lattice (masked grid). These two are quite similar, so much code can be reused in a single interface.
function LatticeReactionSystem(rs, srs, lattice::GridLattice{N, T};
        diagonal_connections = false) where {N, T}
    # Error checks.
    (N > 3) && error("Grids of higher dimension than 3 is currently not supported.")

    # Computes the number of vertices and edges. The two grid types (Cartesian and masked) each
    # uses their own function for this.
    num_verts = count_verts(lattice)
    num_edges = count_edges(lattice; diagonal_connections)

    # Finds all the grid's edges. First computer `flat_to_grid_idx` which is a vector which takes
    # each vertex's flat (scalar) index to its grid index (e.g. (3,5) for a 2d grid). Next compute
    # `grid_to_flat_idx` which is an array (of the same size as the grid) that does the reverse conversion.
    # Especially for masked grids these have non-trivial forms.
    flat_to_grid_idx, grid_to_flat_idx = get_index_converters(lattice, num_verts)
    # Simultaneously iterates through all vertices' flat and grid indices. Finds the (grid) indices
    # of their neighbours (different approaches for the two grid types). Converts these to flat
    # indices and adds the edges to `edge_iterator`.
    cur_vert = 0
    g_size = grid_size(lattice)
    edge_iterator = Vector{Pair{Int64, Int64}}(undef, num_edges)
    for (flat_idx, grid_idx) in enumerate(flat_to_grid_idx)
        for neighbour_grid_idx in get_neighbours(lattice, grid_idx, g_size;
            diagonal_connections)
            cur_vert += 1
            edge_iterator[cur_vert] = flat_idx => grid_to_flat_idx[neighbour_grid_idx...]
        end
    end

    return LatticeReactionSystem(rs, srs, lattice, num_verts, num_edges, edge_iterator)
end

### LatticeReactionSystem Helper Functions ###
# Note, most of these are specifically for (Cartesian or masked) grids, we call them `grid`, not `lattice`.

# Counts the number of vertices on a (Cartesian or masked) grid.
count_verts(grid::CartesianGridRej{N, T}) where {N, T} = prod(grid_size(grid))
count_verts(grid::Array{Bool, N}) where {N} = count(grid)

# Counts and edges on a Cartesian grid. The formula counts the number of internal, side, edge, and
# corner vertices (on the grid). `l,m,n = grid_dims(grid),1,1` ensures that "extra" dimensions get
# length 1. The formula holds even if one or more of l, m, and n are 1.
function count_edges(grid::CartesianGridRej{N, T};
        diagonal_connections = false) where {N, T}
    l, m, n = grid_size(grid)..., 1, 1
    (ni, ns, ne, nc) = diagonal_connections ? (26, 17, 11, 7) : (6, 5, 4, 3)
    num_edges = ni * (l - 2) * (m - 2) * (n - 2) +                            # Edges from internal vertices.
                ns * (2(l - 2) * (m - 2) + 2(l - 2) * (n - 2) + 2(m - 2) * (n - 2)) + # Edges from side vertices.
                ne * (4(l - 2) + 4(m - 2) + 4(n - 2)) +                   # Edges from edge vertices.
                nc * 8                                              # Edges from corner vertices.
    return num_edges
end

# Counts and edges on a masked grid. Does so by looping through all the vertices of the grid, 
# finding their neighbours, and updating the edge count accordingly.
function count_edges(grid::Array{Bool, N}; diagonal_connections = false) where {N}
    g_size = grid_size(grid)
    num_edges = 0
    for grid_idx in get_grid_indices(grid)
        grid[grid_idx] || continue
        num_edges += length(get_neighbours(grid, Tuple(grid_idx), g_size;
            diagonal_connections))
    end
    return num_edges
end

# For a (1d, 2d, or 3d) (Cartesian or masked) grid, returns a vector and an array, permitting the 
# conversion between a vertex's flat (scalar) and grid indices. E.g. for a 2d grid, if grid point (3,2)
# corresponds to the fifth vertex, then `flat_to_grid_idx[5] = (3,2)` and `grid_to_flat_idx[3,2] = 5`.
function get_index_converters(grid::GridLattice{N, T}, num_verts) where {N, T}
    flat_to_grid_idx = Vector{typeof(grid_size(grid))}(undef, num_verts)
    grid_to_flat_idx = Array{Int64}(undef, grid_size(grid))

    # Loops through the flat and grid indices simultaneously, adding them to their respective converters.
    cur_flat_idx = 0
    for grid_idx in get_grid_indices(grid)
        # For a masked grid, grid points with `false` values are skipped.
        (grid isa Array{Bool}) && (!grid[grid_idx]) && continue

        cur_flat_idx += 1
        flat_to_grid_idx[cur_flat_idx] = grid_idx
        grid_to_flat_idx[grid_idx] = cur_flat_idx
    end
    return flat_to_grid_idx, grid_to_flat_idx
end

# For a vertex's grid index, and a lattice, returns the grid indices of all its (valid) neighbours.
function get_neighbours(grid::GridLattice{N, T}, grid_idx, g_size;
        diagonal_connections = false) where {N, T}
    # Depending on the grid's dimension, find all potential neighbours.
    if grid_dims(grid) == 1
        potential_neighbours = [grid_idx .+ (i) for i in -1:1]
    elseif grid_dims(grid) == 2
        potential_neighbours = [grid_idx .+ (i, j) for i in -1:1 for j in -1:1]
    else
        potential_neighbours = [grid_idx .+ (i, j, k) for i in -1:1 for j in -1:1
                                for k in -1:1]
    end

    # Depending on whether diagonal connections are used or not, find valid neighbours.
    if diagonal_connections
        filter!(n_idx -> n_idx !== grid_idx, potential_neighbours)
    else
        filter!(n_idx -> count(n_idx .== grid_idx) == (length(g_size) - 1),
            potential_neighbours)
    end

    # Removes neighbours outside of the grid, and returns the full list.
    return filter(n_idx -> is_valid_grid_point(grid, n_idx, g_size), potential_neighbours)
end

# Checks if a grid index corresponds to a valid grid point. First, check that each dimension of the
# index is within the grid's bounds. Next, perform an extra check for the masked grid.
function is_valid_grid_point(grid::GridLattice{N, T}, grid_idx, g_size) where {N, T}
    if !all(0 < g_idx <= dim_leng for (g_idx, dim_leng) in zip(grid_idx, g_size))
        return false
    end
    return (grid isa Array{Bool}) ? grid[grid_idx...] : true
end

# Gets an iterator over a grid's grid indices. Separate function so we can handle the two grid types
# separately (i.e. not calling `CartesianIndices(ones(grid_size(grid)))` unnecessarily for masked grids).
function get_grid_indices(grid::CartesianGridRej{N, T}) where {N, T}
    CartesianIndices(ones(grid_size(grid)))
end
get_grid_indices(grid::Array{Bool, N}) where {N} = CartesianIndices(grid)

### LatticeReactionSystem-specific Getters ###

# Basic getters (because `LatticeReactionSystem`s are `AbstractSystem`s), normal `lrs.field` does not
# work and these getters must be used throughout all code.
"""
    reactionsystem(lrs::LatticeReactionSystem)

Returns the non-spatial `ReactionSystem` stored in a `LatticeReactionSystem`.
"""
reactionsystem(lrs::LatticeReactionSystem) = getfield(lrs, :reactionsystem)

"""
    spatial_reactions(lrs::LatticeReactionSystem)

Returns a vector with all the spatial reactions stored in a `LatticeReactionSystem`.
"""
spatial_reactions(lrs::LatticeReactionSystem) = getfield(lrs, :spatial_reactions)

"""
    lattice(lrs::LatticeReactionSystem)

Returns the lattice stored in a `LatticeReactionSystem`.
"""
lattice(lrs::LatticeReactionSystem) = getfield(lrs, :lattice)

"""
    num_verts(lrs::LatticeReactionSystem)

Returns the number of vertices (i.e. compartments) in the lattice stored in a `LatticeReactionSystem`.
"""
num_verts(lrs::LatticeReactionSystem) = getfield(lrs, :num_verts)

"""
    num_edges(lrs::LatticeReactionSystem)

Returns the number of edges (i.e. connections between vertices) in the lattice stored in a 
`LatticeReactionSystem`.
"""
num_edges(lrs::LatticeReactionSystem) = getfield(lrs, :num_edges)

"""
    num_species(lrs::LatticeReactionSystem)

Returns the number of species that a `LatticeReactionSystem` contains.
"""
num_species(lrs::LatticeReactionSystem) = getfield(lrs, :num_species)

"""
    spatial_species(lrs::LatticeReactionSystem)

Returns the number of species that can move spatially that a `LatticeReactionSystem` contains.
"""
spatial_species(lrs::LatticeReactionSystem) = getfield(lrs, :spatial_species)

# Returns the parameters in a `LatticeReactionSystem`
MT.parameters(lrs::LatticeReactionSystem) = getfield(lrs, :parameters)

"""
    vertex_parameters(lrs::LatticeReactionSystem)

Returns all the parameters of a `LatticeReactionSystem` whose values are tied to vertices.
"""
vertex_parameters(lrs::LatticeReactionSystem) = getfield(lrs, :vertex_parameters)

"""
    edge_parameters(lrs::LatticeReactionSystem)

Returns all the parameters of a `LatticeReactionSystem` whose values are tied to edges.
"""
edge_parameters(lrs::LatticeReactionSystem) = getfield(lrs, :edge_parameters)

"""
    edge_iterator(lrs::LatticeReactionSystem)

Returns an iterator over all of the edges in the lattice stored in a `LatticeReactionSystem`. Each
edge is a `Pair{Int64, Int64}`, taking the source vertex to the destination vertex.
"""
edge_iterator(lrs::LatticeReactionSystem) = getfield(lrs, :edge_iterator)

# Non-trivial getters.
"""
    is_transport_system(lrs::LatticeReactionSystem)

Returns `true` if all spatial reactions in `lrs` are `TransportReaction`s.
"""
function is_transport_system(lrs::LatticeReactionSystem)
    return all(sr -> sr isa TransportReaction, spatial_reactions(lrs))
end

"""
    has_cartesian_lattice(lrs::LatticeReactionSystem)

Returns `true` if `lrs` was created using a cartesian grid lattice (e.g. created via `CartesianGrid(5,5)`). 
Otherwise, returns `false`.
"""
has_cartesian_lattice(lrs::LatticeReactionSystem) = lattice(lrs) isa
                                                    CartesianGridRej{N, T} where {N, T}

"""
    has_masked_lattice(lrs::LatticeReactionSystem)

Returns `true` if `lrs` was created using a masked grid lattice (e.g. created via `[true true; true false]`). 
Otherwise, returns `false`.
"""
has_masked_lattice(lrs::LatticeReactionSystem) = lattice(lrs) isa Array{Bool, N} where {N}

"""
    has_grid_lattice(lrs::LatticeReactionSystem)

Returns `true` if `lrs` was created using a cartesian or masked grid lattice. Otherwise, returns `false`.
"""
function has_grid_lattice(lrs::LatticeReactionSystem)
    return has_cartesian_lattice(lrs) || has_masked_lattice(lrs)
end

"""
    has_graph_lattice(lrs::LatticeReactionSystem)

Returns `true` if `lrs` was created using a graph grid lattice (e.g. created via `path_graph(5)`). 
Otherwise, returns `false`.
"""
has_graph_lattice(lrs::LatticeReactionSystem) = lattice(lrs) isa SimpleDiGraph

"""
    grid_size(lrs::LatticeReactionSystem)

Returns the size of `lrs`'s lattice (only if it is a cartesian or masked grid lattice). 
E.g. for a lattice `CartesianGrid(4,6)`, `(4,6)` is returned.
"""
grid_size(lrs::LatticeReactionSystem) = grid_size(lattice(lrs))
grid_size(lattice::CartesianGridRej{N, T}) where {N, T} = lattice.dims
grid_size(lattice::Array{Bool, N}) where {N} = size(lattice)
function grid_size(lattice::Graphs.AbstractGraph)
    throw(ArgumentError("Grid size is only defined for LatticeReactionSystems with grid-based lattices (not graph-based)."))
end

"""
    grid_dims(lrs::LatticeReactionSystem)

Returns the number of dimensions of `lrs`'s lattice (only if it is a cartesian or masked grid lattice). 
The output is either `1`, `2`, or `3`.
"""
grid_dims(lrs::LatticeReactionSystem) = grid_dims(lattice(lrs))
grid_dims(lattice::GridLattice{N, T}) where {N, T} = return N
function grid_dims(lattice::Graphs.AbstractGraph)
    throw(ArgumentError("Grid dimensions is only defined for LatticeReactionSystems with grid-based lattices (not graph-based)."))
end

"""
    get_lattice_graph(lrs::LatticeReactionSystem)

Returns lrs's lattice, but in as a graph. Currently does not work for Cartesian lattices.
"""
function get_lattice_graph(lrs::LatticeReactionSystem)
    has_graph_lattice(lrs) && return lattice(lrs)
    return Graphs.SimpleGraphFromIterator(Graphs.SimpleEdge(e[1], e[2])
    for e in edge_iterator(lrs))
end

### Catalyst-based Getters ###

# Get all species.
function species(lrs::LatticeReactionSystem)
    unique([species(reactionsystem(lrs)); spatial_species(lrs)])
end

# Generic ones (simply forwards call to the non-spatial system).
reactions(lrs::LatticeReactionSystem) = reactions(reactionsystem(lrs))

### ModelingToolkit-based Getters ###

# Generic ones (simply forwards call to the non-spatial system)
# The `parameters` MTK getter have a specialised accessor for LatticeReactionSystems.
MT.nameof(lrs::LatticeReactionSystem) = MT.nameof(reactionsystem(lrs))
MT.get_iv(lrs::LatticeReactionSystem) = MT.get_iv(reactionsystem(lrs))
MT.equations(lrs::LatticeReactionSystem) = MT.equations(reactionsystem(lrs))
MT.unknowns(lrs::LatticeReactionSystem) = MT.unknowns(reactionsystem(lrs))
MT.get_metadata(lrs::LatticeReactionSystem) = MT.get_metadata(reactionsystem(lrs))

# Lattice reaction systems should not be combined with compositional modelling.
# Maybe these should be allowed anyway? Still feel a bit weird
function MT.get_eqs(lrs::LatticeReactionSystem)
    MT.get_eqs(reactionsystem(lrs))
end
function MT.get_unknowns(lrs::LatticeReactionSystem)
    MT.get_unknowns(reactionsystem(lrs))
end
function MT.get_ps(lrs::LatticeReactionSystem)
    MT.get_ps(reactionsystem(lrs))
end

# Technically should not be used, but has to be declared for the `show` function to work.
function MT.get_systems(lrs::LatticeReactionSystem)
    return []
end

# Other non-relevant getters.
function MT.independent_variables(lrs::LatticeReactionSystem)
    MT.independent_variables(reactionsystem(lrs))
end

### Edge Parameter Value Generators ###

"""
    make_edge_p_values(lrs::LatticeReactionSystem, make_edge_p_value::Function)

Generates edge parameter values for a lattice reaction system. Only work for (Cartesian or masked) 
grid lattices (without diagonal adjacencies). 

Input:
- `lrs`: The lattice reaction system for which values should be generated.
    - `make_edge_p_value`: a function describing a rule for generating the edge parameter values.

Output:
    - `ep_vals`: A sparse matrix of size (num_verts,num_verts) (where num_verts is the number of 
    vertices in `lrs`). Here, `eps[i,j]` is filled only if there is an edge going from vertex i to 
    vertex j. The value of `eps[i,j]` is determined by `make_edge_p_value`.

Here, `make_edge_p_value` should take two arguments, `src_vert` and `dst_vert`, which correspond to 
the grid indices of an edge's source and destination vertices, respectively. It outputs a single value,
which is the value assigned to that edge.

Example:
    In the following example, we assign the value `0.1` to all edges, except for the one leading from
    vertex (1,1) to vertex (1,2), to which we assign the value `1.0`.
```julia
using Catalyst
rn = @reaction_network begin
    (p,d), 0 <--> X
end
tr = @transport_reaction D X
lattice = CartesianGrid((5,5))
lrs = LatticeReactionSystem(rn, [tr], lattice)

function make_edge_p_value(src_vert, dst_vert)
    if src_vert == (1,1) && dst_vert == (1,2)
        return 1.0
    else
        return 0.1
    end
end

D_vals = make_edge_p_values(lrs, make_edge_p_value)
```
"""
function make_edge_p_values(lrs::LatticeReactionSystem, make_edge_p_value::Function)
    if has_graph_lattice(lrs)
        error("The `make_edge_p_values` function is only meant for lattices with (Cartesian or masked) grid structures. It cannot be applied to graph lattices.")
    end

    # Makes the flat to index grid converts. Predeclared the edge parameter value sparse matrix.
    flat_to_grid_idx = get_index_converters(lattice(lrs), num_verts(lrs))[1]
    values = spzeros(num_verts(lrs), num_verts(lrs))

    # Loops through all edges, and applies the value function to these.
    for e in edge_iterator(lrs)
        # This extra step is needed to ensure that `0` is stored if make_edge_p_value yields a 0.
        # If not, then the sparse matrix simply becomes empty in that position.
        values[e[1], e[2]] = eps()

        values[e[1], e[2]] = make_edge_p_value(flat_to_grid_idx[e[1]],
            flat_to_grid_idx[e[2]])
    end

    return values
end

"""
    make_directed_edge_values(lrs::LatticeReactionSystem, x_vals::Tuple{T,T}, y_vals::Tuple{T,T} = (undef,undef), 
                     z_vals::Tuple{T,T} = (undef,undef)) where {T}

Generates edge parameter values for a lattice reaction system. Only work for (Cartesian or masked) 
grid lattices (without diagonal adjacencies). Each dimension (x, and possibly y and z), and 
direction has assigned its own constant edge parameter value.

Input:
    - `lrs`: The lattice reaction system for which values should be generated.
    - `x_vals::Tuple{T,T}`: The values in the increasing (from a lower x index to a higher x index) 
    and decreasing (from a higher x index to a lower x index) direction along the x dimension.
    - `y_vals::Tuple{T,T}`: The values in the increasing and decreasing direction along the y dimension. 
    Should only be used for 2 and 3-dimensional grids.
    - `z_vals::Tuple{T,T}`: The values in the increasing and decreasing direction along the z dimension. 
    Should only be used for 3-dimensional grids.

Output:
    - `ep_vals`: A sparse matrix of size (num_verts,num_verts) (where num_verts is the number of 
    vertices in `lrs`). Here, `eps[i,j]` is filled only if there is an edge going from vertex i to 
    vertex j. The value of `eps[i,j]` is determined by the `x_vals`, `y_vals`, and `z_vals` Tuples, 
    and vertices i and j's relative position in the grid.

It should be noted that two adjacent vertices will always be different in exactly a single dimension 
(x, y, or z). The corresponding tuple determines which value is assigned.

Example:
    In the following example, we wish to have diffusion in the x dimension, but a constant flow from
    low y values to high y values (so not transportation from high to low y). We achieve it in the 
    following manner:
```julia
using Catalyst
rn = @reaction_network begin
    (p,d), 0 <--> X
end
tr = @transport_reaction D X
lattice = CartesianGrid((5,5))
lrs = LatticeReactionSystem(rn, [tr], lattice)

D_vals = make_directed_edge_values(lrs, (0.1, 0.1), (0.1, 0.0))
```
Here, since we have a 2d grid, we only provide the first two Tuples to `make_directed_edge_values`.
"""
function make_directed_edge_values(lrs::LatticeReactionSystem, x_vals::Tuple{T, T},
        y_vals::Union{Nothing, Tuple{T, T}} = nothing,
        z_vals::Union{Nothing, Tuple{T, T}} = nothing) where {T}
    # Error checks.
    if has_graph_lattice(lrs)
        error("The `make_directed_edge_values` function is only meant for lattices with (Cartesian or masked) grid structures. It cannot be applied to graph lattices.")
    end
    if count(!isnothing(flow) for flow in [x_vals, y_vals, z_vals]) != grid_dims(lrs)
        error("You must provide flows in the same number of dimensions as your lattice has dimensions. The lattice have $(grid_dims(lrs)), and flows where provided for $(count(isnothing(flow) for flow in [x_vals, y_vals, z_vals])) dimensions.")
    end

    # Declares a function that assigns the correct flow value for a given edge.
    function directed_vals_func(src_vert, dst_vert)
        if count(src_vert .== dst_vert) != (grid_dims(lrs) - 1)
            error("The `make_directed_edge_values` function can only be applied to lattices with rigid (non-diagonal) grid structure. It is being evaluated for the edge from $(src_vert) to $(dst_vert), which does not seem directly adjacent on a grid.")
        elseif src_vert[1] != dst_vert[1]
            return (src_vert[1] < dst_vert[1]) ? x_vals[1] : x_vals[2]
        elseif src_vert[2] != dst_vert[2]
            return (src_vert[2] < dst_vert[2]) ? y_vals[1] : y_vals[2]
        elseif src_vert[3] != dst_vert[3]
            return (src_vert[3] < dst_vert[3]) ? z_vals[1] : z_vals[2]
        else
            error("Problem when evaluating adjacency type for the edge from $(src_vert) to $(dst_vert).")
        end
    end

    # Uses the make_edge_p_values function to compute the output.
    return make_edge_p_values(lrs, directed_vals_func)
end
