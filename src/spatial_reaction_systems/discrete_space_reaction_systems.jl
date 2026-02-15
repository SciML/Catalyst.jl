### New Type Unions ###

# Cartesian and masked grids share several traits, hence we declare a common (union) type for them.
const GridLattice{N, T} = Union{Array{Bool, N}, CartesianGridRej{N, T}}

### Lattice Reaction Network Structure ###

"""
$(TYPEDEF)

A representation of a spatial system of chemical reactions on a discrete space.

# Fields
$(FIELDS)

Arguments:
- `rs`: The non-spatial [`ReactionSystem`](@ref) model that is expanded to a spatial model.
- `srs`: A vector of spatial reactions. These provide the rules for how species may move spatially.
- `dspace`: Either a Cartesian grid, a masked grid, or a graph. This describes the discrete space
to which the non-spatial model is expanded.

Keyword Arguments:
- `diagonal_connections = false`: Only relevant for Cartesian and masked spaces. If `true`,
diagonally adjacent compartments are considered adjacent, and spatial reactions in between these
are possible.

Example:
```julia
# Fetch packages.
using Catalyst, OrdinaryDiffEqDefault
import CairoMakie

# Creates the `DiscreteSpaceReactionSystem` model.
rs = @reaction_network begin
    (p,d), 0 <--> X
end
diffusion_rx = @transport_reaction D X
space = CartesianGrid((5,5))
dsrs = DiscreteSpaceReactionSystem(rs, [diffusion_rx], space)

# Simulates the model (using ODE and jumps).
u0 = [:X => rand(5,5)]
tspan = (0.0, 1.0)
ps = [:p => 1.0, :d => 0.5, :D => 0.1]
oprob = ODEProblem(dsrs, u0, tspan, ps)
osol = solve(oprob)

# Saves an animation of the solution to the file "dspace_animation.mp4".
dspace_animation(osol, :X, dsrs, "dspace_animation.mp4")
```

Notes:
- Spatial modelling in Catalyst is still a work in progress, any feedback (or contributions) to this
is highly welcome.
- `DiscreteSpaceReactionSystem`s are primarily intended to model systems in discrete space. Modelling
continuous space systems with them is possible, but requires the user to determine the discretisation
(the space). Better support for continuous space models is a work in progress.
- Catalyst contains extensive documentation on spatial modelling, which can be found [here](https://docs.sciml.ai/Catalyst/stable/spatial_modelling/dspace_reaction_systems/).
"""
struct DiscreteSpaceReactionSystem{Q, R, S, T} <: MT.AbstractSystem
    # Input values.
    """The (non-spatial) reaction system within each vertex."""
    reactionsystem::ReactionSystem{Q}
    """The spatial reactions defined between individual vertices."""
    spatial_reactions::Vector{R}
    """The discrete space on which the (discrete) spatial system is defined."""
    dspace::S

    # Derived values.
    """The number of vertices (compartments)."""
    num_verts::Int64
    """The number of edges."""
    num_edges::Int64
    """The number of species."""
    num_species::Int64

    """List of species that may move spatially."""
    spatial_species::Vector{SymbolicT}
    """
    All parameters related to the discrete space reaction system
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
    An iterator over all the discrete space's edges. Currently, the format is always a Vector{Pair{Int64,Int64}}.
    However, in the future, different types could potentially be used for different types of discrete space
    (E.g. for a Cartesian grid, we do not technically need to enumerate each edge)
    """
    edge_iterator::T

    # Other values.
    "The name of the discrete space reaction system. Typically taken directly from the base `ReactionSystem`."
    name::Symbol

    function DiscreteSpaceReactionSystem(rs::ReactionSystem{Q}, spatial_reactions::Vector{R},
            dspace::S, num_verts::Int64, num_edges::Int64, edge_iterator::T;
            name::Symbol = MT.nameof(rs)) where {Q, R, S, T}
        # Error checks.
        if !(R <: AbstractSpatialReaction)
            throw(ArgumentError("The second argument must be a vector of AbstractSpatialReaction subtypes."))
        end
        if !iscomplete(rs)
            throw(ArgumentError("A non-complete `ReactionSystem` was used as input, this is not permitted."))
        end
        if length(reactions(rs)) != length(equations(rs))
            throw(ArgumentError("The `ReactionSystem` used as input contain equations (in addition to reactions). This is not permitted."))
        end
        if length(species(rs)) != length(unknowns(rs))
            throw(ArgumentError("The `ReactionSystem` used as input contain variable unknowns (in addition to species unknowns). This is not permitted (the input `ReactionSystem` must contain species unknowns only)."))
        end
        if !isempty(MT.continuous_events(rs)) || !isempty(MT.discrete_events(rs))
            throw(ArgumentError("The `ReactionSystem` used as input to `DiscreteSpaceReactionSystem` contain events. These will be ignored in any simulations based on the created `DiscreteSpaceReactionSystem`."))
        end
        if !isnothing(MT.get_parent(rs)) && !isempty(MT.get_systems(MT.get_parent(rs)))
            @warn "The `ReactionSystem` used as input to `DiscreteSpaceReactionSystem` was originally created as a hierarchical model. While this won't necessarily result in errors, it has not been well-tested, and is not recommended."
        end
        if !isempty(MT.observed(rs))
            @warn "The `ReactionSystem` used as input to `DiscreteSpaceReactionSystem` contain observables. It will not be possible to access these from the created `DiscreteSpaceReactionSystem`."
        end

        # Computes the species which are parts of spatial reactions. Also counts the total number of
        # species types.
        if isempty(spatial_reactions)
            spat_species = SymbolicT[]
        else
            spat_species = unique(reduce(vcat,
                [spatial_species(sr) for sr in spatial_reactions]))
        end
        num_species = length(unique([species(rs); spat_species]))

        # Computes the sets of vertex, edge, and all, parameters.
        rs_edge_parameters = filter(isedgeparameter, parameters(rs))
        if isempty(spatial_reactions)
            srs_edge_parameters = SymbolicT[]
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

        # Additional error checks.
        if any(is_array_symvar(symvar) for symvar in [ps; species(rs)])
            throw(ArgumentError("Some species and/or parameters used to create the `DiscreteSpaceReactionSystem` are array variables ($(filter(is_array_symvar, [ps; species(rs)]))). This is currently not supported."))
        end

        return new{Q, R, S, T}(
            rs, spatial_reactions, dspace, num_verts, num_edges, num_species,
            spat_species, ps, vertex_parameters, edge_parameters, edge_iterator)
    end
end

# Checks if a variable is a (non-/)scalarised array symbolic variable.
is_array_symvar(sym) = SymbolicUtils.is_array_shape(SymbolicUtils.shape(sym)) || (iscall(sym) && operation(sym) === getindex)

# Creates a DiscreteSpaceReactionSystem from a (directed) Graph dspace (graph grid).
function DiscreteSpaceReactionSystem(rs, srs, dspace::DiGraph; kwargs...)
    num_verts = nv(dspace)
    num_edges = ne(dspace)
    edge_iterator = [e.src => e.dst for e in edges(dspace)]
    return DiscreteSpaceReactionSystem(rs, srs, dspace, num_verts, num_edges, edge_iterator; kwargs...)
end
# Creates a DiscreteSpaceReactionSystem from a (undirected) Graph dspace (graph grid).
function DiscreteSpaceReactionSystem(rs, srs, dspace::SimpleGraph; kwargs...)
    DiscreteSpaceReactionSystem(rs, srs, DiGraph(dspace); kwargs...)
end

# Creates a DiscreteSpaceReactionSystem from a CartesianGrid dspace (cartesian grid) or a Boolean Array
# dspace (masked grid). These two are quite similar, so much code can be reused in a single interface.
function DiscreteSpaceReactionSystem(rs, srs, dspace::GridLattice{N, T};
        diagonal_connections = false, kwargs...) where {N, T}
    # Error checks.
    (N > 3) && error("Grids of higher dimension than 3 is currently not supported.")

    # Computes the number of vertices and edges. The two grid types (Cartesian and masked) each
    # uses their own function for this.
    num_verts = count_verts(dspace)
    num_edges = count_edges(dspace; diagonal_connections)

    # Finds all the grid's edges. First computer `fspat_to_grid_idx` which is a vector which takes
    # each vertex's flat (scalar) index to its grid index (e.g. (3,5) for a 2d grid). Next compute
    # `grid_to_fspat_idx` which is an array (of the same size as the grid) that does the reverse conversion.
    # Especially for masked grids these have non-trivial forms.
    fspat_to_grid_idx, grid_to_fspat_idx = get_index_converters(dspace, num_verts)
    # Simultaneously iterates through all vertices' flat and grid indices. Finds the (grid) indices
    # of their neighbours (different approaches for the two grid types). Converts these to flat
    # indices and adds the edges to `edge_iterator`.
    cur_vert = 0
    g_size = grid_size(dspace)
    edge_iterator = Vector{Pair{Int64, Int64}}(undef, num_edges)
    for (fspat_idx, grid_idx) in enumerate(fspat_to_grid_idx)
        for neighbour_grid_idx in get_neighbours(dspace, grid_idx, g_size;
            diagonal_connections)
            cur_vert += 1
            edge_iterator[cur_vert] = fspat_idx => grid_to_fspat_idx[neighbour_grid_idx...]
        end
    end

    return DiscreteSpaceReactionSystem(rs, srs, dspace, num_verts, num_edges, edge_iterator; kwargs...)
end

### DiscreteSpaceReactionSystem Helper Functions ###
# Note, most of these are specifically for (Cartesian or masked) grids, we call them `grid`, not `dspace`.

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
# corresponds to the fifth vertex, then `fspat_to_grid_idx[5] = (3,2)` and `grid_to_fspat_idx[3,2] = 5`.
function get_index_converters(grid::GridLattice{N, T}, num_verts) where {N, T}
    fspat_to_grid_idx = Vector{typeof(grid_size(grid))}(undef, num_verts)
    grid_to_fspat_idx = Array{Int64}(undef, grid_size(grid))

    # Loops through the flat and grid indices simultaneously, adding them to their respective converters.
    cur_fspat_idx = 0
    for grid_idx in get_grid_indices(grid)
        # For a masked grid, grid points with `false` values are skipped.
        (grid isa Array{Bool}) && (!grid[grid_idx]) && continue

        cur_fspat_idx += 1
        fspat_to_grid_idx[cur_fspat_idx] = grid_idx
        grid_to_fspat_idx[grid_idx] = cur_fspat_idx
    end
    return fspat_to_grid_idx, grid_to_fspat_idx
end

# For a vertex's grid index, and a dspace, returns the grid indices of all its (valid) neighbours.
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

### DiscreteSpaceReactionSystem-specific Getters ###

# Basic getters (because `DiscreteSpaceReactionSystem`s are `AbstractSystem`s), normal `dsrs.field` does not
# work and these getters must be used throughout all code.
"""
    reactionsystem(dsrs::DiscreteSpaceReactionSystem)

Returns the non-spatial `ReactionSystem` stored in a `DiscreteSpaceReactionSystem`.
"""
reactionsystem(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :reactionsystem)

"""
    spatial_reactions(dsrs::DiscreteSpaceReactionSystem)

Returns a vector with all the spatial reactions stored in a `DiscreteSpaceReactionSystem`.
"""
spatial_reactions(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :spatial_reactions)

"""
    dspace(dsrs::DiscreteSpaceReactionSystem)

Returns the dspace (i.e. discrete space) stored in a `DiscreteSpaceReactionSystem`.
"""
dspace(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :dspace)

"""
    num_verts(dsrs::DiscreteSpaceReactionSystem)

Returns the number of vertices (i.e. compartments) in the discrete space stored in a `DiscreteSpaceReactionSystem`.
"""
num_verts(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :num_verts)

"""
    num_edges(dsrs::DiscreteSpaceReactionSystem)

Returns the number of edges (i.e. connections between vertices) in the discrete space stored in a
`DiscreteSpaceReactionSystem`.
"""
num_edges(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :num_edges)

"""
    num_species(dsrs::DiscreteSpaceReactionSystem)

Returns the number of species that a `DiscreteSpaceReactionSystem` contains.
"""
num_species(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :num_species)

"""
    spatial_species(dsrs::DiscreteSpaceReactionSystem)

Returns the species that can move spatially in a `DiscreteSpaceReactionSystem`.
"""
spatial_species(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :spatial_species)

# Returns the parameters in a `DiscreteSpaceReactionSystem`
MT.parameters(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :parameters)

"""
    vertex_parameters(dsrs::DiscreteSpaceReactionSystem)

Returns all the parameters of a `DiscreteSpaceReactionSystem` whose values are tied to vertices.
"""
vertex_parameters(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :vertex_parameters)

"""
    edge_parameters(dsrs::DiscreteSpaceReactionSystem)

Returns all the parameters of a `DiscreteSpaceReactionSystem` whose values are tied to edges.
"""
edge_parameters(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :edge_parameters)

"""
    edge_iterator(dsrs::DiscreteSpaceReactionSystem)

Returns an iterator over all of the edges in the discrete space stored in a `DiscreteSpaceReactionSystem`. Each
edge is a `Pair{Int64, Int64}`, taking the source vertex to the destination vertex.
"""
edge_iterator(dsrs::DiscreteSpaceReactionSystem) = getfield(dsrs, :edge_iterator)

# Non-trivial getters.
"""
    is_transport_system(dsrs::DiscreteSpaceReactionSystem)

Returns `true` if all spatial reactions in `dsrs` are `TransportReaction`s.
"""
function is_transport_system(dsrs::DiscreteSpaceReactionSystem)
    return all(sr -> sr isa TransportReaction, spatial_reactions(dsrs))
end

"""
    has_cartesian_dspace(dsrs::DiscreteSpaceReactionSystem)

Returns `true` if `dsrs` was created using a cartesian grid discrete space (e.g. created via `CartesianGrid(5,5)`).
Otherwise, returns `false`.
"""
function has_cartesian_dspace(dsrs::DiscreteSpaceReactionSystem)
    dspace(dsrs) isa CartesianGridRej{N, T} where {N, T}
end

"""
    has_masked_dspace(dsrs::DiscreteSpaceReactionSystem)

Returns `true` if `dsrs` was created using a masked grid discrete space (e.g. created via `[true true; true false]`).
Otherwise, returns `false`.
"""
has_masked_dspace(dsrs::DiscreteSpaceReactionSystem) = dspace(dsrs) isa Array{Bool, N} where {N}

"""
    has_grid_dspace(dsrs::DiscreteSpaceReactionSystem)

Returns `true` if `dsrs` was created using a cartesian or masked grid discrete space. Otherwise, returns `false`.
"""
function has_grid_dspace(dsrs::DiscreteSpaceReactionSystem)
    return has_cartesian_dspace(dsrs) || has_masked_dspace(dsrs)
end

"""
    has_graph_dspace(dsrs::DiscreteSpaceReactionSystem)

Returns `true` if `dsrs` was created using a graph grid discrete space (e.g. created via `path_graph(5)`).
Otherwise, returns `false`.
"""
has_graph_dspace(dsrs::DiscreteSpaceReactionSystem) = dspace(dsrs) isa SimpleDiGraph

"""
    grid_size(dsrs::DiscreteSpaceReactionSystem)

Returns the size of `dsrs`'s discrete space (only if it is a cartesian or masked grid discrete space).
E.g. for a discrete space `CartesianGrid(4,6)`, `(4,6)` is returned.
"""
grid_size(dsrs::DiscreteSpaceReactionSystem) = grid_size(dspace(dsrs))
grid_size(dspace::CartesianGridRej{N, T}) where {N, T} = dspace.dims
grid_size(dspace::Array{Bool, N}) where {N} = size(dspace)
function grid_size(dspace::Graphs.AbstractGraph)
    throw(ArgumentError("Grid size is only defined for DiscreteSpaceReactionSystems with grid-based discrete spaces (not graph-based)."))
end

"""
    grid_dims(dsrs::DiscreteSpaceReactionSystem)

Returns the number of dimensions of `dsrs`'s discrete space (only if it is a cartesian or masked grid discrete space).
The output is either `1`, `2`, or `3`.
"""
grid_dims(dsrs::DiscreteSpaceReactionSystem) = grid_dims(dspace(dsrs))
grid_dims(dspace::GridLattice{N, T}) where {N, T} = return N
function grid_dims(dspace::Graphs.AbstractGraph)
    throw(ArgumentError("Grid dimensions is only defined for DiscreteSpaceReactionSystems with grid-based discrete spaces (not graph-based)."))
end

"""
    get_dspace_graph(dsrs::DiscreteSpaceReactionSystem)

Returns dsrs's discrete space as a graph. Currently does not work for Cartesian discrete spaces.
"""
function get_dspace_graph(dsrs::DiscreteSpaceReactionSystem)
    has_graph_dspace(dsrs) && return dspace(dsrs)
    return Graphs.SimpleGraphFromIterator(Graphs.SimpleEdge(e[1], e[2])
    for e in edge_iterator(dsrs))
end

### Catalyst-based Getters ###

# Get all species.
function species(dsrs::DiscreteSpaceReactionSystem)
    unique([species(reactionsystem(dsrs)); spatial_species(dsrs)])
end

# Generic ones (simply forwards call to the non-spatial system).
reactions(dsrs::DiscreteSpaceReactionSystem) = reactions(reactionsystem(dsrs))

### ModelingToolkit-based Getters ###

# Generic ones (simply forwards call to the non-spatial system)
# The `parameters` MTK getter have a specialised accessor for DiscreteSpaceReactionSystems.
MT.nameof(lrs::LatticeReactionSystem) = MT.nameof(reactionsystem(lrs))
MT.get_iv(dsrs::DiscreteSpaceReactionSystem) = MT.get_iv(reactionsystem(dsrs))
MT.equations(dsrs::DiscreteSpaceReactionSystem) = MT.equations(reactionsystem(dsrs))
MT.unknowns(dsrs::DiscreteSpaceReactionSystem) = MT.unknowns(reactionsystem(dsrs))
MT.get_metadata(dsrs::DiscreteSpaceReactionSystem) = MT.get_metadata(reactionsystem(dsrs))

# Lattice reaction systems should not be combined with compositional modelling.
# Maybe these should be allowed anyway? Still feel a bit weird
function MT.get_eqs(dsrs::DiscreteSpaceReactionSystem)
    MT.get_eqs(reactionsystem(dsrs))
end
function MT.get_unknowns(dsrs::DiscreteSpaceReactionSystem)
    MT.get_unknowns(reactionsystem(dsrs))
end
function MT.get_ps(dsrs::DiscreteSpaceReactionSystem)
    MT.get_ps(reactionsystem(dsrs))
end

# Technically should not be used, but has to be declared for the `show` function to work.
function MT.get_systems(dsrs::DiscreteSpaceReactionSystem)
    return []
end

# Other non-relevant getters.
function MT.independent_variables(dsrs::DiscreteSpaceReactionSystem)
    MT.independent_variables(reactionsystem(dsrs))
end

### Edge Parameter Value Generators ###

"""
    make_edge_p_values(dsrs::DiscreteSpaceReactionSystem, make_edge_p_value::Function)

Generates edge parameter values for a discrete space reaction system. Only works for (Cartesian or masked)
grid discrete spaces (without diagonal adjacencies).

Input:
- `dsrs`: The discrete space reaction system for which values should be generated.
- `make_edge_p_value`: a function describing a rule for generating the edge parameter values.

Output:
    - `ep_vals`: A sparse matrix of size (num_verts,num_verts) (where num_verts is the number of
    vertices in `dsrs`). Here, `eps[i,j]` is filled only if there is an edge going from vertex i to
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
space = CartesianGrid((5,5))
dsrs = DiscreteSpaceReactionSystem(rn, [tr], space)

function make_edge_p_value(src_vert, dst_vert)
    if src_vert == (1,1) && dst_vert == (1,2)
        return 1.0
    else
        return 0.1
    end
end

D_vals = make_edge_p_values(dsrs, make_edge_p_value)
```
"""
function make_edge_p_values(dsrs::DiscreteSpaceReactionSystem, make_edge_p_value::Function)
    if has_graph_dspace(dsrs)
        error("The `make_edge_p_values` function is only meant for discrete spaces with (Cartesian or masked) grid structures. It cannot be applied to graph discrete spaces.")
    end

    # Makes the flat to index grid converts. Predeclared the edge parameter value sparse matrix.
    fspat_to_grid_idx = get_index_converters(dspace(dsrs), num_verts(dsrs))[1]
    values = spzeros(num_verts(dsrs), num_verts(dsrs))

    # Loops through all edges, and applies the value function to these.
    for e in edge_iterator(dsrs)
        # This extra step is needed to ensure that `0` is stored if make_edge_p_value yields a 0.
        # If not, then the sparse matrix simply becomes empty in that position.
        values[e[1], e[2]] = eps()

        values[e[1], e[2]] = make_edge_p_value(fspat_to_grid_idx[e[1]],
            fspat_to_grid_idx[e[2]])
    end

    return values
end

"""
    make_directed_edge_values(dsrs::DiscreteSpaceReactionSystem, x_vals::Tuple{T,T}, y_vals::Tuple{T,T} = (undef,undef),
                     z_vals::Tuple{T,T} = (undef,undef)) where {T}

Generates edge parameter values for a discrete space reaction system. Only works for (Cartesian or masked)
grid discrete spaces (without diagonal adjacencies). Each dimension (x, and possibly y and z), and
direction has assigned its own constant edge parameter value.

Input:
    - `dsrs`: The discrete space reaction system for which values should be generated.
    - `x_vals::Tuple{T,T}`: The values in the increasing (from a lower x index to a higher x index)
    and decreasing (from a higher x index to a lower x index) direction along the x dimension.
    - `y_vals::Tuple{T,T}`: The values in the increasing and decreasing direction along the y dimension.
    Should only be used for 2 and 3-dimensional grids.
    - `z_vals::Tuple{T,T}`: The values in the increasing and decreasing direction along the z dimension.
    Should only be used for 3-dimensional grids.

Output:
    - `ep_vals`: A sparse matrix of size (num_verts,num_verts) (where num_verts is the number of
    vertices in `dsrs`). Here, `eps[i,j]` is filled only if there is an edge going from vertex i to
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
space = CartesianGrid((5,5))
dsrs = DiscreteSpaceReactionSystem(rn, [tr], space)

D_vals = make_directed_edge_values(dsrs, (0.1, 0.1), (0.1, 0.0))
```
Here, since we have a 2d grid, we only provide the first two Tuples to `make_directed_edge_values`.
"""
function make_directed_edge_values(dsrs::DiscreteSpaceReactionSystem, x_vals::Tuple{T, T},
        y_vals::Union{Nothing, Tuple{T, T}} = nothing,
        z_vals::Union{Nothing, Tuple{T, T}} = nothing) where {T}
    # Error checks.
    if has_graph_dspace(dsrs)
        error("The `make_directed_edge_values` function is only meant for discrete spaces with (Cartesian or masked) grid structures. It cannot be applied to graph discrete spaces.")
    end
    if count(!isnothing(flow) for flow in [x_vals, y_vals, z_vals]) != grid_dims(dsrs)
        error("You must provide flows in the same number of dimensions as your discrete space has dimensions. The discrete space have $(grid_dims(dsrs)), and flows where provided for $(count(isnothing(flow) for flow in [x_vals, y_vals, z_vals])) dimensions.")
    end

    # Declares a function that assigns the correct flow value for a given edge.
    function directed_vals_func(src_vert, dst_vert)
        if count(src_vert .== dst_vert) != (grid_dims(dsrs) - 1)
            error("The `make_directed_edge_values` function can only be applied to discrete spaces with rigid (non-diagonal) grid structure. It is being evaluated for the edge from $(src_vert) to $(dst_vert), which does not seem directly adjacent on a grid.")
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
    return make_edge_p_values(dsrs, directed_vals_func)
end
