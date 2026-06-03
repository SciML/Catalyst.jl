### Preparations ###

# Fetch packages.
using Catalyst, Graphs, OrdinaryDiffEqTsit5, Test

# Fetch test networks.
include("../spatial_test_networks.jl")

# Pre-declares a set of grids.
grids = [very_small_2d_cartesian_grid, very_small_2d_masked_grid, very_small_2d_graph_grid]


### Tests DiscreteSpaceReactionSystem Getters Correctness ###

# Test case 1.
let
    rs = @reaction_network begin
        (p, 1), 0 <--> X
    end
    tr = @transport_reaction d X
    for grid in grids
        dsrs = DiscreteSpaceReactionSystem(rs, [tr], grid)

        @unpack X, p = rs
        d = edge_parameters(dsrs)[1]
        @test issetequal(species(dsrs), [X])
        @test issetequal(spatial_species(dsrs), [X])
        @test issetequal(parameters(dsrs), [p, d])
        @test issetequal(vertex_parameters(dsrs), [p])
        @test issetequal(edge_parameters(dsrs), [d])
    end
end

# Test case 2.
let
    rs = @reaction_network begin
        @parameters p1 p2 [edgeparameter=true]
    end
    @unpack p1, p2 = rs

    @test !isedgeparameter(p1)
    @test isedgeparameter(p2)
end

# Test case 3.
let
    rs = @reaction_network begin
        @parameters pX pY dX [edgeparameter=true] dY
        (pX, 1), 0 <--> X
        (pY, 1), 0 <--> Y
    end
    tr_1 = @transport_reaction dX X
    tr_2 = @transport_reaction dY Y
    for grid in grids
        dsrs = DiscreteSpaceReactionSystem(rs, [tr_1, tr_2], grid)

        @unpack X, Y, pX, pY, dX, dY = rs
        @test issetequal(species(dsrs), [X, Y])
        @test issetequal(spatial_species(dsrs), [X, Y])
        @test issetequal(parameters(dsrs), [pX, pY, dX, dY])
        @test issetequal(vertex_parameters(dsrs), [pX, pY, dY])
        @test issetequal(edge_parameters(dsrs), [dX])
    end
end

# Test case 4.
let
    rs = @reaction_network begin
        @parameters dX p
        (pX, 1), 0 <--> X
        (pY, 1), 0 <--> Y
    end
    tr_1 = @transport_reaction dX X
    for grid in grids
        dsrs = DiscreteSpaceReactionSystem(rs, [tr_1], grid)

        @unpack dX, p, X, Y, pX, pY = rs
        @test issetequal(species(dsrs), [X, Y])
        @test issetequal(spatial_species(dsrs), [X])
        @test issetequal(parameters(dsrs), [dX, p, pX, pY])
        @test issetequal(vertex_parameters(dsrs), [dX, p, pX, pY])
        @test issetequal(edge_parameters(dsrs), [])
    end
end

# Test case 5.
let
    rs = @reaction_network begin
        @species W(t)
        @parameters pX pY dX [edgeparameter=true] dY
        (pX, 1), 0 <--> X
        (pY, 1), 0 <--> Y
        (pZ, 1), 0 <--> Z
        (pV, 1), 0 <--> V
    end
    @unpack dX, X, V = rs
    @parameters dV dW
    t = default_t()
    @species W(t)
    tr_1 = TransportReaction(dX, X)
    tr_2 = @transport_reaction dY Y
    tr_3 = @transport_reaction dZ Z
    tr_4 = TransportReaction(dV, V)
    tr_5 = TransportReaction(dW, W)
    for grid in grids
        dsrs = DiscreteSpaceReactionSystem(rs, [tr_1, tr_2, tr_3, tr_4, tr_5], grid)

        @unpack pX, pY, pZ, pV, dX, dY, X, Y, Z, V = rs
        dZ, dV, dW = edge_parameters(dsrs)[2:end]
        @test issetequal(species(dsrs), [W, X, Y, Z, V])
        @test issetequal(spatial_species(dsrs), [X, Y, Z, V, W])
        @test issetequal(parameters(dsrs), [pX, pY, dX, dY, pZ, pV, dZ, dV, dW])
        @test issetequal(vertex_parameters(dsrs), [pX, pY, dY, pZ, pV])
        @test issetequal(edge_parameters(dsrs), [dX, dZ, dV, dW])
    end
end

# Test case 6.
let
    rs = @reaction_network customname begin
        (p, 1), 0 <--> X
    end
    tr = @transport_reaction d X
    for grid in grids
        dsrs = DiscreteSpaceReactionSystem(rs, [tr], grid)

        @test nameof(dsrs) == :customname
    end
end

# Tests using various more obscure types of getters.
let
    # Create DiscreteSpaceReactionSystems.
    t = default_t()
    @parameters p d kB kD
    @species X(t) X2(t)
    rxs = [
        Reaction(p, [], [X])
        Reaction(d, [X], [])
        Reaction(kB, [X], [X2], [2], [1])
        Reaction(kD, [X2], [X], [1], [2])
    ]
    @named rs = ReactionSystem(rxs, t; metadata = [MiscSystemData => "Metadata string"])
    rs = complete(rs)
    tr = @transport_reaction D X2
    dsrs = DiscreteSpaceReactionSystem(rs, [tr], small_2d_cartesian_grid)

    # Generic ones (simply forwards call to the non-spatial system).
    @test isequal(reactions(dsrs), rxs)
    @test isequal(nameof(dsrs), :rs)
    @test isequal(ModelingToolkitBase.get_iv(dsrs), t)
    @test isequal(equations(dsrs), rxs)
    @test isequal(unknowns(dsrs), [X, X2])
    @test isequal(ModelingToolkitBase.getmetadata(dsrs, MiscSystemData, nothing), "Metadata string")
    @test isequal(ModelingToolkitBase.get_eqs(dsrs), rxs)
    @test isequal(ModelingToolkitBase.get_unknowns(dsrs), [X, X2])
    @test isequal(ModelingToolkitBase.get_ps(dsrs), [p, d, kB, kD])
    @test isequal(ModelingToolkitBase.get_systems(dsrs), [])
    @test isequal(independent_variables(dsrs), [t])
end

### Tests Error generation ###

# Network where diffusion species is not declared in non-spatial network.
let
    rs = @reaction_network begin
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D Y
    for grid in grids
        @test_throws ErrorException DiscreteSpaceReactionSystem(rs, [tr], grid)
    end
end

# Network where the rate depend on a species
let
    rs = @reaction_network begin
        @species Y(t)
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D*Y X
    for grid in grids
        @test_throws ErrorException DiscreteSpaceReactionSystem(rs, [tr], grid)
    end
end

# Network with edge parameter in non-spatial reaction rate.
let
    rs = @reaction_network begin
        @parameters p [edgeparameter=true]
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D X
    for grid in grids
        @test_throws ErrorException DiscreteSpaceReactionSystem(rs, [tr], grid)
    end
end

# Network where metadata has been added in rs (which is not seen in transport reaction).
let
    rs = @reaction_network begin
        @species X(t) [description="Species with added metadata"]
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D X
    for grid in grids
        @test_throws ErrorException DiscreteSpaceReactionSystem(rs, [tr], grid)

        rs = @reaction_network begin
            @parameters D [description="Parameter with added metadata"]
            (p, d), 0 <--> X
        end
        tr = @transport_reaction D X
        @test_throws ErrorException DiscreteSpaceReactionSystem(rs, [tr], grid)
    end
end

# The second argument must be a vector of AbstractSpatialReaction subtypes.
let
    # Define an invalid spatial reaction type (not a subtype of AbstractSpatialReaction)
    struct InvalidSpatialReactionType end

    # Attempt to create the DiscreteSpaceReactionSystem with InvalidSpatialReactionType
    for grid in grids
        @test_throws ArgumentError DiscreteSpaceReactionSystem(binding_system, [], grid)
        @test_throws ArgumentError DiscreteSpaceReactionSystem(binding_system, [InvalidSpatialReactionType()], grid)
    end
end

# Tests various networks with non-permitted content.
let
    tr = @transport_reaction D X

    # Variable unknowns.
    rs1 = @reaction_network begin
        @variables V(t)
        (p,d), 0 <--> X
    end
    @test_throws ArgumentError DiscreteSpaceReactionSystem(rs1, [tr], short_path)

    # Non-reaction equations.
    rs2 = @reaction_network begin
        @equations D(V) ~ X - V
        (p,d), 0 <--> X
    end
    @test_throws ArgumentError DiscreteSpaceReactionSystem(rs2, [tr], short_path)

    # Events.
    rs3 = @reaction_network begin
        @discrete_events [1.0] => [p => p + 1]
        (p,d), 0 <--> X
    end
    @test_throws ArgumentError DiscreteSpaceReactionSystem(rs3, [tr], short_path)

    # Observables (only generates a warning).
    rs4 = @reaction_network begin
        @observables X2 ~ 2X
        (p,d), 0 <--> X
    end
    @test_logs (:warn, r"The `ReactionSystem` used as input to `DiscreteSpaceReactionSystem` contain observables. It *") match_mode=:any DiscreteSpaceReactionSystem(rs4, [tr], short_path)
end

# Tests for hierarchical input system (should yield a warning).
let
    t = default_t()
    @parameters d D
    @species X(t)
    rxs = [Reaction(d, [X], [])]
    @named rs1 = ReactionSystem(rxs, t)
    @named rs2 = ReactionSystem(rxs, t; systems = [rs1])
    rs2 = complete(rs2)
    @test_logs (:warn, r"The `ReactionSystem` used as input to `DiscreteSpaceReactionSystem` was originally created as a hierarchical model. While *") match_mode=:any DiscreteSpaceReactionSystem(rs2, [TransportReaction(D, X)], CartesianGrid((2,2)))
end

# Tests for non-complete input `ReactionSystem`.
let
    tr = @transport_reaction D X
    rs = @network_component begin
        (p,d), 0 <--> X
    end
    @test_throws ArgumentError DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,2)))
end

# Tests for array parameters/species.
let
    tr = @transport_reaction D Y

    rs1 = @reaction_network begin
        @species X(t)[1:2] Y(t)
        (k1,k2), X[1] <--> X[2]
    end
    @test_throws ArgumentError DiscreteSpaceReactionSystem(rs1, [tr], CartesianGrid((2,2)))

    rs2 =  @reaction_network begin
        @species Y(t)
        @parameters k[1:2,1:2]
        (k[1,1],k[1,2]), X11 <--> X12
        (k[2,1],k[2,2]), X21 <--> X22
    end
    @test_throws ArgumentError DiscreteSpaceReactionSystem(rs2, [tr], CartesianGrid((2,2)))
end

### Tests Grid Vertex and Edge Number Computation ###

# Tests that the correct numbers are computed for num_edges.
let
    # Function counting the values in an iterator by stepping through it.
    function iterator_count(iterator)
        count = 0
        foreach(e -> count+=1, iterator)
        return count
    end

    # Cartesian and masked grid (test diagonal edges as well).
    for space in [small_1d_cartesian_grid, small_2d_cartesian_grid, small_3d_cartesian_grid,
                random_1d_masked_grid, random_2d_masked_grid, random_3d_masked_grid]
        dsrs1 = DiscreteSpaceReactionSystem(SIR_system, SIR_srs_1, space)
        dsrs2 = DiscreteSpaceReactionSystem(SIR_system, SIR_srs_1, space; diagonal_connections=true)
        @test num_edges(dsrs1) == iterator_count(edge_iterator(dsrs1))
        @test num_edges(dsrs2) == iterator_count(edge_iterator(dsrs2))
    end

    # Graph grids (cannot test diagonal connections).
    for space in [small_2d_graph_grid, small_3d_graph_grid, undirected_cycle, small_directed_cycle, unconnected_graph]
        dsrs1 = DiscreteSpaceReactionSystem(SIR_system, SIR_srs_1, space)
        @test num_edges(dsrs1) == iterator_count(edge_iterator(dsrs1))
    end
end

### Tests Edge Value Computation Helper Functions ###

# Checks that we compute the correct values across various types of grids.
let
    # Prepares the model and the function that determines the edge values.
    rn = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    function make_edge_p_value(src_vert, dst_vert)
        return prod(src_vert) + prod(dst_vert)
    end

    # Loops through a variety of grids, checks that `make_edge_p_values` yields the correct values.
    for space in [small_1d_cartesian_grid, small_2d_cartesian_grid, small_3d_cartesian_grid,
                 small_1d_masked_grid, small_2d_masked_grid, small_3d_masked_grid,
                 random_1d_masked_grid, random_2d_masked_grid, random_3d_masked_grid]
        dsrs = DiscreteSpaceReactionSystem(rn, [tr], space)
        fspat_to_grid_idx = Catalyst.get_index_converters(dspace(dsrs), num_verts(dsrs))[1]
        edge_values = make_edge_p_values(dsrs, make_edge_p_value)

        for e in edge_iterator(dsrs)
            @test edge_values[e[1], e[2]] == make_edge_p_value(fspat_to_grid_idx[e[1]], fspat_to_grid_idx[e[2]])
        end
    end
end

# Checks that all species end up in the correct place in a pure flow system (checking various dimensions).
let
    # Prepares a system with a single species which is transported only.
    rn = @reaction_network begin
        @species X(t)
    end
    n = 5
    tr = @transport_reaction D X
    tspan = (0.0, 1000.0)
    u0 = [:X => 1.0]

    # Checks the 1d case.
    dsrs = DiscreteSpaceReactionSystem(rn, [tr], CartesianGrid(n))
    ps = [:D => make_directed_edge_values(dsrs, (10.0, 0.0))]
    oprob = ODEProblem(dsrs, u0, tspan, ps)
    @test isapprox(solve(oprob, Tsit5()).u[end][5], n, rtol=1e-6)

    # Checks the 2d case (both with 1d and 2d flow).
    dsrs = DiscreteSpaceReactionSystem(rn, [tr], CartesianGrid((n,n)))

    ps = [:D => make_directed_edge_values(dsrs, (1.0, 0.0), (0.0, 0.0))]
    oprob = ODEProblem(dsrs, u0, tspan, ps)
    @test all(isapprox.(solve(oprob, Tsit5()).u[end][5:5:25], n, rtol=1e-6))

    ps = [:D => make_directed_edge_values(dsrs, (1.0, 0.0), (1.0, 0.0))]
    oprob = ODEProblem(dsrs, u0, tspan, ps)
    @test isapprox(solve(oprob, Tsit5()).u[end][25], n^2, rtol=1e-6)

    # Checks the 3d case (both with 1d and 2d flow).
    dsrs = DiscreteSpaceReactionSystem(rn, [tr], CartesianGrid((n,n,n)))

    ps = [:D => make_directed_edge_values(dsrs, (1.0, 0.0), (0.0, 0.0), (0.0, 0.0))]
    oprob = ODEProblem(dsrs, u0, tspan, ps)
    @test all(isapprox.(solve(oprob, Tsit5()).u[end][5:5:125], n, rtol=1e-6))

    ps = [:D => make_directed_edge_values(dsrs, (1.0, 0.0), (1.0, 0.0), (0.0, 0.0))]
    oprob = ODEProblem(dsrs, u0, tspan, ps)
    @test all(isapprox.(solve(oprob, Tsit5()).u[end][25:25:125], n^2, rtol=1e-6))

    ps = [:D => make_directed_edge_values(dsrs, (1.0, 0.0), (1.0, 0.0), (1.0, 0.0))]
    oprob = ODEProblem(dsrs, u0, tspan, ps)
    @test isapprox(solve(oprob, Tsit5()).u[end][125], n^3, rtol=1e-6)
end

# Checks that erroneous input yields errors.
let
    rn = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    tspan = (0.0, 10000.0)
    make_edge_p_value(src_vert, dst_vert) = rand()

    # Graph grids.
    dsrs = DiscreteSpaceReactionSystem(rn, [tr], path_graph(5))
    @test_throws Exception make_edge_p_values(dsrs, make_edge_p_value,)
    @test_throws Exception make_directed_edge_values(dsrs, (1.0, 0.0))

    # Wrong dimensions to `make_directed_edge_values`.
    dsrs_1d = DiscreteSpaceReactionSystem(rn, [tr], CartesianGrid(5))
    dsrs_2d = DiscreteSpaceReactionSystem(rn, [tr], fill(true,5,5))
    dsrs_3d = DiscreteSpaceReactionSystem(rn, [tr], CartesianGrid((5,5,5)))

    @test_throws Exception make_directed_edge_values(dsrs_1d, (1.0, 0.0), (1.0, 0.0))
    @test_throws Exception make_directed_edge_values(dsrs_1d, (1.0, 0.0), (1.0, 0.0), (1.0, 0.0))
    @test_throws Exception make_directed_edge_values(dsrs_2d, (1.0, 0.0))
    @test_throws Exception make_directed_edge_values(dsrs_2d, (1.0, 0.0), (1.0, 0.0), (1.0, 0.0))
    @test_throws Exception make_directed_edge_values(dsrs_3d, (1.0, 0.0))
    @test_throws Exception make_directed_edge_values(dsrs_3d, (1.0, 0.0), (1.0, 0.0))
end
