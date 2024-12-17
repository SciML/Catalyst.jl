### Preparations ###

# Fetch packages.
using Catalyst, Graphs, OrdinaryDiffEqTsit5, Test

# Fetch test networks.
include("../spatial_test_networks.jl")

# Pre-declares a set of grids.
grids = [very_small_2d_cartesian_grid, very_small_2d_masked_grid, very_small_2d_graph_grid]


### Tests LatticeReactionSystem Getters Correctness ###

# Test case 1.
let
    rs = @reaction_network begin
        (p, 1), 0 <--> X
    end
    tr = @transport_reaction d X
    for grid in grids
        lrs = LatticeReactionSystem(rs, [tr], grid)

        @unpack X, p = rs
        d = edge_parameters(lrs)[1]
        @test issetequal(species(lrs), [X])
        @test issetequal(spatial_species(lrs), [X])
        @test issetequal(parameters(lrs), [p, d])
        @test issetequal(vertex_parameters(lrs), [p])
        @test issetequal(edge_parameters(lrs), [d])
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
        lrs = LatticeReactionSystem(rs, [tr_1, tr_2], grid)

        @unpack X, Y, pX, pY, dX, dY = rs
        @test issetequal(species(lrs), [X, Y])
        @test issetequal(spatial_species(lrs), [X, Y])
        @test issetequal(parameters(lrs), [pX, pY, dX, dY])
        @test issetequal(vertex_parameters(lrs), [pX, pY, dY])
        @test issetequal(edge_parameters(lrs), [dX])
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
        lrs = LatticeReactionSystem(rs, [tr_1], grid)

        @unpack dX, p, X, Y, pX, pY = rs
        @test issetequal(species(lrs), [X, Y])
        @test issetequal(spatial_species(lrs), [X])
        @test issetequal(parameters(lrs), [dX, p, pX, pY])
        @test issetequal(vertex_parameters(lrs), [dX, p, pX, pY])
        @test issetequal(edge_parameters(lrs), [])
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
        lrs = LatticeReactionSystem(rs, [tr_1, tr_2, tr_3, tr_4, tr_5], grid)

        @unpack pX, pY, pZ, pV, dX, dY, X, Y, Z, V = rs
        dZ, dV, dW = edge_parameters(lrs)[2:end]
        @test issetequal(species(lrs), [W, X, Y, Z, V])
        @test issetequal(spatial_species(lrs), [X, Y, Z, V, W])
        @test issetequal(parameters(lrs), [pX, pY, dX, dY, pZ, pV, dZ, dV, dW])
        @test issetequal(vertex_parameters(lrs), [pX, pY, dY, pZ, pV])
        @test issetequal(edge_parameters(lrs), [dX, dZ, dV, dW])
    end
end

# Test case 6.
let
    rs = @reaction_network customname begin
        (p, 1), 0 <--> X
    end
    tr = @transport_reaction d X
    for grid in grids
        lrs = LatticeReactionSystem(rs, [tr], grid)

        @test nameof(lrs) == :customname
    end
end

# Tests using various more obscure types of getters.
let
    # Create LatticeReactionsSystems.
    t = default_t()
    @parameters p d kB kD
    @species X(t) X2(t)
    rxs = [
        Reaction(p, [], [X])
        Reaction(d, [X], [])
        Reaction(kB, [X], [X2], [2], [1])
        Reaction(kD, [X2], [X], [1], [2])
    ]
    @named rs = ReactionSystem(rxs, t; metadata = "Metadata string")
    rs = complete(rs)
    tr = @transport_reaction D X2
    lrs = LatticeReactionSystem(rs, [tr], small_2d_cartesian_grid)

    # Generic ones (simply forwards call to the non-spatial system).
    @test isequal(reactions(lrs), rxs)
    @test isequal(nameof(lrs), :rs)
    @test isequal(ModelingToolkit.get_iv(lrs), t)
    @test isequal(equations(lrs), rxs)
    @test isequal(unknowns(lrs), [X, X2])
    @test isequal(ModelingToolkit.get_metadata(lrs), "Metadata string")
    @test isequal(ModelingToolkit.get_eqs(lrs), rxs)
    @test isequal(ModelingToolkit.get_unknowns(lrs), [X, X2])
    @test isequal(ModelingToolkit.get_ps(lrs), [p, d, kB, kD])
    @test isequal(ModelingToolkit.get_systems(lrs), [])
    @test isequal(independent_variables(lrs), [t])
end

### Tests Error generation ###

# Network where diffusion species is not declared in non-spatial network.
let
    rs = @reaction_network begin
        (p, d), 0 <--> X
    end
    tr = @transport_reaction D Y
    for grid in grids
        @test_throws ErrorException LatticeReactionSystem(rs, [tr], grid)
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
        @test_throws ErrorException LatticeReactionSystem(rs, [tr], grid)
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
        @test_throws ErrorException LatticeReactionSystem(rs, [tr], grid)
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
        @test_throws ErrorException LatticeReactionSystem(rs, [tr], grid)

        rs = @reaction_network begin
            @parameters D [description="Parameter with added metadata"]
            (p, d), 0 <--> X
        end
        tr = @transport_reaction D X
        @test_throws ErrorException LatticeReactionSystem(rs, [tr], grid)
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
    @test_throws ArgumentError LatticeReactionSystem(rs1, [tr], short_path)

    # Non-reaction equations.
    rs2 = @reaction_network begin
        @equations D(V) ~ X - V
        (p,d), 0 <--> X
    end
    @test_throws ArgumentError LatticeReactionSystem(rs2, [tr], short_path)

    # Events.
    rs3 = @reaction_network begin
        @discrete_events [1.0] => [p ~ p + 1]
        (p,d), 0 <--> X
    end
    @test_throws ArgumentError LatticeReactionSystem(rs3, [tr], short_path)

    # Observables (only generates a warning).
    rs4 = @reaction_network begin
        @observables X2 ~ 2X
        (p,d), 0 <--> X
    end
    @test_logs (:warn, r"The `ReactionSystem` used as input to `LatticeReactionSystem` contain observables. It *") match_mode=:any LatticeReactionSystem(rs4, [tr], short_path)
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
    @test_logs (:warn, r"The `ReactionSystem` used as input to `LatticeReactionSystem` was originally created as a hierarchical model. While *") match_mode=:any LatticeReactionSystem(rs2, [TransportReaction(D, X)], CartesianGrid((2,2)))
end

# Tests for non-complete input `ReactionSystem`.
let
    tr = @transport_reaction D X
    rs = @network_component begin
        (p,d), 0 <--> X
    end
    @test_throws ArgumentError LatticeReactionSystem(rs, [tr], CartesianGrid((2,2)))
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
    for lattice in [small_1d_cartesian_grid, small_2d_cartesian_grid, small_3d_cartesian_grid,
                random_1d_masked_grid, random_2d_masked_grid, random_3d_masked_grid]
        lrs1 = LatticeReactionSystem(SIR_system, SIR_srs_1, lattice)
        lrs2 = LatticeReactionSystem(SIR_system, SIR_srs_1, lattice; diagonal_connections=true)
        @test num_edges(lrs1) == iterator_count(edge_iterator(lrs1))
        @test num_edges(lrs2) == iterator_count(edge_iterator(lrs2))
    end

    # Graph grids (cannot test diagonal connections).
    for lattice in [small_2d_graph_grid, small_3d_graph_grid, undirected_cycle, small_directed_cycle, unconnected_graph]
        lrs1 = LatticeReactionSystem(SIR_system, SIR_srs_1, lattice)
        @test num_edges(lrs1) == iterator_count(edge_iterator(lrs1))
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
    for grid in [small_1d_cartesian_grid, small_2d_cartesian_grid, small_3d_cartesian_grid,
                 small_1d_masked_grid, small_2d_masked_grid, small_3d_masked_grid,
                 random_1d_masked_grid, random_2d_masked_grid, random_3d_masked_grid]
        lrs = LatticeReactionSystem(rn, [tr], grid)
        flat_to_grid_idx = Catalyst.get_index_converters(lattice(lrs), num_verts(lrs))[1]
        edge_values = make_edge_p_values(lrs, make_edge_p_value)

        for e in edge_iterator(lrs)
            @test edge_values[e[1], e[2]] == make_edge_p_value(flat_to_grid_idx[e[1]], flat_to_grid_idx[e[2]])
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
    lrs = LatticeReactionSystem(rn, [tr], CartesianGrid(n))
    ps = [:D => make_directed_edge_values(lrs, (10.0, 0.0))]
    oprob = ODEProblem(lrs, u0, tspan, ps)
    @test isapprox(solve(oprob, Tsit5()).u[end][5], n, rtol=1e-6)

    # Checks the 2d case (both with 1d and 2d flow).
    lrs = LatticeReactionSystem(rn, [tr], CartesianGrid((n,n)))

    ps = [:D => make_directed_edge_values(lrs, (1.0, 0.0), (0.0, 0.0))]
    oprob = ODEProblem(lrs, u0, tspan, ps)
    @test all(isapprox.(solve(oprob, Tsit5()).u[end][5:5:25], n, rtol=1e-6))

    ps = [:D => make_directed_edge_values(lrs, (1.0, 0.0), (1.0, 0.0))]
    oprob = ODEProblem(lrs, u0, tspan, ps)
    @test isapprox(solve(oprob, Tsit5()).u[end][25], n^2, rtol=1e-6)

    # Checks the 3d case (both with 1d and 2d flow).
    lrs = LatticeReactionSystem(rn, [tr], CartesianGrid((n,n,n)))

    ps = [:D => make_directed_edge_values(lrs, (1.0, 0.0), (0.0, 0.0), (0.0, 0.0))]
    oprob = ODEProblem(lrs, u0, tspan, ps)
    @test all(isapprox.(solve(oprob, Tsit5()).u[end][5:5:125], n, rtol=1e-6))

    ps = [:D => make_directed_edge_values(lrs, (1.0, 0.0), (1.0, 0.0), (0.0, 0.0))]
    oprob = ODEProblem(lrs, u0, tspan, ps)
    @test all(isapprox.(solve(oprob, Tsit5()).u[end][25:25:125], n^2, rtol=1e-6))

    ps = [:D => make_directed_edge_values(lrs, (1.0, 0.0), (1.0, 0.0), (1.0, 0.0))]
    oprob = ODEProblem(lrs, u0, tspan, ps)
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
    lrs = LatticeReactionSystem(rn, [tr], path_graph(5))
    @test_throws Exception make_edge_p_values(lrs, make_edge_p_value,)
    @test_throws Exception make_directed_edge_values(lrs, (1.0, 0.0))

    # Wrong dimensions to `make_directed_edge_values`.
    lrs_1d = LatticeReactionSystem(rn, [tr], CartesianGrid(5))
    lrs_2d = LatticeReactionSystem(rn, [tr], fill(true,5,5))
    lrs_3d = LatticeReactionSystem(rn, [tr], CartesianGrid((5,5,5)))

    @test_throws Exception make_directed_edge_values(lrs_1d, (1.0, 0.0), (1.0, 0.0))
    @test_throws Exception make_directed_edge_values(lrs_1d, (1.0, 0.0), (1.0, 0.0), (1.0, 0.0))
    @test_throws Exception make_directed_edge_values(lrs_2d, (1.0, 0.0))
    @test_throws Exception make_directed_edge_values(lrs_2d, (1.0, 0.0), (1.0, 0.0), (1.0, 0.0))
    @test_throws Exception make_directed_edge_values(lrs_3d, (1.0, 0.0))
    @test_throws Exception make_directed_edge_values(lrs_3d, (1.0, 0.0), (1.0, 0.0))
end
