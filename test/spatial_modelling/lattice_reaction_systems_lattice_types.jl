### Preparations ###

# Fetch packages.
using Catalyst, Graphs, OrdinaryDiffEqTsit5, OrdinaryDiffEqBDF, Test

# Fetch test networks.
include("../spatial_test_networks.jl")


### Run Tests ###

# Test errors when attempting to create networks with dimensions > 3.
let 
    @test_throws Exception LatticeReactionSystem(brusselator_system, brusselator_srs_1, CartesianGrid((5, 5, 5, 5)))
    @test_throws Exception LatticeReactionSystem(brusselator_system, brusselator_srs_1, fill(true, 5, 5, 5, 5))
end

# Checks that getter functions give the correct output.
let 
    # Create LatticeReactionsSystems.
    cartesian_1d_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_1d_cartesian_grid)
    cartesian_2d_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_2d_cartesian_grid)
    cartesian_3d_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_3d_cartesian_grid)
    masked_1d_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_1d_masked_grid)
    masked_2d_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_2d_masked_grid)
    masked_3d_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, small_3d_masked_grid)
    graph_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, medium_2d_graph_grid)

    # Test lattice type getters.
    @test has_cartesian_lattice(cartesian_2d_lrs)
    @test !has_cartesian_lattice(masked_2d_lrs)
    @test !has_cartesian_lattice(graph_lrs)

    @test !has_masked_lattice(cartesian_2d_lrs)
    @test has_masked_lattice(masked_2d_lrs)
    @test !has_masked_lattice(graph_lrs)

    @test has_grid_lattice(cartesian_2d_lrs)
    @test has_grid_lattice(masked_2d_lrs)
    @test !has_grid_lattice(graph_lrs)

    @test !has_graph_lattice(cartesian_2d_lrs)
    @test !has_graph_lattice(masked_2d_lrs)
    @test has_graph_lattice(graph_lrs)

    # Checks grid dimensions.
    @test grid_dims(cartesian_1d_lrs) == 1
    @test grid_dims(cartesian_2d_lrs) == 2
    @test grid_dims(cartesian_3d_lrs) == 3
    @test grid_dims(masked_1d_lrs) == 1
    @test grid_dims(masked_2d_lrs) == 2
    @test grid_dims(masked_3d_lrs) == 3
    @test_throws ArgumentError grid_dims(graph_lrs)

    # Checks grid sizes.
    @test grid_size(cartesian_1d_lrs) == (5,)
    @test grid_size(cartesian_2d_lrs) == (5,5)
    @test grid_size(cartesian_3d_lrs) == (5,5,5)
    @test grid_size(masked_1d_lrs) == (5,)
    @test grid_size(masked_2d_lrs) == (5,5)
    @test grid_size(masked_3d_lrs) == (5,5,5)
    @test_throws ArgumentError grid_size(graph_lrs)
end

# Checks grid dimensions for 2d and 3d grids where some dimension is equal to 1.
let
    # Creates LatticeReactionSystems
    cartesian_2d_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, CartesianGrid((1,5)))
    cartesian_3d_lrs_1 = LatticeReactionSystem(brusselator_system, brusselator_srs_1, CartesianGrid((1,5,5)))
    cartesian_3d_lrs_2 = LatticeReactionSystem(brusselator_system, brusselator_srs_1, CartesianGrid((1,1,5)))
    masked_2d_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, fill(true, 1, 5))
    masked_3d_lrs_1 = LatticeReactionSystem(brusselator_system, brusselator_srs_1, fill(true, 1, 5,5))
    masked_3d_lrs_2 = LatticeReactionSystem(brusselator_system, brusselator_srs_1, fill(true, 1, 1,5))

    # Check grid dimensions.
    @test grid_dims(cartesian_2d_lrs) == 2
    @test grid_dims(cartesian_3d_lrs_1) == 3
    @test grid_dims(cartesian_3d_lrs_2) == 3
    @test grid_dims(masked_2d_lrs) == 2
    @test grid_dims(masked_3d_lrs_1) == 3
    @test grid_dims(masked_3d_lrs_2) == 3
end

# Checks that some grids, created using different approaches, generates the same spatial structures.
# Checks that some grids, created using different approaches, generates the same simulation output.
let 
    # Create LatticeReactionsSystems.
    cartesian_grid = CartesianGrid((5, 5))
    masked_grid = fill(true, 5, 5)
    graph_grid = Graphs.grid([5, 5])

    cartesian_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, cartesian_grid)
    masked_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, masked_grid)
    graph_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, graph_grid)

    # Check internal structures.
    @test Catalyst.isequivalent(reactionsystem(cartesian_lrs), reactionsystem(masked_lrs))
    @test Catalyst.isequivalent(reactionsystem(masked_lrs), reactionsystem(graph_lrs))
    @test spatial_reactions(cartesian_lrs) == spatial_reactions(masked_lrs) == spatial_reactions(graph_lrs)
    @test num_verts(cartesian_lrs) == num_verts(masked_lrs) == num_verts(graph_lrs)
    @test num_edges(cartesian_lrs) == num_edges(masked_lrs) == num_edges(graph_lrs)
    @test num_species(cartesian_lrs) == num_species(masked_lrs) == num_species(graph_lrs)
    @test isequal(spatial_species(cartesian_lrs), spatial_species(masked_lrs)) 
    @test isequal(spatial_species(masked_lrs), spatial_species(graph_lrs)) 
    @test isequal(parameters(cartesian_lrs), parameters(masked_lrs))
    @test isequal(parameters(masked_lrs), parameters(graph_lrs))
    @test isequal(vertex_parameters(cartesian_lrs), vertex_parameters(masked_lrs))
    @test isequal(edge_parameters(masked_lrs), edge_parameters(graph_lrs))
    @test issetequal(edge_iterator(cartesian_lrs), edge_iterator(masked_lrs))
    @test issetequal(edge_iterator(masked_lrs), edge_iterator(graph_lrs))

    # Checks that simulations yields the same output.
    X_vals = rand(num_verts(cartesian_lrs))
    u0_cartesian = [:X => reshape(X_vals, 5, 5), :Y => 2.0]
    u0_masked = [:X => reshape(X_vals, 5, 5), :Y => 2.0]
    u0_graph = [:X => X_vals, :Y => 2.0]
    B_vals = rand(num_verts(cartesian_lrs))
    pV_cartesian = [:A => 0.5 .+ reshape(B_vals, 5, 5), :B => 4.0]
    pV_masked = [:A => 0.5 .+ reshape(B_vals, 5, 5), :B => 4.0]
    pV_graph = [:A => 0.5 .+ B_vals, :B => 4.0]
    pE = [:dX => 0.2]

    cartesian_oprob = ODEProblem(cartesian_lrs, u0_cartesian, (0.0, 100.0), [pV_cartesian; pE])
    masked_oprob = ODEProblem(masked_lrs, u0_masked, (0.0, 100.0), [pV_masked; pE])
    graph_oprob = ODEProblem(graph_lrs, u0_graph, (0.0, 100.0), [pV_graph; pE])

    cartesian_sol = solve(cartesian_oprob, FBDF(); saveat=0.1)
    masked_sol = solve(masked_oprob, FBDF(); saveat=0.1)
    graph_sol = solve(graph_oprob, FBDF(); saveat=0.1)

    @test cartesian_sol.u == masked_sol.u == graph_sol.u
end

# Checks that a regular grid with absent vertices generate the same output as corresponding graph.
let
    # Create LatticeReactionsSystems.
    masked_grid = [true true true; true false true; true true true]
    graph_grid = SimpleGraph(8)
    add_edge!(graph_grid, 1, 2); add_edge!(graph_grid, 2, 1);
    add_edge!(graph_grid, 2, 3); add_edge!(graph_grid, 3, 2);
    add_edge!(graph_grid, 3, 5); add_edge!(graph_grid, 5, 3);
    add_edge!(graph_grid, 5, 8); add_edge!(graph_grid, 8, 5);
    add_edge!(graph_grid, 8, 7); add_edge!(graph_grid, 7, 8);
    add_edge!(graph_grid, 7, 6); add_edge!(graph_grid, 6, 7);
    add_edge!(graph_grid, 6, 4); add_edge!(graph_grid, 4, 6);
    add_edge!(graph_grid, 4, 1); add_edge!(graph_grid, 1, 4);

    masked_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, masked_grid)
    graph_lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, graph_grid)
    
    # Check internal structures.
    @test num_verts(masked_lrs) == num_verts(graph_lrs)
    @test num_edges(masked_lrs) == num_edges(graph_lrs)
    @test issetequal(edge_iterator(masked_lrs), edge_iterator(graph_lrs))
    
    # Checks that simulations yields the same output.
    u0_masked_grid = [:X => [1. 4. 6.; 2. 0. 7.; 3. 5. 8.], :Y => 2.0]
    u0_graph_grid = [:X => [1., 2., 3., 4., 5., 6., 7., 8.], :Y => 2.0]
    pV_masked_grid = [:A => 0.5 .+ [1. 4. 6.; 2. 0. 7.; 3. 5. 8.], :B => 4.0]
    pV_graph_grid = [:A => 0.5 .+ [1., 2., 3., 4., 5., 6., 7., 8.], :B => 4.0]
    pE = [:dX => 0.2]
    
    base_oprob = ODEProblem(masked_lrs, u0_masked_grid, (0.0, 100.0), [pV_masked_grid; pE])
    base_osol = solve(base_oprob, QNDF(); saveat=0.1, abstol=1e-9, reltol=1e-9)

    for jac in [false, true], sparse in [false, true]
        masked_oprob = ODEProblem(masked_lrs, u0_masked_grid, (0.0, 100.0), [pV_masked_grid; pE]; jac, sparse)
        graph_oprob = ODEProblem(graph_lrs, u0_graph_grid, (0.0, 100.0), [pV_graph_grid; pE]; jac, sparse)
        masked_sol = solve(masked_oprob, QNDF(); saveat=0.1, abstol=1e-9, reltol=1e-9)
        graph_sol = solve(graph_oprob, QNDF(); saveat=0.1, abstol=1e-9, reltol=1e-9)
        @test base_osol ≈ masked_sol atol = 1e-6 rtol = 1e-6
        @test base_osol ≈  graph_sol atol = 1e-6 rtol = 1e-6
        @test masked_sol ≈ graph_sol atol = 1e-6 rtol = 1e-6
    end
end

# For a system which is a single line of vertices: (O-O-O-O-X-O-O-O), ensures that different simulations
# approach yield the same result. Checks for both masked and Cartesian grid. For both, simulates where
# initial conditions/vertex parameters are either a vector of the same length as the number of vertices (7),
# Or as the grid. Here, we try grid sizes (n), (1,n), and (1,n,1) (so the same grid, but in 1d, 2d, and 3d).
# For the Cartesian grid, we cannot represent the gap, so we make simulations both for length-4 and
# length-3 grids.
let
    # Declares the initial condition/parameter values.
    S_vals = [500.0, 600.0, 700.0, 800.0, 0.0, 900.0, 1000.0, 1100.0]
    I_val = 1.0
    R_val = 1.0
    α_vals = [0.1, 0.11, 0.12, 0.13, 0.0, 0.14, 0.15, 0.16]
    β_val = 0.01
    dS_val = 0.05
    SIR_p = [:α => 0.1 / 1000, :β => 0.01]

    # Declares the grids (1d, 2d, and 3d). For each dimension, there are a 2 Cartesian grids (length 4 and 3).
    cart_grid_1d_1 = CartesianGrid(4)
    cart_grid_1d_2 = CartesianGrid(3)
    cart_grid_2d_1 = CartesianGrid((4,1))
    cart_grid_2d_2 = CartesianGrid((3,1))
    cart_grid_3d_1 = CartesianGrid((1,4,1))
    cart_grid_3d_2 = CartesianGrid((1,3,1))

    masked_grid_1d = [true, true, true, true, false, true, true, true]
    masked_grid_2d = reshape(masked_grid_1d,8,1)
    masked_grid_3d = reshape(masked_grid_1d,1,8,1)

    # Creaets a base solution to which we will compare all simulations.
    lrs_base = LatticeReactionSystem(SIR_system, SIR_srs_1, masked_grid_1d)
    oprob_base = ODEProblem(lrs_base, [:S => S_vals, :I => I_val, :R => R_val], (0.0, 100.0), [:α => α_vals, :β => β_val, :dS => dS_val])
    sol_base = solve(oprob_base, Tsit5(); saveat = 1.0, abstol = 1e-9, reltol = 1e-9)

    # Checks simulations for the masked grid (covering all 7 vertices, with a gap in the middle).
    for grid in [masked_grid_1d, masked_grid_2d, masked_grid_3d]
        # Checks where the values are vectors of length equal to the number of vertices.
        lrs = LatticeReactionSystem(SIR_system, SIR_srs_1, grid)
        u0 = [:S => [S_vals[1:4]; S_vals[6:8]], :I => I_val, :R => R_val]
        ps = [:α => [α_vals[1:4]; α_vals[6:8]], :β => β_val, :dS => dS_val]
        oprob = ODEProblem(lrs, u0, (0.0, 100.0), ps)
        sol = solve(oprob, Tsit5(); saveat = 1.0, abstol = 1e-9, reltol = 1e-9)
        @test sol ≈ sol_base

        # Checks where the values are arrays of size equal to the grid.
        u0 = [:S => reshape(S_vals, grid_size(lrs)), :I => I_val, :R => R_val]
        ps = [:α => reshape(α_vals, grid_size(lrs)), :β => β_val, :dS => dS_val]
        oprob = ODEProblem(lrs, u0, (0.0, 100.0), ps)
        sol = solve(oprob, Tsit5(); saveat = 1.0, abstol = 1e-9, reltol = 1e-9)
        @test sol ≈ sol_base
    end

    # Checks simulations for the first Cartesian grids (covering vertices 1 to 4).
    for grid in [cart_grid_1d_1, cart_grid_2d_1, cart_grid_3d_1]
        # Checks where the values are vectors of length equal to the number of vertices.
        lrs = LatticeReactionSystem(SIR_system, SIR_srs_1, grid)
        u0 = [:S => S_vals[1:4], :I => I_val, :R => R_val]
        ps = [:α => α_vals[1:4], :β => β_val, :dS => dS_val]
        oprob = ODEProblem(lrs, u0, (0.0, 100.0), ps)
        sol = solve(oprob, Tsit5(); saveat = 1.0, abstol = 1e-9, reltol = 1e-9)
        @test hcat(sol.u...) ≈ sol_base[1:12,:]

        # Checks where the values are arrays of size equal to the grid.
        u0 = [:S => reshape(S_vals[1:4], grid_size(lrs)), :I => I_val, :R => R_val]
        ps = [:α => reshape(α_vals[1:4], grid_size(lrs)), :β => β_val, :dS => dS_val]
        oprob = ODEProblem(lrs, u0, (0.0, 100.0), ps)
        sol = solve(oprob, Tsit5(); saveat = 1.0, abstol = 1e-9, reltol = 1e-9)
        @test hcat(sol.u...) ≈ sol_base[1:12,:]
    end

    # Checks simulations for the second Cartesian grids (covering vertices 6 to 8).
    for grid in [cart_grid_1d_2, cart_grid_2d_2, cart_grid_3d_2]
        # Checks where the values are vectors of length equal to the number of vertices.
        lrs = LatticeReactionSystem(SIR_system, SIR_srs_1, grid)
        u0 = [:S => S_vals[6:8], :I => I_val, :R => R_val]
        ps = [:α => α_vals[6:8], :β => β_val, :dS => dS_val]
        oprob = ODEProblem(lrs, u0, (0.0, 100.0), ps)
        sol = solve(oprob, Tsit5(); saveat = 1.0, abstol = 1e-9, reltol = 1e-9)
        @test hcat(sol.u...) ≈ sol_base[13:end,:]

        # Checks where the values are arrays of size equal to the grid.
        u0 = [:S => reshape(S_vals[6:8], grid_size(lrs)), :I => I_val, :R => R_val]
        ps = [:α => reshape(α_vals[6:8], grid_size(lrs)), :β => β_val, :dS => dS_val]
        oprob = ODEProblem(lrs, u0, (0.0, 100.0), ps)
        sol = solve(oprob, Tsit5(); saveat = 1.0, abstol = 1e-9, reltol = 1e-9)
        @test hcat(sol.u...) ≈ sol_base[13:end,:]
    end
end
