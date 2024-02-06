### Preparations ###

# Fetch packages.
using Catalyst, Graphs, OrdinaryDiffEq, Test

# Fetch test networks.
include("../spatial_test_networks.jl")

### Run Tests ###

# Test errors when attempting to create networks with dimension > 3.
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
    @test_throws Exception grid_dims(graph_lrs)
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
    @test cartesian_lrs.rs == masked_lrs.rs == graph_lrs.rs
    @test cartesian_lrs.spatial_reactions == masked_lrs.spatial_reactions == graph_lrs.spatial_reactions
    @test cartesian_lrs.num_verts == masked_lrs.num_verts == graph_lrs.num_verts
    @test cartesian_lrs.num_edges == masked_lrs.num_edges == graph_lrs.num_edges
    @test cartesian_lrs.num_species == masked_lrs.num_species == graph_lrs.num_species
    @test isequal(cartesian_lrs.spat_species, masked_lrs.spat_species) 
    @test isequal(masked_lrs.spat_species, graph_lrs.spat_species) 
    @test isequal(cartesian_lrs.parameters, masked_lrs.parameters)
    @test isequal(masked_lrs.parameters, graph_lrs.parameters)
    @test isequal(cartesian_lrs.vertex_parameters, masked_lrs.vertex_parameters)
    @test isequal(masked_lrs.edge_parameters, graph_lrs.edge_parameters)
    @test issetequal(cartesian_lrs.edge_iterator, masked_lrs.edge_iterator)
    @test issetequal(masked_lrs.edge_iterator, graph_lrs.edge_iterator)

    # Checks that simulations yields the same output.
    X_vals = rand(cartesian_lrs.num_verts)
    u0_cartesian = [:X => reshape(X_vals, 5, 5), :Y => 2.0]
    u0_masked = [:X => reshape(X_vals, 5, 5), :Y => 2.0]
    u0_graph = [:X => X_vals, :Y => 2.0]
    B_vals = rand(cartesian_lrs.num_verts)
    pV_cartesian = [:A => 0.5 .+ reshape(B_vals, 5, 5), :B => 4.0]
    pV_masked = [:A => 0.5 .+ reshape(B_vals, 5, 5), :B => 4.0]
    pV_graph = [:A => 0.5 .+ B_vals, :B => 4.0]
    pE = [:dX => 0.2]

    cartesian_oprob = ODEProblem(cartesian_lrs, u0_cartesian, (0.0, 100.0), [pV_cartesian; pE])
    masked_oprob = ODEProblem(masked_lrs, u0_masked, (0.0, 100.0), [pV_masked; pE])
    graph_oprob = ODEProblem(graph_lrs, u0_graph, (0.0, 100.0), [pV_graph; pE])

    cartesian_sol = solve(cartesian_oprob, QNDF(); saveat=0.1)
    masked_sol = solve(masked_oprob, QNDF(); saveat=0.1)
    graph_sol = solve(graph_oprob, QNDF(); saveat=0.1)

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
    @test masked_lrs.num_verts == graph_lrs.num_verts
    @test masked_lrs.num_edges == graph_lrs.num_edges
    @test issetequal(masked_lrs.edge_iterator, graph_lrs.edge_iterator)
    
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
        @test base_osol ≈ masked_sol ≈ graph_sol
    end
end