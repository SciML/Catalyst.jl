### Preparations ###

# Fetch packages.
using Catalyst, Graphs, OrdinaryDiffEq, Test

### Run Tests ###

# Test errors when attempting to create networks with dimension > 3.
let 
    @test_throws Exception LatticeReactionSystem(brusselator_system, srs, CartesianGrid((5, 5, 5, 5)))
    @test_throws Exception LatticeReactionSystem(brusselator_system, srs, fill(true, 5, 5, 5, 5))
end

# Checks that getter functions give the correct output.
let 
    # Create LatticeReactionsSystems.
    cartesian_1d_lrs = LatticeReactionSystem(brusselator_system, srs, small_1d_cartesian_grid)
    cartesian_2d_lrs = LatticeReactionSystem(brusselator_system, srs, small_2d_cartesian_grid)
    cartesian_3d_lrs = LatticeReactionSystem(brusselator_system, srs, small_3d_cartesian_grid)
    masked_1d_lrs = LatticeReactionSystem(brusselator_system, srs, small_1d_masked_grid)
    masked_2d_lrs = LatticeReactionSystem(brusselator_system, srs, small_2d_masked_grid)
    masked_3d_lrs = LatticeReactionSystem(brusselator_system, srs, small_3d_masked_grid)
    graph_lrs = LatticeReactionSystem(brusselator_system, srs, small_2d_grid)

    # Test lattice type getters.
    @test has_cartesian_lattice(cartesian_2d_lrs)
    @test has_cartesian_lattice(masked_2d_lrs)
    @test has_cartesian_lattice(graph_lrs)

    @test !has_masked_lattice(cartesian_2d_lrs)
    @test !has_masked_lattice(masked_2d_lrs)
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

# Checks that some grids, created using different approaches, generates the same spatial structures.
# Checks that some grids, created using different approaches, generates the same simulation output.
let 
    # Create LatticeReactionsSystems.
    cartesian_grid = Graphs.grid([5, 5])
    masked_grid = fill(true, 5, 5)
    graph_grid = Graphs.grid([5, 5])

    cartesian_lrs = LatticeReactionSystem(brusselator_system, srs, cartesian_grid)
    masked_lrs = LatticeReactionSystem(brusselator_system, srs, masked_grid)
    graph_lrs = LatticeReactionSystem(brusselator_system, srs, graph_grid)

    # Check internal structures.
    @test cartesian_lrs.rs == masked_lrs.rs == graph_lrs.rs
    @test cartesian_lrs.spatial_reactions == masked_lrs.spatial_reactions == graph_lrs.spatial_reactions
    @test cartesian_lrs.num_verts == masked_lrs.num_verts == graph_lrs.num_verts
    @test cartesian_lrs.num_edges == masked_lrs.num_edges == graph_lrs.num_edges
    @test cartesian_lrs.num_species == masked_lrs.num_species == graph_lrs.num_species
    @test cartesian_lrs.spat_species == masked_lrs.spat_species == graph_lrs.spat_species
    @test cartesian_lrs.parameters == masked_lrs.parameters == graph_lrs.parameters
    @test cartesian_lrs.vertex_parameters == masked_lrs.vertex_parameters == graph_lrs.vertex_parameters
    @test cartesian_lrs.edge_parameters == masked_lrs.edge_parameters == graph_lrs.edge_parameters
    @test cartesian_lrs.directed_edges == masked_lrs.directed_edges == graph_lrs.directed_edges
    @test cartesian_lrs.edge_list == masked_lrs.edge_list == graph_lrs.edge_list

    # Checks that simulations yields the same output.
    u0 = [:X => rand_v_vals(graph_lrs.lattice, 10.0), :Y => 2.0]
    pV = [:A => 0.5 .+ rand_v_vals(graph_lrs.lattice, 0.5), :B => 4.0]
    pE= map(sp -> sp => 0.2, spatial_param_syms(lrs))

    cartesian_oprob = ODEProblem(cartesian_lrs, u0, (0.0, 100.0), (pV, pE))
    masked_oprob = ODEProblem(masked_lrs, u0, (0.0, 100.0), (pV, pE))
    graph_oprob = ODEProblem(graph_lrs, u0, (0.0, 100.0), (pV, pE))

    cartesian_sol = solve(cartesian_oprob, QNDF(); saveat=0.1)
    masked_sol = solve(masked_oprob, QNDF(); saveat=0.1)
    graph_sol = solvegraph(_oprob, QNDF(); saveat=0.1)

    @test cartesian_sol.u ≈ masked_sol.u ≈ graph_sol
end

# Checks that a regular grid with absent vertices generate the same output as corresponding graph.
let
    # Create LatticeReactionsSystems.
    masked_grid = [true true true; true false true; true true true]
    graph_grid = cycle_graph(8)
    
    masked_lrs = LatticeReactionSystem(brusselator_system, srs, masked_grid)
    graph_lrs = LatticeReactionSystem(brusselator_system, srs, graph_grid)

    # Check internal structures.
    @test masked_lrs.num_verts == graph_lrs.num_verts
    @test masked_lrs.num_edges == graph_lrs.num_edges
    @test masked_lrs.edge_list == graph_lrs.edge_list

    # Checks that simulations yields the same output.
    u0 = [:X => rand_v_vals(graph_lrs.lattice, 10.0), :Y => 2.0]
    pV = [:A => 0.5 .+ rand_v_vals(graph_lrs.lattice, 0.5), :B => 4.0]
    pE= map(sp -> sp => 0.2, spatial_param_syms(lrs))

    masked_oprob = ODEProblem(masked_lrs, u0, (0.0, 100.0), (pV, pE))
    graph_oprob = ODEProblem(graph_lrs, u0, (0.0, 100.0), (pV, pE))

    masked_sol = solve(masked_oprob, QNDF(); saveat=0.1)
    graph_sol = solvegraph(_oprob, QNDF(); saveat=0.1)

    @test masked_sol.u ≈ graph_sol
end