### Preparations ###

# Fetch packages.
using Catalyst, Graphs, OrdinaryDiffEq, Test

### Run Tests ###

# Checks that some grids, created using different approaches, generates the same spatial structures.
# Checks that some grids, created using different approaches, generates the same simulation output.
let 
    # Create LatticeReactionsSystems
    cartesian_grid = Graphs.grid([5, 5])
    regular_grid = fill(true, 5, 5)
    graph_grid = Graphs.grid([5, 5])

    cartesian_lrs = LatticeReactionSystem(brusselator_system, srs, cartesian_grid)
    regular_lrs = LatticeReactionSystem(brusselator_system, srs, regular_grid)
    graph_lrs = LatticeReactionSystem(brusselator_system, srs, graph_grid)

    # Check internal structures.
    @test cartesian_lrs.rs == regular_lrs.rs == graph_lrs.rs
    @test cartesian_lrs.spatial_reactions == regular_lrs.spatial_reactions == graph_lrs.spatial_reactions
    @test cartesian_lrs.num_verts == regular_lrs.num_verts == graph_lrs.num_verts
    @test cartesian_lrs.num_edges == regular_lrs.num_edges == graph_lrs.num_edges
    @test cartesian_lrs.num_species == regular_lrs.num_species == graph_lrs.num_species
    @test cartesian_lrs.spat_species == regular_lrs.spat_species == graph_lrs.spat_species
    @test cartesian_lrs.parameters == regular_lrs.parameters == graph_lrs.parameters
    @test cartesian_lrs.vertex_parameters == regular_lrs.vertex_parameters == graph_lrs.vertex_parameters
    @test cartesian_lrs.edge_parameters == regular_lrs.edge_parameters == graph_lrs.edge_parameters
    @test cartesian_lrs.directed_edges == regular_lrs.directed_edges == graph_lrs.directed_edges
    @test cartesian_lrs.edge_list == regular_lrs.edge_list == graph_lrs.edge_list

    # Checks that simulations yields the same output.
    u0 = [:X => rand_v_vals(graph_lrs.lattice, 10.0), :Y => 2.0]
    pV = [:A => 0.5 .+ rand_v_vals(graph_lrs.lattice, 0.5), :B => 4.0]
    pE= map(sp -> sp => 0.2, spatial_param_syms(lrs))

    cartesian_oprob = ODEProblem(cartesian_lrs, u0, (0.0, 100.0), (pV, pE))
    regular_oprob = ODEProblem(regular_lrs, u0, (0.0, 100.0), (pV, pE))
    graph_oprob = ODEProblem(graph_lrs, u0, (0.0, 100.0), (pV, pE))

    cartesian_sol = solve(cartesian_oprob, QNDF(); saveat=0.1)
    regular_sol = solve(regular_oprob, QNDF(); saveat=0.1)
    graph_sol = solvegraph(_oprob, QNDF(); saveat=0.1)

    @test cartesian_sol.u ≈ regular_sol.u ≈ graph_sol
end

# Checks that a regular grid with absent vertices generate the same output as corresponding graph.
let
    # Create LatticeReactionsSystems
    regular_grid = [true true true; true false true; true true true]
    graph_grid = cycle_graph(8)
    
    regular_lrs = LatticeReactionSystem(brusselator_system, srs, regular_grid)
    graph_lrs = LatticeReactionSystem(brusselator_system, srs, graph_grid)

    # Check internal structures.
    @test regular_lrs.num_verts == graph_lrs.num_verts
    @test regular_lrs.num_edges == graph_lrs.num_edges
    @test regular_lrs.edge_list == graph_lrs.edge_list

    # Checks that simulations yields the same output.
    u0 = [:X => rand_v_vals(graph_lrs.lattice, 10.0), :Y => 2.0]
    pV = [:A => 0.5 .+ rand_v_vals(graph_lrs.lattice, 0.5), :B => 4.0]
    pE= map(sp -> sp => 0.2, spatial_param_syms(lrs))

    regular_oprob = ODEProblem(regular_lrs, u0, (0.0, 100.0), (pV, pE))
    graph_oprob = ODEProblem(graph_lrs, u0, (0.0, 100.0), (pV, pE))

    regular_sol = solve(regular_oprob, QNDF(); saveat=0.1)
    graph_sol = solvegraph(_oprob, QNDF(); saveat=0.1)

    @test regular_sol.u ≈ graph_sol
end