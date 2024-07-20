### Preparations ###

# Fetch packages.
using Catalyst, Graphs, JumpProcesses, OrdinaryDiffEq, SparseArrays, Test

# Fetch test networks.
include("../spatial_test_networks.jl")

### Problem & Integrator Interfacing Function Tests ###

# Checks `lat_getu` for ODE and Jump problem and integrators.
# Checks `lat_setu!` for ODE and Jump problem and integrators.
# Checks for all types of lattices.
# Checks for symbol and symbolic variables input.
let
    # Declares various types of lattices and corresponding initial values of `X`.
    lattice_cartesian = CartesianGrid((2,2,2))
    lattice_masked = [true true; false true]
    lattice_graph = cycle_graph(5)
    val0_cartesian = fill(1.0, 2, 2, 2)
    val0_masked = sparse([1.0 2.0; 0.0 3.0])
    val0_graph = [1.0, 2.0, 3.0, 4.0, 5.0]

    # Unpacks the `X`  and `Y` symbolic variable (so that indexing using it can be tested).
    @unpack X, Y = brusselator_system

    # Loops through all alternative lattices and values. Checks that `lat_getu` works in all cases.
    for (lattice, val0) in zip([lattice_cartesian, lattice_masked, lattice_graph], [val0_cartesian, val0_masked, val0_graph])
        # Prepares various problems and integrators. Uses `deepcopy` to ensure there is no cross-talk
        # between the different u vectors as they get updated.
        lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, lattice)
        u0 = [:X => val0, :Y => 0.5]
        ps = [:A => 1.0, :B => 2.0, :dX => 0.1]
        oprob = ODEProblem(lrs, deepcopy(u0), (0.0, 1.0), ps)
        dprob = DiscreteProblem(lrs, deepcopy(u0), (0.0, 1.0), ps)
        jprob = JumpProblem(lrs, dprob, NSM())
        oint = init(deepcopy(oprob), Tsit5())
        jint = init(deepcopy(jprob), SSAStepper())
        
        # Check that `lat_getu` retrieves the correct values.
        @test lat_getu(oprob, :X, lrs) == lat_getu(oprob, X, lrs) == lat_getu(oprob, brusselator_system.X, lrs) == val0
        @test lat_getu(oint, :X, lrs) == lat_getu(oint, X, lrs) == lat_getu(oint, brusselator_system.X, lrs) == val0
        @test lat_getu(jprob, :X, lrs) == lat_getu(jprob, X, lrs) == lat_getu(jprob, brusselator_system.X, lrs) == val0
        @test lat_getu(jint, :X, lrs) == lat_getu(jint, X, lrs) == lat_getu(jint, brusselator_system.X, lrs) == val0
        
        # Updates Y and checks its content.
        lat_setu!(oprob, :Y, lrs, val0)
        @test lat_getu(oprob, :Y, lrs) == lat_getu(oprob, Y, lrs) == lat_getu(oprob, brusselator_system.Y, lrs) == val0
        lat_setu!(oint, :Y, lrs, val0)
        @test lat_getu(oint, :Y, lrs) == lat_getu(oint, Y, lrs) == lat_getu(oint, brusselator_system.Y, lrs) == val0
        lat_setu!(jprob, :Y, lrs, val0)
        @test lat_getu(jprob, :Y, lrs) == lat_getu(jprob, Y, lrs) == lat_getu(jprob, brusselator_system.Y, lrs) == val0
        lat_setu!(jint, :Y, lrs, val0)
        @test lat_getu(jint, :Y, lrs) == lat_getu(jint, Y, lrs) == lat_getu(jint, brusselator_system.Y, lrs) == val0

        # Tries where we change a spatially non-uniform variable to spatially uniform.
        lat_setu!(oprob, X, lrs, 0.0)
        @test all(isequal(0.0), lat_getu(oprob, X, lrs))
        lat_setu!(oint, X, lrs, 0.0)
        @test all(isequal(0.0), lat_getu(oint, X, lrs))
        lat_setu!(jprob, X, lrs, 0.0)
        @test all(isequal(0.0), lat_getu(jprob, X, lrs))
        lat_setu!(jint, X, lrs, 0.0)
        @test all(isequal(0.0), lat_getu(jint, X, lrs))
    end
end

# Checks `lat_getp` for ODEproblem and integrators.
# Checks `lat_setp!` for ODE problem and integrators.
# Checks for all types of lattices.
# Checks for symbol and symbolic variables input.
let
    # Declares various types of lattices and corresponding initial values of `A`.
    lattice_cartesian = CartesianGrid((2,2,2))
    lattice_masked = [true true; false true]
    lattice_graph = cycle_graph(5)
    val0_cartesian = fill(1.0, 2, 2, 2)
    val0_masked = sparse([1.0 2.0; 0.0 3.0])
    val0_graph = [1.0, 2.0, 3.0, 4.0, 5.0]

    # Unpacks the `A` and `B` symbolic variable (so that indexing using it can be tested).
    @unpack A, B = brusselator_system

    # Loops through all alternative lattices and values. Checks that `lat_getp` works in all cases.
    for (lattice, val0) in zip([lattice_cartesian, lattice_masked, lattice_graph], [val0_cartesian, val0_masked, val0_graph])
        # Prepares various problems and integrators. Uses `deepcopy` to ensure there is no cross-talk
        # between the different p vectors as they get updated.
        lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, lattice)
        u0 = [:X => 1.0, :Y => 0.5]
        ps = [:A => val0, :B => 2.0, :dX => 0.1]
        oprob = ODEProblem(lrs, u0, (0.0, 1.0), deepcopy(ps))
        oint = init(deepcopy(oprob), Tsit5())
        
        # Check that `lat_getp` retrieves the correct values.
        @test lat_getp(oprob, :A, lrs) == lat_getp(oprob, A, lrs) == lat_getp(oprob, brusselator_system.A, lrs) == val0
        @test lat_getp(oint, :A, lrs) == lat_getp(oint, A, lrs) == lat_getp(oint, brusselator_system.A, lrs) == val0
        
        # Updates Y and checks its content.
        lat_setp!(oprob, :B, lrs, val0)
        @test lat_getp(oprob, :B, lrs) == lat_getp(oprob, B, lrs) == lat_getp(oprob, brusselator_system.B, lrs) == val0
        lat_setp!(oint, :B, lrs, val0)
        @test lat_getp(oint, :B, lrs) == lat_getp(oint, B, lrs) == lat_getp(oint, brusselator_system.B, lrs) == val0

        # Tries where we change a spatially non-uniform variable to spatially uniform.
        lat_setp!(oprob, A, lrs, 0.0)
        @test all(isequal(0.0), lat_getp(oprob, A, lrs))
        lat_setp!(oint, A, lrs, 0.0)
        @test all(isequal(0.0), lat_getp(oint, A, lrs))
    end
end

# Checks that `lat_getp` and `lat_setp!` generates errors when applied to `JumpProblem`s and their integrators.
let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, very_small_1d_cartesian_grid)
    u0 = [:X => 1, :Y => 0]
    ps = [:A => 3.0, :B => 2.0, :dX => 0.1]
    dprob = DiscreteProblem(lrs, u0, (0.0, 1.0), ps)
    jprob = JumpProblem(lrs, dprob, NSM())
    jint = init(jprob, SSAStepper())

    @test_throws Exception lat_getp(jprob, :A, lrs)
    @test_throws Exception lat_getp(jprob, :A, lrs)
    @test_throws Exception lat_setp!(jint, :A, lrs, 0.0)
    @test_throws Exception lat_setp!(jint, :A, lrs, 0.0)
end

# Checks that `lat_getp` and `lat_setp!` generates errors when applied to edge parameters.
let
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_1, very_small_1d_cartesian_grid)
    u0 = [:X => 1.0, :Y => 0.0]
    ps = [:A => 3.0, :B => 2.0, :dX => 0.1]
    oprob = ODEProblem(lrs, u0, (0.0, 1.0), ps)
    oint = init(deepcopy(oprob), Tsit5())

    @test_throws ArgumentError lat_getp(oprob, :dX, lrs)
    @test_throws ArgumentError lat_setp!(oprob, :dX, lrs, 0.0)
end

### Simulation `lat_getu` Tests ###

# Basic test. For simulations without change in system, check that the solution corresponds to known
# initial condition throughout the solution. 
# Checks using both `t` sampling` and normal time step sampling.
# Checks for both ODE and jump simulations.
# Checks for all lattice types.
let 
    # Prepare `LatticeReactionSystem`s.
    rs = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    tr = @transport_reaction D X1
    lrs1 = LatticeReactionSystem(rs, [tr], CartesianGrid((2,)))
    lrs2 = LatticeReactionSystem(rs, [tr], CartesianGrid((2,3)))
    lrs3 = LatticeReactionSystem(rs, [tr], CartesianGrid((2,3,2)))
    lrs4 = LatticeReactionSystem(rs, [tr], [true, true, false, true])
    lrs5 = LatticeReactionSystem(rs, [tr], [true false; true true])
    lrs6 = LatticeReactionSystem(rs, [tr], cycle_graph(4))

    # Create problem inputs.
    u0_1 = Dict([:X1 => 0, :X2 => [1, 2]])
    u0_2 = Dict([:X1 => 0, :X2 => [1 2 3; 4 5 6]])
    u0_3 = Dict([:X1 => 0, :X2 => fill(1, 2, 3, 2)])
    u0_4 = Dict([:X1 => 0, :X2 => sparse([1, 2, 0, 3])])
    u0_5 = Dict([:X1 => 0, :X2 => sparse([1 0; 2 3])])
    u0_6 = Dict([:X1 => 0, :X2 => [1, 2, 3, 4]])
    tspan = (0.0, 1.0)
    ps = [:k1 => 0.0, :k2 => 0.0, :D => 0.0]

    # Loops through all lattice cases and check that they are correct.
    for (u0,lrs) in zip([u0_1, u0_2, u0_3, u0_4, u0_5, u0_6], [lrs1, lrs2, lrs3, lrs4, lrs5, lrs6])
        # Simulates ODE version and checks `lat_getu` on its solution.
        oprob = ODEProblem(lrs, u0, tspan, ps)
        osol = solve(oprob, Tsit5(), saveat = 0.5)
        @test lat_getu(osol, :X1, lrs) == lat_getu(osol, :X1, lrs; t = 0.0:0.5:1.0)
        @test all(all(val == Float64(u0[:X1]) for val in vals) for vals in lat_getu(osol, :X1, lrs))
        @test lat_getu(osol, :X2, lrs) == lat_getu(osol, :X2, lrs; t = 0.0:0.5:1.0) == fill(u0[:X2], 3)

        # Simulates jump version and checks `lat_getu` on its solution.
        dprob = DiscreteProblem(lrs, u0, tspan, ps)
        jprob = JumpProblem(lrs, dprob, NSM())
        jsol = solve(jprob, SSAStepper(), saveat = 0.5)
        @test lat_getu(jsol, :X1, lrs) == lat_getu(jsol, :X1, lrs; t = 0.0:0.5:1.0)
        @test all(all(val == Float64(u0[:X1]) for val in vals) for vals in lat_getu(jsol, :X1, lrs))
        @test lat_getu(jsol, :X2, lrs) == lat_getu(jsol, :X2, lrs; t = 0.0:0.5:1.0) == fill(u0[:X2], 3)
    end
end

# Checks on simulations where the system changes in time.
# Checks that a solution has correct initial condition and end point (steady state).
# Checks that solution is monotonously increasing/decreasing (it should be for this problem).
let
    # Prepare `LatticeReactionSystem`s.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    lrs = LatticeReactionSystem(rs, [tr], CartesianGrid((2,)))

    # Prepares a corresponding ODEProblem.
    u0 = [:X => [1.0, 3.0]]
    tspan = (0.0, 50.0)
    ps = [:p => 2.0, :d => 1.0, :D => 0.01]
    oprob = ODEProblem(lrs, u0, tspan, ps)

    # Simulates the ODE. Checks that the start/end points are correct.
    # Check that the first vertex is monotonously increasing in values, and that the second one is 
    # monotonously decreasing. The non evenly spaced `saveat` is so that non-monotonicity is
    # not produced due to numeric errors.
    saveat = [0.0, 1.0, 5.0, 10.0, 50.0]
    sol = solve(oprob, Vern7(); abstol = 1e-8, reltol = 1e-8)
    vals = lat_getu(sol, :X, lrs)
    @test vals[1] == [1.0, 3.0]
    @test vals[end] ≈ [2.0, 2.0]
    for i = 1:(length(saveat) - 1)
        @test vals[i][1] < vals[i + 1][1]
        @test vals[i][2] > vals[i + 1][2]
    end
end

# Checks interpolation when sampling at time point. Check that values at `t` is in between the 
# sample points. Does so by checking that in simulation which is monotonously decreasing/increasing.
let
    # Prepare `LatticeReactionSystem`s.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    lrs = LatticeReactionSystem(rs, [tr], CartesianGrid((2,)))

    # Solved a corresponding ODEProblem.
    u0 = [:X => [1.0, 3.0]]
    tspan = (0.0, 1.0)
    ps = [:p => 2.0, :d => 1.0, :D => 0.0]
    oprob = ODEProblem(lrs, u0, tspan, ps)

    # Solves and check the interpolation of t.
    sol = solve(oprob, Tsit5(); saveat = 1.0)
    t5_vals = lat_getu(sol, :X, lrs; t = [0.5])[1]
    @test sol.u[1][1] < t5_vals[1] < sol.u[2][1]
    @test sol.u[1][2] > t5_vals[2] > sol.u[2][2]
end

# Checks that attempting to sample `t` outside tspan range yields an error.
let
    # Prepare `LatticeReactionSystem`s.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    lrs = LatticeReactionSystem(rs, [tr], CartesianGrid((2,)))

    # Solved a corresponding ODEProblem.
    u0 = [:X => 1.0]
    tspan = (1.0, 2.0)
    ps = [:p => 2.0, :d => 1.0, :D => 1.0]
    oprob = ODEProblem(lrs, u0, tspan, ps)

    # Solves and check the interpolation of t.
    sol = solve(oprob, Tsit5(); saveat = 1.0)
    @test_throws Exception lat_getu(sol, :X, lrs; t = [0.0])
    @test_throws Exception lat_getu(sol, :X, lrs; t = [3.0])
end

# Checks that attempting to sample `t` outside tspan range yields an error.
let
    # Prepare `LatticeReactionSystem`s.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    lrs = LatticeReactionSystem(rs, [tr], CartesianGrid((2,)))

    # Solved a corresponding ODEProblem.
    u0 = [:X => 1.0]
    tspan = (1.0, 2.0)
    ps = [:p => 2.0, :d => 1.0, :D => 1.0]
    oprob = ODEProblem(lrs, u0, tspan, ps)

    # Solves and check the interpolation of t.
    sol = solve(oprob, Tsit5(); saveat = 1.0)
    @test_throws Exception lat_getu(sol, :X, lrs; t = [0.0])
    @test_throws Exception lat_getu(sol, :X, lrs; t = [3.0])
end

# Checks that applying `lat_getu` to a 3d masked lattice yields an error.
let
    # Prepare `LatticeReactionSystem`s.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    lrs = LatticeReactionSystem(rs, [tr], rand([false, true], 2, 3, 4))

    # Solved a corresponding ODEProblem.
    u0 = [:X => 1.0]
    tspan = (1.0, 2.0)
    ps = [:p => 2.0, :d => 1.0, :D => 1.0]
    oprob = ODEProblem(lrs, u0, tspan, ps)

    # Solves and check the interpolation of t.
    sol = solve(oprob, Tsit5(); saveat = 1.0)
    @test_throws Exception lat_getu(sol, :X, lrs)
end

# Checks that `lat_getu` works for all types of symbols.
let
    t = default_t()
    @species X(t)
    @parameters d
    @named rs = ReactionSystem([Reaction(d, [X], [])], t)
    rs = complete(rs)
    tr = @transport_reaction D X
    lrs = LatticeReactionSystem(rs, [tr], CartesianGrid(2,))

    # Solved a corresponding ODEProblem.
    u0 = [:X => 1.0]
    tspan = (0.0, 1.0)
    ps = [:d => 1.0, :D => 0.1]
    oprob = ODEProblem(lrs, u0, tspan, ps)

    # Solves and check the interpolation of t.
    sol = solve(oprob, Tsit5(); saveat = 1.0)
    @test lat_getu(sol, X, lrs) == lat_getu(sol, rs.X, lrs) == lat_getu(sol, :X, lrs)
    @test lat_getu(sol, X, lrs; t = 0.0:0.5:1.0) == lat_getu(sol, rs.X, lrs; t = 0.0:0.5:1.0) == lat_getu(sol, :X, lrs; t = 0.0:0.5:1.0)
end

### ODEProblem & Integrator Rebuilding ###

# Checks that the `rebuild_lat_internals!` function is correctly applied to an ODEProblem.
let
    # Creates a Brusselator `LatticeReactionSystem`.
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_2, very_small_2d_cartesian_grid)

    # Checks for all combinations of Jacobian and sparsity.
    for jac in [false, true], sparse in [false, true]
        # Creates an initial ODEProblem.
        u0 = [:X => 1.0, :Y => [1.0 2.0; 3.0 4.0]]
        dY_vals = spzeros(4,4)
        dY_vals[1,2] = 0.1; dY_vals[2,1] = 0.1; 
        dY_vals[1,3] = 0.2; dY_vals[3,1] = 0.2; 
        dY_vals[2,4] = 0.3; dY_vals[4,2] = 0.3; 
        dY_vals[3,4] = 0.4; dY_vals[4,3] = 0.4; 
        ps = [:A => 1.0, :B => [4.0 5.0; 6.0 7.0], :dX => 0.1, :dY => dY_vals]
        oprob_1 = ODEProblem(lrs, u0, (0.0, 10.0), ps; jac, sparse)

        # Creates an alternative version of the ODEProblem.
        dX_vals = spzeros(4,4)
        dX_vals[1,2] = 0.01; dX_vals[2,1] = 0.01; 
        dX_vals[1,3] = 0.02; dX_vals[3,1] = 0.02; 
        dX_vals[2,4] = 0.03; dX_vals[4,2] = 0.03; 
        dX_vals[3,4] = 0.04; dX_vals[4,3] = 0.04; 
        ps = [:A => [1.1 1.2; 1.3 1.4], :B => 5.0, :dX => dX_vals, :dY => 0.01]
        oprob_2 = ODEProblem(lrs, u0, (0.0, 10.0), ps; jac, sparse)

        # Modifies the initial ODEProblem to be identical to the new one.
        lat_setp!(oprob_1, :A, lrs, [1.1 1.2; 1.3 1.4])
        lat_setp!(oprob_1, :B, lrs, 5.0)
        oprob_1.ps[:dX] = dX_vals
        oprob_1.ps[:dY] = [0.01]
        rebuild_lat_internals!(oprob_1)

        # Checks that simulations of the two `ODEProblem`s are identical.
        @test solve(oprob_1, Rodas5P()) ≈ solve(oprob_2, Rodas5P())
    end
end

# Checks that the `rebuild_lat_internals!` function is correctly applied to an integrator.
# Does through by applying it within a callback, and compare to simulations without callback.
# To keep test faster, only check for `jac = sparse = true` only.
let
    # Prepares problem inputs.
    lrs = LatticeReactionSystem(brusselator_system, brusselator_srs_2, very_small_2d_cartesian_grid)
    u0 = [:X => 1.0, :Y => [1.0 2.0; 3.0 4.0]]
    A1 = 1.0
    B1 = [4.0 5.0; 6.0 7.0]
    A2 = [1.1 1.2; 1.3 1.4]
    B2 = 5.0
    dY_vals = spzeros(4,4)
    dY_vals[1,2] = 0.1; dY_vals[2,1] = 0.1; 
    dY_vals[1,3] = 0.2; dY_vals[3,1] = 0.2; 
    dY_vals[2,4] = 0.3; dY_vals[4,2] = 0.3; 
    dY_vals[3,4] = 0.4; dY_vals[4,3] = 0.4; 
    dX_vals = spzeros(4,4)
    dX_vals[1,2] = 0.01; dX_vals[2,1] = 0.01; 
    dX_vals[1,3] = 0.02; dX_vals[3,1] = 0.02; 
    dX_vals[2,4] = 0.03; dX_vals[4,2] = 0.03; 
    dX_vals[3,4] = 0.04; dX_vals[4,3] = 0.04; 
    dX1 = 0.1
    dY1 = dY_vals
    dX2 = dX_vals
    dY2 = 0.01
    ps_1 = [:A => A1, :B => B1, :dX => dX1, :dY => dY1]
    ps_2 = [:A => A2, :B => B2, :dX => dX2, :dY => dY2]

    # Creates simulation through two different separate simulations.
    oprob_1_1 = ODEProblem(lrs, u0, (0.0, 5.0), ps_1; jac = true, sparse = true)
    sol_1_1 = solve(oprob_1_1, Rosenbrock23(); saveat = 1.0, abstol = 1e-8, reltol = 1e-8)
    u0_1_2 = [:X => sol_1_1.u[end][1:2:end], :Y => sol_1_1.u[end][2:2:end]]
    oprob_1_2 = ODEProblem(lrs, u0_1_2, (0.0, 5.0), ps_2; jac = true, sparse = true)
    sol_1_2 = solve(oprob_1_2, Rosenbrock23(); saveat = 1.0, abstol = 1e-8, reltol = 1e-8)

    # Creates simulation through a single simulation with a callback
    oprob_2 = ODEProblem(lrs, u0, (0.0, 10.0), ps_1; jac = true, sparse = true)
    condition(u, t, integrator) = (t == 5.0)
    function affect!(integrator)
        lat_setp!(integrator, :A, lrs, A2)
        lat_setp!(integrator, :B, lrs, B2)
        integrator.ps[:dX] = dX2
        integrator.ps[:dY] = [dY2]
        rebuild_lat_internals!(integrator)
    end
    callback = DiscreteCallback(condition, affect!)
    sol_2 = solve(oprob_2, Rosenbrock23(); saveat = 1.0, tstops = [5.0], callback, abstol = 1e-8, reltol = 1e-8)

    # Check that trajectories are equivalent.
    @test [sol_1_1.u; sol_1_2.u] ≈ sol_2.u
end

# Currently not supported for jump stuff, check that corresponding functions yield errors.
let
    # Prepare `LatticeReactionSystem`.
    rs = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    tr = @transport_reaction D X1 
    grid = CartesianGrid((2,2))
    lrs = LatticeReactionSystem(rs, [tr], grid)

    # Create problems.
    u0 = [:X1 => 2, :X2 => [5 6; 7 8]]
    tspan = (0.0, 10.0)
    ps = [:k1 => 1.5, :k2 => [1.0 1.5; 2.0 3.5], :D => 0.1]
    dprob = DiscreteProblem(lrs, u0, tspan, ps)
    jprob = JumpProblem(lrs, dprob, NSM())
    jinit = init(jprob, SSAStepper())

    # Checks that rebuilding errors.
    @test_throws Exception rebuild_lat_internals!(dprob)
    @test_throws Exception rebuild_lat_internals!(jprob)
    @test_throws Exception rebuild_lat_internals!(jinit)
end
