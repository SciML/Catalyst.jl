### Preparations ###

# Fetch packages.
using Catalyst, Graphs, JumpProcesses, OrdinaryDiffEqVerner, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, SparseArrays, Test

# Fetch test networks.
include("../spatial_test_networks.jl")

### Problem & Integrator Interfacing Function Tests ###

# Checks `spat_getu` for ODE and Jump problem and integrators.
# Checks `spat_setu!` for ODE and Jump problem and integrators.
# Checks for all types of spaces.
# Checks for symbol and symbolic variables input.
let
    # Declares various types of spaces and corresponding initial values of `X`.
    dspace_cartesian = CartesianGrid((2,2,2))
    dspace_masked = [true true; false true]
    dspace_graph = cycle_graph(5)
    val0_cartesian = fill(1.0, 2, 2, 2)
    val0_masked = sparse([1.0 2.0; 0.0 3.0])
    val0_graph = [1.0, 2.0, 3.0, 4.0, 5.0]

    # Unpacks the `X`  and `Y` symbolic variable (so that indexing using it can be tested).
    @unpack X, Y = brusselator_system

    # Loops through all alternative spaces and values. Checks that `spat_getu` works in all cases.
    for (space, val0) in zip([dspace_cartesian, dspace_masked, dspace_graph], [val0_cartesian, val0_masked, val0_graph])
        # Prepares various problems and integrators. Uses `deepcopy` to ensure there is no cross-talk
        # between the different u vectors as they get updated.
        dsrs = DiscreteSpaceReactionSystem(brusselator_system, brusselator_srs_1, space)
        u0 = [:X => val0, :Y => 0.5]
        ps = [:A => 1.0, :B => 2.0, :dX => 0.1]
        oprob = ODEProblem(dsrs, deepcopy(u0), (0.0, 1.0), ps)
        jprob = JumpProblem(dsrs, deepcopy(u0), (0.0, 1.0), ps)
        oint = init(deepcopy(oprob), Tsit5())
        jint = init(deepcopy(jprob), SSAStepper())
        
        # Check that `spat_getu` retrieves the correct values.
        @test spat_getu(oprob, :X, dsrs) == spat_getu(oprob, X, dsrs) == spat_getu(oprob, brusselator_system.X, dsrs) == val0
        @test spat_getu(oint, :X, dsrs) == spat_getu(oint, X, dsrs) == spat_getu(oint, brusselator_system.X, dsrs) == val0
        @test spat_getu(jprob, :X, dsrs) == spat_getu(jprob, X, dsrs) == spat_getu(jprob, brusselator_system.X, dsrs) == val0
        @test spat_getu(jint, :X, dsrs) == spat_getu(jint, X, dsrs) == spat_getu(jint, brusselator_system.X, dsrs) == val0
        
        # Updates Y and checks its content.
        spat_setu!(oprob, :Y, dsrs, val0)
        @test spat_getu(oprob, :Y, dsrs) == spat_getu(oprob, Y, dsrs) == spat_getu(oprob, brusselator_system.Y, dsrs) == val0
        spat_setu!(oint, :Y, dsrs, val0)
        @test spat_getu(oint, :Y, dsrs) == spat_getu(oint, Y, dsrs) == spat_getu(oint, brusselator_system.Y, dsrs) == val0
        spat_setu!(jprob, :Y, dsrs, val0)
        @test spat_getu(jprob, :Y, dsrs) == spat_getu(jprob, Y, dsrs) == spat_getu(jprob, brusselator_system.Y, dsrs) == val0
        spat_setu!(jint, :Y, dsrs, val0)
        @test spat_getu(jint, :Y, dsrs) == spat_getu(jint, Y, dsrs) == spat_getu(jint, brusselator_system.Y, dsrs) == val0

        # Tries where we change a spatially non-uniform variable to spatially uniform.
        spat_setu!(oprob, X, dsrs, 0.0)
        @test all(isequal(0.0), spat_getu(oprob, X, dsrs))
        spat_setu!(oint, X, dsrs, 0.0)
        @test all(isequal(0.0), spat_getu(oint, X, dsrs))
        spat_setu!(jprob, X, dsrs, 0.0)
        @test all(isequal(0.0), spat_getu(jprob, X, dsrs))
        spat_setu!(jint, X, dsrs, 0.0)
        @test all(isequal(0.0), spat_getu(jint, X, dsrs))
    end
end

# Checks `spat_getp` for ODEproblem and integrators.
# Checks `spat_setp!` for ODE problem and integrators.
# Checks for all types of spaces.
# Checks for symbol and symbolic variables input.
let
    # Declares various types of spaces and corresponding initial values of `A`.
    dspace_cartesian = CartesianGrid((2,2,2))
    dspace_masked = [true true; false true]
    dspace_graph = cycle_graph(5)
    val0_cartesian = fill(1.0, 2, 2, 2)
    val0_masked = sparse([1.0 2.0; 0.0 3.0])
    val0_graph = [1.0, 2.0, 3.0, 4.0, 5.0]

    # Unpacks the `A` and `B` symbolic variable (so that indexing using it can be tested).
    @unpack A, B = brusselator_system

    # Loops through all alternative spaces and values. Checks that `spat_getp` works in all cases.
    for (space, val0) in zip([dspace_cartesian, dspace_masked, dspace_graph], [val0_cartesian, val0_masked, val0_graph])
        # Prepares various problems and integrators. Uses `deepcopy` to ensure there is no cross-talk
        # between the different p vectors as they get updated.
        dsrs = DiscreteSpaceReactionSystem(brusselator_system, brusselator_srs_1, space)
        u0 = [:X => 1.0, :Y => 0.5]
        ps = [:A => val0, :B => 2.0, :dX => 0.1]
        oprob = ODEProblem(dsrs, u0, (0.0, 1.0), deepcopy(ps))
        oint = init(deepcopy(oprob), Tsit5())
        
        # Check that `spat_getp` retrieves the correct values.
        @test spat_getp(oprob, :A, dsrs) == spat_getp(oprob, A, dsrs) == spat_getp(oprob, brusselator_system.A, dsrs) == val0
        @test spat_getp(oint, :A, dsrs) == spat_getp(oint, A, dsrs) == spat_getp(oint, brusselator_system.A, dsrs) == val0
        
        # Updates Y and checks its content.
        spat_setp!(oprob, :B, dsrs, val0)
        @test spat_getp(oprob, :B, dsrs) == spat_getp(oprob, B, dsrs) == spat_getp(oprob, brusselator_system.B, dsrs) == val0
        spat_setp!(oint, :B, dsrs, val0)
        @test spat_getp(oint, :B, dsrs) == spat_getp(oint, B, dsrs) == spat_getp(oint, brusselator_system.B, dsrs) == val0

        # Tries where we change a spatially non-uniform variable to spatially uniform.
        spat_setp!(oprob, A, dsrs, 0.0)
        @test all(isequal(0.0), spat_getp(oprob, A, dsrs))
        spat_setp!(oint, A, dsrs, 0.0)
        @test all(isequal(0.0), spat_getp(oint, A, dsrs))
    end
end

# Checks that `spat_getp` and `spat_setp!` generates errors when applied to `JumpProblem`s and their integrators.
let
    dsrs = DiscreteSpaceReactionSystem(brusselator_system, brusselator_srs_1, very_small_1d_cartesian_grid)
    u0 = [:X => 1, :Y => 0]
    ps = [:A => 3.0, :B => 2.0, :dX => 0.1]
    jprob = JumpProblem(dsrs, u0, (0.0, 1.0), ps)
    jint = init(jprob, SSAStepper())

    @test_throws Exception spat_getp(jprob, :A, dsrs)
    @test_throws Exception spat_getp(jprob, :A, dsrs)
    @test_throws Exception spat_setp!(jint, :A, dsrs, 0.0)
    @test_throws Exception spat_setp!(jint, :A, dsrs, 0.0)
end

# Checks that `spat_getp` and `spat_setp!` generates errors when applied to edge parameters.
let
    dsrs = DiscreteSpaceReactionSystem(brusselator_system, brusselator_srs_1, very_small_1d_cartesian_grid)
    u0 = [:X => 1.0, :Y => 0.0]
    ps = [:A => 3.0, :B => 2.0, :dX => 0.1]
    oprob = ODEProblem(dsrs, u0, (0.0, 1.0), ps)
    oint = init(deepcopy(oprob), Tsit5())

    @test_throws ArgumentError spat_getp(oprob, :dX, dsrs)
    @test_throws ArgumentError spat_setp!(oprob, :dX, dsrs, 0.0)
end

### Simulation `spat_getu` Tests ###

# Basic test. For simulations without change in system, check that the solution corresponds to known
# initial condition throughout the solution. 
# Checks using both `t` sampling` and normal time step sampling.
# Checks for both ODE and jump simulations.
# Checks for all discrete space types.
let 
    # Prepare `DiscreteSpaceReactionSystem`s.
    rs = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    tr = @transport_reaction D X1
    dsrs1 = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,)))
    dsrs2 = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,3)))
    dsrs3 = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,3,2)))
    dsrs4 = DiscreteSpaceReactionSystem(rs, [tr], [true, true, false, true])
    dsrs5 = DiscreteSpaceReactionSystem(rs, [tr], [true false; true true])
    dsrs6 = DiscreteSpaceReactionSystem(rs, [tr], cycle_graph(4))

    # Create problem inputs.
    u0_1 = Dict([:X1 => 0, :X2 => [1, 2]])
    u0_2 = Dict([:X1 => 0, :X2 => [1 2 3; 4 5 6]])
    u0_3 = Dict([:X1 => 0, :X2 => fill(1, 2, 3, 2)])
    u0_4 = Dict([:X1 => 0, :X2 => sparse([1, 2, 0, 3])])
    u0_5 = Dict([:X1 => 0, :X2 => sparse([1 0; 2 3])])
    u0_6 = Dict([:X1 => 0, :X2 => [1, 2, 3, 4]])
    tspan = (0.0, 1.0)
    ps = [:k1 => 0.0, :k2 => 0.0, :D => 0.0]

    # Loops through all discrete space cases and check that they are correct.
    for (u0,dsrs) in zip([u0_1, u0_2, u0_3, u0_4, u0_5, u0_6], [dsrs1, dsrs2, dsrs3, dsrs4, dsrs5, dsrs6])
        # Simulates ODE version and checks `spat_getu` on its solution.
        oprob = ODEProblem(dsrs, u0, tspan, ps)
        osol = solve(oprob, Tsit5(), saveat = 0.5)
        @test spat_getu(osol, :X1, dsrs) == spat_getu(osol, :X1, dsrs; t = 0.0:0.5:1.0)
        @test all(all(val == Float64(u0[:X1]) for val in vals) for vals in spat_getu(osol, :X1, dsrs))
        @test spat_getu(osol, :X2, dsrs) == spat_getu(osol, :X2, dsrs; t = 0.0:0.5:1.0) == fill(u0[:X2], 3)

        # Simulates jump version and checks `spat_getu` on its solution.
        jprob = JumpProblem(dsrs, u0, tspan, ps)
        jsol = solve(jprob, SSAStepper(), saveat = 0.5)
        @test spat_getu(jsol, :X1, dsrs) == spat_getu(jsol, :X1, dsrs; t = 0.0:0.5:1.0)
        @test all(all(val == Float64(u0[:X1]) for val in vals) for vals in spat_getu(jsol, :X1, dsrs))
        @test spat_getu(jsol, :X2, dsrs) == spat_getu(jsol, :X2, dsrs; t = 0.0:0.5:1.0) == fill(u0[:X2], 3)
    end
end

# Checks on simulations where the system changes in time.
# Checks that a solution has correct initial condition and end point (steady state).
# Checks that solution is monotonously increasing/decreasing (it should be for this problem).
let
    # Prepare `DiscreteSpaceReactionSystem`s.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    dsrs = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,)))

    # Prepares a corresponding ODEProblem.
    u0 = [:X => [1.0, 3.0]]
    tspan = (0.0, 50.0)
    ps = [:p => 2.0, :d => 1.0, :D => 0.01]
    oprob = ODEProblem(dsrs, u0, tspan, ps)

    # Simulates the ODE. Checks that the start/end points are correct.
    # Check that the first vertex is monotonously increasing in values, and that the second one is 
    # monotonously decreasing. The non evenly spaced `saveat` is so that non-monotonicity is
    # not produced due to numeric errors.
    saveat = [0.0, 1.0, 5.0, 10.0, 50.0]
    sol = solve(oprob, Vern7(); abstol = 1e-8, reltol = 1e-8)
    vals = spat_getu(sol, :X, dsrs)
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
    # Prepare `DiscreteSpaceReactionSystem`s.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    dsrs = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,)))

    # Solved a corresponding ODEProblem.
    u0 = [:X => [1.0, 3.0]]
    tspan = (0.0, 1.0)
    ps = [:p => 2.0, :d => 1.0, :D => 0.0]
    oprob = ODEProblem(dsrs, u0, tspan, ps)

    # Solves and check the interpolation of t.
    sol = solve(oprob, Tsit5(); saveat = 1.0)
    t5_vals = spat_getu(sol, :X, dsrs; t = [0.5])[1]
    @test sol.u[1][1] < t5_vals[1] < sol.u[2][1]
    @test sol.u[1][2] > t5_vals[2] > sol.u[2][2]
end

# Checks that attempting to sample `t` outside tspan range yields an error.
let
    # Prepare `DiscreteSpaceReactionSystem`s.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    dsrs = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,)))

    # Solved a corresponding ODEProblem.
    u0 = [:X => 1.0]
    tspan = (1.0, 2.0)
    ps = [:p => 2.0, :d => 1.0, :D => 1.0]
    oprob = ODEProblem(dsrs, u0, tspan, ps)

    # Solves and check the interpolation of t.
    sol = solve(oprob, Tsit5(); saveat = 1.0)
    @test_throws Exception spat_getu(sol, :X, dsrs; t = [0.0])
    @test_throws Exception spat_getu(sol, :X, dsrs; t = [3.0])
end

# Checks that attempting to sample `t` outside tspan range yields an error.
let
    # Prepare `DiscreteSpaceReactionSystem`s.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    dsrs = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,)))

    # Solved a corresponding ODEProblem.
    u0 = [:X => 1.0]
    tspan = (1.0, 2.0)
    ps = [:p => 2.0, :d => 1.0, :D => 1.0]
    oprob = ODEProblem(dsrs, u0, tspan, ps)

    # Solves and check the interpolation of t.
    sol = solve(oprob, Tsit5(); saveat = 1.0)
    @test_throws Exception spat_getu(sol, :X, dsrs; t = [0.0])
    @test_throws Exception spat_getu(sol, :X, dsrs; t = [3.0])
end

# Checks that applying `spat_getu` to a 3d masked space yields an error.
let
    # Prepare `DiscreteSpaceReactionSystem`s.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    tr = @transport_reaction D X
    dsrs = DiscreteSpaceReactionSystem(rs, [tr], rand([false, true], 2, 3, 4))

    # Solved a corresponding ODEProblem.
    u0 = [:X => 1.0]
    tspan = (1.0, 2.0)
    ps = [:p => 2.0, :d => 1.0, :D => 1.0]
    oprob = ODEProblem(dsrs, u0, tspan, ps)

    # Solves and check the interpolation of t.
    sol = solve(oprob, Tsit5(); saveat = 1.0)
    @test_throws Exception spat_getu(sol, :X, dsrs)
end

# Checks that `spat_getu` works for all types of symbols.
let
    t = default_t()
    @species X(t)
    @parameters d
    @named rs = ReactionSystem([Reaction(d, [X], [])], t)
    rs = complete(rs)
    tr = @transport_reaction D X
    dsrs = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid(2,))

    # Solved a corresponding ODEProblem.
    u0 = [:X => 1.0]
    tspan = (0.0, 1.0)
    ps = [:d => 1.0, :D => 0.1]
    oprob = ODEProblem(dsrs, u0, tspan, ps)

    # Solves and check the interpolation of t.
    sol = solve(oprob, Tsit5(); saveat = 1.0)
    @test spat_getu(sol, X, dsrs) == spat_getu(sol, rs.X, dsrs) == spat_getu(sol, :X, dsrs)
    @test spat_getu(sol, X, dsrs; t = 0.0:0.5:1.0) == spat_getu(sol, rs.X, dsrs; t = 0.0:0.5:1.0) == spat_getu(sol, :X, dsrs; t = 0.0:0.5:1.0)
end

### ODEProblem & Integrator Rebuilding ###

# Checks that the `rebuild_spat_internals!` function is correctly applied to an ODEProblem.
let
    # Creates a Brusselator `DiscreteSpaceReactionSystem`.
    dsrs = DiscreteSpaceReactionSystem(brusselator_system, brusselator_srs_2, very_small_2d_cartesian_grid)

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
        oprob_1 = ODEProblem(dsrs, u0, (0.0, 10.0), ps; jac, sparse)

        # Creates an alternative version of the ODEProblem.
        dX_vals = spzeros(4,4)
        dX_vals[1,2] = 0.01; dX_vals[2,1] = 0.01; 
        dX_vals[1,3] = 0.02; dX_vals[3,1] = 0.02; 
        dX_vals[2,4] = 0.03; dX_vals[4,2] = 0.03; 
        dX_vals[3,4] = 0.04; dX_vals[4,3] = 0.04; 
        ps = [:A => [1.1 1.2; 1.3 1.4], :B => 5.0, :dX => dX_vals, :dY => 0.01]
        oprob_2 = ODEProblem(dsrs, u0, (0.0, 10.0), ps; jac, sparse)

        # Modifies the initial ODEProblem to be identical to the new one.
        spat_setp!(oprob_1, :A, dsrs, [1.1 1.2; 1.3 1.4])
        spat_setp!(oprob_1, :B, dsrs, 5.0)
        oprob_1.ps[:dX] = dX_vals
        oprob_1.ps[:dY] = [0.01]
        rebuild_spat_internals!(oprob_1)

        # Checks that simulations of the two `ODEProblem`s are identical.
        @test solve(oprob_1, Rodas5P()) ≈ solve(oprob_2, Rodas5P())
    end
end

# Checks that the `rebuild_spat_internals!` function is correctly applied to an integrator.
# Does through by applying it within a callback, and compare to simulations without callback.
# To keep test faster, only check for `jac = sparse = true` only.
let
    # Prepares problem inputs.
    dsrs = DiscreteSpaceReactionSystem(brusselator_system, brusselator_srs_2, very_small_2d_cartesian_grid)
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
    oprob_1_1 = ODEProblem(dsrs, u0, (0.0, 5.0), ps_1; jac = true, sparse = true)
    sol_1_1 = solve(oprob_1_1, Rosenbrock23(); saveat = 1.0, abstol = 1e-8, reltol = 1e-8)
    u0_1_2 = [:X => sol_1_1.u[end][1:2:end], :Y => sol_1_1.u[end][2:2:end]]
    oprob_1_2 = ODEProblem(dsrs, u0_1_2, (0.0, 5.0), ps_2; jac = true, sparse = true)
    sol_1_2 = solve(oprob_1_2, Rosenbrock23(); saveat = 1.0, abstol = 1e-8, reltol = 1e-8)

    # Creates simulation through a single simulation with a callback
    oprob_2 = ODEProblem(dsrs, u0, (0.0, 10.0), ps_1; jac = true, sparse = true)
    condition(u, t, integrator) = (t == 5.0)
    function affect!(integrator)
        spat_setp!(integrator, :A, dsrs, A2)
        spat_setp!(integrator, :B, dsrs, B2)
        integrator.ps[:dX] = dX2
        integrator.ps[:dY] = [dY2]
        rebuild_spat_internals!(integrator)
    end
    callback = DiscreteCallback(condition, affect!)
    sol_2 = solve(oprob_2, Rosenbrock23(); saveat = 1.0, tstops = [5.0], callback, abstol = 1e-8, reltol = 1e-8)

    # Check that trajectories are equivalent.
    @test [sol_1_1.u; sol_1_2.u] ≈ sol_2.u
end

# Currently not supported for jump stuff, check that corresponding functions yield errors.
let
    # Prepare `DiscreteSpaceReactionSystem`.
    rs = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    tr = @transport_reaction D X1 
    grid = CartesianGrid((2,2))
    dsrs = DiscreteSpaceReactionSystem(rs, [tr], grid)

    # Create problems.
    u0 = [:X1 => 2, :X2 => [5 6; 7 8]]
    tspan = (0.0, 10.0)
    ps = [:k1 => 1.5, :k2 => [1.0 1.5; 2.0 3.5], :D => 0.1]
    jprob = JumpProblem(dsrs, u0, tspan, ps)
    jinit = init(jprob, SSAStepper())

    # Checks that rebuilding errors.
    @test_throws Exception rebuild_spat_internals!(dprob)
    @test_throws Exception rebuild_spat_internals!(jprob)
    @test_throws Exception rebuild_spat_internals!(jinit)
end
