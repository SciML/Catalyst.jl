### Preparations ###

# Fetch packages.
using Catalyst, Graphs, JumpProcesses, OrdinaryDiffEq, SparseArrays, Test

### `get_lrs_vals` Tests ###

# Basic test. For simulations without effect of system, check that solution correspond to known
# initial condition throughout solution. 
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

    # Loops through a lattice cases and check that they are correct.
    for (u0,lrs) in zip([u0_1, u0_2, u0_3, u0_4, u0_5, u0_6], [lrs1, lrs2, lrs3, lrs4, lrs5, lrs6])
        # Simulates ODE version and checks `get_lrs_vals` on its solution.
        oprob = ODEProblem(lrs, u0, tspan, ps)
        osol = solve(oprob, Tsit5(), saveat = 0.5)
        @test get_lrs_vals(osol, :X1, lrs) == get_lrs_vals(osol, :X1, lrs; t = 0.0:0.5:1.0)
        @test all(all(val == Float64(u0[:X1]) for val in vals) for vals in get_lrs_vals(osol, :X1, lrs))
        @test get_lrs_vals(osol, :X2, lrs) == get_lrs_vals(osol, :X2, lrs; t = 0.0:0.5:1.0) == fill(u0[:X2], 3)

        # Simulates jump version and checks `get_lrs_vals` on its solution.
        dprob = DiscreteProblem(lrs, u0, tspan, ps)
        jprob = JumpProblem(lrs, dprob, NSM())
        jsol = solve(jprob, SSAStepper(), saveat = 0.5)
        @test get_lrs_vals(jsol, :X1, lrs) == get_lrs_vals(jsol, :X1, lrs; t = 0.0:0.5:1.0)
        @test all(all(val == Float64(u0[:X1]) for val in vals) for vals in get_lrs_vals(jsol, :X1, lrs))
        @test get_lrs_vals(jsol, :X2, lrs) == get_lrs_vals(jsol, :X2, lrs; t = 0.0:0.5:1.0) == fill(u0[:X2], 3)
    end
end

# Checks on simulations where the system changes in time.
# Checks that solution have correct initial condition and end point (steady state).
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
    vals = get_lrs_vals(sol, :X, lrs)
    @test vals[1] == [1.0, 3.0]
    @test vals[end] ≈ [2.0, 2.0]
    for i = 1:(length(saveat) - 1)
        @test vals[i][1] < vals[i + 1][1]
        @test vals[i][2] > vals[i + 1][2]
    end
end

# Checks interpolation when sampling at time point. Check that value at `t` is inbetween the 
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
    t5_vals = get_lrs_vals(sol, :X, lrs; t = [0.5])[1]
    @test sol.u[1][1] < t5_vals[1] < sol.u[2][1]
    @test sol.u[1][2] > t5_vals[2] > sol.u[2][2]
end

### Error Tests ###

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
    @test_throws Exception get_lrs_vals(sol, :X, lrs; t = [0.0])
    @test_throws Exception get_lrs_vals(sol, :X, lrs; t = [3.0])
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
    @test_throws Exception get_lrs_vals(sol, :X, lrs; t = [0.0])
    @test_throws Exception get_lrs_vals(sol, :X, lrs; t = [3.0])
end

# Checks that applying `get_lrs_vals` to a 3d masked lattice yields an error.
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
    @test_throws Exception get_lrs_vals(sol, :X, lrs)
end

### Other Tests ###

# Checks that `get_lrs_vals` works for all types of symbols.
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
    @test get_lrs_vals(sol, X, lrs) == get_lrs_vals(sol, rs.X, lrs) == get_lrs_vals(sol, :X, lrs)
    @test get_lrs_vals(sol, X, lrs; t = 0.0:0.5:1.0) == get_lrs_vals(sol, rs.X, lrs; t = 0.0:0.5:1.0) == get_lrs_vals(sol, :X, lrs; t = 0.0:0.5:1.0)
end