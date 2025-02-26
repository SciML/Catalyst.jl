
# Fetch packages.
using Catalyst, DataInterpolations, JumpProcesses, OrdinaryDiffEq, StochasticDiffEq 

# Sets the default `t` to use.
t = default_t()

### Basic Checks ###

# Checks with oscillating time function.
# Checks that interpolating into DSL yields identical result (both simulation and model).
# Checks that ODE simulations are correct.
let
    # Defines an input process (modified sinus wave).
    tend = 5.0
    ts = collect(0.0:0.0001:tend)
    A = 2.0; f = 1.6; ϕ = 0.5;
    spline = LinearInterpolation(A .* (sin.(2π * f .* ts .- ϕ) .+ 1) /2, ts)
    @parameters (pIn::typeof(spline))(..)

    # Defines a `ReactionSystem` without the input parameter.
    @species X(t) Y1(t) Y2(t)
    @parameters p d k1 k2
    rxs_base = [
        Reaction(A * (sin(2π * f * t .- ϕ) + 1) /2, [], [X])
        Reaction(d, [X], [])
        Reaction(k1*X, [Y1], [Y2])
        Reaction(k2, [Y2], [Y1])
    ]
    @named rs_base = ReactionSystem(rxs_base, t)
    rs_base = complete(rs_base)

    # Defines a `ReactionSystem` with the input parameter (Programmatically).
    rxs_pIn = [
        Reaction(pIn(t), [], [X])
        Reaction(d, [X], [])
        Reaction(k1*X, [Y1], [Y2])
        Reaction(k2, [Y2], [Y1])
    ]
    @named rs_pIn = ReactionSystem(rxs_pIn, t)
    rs_pIn = complete(rs_pIn)

    # Defines a `ReactionSystem` with the input parameter (DSL).
    input = pIn(t)
    rs_pIn_dsl = @reaction_network rs_pIn begin
        ($input,d), 0 <--> X
        (k1*X,k2), Y1 <--> Y2
    end

    # Checks that the two model declarations are identical.
    @test isequal(rs_pIn, rs_pIn_dsl)

    # Makes ODE simulation using different approaches.
    u0 = [X => 1.0, Y1 => 2.0, Y2 => 3.0]
    ps_base = [d => 0.5, k1 => 0.1, k2 => 0.4]
    ps = [ps_base; pIn => spline]
    oprob_base = ODEProblem(rs_base, u0, tend, ps_base)
    oprob_prog = ODEProblem(rs_pIn, u0, tend, ps)
    oprob_dsl = ODEProblem(rs_pIn_dsl, u0, tend, ps)
    sol_base = solve(oprob_base; saveat = 0.1, abstol = 1e-8, reltol = 1e-8)
    sol_prog = solve(oprob_prog; saveat = 0.1, abstol = 1e-8, reltol = 1e-8)
    sol_dsl = solve(oprob_dsl; saveat = 0.1, abstol = 1e-8, reltol = 1e-8)

    # Checks that simulations are identical.
    @test sol_base[X] ≈ sol_prog[X] atol = 1e-8 rtol = 1e-8
    @test sol_base[X] ≈ sol_dsl[X] atol = 1e-8 rtol = 1e-8
end

# Checks a simple model for different problem types.
# Currently works for ODE and SDE, not for Jump.
let
    # Defines an input process (a decaying function).
    tend = 10.0
    ts = collect(0.0:0.001:tend)
    spline = LinearInterpolation((10 .+ ts) ./ (2 .+ ts), ts)
    @parameters (pIn::typeof(spline))(..)

    # Defines the `ReactionSystem` with/without the temporal parameter.
    @species X(t)
    @parameters d
    @named rs_base = ReactionSystem([
        Reaction((10 + t) / (2 + t), [], [X]),
        Reaction(d, [X], [])
    ])
    @named rs_pIn = ReactionSystem([
        Reaction(pIn(t), [], [X]),
        Reaction(d, [X], [])
    ])
    rs_base = complete(rs_base)
    rs_pIn = complete(rs_pIn)

    # Sets simulation conditions.
    u0 = [X => 50]
    ps_base = [d => 0.1]
    ps = [ps_base; pIn => spline]

    # Checks ODE simulations.
    ode_sol_base = solve(ODEProblem(rs_base, u0, tend, ps_base); reltol = 1e-8, abstol = 1e-8, saveat = 0.1)
    ode_sol_pIn = solve(ODEProblem(rs_pIn, u0, tend, ps); reltol = 1e-8, abstol = 1e-8, saveat = 0.1)
    @test ode_sol_base ≈ ode_sol_pIn

    # Checks SDE simulations.
    esde_prob_base = EnsembleProblem(SDEProblem(rs_base, u0, tend, ps_base))
    esde_prob_pIn = EnsembleProblem(SDEProblem(rs_pIn, u0, tend, ps))
    sde_sol_mean_base = timeseries_point_mean(solve(esde_prob_base, ImplicitEM(); trajectories = 1000), 1.0:10.0).u
    sde_sol_mean_pIn = timeseries_point_mean(solve(esde_prob_pIn, ImplicitEM(); trajectories = 1000), 1.0:10.0).u
    @test sde_sol_mean_base ≈ sde_sol_mean_pIn rtol = 1e-1 atol = 1e-1

    # Checks Jump simulations (cannot currently be created).
    # ejmp_prob_base = EnsembleProblem(JumpProblem(JumpInputs(rs_base, u0, tend, ps_base)))
    # ejmp_prob_pIn = EnsembleProblem(JumpProblem(JumpInputs(rs_pIn, u0, tend, ps)))
    # jmp_sol_mean_base = timeseries_point_mean(solve(ejmp_prob_base; trajectories = 1000), 1.0:10.0).u
    # jmp_sol_mean_pIn = timeseries_point_mean(solve(ejmp_prob_pIn; trajectories = 1000), 1.0:10.0).u
    @test_broken jmp_sol_mean_base ≈ jmp_sol_meanpIn rtol = 1e-1 atol = 1e-1
end