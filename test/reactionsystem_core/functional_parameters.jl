
# Fetch packages.
using Catalyst, DataInterpolations, DynamicQuantities, JumpProcesses, OrdinaryDiffEqDefault, OrdinaryDiffEqTsit5, StochasticDiffEq, Test

# Sets the default `t` to use.
t = default_t()

### Basic Checks ###

# Checks with oscillating time function.
# Checks that interpolating into DSL yields identical results (both simulation and model).
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
    rs_pIn_dsl = @reaction_network rs_pIn begin
        ($pIn(t),d), 0 <--> X
        (k1*X,k2), Y1 <--> Y2
    end

    # Checks that the two model declarations are equivalent.
    @test Catalyst.isequivalent(rs_pIn, rs_pIn_dsl)

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
    sde_sol_mean_base = EnsembleAnalysis.timeseries_point_mean(solve(esde_prob_base, ImplicitEM(); trajectories = 1000), 1.0:10.0).u
    sde_sol_mean_pIn = EnsembleAnalysis.timeseries_point_mean(solve(esde_prob_pIn, ImplicitEM(); trajectories = 1000), 1.0:10.0).u
    @test sde_sol_mean_base ≈ sde_sol_mean_pIn rtol = 1e-1 atol = 1e-1

    # Checks Jump simulations (uses Tsit5 because model has time-dependent rates -> VariableRateJumps).
    ejmp_prob_base = EnsembleProblem(JumpProblem(rs_base, u0, (0.0, tend), ps_base))
    ejmp_prob_pIn = EnsembleProblem(JumpProblem(rs_pIn, u0, (0.0, tend), ps))
    jmp_sol_mean_base = EnsembleAnalysis.timeseries_point_mean(solve(ejmp_prob_base, Tsit5(); trajectories = 1000), 1.0:10.0).u
    jmp_sol_mean_pIn = EnsembleAnalysis.timeseries_point_mean(solve(ejmp_prob_pIn, Tsit5(); trajectories = 1000), 1.0:10.0).u
    @test jmp_sol_mean_base ≈ jmp_sol_mean_pIn rtol = 1e-1 atol = 1e-1
end

# Checks correctness for non-time functions.
# Uses an SIR model where we have modified the infection step to not scale linearly with I.
# Intended to add cases with 2d functions here as well, but this is currently not supported.
let
    # Declares the functional parameter.
    Is = collect(0.0:0.0001:200.0)
    spline1d = LinearInterpolation(100.0*Is ./ (100.0 .+ Is), Is)
    @parameters (i_rate::typeof(spline1d))(..)

    # Declares the models.
    sir = @reaction_network begin
        k1*100*I/(100 + I), S --> I
        k2, I --> R
    end
    sir_funcp = @reaction_network rs begin
        k1*$i_rate(I), S --> I
        k2, I --> R
    end

    # Simulates the models for the same conditions. Check that simulations are identical.
    u0 = [:S => 99.0, :I => 10.0, :R => 0.0]
    ps = [:k1 => 0.01, :k2 => 0.01]
    oprob = ODEProblem(sir, u0, 200.0, ps)
    oprob_funcp = ODEProblem(sir_funcp, u0, 200.0, [ps; i_rate => spline1d])
    sol = solve(oprob)
    sol_funcp = solve(oprob_funcp)
    @test sol.u ≈ sol_funcp.u atol = 1e-6 rtol = 1e-6
end

### Other tests ###

# Checks the combination of functional parameters and units.
# Unclear what the final expected behaviour should be here (read https://github.com/SciML/ModelingToolkit.jl/issues/3420).
# `get_unit` works when called on the functional parameters as a call (`get_unit(pIn(t))`),
# but not just on itself (`get_unit(pIn)`). This is currently how things work, and if MTK
# changes how things work we will at least detect it in this test.
let
    # Defines an input process (just some decaying function of t).
    ts = collect(0.0:0.01:1.0)
    spline = LinearInterpolation(1 ./ (1 .+ ts), ts)
    @parameters (pIn::typeof(spline))(..) [unit=u"1/(s*m^3)"]

    # Defines a `ReactionSystem` using the input parameter (birth/death process, birth split in two parameters).
    # Checks that the units of the reaction rates are correct.
    @independent_variables t [unit=u"s"]
    @species X(t) [unit=u"mol/m^3"]
    @parameters p_base [unit=u"mol"] d [unit=u"1/s"]
    rxs = [
        Reaction(p_base*pIn(t), [], [X])
        Reaction(d, [X], [])
    ]
    @named rs = ReactionSystem(rxs, t)
    @test issetequal([Catalyst.get_unit(rx.rate) for rx in reactions(rs)], [u"mol/(s*m^3)", u"1/s"])
end

# Tests that a functional parameter can be interpolated as the function only, as a function of a 
# symbolic variable, or interpolated with the functional parameter and argument separately.
let 
    # Prepares the functional parameter.
    ts = collect(0.0:0.1:1.0)
    spline = LinearInterpolation(log.(ts), ts)
    t_var = default_t()
    @parameters (pIn::typeof(spline))(..)
    pIn_func = pIn(t_var)

    # Creates models using different approaches, check that all yield the same model.
    rn1 = @reaction_network rn begin
        p, 0 --> X
        $pIn_func, 0 --> X
    end
    rn2 = @reaction_network rn begin
        p, 0 --> X
        $pIn(t), 0 --> X
    end
    rn3 = @reaction_network rn begin
        p, 0 --> X
        $pIn($t_var), 0 --> X
    end
    @test all(rn -> Catalyst.isequivalent(rn1, rn), (rn2, rn3))
end
