### Prepares Tests ###

# Fetch packages.
using Catalyst, JumpProcesses, OrdinaryDiffEqTsit5, StochasticDiffEq, Statistics, Test
using Symbolics: unwrap

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# Sets the default `t` to use.
t = default_t()

# Fetch test functions.
include("../test_functions.jl")


### Base Tests ###

# Checks that systems with symbolic stoichiometries, created using different approaches, are identical.
let
    @parameters p k d::Float64 n1::Int64 n2 n3
    @species X(t) Y(t)
    rxs1 = [
        Reaction(p, nothing, [X], nothing, [n1])
        Reaction(k, [X], [Y], [n2], [n3])
        Reaction(d, [Y], nothing)
    ]
    rs1 = complete(ReactionSystem(rxs1, t; name = :rs))

    rxs2 = [
        @reaction p, 0 --> $n1*X
        @reaction k, n2*X --> n3*Y
        @reaction $d, Y --> 0
    ]
    rs2 = complete(ReactionSystem(rxs2, t; name = :rs))

    rs3 = @reaction_network rs begin
        @parameters d::Float64 n1::Int64
        p, 0 --> n1*X
        k, n2*X --> n3*Y
        d, Y --> 0
    end

    @test Catalyst.isequivalent(rs1, rs2)
    @test Catalyst.isequivalent(rs2, rs3)
    @test issetequal(unknowns(rs1), [X, Y])
    @test issetequal(parameters(rs1), [p, k, d, n1, n2, n3])
    @test SymbolicUtils.symtype(d) == Float64
    @test SymbolicUtils.symtype(n1) == Int64
end

# Declares a network, parameter values, and initial conditions, to be used for the next couple of tests.
begin
    # Prepares parameters and species and their values.
    @parameters k α::Int64
    @species A(t), B(t), C(t), D(t)
    u0_1 = (A => 3.0, B => 2.0, C => 3.0, D => 5.0)
    ps_1 = (k => 5.0, α => 2)
    u0_2 = [u[2] for u in u0_1]
    ps_2 = Tuple(p[2] for p in ps_1)
    τ = 1.5

    # Creates `ReactionSystem` model.
    g = k + α * C
    rs = @reaction_network rs begin
        @parameters k α::Int64
        t*k, 2*(α^2)*A --> $g*B
        1.0, α*A + 2*B --> k*C + α*D
    end
end

# Compares the Catalyst-generated ODE function to a manually computed ODE function.
let
    # With combinatoric ratelaws.
    function oderhs(u, p, t)
        k,α = p
        A,B,C,D = u
        n = 2 * α^2
        rl = t * k / factorial(n) * A^n
        rl2 = A^α * B^2 / (2 * factorial(α))

        du = zeros(4)
        du[1] = -n * rl - α * rl2
        du[2] = (k + α * C) * rl - 2 * rl2
        du[3] = k * rl2
        du[4] = α * rl2
        return du
    end
    @test f_eval(rs, u0_1, ps_1, τ) ≈ oderhs(u0_2, ps_2, τ)

    # Without combinatoric ratelaws.
    function oderhs_no_crl(u, p, t)
        k,α = p
        A,B,C,D = u
        n = 2 * α^2
        rl = t * k * A^n
        rl2 = A^α * B^2

        du = zeros(4)
        du[1] = -n * rl - α * rl2
        du[2] = (k + α * C) * rl - 2 * rl2
        du[3] = k * rl2
        du[4] = α * rl2
        return du
    end
    @test f_eval(rs, u0_1, ps_1, τ; combinatoric_ratelaws = false) ≈ oderhs_no_crl(u0_2, ps_2, τ)
end

# Compares the Catalyst-generated SDE noise function to a manually computed SDE noise function.
let
    # With combinatoric ratelaws.
    function sdenoise(u, p, t)
        k,α = p
        A,B,C,D = u
        n = 2 * α^2
        rl = sqrt(t * k / factorial(n) * A^n)
        rl2 = sqrt(A^α * B^2 / (2 * factorial(α)))

        du = zeros(4,2)
        du = [-n*rl (-α*rl2);
             (k + α * C)*rl (-2*rl2);
             0.0 k*rl2;
             0.0 α*rl2]
        return du
    end
    @test g_eval(rs, u0_1, ps_1, τ) ≈ sdenoise(u0_2, ps_2, τ)

    # Without combinatoric ratelaws.
    function sdenoise_no_crl(u, p, t)
        k,α = p
        A,B,C,D = u
        n = 2 * α^2
        rl = sqrt(t * k * A^n)
        rl2 = sqrt(A^α * B^2)

        du = zeros(4,2)
        du = [-n*rl (-α*rl2);
             (k + α * C)*rl (-2*rl2);
             0.0 k*rl2;
             0.0 α*rl2]
        return du
    end
    @test g_eval(rs, u0_1, ps_1, τ; combinatoric_ratelaws = false) ≈ sdenoise_no_crl(u0_2, ps_2, τ)
end

# Compares the Catalyst-generated jump process jumps to manually computed jumps.
let
    # Manually declares the systems two jumps.
    function r1(u, p, t)
        k, α = p
        A = u[1]
        t * k * binomial(A, 2 * α^2)
    end
    function affect1!(integrator)
        k, α = integrator.p
        C = integrator.u[3]
        integrator.u[1] -= 2 * α^2
        integrator.u[2] += (k + α * C)
        nothing
    end
    function r2(u, p, t)
        k, α = p
        A,B = u[1:2]
        binomial(Int64(A), Int64(α)) * B * (B - 1) / 2
    end
    function affect2!(integrator)
        k, α = integrator.p
        integrator.u[1] -= α
        integrator.u[2] -= 2
        integrator.u[3] += k
        integrator.u[4] += α
        nothing
    end
    jumps = [VariableRateJump(r1, affect1!), VariableRateJump(r2, affect2!)]

    # Checks that the Catalyst-generated functions are equal to the manually declared ones.
    @test_broken for i in 1:2 # @Sam: something internal jump-related is breaking, can you have a look so that it is fixed right?
        catalyst_jsys = make_sck_jump(rs)
        unknownoid = Dict(unknown => i for (i, unknown) in enumerate(unknowns(catalyst_jsys)))
        catalyst_vrj = ModelingToolkitBase.assemble_vrj(catalyst_jsys, ModelingToolkitBase.jumps(catalyst_jsys)[i], unknownoid)
        @test isapprox(catalyst_vrj.rate(u0_2, ps_2, τ), jumps[i].rate(u0_2, ps_2, τ))

        fake_integrator1 = (u = copy(u0_2), p = ps_2, t = τ)
        fake_integrator2 = deepcopy(fake_integrator1)
        catalyst_vrj.affect!(fake_integrator1)
        jumps[i].affect!(fake_integrator2)
        @test fake_integrator1 == fake_integrator2
    end
end

# Tests symbolic stoichiometries in simulations.
# Tests for decimal numbered symbolic stoichiometries.
let
    # Declares models. The references models have the `n` parameters so they can use the same
    # parameter vectors as the non-reference ones.
    rs_int = @reaction_network begin
        @parameters n::Int64
        (k1, k2), n*X1 <--> X2
    end
    rs_dec = @reaction_network begin
        @parameters n::Float64
        (k1, k2), n*X1 <--> X2
    end
    rs_ref_int = @reaction_network begin
        @parameters n::Int64
        (k1, k2), 3*X1 <--> X2
    end
    rs_ref_dec = @reaction_network begin
        @parameters n::Float64
        (k1, k2), 2.5*X1 <--> X2
    end

    # Set simulation settings. Initial conditions are design to start, more or less, at
    # steady state concentrations.
    # Values are selected so that stochastic tests should always pass within the bounds (independent
    # of seed).
    u0_int = [:X1 => 150, :X2 => 600]
    u0_dec = [:X1 => 100, :X2 => 600]
    tspan_det = (0.0, 1.0)
    tspan_stoch = (0.0, 10000.0)
    ps_int = (:k1 => 0.00001, :k2 => 0.01, :n => 3)
    ps_dec = (:k1 => 0.00001, :k2 => 0.01, :n => 2.5)

    # Test ODE simulations with integer coefficients.
    oprob_int = ODEProblem(rs_int, u0_int, tspan_det, ps_int)
    oprob_int_ref = ODEProblem(rs_ref_int, u0_int, tspan_det, ps_int)
    @test solve(oprob_int, Tsit5())[:X1][end] ≈ solve(oprob_int_ref, Tsit5())[:X1][end]

    # Test ODE simulations with decimal coefficients.
    oprob_dec = ODEProblem(rs_dec, u0_dec, tspan_det, ps_dec; combinatoric_ratelaws = false)
    oprob_dec_ref = ODEProblem(rs_ref_dec, u0_dec, tspan_det, ps_dec; combinatoric_ratelaws = false)
    @test solve(oprob_dec, Tsit5())[:X1][end] ≈ solve(oprob_dec_ref, Tsit5())[:X1][end]

    # Test SDE simulations with integer coefficients.
    sprob_int = SDEProblem(rs_int, u0_int, tspan_stoch, ps_int)
    sprob_int_ref = SDEProblem(rs_ref_int, u0_dec, tspan_stoch, ps_int)
    ssol_int = solve(sprob_int, ImplicitEM(); seed)
    ssol_int_ref = solve(sprob_int_ref, ImplicitEM(); seed)
    @test mean(ssol_int[:X1]) ≈ mean(ssol_int_ref[:X1]) atol = 2*1e0

    # Test SDE simulations with decimal coefficients.
    sprob_dec = SDEProblem(rs_dec, u0_dec, tspan_stoch, ps_dec; combinatoric_ratelaws = false)
    sprob_dec_ref = SDEProblem(rs_ref_dec, u0_dec, tspan_stoch, ps_dec; combinatoric_ratelaws = false)
    ssol_dec = solve(sprob_dec, ImplicitEM(); seed)
    ssol_dec_ref = solve(sprob_dec_ref, ImplicitEM(); seed)
    @test mean(ssol_dec[:X1]) ≈ mean(ssol_dec_ref[:X1]) atol = 2*1e0

    # Test Jump simulations with integer coefficients.
    jprob_int = JumpProblem(rs_int, u0_int, tspan_stoch, ps_int; rng, save_positions = (false, false))
    jprob_int_ref = JumpProblem(rs_ref_int, u0_int, tspan_stoch, ps_int; rng, save_positions = (false, false))
    jsol_int = solve(jprob_int, SSAStepper(); seed, saveat = 1.0)
    jsol_int_ref = solve(jprob_int_ref, SSAStepper(); seed, saveat = 1.0)
    @test mean(jsol_int[:X1]) ≈ mean(jsol_int_ref[:X1]) atol = 1e-2 rtol = 1e-2
end

# Check that jump simulations (implemented with and without symbolic stoichiometries) yield simulations
# with identical mean number of species at the end of the simulations.
# Also checks that ODE simulations are identical.
let
    # Creates the models.
    sir = @reaction_network begin
        @parameters n::Int64 k::Int64
        i, S + n*I --> k*I
        r, n*I --> n*R
    end
    sir_ref = @reaction_network begin
        i, S + I --> 2I
        r, I --> R
    end

    ps = [:i => 1e-4, :r => 1e-2, :n => 1.0, :k => 2.0]
    ps_ref = [:i => 1e-4, :r => 1e-2]
    tspan = (0.0, 250.0) # tspan[2] is selected so that it is in the middle of the outbreak peak.
    u0 = [:S => 999.0, :I => 1.0, :R => 0.0]
    @test issetequal(unknowns(sir), unknowns(sir_ref))

    # Checks that ODE simulations are identical.
    oprob = ODEProblem(sir, u0, tspan, ps)
    oprob_ref = ODEProblem(sir_ref, u0, tspan, ps_ref)
    @test solve(oprob, Tsit5()) ≈ solve(oprob_ref, Tsit5())

    # Jumps. First ensemble problems for each systems is created.
    jprob = JumpProblem(sir, u0, tspan, ps; rng, save_positions = (false, false))
    jprob_ref = JumpProblem(sir_ref, u0, tspan, ps_ref; rng, save_positions = (false, false))
    eprob = EnsembleProblem(jprob)
    eprob_ref = EnsembleProblem(jprob_ref)

    # Simulates both ensemble problems. Checks that the distribution of values at the simulation ends is similar.
    sols = solve(eprob, SSAStepper(); trajectories = 10000)
    sols_ref = solve(eprob_ref, SSAStepper(); trajectories = 10000)
    end_vals = [[sol[s][end] for sol in sols.u] for s in species(sir)]
    end_vals_ref = [[sol[s][end] for sol in sols_ref.u] for s in species(sir_ref)]
    @test mean.(end_vals_ref) ≈ mean.(end_vals) atol=1e-1 rtol = 1e-1
    @test var.(end_vals_ref) ≈ var.(end_vals) atol=1e-1 rtol = 1e-1
end
