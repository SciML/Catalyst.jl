### Prepares Tests ###

# Fetch packages.
using Catalyst, JumpProcesses, OrdinaryDiffEq, StochasticDiffEq, Statistics, Test
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
    rs1 = ReactionSystem(rxs1, t; name = :rs)

    rxs2 = [
        @reaction p, 0 --> $n1*X
        @reaction k, n2*X --> n3*Y
        @reaction $d, Y --> 0
    ]
    rs2 = ReactionSystem(rxs2, t; name = :rs)

    rs3 = @reaction_network rs begin
        @parameters d::Float64 n1::Int64
        p, 0 --> n1*X
        k, n2*X --> n3*Y
        d, Y --> 0
    end

    @test rs1 == rs2 == rs3
    @test issetequal(unknowns(rs1), [X, Y])
    @test issetequal(parameters(rs1), [p, k, d, n1, n2, n3])
    @test unwrap(d) isa SymbolicUtils.BasicSymbolic{Float64}
    @test unwrap(n1) isa SymbolicUtils.BasicSymbolic{Int64}
end

# Declares a network, parameter values, and initial conditions, to be used for the next couple of tests.
begin
    # Prepares parameters and species and their values.
    @parameters k α::Int64
    @species A(t), B(t), C(t), D(t)
    u0_1 = (A => 3.0, B => 2.0, C => 3.0, D => 5.0)
    ps_1 = (k => 5.0, α => 2)
    u0_2 = [u[2] for u in u0_1]
    ps_2 = [p[2] for p in ps_1]
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
    for i in 1:2
        catalyst_jsys = convert(JumpSystem, rs)
        unknownoid = Dict(unknown => i for (i, unknown) in enumerate(unknowns(catalyst_jsys)))
        catalyst_vrj = ModelingToolkit.assemble_vrj(catalyst_jsys, equations(catalyst_jsys)[i], unknownoid)
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
    rs_int = @reaction_network begin
        @parameters n1::Int64 n2::Int64
        p, 0 -->X
        k, n1*X --> n2*Y
        d, Y --> 0
    end
    rs_dec = @reaction_network begin
        @parameters n1::Float64 n2::Float64
        p, 0 -->X
        k, n1*X --> n2*Y
        d, Y --> 0
    end
    rs_ref_int = @reaction_network begin
        @parameters n1::Int64 n2::Int64
        p, 0 -->X
        k, 3*X --> 4*Y
        d, Y --> 0
    end
    rs_ref_dec = @reaction_network begin
        @parameters n1::Float64 n2::Float64
        p, 0 -->X
        k, 0.5*X --> 2.5*Y
        d, Y --> 0
    end

    u0 = [:X => 1, :Y => 1]
    ps_int = [:p => 1.0, :n1 => 3, :n2 => 4, :d => 0.5]
    ps_dec = [:p => 1.0, :n1 => 0.5, :n2 => 2.5, :d => 0.5]
    tspan = (0.0, 10.0)

    # Test ODE simulations with integer coefficients.
    oprob_int = ODEProblem(rs_int, u0, tspan, ps_int)
    oprob_int_ref = ODEProblem(rs_ref_int, u0, tspan, ps_int)
    @test solve(oprob_int, Tsit5()) == solve(oprob_int_ref, Tsit5())

    # Test ODE simulations with decimal coefficients.
    oprob_dec = ODEProblem(rs_dec, u0, tspan, ps_dec; combinatoric_ratelaws = false)
    oprob_dec_ref = ODEProblem(rs_ref_dec, u0, tspan, ps_dec; combinatoric_ratelaws = false)
    @test solve(oprob_dec, Tsit5()) == solve(oprob_dec_ref, Tsit5())

    # Test SDE simulations with integer coefficients.
    sprob_int = SDEProblem(rs_int, u0, tspan, ps_int)
    sprob_int_ref = SDEProblem(rs_ref_int, u0, tspan, ps_int)
    @test solve(sprob_int, ImplicitEM(); seed) == solve(sprob_int_ref, ImplicitEM(); seed)

    # Test SDE simulations with decimal coefficients.
    sprob_dec = SDEProblem(rs_dec, u0, tspan, ps_dec; combinatoric_ratelaws = false)
    sprob_dec_ref = SDEProblem(rs_ref_dec, u0, tspan, ps_dec; combinatoric_ratelaws = false)
    @test solve(sprob_dec, ImplicitEM(); seed) == solve(sprob_dec_ref, ImplicitEM(); seed)

    # Tests jump simulations with integer coefficients.
    dprob_int = DiscreteProblem(rs_int, u0, tspan, ps_int)
    dprob_int_ref = DiscreteProblem(rs_ref_int, u0, tspan, ps_int)
    jprob_int = JumpProblem(rs_int, dprob_int, Direct(); rng)
    jprob_int_ref = JumpProblem(rs_ref_int, dprob_int_ref, Direct(); rng)
    @test solve(jprob_int, SSAStepper(); seed) == solve(jprob_int_ref, SSAStepper(); seed)
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
    dprob = DiscreteProblem(sir, u0, tspan, ps) 
    dprob_ref = DiscreteProblem(sir_ref, u0, tspan, ps_ref) 
    jprob = JumpProblem(sir, dprob, Direct(); rng, save_positions = (false, false)) 
    jprob_ref = JumpProblem(sir_ref, dprob_ref, Direct(); rng, save_positions = (false, false)) 
    eprob = EnsembleProblem(jprob)
    eprob_ref = EnsembleProblem(jprob_ref)

    # Simulates both ensemble problems. Checks that the distribution of values at the simulation ends is similar.
    sols = solve(eprob, SSAStepper(); trajectories = 10000)
    sols_ref = solve(eprob_ref, SSAStepper(); trajectories = 10000)
    end_vals = [[sol[s][end] for sol in sols.u] for s in species(sir)]
    end_vals_ref = [[sol[s][end] for sol in sols_ref.u] for s in species(sir_ref)]
    @test mean.(end_vals_ref) ≈ mean.(end_vals) atol=1e-2 rtol = 1e-2
    @test var.(end_vals_ref) ≈ var.(end_vals) atol=1e-1 rtol = 1e-1
end