# Fetch packages.
using Catalyst, JumpProcesses, OrdinaryDiffEqTsit5, Statistics, Test, Random

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# test JumpInputs function auto problem selection
let
    rn = @reaction_network begin
        k*(1 + sin(t)), 0 --> A
    end
    jinput = JumpInputs(rn, [:A => 0], (0.0, 10.0), [:k => .5], remake_warn = false)
    @test jinput.prob isa ODEProblem
    jprob = JumpProblem(jinput; rng)
    sol = solve(jprob, Tsit5())
    @test sol(10.0; idxs = :A) > 0

    rn = @reaction_network begin
        k, 0 --> A
    end
    jinput = JumpInputs(rn, [:A => 0], (0.0, 10.0), [:k => .5])
    @test jinput.prob isa DiscreteProblem
    jprob = JumpProblem(jinput; rng)
    sol = solve(jprob)
    @test sol(10.0; idxs = :A) > 0

    rn = @reaction_network begin
        @parameters λ
        k*V, 0 --> A
        @equations D(V) ~ λ*V
        @continuous_events begin
            [V ~ 2.0] => [V ~ V/2, A ~ A/2]
        end
    end
    jinput = JumpInputs(rn, [:A => 0, :V => 1.0], (0.0, 10.0), [:k => 1.0, :λ => .4],
        remake_warn = false)
    @test jinput.prob isa ODEProblem
    jprob = JumpProblem(jinput; rng)
    sol = solve(jprob, Tsit5())
end

# solution correctness tests
let
    seed = 1111
    Random.seed!(rng, seed)
    rn = @reaction_network begin
        β, X --> 0
        β, Y --> 0
        α, 0 --> Y
        (α * (1 + Y)), 0 --> X, [physical_scale = PhysicalScale.ODE]
    end
    p = (α = 6.0, β = 2.0, X₀ = 2.0, Y₀ = 1.0)
    u0map = [:X => p.X₀, :Y => p.Y₀]
    pmap = [:α => p.α, :β => p.β]
    tspan = (0.0, 20.0)
    jinputs = JumpInputs(rn, u0map, tspan, pmap; save_positions = (false, false))
    jprob = JumpProblem(jinputs; rng, save_positions = (false, false))
    times = range(0.0, tspan[2], length = 100)
    Nsims = 4000
    Xv = zeros(length(times))
    Yv = zeros(length(times))
    for n in 1:Nsims
        sol = solve(jprob, Tsit5(); saveat = times, seed)
        Xv .+= sol[:X]
        Yv .+= sol[:Y]
        seed += 1
    end
    Xv ./= Nsims
    Yv ./= Nsims

    function Yf(t, p)
        local α, β, X₀, Y₀ = p
        return (α / β) + (Y₀ - α / β) * exp(-β * t)
    end
    function Xf(t, p)
        local α, β, X₀, Y₀ = p
        return (α / β) + (α^2 / β^2) + α * (Y₀ - α / β) * t * exp(-β * t) +
               (X₀ - α / β - α^2 / β^2) * exp(-β * t)
    end
    Xact = [Xf(t, p) for t in times]
    Yact = [Yf(t, p) for t in times]
    @test all(abs.(Xv .- Xact) .<= 0.05 .* Xv)
    @test all(abs.(Yv .- Yact) .<= 0.05 .* Yv)

    function affect!(integ, u, p, ctx)
        savevalues!(integ, true)
        terminate!(integ)
        nothing
    end
    rn = @network_component begin
        β, X --> 0
        β, Y --> 0
        α, 0 --> Y
        (α * (1 + Y)), 0 --> X, [physical_scale = PhysicalScale.ODE]
    end
    t = Catalyst.default_t()
    cevents = [t ~ 0.2] => (affect!, [], [], [], nothing)
    @named rn2 = ReactionSystem([], t; continuous_events = cevents)
    rn = complete(extend(rn2, rn))
    tspan = (0.0, 200.0)
    jinputs = JumpInputs(rn, u0map, tspan, pmap; save_positions = (false, false))
    jprob = JumpProblem(jinputs; rng, save_positions = (false, false))
    Xsamp = 0.0
    Nsims = 4000
    for n in 1:Nsims
        sol = solve(jprob, Tsit5(); saveat = tspan[2], seed)
        @test SciMLBase.successful_retcode(sol)
        Xsamp += sol[1, end]
        seed += 1
    end
    Xsamp /= Nsims
    @test_broken abs(Xsamp - Xf(0.2, p)) < 0.05 * Xf(0.2, p)
end

# Checks that a disjoint hybrid model (i.e. where the Jump and ODE parts do not interact) gives the
# same results as simulating the two parts separately.
let
    # Creates the disjoint ODE/Jump models, and the hybrid model combining both.
    rn_ode = @reaction_network begin
        A, ∅ → X
        1, 2X + Y → 3X
        B, X → Y
        1, X → ∅
    end
    rn_jump = @reaction_network begin
        (p,d), 0 <--> V
        (k1,k2), V + W <--> VW
    end
    rn_hybrid = @reaction_network begin
        A, ∅ → X, [physical_scale = PhysicalScale.ODE]
        1, 2X + Y → 3X, [physical_scale = PhysicalScale.ODE]
        B, X → Y, [physical_scale = PhysicalScale.ODE]
        1, X → ∅, [physical_scale = PhysicalScale.ODE]
        (p,d), 0 <--> V
        (k1,k2), V + W <--> VW
    end

    # Sets simulation conditions and creates problems corresponding to the different models.
    u0_ode = [:X => 1.0, :Y => 0.0]
    u0_jump = [:V => 1, :W => 20, :VW => 0.0]
    u0_hybrid = (u0_ode..., u0_jump...)
    ps_ode = [:A => 1.0, :B => 4.0]
    ps_jump = [:p => 2.0, :d => 0.1, :k1 => 0.2, :k2 => 1.0]
    ps_hybrid = [ps_ode; ps_jump]
    tspan = (0.0, 10000.0)
    ode_prob = ODEProblem(rn_ode, u0_ode, tspan, ps_ode)
    jump_prob = JumpProblem(JumpInputs(rn_jump, u0_jump, tspan, ps_jump; remake_warn = false); save_positions = (false,false), rng)
    hybrid_prob = JumpProblem(JumpInputs(rn_hybrid, u0_hybrid, tspan, ps_hybrid; remake_warn = false); save_positions = (false,false), rng)

    # Performs simulations. Checks that ODE parts are identical. Check that jump parts have similar statistics.
    ode_sol = solve(ode_prob, Tsit5(); saveat = 1.0, abstol = 1e-10, reltol = 1e-10)
    jump_sol = solve(jump_prob, SSAStepper(); saveat = 1.0)
    hybrid_sol = solve(hybrid_prob, Tsit5(); saveat = 1.0, abstol = 1e-10, reltol = 1e-10)
    @test ode_sol[:Y] ≈ hybrid_sol[:Y] atol = 1e-4 rtol = 1e-4
    @test mean(jump_sol[:V]) ≈ mean(hybrid_sol[:V]) atol = 1e-1 rtol = 1e-1
    @test mean(jump_sol[:W]) ≈ mean(hybrid_sol[:W]) atol = 1e-1 rtol = 1e-1
    @test mean(jump_sol[:VW]) ≈ mean(hybrid_sol[:VW]) atol = 1e-1 rtol = 1e-1
end

### Other Tests ###

# Checks that the full pipeline of symbolic accessing/updating problems/integrators/solutions
# of hybrid models work.
# CURRENTLY NOT WORKING DUE TO `remake` ISSUE WITH HYBRID MODELS.
let
    # Creates model and initial problem (with defaults for one species and one parameter).
    rn_hybrid = @reaction_network begin
        @species X(t) = 1.0
        @parameters p = 1.0
        p, 0 --> X, [physical_scale = PhysicalScale.ODE]
        (kB,kD), 2X <--> X2
        k, X2 --> Y2, [physical_scale = PhysicalScale.ODE]
        (kB,kD), 2Y <--> Y2
        d, Y --> 0, [physical_scale = PhysicalScale.ODE]
    end
    u0 = [:X2 => 0.0, :Y => 0.0, :Y2 => 0.0]
    ps = [:kB => 2.0, :kD => 0.4, :k => 1.0, :d => 0.5]
    prob = JumpProblem(JumpInputs(rn_hybrid, u0, (0.0, 5.0), ps; remake_warn = false))

    # Check problem content, updates problem with remake, checks again.
    @test prob[:X] == 1.0
    @test prob[:X2] == 0.0
    @test prob[:Y] == 0.0
    @test prob[:Y2] == 0.0
    @test prob.ps[:p] == 1.0
    @test prob.ps[:kB] == 2.0
    @test prob.ps[:kD] == 0.4
    @test prob.ps[:k] == 1.0
    @test prob.ps[:d] == 0.5

    # AS OF MTK 9.72 THESE ARE BROKEN AND GIVE ERRORS
    # prob = remake(prob; u0 = [:X => 3.0, :X2 => 2.0, :Y => 1.0], p = [:p => 5.0, :dD => 0.8, :k => 2.0])
    # @test prob[:X] == 3.0
    # @test prob[:X2] == 2.0
    # @test prob[:Y] == 1.0
    # @test prob[:Y2] == 0.0
    # @test prob.ps[:p] == 5.0
    # @test prob.ps[:kB] == 2.0
    # @test prob.ps[:kD] == 0.8
    # @test prob.ps[:k] == 2.0
    # @test prob.ps[:d] == 0.5

    # # Creates, checks, updates, and checks an integrator.
    # int = init(prob, Tsit5())
    # @test int[:X] == 3.0
    # @test int[:X2] == 2.0
    # @test int[:Y] == 1.0
    # @test int[:Y2] == 0.0
    # @test int.ps[:p] == 5.0
    # @test int.ps[:kB] == 2.0
    # @test int.ps[:kD] == 0.8
    # @test int.ps[:k] == 2.0
    # @test int.ps[:d] == 0.5
    # @test int[:X] = 4.0
    # @test int[:X2] = 3.0
    # @test int[:Y] = 2.0
    # @test int.ps[:p] = 6.0
    # @test int.ps[:kB] = 3.0
    # @test int.ps[:kD] = 0.6
    # @test int.ps[:k] = 3.0
    # @test int[:X] == 4.0
    # @test int[:X2] == 3.0
    # @test int[:Y] == 2.0
    # @test int[:Y2] == 0.0
    # @test int.ps[:p] == 6.0
    # @test int.ps[:kB] == 3.0
    # @test int.ps[:kD] == 0.5
    # @test int.ps[:k] == 3.0
    # @test int.ps[:d] == 0.5

    # # Solves and checks values of solution.
    # sol = solve(prob, Tsit5(); maxiters = 1, verbose = false)
    # @test sol[:X][1] == 3.0
    # @test sol[:X2][1] == 2.0
    # @test sol[:Y][1] == 1.0
    # @test sol[:Y2][1] == 0.0
    # @test sol.ps[:p] == 5.0
    # @test sol.ps[:kB] == 2.0
    # @test sol.ps[:kD] == 0.8
    # @test sol.ps[:k] == 2.0
    # @test sol.ps[:d] == 0.5
end

# Checks that various model options (observables, events, defaults and metadata, differential equations,
# non-default_iv) works for hybrid models.
let
    # Creates the model (X species is pure jump, Y is pure ODE, and  Z1,Z2 are mixed).
    # Hybrid species have non-constant rates containing the two other species.
    rn = @reaction_network begin
        @ivs τ
        @differentials Δ = Differential(τ)
        @parameters X0 Y0
        @species X(τ) = X0 [description = "A jump only species"] Y(τ) = Y0 [description = "An ODE only species"]
        @observables Ztot ~ Z1 + Z2
        @discrete_events [1.0] => [Z1 ~ Z1 + 1.0]
        @continuous_events [Y ~ 1.0] => [Y ~ 5.0]
        @equations Δ(V) ~ Z1 + X^2 - V
        (p,d), 0 <--> X
        d, Y --> 0, [physical_scale = PhysicalScale.ODE]
        (k*X, k*Y), Z1 <--> Z2, ([physical_scale = PhysicalScale.ODE], [physical_scale = PhysicalScale.Jump])
    end

    # Simulates the model.
    u0 = [:Z1 => 1.0, :Z2 => 2.0, :V => 1.5]
    ps = [:p => 2.0, :d => 0.5, :k => 1.0, :X0 => 0.1, :Y0 => 4.0]
    prob = JumpProblem(JumpInputs(rn, u0, (0.0, 10.0), ps; remake_warn = false))
    sol = solve(prob, Tsit5(); tstops = [1.0 - eps(), 1.0 + eps()])

    # Tests that the model contain the correct stuff.
    @test ModelingToolkit.getdescription(rn.X) == "A jump only species"
    @test ModelingToolkit.getdescription(rn.Y) == "An ODE only species"
    @test sol[:X][1] == sol.ps[:X0] == 0.1
    @test sol[:Y][1] == sol.ps[:Y0] == 4.0
    @test sol[:V][1] == 1.5
    @test sol[:Ztot] ≈ sol[rn.Z1 + rn.Z2]
    @test minimum(sol[:Y]) ≈ 1.0
    @test maximum(sol[:Y]) ≈ 5.0 atol = 1e-1 rtol = 1e-1
    @test all(isequal([rn.τ], Symbolics.arguments(Symbolics.value(u))) for u in unknowns(rn))
    @test sol(1.0 - eps(); idxs = :Z1) + 1 ≈ sol(1.0 + eps(); idxs = :Z1)
end

# Checks the types of species when various combinations of default/non-default types are used.
let
    # Creates model and parameter set. Model have one pure ODE and one pure Jump species.
    rn = @reaction_network begin
        d, X --> 0
        d, Y --> 0, [physical_scale = PhysicalScale.ODE]
    end
    ps = [:d => 0.5]

    # Checks case: X Int64 (default), Y Float64 (default). Everything is converted to Float64.
    u0 = (:X => 1, :Y => 1.0)
    prob = JumpProblem(JumpInputs(rn, u0, (0.0, 1.0), ps; remake_warn = false))
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Int32 (non-default), Y Float64 (default). Everything is converted to Float64.
    u0 = (:X => Int32(1), :Y => 1.0)
    prob = JumpProblem(JumpInputs(rn, u0, (0.0, 1.0), ps; remake_warn = false))
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Int64 (default), Y Float32 (non-default). Everything is converted to Float64.
    u0 = (:X => 1, :Y => 1.0f0)
    prob = JumpProblem(JumpInputs(rn, u0, (0.0, 1.0), ps; remake_warn = false))
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Int32 (non-default), Y Float32 (non-default). Everything is converted to Float64.
    u0 = (:X => Int32(1), :Y => 1.0f0)
    prob = JumpProblem(JumpInputs(rn, u0, (0.0, 1.0), ps; remake_warn = false))
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Float32 (non-default float), Y Float32 (non-default). Everything is converted to Float32.
    u0 = (:X => 1.0f0, :Y => 1.0f0)
    ps = [:d => 0.5f0]
    prob = JumpProblem(JumpInputs(rn, u0, (0.0, 1.0), ps; remake_warn = false))
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float32
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float32
end
