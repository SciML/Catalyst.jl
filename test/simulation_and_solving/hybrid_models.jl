# Fetch packages.
using Catalyst, JumpProcesses, ModelingToolkitBase, OrdinaryDiffEqTsit5, Statistics, Test, Random

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# test JumpProblem from ReactionSystem with different reaction types
let
    # Time-dependent rate -> VariableRateJump (needs ODE solver)
    rn = @reaction_network begin
        k*(1 + sin(t)), 0 --> A
    end
    jprob = JumpProblem(rn, [:A => 0], (0.0, 10.0), [:k => .5]; rng)
    sol = solve(jprob, Tsit5())
    @test sol(10.0; idxs = :A) > 0

    # Constant rate -> MassActionJump (can use SSAStepper)
    rn = @reaction_network begin
        k, 0 --> A
    end
    jprob = JumpProblem(rn, [:A => 0], (0.0, 10.0), [:k => .5]; rng)
    sol = solve(jprob, SSAStepper())
    @test sol(10.0; idxs = :A) > 0

    # Hybrid model with ODE equations and events
    rn = @reaction_network begin
        @parameters λ
        k*V, 0 --> A
        @equations D(V) ~ λ*V
        @continuous_events begin
            [V ~ 2.0] => [V ~ V/2, A ~ A/2]
        end
    end
    jprob = JumpProblem(rn, [:A => 0, :V => 1.0], (0.0, 10.0), [:k => 1.0, :λ => .4]; rng)
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
    jprob = JumpProblem(rn, u0map, tspan, pmap; save_positions = (false, false), rng)
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

    function affect!(m, o, ctx, integ)
        savevalues!(integ, true)
        terminate!(integ)
        return (;)
    end
    rn = @network_component begin
        β, X --> 0
        β, Y --> 0
        α, 0 --> Y
        (α * (1 + Y)), 0 --> X, [physical_scale = PhysicalScale.ODE]
    end
    t = Catalyst.default_t()
    cevents = [t ~ 0.2] => ModelingToolkitBase.ImperativeAffect(affect!)
    @named rn2 = ReactionSystem([], t; continuous_events = cevents)
    rn = complete(extend(rn2, rn))
    tspan = (0.0, 200.0)
    jprob = JumpProblem(rn, u0map, tspan, pmap; save_positions = (false, false), rng)
    Xsamp = 0.0
    Nsims = 4000
    for n in 1:Nsims
        sol = solve(jprob, Tsit5(); saveat = tspan[2], seed)
        @test SciMLBase.successful_retcode(sol)
        Xsamp += sol[1, end]
        seed += 1
    end
    Xsamp /= Nsims
    @test abs(Xsamp - Xf(0.2, p)) < 0.05 * Xf(0.2, p)
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
    jump_prob = JumpProblem(rn_jump, u0_jump, tspan, ps_jump; save_positions = (false, false), rng)
    hybrid_prob = JumpProblem(rn_hybrid, u0_hybrid, tspan, ps_hybrid; save_positions = (false, false), rng)

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
    prob = JumpProblem(rn_hybrid, u0, (0.0, 5.0), ps)

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
        @variables Ztot(τ) = 0.0  # Default needed for MTK v10+ callback compilation
        @observables Ztot ~ Z1 + Z2
        @discrete_events [1.0] => [Z1 ~ Pre(Z1) + 1.0]
        @continuous_events [Y ~ 1.0] => [Y ~ 5.0]
        @equations Δ(V) ~ Z1 + X^2 - V
        (p,d), 0 <--> X
        d, Y --> 0, [physical_scale = PhysicalScale.ODE]
        (k*X, k*Y), Z1 <--> Z2, ([physical_scale = PhysicalScale.ODE], [physical_scale = PhysicalScale.Jump])
    end

    # Simulates the model.
    u0 = [:Z1 => 1.0, :Z2 => 2.0, :V => 1.5]
    ps = [:p => 2.0, :d => 0.5, :k => 1.0, :X0 => 0.1, :Y0 => 4.0]
    prob = JumpProblem(rn, u0, (0.0, 10.0), ps)
    sol = solve(prob, Tsit5(); tstops = [1.0 - eps(), 1.0 + eps()])

    # Tests that the model contain the correct stuff.
    @test ModelingToolkitBase.getdescription(rn.X) == "A jump only species"
    @test ModelingToolkitBase.getdescription(rn.Y) == "An ODE only species"
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
    prob = JumpProblem(rn, u0, (0.0, 1.0), ps)
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Int32 (non-default), Y Float64 (default). Everything is converted to Float64.
    u0 = (:X => Int32(1), :Y => 1.0)
    prob = JumpProblem(rn, u0, (0.0, 1.0), ps)
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Int64 (default), Y Float32 (non-default). Everything is converted to Float64.
    u0 = (:X => 1, :Y => 1.0f0)
    prob = JumpProblem(rn, u0, (0.0, 1.0), ps)
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Int32 (non-default), Y Float32 (non-default). Everything is converted to Float64.
    u0 = (:X => Int32(1), :Y => 1.0f0)
    prob = JumpProblem(rn, u0, (0.0, 1.0), ps)
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Float32 (non-default float), Y Float32 (non-default). Everything is converted to Float32.
    u0 = (:X => 1.0f0, :Y => 1.0f0)
    ps = [:d => 0.5f0]
    prob = JumpProblem(rn, u0, (0.0, 1.0), ps)
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float32
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float32
end

### save_positions Tests ###

# Tests that save_positions keyword is correctly forwarded for pure SSA (DiscreteProblem) JumpProblems.
# With save_positions=(false, false), saveat should control output times exactly.
let
    # Pure SSA model with mass action jumps only.
    rn = @reaction_network begin
        (p, d), 0 <--> A
    end
    u0 = [:A => 10]
    ps = [:p => 1.0, :d => 0.1]
    tspan = (0.0, 10.0)
    times = 0.0:1.0:10.0

    # With save_positions=(false, false), saveat should give exactly the specified times.
    prob = JumpProblem(rn, u0, tspan, ps; save_positions = (false, false), rng)
    sol = solve(prob, SSAStepper(); saveat = times)
    @test length(sol.t) == length(times)
    @test sol.t ≈ collect(times)

    # With default save_positions (true, true), there will be extra points from jump times.
    prob_default = JumpProblem(rn, u0, tspan, ps; rng)
    sol_default = solve(prob_default, SSAStepper(); saveat = times)
    @test length(sol_default.t) > length(times)

    # Verify the aggregator has correct save_positions setting.
    disc_agg = prob.discrete_jump_aggregation
    @test disc_agg.save_positions == (false, false)
    disc_agg_default = prob_default.discrete_jump_aggregation
    @test disc_agg_default.save_positions == (true, true)
end

# Tests that save_positions keyword is correctly forwarded for hybrid (ODEProblem) JumpProblems
# with mass action jumps.
let
    # Hybrid model with ODE species and mass action jumps.
    rn = @reaction_network begin
        d, X --> 0, [physical_scale = PhysicalScale.ODE]
        (p, d), 0 <--> A
    end
    u0 = [:X => 10.0, :A => 5]
    ps = [:p => 1.0, :d => 0.1]
    tspan = (0.0, 10.0)
    times = 0.0:1.0:10.0

    # With save_positions=(false, false), saveat should control the number of output points.
    prob = JumpProblem(rn, u0, tspan, ps; save_positions = (false, false), rng)
    sol = solve(prob, Tsit5(); saveat = times)
    @test length(sol.t) == length(times)
    @test sol.t ≈ collect(times)

    # With default save_positions (true, true), there will be extra points from jump times.
    prob_default = JumpProblem(rn, u0, tspan, ps; rng)
    sol_default = solve(prob_default, Tsit5(); saveat = times)
    @test length(sol_default.t) > length(times)

    # Verify the aggregator has correct save_positions setting.
    disc_agg = prob.discrete_jump_aggregation
    @test disc_agg.save_positions == (false, false)
    disc_agg_default = prob_default.discrete_jump_aggregation
    @test disc_agg_default.save_positions == (true, true)
end

# Tests that save_positions keyword is correctly forwarded for hybrid JumpProblems
# with variable rate jumps (created via time-dependent rates).
let
    # Hybrid model with variable rate jump (time-dependent rate).
    rn = @reaction_network begin
        k * (1 + sin(t)), 0 --> A
        d, A --> 0
    end
    u0 = [:A => 0]
    ps = [:k => 0.5, :d => 0.1]
    tspan = (0.0, 10.0)
    times = 0.0:1.0:10.0

    # With save_positions=(false, false), the ContinuousCallback for VRJ should not save.
    prob = JumpProblem(rn, u0, tspan, ps; save_positions = (false, false), rng)
    sol = solve(prob, Tsit5(); saveat = times)

    # Should get approximately the expected number of times (may have small variations
    # due to event handling, but shouldn't have a point for every jump).
    @test length(sol.t) <= length(times) + 5  # Allow small tolerance for events

    # Verify the ContinuousCallback for VRJ has correct save_positions.
    cc = prob.jump_callback.continuous_callbacks[1]
    @test cc.save_positions == [false, false]

    # With default save_positions, ContinuousCallback should save.
    prob_default = JumpProblem(rn, u0, tspan, ps; rng)
    cc_default = prob_default.jump_callback.continuous_callbacks[1]
    @test cc_default.save_positions == [true, true]
end

# Tests that save_positions keyword works correctly for hybrid JumpProblems with
# a mix of mass action jumps and variable rate jumps.
let
    # Hybrid model with both mass action and variable rate jumps.
    rn = @reaction_network begin
        k * (1 + sin(t)), 0 --> A           # Variable rate jump
        d, A --> 0                           # Mass action jump
        (p, d2), 0 <--> B                    # Mass action jumps
    end
    u0 = [:A => 0, :B => 5]
    ps = [:k => 0.5, :d => 0.1, :p => 1.0, :d2 => 0.2]
    tspan = (0.0, 10.0)
    times = 0.0:1.0:10.0

    # With save_positions=(false, false), both mass action and VRJ should not save at jump times.
    prob = JumpProblem(rn, u0, tspan, ps; save_positions = (false, false), rng)
    sol = solve(prob, Tsit5(); saveat = times)

    # Should get approximately the expected number of times.
    @test length(sol.t) <= length(times) + 5

    # Verify both aggregator and ContinuousCallback have correct save_positions.
    disc_agg = prob.discrete_jump_aggregation
    @test disc_agg.save_positions == (false, false)
    cc = prob.jump_callback.continuous_callbacks[1]
    @test cc.save_positions == [false, false]
end

# Tests that saveat with save_positions=(false, false) gives correct solution values
# by comparing against interpolation-based access.
let
    rn = @reaction_network begin
        d, X --> 0, [physical_scale = PhysicalScale.ODE]
        (p, d), 0 <--> A
    end
    u0 = [:X => 10.0, :A => 5]
    ps = [:p => 1.0, :d => 0.1]
    tspan = (0.0, 10.0)
    times = 0.0:1.0:10.0

    # Create problem with save_positions=(false, false).
    prob = JumpProblem(rn, u0, tspan, ps; save_positions = (false, false), rng)

    # Solve with saveat.
    sol_saveat = solve(prob, Tsit5(); saveat = times)

    # Verify we can access solution values correctly.
    @test length(sol_saveat[:X]) == length(times)
    @test length(sol_saveat[:A]) == length(times)

    # The ODE species X should follow exponential decay (ignoring stochastic effects).
    # Just verify the values are reasonable.
    @test sol_saveat[:X][1] == 10.0
    @test sol_saveat[:X][end] < 10.0  # Should have decayed
end

### make_hybrid_model Tests ###

# Tests that make_hybrid_model produces correct structure for pure ODE, SDE, and Jump cases.
let
    rn = @reaction_network begin
        k1, S --> P
        k2, P --> S
    end

    # Pure ODE: should have ODE equations, no brownians, no jumps.
    sys_ode = make_hybrid_model(rn; default_scale = PhysicalScale.ODE)
    @test length(equations(sys_ode)) == 2
    @test isempty(ModelingToolkitBase.get_brownians(sys_ode))
    @test isempty(ModelingToolkitBase.get_jumps(sys_ode))

    # Pure SDE: should have ODE+noise equations with brownians, no jumps.
    sys_sde = make_hybrid_model(rn; default_scale = PhysicalScale.SDE)
    @test length(equations(sys_sde)) == 2
    @test length(ModelingToolkitBase.get_brownians(sys_sde)) == 2
    @test isempty(ModelingToolkitBase.get_jumps(sys_sde))

    # Pure Jump: should have no ODE equations, no brownians, only jumps.
    sys_jump = make_hybrid_model(rn; default_scale = PhysicalScale.Jump)
    @test isempty(equations(sys_jump))
    @test isempty(ModelingToolkitBase.get_brownians(sys_jump))
    @test length(ModelingToolkitBase.get_jumps(sys_jump)) == 2
end

# Tests that make_rre_ode and make_cle_sde thin wrappers produce equivalent results
# to direct make_hybrid_model calls.
let
    rn = @reaction_network begin
        k1, S --> P
        k2, P --> S
    end

    # make_rre_ode should produce same equations as make_hybrid_model with ODE override.
    sys_ode = make_rre_ode(rn)
    sys_hybrid_ode = make_hybrid_model(rn; default_scale = PhysicalScale.ODE)
    @test length(equations(sys_ode)) == length(equations(sys_hybrid_ode))
    for (eq1, eq2) in zip(equations(sys_ode), equations(sys_hybrid_ode))
        @test isequal(eq1.lhs, eq2.lhs)
        @test isequal(eq1.rhs, eq2.rhs)
    end

    # make_cle_sde should produce same structure as make_hybrid_model with SDE override.
    sys_sde = make_cle_sde(rn)
    sys_hybrid_sde = make_hybrid_model(rn; default_scale = PhysicalScale.SDE)
    @test length(equations(sys_sde)) == length(equations(sys_hybrid_sde))
    @test length(ModelingToolkitBase.get_brownians(sys_sde)) ==
          length(ModelingToolkitBase.get_brownians(sys_hybrid_sde))

    # make_sck_jump should produce same number of jumps as make_hybrid_model with Jump override.
    sys_jump = make_sck_jump(rn)
    sys_hybrid_jump = make_hybrid_model(rn; default_scale = PhysicalScale.Jump)
    @test length(ModelingToolkitBase.get_jumps(sys_jump)) ==
          length(ModelingToolkitBase.get_jumps(sys_hybrid_jump))
end

# Tests that the Brownian noise matrix extracted by mtkcompile matches assemble_diffusion.
let
    rn = @reaction_network begin
        k1, S --> P
        k2, P --> S
    end

    # Build via make_cle_sde (Brownian-based) and compile.
    sys_sde = make_cle_sde(rn)
    compiled = ModelingToolkitBase.mtkcompile(sys_sde)
    noise_matrix_brownian = ModelingToolkitBase.get_noise_eqs(compiled)

    # Build the old-style noise matrix via assemble_diffusion for comparison.
    flatrs = Catalyst.flatten(rn)
    ists, ispcs = Catalyst.get_indep_sts(flatrs, false)
    noise_matrix_old = Catalyst.assemble_diffusion(flatrs, ists, ispcs;
        combinatoric_ratelaws = true, remove_conserved = false, expand_catalyst_funs = true)

    # Both should be 2×2 matrices. Verify they are symbolically equivalent.
    @test size(noise_matrix_brownian) == size(noise_matrix_old)
    for i in axes(noise_matrix_brownian, 1), j in axes(noise_matrix_brownian, 2)
        @test Symbolics._iszero(Symbolics.simplify(noise_matrix_brownian[i, j] - noise_matrix_old[i, j]))
    end
end

# Tests make_hybrid_model with per-reaction physical_scales override.
let
    rn = @reaction_network begin
        k1, S --> P
        k2, P --> S
    end

    # Override reaction 1 to ODE, reaction 2 to Jump.
    sys = make_hybrid_model(rn; physical_scales = [1 => PhysicalScale.ODE, 2 => PhysicalScale.Jump])
    @test length(equations(sys)) == 2  # ODE equations present (for the ODE species)
    @test isempty(ModelingToolkitBase.get_brownians(sys))
    @test length(ModelingToolkitBase.get_jumps(sys)) == 1  # One jump (reaction 2)
end

# Tests that unresolved PhysicalScale.Auto throws an error.
let
    rn = @reaction_network begin
        k1, S --> P
        k2, P --> S
    end

    @test_throws ErrorException make_hybrid_model(rn; default_scale = PhysicalScale.Auto)
end

# Tests that remove_conserved with jump-scale reactions throws an error.
let
    rn = @reaction_network begin
        k1, S --> P
        k2, P --> S
    end

    @test_throws ArgumentError make_hybrid_model(rn;
        default_scale = PhysicalScale.Jump, remove_conserved = true)
end

# Tests that events pass through to the hybrid system.
let
    rn = @reaction_network begin
        @continuous_events begin
            [S ~ 50.0] => [S ~ 25.0]
        end
        @discrete_events [1.0] => [P ~ Pre(P) + 1.0]
        k1, S --> P
        k2, P --> S
    end

    sys = make_hybrid_model(rn; default_scale = PhysicalScale.ODE)
    @test length(ModelingToolkitBase.get_continuous_events(sys)) == 1
    @test length(ModelingToolkitBase.get_discrete_events(sys)) == 1
end

# Tests that make_hybrid_model works for a mixed ODE+SDE+Jump hybrid system.
let
    rn = @reaction_network begin
        k1, S --> P, [physical_scale = PhysicalScale.ODE]
        k2, P --> S, [physical_scale = PhysicalScale.SDE]
        k3, S --> 0, [physical_scale = PhysicalScale.Jump]
    end

    sys = make_hybrid_model(rn)
    @test length(equations(sys)) == 2  # ODE+SDE contribute drift equations
    @test length(ModelingToolkitBase.get_brownians(sys)) == 1  # One SDE reaction → one Brownian
    @test length(ModelingToolkitBase.get_jumps(sys)) == 1  # One jump reaction
end

# Tests that SDEProblem construction works end-to-end via the refactored path.
let
    rn = @reaction_network begin
        k1, S --> P
        k2, P --> S
    end

    prob = SDEProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 1.0), [:k1 => 1.0, :k2 => 0.5])
    @test prob.u0 == [100.0, 0.0]
    @test prob.tspan == (0.0, 1.0)
end

# Tests that make_sck_jump preserves VariableRateJump metadata.
# Uses independent species so VRJ classification doesn't propagate via dependency graph.
let
    rn = @reaction_network begin
        k1, A --> B, [physical_scale = PhysicalScale.VariableRateJump]
        k2, S --> P
    end

    sys = make_sck_jump(rn)
    jumps = ModelingToolkitBase.get_jumps(sys)
    # Reaction 1 should be a VariableRateJump (from metadata).
    # Reaction 2 should be a MassActionJump (independent species, default for make_sck_jump).
    has_vrj = any(j -> j isa JumpProcesses.VariableRateJump, jumps)
    has_maj = any(j -> j isa JumpProcesses.MassActionJump, jumps)
    @test has_vrj
    @test has_maj
end
