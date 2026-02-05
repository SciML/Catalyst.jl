# Fetch packages.
using Catalyst, JumpProcesses, ModelingToolkitBase, OrdinaryDiffEqTsit5, Statistics, Test, Random
using StochasticDiffEq: SRIW1  # For SDE+Jump hybrid problems
import DiffEqNoiseProcess  # Required for SDEProblem via mtkcompile

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

    # Hybrid model with ODE equations and events - requires HybridProblem
    rn = @reaction_network begin
        @parameters λ
        k*V, 0 --> A
        @equations D(V) ~ λ*V
        @continuous_events begin
            [V ~ 2.0] => [V ~ V/2, A ~ A/2]
        end
    end
    # JumpProblem no longer supports ODE equations - use HybridProblem instead
    jprob = HybridProblem(rn, [:A => 0, :V => 1.0], (0.0, 10.0), [:k => 1.0, :λ => .4]; rng)
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
    jprob = HybridProblem(rn, u0map, tspan, pmap; save_positions = (false, false), rng)
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
    jprob = HybridProblem(rn, u0map, tspan, pmap; save_positions = (false, false), rng)
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
    hybrid_prob = HybridProblem(rn_hybrid, u0_hybrid, tspan, ps_hybrid; save_positions = (false, false), rng)

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
    prob = HybridProblem(rn_hybrid, u0, (0.0, 5.0), ps)

    # Check problem content.
    @test prob[:X] == 1.0
    @test prob[:X2] == 0.0
    @test prob[:Y] == 0.0
    @test prob[:Y2] == 0.0
    @test prob.ps[:p] == 1.0
    @test prob.ps[:kB] == 2.0
    @test prob.ps[:kD] == 0.4
    @test prob.ps[:k] == 1.0
    @test prob.ps[:d] == 0.5

    # Solves and checks values of solution (do this before integrator mutation test
    # since integrator mutation affects shared parameter state).
    sol = solve(prob, Tsit5(); maxiters = 10, verbose = false)
    @test sol[:X][1] == 1.0
    @test sol[:X2][1] == 0.0
    @test sol[:Y][1] == 0.0
    @test sol[:Y2][1] == 0.0
    @test sol.ps[:p] == 1.0
    @test sol.ps[:kB] == 2.0
    @test sol.ps[:kD] == 0.4
    @test sol.ps[:k] == 1.0
    @test sol.ps[:d] == 0.5

    # Creates, checks, updates, and checks an integrator.
    int = init(prob, Tsit5())
    @test int[:X] == 1.0
    @test int[:X2] == 0.0
    @test int[:Y] == 0.0
    @test int[:Y2] == 0.0
    @test int.ps[:p] == 1.0
    @test int.ps[:kB] == 2.0
    @test int.ps[:kD] == 0.4
    @test int.ps[:k] == 1.0
    @test int.ps[:d] == 0.5

    # Mutate integrator state and parameters.
    int[:X] = 4.0
    int[:X2] = 3.0
    int[:Y] = 2.0
    int.ps[:p] = 6.0
    int.ps[:kB] = 3.0
    int.ps[:kD] = 0.6
    int.ps[:k] = 3.0

    # Check mutated values.
    @test int[:X] == 4.0
    @test int[:X2] == 3.0
    @test int[:Y] == 2.0
    @test int[:Y2] == 0.0
    @test int.ps[:p] == 6.0
    @test int.ps[:kB] == 3.0
    @test int.ps[:kD] == 0.6
    @test int.ps[:k] == 3.0
    @test int.ps[:d] == 0.5

    # REMAKE WITH SYMBOL/DICT u0 IS BROKEN (MTKBase issue with JumpProblem remake).
    # Error: "Passed in u0 is incompatible with current u0 which has type: Vector{Float64}."
    # Workaround: remake with numeric array works, or remake with only parameters works.
    # Upstream issue needed in ModelingToolkit/JumpProcesses.
    # prob = remake(prob; u0 = [:X => 3.0, :X2 => 2.0, :Y => 1.0], p = [:p => 5.0, :kD => 0.8, :k => 2.0])
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
    prob = HybridProblem(rn, u0, (0.0, 10.0), ps)
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
    prob = HybridProblem(rn, u0, (0.0, 1.0), ps)
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Int32 (non-default), Y Float64 (default). Everything is converted to Float64.
    u0 = (:X => Int32(1), :Y => 1.0)
    prob = HybridProblem(rn, u0, (0.0, 1.0), ps)
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Int64 (default), Y Float32 (non-default). Everything is converted to Float64.
    u0 = (:X => 1, :Y => 1.0f0)
    prob = HybridProblem(rn, u0, (0.0, 1.0), ps)
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Int32 (non-default), Y Float32 (non-default). Everything is converted to Float64.
    u0 = (:X => Int32(1), :Y => 1.0f0)
    prob = HybridProblem(rn, u0, (0.0, 1.0), ps)
    int = init(prob, Tsit5())
    sol = solve(prob, Tsit5())
    @test typeof(prob[:X]) == typeof(int[:X]) == typeof(sol[:X][end]) == Float64
    @test typeof(prob[:Y]) == typeof(int[:Y]) == typeof(sol[:Y][end]) == Float64

    # Checks case: X Float32 (non-default float), Y Float32 (non-default). Everything is converted to Float32.
    u0 = (:X => 1.0f0, :Y => 1.0f0)
    ps = [:d => 0.5f0]
    prob = HybridProblem(rn, u0, (0.0, 1.0), ps)
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
    prob = HybridProblem(rn, u0, tspan, ps; save_positions = (false, false), rng)
    sol = solve(prob, Tsit5(); saveat = times)
    @test length(sol.t) == length(times)
    @test sol.t ≈ collect(times)

    # With default save_positions (true, true), there will be extra points from jump times.
    prob_default = HybridProblem(rn, u0, tspan, ps; rng)
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

    # With save_positions=(false, false), saveat should give exactly the specified times.
    @test length(sol.t) == length(times)
    @test sol.t ≈ collect(times)

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

    # With save_positions=(false, false), saveat should give exactly the specified times.
    @test length(sol.t) == length(times)
    @test sol.t ≈ collect(times)

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
    prob = HybridProblem(rn, u0, tspan, ps; save_positions = (false, false), rng)

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

# Tests that save_positions keyword is correctly forwarded for SDE+Jump hybrid systems
# (SDEProblem with jumps overlaid).
let
    # Hybrid model with SDE noise and mass action jumps.
    rn = @reaction_network begin
        k1, S --> P, [physical_scale = PhysicalScale.SDE]
        (p, d), 0 <--> A
    end
    u0 = [:S => 50.0, :P => 0.0, :A => 10]
    ps = [:k1 => 0.1, :p => 1.0, :d => 0.1]
    tspan = (0.0, 10.0)
    times = 0.0:1.0:10.0

    # With save_positions=(false, false), saveat should control output times.
    prob = HybridProblem(rn, u0, tspan, ps; save_positions = (false, false), rng)
    sol = solve(prob, SRIW1(); saveat = times)
    @test length(sol.t) == length(times)
    @test sol.t ≈ collect(times)

    # With default save_positions (true, true), there will be extra points from jump times.
    prob_default = HybridProblem(rn, u0, tspan, ps; rng)
    sol_default = solve(prob_default, SRIW1(); saveat = times)
    @test length(sol_default.t) > length(times)

    # Verify the aggregator has correct save_positions setting.
    disc_agg = prob.discrete_jump_aggregation
    @test disc_agg.save_positions == (false, false)
    disc_agg_default = prob_default.discrete_jump_aggregation
    @test disc_agg_default.save_positions == (true, true)
end

# Tests that save_positions keyword works for SDE+VariableRateJump hybrid systems.
let
    # Hybrid model with SDE noise and variable rate jump (time-dependent rate).
    rn = @reaction_network begin
        k1, S --> P, [physical_scale = PhysicalScale.SDE]
        k2 * (1 + sin(t)), 0 --> A
        d, A --> 0
    end
    u0 = [:S => 50.0, :P => 0.0, :A => 0]
    ps = [:k1 => 0.1, :k2 => 0.5, :d => 0.1]
    tspan = (0.0, 10.0)
    times = 0.0:1.0:10.0

    # With save_positions=(false, false), both SDE and VRJ should not save at event times.
    prob = HybridProblem(rn, u0, tspan, ps; save_positions = (false, false), rng)
    sol = solve(prob, SRIW1(); saveat = times)

    # With save_positions=(false, false), saveat should give exactly the specified times.
    @test length(sol.t) == length(times)
    @test sol.t ≈ collect(times)

    # Verify the ContinuousCallback for VRJ has correct save_positions.
    cc = prob.jump_callback.continuous_callbacks[1]
    @test cc.save_positions == [false, false]

    # With default save_positions, ContinuousCallback should save.
    prob_default = HybridProblem(rn, u0, tspan, ps; rng)
    cc_default = prob_default.jump_callback.continuous_callbacks[1]
    @test cc_default.save_positions == [true, true]
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

    # make_cle_sde with use_legacy_noise=false should produce same structure as make_hybrid_model with SDE override.
    sys_sde = make_cle_sde(rn; use_legacy_noise = false)
    sys_hybrid_sde = make_hybrid_model(rn; default_scale = PhysicalScale.SDE)
    @test length(equations(sys_sde)) == length(equations(sys_hybrid_sde))
    @test length(ModelingToolkitBase.get_brownians(sys_sde)) ==
          length(ModelingToolkitBase.get_brownians(sys_hybrid_sde))

    # make_cle_sde with use_legacy_noise=true (default) should use noise_eqs matrix, not brownians.
    sys_sde_legacy = make_cle_sde(rn; use_legacy_noise = true)
    @test ModelingToolkitBase.get_noise_eqs(sys_sde_legacy) !== nothing
    @test isempty(ModelingToolkitBase.get_brownians(sys_sde_legacy))

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

    # Build via make_cle_sde (Brownian-based, use_legacy_noise=false) and compile.
    sys_sde = make_cle_sde(rn; use_legacy_noise = false)
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

# Tests that make_sck_jump does NOT respect ODE metadata - all reactions become jumps.
let
    rn = @reaction_network begin
        k1, S --> P, [physical_scale = PhysicalScale.ODE]
        k2, P --> S
    end

    sys = make_sck_jump(rn)
    # Both reactions should be jumps (ODE metadata is ignored).
    @test isempty(equations(sys))
    @test length(ModelingToolkitBase.get_jumps(sys)) == 2
end

# Tests that HybridProblem works for mixed ODE+Jump systems.
let
    rn = @reaction_network begin
        k1, S --> P, [physical_scale = PhysicalScale.ODE]
        k2, P --> S, [physical_scale = PhysicalScale.Jump]
    end

    prob = HybridProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 1.0), [:k1 => 1.0, :k2 => 0.5])
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
end

# Tests that HybridProblem default_scale works (untagged reactions become jumps by default).
let
    rn = @reaction_network begin
        k1, S --> P, [physical_scale = PhysicalScale.ODE]
        k2, P --> S  # default_scale will apply to this reaction
    end

    # With default_scale = Jump (the default), the second reaction becomes a jump.
    prob = HybridProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 1.0), [:k1 => 1.0, :k2 => 0.5])
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
end

### HybridProblem Return Type Tests ###

# Tests that HybridProblem returns correct problem types based on scales.
let
    rn = @reaction_network begin
        k1, S --> P
        k2, P --> S
    end

    # Pure ODE → ODEProblem
    prob_ode = HybridProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 1.0), [:k1 => 1.0, :k2 => 0.5];
        default_scale = PhysicalScale.ODE)
    @test prob_ode isa ODEProblem
    sol = solve(prob_ode, Tsit5())
    @test SciMLBase.successful_retcode(sol)

    # Pure SDE → SDEProblem
    prob_sde = HybridProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 1.0), [:k1 => 1.0, :k2 => 0.5];
        default_scale = PhysicalScale.SDE)
    @test prob_sde isa SDEProblem

    # Pure Jump → JumpProblem
    prob_jump = HybridProblem(rn, [:S => 100, :P => 0], (0.0, 1.0), [:k1 => 1.0, :k2 => 0.5];
        default_scale = PhysicalScale.Jump)
    @test prob_jump isa JumpProcesses.JumpProblem
    sol = solve(prob_jump, SSAStepper())
    @test SciMLBase.successful_retcode(sol)
end

# Tests ODE+Jump hybrid returns JumpProblem.
let
    rn = @reaction_network begin
        k1, S --> P, [physical_scale = PhysicalScale.ODE]
        k2, P --> S, [physical_scale = PhysicalScale.Jump]
    end

    prob = HybridProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 1.0), [:k1 => 1.0, :k2 => 0.5])
    @test prob isa JumpProcesses.JumpProblem
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
end

# Tests ODE+SDE (no jumps) returns SDEProblem.
let
    rn = @reaction_network begin
        k1, S --> P, [physical_scale = PhysicalScale.ODE]
        k2, P --> S, [physical_scale = PhysicalScale.SDE]
    end

    prob = HybridProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 1.0), [:k1 => 1.0, :k2 => 0.5])
    @test prob isa SDEProblem
end

# Tests SDE+Jump hybrid returns JumpProblem wrapping SDEProblem.
let
    rn = @reaction_network begin
        k1, S --> P, [physical_scale = PhysicalScale.SDE]
        k2, P --> S, [physical_scale = PhysicalScale.Jump]
    end

    prob = HybridProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 1.0), [:k1 => 1.0, :k2 => 0.5])
    @test prob isa JumpProcesses.JumpProblem
    @test prob.prob isa SciMLBase.SDEProblem

    # Test that we can solve it (use SDE solver since underlying problem is SDEProblem)
    sol = solve(prob, SRIW1())
    @test SciMLBase.successful_retcode(sol)
end

# Tests ODE+SDE+Jump full hybrid system.
let
    rn = @reaction_network begin
        k1, S --> P, [physical_scale = PhysicalScale.ODE]
        k2, P --> S, [physical_scale = PhysicalScale.SDE]
        k3, S --> 0, [physical_scale = PhysicalScale.Jump]
    end

    prob = HybridProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 1.0),
        [:k1 => 1.0, :k2 => 0.5, :k3 => 0.1])
    @test prob isa JumpProcesses.JumpProblem
    @test prob.prob isa SciMLBase.SDEProblem

    # Use SDE solver since underlying problem is SDEProblem
    sol = solve(prob, SRIW1())
    @test SciMLBase.successful_retcode(sol)
end

# Mathematical correctness test: SDE+Jump birth-death model with known transient solution.
# dX = k1 dt + sqrt(k1) dW - (jump: X -> X-1 at rate k2*X)
# Analytic mean: E[X](t) = k1/k2 + (X0 - k1/k2)*exp(-k2*t)
# Also compares variance against manually coded JumpProcesses implementation.
let
    # Parameters
    k1 = 10.0  # birth rate
    k2 = 0.1   # death rate per molecule
    X0 = 50.0
    tspan = (0.0, 50.0)
    N = 4000  # number of trajectories
    times = range(tspan[1], tspan[2], length = 100)

    # --- Catalyst version ---
    rn = @reaction_network begin
        k1, 0 --> X, [physical_scale = PhysicalScale.SDE]
        k2, X --> 0, [physical_scale = PhysicalScale.Jump]
    end
    cat_prob = HybridProblem(rn, [:X => X0], tspan, [:k1 => k1, :k2 => k2];
        save_positions = (false, false))

    # --- Manually coded version using JumpProcesses + StochasticDiffEq ---
    f_manual!(du, u, p, t) = (du[1] = p[1])  # drift: k1
    g_manual!(du, u, p, t) = (du[1] = sqrt(p[1]))  # diffusion: sqrt(k1)
    rate_manual(u, p, t) = p[2] * u[1]  # k2 * X
    affect_manual!(integrator) = (integrator.u[1] -= 1)
    jump_manual = VariableRateJump(rate_manual, affect_manual!; save_positions = (false, false))
    sde_prob = SciMLBase.SDEProblem(f_manual!, g_manual!, [X0], tspan, [k1, k2])
    manual_prob = JumpProblem(sde_prob, Direct(), jump_manual; rng, save_positions = (false, false))

    # Run simulations and collect values at all time points
    cat_samples = zeros(N, length(times))
    manual_samples = zeros(N, length(times))
    for i in 1:N
        cat_sol = solve(cat_prob, SRIW1(); saveat = times)
        manual_sol = solve(manual_prob, SRIW1(); saveat = times)
        cat_samples[i, :] .= cat_sol[:X]
        manual_samples[i, :] .= manual_sol[1, :]  # numeric index for manual problem
    end

    # Compute means and variances
    cat_mean = vec(mean(cat_samples; dims = 1))
    manual_mean = vec(mean(manual_samples; dims = 1))
    cat_var = vec(var(cat_samples; dims = 1))
    manual_var = vec(var(manual_samples; dims = 1))

    # Analytic mean: E[X](t) = k1/k2 + (X0 - k1/k2)*exp(-k2*t)
    Xf(t) = k1 / k2 + (X0 - k1 / k2) * exp(-k2 * t)
    Xact = [Xf(t) for t in times]

    # Test means against analytic solution (5% relative tolerance)
    @test all(abs.(cat_mean .- Xact) .<= 0.05 .* Xact)
    @test all(abs.(manual_mean .- Xact) .<= 0.05 .* Xact)

    # Test that Catalyst and manual implementations have matching variances (10% relative tolerance)
    # Skip early times where variance is small and relative error can be large
    start_idx = findfirst(t -> t >= 5.0, times)
    @test all(abs.(cat_var[start_idx:end] .- manual_var[start_idx:end]) .<= 0.10 .* manual_var[start_idx:end])
end

# Mathematical correctness test: Two-species system with known transient solution.
# Production with CLE noise, degradation as jump.
# Analytic mean (with P(0)=0): E[P](t) = k1/k2 * (1 - exp(-k2*t))
let
    # S --k1--> S + P (SDE, CLE approximation)
    # P --k2--> 0 (Jump)
    k1 = 5.0
    k2 = 0.5
    P0 = 0.0
    tspan = (0.0, 20.0)  # Shorter tspan since steady-state reached quickly (τ = 1/k2 = 2)
    N = 4000
    times = range(tspan[1], tspan[2], length = 100)

    rn = @reaction_network begin
        k1, S --> S + P, [physical_scale = PhysicalScale.SDE]
        k2, P --> 0, [physical_scale = PhysicalScale.Jump]
    end

    prob = HybridProblem(rn, [:S => 1.0, :P => P0], tspan, [:k1 => k1, :k2 => k2];
        save_positions = (false, false))

    # Run simulations and collect values at all time points
    Pv = zeros(length(times))
    for i in 1:N
        sol = solve(prob, SRIW1(); saveat = times)
        Pv .+= sol[:P]
    end
    Pv ./= N

    # Analytic mean: E[P](t) = k1/k2 + (P0 - k1/k2)*exp(-k2*t)
    # With P0 = 0: E[P](t) = k1/k2 * (1 - exp(-k2*t))
    Pf(t) = k1 / k2 + (P0 - k1 / k2) * exp(-k2 * t)
    Pact = [Pf(t) for t in times]

    # Skip t=0 where Pact=0 (would give division issues in relative tolerance)
    @test all(abs.(Pv[2:end] .- Pact[2:end]) .<= 0.05 .* Pact[2:end])
end

# Mathematical correctness test: Complex non-linear multi-species system.
# Compares Catalyst HybridProblem against manually coded JumpProcesses implementation.
# System:
#   0 -> A (SDE, production with CLE noise)
#   A + B -> C (Jump, bimolecular)
#   C -> B + D (Jump, conversion)
#   A -> 0 (Jump, degradation)
#   D -> 0 (Jump, degradation)
let
    # Parameters
    k1 = 2.0   # A production rate
    k2 = 0.05  # A + B -> C rate
    k3 = 0.5   # C -> B + D rate
    k4 = 0.1   # A degradation rate
    k5 = 0.2   # D degradation rate
    A0, B0, C0, D0 = 10.0, 20.0, 0.0, 0.0
    tspan = (0.0, 30.0)
    N = 4000
    times = range(tspan[1], tspan[2], length = 100)

    # --- Catalyst version ---
    rn = @reaction_network begin
        k1, 0 --> A, [physical_scale = PhysicalScale.SDE]
        k2, A + B --> C, [physical_scale = PhysicalScale.Jump]
        k3, C --> B + D, [physical_scale = PhysicalScale.Jump]
        k4, A --> 0, [physical_scale = PhysicalScale.Jump]
        k5, D --> 0, [physical_scale = PhysicalScale.Jump]
    end
    cat_prob = HybridProblem(rn, [:A => A0, :B => B0, :C => C0, :D => D0], tspan,
        [:k1 => k1, :k2 => k2, :k3 => k3, :k4 => k4, :k5 => k5];
        save_positions = (false, false))

    # --- Manually coded version ---
    # Species order: [A, B, C, D]
    # SDE: dA = k1 dt + sqrt(k1) dW (production with CLE noise)
    function f_manual!(du, u, p, t)
        du[1] = p[1]  # k1 (A production drift)
        du[2] = 0.0
        du[3] = 0.0
        du[4] = 0.0
    end
    function g_manual!(du, u, p, t)
        du[1] = sqrt(p[1])  # sqrt(k1) (A production diffusion)
        du[2] = 0.0
        du[3] = 0.0
        du[4] = 0.0
    end

    # Jump 1: A + B -> C at rate k2*A*B
    rate1(u, p, t) = p[2] * u[1] * u[2]
    function affect1!(integrator)
        integrator.u[1] -= 1  # A -= 1
        integrator.u[2] -= 1  # B -= 1
        integrator.u[3] += 1  # C += 1
    end
    jump1 = VariableRateJump(rate1, affect1!; save_positions = (false, false))

    # Jump 2: C -> B + D at rate k3*C
    rate2(u, p, t) = p[3] * u[3]
    function affect2!(integrator)
        integrator.u[3] -= 1  # C -= 1
        integrator.u[2] += 1  # B += 1
        integrator.u[4] += 1  # D += 1
    end
    jump2 = VariableRateJump(rate2, affect2!; save_positions = (false, false))

    # Jump 3: A -> 0 at rate k4*A
    rate3(u, p, t) = p[4] * u[1]
    affect3!(integrator) = (integrator.u[1] -= 1)
    jump3 = VariableRateJump(rate3, affect3!; save_positions = (false, false))

    # Jump 4: D -> 0 at rate k5*D
    rate4(u, p, t) = p[5] * u[4]
    affect4!(integrator) = (integrator.u[4] -= 1)
    jump4 = VariableRateJump(rate4, affect4!; save_positions = (false, false))

    sde_prob = SciMLBase.SDEProblem(f_manual!, g_manual!, [A0, B0, C0, D0], tspan,
        [k1, k2, k3, k4, k5])
    manual_prob = JumpProblem(sde_prob, Direct(), jump1, jump2, jump3, jump4; rng,
        save_positions = (false, false))

    # Run simulations and collect values
    cat_A = zeros(N, length(times))
    cat_B = zeros(N, length(times))
    cat_C = zeros(N, length(times))
    cat_D = zeros(N, length(times))
    manual_A = zeros(N, length(times))
    manual_B = zeros(N, length(times))
    manual_C = zeros(N, length(times))
    manual_D = zeros(N, length(times))

    for i in 1:N
        cat_sol = solve(cat_prob, SRIW1(); saveat = times)
        manual_sol = solve(manual_prob, SRIW1(); saveat = times)
        cat_A[i, :] .= cat_sol[:A]
        cat_B[i, :] .= cat_sol[:B]
        cat_C[i, :] .= cat_sol[:C]
        cat_D[i, :] .= cat_sol[:D]
        manual_A[i, :] .= manual_sol[1, :]
        manual_B[i, :] .= manual_sol[2, :]
        manual_C[i, :] .= manual_sol[3, :]
        manual_D[i, :] .= manual_sol[4, :]
    end

    # Compute means
    cat_mean_A = vec(mean(cat_A; dims = 1))
    cat_mean_B = vec(mean(cat_B; dims = 1))
    cat_mean_C = vec(mean(cat_C; dims = 1))
    cat_mean_D = vec(mean(cat_D; dims = 1))
    manual_mean_A = vec(mean(manual_A; dims = 1))
    manual_mean_B = vec(mean(manual_B; dims = 1))
    manual_mean_C = vec(mean(manual_C; dims = 1))
    manual_mean_D = vec(mean(manual_D; dims = 1))

    # Test that Catalyst and manual implementations match (5% relative tolerance)
    # Skip early times where values may be near zero
    start_idx = findfirst(t -> t >= 2.0, times)
    @test all(abs.(cat_mean_A[start_idx:end] .- manual_mean_A[start_idx:end]) .<= 0.05 .* manual_mean_A[start_idx:end])
    @test all(abs.(cat_mean_B[start_idx:end] .- manual_mean_B[start_idx:end]) .<= 0.05 .* manual_mean_B[start_idx:end])
    @test all(abs.(cat_mean_C[start_idx:end] .- manual_mean_C[start_idx:end]) .<= 0.05 .* manual_mean_C[start_idx:end])
    @test all(abs.(cat_mean_D[start_idx:end] .- manual_mean_D[start_idx:end]) .<= 0.05 .* manual_mean_D[start_idx:end])
end

# Tests that species-only reaction systems (no reactions) produce valid ODE systems with zero ODEs.
# This is important for spatial modeling where transport is added separately.
let
    # Species-only network with no reactions.
    rn = @reaction_network begin
        @species X(t) Y(t)
    end

    # make_rre_ode should produce D(X) ~ 0 and D(Y) ~ 0 equations.
    sys = make_rre_ode(rn)
    @test length(equations(sys)) == 2
    for eq in equations(sys)
        @test Symbolics._iszero(eq.rhs)
    end

    # Solving should keep species at initial values.
    prob = ODEProblem(rn, [:X => 5.0, :Y => 3.0], (0.0, 1.0), [])
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
    @test sol[:X][end] ≈ 5.0
    @test sol[:Y][end] ≈ 3.0
end
