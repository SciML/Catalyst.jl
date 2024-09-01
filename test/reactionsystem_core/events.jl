### Prepares Tests ###

# Fetch packages.
using Catalyst, DiffEqCallbacks, JumpProcesses, OrdinaryDiffEq, StochasticDiffEq, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# Sets the default `t` and `D` to use.
t = default_t()
D = default_time_deriv()


### Basic Tests ###

# Test discrete event is propagated to ODE solver correctly.
let
    # Creates model (essentially a jagged oscillation, where `V` is reset to 1.0 every 1.0 time units).
    @variables V(t)=1.0
    eqs = [D(V) ~ V]
    discrete_events = [1.0 => [V ~ 1.0]]
    rxs = [
        @reaction $V, 0 --> A
        @reaction 1.0, A --> 0
    ]
    @named rs = ReactionSystem([rxs; eqs], t; discrete_events)
    @test length(ModelingToolkit.continuous_events(rs)) == 0
    @test length(ModelingToolkit.discrete_events(rs)) == 1

    # Tests in simulation.
    osys = complete(convert(ODESystem, complete(rs)))
    @test length(ModelingToolkit.continuous_events(osys)) == 0
    @test length(ModelingToolkit.discrete_events(osys)) == 1
    oprob = ODEProblem(osys, [osys.A => 0.0], (0.0, 20.0))
    sol = solve(oprob, Tsit5())
    @test sol(10 + 10*eps(), idxs = V) ≈ 1.0
end

# Test continuous event is propagated to the ODE solver.
let
    # Creates model (a production/degradation system, but both reactions stop at `t=2.5`).
    @parameters α=5.0 β=1.0
    @species V(t)=0.0
    rxs = [
        Reaction(α, nothing, [V]),
        Reaction(β, [V], nothing)
    ]
    continuous_events = [V ~ 2.5] => [α ~ 0, β ~ 0]
    @named rs = ReactionSystem(rxs, t; continuous_events)
    @test length(ModelingToolkit.continuous_events(rs)) == 1
    @test length(ModelingToolkit.discrete_events(rs)) == 0

    # Tests in simulation.
    osys = complete(convert(ODESystem, complete(rs)))
    @test length(ModelingToolkit.continuous_events(osys)) == 1
    @test length(ModelingToolkit.discrete_events(osys)) == 0
    oprob = ODEProblem(osys, [], (0.0, 20.0))
    sol = solve(oprob, Tsit5())
    @test sol(20.0, idxs = V) ≈ 2.5
end

# Tests that species/variables/parameters only encountered in events are added to `ReactionSystem`s properly.
# Tests for both discrete and continuous events. Tests that these quantities can be accessed in Problems.
# Tests that metadata for these quantities are saved properly
let
    # Creates model.
    @parameters p d α::Int64 = 1
    @species X(t) A(t) = 2 [description="A species"]
    @variables a(t) = 3
    rxs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing)
    ]
    continuous_events = [α ~ t] => [A ~ A + a]
    discrete_events = [2.0 => [A ~ α + a]]
    @named rs_ce = ReactionSystem(rxs, t; continuous_events)
    @named rs_de = ReactionSystem(rxs, t; discrete_events)
    continuous_events = [[α ~ t] => [A ~ A + α]]
    discrete_events = [2.0 => [A ~ a]]
    @named rs_ce_de = ReactionSystem(rxs, t; continuous_events, discrete_events)
    rs_ce = complete(rs_ce)
    rs_de = complete(rs_de)
    rs_ce_de = complete(rs_ce_de)

    # Tests model content.
    issetequal(species(rs_ce), [X, A])
    issetequal(species(rs_de), [X, A])
    issetequal(species(rs_ce_de), [X, A])
    issetequal(unknowns(rs_ce), [X, A, a])
    issetequal(unknowns(rs_de), [X, A, a])
    issetequal(unknowns(rs_ce_de), [X, A, a])
    issetequal(parameters(rs_ce), [p, d, α])
    issetequal(parameters(rs_de), [p, d, α])
    issetequal(parameters(rs_ce_de), [p, d, α])
    @test Symbolics.unwrap(rs_ce_de.α) isa Symbolics.BasicSymbolic{Int64}
    @test Symbolics.unwrap(rs_de.α) isa Symbolics.BasicSymbolic{Int64}
    @test Symbolics.unwrap(rs_ce_de.α) isa Symbolics.BasicSymbolic{Int64}
    @test ModelingToolkit.getdescription(rs_ce_de.A) == "A species"
    @test ModelingToolkit.getdescription(rs_de.A) == "A species"
    @test ModelingToolkit.getdescription(rs_ce_de.A) == "A species"

    # Tests that species/variables/parameters can be accessed correctly one a MTK problem have been created.
    u0 = [X => 1]
    tspan = (0.0, 10.0)
    ps = [p => 10.0, d => 0.2]
    for XProb in [ODEProblem, SDEProblem], rs in [rs_ce, rs_de, rs_ce_de]
        prob = XProb(rs, u0, (0.0, 10.0), ps)
        @test prob[A] == 2
        @test prob[a] == 3
        @test prob.ps[α] == 1
        @test prob.ps[α] isa Int64
    end

    # Handles `DiscreteProblem`s and `JumpProblem`s (these cannot contain continuous events or variables).
    discrete_events = [2.0 => [A ~ A + α]]
    @named rs_de_2 = ReactionSystem(rxs, t; discrete_events)
    rs_de_2 = complete(rs_de_2)
    dprob = DiscreteProblem(rs_de_2, u0, (0.0, 10.0), ps)
    jprob = JumpProblem(rs_de_2, dprob, Direct())
    for prob in [dprob, jprob]
        @test dprob[A] == 2
        @test dprob.ps[α] == 1
        @test dprob.ps[α] isa Int64
    end
end


### Event Input Checks ###

# Checks that singular events can be provided as vectors/not as vectors and this does not matter.
let
    @parameters p d
    @species X(t)
    rxs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing)
    ]
    ce = [X ~ 1.0] => [X ~ 0.5]
    de = [2.0] => [p ~ 1.0]
    rs1 = ReactionSystem(rxs, t; continuous_events = ce, discrete_events = de, name = :rs)
    rs2 = ReactionSystem(rxs, t; continuous_events = [ce], discrete_events = de, name = :rs)
    rs3 = ReactionSystem(rxs, t; continuous_events = ce, discrete_events = [de], name = :rs)
    rs4 = ReactionSystem(rxs, t; continuous_events = [ce], discrete_events = [de], name = :rs)
    @test rs1 == rs2 == rs3 == rs4
end

# Checks that various various erroneous forms yield errors.
# I.e. ensures affects/conditions requires vector forms in the right cases.
let
    # Prepares the model reaction.
    @parameters p d
    @species X(t)
    rxs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing)
    ]

    # Declares various misformatted events .
    @test_broken false # Some misformatted tests should throw error at this stage, but does not (https://github.com/SciML/ModelingToolkit.jl/issues/2612).
    continuous_events_bad = [
        X ~ 1.0 => [X ~ 0.5],       # Scalar condition.
        [X ~ 1.0] => X ~ 0.5,       # Scalar affect.
        (X ~ 1.0,) => [X ~ 0.5],    # Tuple condition.
        [X ~ 1.0] => (X ~ 0.5,),    # Tuple affect.
        [X - 1.0] => [X ~ 0.5],     # Non-equation condition (1).
        [X == 1.0] => [X ~ 0.5],    # Non-equation condition (2).
        # [X ~ 1.0] => [X ~ 0.5, p ~ 0.5], # No error on system creation, but permitted. Should probably throw an early error.
    ]
    discrete_events_bad = [
        [2.0] => p ~ 1.0,       # Scalar affect.
        [2.0] => (p ~ 1.0, ),   # Tuple affect.
        #[X > 2.0] => [p ~ 1.0], # Vector conditions. Should probably throw an error here already, currently does not.
        #(1.0, 2.0) => [p ~ 1.0] # Tuple condition. Should probably throw an error here already, currently does not.
    ]

    # Checks that errors are produced.
    for continuous_events in continuous_events_bad
        @test_throws Exception @named rs = ReactionSystem(rxs, t; continuous_events)
    end
    for discrete_events in discrete_events_bad
        @test_throws Exception @named rs = ReactionSystem(rxs, t; discrete_events)
    end
end


### DSL-based Tests ###

# Compares models with complicated events that are created programmatically/with the DSL.
# Checks that simulations are correct.
# Checks that various simulation inputs works.
# Checks continuous, discrete, preset time, and periodic events.
# Tests event affecting non-species components.
let
    # Creates model via DSL.
    rn_dsl = @reaction_network rn begin
        @parameters thres=7.0 dY_up
        @variables Z(t)
        @continuous_events begin
            [t ~ 2.5] => [p ~ p + 0.2]
            [X ~ thres, Y ~ X] => [X ~ X - 0.5, Z ~ Z + 0.1]
        end
        @discrete_events begin
            2.0 => [dX ~ dX + 0.01, dY ~ dY + dY_up]
            [1.0, 5.0] => [p ~ p - 0.1]
            (Z > Y) => [Z ~ Z - 0.1]
        end

        (p, dX), 0 <--> X
        (1.1*p, dY), 0 <--> Y
    end

    # Creates model programmatically.
    t = default_t()
    @variables Z(t)
    @species X(t) Y(t)
    @parameters p dX dY thres=7.0 dY_up
    rxs = [
        Reaction(p, nothing, [X], nothing, [1])
        Reaction(dX, [X], nothing, [1], nothing)
        Reaction(1.1*p, nothing, [Y], nothing, [1])
        Reaction(dY, [Y], nothing, [1], nothing)
    ]
    continuous_events = [
        [t ~ 2.5] => [p ~ p + 0.2]
        [X ~ thres, Y ~ X] => [X ~ X - 0.5, Z ~ Z + 0.1]
    ]
    discrete_events = [
        2.0 => [dX ~ dX + 0.01, dY ~ dY + dY_up]
        [1.0, 5.0] => [p ~ p - 0.1]
        (Z > Y) => [Z ~ Z - 0.1]
    ]
    rn_prog = ReactionSystem(rxs, t; continuous_events, discrete_events, name = :rn)
    rn_prog = complete(rn_prog)

    # Tests that approaches yield identical results.
    @test isequal(rn_dsl, rn_prog)

    u0 = [X => 6.0, Y => 4.5, Z => 5.5]
    tspan = (0.0, 20.0)
    ps = [p => 0.5, dX => 0.025, dY => 0.025, dY_up => 0.01]

    sol_dsl = solve(ODEProblem(rn_dsl, u0, tspan, ps), Tsit5())
    sol_prog = solve(ODEProblem(rn_prog, u0, tspan, ps), Tsit5())
    @test sol_dsl == sol_prog
end

# Checks that misformatted events yields errors in the DSL.
let
    # Quantity in event not declared elsewhere (continuous events).
    @test_throws Exception @eval @reaction_network begin
        @continuous_events X ~ 2.0 => [X ~ X + 1]
    end

    # Scalar condition (continuous events).
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @continuous_events X ~ 2.0 => [X ~ X + 1]
    end

    # Scalar affect (continuous events).
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @continuous_events [X ~ 2.0] => X ~ X + 1
    end

    # Tuple condition (continuous events).
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @continuous_events (X ~ 2.0,) => [X ~ X + 1]
    end

    # Tuple affect (continuous events).
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @continuous_events [X ~ 2.0] => (X ~ X + 1,)
    end

    # Non-equation condition (continuous events).
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @continuous_events [X - 2.0] => [X ~ X + 1]
    end

    # Quantity in event not declared elsewhere (discrete events).
    @test_throws Exception @eval @reaction_network begin
        @discrete_events X ~ 2.0 => [X ~ X + 1]
    end

    # Scalar affect (discrete events).
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @discrete_events 1.0 => X ~ X + 1
    end

    # Tuple affect (discrete events).
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @discrete_events 1.0 => (X ~ X + 1, )
    end

    # Equation condition (discrete events).
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @discrete_events X ~ 1.0 => [X ~ X + 1]
    end
end


### Additional Correctness Tests ###

# Tests that events are properly triggered for SDEs.
# Tests for continuous events, and all three types of discrete events.
let
    # Creates model with all types of events. The `e` parameters track whether events are triggered.
    rn = @reaction_network begin
        @parameters e1=0 e2=0 e3=0 e4=0
        @continuous_events begin
            [X ~ 1000.0] => [e1 ~ 1]
        end
        @discrete_events begin
            [1.0] => [e2 ~ 1]
            1.0 => [e3 ~ 1]
            (Y > 1000.0) & (e4==0) => [e4 ~ 1]
        end
        (p,d), 0 <--> X
        (p,d), 0 <--> Y
    end

    # Simulates the model for conditions where it *definitely* will cross `X = 1000.0`
    u0 = [:X => 999.9, :Y => 999.9]
    ps = [:p => 10.0, :d => 0.001]
    sprob = SDEProblem(rn, u0, (0.0, 2.0), ps)
    sol = solve(sprob, ImplicitEM(); seed)

    # Checks that all `e` parameters have been updated properly.
    @test sol.ps[:e1] == 1
    @test sol.ps[:e2] == 1
    @test sol.ps[:e3] == 1
    @test sol.ps[:e4] == 1
end

# Tests that events are properly triggered for Jump simulations.
# Tests for all three types of discrete events.
let
    # Creates model with all types of events. The `e` parameters track whether events are triggered.
    rn = @reaction_network begin
        @parameters e1=0 e2=0 e3=0
        @discrete_events begin
            [1.0] => [e1 ~ 1]
            1.0 => [e2 ~ 1]
            (X > 1000.0) & (e3==0) => [e3 ~ 1]
        end
        (p,d), 0 <--> X
    end

    # Simulates the model for conditions where it *definitely* will cross `X = 1000.0`
    u0 = [:X => 999]
    ps = [:p => 10.0, :d => 0.001]
    dprob = DiscreteProblem(rn, u0, (0.0, 2.0), ps)
    jprob = JumpProblem(rn, dprob, Direct(); rng)
    sol = solve(jprob, SSAStepper(); seed)

    # Checks that all `e` parameters have been updated properly.
    @test sol.ps[:e1] == 1
    @test sol.ps[:e2] == 1 
    @test sol.ps[:e3] == 1
end

# Compares simulations using MTK type events with those generated through callbacks.
# Jump simulations must be handles differently (since these only accepts discrete callbacks).
# Checks for all types of discrete callbacks, and for continuous callbacks.
# Turns of noise for SDE simulations (not sure seeding works when callbacks/events declared differently).
let
    # Creates models. Jump simulations requires one with discrete events only.
    rn = @reaction_network begin
        @default_noise_scaling 0.0
        @parameters add::Int64
        (p,d), 0 <--> X
        (p,d), 0 <--> Y
    end
    rn_events = @reaction_network begin
        @default_noise_scaling 0.0
        @parameters add::Int64
        @continuous_events begin
            [X ~ 90.0] => [X ~ X + 10.0]
        end
        @discrete_events begin
            [5.0, 10.0] => [X ~ X + add, Y ~ Y + add]
            20.0 => [X ~ X + add]
            (Y < X) => [Y ~ Y + add]
        end
        (p,d), 0 <--> X
        (p,d), 0 <--> Y
    end
    rn_dics_events = @reaction_network begin
        @parameters add::Int64
        @discrete_events begin
            [5.0, 10.0] => [X ~ X + add, Y ~ Y + add]
            20.0 => [X ~ X + add]
            (Y < X) => [Y ~ Y + add]
        end
        (p,d), 0 <--> X
        (p,d), 0 <--> Y
    end

    # Sets simulation inputs.
    u0 = [:X => 100, :Y => 150]
    tspan = (0.0, 90.0)
    ps = [:p => 100.0, :d => 1.0, :add => 100]

    # Create callbacks
    cb_cont = ContinuousCallback((u, t, int) -> (int[:X] - 90.0), int -> (int[:X] += 10.0))
    cb_disc_1 = PresetTimeCallback([5.0, 10.0], int -> (int[:X] += int.ps[:add]; int[:Y] += int.ps[:add];))
    cb_disc_2 = PresetTimeCallback(20.0:20.0:tspan[end], int -> (int[:X] += int.ps[:add]))
    cb_disc_3 = DiscreteCallback((u,t,i) -> i[:Y] < i[:X], int -> (int[:Y] += int.ps[:add]))
    callback = CallbackSet(cb_cont, cb_disc_1, cb_disc_2, cb_disc_3)

    # Checks for ODE simulations.
    oprob = ODEProblem(rn, u0, tspan, ps)
    oprob_events = ODEProblem(rn_events, u0, tspan, ps)
    osol = solve(oprob, Tsit5(); callback)
    osol_events = solve(oprob_events, Tsit5())
    @test osol == osol_events

    # Checks for SDE simulations (note, non-seed dependant test should be created instead).
    sprob = SDEProblem(rn, u0, tspan, ps)
    sprob_events = SDEProblem(rn_events, u0, tspan, ps)
    ssol = solve(sprob, ImplicitEM(); seed, callback)
    ssol_events = solve(sprob_events, ImplicitEM(); seed)
    @test ssol == ssol_events

    # Checks for Jump simulations. (note, non-seed dependant test should be created instead)
    # Note that periodic discrete events are currently broken for jump processes (and unlikely to be fixed soon due to have events are implemented).
    callback = CallbackSet(cb_disc_1, cb_disc_2, cb_disc_3)
    dprob = DiscreteProblem(rn, u0, tspan, ps)
    dprob_events = DiscreteProblem(rn_dics_events, u0, tspan, ps)
    jprob = JumpProblem(rn, dprob, Direct(); rng)
    jprob_events = JumpProblem(rn_dics_events, dprob_events, Direct(); rng)
    sol = solve(jprob, SSAStepper(); seed, callback)
    sol_events = solve(jprob_events, SSAStepper(); seed)
    @test_broken sol == sol_events  # seems to be not identical in the sample paths
end
