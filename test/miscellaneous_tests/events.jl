### Prepares Tests ###

# Fetch packages.
using Catalyst, OrdinaryDiffEq, Test

# Sets the default `t` and `D` to use.
t = default_t()
D = default_time_deriv()

### Basic Tests ### 

# Test discrete event is propagated to ODE solver correctly.
let
    # Creates model (essentially a jagged oscillation, where `V` is reset to 1.0 every `1.0` time units).
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
    oprob = ODEProblem(rs, [], (0.0, 20.0))
    sol = solve(oprob, Tsit5())
    @test sol(20.0, idxs = V) ≈ 2.5
end

let
    # Creates model.
    @parameters p d α = 1.0
    @species X(t) A(t) = 2
    @variables a(t) = 3
    rxs = [
        Reaction(p, nothing, [X]),
        Reaction(d, [X], nothing)
    ]
    continuous_events = [(a == t) => A ~ A + α]
    discrete_events = [2.0 => A ~ α + a]
    @named rs_ce = ReactionSystem(rxs, t; continuous_events)
    @named rs_de = ReactionSystem(rxs, t; discrete_events)
    continuous_events = [(a == t) => A ~ A + a]
    discrete_events = [2.0 => A ~ α + a]
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

    # Tests that species/variables/parameters can be accessed correctly one a MTK structure have been created.
    u0 = [X => 1]
    tspan = (0.0, 10.0)
    ps = [p => 10.0, d => 0.2]
    for XProb in [ODEProblem, SDEProblem, DiscreteProblem], rs in [rs_ce, rs_de, rs_ce_de]
        prob = XProb(rs, u0, (0.0, 10.0), ps)
        @test prob[A] == 2
        @test prob[a] == 3
        @test prob.ps[α] == 1.0
    end
end