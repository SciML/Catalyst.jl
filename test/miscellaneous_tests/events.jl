### Prepares Tests ###

# Fetch packages.
using Catalyst, OrdinaryDiffEq, Test

# Sets the default `t` and `D` to use.
t = default_t()
D = default_time_deriv()

### Basic Tests ### 

# Test discrete event is propagated to ODE solver correctly.
let
    # Creates model.
    @variables V(t)=1.0
    eqs = [D(V) ~ V]
    discrete_events = [1.0 => [V ~ 1.0]]
    rxs = [(@reaction $V, 0 --> A), (@reaction 1.0, A --> 0)]
    @named rs = ReactionSystem([rxs; eqs], t; discrete_events)
    rs = complete(rs)
    @test length(ModelingToolkit.discrete_events(rs)) == 1
    @test length(ModelingToolkit.continuous_events(rs)) == 0

    # Tests in simulation.
    setdefaults!(rs, [:A => 0.0])
    osys = complete(convert(ODESystem, rs))
    @test length(ModelingToolkit.discrete_events(osys)) == 1
    oprob = ODEProblem(osys, [], (0.0, 20.0))
    sol = solve(oprob, Tsit5())
    @test sol(10 + 10 * eps(), idxs = V) ≈ 1.0
end

# Test continuous event is propagated to the ODE solver.
let
    # Creates model.
    @parameters α=5.0 β=1.0
    @species V(t) = 0.0
    rxs = [Reaction(α, nothing, [V]), Reaction(β, [V], nothing)]
    continuous_events = [V ~ 2.5] => [α ~ 0, β ~ 0]
    @named rs = ReactionSystem(rxs, t; continuous_events)
    rs = complete(rs)
    @test length(ModelingToolkit.discrete_events(rs)) == 0
    @test length(ModelingToolkit.continuous_events(rs)) == 1

    # Tests in simulation.
    osys = convert(ODESystem, rs)
    @test length(ModelingToolkit.continuous_events(osys)) == 1
    oprob = ODEProblem(rs, [], (0.0, 20.0))
    sol = solve(oprob, Tsit5())
    @test sol(20.0, idxs = V) ≈ 2.5
end
