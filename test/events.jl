using Test, Catalyst, ModelingToolkit, OrdinaryDiffEq

# test discrete event is propagated to ODE solver correctly
@variables t V(t)=1.0
D = Differential(t)
eqs = [D(V) ~ V]
discrete_events = [1.0 => [V ~ 1.0]]
rxs = [(@reaction $V, 0 --> A), (@reaction 1.0, A --> 0)]
@named rs = ReactionSystem([rxs; eqs], t; discrete_events)
setdefaults!(rs, [:A => 0.0])
oprob = ODEProblem(rs, [], (0.0, 20.0))
sol = solve(oprob, Tsit5())
@test isapprox(sol(10+10*eps(), idxs = V), 1.0)
