using Catalyst, LinearAlgebra, OrdinaryDiffEq, Test, NonlinearSolve

# Repressilator model
@parameters t α₀ α K n δ β μ
@variables m(t) P(t) R(t)
rxs = [
       Reaction(α₀, nothing, [m]),
       Reaction(α / (1 + (R/K)^n), nothing, [m]),
       Reaction(δ, [m], nothing),
       Reaction(β, [m], [m,P]),
       Reaction(μ, [P], nothing)
      ]

specs = [m,P,R]
pars  = [α₀,α,K,n,δ,β,μ]
@named rs = ReactionSystem(rxs, t, specs, pars)

# using ODESystem components
@named sys₁ = convert(ODESystem, rs; include_zero_odes=false)
@named sys₂ = convert(ODESystem, rs; include_zero_odes=false)
@named sys₃ = convert(ODESystem, rs; include_zero_odes=false)
connections = [sys₁.R ~ sys₃.P,
               sys₂.R ~ sys₁.P,
               sys₃.R ~ sys₂.P]
@named connected = ODESystem(connections, t, [], [], systems=[sys₁,sys₂,sys₃])
oderepressilator = structural_simplify(connected)

pvals = [sys₁.α₀ => 5e-4,
         sys₁.α => .5,
         sys₁.K => 40.0,
         sys₁.n => 2,
         sys₁.δ => (log(2)/120),
         sys₁.β => (20*log(2)/120),
         sys₁.μ => (log(2)/600),
         sys₂.α₀ => 5e-4,
         sys₂.α => .5,
         sys₂.K => 40.0,
         sys₂.n => 2,
         sys₂.δ => (log(2)/120),
         sys₂.β => (20*log(2)/120),
         sys₂.μ => (log(2)/600),
         sys₃.α₀ => 5e-4,
         sys₃.α => .5,
         sys₃.K => 40.0,
         sys₃.n => 2,
         sys₃.δ => (log(2)/120),
         sys₃.β => (20*log(2)/120),
         sys₃.μ => (log(2)/600)]
u₀    = [sys₁.m => 0.0, sys₁.P => 20.0, sys₂.m => 0.0, sys₂.P => 0.0, sys₃.m => 0.0, sys₃.P => 0.0]
tspan = (0.0, 100000.0)
oprob = ODEProblem(oderepressilator, u₀, tspan, pvals)
sol = solve(oprob, Tsit5())

# hardcoded network
function repress!(f, y, p, t)
    α = p.α; α₀ = p.α₀; β = p.β; δ = p.δ; μ = p.μ; K = p.K; n = p.n
    f[1] = α / (1 + (y[6] / K)^n) - δ * y[1] + α₀
    f[2] = α / (1 + (y[4] / K)^n) - δ * y[2] + α₀
    f[3] = α / (1 + (y[5] / K)^n) - δ * y[3] + α₀
    f[4] = β * y[1] - μ * y[4]
    f[5] = β * y[2] - μ * y[5]
    f[6] = β * y[3] - μ * y[6]
    nothing
end
ps = (α₀=5e-4, α=.5, K=40.0, n=2, δ=(log(2)/120), β=(20*log(2)/120), μ=(log(2)/600))
u0 = [0.0,0.0,0.0,20.0,0.0,0.0]
oprob2 = ODEProblem(repress!, u0, tspan, ps)
sol2 = solve(oprob2, Tsit5())
tvs = 0:1:tspan[end]
@test all(isapprox.(sol(tvs,idxs=sys₁.P),sol2(tvs, idxs=4),atol=1e-4))

# using ReactionSystem components
@named sys₁ = ReactionSystem(rxs, t, specs, pars)
@named sys₂ = ReactionSystem(rxs, t, specs, pars)
@named sys₃ = ReactionSystem(rxs, t, specs, pars)
connections = [ParentScope(sys₁.R) ~ ParentScope(sys₃.P),
               ParentScope(sys₂.R) ~ ParentScope(sys₁.P),
               ParentScope(sys₃.R) ~ ParentScope(sys₂.P)]
@named csys = ODESystem(connections, t, [], [])
@named repressilator = ReactionSystem(t; systems=[csys,sys₁,sys₂,sys₃])
@named oderepressilator2 = convert(ODESystem, repressilator, include_zero_odes=false)
sys2 = structural_simplify(oderepressilator2)  # FAILS currently
oprob = ODEProblem(sys2, u₀, tspan, pvals)
sol = solve(oprob, Tsit5())
@test all(isapprox.(sol(tvs,idxs=sys₁.P),sol2(tvs, idxs=4),atol=1e-4))

# test conversion to nonlinear system 
@named nsys = NonlinearSystem(connections, [], [])
@named ssrepressilator = ReactionSystem(t; systems=[nsys,sys₁,sys₂,sys₃])
@named nlrepressilator = convert(NonlinearSystem, ssrepressilator, include_zero_odes=false)
sys2 = structural_simplify(nlrepressilator)
@test length(equations(sys2)) == 6
nlprob = NonlinearProblem(sys2, u₀, pvals)
sol = solve(nlprob, NewtonRaphson(), tol=1e-9)
@test sol[sys₁.P] ≈ sol[sys₂.P] ≈ sol[sys₃.P]
@test sol[sys₁.m] ≈ sol[sys₂.m] ≈ sol[sys₃.m]
@test sol[sys₁.R] ≈ sol[sys₂.R] ≈ sol[sys₃.R]

# TODO add conversion to SDE and JumpSystems once supported

# adding algebraic constraints
@parameters t, r₊, r₋, β
@variables A(t), B(t), C(t), D(t)
rxs1 = [Reaction(r₊, [A,B], [C])]
rxs2 = [Reaction(r₋, [C], [A,B])]
@named rs1 = ReactionSystem(rxs1, t, [A,B,C], [r₊])
@named rs2 = ReactionSystem(rxs2, t, [A,B,C], [r₋])
@named rs  = extend(rs1, rs2)
@test issetequal(states(rs), [A,B,C])
@test issetequal(parameters(rs), [r₊,r₋])
@test issetequal(equations(rs), union(rxs1,rxs2))
A2 = ModelingToolkit.ParentScope(A)
B2 = ModelingToolkit.ParentScope(B)
nseqs = [D ~ 2*A2 + β*B2]
@named ns = ODESystem(nseqs, t, [A2,B2,D], [β])
rs = compose(rs, [ns])
osys = convert(ODESystem, rs; include_zero_odes=false)
p  = [r₊ => 1.0, r₋ => 2.0, ns.β => 3.0]
u₀ = [A => 1.0, B => 2.0, C => 0.0]
oprob = ODEProblem(structural_simplify(osys), u₀, (0.0,10.0), p)
sol = solve(oprob, Tsit5())
@test isapprox(0, norm(sol[ns.D] .- 2*sol[A] - 3*sol[B]), atol=(100*eps())) 

# test API functions for composed model
@test issetequal(species(rs), [A,B,C])
@test issetequal(states(rs), [A,B,C,ns.D])
@test issetequal(reactionparams(rs), [r₊,r₋])
@test issetequal(parameters(rs), [r₊,r₋,ns.β])
@test issetequal(reactions(rs), union(rxs1,rxs2))
@test issetequal(filter(eq -> eq isa Reaction, equations(rs)), union(rxs1,rxs2))
@test issetequal(filter(eq -> eq isa Equation, equations(rs)), [ModelingToolkit.namespace_equation(nseqs[1],ns)])