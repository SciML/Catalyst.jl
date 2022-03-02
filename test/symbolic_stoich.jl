using Catalyst, ModelingToolkit, OrdinaryDiffEq, Test, LinearAlgebra

@parameters t k α
@variables t, A(t), B(t), C(t), D(t)
rxs = [Reaction(t*k, [A], [B], [2*α^2], [k+α*C])
       Reaction(1.0, [A,B], [C,D], [α,2], [k,α])
]
@named rs = ReactionSystem(rxs, t)
@test issetequal(states(rs), [A,B,C,D])
@test issetequal(parameters(rs), [k,α])
osys = convert(ODESystem, rs)

u0map = [A => 3.0, B => 2.0, C => 3.0, D => 0.0]
pmap  = (k => 2.5, α => 2)
tspan = (0.0,5.0)
oprob = ODEProblem(osys, u0map, tspan, pmap)
oprob = remake(oprob, p=Tuple(pv[2] for pv in pmap))
sol   = solve(oprob, Tsit5())

function oderhs(du,u,p,t)
    k = p[1]; α = p[2]
    A = u[1]; B = u[2]; C = u[3]; D = u[4]
    n = 2*α^2
    @show n
    rl = t*k/factorial(n) * A^n
    rl2 = A^α*B^2 / (2 * factorial(α))
    du[1] = -n*rl - α*rl2
    du[2] = (k+α*C)*rl - 2*rl2
    du[3] = k*rl2
    du[4] = α*rl2
end
oprob2 = ODEProblem(oderhs, [uv[2] for uv in u0map], tspan, Tuple(pv[2] for pv in pmap))
sol2 = solve(oprob2, Tsit5())
@test norm(sol - sol2(sol.t)) < 100*eps()