using Catalyst, ModelingToolkit, OrdinaryDiffEq, Test, LinearAlgebra, DiffEqJump

@parameters k α
@variables t, A(t), B(t), C(t), D(t)
rxs = [Reaction(t * k, [A], [B], [2 * α^2], [k + α * C])
       Reaction(1.0, [A, B], [C, D], [α, 2], [k, α])]
@named rs = ReactionSystem(rxs, t)
@test issetequal(states(rs), [A, B, C, D])
@test issetequal(parameters(rs), [k, α])
osys = convert(ODESystem, rs)

g = (k + α * C)
rs2 = @reaction_network rs begin
    t * k, 2 * α^2 * A --> $g * B
    1.0, α * A + 2 * B --> k * C + α * D
end k α
@test rs2 == rs

rxs2 = [(@reaction t * k, 2 * α^2 * A --> $g * B),
        (@reaction 1.0, α * A + 2 * B --> k * C + α * D)]
rs3 = ReactionSystem(rxs2, t; name = :rs)
@test rs3 == rs

u0map = [A => 3.0, B => 2.0, C => 3.0, D => 1.5]
pmap = (k => 2.5, α => 2)
tspan = (0.0, 5.0)
oprob = ODEProblem(osys, u0map, tspan, pmap)
# this is a hack because of https://github.com/SciML/ModelingToolkit.jl/issues/1475
oprob = remake(oprob; p = Tuple(pv[2] for pv in pmap))
du1 = zeros(size(oprob.u0))
oprob.f(du1, oprob.u0, oprob.p, 1.5)

function oderhs(du, u, p, t)
    k = p[1]
    α = p[2]
    A = u[1]
    B = u[2]
    C = u[3]
    D = u[4]
    n = 2 * α^2
    rl = t * k / factorial(n) * A^n
    rl2 = A^α * B^2 / (2 * factorial(α))
    du[1] = -n * rl - α * rl2
    du[2] = (k + α * C) * rl - 2 * rl2
    du[3] = k * rl2
    return du[4] = α * rl2
end
u0 = [uv[2] for uv in u0map]
p = Tuple(pv[2] for pv in pmap)
oprob2 = ODEProblem(oderhs, u0, tspan, p)
du2 = copy(du1)
oprob2.f(du2, oprob2.u0, oprob2.p, 1.5)
@test norm(du1 .- du2) < 100 * eps()

# test without rate law scalings
osys = convert(ODESystem, rs; combinatoric_ratelaws = false)
oprob = ODEProblem(osys, u0map, tspan, pmap)
function oderhs(du, u, p, t)
    k = p[1]
    α = p[2]
    A = u[1]
    B = u[2]
    C = u[3]
    D = u[4]
    n = 2 * α^2
    rl = t * k * A^n
    rl2 = A^α * B^2
    du[1] = -n * rl - α * rl2
    du[2] = (k + α * C) * rl - 2 * rl2
    du[3] = k * rl2
    return du[4] = α * rl2
end
oprob2 = ODEProblem(oderhs, [uv[2] for uv in u0map], tspan, oprob.p)
du1 .= 0;
du2 .= 0;
oprob.f(du1, oprob.u0, oprob.p, 1.5)
oprob2.f(du2, oprob2.u0, oprob2.p, 1.5)
@test norm(du1 .- du2) < 100 * eps()

# SDESystem test
ssys = convert(SDESystem, rs)
sf = SDEFunction{false}(ssys, states(ssys), parameters(ssys))
G = sf.g(u0, p, 1.0)
function sdenoise(u, p, t)
    k = p[1]
    α = p[2]
    A = u[1]
    B = u[2]
    C = u[3]
    D = u[4]
    n = 2 * α^2
    rl = sqrt(t * k / factorial(n) * A^n)
    rl2 = sqrt(A^α * B^2 / (2 * factorial(α)))
    return G = [-n*rl (-α*rl2);
                (k + α * C)*rl (-2*rl2);
                0.0 k*rl2;
                0.0 α*rl2]
end
G2 = sdenoise(u0, p, 1.0)
@test norm(G - G2) < 100 * eps()

# SDESystem test with no combinatoric rate laws
ssys = convert(SDESystem, rs; combinatoric_ratelaws = false)
sf = SDEFunction{false}(ssys, states(ssys), parameters(ssys))
G = sf.g(u0, p, 1.0)
function sdenoise(u, p, t)
    k = p[1]
    α = p[2]
    A = u[1]
    B = u[2]
    C = u[3]
    D = u[4]
    n = 2 * α^2
    rl = sqrt(t * k * A^n)
    rl2 = sqrt(A^α * B^2)
    return G = [-n*rl (-α*rl2);
                (k + α * C)*rl (-2*rl2);
                0.0 k*rl2;
                0.0 α*rl2]
end
G2 = sdenoise(u0, p, 1.0)
@test norm(G - G2) < 100 * eps()

# JumpSystem test
js = convert(JumpSystem, rs)
u0map = [A => 3, B => 2, C => 3, D => 5]
u0 = [uv[2] for uv in u0map]
pmap = (k => 5, α => 2)
p = Tuple(pv[2] for pv in pmap)
statetoid = Dict(state => i for (i, state) in enumerate(states(js)))
function r1(u, p, t)
    α = p[2]
    A = u[1]
    return t * k * binomial(A, 2 * α^2)
end
function affect1!(integrator)
    A = integrator.u[1]
    B = integrator.u[2]
    C = integrator.u[3]
    k = p[1]
    α = p[2]
    integrator.u[1] -= 2 * α^2
    integrator.u[2] += (k + α * C)
    return nothing
end
ttt = 1.5
j1 = VariableRateJump(r1, affect1!)
vrj = ModelingToolkit.assemble_vrj(js, equations(js)[2], statetoid)
@test isapprox(vrj.rate(u0, p, ttt), r1(u0, p, ttt))
fake_integrator1 = (u = copy(u0), p = p, t = ttt);
fake_integrator2 = deepcopy(fake_integrator1);
vrj.affect!(fake_integrator1);
affect1!(fake_integrator2);
@test fake_integrator1 == fake_integrator2

function r2(u, p, t)
    α = p[2]
    A = u[1]
    B = u[2]
    return binomial(A, α) * B * (B - 1) / 2
end
function affect2!(integrator)
    A = integrator.u[1]
    B = integrator.u[2]
    C = integrator.u[3]
    D = integrator.u[4]
    k = p[1]
    α = p[2]
    integrator.u[1] -= α
    integrator.u[2] -= 2
    integrator.u[3] += k
    integrator.u[4] += α
    return nothing
end
j2 = ConstantRateJump(r2, affect2!)
crj = ModelingToolkit.assemble_crj(js, equations(js)[1], statetoid)
@test isapprox(crj.rate(u0, p, ttt), r2(u0, p, ttt))
fake_integrator1 = (u = copy(u0), p = p, t = ttt);
fake_integrator2 = deepcopy(fake_integrator1);
crj.affect!(fake_integrator1);
affect2!(fake_integrator2);
@test fake_integrator1 == fake_integrator2

# a few simple solving tests via the SIR Model
@parameters α β γ k
@variables t, S(t), I(t), R(t)
rxs = [Reaction(α, [S, I], [I], [1, 1], [2]),
       Reaction(β, [I], [R], [1], [1])]
@named sir_ref = ReactionSystem(rxs, t)

rxs2 = [Reaction(α, [S, I], [I], [γ, 1], [k]),
        Reaction(β, [I], [R], [γ], [γ])]
@named sir = ReactionSystem(rxs2, t)

@test issetequal(states(sir_ref), states(sir))

# ODEs
p1 = (α => 0.1 / 1000, β => 0.01)
p2 = (α => 0.1 / 1000, β => 0.01, γ => 1, k => 2)
tspan = (0.0, 250.0)
u0 = [S => 999.0, I => 1.0, R => 0.0]
oprob = ODEProblem(sir_ref, u0, tspan, p1)
sol = solve(oprob, Tsit5())
# here we hack around https://github.com/SciML/ModelingToolkit.jl/issues/1475
p2dict = Dict(p2)
pvs = Tuple(p2dict[s] for s in parameters(sir))
@test all(isequal.(parameters(sir), parameters(convert(ODESystem, sir))))
oprob2 = ODEProblem(sir, u0, tspan, pvs)
sol2 = solve(oprob2, Tsit5())
@test norm(sol - sol2(sol.t)) < 1e-10

# jumps
Nsims = 10000
u0 = [S => 999, I => 1, R => 0]
jsys = convert(JumpSystem, sir_ref)
dprob = DiscreteProblem(jsys, u0, tspan, p1)
jprob = JumpProblem(jsys, dprob, Direct(); save_positions = (false, false))
function getmean(jprob, tf)
    m = zeros(3)
    for i in 1:Nsims
        sol = solve(jprob, SSAStepper())
        m .+= sol(tspan[2])
    end
    m /= Nsims
    return m
end
m1 = getmean(jprob, tspan[2])
jsys2 = convert(JumpSystem, sir)
p2dict = Dict(p2)
pvs = Tuple(p2dict[s] for s in parameters(sir))
@test all(isequal.(parameters(sir), parameters(jsys2)))
dprob2 = DiscreteProblem(jsys2, u0, tspan, pvs)
jprob2 = JumpProblem(jsys2, dprob2, Direct(); save_positions = (false, false))
m2 = getmean(jprob2, tspan[2])
@test maximum(abs.(m1[2:3] .- m2[2:3]) ./ m1[2:3]) < 0.05
