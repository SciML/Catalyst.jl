using Catalyst, LinearAlgebra, DiffEqJump, Test, OrdinaryDiffEq, StochasticDiffEq

const MT = ModelingToolkit

@parameters t k[1:20]
@variables A(t) B(t) C(t) D(t)
rxs = [Reaction(k[1], nothing, [A]),            # 0 -> A
       Reaction(k[2], [B], nothing),            # B -> 0
       Reaction(k[3],[A],[C]),                  # A -> C
       Reaction(k[4], [C], [A,B]),              # C -> A + B
       Reaction(k[5], [C], [A], [1], [2]),      # C -> A + A
       Reaction(k[6], [A,B], [C]),              # A + B -> C
       Reaction(k[7], [B], [A], [2], [1]),      # 2B -> A
       Reaction(k[8], [A,B], [A,C]),            # A + B -> A + C
       Reaction(k[9], [A,B], [C,D]),            # A + B -> C + D
       Reaction(k[10], [A], [C,D], [2], [1,1]), # 2A -> C + D
       Reaction(k[11], [A], [A,B], [2], [1,1]), # 2A -> A + B
       Reaction(k[12], [A,B,C], [C,D], [1,3,4], [2, 3]),          # A+3B+4C -> 2C + 3D
       Reaction(k[13], [A,B], nothing, [3,1], nothing),           # 3A+B -> 0
       Reaction(k[14], nothing, [A], nothing, [2]),               # 0 -> 2A
       Reaction(k[15]*A/(2+A), [A], nothing; only_use_rate=true), # A -> 0 with custom rate
       Reaction(k[16], [A], [B]; only_use_rate=true),             # A -> B with custom rate.
       Reaction(k[17]*A*exp(B), [C], [D], [2], [1]),              # 2C -> D with non constant rate.
       Reaction(k[18]*B, nothing, [B], nothing, [2]),             # 0 -> 2B with non constant rate.
       Reaction(k[19]*t, [A], [B]),                                # A -> B with non constant rate.
       Reaction(k[20]*t*A, [B,C], [D],[2,1],[2])                  # 2A +B -> 2C with non constant rate.
  ]
@named rs = ReactionSystem(rxs,t,[A,B,C,D],k)
odesys = convert(ODESystem,rs)
sdesys = convert(SDESystem,rs)

# test equation only constructor
@named rs2 = ReactionSystem(rxs,t)
@test Catalyst.isequal_ignore_names(rs,rs2)

# test show
io = IOBuffer()
show(io, rs)
str = String(take!(io))
@test count(isequal('\n'), str) < 30

# defaults test
def_p = [ki => float(i) for (i, ki) in enumerate(k)]
def_u0 = [A => 0.5, B => 1., C=> 1.5, D => 2.0]
defs = merge(Dict(def_p), Dict(def_u0))

@named rs = ReactionSystem(rxs,t,[A,B,C,D],k; defaults=defs)
odesys = convert(ODESystem,rs)
sdesys = convert(SDESystem,rs)
js = convert(JumpSystem,rs)
nlsys = convert(NonlinearSystem,rs)

@test ModelingToolkit.get_defaults(rs) ==
       ModelingToolkit.get_defaults(odesys) ==
       ModelingToolkit.get_defaults(sdesys) ==
       ModelingToolkit.get_defaults(js) ==
       ModelingToolkit.get_defaults(nlsys) ==
       defs

u0map = [A=>5.] # was 0.5
pmap = [k[1]=>5.] # was 1.
prob = ODEProblem(rs, u0map, (0,10.), pmap)
@test prob.p[1] == 5.
@test prob.u0[1] == 5.
u0 = [10., 11., 12., 13.]
ps = [float(x) for x in 100:119]
prob = ODEProblem(rs, u0, (0,10.), ps)
@test prob.p == ps
@test prob.u0 == u0

# hard coded ODE rhs
function oderhs(u,k,t)
       A = u[1]; B = u[2]; C = u[3]; D = u[4];
       du = zeros(eltype(u),4)
       du[1] = k[1] - k[3]*A + k[4]*C + 2*k[5]*C - k[6]*A*B + k[7]*B^2/2 - k[9]*A*B - k[10]*A^2 - k[11]*A^2/2 - k[12]*A*B^3*C^4/144 - 3*k[13]*A^3*B/6 + 2*k[14] - k[15]*A/(2+A) - k[16] - k[19]*t*A
       du[2] = -k[2]*B + k[4]*C - k[6]*A*B - k[7]*B^2 - k[8]*A*B - k[9]*A*B + k[11]*A^2/2 - 3*k[12]*A*B^3*C^4/144 - k[13]*A^3*B/6 + k[16] + 2*k[18]*B + k[19]*t*A - 2*k[20]*t*A*B^2*C
       du[3] = k[3]*A - k[4]*C - k[5]*C  + k[6]*A*B + k[8]*A*B + k[9]*A*B + k[10]*A^2/2 - 2*k[12]*A*B^3*C^4/144 - 2*k[17]*A*exp(B)*C^2/2 - k[20]*t*A*B^2*C
       du[4] = k[9]*A*B + k[10]*A^2/2 + 3*k[12]*A*B^3*C^4/144 + k[17]*A*exp(B)*C^2/2 + 2*k[20]*t*A*B^2*C
       du
end

# sde noise coefs
function sdenoise(u,k,t)
       A = u[1]; B = u[2]; C = u[3]; D = u[4];
       G = zeros(eltype(u),length(k),length(u))
       z = zero(eltype(u))

       G = [sqrt(k[1]) z  z z;
            z -sqrt(k[2]*B) z z;
            -sqrt(k[3]*A) z sqrt(k[3]*A) z;
            sqrt(k[4]*C) sqrt(k[4]*C) -sqrt(k[4]*C) z;
            2*sqrt(k[5]*C) z -sqrt(k[5]*C) z;
            -sqrt(k[6]*A*B) -sqrt(k[6]*A*B) sqrt(k[6]*A*B) z;
            sqrt(k[7]*B^2/2) -2*sqrt(k[7]*B^2/2) z z;
            z -sqrt(k[8]*A*B) sqrt(k[8]*A*B) z;
            -sqrt(k[9]*A*B) -sqrt(k[9]*A*B) sqrt(k[9]*A*B) sqrt(k[9]*A*B);
            -2*sqrt(k[10]*A^2/2) z sqrt(k[10]*A^2/2) sqrt(k[10]*A^2/2);
            -sqrt(k[11]*A^2/2) sqrt(k[11]*A^2/2) z z;
            -sqrt(k[12]*A*B^3*C^4/144) -3*sqrt(k[12]*A*B^3*C^4/144) -2*sqrt(k[12]*A*B^3*C^4/144) 3*sqrt(k[12]*A*B^3*C^4/144);
            -3*sqrt(k[13]*A^3*B/6) -sqrt(k[13]*A^3*B/6) z z;
            2*sqrt(k[14]) z z z;
            -sqrt(k[15]*A/(2+A)) z z z;
            -sqrt(k[16]) sqrt(k[16]) z z;
            z z -2*sqrt(k[17]*A*exp(B)*C^2/2) sqrt(k[17]*A*exp(B)*C^2/2);
            z 2*sqrt(k[18]*B) z z;
            -sqrt(k[19]*t*A) sqrt(k[19]*t*A) z z;
            z -2*sqrt(k[20]*t*A*B^2*C) -sqrt(k[20]*t*A*B^2*C) +2*sqrt(k[20]*t*A*B^2*C)]'
       return G
end

# test by evaluating drift and diffusion terms
p  = rand(length(k))
u  = rand(length(k))
t  = 0.
du = oderhs(u,p,t)
G  = sdenoise(u,p,t)
sdesys = convert(SDESystem,rs)
sf = SDEFunction{false}(sdesys, states(rs), parameters(rs))
du2 = sf.f(u,p,t)
@test norm(du-du2) < 100*eps()
G2 = sf.g(u,p,t)
@test norm(G-G2) < 100*eps()

# test conversion to NonlinearSystem
ns = convert(NonlinearSystem,rs)
fnl = eval(generate_function(ns)[2])
dunl = similar(du)
fnl(dunl,u,p)
@test norm(du-dunl) < 100*eps()

# tests the noise_scaling argument.
p  = rand(length(k)+1)
u  = rand(length(k))
t  = 0.
G  = p[21]*sdenoise(u,p,t)
@variables η
sdesys_noise_scaling = convert(SDESystem,rs;noise_scaling=η)
sf = SDEFunction{false}(sdesys_noise_scaling, states(rs), parameters(sdesys_noise_scaling))
G2 = sf.g(u,p,t)
@test norm(G-G2) < 100*eps()

# tests the noise_scaling vector argument.
p  = rand(length(k)+3)
u  = rand(length(k))
t  = 0.
G  = vcat(fill(p[21],8),fill(p[22],3),fill(p[23],9))' .* sdenoise(u,p,t)
@variables η[1:3]
sdesys_noise_scaling = convert(SDESystem,rs;noise_scaling=vcat(fill(η[1],8),fill(η[2],3),fill(η[3],9)))
sf = SDEFunction{false}(sdesys_noise_scaling, states(rs), parameters(sdesys_noise_scaling))
G2 = sf.g(u,p,t)
@test norm(G-G2) < 100*eps()

# tests using previous parameter for noise scaling
p  = rand(length(k))
u  = rand(length(k))
t  = 0.
G  = [p p p p]' .* sdenoise(u,p,t)
sdesys_noise_scaling = convert(SDESystem,rs;noise_scaling=k)
sf = SDEFunction{false}(sdesys_noise_scaling, states(rs), parameters(sdesys_noise_scaling))
G2 = sf.g(u,p,t)
@test norm(G-G2) < 100*eps()

# test with JumpSystem
js = convert(JumpSystem, rs)

midxs = 1:14
cidxs = 15:18
vidxs = 19:20
@test all(map(i -> typeof(equations(js)[i]) <: DiffEqJump.MassActionJump, midxs))
@test all(map(i -> typeof(equations(js)[i]) <: DiffEqJump.ConstantRateJump, cidxs))
@test all(map(i -> typeof(equations(js)[i]) <: DiffEqJump.VariableRateJump, vidxs))

pars = rand(length(k)); u0 = rand(1:10,4); ttt = rand();
jumps = Vector{Union{ConstantRateJump, MassActionJump, VariableRateJump}}(undef,length(rxs))

jumps[1] = MassActionJump(pars[1], Vector{Pair{Int,Int}}(), [1 => 1]);
jumps[2] = MassActionJump(pars[2], [2 => 1], [2 => -1]);
jumps[3] = MassActionJump(pars[3], [1 => 1], [1 => -1, 3 => 1]);
jumps[4] = MassActionJump(pars[4], [3 => 1], [1 => 1, 2 => 1, 3 => -1]);
jumps[5] = MassActionJump(pars[5], [3 => 1], [1 => 2, 3 => -1]);
jumps[6] = MassActionJump(pars[6], [1 => 1, 2 => 1], [1 => -1, 2 => -1, 3 => 1]);
jumps[7] = MassActionJump(pars[7], [2 => 2], [1 => 1, 2 => -2]);
jumps[8] = MassActionJump(pars[8], [1 => 1, 2 => 1], [2 => -1, 3 => 1]);
jumps[9] = MassActionJump(pars[9], [1 => 1, 2 => 1], [1 => -1, 2 => -1, 3 => 1, 4 => 1]);
jumps[10] = MassActionJump(pars[10], [1 => 2], [1 => -2, 3 => 1, 4 => 1]);
jumps[11] = MassActionJump(pars[11], [1 => 2], [1 => -1, 2 => 1]);
jumps[12] = MassActionJump(pars[12], [1 => 1, 2 => 3, 3 => 4], [1 => -1, 2 => -3, 3 => -2, 4 => 3]);
jumps[13] = MassActionJump(pars[13], [1 => 3, 2 => 1], [1 => -3, 2 => -1]);
jumps[14] = MassActionJump(pars[14], Vector{Pair{Int,Int}}(), [1 => 2]);

jumps[15] = ConstantRateJump((u,p,t) -> p[15]*u[1]/(2+u[1]), integrator -> (integrator.u[1] -= 1))
jumps[16] = ConstantRateJump((u,p,t) -> p[16], integrator -> (integrator.u[1] -= 1; integrator.u[2] += 1;))
jumps[17] = ConstantRateJump((u,p,t) -> p[17]*u[1]*exp(u[2])*binomial(u[3],2), integrator -> (integrator.u[3] -= 2; integrator.u[4] += 1))
jumps[18] = ConstantRateJump((u,p,t) -> p[18]*u[2], integrator -> (integrator.u[2] += 2))

jumps[19] = VariableRateJump((u,p,t) -> p[19]*u[1]*t, integrator -> (integrator.u[1] -= 1; integrator.u[2] += 1))
jumps[20] = VariableRateJump((u,p,t) -> p[20]*t*u[1]*binomial(u[2],2)*u[3], integrator -> (integrator.u[2] -= 2; integrator.u[3] -= 1; integrator.u[4] += 2))

statetoid = Dict(state => i for (i,state) in enumerate(states(js)))
jspmapper = ModelingToolkit.JumpSysMajParamMapper(js, pars)
symmaj = ModelingToolkit.assemble_maj(equations(js).x[1], statetoid, jspmapper)
maj    = MassActionJump(symmaj.param_mapper(pars), symmaj.reactant_stoch, symmaj.net_stoch, symmaj.param_mapper, scale_rates=false)
for i in midxs
  @test abs(jumps[i].scaled_rates - maj.scaled_rates[i]) < 100*eps()
  @test jumps[i].reactant_stoch == maj.reactant_stoch[i]
  @test jumps[i].net_stoch == maj.net_stoch[i]
end
for i in cidxs
  crj = ModelingToolkit.assemble_crj(js, equations(js)[i], statetoid)
  @test isapprox(crj.rate(u0,p,ttt), jumps[i].rate(u0,p,ttt))
  fake_integrator1 = (u=zeros(4),p=p,t=0); fake_integrator2 = deepcopy(fake_integrator1);
  crj.affect!(fake_integrator1);
  jumps[i].affect!(fake_integrator2);
  @test fake_integrator1 == fake_integrator2
end
for i in vidxs
  crj = ModelingToolkit.assemble_vrj(js, equations(js)[i], statetoid)
  @test isapprox(crj.rate(u0,p,ttt), jumps[i].rate(u0,p,ttt))
  fake_integrator1 = (u=zeros(4),p=p,t=0.); fake_integrator2 = deepcopy(fake_integrator1);
  crj.affect!(fake_integrator1); jumps[i].affect!(fake_integrator2);
  @test fake_integrator1 == fake_integrator2
end


# test for https://github.com/SciML/ModelingToolkit.jl/issues/436
@parameters t
@variables S(t) I(t)
rxs = [Reaction(1,[S],[I]), Reaction(1.1,[S],[I])]
@named rs = ReactionSystem(rxs, t, [S,I], [])
js = convert(JumpSystem, rs)
dprob = DiscreteProblem(js, [S => 1, I => 1], (0.0,10.0))
jprob = JumpProblem(js, dprob, Direct())
sol = solve(jprob, SSAStepper())

# test for https://github.com/SciML/ModelingToolkit.jl/issues/1042
jprob = JumpProblem(rs, dprob, Direct(), save_positions=(false,false))


@parameters k1 k2
@variables R(t)
rxs = [Reaction(k1*S, [S,I], [I], [2,3], [2]),
       Reaction(k2*R, [I], [R]) ]
@named rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])
@test isequal(oderatelaw(equations(rs)[1]), k1*S*S^2*I^3/(factorial(2)*factorial(3)))
@test_skip isequal(jumpratelaw(equations(eqs)[1]), k1*S*binomial(S,2)*binomial(I,3))
dep = Set()
ModelingToolkit.get_variables!(dep, rxs[2], Set(states(rs)))
dep2  = Set([R,I])
@test dep == dep2
dep = Set()
ModelingToolkit.modified_states!(dep, rxs[2], Set(states(rs)))
@test dep == Set([R,I])

isequal2(a,b) = isequal(simplify(a), simplify(b))

@test isequal2(jumpratelaw(rxs[1]), k1*S*S*(S-1)*I*(I-1)*(I-2)/12)
@test isequal2(jumpratelaw(rxs[1]; combinatoric_ratelaw=false), k1*S*S*(S-1)*I*(I-1)*(I-2))
@test isequal2(oderatelaw(rxs[1]), k1*S*S^2*I^3/12)
@test isequal2(oderatelaw(rxs[1]; combinatoric_ratelaw=false), k1*S*S^2*I^3)

@named rs2 = ReactionSystem(rxs, t, [S,I,R], [k1,k2]; combinatoric_ratelaws=false)

#test ODE scaling:
os = convert(ODESystem,rs)
@test isequal2(equations(os)[1].rhs, -2*k1*S*S^2*I^3/12)
os = convert(ODESystem,rs; combinatoric_ratelaws=false)
@test isequal2(equations(os)[1].rhs, -2*k1*S*S^2*I^3)
os2 = convert(ODESystem, rs2)
@test isequal2(equations(os2)[1].rhs, -2*k1*S*S^2*I^3)
os3 = convert(ODESystem, rs2; combinatoric_ratelaws=true)
@test isequal2(equations(os3)[1].rhs, -2*k1*S*S^2*I^3/12)

# test ConstantRateJump rate scaling
js = convert(JumpSystem,rs)
@test isequal2(equations(js)[1].rate, k1*S*S*(S-1)*I*(I-1)*(I-2)/12)
js = convert(JumpSystem,rs;combinatoric_ratelaws=false)
@test isequal2(equations(js)[1].rate, k1*S*S*(S-1)*I*(I-1)*(I-2))
js2 = convert(JumpSystem,rs2)
@test isequal2(equations(js2)[1].rate, k1*S*S*(S-1)*I*(I-1)*(I-2))
js3 = convert(JumpSystem,rs2; combinatoric_ratelaws=true)
@test isequal2(equations(js3)[1].rate, k1*S*S*(S-1)*I*(I-1)*(I-2)/12)

# test MassActionJump rate scaling
rxs = [Reaction(k1, [S,I], [I], [2,3], [2]),
       Reaction(k2, [I], [R]) ]
@named rs = ReactionSystem(rxs, t, [S,I,R], [k1,k2])
js = convert(JumpSystem, rs)
@test isequal2(equations(js)[1].scaled_rates, k1/12)
js = convert(JumpSystem,rs; combinatoric_ratelaws=false)
@test isequal2(equations(js)[1].scaled_rates, k1)

# test building directly from rxs
@parameters x,y
rxs = [Reaction(x*t*A*B+y, [A], nothing)]
@named rs1 = ReactionSystem(rxs, t, [A,B], [x,y])
@named rs2 = ReactionSystem(rxs, t)
@test Catalyst.isequal_ignore_names(rs1,rs2)

@variables t, L(t), H(t)
obs = [Equation(L, 2*x + y)]
@named rs3 = ReactionSystem(rxs, t; observed=obs)
L2 = L
@unpack L = rs3
@test isequal(L,L2)

# test non-integer stoichiometry goes through
@parameters k b
@variables t A(t) B(t) C(t) D(t)
rx1 = Reaction(k,[B,C],[B,D], [2.5,1],[3.5, 2.5])
rx2 = Reaction(2*k, [B], [D], [1], [2.5])
rx3 = Reaction(2*k, [B], [D], [2.5], [2])
@named mixedsys = ReactionSystem([rx1,rx2,rx3],t,[A,B,C,D],[k,b])
osys = convert(ODESystem, mixedsys; combinatoric_ratelaws=false)

# test for constant and boundary condition species
function f!(du,u,p,t)
    k1 = p[1]; k2 = p[2]; A = p[3]
    B = u[1]; D = u[2]; E = u[3]; C = u[4]
    du[1] = k1*A - k2*B
    du[2] = -k1*C*D + k2*E
    du[3] = k1*C*D - k2*E
    du[4] = -C
    nothing
end
function fs!(du,u,p,t)
    k1 = p[1]; k2 = p[2]; A = p[3]
    B = u[1]; D = u[2]; E = u[3]; C = u[4]
    du[1] = k1*A - k2*B
    du[2] = -k1*C*D + k2*E
    du[3] = k1*C*D - k2*E
    nothing
end
function gs!(dg,u,p,t)
    k1 = p[1]; k2 = p[2]; A = p[3]
    B = u[1]; D = u[2]; E = u[3]; C = u[4]
    dg .= 0.0
    dg[1,1] = sqrt(k1*A);    dg[1,2] = - sqrt(k2*B)
    dg[2,3] = -sqrt(k1*C*D); dg[2,4] = sqrt(k2*E)
    dg[3,3] = -dg[2,3];       dg[3,4] = -dg[2,4]
    nothing
end
let
    @parameters k1 k2
    @variables t A(t) [isconstant=true] B(t) C(t) [isbc=true] D(t) E(t)
    rxs = [(@reaction k1, $A --> B),
           (@reaction k2, B --> $A),
           (@reaction k1, $C + D --> E),
           (@reaction k2, E --> $C + D)]
    Dt = Differential(t)
    csys = ODESystem(Equation[Dt(C) ~ -C], t; name=:rs)
    @named rs = ReactionSystem(rxs, t; constraints=csys)
    osys = convert(ODESystem, rs)
    @test issetequal(MT.get_states(osys), [B, C, D, E])
    @test issetequal(MT.get_ps(osys), [k1, k2, A])

    # test nonlinear systems
    u0 = [1.0, 2.0, 3.0, 4.0]
    p = [2.5, 3.5, 2.0]
    f!(rand(4),u0,p,1.0)
    u0map = [B,D,E,C] .=> u0
    pmap = [k1,k2,A] .=> p
    tspan = (0.0,5.0)
    oprob1 = ODEProblem(osys, u0map, tspan, pmap)
    sts = [B,D,E,C]
    syms = [:B,:D,:E,:C]
    ofun = ODEFunction(f!; syms)
    oprob2 = ODEProblem(ofun, u0, tspan, p)
    saveat = tspan[2] / 50
    abstol = 1e-10
    reltol = 1e-10
    sol1 = solve(oprob1, Tsit5(); saveat, abstol, reltol)
    sol2 = solve(oprob2, Tsit5(); saveat, abstol, reltol)
    for i in eachindex(sts)
        @test isapprox(sol1[sts[i]], sol2[syms[i]])
    end

    # test sde systems
    @named rs = ReactionSystem(rxs, t)   # add constraint csys when supported!
    ssys = convert(SDESystem, rs)
    @test issetequal(MT.get_states(ssys), [B, C, D, E])
    @test issetequal(MT.get_ps(ssys), [k1, k2, A])
    du1 = zeros(4); du2 = zeros(4)
    sprob = SDEProblem(ssys, u0map, tspan, pmap; check_length=false)
    sprob.f(du1, u0, p, 1.0)
    fs!(du2, u0, p, 1.0)
    @test isapprox(du1, du2)
    dg1 = zeros(4,4); dg2 = zeros(4,4)
    sprob.g(dg1, u0, p, 1.0)
    gs!(dg2, u0, p, t)
    @test isapprox(dg1, dg2)

    # test jump systems
    rxs = [(@reaction k1, $A --> B),
           (@reaction k2, B --> $A),
           (@reaction k1, $C + D --> E),
           (@reaction k2, E --> $C + D),
           (@reaction k1*t, $A + $C--> B),
           (@reaction k1*B, 2*$A --> $C + B)]
    @named rs = ReactionSystem(rxs, t)
    jsys = convert(JumpSystem, rs)
    @test issetequal(states(jsys), [B,C,D,E])
    @test issetequal(parameters(jsys), [k1, k2, A])
    majrates = [k1*A, k2, k1*C, k2]
    majrs = [[],[B => 1],[D => 1],[E => 1]]
    majns = [[B => 1],[B => -1],[D => -1, E => 1],[D => 1, E => -1]]
    for (i,maj) in enumerate(equations(jsys).x[1])
        @test isequal(maj.scaled_rates, majrates[i])
        @test issetequal(maj.reactant_stoch, majrs[i])
        @test issetequal(maj.net_stoch, majns[i])
    end
    crj = equations(jsys).x[2][1]
    @test isequal(crj.rate, k1*B*A*(A-1)/2)
    @test issetequal(crj.affect!, [B ~ B + 1])
    vrj = equations(jsys).x[3][1]
    @test isequal(vrj.rate, k1*t*A*C)
    @test issetequal(vrj.affect!, [B ~ B + 1])
end
