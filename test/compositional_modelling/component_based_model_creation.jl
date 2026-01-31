#! format: off

### Prepares Tests ###

# Fetch packages.
using Catalyst, LinearAlgebra, OrdinaryDiffEqTsit5, NonlinearSolve, Test
using ModelingToolkitBase: nameof, getname
const MT = ModelingToolkitBase

# Sets the default `t` to use.
t = default_t()


### Run Tests ###

# Repressilator model.
let
    @parameters α₀ α K n δ β μ
    @species m(t) P(t) R(t)
    rxs = [
        Reaction(α₀, nothing, [m]),
        Reaction(α / (1 + (R / K)^n), nothing, [m]),
        Reaction(δ, [m], nothing),
        Reaction(β, [m], [m, P]),
        Reaction(μ, [P], nothing),
    ]

    specs = [m, P, R]
    pars = [α₀, α, K, n, δ, β, μ]
    @named rs = ReactionSystem(rxs, t, specs, pars)
    rs = complete(rs)

    # Using ODESystem components.
    @named sys₁ = make_rre_ode(rs; include_zero_odes = false)
    @named sys₂ = make_rre_ode(rs; include_zero_odes = false)
    @named sys₃ = make_rre_ode(rs; include_zero_odes = false)
    connections = [sys₁.R ~ sys₃.P,
        sys₂.R ~ sys₁.P,
        sys₃.R ~ sys₂.P]
    @named connected = System(connections, t, [], [], systems = [sys₁, sys₂, sys₃])
    oderepressilator = mtkcompile(connected)

    pvals = [sys₁.α₀ => 5e-4, sys₁.α => 0.5, sys₁.K => 40.0, sys₁.n => 2,
        sys₁.δ => (log(2) / 120), sys₁.β => (20 * log(2) / 120),
        sys₁.μ => (log(2) / 600), sys₂.α₀ => 5e-4, sys₂.α => 0.5, sys₂.K => 40.0,
        sys₂.n => 2, sys₂.δ => (log(2) / 120), sys₂.β => (20 * log(2) / 120),
        sys₂.μ => (log(2) / 600), sys₃.α₀ => 5e-4, sys₃.α => 0.5, sys₃.K => 40.0,
        sys₃.n => 2, sys₃.δ => (log(2) / 120), sys₃.β => (20 * log(2) / 120),
        sys₃.μ => (log(2) / 600)]
    u₀ = [sys₁.m => 0.0, sys₁.P => 20.0, sys₂.m => 0.0, sys₂.P => 0.0,
        sys₃.m => 0.0, sys₃.P => 0.0]
    tspan = (0.0, 100000.0)
    oprob = ODEProblem(oderepressilator, u₀, tspan, pvals)
    sol = solve(oprob, Tsit5())

    # Hardcoded network.
    function repress!(f, y, p, t)
        α = p.α
        α₀ = p.α₀
        β = p.β
        δ = p.δ
        μ = p.μ
        K = p.K
        n = p.n
        f[1] = α / (1 + (y[6] / K)^n) - δ * y[1] + α₀
        f[2] = α / (1 + (y[4] / K)^n) - δ * y[2] + α₀
        f[3] = α / (1 + (y[5] / K)^n) - δ * y[3] + α₀
        f[4] = β * y[1] - μ * y[4]
        f[5] = β * y[2] - μ * y[5]
        f[6] = β * y[3] - μ * y[6]
        nothing
    end
    ps = (α₀ = 5e-4, α = 0.5, K = 40.0, n = 2, δ = (log(2) / 120), β = (20 * log(2) / 120),
        μ = (log(2) / 600))
    u0 = [0.0, 0.0, 0.0, 20.0, 0.0, 0.0]
    oprob2 = ODEProblem(repress!, u0, tspan, ps)
    sol2 = solve(oprob2, Tsit5())
    tvs = 0:1:tspan[end]
    @test all(isapprox.(sol(tvs, idxs = sys₁.P), sol2(tvs, idxs = 4), atol = 1e-4))

    # Using ReactionSystem components.
    @named sys₁ = ReactionSystem(rxs, t, specs, pars)
    @named sys₂ = ReactionSystem(rxs, t, specs, pars)
    @named sys₃ = ReactionSystem(rxs, t, specs, pars)
    connections = [ParentScope(sys₁.R) ~ ParentScope(sys₃.P),
        ParentScope(sys₂.R) ~ ParentScope(sys₁.P),
        ParentScope(sys₃.R) ~ ParentScope(sys₂.P)]
    @named csys = ReactionSystem(connections, t; observed = Equation[])
    @named repressilator = ReactionSystem(t; systems = [csys, sys₁, sys₂, sys₃])
    repressilator = complete(repressilator)
    @named oderepressilator2 = make_rre_ode(repressilator, include_zero_odes = false)
    sys2 = mtkcompile(oderepressilator2)  # FAILS currently
    oprob = ODEProblem(sys2, [u₀; pvals], tspan)
    sol = solve(oprob, Tsit5())
    @test all(isapprox.(sol(tvs, idxs = sys₁.P), sol2(tvs, idxs = 4), atol = 1e-4))

    # Test conversion to nonlinear system.
    u₀_nl = [sys₁.m => 0.0, sys₁.P => 20.0, sys₁.R => 0.0, sys₂.m => 0.0, sys₂.P => 0.0,
        sys₂.R => 0.0, sys₃.m => 0.0, sys₃.P => 0.0, sys₃.R => 0.0]
    @named nsys = ReactionSystem(connections, t; observed = Equation[])
    @named ssrepressilator = ReactionSystem(t; systems = [nsys, sys₁, sys₂, sys₃])
    ssrepressilator = complete(ssrepressilator)
    @named nlrepressilator = make_rre_algeqs(ssrepressilator)
    sys2 = mtkcompile(nlrepressilator)
    @test_broken length(equations(sys2)) <= 6  # BUG in MTK: mtkcompile does not eliminate these like it does for ODEs
    nlprob = NonlinearProblem(sys2, [u₀_nl; pvals])
    sol = solve(nlprob; abstol = 1e-9)
    @test sol[sys₁.P] ≈ sol[sys₂.P] ≈ sol[sys₃.P]
    @test sol[sys₁.m] ≈ sol[sys₂.m] atol=1e-7
    @test sol[sys₁.m] ≈ sol[sys₃.m] atol=1e-7
    @test sol[sys₁.R] ≈ sol[sys₂.R] ≈ sol[sys₃.R]

    # Flattening.
    fsys = Catalyst.flatten(ssrepressilator)
    fsys = complete(fsys)
    @named nlrepressilator = make_rre_algeqs(fsys)
    sys2 = mtkcompile(nlrepressilator)
@test_broken length(equations(sys2)) <= 6  # BUG in MTK: mtkcompile does not eliminate these like it does for ODEs
    nlprob = NonlinearProblem(sys2, [u₀_nl; pvals])
    sol = solve(nlprob; abstol = 1e-9)
    @test sol[sys₁.P] ≈ sol[sys₂.P] ≈ sol[sys₃.P]
    @test sol[sys₁.m] ≈ sol[sys₂.m] atol=1e-7
    @test sol[sys₁.m] ≈ sol[sys₃.m] atol=1e-7
    @test sol[sys₁.R] ≈ sol[sys₂.R] ≈ sol[sys₃.R]

    # Test constraints.
    connections = [sys₁.R ~ sys₃.P,
        sys₂.R ~ sys₁.P,
        sys₃.R ~ sys₂.P]
    @named csys = ReactionSystem(connections, t; observed = Equation[])
    @named repressilator2 = ReactionSystem(connections, t; systems = [sys₁, sys₂, sys₃])
    repressilator2 = complete(repressilator2)
    @named nlrepressilator = make_rre_algeqs(repressilator2)
    sys2 = mtkcompile(nlrepressilator)
    @test_broken length(equations(sys2)) <= 6  # BUG in MTK: mtkcompile does not eliminate these like it does for ODEs
    nlprob = NonlinearProblem(sys2, [u₀_nl; pvals])
    sol = solve(nlprob; abstol = 1e-9)
    @test sol[sys₁.P] ≈ sol[sys₂.P] ≈ sol[sys₃.P]
    @test sol[sys₁.m] ≈ sol[sys₂.m] atol=1e-7
    @test sol[sys₁.m] ≈ sol[sys₃.m] atol=1e-7
    @test sol[sys₁.R] ≈ sol[sys₂.R] ≈ sol[sys₃.R]

    # Test constraint system variables are accessible through Base.getproperty
    # even if they do not appear in the original ReactionSystem.
    network = @network_component
    @parameters a
    @variables x(t)
    @named constraints = ReactionSystem([x ~ a], t; observed = Equation[])
    extended = extend(constraints, network)
    @test isequal(extended.a, MT.namespace_expr(a, extended))
    @test isequal(extended.x, MT.namespace_expr(x, extended))
    # and after conversion to an AbstractSystem
    extended = complete(extended)
    system = make_rre_algeqs(extended)
    @test isequal(system.a, MT.namespace_expr(a, system))
    @test isequal(system.x, MT.namespace_expr(x, system; ivs = independent_variables(extended)))
    @test length(equations(system)) == 1
    @test equations(system) == equations(constraints)

    # Test that the namespacing still works if the extended system takes the name
    # of the ReactionSystem.
    @named constraints2 = ReactionSystem([x ~ a], t; observed = Equation[])
    extended = extend(constraints2, network; name = nameof(network))
    @test isequal(extended.a, MT.namespace_expr(a, extended))
    @test isequal(extended.x, MT.namespace_expr(x, extended))
    # and after conversion to an AbstractSystem.
    extended = complete(extended)
    system = make_rre_algeqs(extended)
    @test isequal(system.a, MT.namespace_expr(a, system))
    @test isequal(system.x, MT.namespace_expr(x, system; ivs = independent_variables(extended)))
    @test length(equations(system)) == 1
    @test Set(equations(system)) == Set(equations(constraints2))

    # Test that extending a system with constraints correctly handles default values.
    network = @network_component
    subnetwork = @network_component
    @parameters a=1 b=2
    @variables x(t)=a y(t)=b
    @named constraints = ReactionSystem([x ~ a], t; observed = Equation[])
    @named subsystemconstraints = ReactionSystem([y ~ b], t; observed = Equation[])
    extended = extend(constraints, network)
    subextended = extend(subsystemconstraints, subnetwork)
    extended = compose(extended, subextended)
    defs = MT.initial_conditions(extended)
    @test Symbolics.value(get(defs, a, nothing)) == 1
    @test isequal(Symbolics.getdefaultval(x), a)
    @test Symbolics.value(get(defs, subextended.b, nothing)) == 2
    # Note: default values are stored with original (non-namespaced) symbols
    @test isequal(Symbolics.getdefaultval(subextended.y), b)

    extended = extend(constraints, network; name = nameof(network))
    subextended = extend(subsystemconstraints, subnetwork, name = nameof(subnetwork))
    extended = compose(extended, subextended)
    defs = MT.initial_conditions(extended)
    @test Symbolics.value(get(defs, a, nothing)) == 1
    @test isequal(Symbolics.getdefaultval(x), a)
    @test Symbolics.value(get(defs, subextended.b, nothing)) == 2
    @test isequal(Symbolics.getdefaultval(subextended.y), b)

    # Test that the observables of constraint systems are accessible after
    # extending a ReactionSystem.
    network = @network_component
    subnetwork = @network_component
    @parameters a b
    @variables x(t) y(t)
    @named constraints = ReactionSystem([x ~ a], t; observed = Equation[])
    @named subconstraints = ReactionSystem([y ~ b], t; observed = Equation[])

    extended = extend(constraints, network; name = nameof(network))
    subextended = extend(subconstraints, subnetwork, name = nameof(subnetwork))
    extended = compose(extended, subextended)
    @test isequal(extended.a, MT.namespace_expr(a, extended))
    @test isequal(extended.x, MT.namespace_expr(x, extended))
    extended = complete(extended)
    odesystem = complete(make_rre_ode(extended))
    nlsystem = complete(make_rre_algeqs(extended))

    obs = Set([MT.observed(constraints);
            [MT.namespace_equation(o, subextended)
                for o in MT.observed(subconstraints)]])
    @test Set(MT.observed(extended)) == obs
    @test Set(MT.observed(odesystem)) == obs
    @test Set(MT.observed(nlsystem)) == obs

    extended = extend(constraints, network)
    subextended = extend(subconstraints, subnetwork)
    extended = compose(extended, subextended)
    @test isequal(extended.a, MT.namespace_expr(a, extended))
    @test isequal(extended.x, MT.namespace_expr(x, extended))
    extended = complete(extended)
    odesystem = complete(make_rre_ode(extended))
    nlsystem = complete(make_rre_algeqs(extended))

    obs = Set([MT.observed(constraints);
            [MT.namespace_equation(o, subextended)
                for o in MT.observed(subconstraints)]])
    @test Set(MT.observed(extended)) == obs
    @test Set(MT.observed(odesystem)) == obs
    @test Set(MT.observed(nlsystem)) == obs

    # Test can make ODESystem.
    @named oderepressilator = make_rre_ode(repressilator2, include_zero_odes = false)
    sys2 = mtkcompile(oderepressilator)  # FAILS currently
    oprob = ODEProblem(sys2, [u₀; pvals], tspan)
    sol = solve(oprob, Tsit5())
    @test all(isapprox.(sol(tvs, idxs = sys₁.P), sol2(tvs, idxs = 4), atol = 1e-4))

    # Test extending with NonlinearSystem.
    @named repressilator2 = ReactionSystem(t; systems = [sys₁, sys₂, sys₃])
    repressilator2 = Catalyst.flatten(repressilator2)
    repressilator2 = extend(csys, repressilator2)
    repressilator2 = complete(repressilator2)
    @named nlrepressilator = make_rre_algeqs(repressilator2)
    sys2 = mtkcompile(nlrepressilator)
    @test_broken length(equations(sys2)) <= 6  # BUG in MTK: mtkcompile does not eliminate these like it does for ODEs
    nlprob = NonlinearProblem(sys2, [u₀_nl; pvals])
    sol = solve(nlprob; abstol = 1e-9)
    @test sol[sys₁.P] ≈ sol[sys₂.P] ≈ sol[sys₃.P]
    @test sol[sys₁.m] ≈ sol[sys₂.m] atol=1e-7
    @test sol[sys₁.m] ≈ sol[sys₃.m] atol=1e-7
    @test sol[sys₁.R] ≈ sol[sys₂.R] ≈ sol[sys₃.R]
end

# TODO add conversion to SDE and JumpSystems once supported.

# Adding algebraic constraints.
let
    @parameters r₊, r₋, β
    @species A(t), B(t), C(t)
    @variables D(t)
    rxs1 = [Reaction(r₊, [A, B], [C])]
    rxs2 = [Reaction(r₋, [C], [A, B])]
    @named rs1 = ReactionSystem(rxs1, t, [A, B, C], [r₊])
    @named rs2 = ReactionSystem(rxs2, t, [A, B, C], [r₋])
    @named rs = extend(rs1, rs2)
    @test issetequal(unknowns(rs), [A, B, C])
    @test issetequal(parameters(rs), [r₊, r₋])
    @test issetequal(equations(rs), union(rxs1, rxs2))
    A2 = MT.ParentScope(A)
    B2 = MT.ParentScope(B)
    nseqs = [D ~ 2 * A2 + β * B2]
    @named ns = ReactionSystem(nseqs, t; observed = Equation[])
    rs = compose(rs, [ns])
    rs = complete(rs)
    osys = make_rre_ode(rs; include_zero_odes = false)
    p = [r₊ => 1.0, r₋ => 2.0, ns.β => 3.0]
    u₀ = [A => 1.0, B => 2.0, C => 0.0]
    oprob = ODEProblem(mtkcompile(osys), u₀, (0.0, 10.0), p)
    sol = solve(oprob, Tsit5())
    @test isapprox(0, norm(sol[ns.D] .- 2 * sol[A] - 3 * sol[B]), atol = (100 * eps()))
    psyms = [:r₊ => 1.0, :r₋ => 2.0, :ns₊β => 3.0]
    u₀syms = [:A => 1.0, :B => 2.0, :C => 0.0]
    p = symmap_to_varmap(osys, psyms)
    u₀ = symmap_to_varmap(osys, u₀syms)
    oprob = ODEProblem(mtkcompile(osys), u₀, (0.0, 10.0), p)
    sol = solve(oprob, Tsit5())
    @test isapprox(0, norm(sol[ns.D] .- 2 * sol[A] - 3 * sol[B]), atol = (100 * eps()))

    # Test API functions for composed model.
    @test issetequal(species(rs), [A, B, C])
    @test issetequal(unknowns(rs), [A, B, C, ns.D])
    @test issetequal(parameters(rs), [r₊, r₋, ns.β])
    @test issetequal(reactions(rs), union(rxs1, rxs2))
    @test issetequal(filter(eq -> eq isa Reaction, equations(rs)), union(rxs1, rxs2))
    @test issetequal(filter(eq -> eq isa Equation, equations(rs)),
                    [MT.namespace_equation(nseqs[1], ns)])

    # Check several levels of nesting namespace and filter ok for the API functions.
    @parameters p1, p2a, p2b, p3a, p3b
    @species A1(t), A2a(t), A3a(t)
    @variables A2b(t), A3b(t)  # Variables for algebraic equations (not species)
    rxs1 = [Reaction(p1, [A1], nothing)]
    rxs2 = [Reaction(p2a, [A2a], nothing), Reaction(p2b, [ParentScope(A1)], nothing)]
    eqs2 = [ParentScope(A1) ~ ParentScope(p1) * A2b]
    rxs3 = [
        Reaction(p3a, [A3a], nothing),
        Reaction(ParentScope(p2a), nothing, [ParentScope(A2a)]),
    ]
    eqs3 = [ParentScope(A2a) ~ p3b * A3b]

    @named rs3 = ReactionSystem(rxs3, t, [A3a, ParentScope(A2a)], [p3a, ParentScope(p2a)];
                                combinatoric_ratelaws = false)
    @named ns3 = ReactionSystem(eqs3, t; observed = Equation[])
    @named rs2 = ReactionSystem(rxs2, t, [A2a, ParentScope(A1)], [p2a, p2b],
                                systems = [rs3, ns3]; combinatoric_ratelaws = true)
    @named ns2 = ReactionSystem(eqs2, t; observed = Equation[])
    @named rs1 = ReactionSystem(rxs1, t, [A1], [p1], systems = [rs2, ns2];
                                combinatoric_ratelaws = false)

    # Namespaced reactions.
    nrxs1 = [Reaction(p1, [A1], nothing)]
    nrxs2 = [Reaction(rs2.p2a, [rs2.A2a], nothing), Reaction(rs2.p2b, [A1], nothing)]
    neqs2 = [A1 ~ p1 * ns2.A2b]
    nrxs3 = [
        Reaction(rs2.rs3.p3a, [rs2.rs3.A3a], nothing),
        Reaction(rs2.p2a, nothing, [rs2.A2a]),
    ]
    neqs3 = [rs2.A2a ~ rs2.ns3.p3b * rs2.ns3.A3b]
    rxs = vcat(nrxs1, nrxs2, nrxs3)
    eqs = vcat(nrxs1, nrxs2, neqs2, nrxs3, neqs3)

    @test issetequal(unknowns(rs1), [A1, rs2.A2a, ns2.A2b, rs2.rs3.A3a, rs2.ns3.A3b]) 
    @test issetequal(species(rs1), [A1, rs2.A2a, rs2.rs3.A3a])
    @test issetequal(parameters(rs1), [p1, rs2.p2a, rs2.p2b, rs2.rs3.p3a, rs2.ns3.p3b])
    @test issetequal(rxs, reactions(rs1))
    @test issetequal(eqs, equations(rs1))
    @test Catalyst.combinatoric_ratelaws(rs1)
    @test Catalyst.combinatoric_ratelaws(Catalyst.flatten(rs1))
end

# Test ODE and algebraic constraints with extend and compose.
let
    rn = @network_component rn begin
        @parameters k1 k2
        (k1, k2), A <--> B
    end
    @parameters a, b
    @unpack A = rn
    @variables C(t)
    D = default_time_deriv()
    eqs = [D(C) ~ -b * C + a * A]
    @named osys = ReactionSystem(eqs, t; observed = Equation[])
    rn2 = extend(osys, rn)
    rnodes = make_rre_ode(complete(rn2))

    # Ensure right number of equations are generated.
    @variables G(t)
    eqs = [D(G) ~ -G]
    @named osys2 = ReactionSystem(eqs, t; observed = Equation[])
    rn3 = compose(rn2, osys2)
    @test length(equations(rn3)) == 4

    # Check conversions work with algebraic constraints.
    # Note: ReactionSystems must be created WITH an iv to extend with other ReactionSystems.
    eqs = [0 ~ -a * A + C, 0 ~ -b * C + a * A]
    @named nlsys = ReactionSystem(eqs, t; observed = Equation[])
    rn2 = complete(extend(nlsys, rn))
    rnodes = make_rre_ode(rn2)
    rnnlsys = make_rre_algeqs(rn2)
end

# https://github.com/SciML/ModelingToolkit.jl/issues/1274
let
    @parameters p1 p2
    @species A(t)
    rxs1 = [Reaction(p1, [A], nothing)]
    rxs2 = [Reaction(p2, [ParentScope(A)], nothing)]
    @named rs1 = ReactionSystem(rxs1, t)
    @named rs2 = ReactionSystem(rxs2, t)
    rsc = compose(rs1, [rs2])
    rsc = complete(rsc)
    orsc = make_rre_ode(rsc)
    @test length(equations(orsc)) == 1
end


# https://github.com/SciML/Catalyst.jl/issues/545
let
    rn_AB = @reaction_network AB begin
        @parameters k1 n
        k1, A --> n*B
    end

    rn_BC = @reaction_network BC begin
        @parameters k2
        k2, B --> C
    end

    @named rs = ReactionSystem(t; systems = [rn_AB, rn_BC])
    sts = unknowns(rs)
    @test issetequal(sts, (@species AB₊A(t) AB₊B(t) BC₊B(t) BC₊C(t)))
    ps = parameters(rs)
    @test issetequal(ps, (@parameters AB₊k1 AB₊n BC₊k2))
    rxs = reactions(rs)
    @parameters AB₊n
    rxs2 = Reaction[(@reaction AB₊k1, AB₊A --> $(AB₊n)*AB₊B), (@reaction BC₊k2, BC₊B --> BC₊C)]
    @test (length(rxs) == length(rxs2)) && issubset(rxs, rxs2)
end

# Test ordering of unknowns and equations.
let
    @parameters k1 k2 k3
    @variables V1(t) V2(t) V3(t)
    @species A1(t) A2(t) A3(t) B1(t) B2(t) B3(t)
    D = default_time_deriv()
    rx1 = Reaction(k1*V1, [A1], [B1])
    eq1 = D(V1) ~ -V1
    @named rs1 = ReactionSystem([rx1, eq1], t)
    rx2 = Reaction(k2*V2, [A2], [B2])
    eq2 = D(V2) ~ -V2
    @named rs2 = ReactionSystem([rx2, eq2], t)
    rx3 = Reaction(k3*V3, [A3], [B3])
    eq3 = D(V3) ~ -V3
    @named rs3 = ReactionSystem([rx3, eq3], t)
    @named rs23 = compose(rs2, [rs3])
    @test length(unknowns(rs23)) == 6
    @test all(p -> isequal(p[1], p[2]), zip(unknowns(rs23)[1:4], species(rs23)))
    @test length(equations(rs23)) == 4
    @test all(p -> isequal(p[1], p[2]), zip(equations(rs23)[1:2], reactions(rs23)))
    @named rs123 = compose(rs1, [rs23])
    @test length(unknowns(rs123)) == 9
    @test all(p -> isequal(p[1], p[2]), zip(unknowns(rs123)[1:6], species(rs123)))
    @test length(equations(rs123)) == 6
    @test length(reactions(rs123)) == 3
    @test all(p -> isequal(p[1], p[2]), zip(equations(rs123)[1:3], reactions(rs123)))

    @test numspecies(rs123) == 6
    @test issetequal(nonspecies(rs123), [V1, rs23.V2, rs23.rs3.V3])
end

# Tests that conversion with defaults works for a composed model.
let
    rn1 = @network_component rn1 begin
        @parameters p=1.0 r=2.0
        @species X(t) = 3.0 Y(t) = 4.0
        (p1, d), 0 <--> X
        (p2, r), 0 <--> Z
    end
    rn2 = @network_component rn1 begin
        @parameters p=10. q=20.0
        @species X(t) = 30.0 Z(t) = 40.0
        (p1, d), 0 <--> X
        (p2, q), 0 <--> Z
    end
    composed_reaction_system = compose(rn1, [rn2])
    composed_reaction_system = complete(composed_reaction_system)
    osys = make_rre_ode(composed_reaction_system)
    parameters(osys)[1].metadata

    defs = MT.initial_conditions(osys)
    @unpack p, r, X, Y = rn1
    @test Symbolics.value(defs[p]) == 1.0
    @test Symbolics.value(defs[r]) == 2.0
    @test Symbolics.value(defs[X]) == 3.0
    @test Symbolics.value(defs[Y]) == 4.0
    @test Symbolics.value(defs[rn2.p]) == 10.0
    @test Symbolics.value(defs[rn2.q]) == 20.0
    @test Symbolics.value(defs[rn2.X]) == 30.0
    @test Symbolics.value(defs[rn2.Z]) == 40.0
end

# test scoping in compose
# code adapted from ModelingToolkit.jl tests
# Note: DelayParentScope was removed in MTK v10/v11
let
    t = default_t()
    D = default_time_deriv()
    @species x1(t) x2(t)
    @variables x3(t) x5(t)
    x2 = ParentScope(x2)
    x3 = ParentScope(ParentScope(x3))
    x5 = GlobalScope(x5)
    @parameters p1 p2 p3 p5
    p2 = ParentScope(p2)
    p3 = ParentScope(ParentScope(p3))
    p5 = GlobalScope(p5)
    rxs = [Reaction(p1, nothing, [x1]), Reaction(p2, [x2], nothing),
           D(x3) ~ p3, D(x5) ~ p5]
    @named sys1 = ReactionSystem(rxs, t)
    @test isequal(x1, only(unknowns(sys1)))
    @test isequal(x1, only(species(sys1)))
    @test isequal(p1, only(parameters(sys1)))
    @named sys2 = ReactionSystem([], t; systems = [sys1])
    @test length(unknowns(sys2)) == 2
    @test any(isequal(x2), unknowns(sys2))
    @test any(isequal(x2), species(sys2))
    @test length(parameters(sys2)) == 2
    @test any(isequal(p2), parameters(sys2))
    @named sys3 = ReactionSystem(Equation[], t)
    sys3 = sys3 ∘ sys2
    @test length(unknowns(sys3)) == 3
    @test any(isequal(x3), unknowns(sys3))
    @test length(species(sys3)) == 2
    @test length(parameters(sys3)) == 3
    @test any(isequal(p3), parameters(sys3))
    sys4 = complete(sys3)
    @test length(unknowns(sys3)) == 3
    @test length(parameters(sys4)) == 4
    @test any(isequal(p5), parameters(sys4))
end
