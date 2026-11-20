### Prepares Tests ###

# Fetch packages.
using Catalyst, DynamicQuantities, JumpProcesses, OrdinaryDiffEq, SciMLBase, StaticArrays,
    StochasticDiffEq, Test
using ModelingToolkitBase: Pre
const MT = ModelingToolkitBase

# Sets the default `t` and `D` to use, and a stable rng seed.
t = default_t()
D = default_time_deriv()
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# The ODE drift of a problem at its initial state, keyed by unknown.
function drift_by_unknown(prob)
    du = prob.f(prob.u0, prob.p, 0.0)
    return Dict(zip(unknowns(prob.f.sys), du))
end

# The (legacy, matrix) SDE noise of a problem at its initial state, one row per unknown.
function noise_by_unknown(prob)
    g = prob.g(prob.u0, prob.p, 0.0)
    return Dict(u => g[i, :] for (i, u) in enumerate(unknowns(prob.f.sys)))
end

# Whether two expressions agree at the given values.
same_value(ex1, ex2, vals) = Symbolics.value(Symbolics.substitute(ex1 - ex2, vals)) ≈ 0


### Comparisons With Hand-Merged Models ###

# An aliased network gives the same ODE, SDE, jump and hybrid models as the hand-merged one.
# `2A + B --> E` with `A ~ B` is `3B --> E`, with that reaction's own mass action law.
let
    rn = @reaction_network begin
        @aliases begin
            A ~ B
            k1 ~ k2
        end
        k1, A --> C
        k2, B --> D
        k3, 2A + B --> E
        d, E --> 0
    end
    rn_merged = @reaction_network begin
        k2, B --> C
        k2, B --> D
        k3, 3B --> E
        d, E --> 0
    end
    @test issetequal(species(rn), species(rn_merged))
    @test isequal(reactions(rn), reactions(rn_merged))
    @test isequal(oderatelaw.(reactions(rn)), oderatelaw.(reactions(rn_merged)))
    @test isequal(jumpratelaw.(reactions(rn)), jumpratelaw.(reactions(rn_merged)))

    # The eliminated names can be used when building problems.
    u0 = [:A => 10, :C => 0, :D => 0, :E => 0]
    ps = [:k1 => 0.5, :k3 => 0.01, :d => 0.2]
    u0_merged = [:B => 10, :C => 0, :D => 0, :E => 0]
    ps_merged = [:k2 => 0.5, :k3 => 0.01, :d => 0.2]

    # ODE drift, SDE drift and noise.
    oprob = ODEProblem(rn, u0, (0.0, 10.0), ps)
    oprob_merged = ODEProblem(rn_merged, u0_merged, (0.0, 10.0), ps_merged)
    @test drift_by_unknown(oprob) == drift_by_unknown(oprob_merged)
    sprob = SDEProblem(rn, u0, (0.0, 10.0), ps)
    sprob_merged = SDEProblem(rn_merged, u0_merged, (0.0, 10.0), ps_merged)
    @test drift_by_unknown(sprob) == drift_by_unknown(sprob_merged)
    @test noise_by_unknown(sprob) == noise_by_unknown(sprob_merged)

    # Jump propensities and updates: identical SSA trajectories for the same seed.
    jprob = JumpProblem(rn, u0, (0.0, 10.0), ps)
    jprob_merged = JumpProblem(rn_merged, u0_merged, (0.0, 10.0), ps_merged)
    @test eltype(jprob.prob.u0) == Int
    jsol = solve(jprob, SSAStepper(); rng = StableRNG(seed))
    jsol_merged = solve(jprob_merged, SSAStepper(); rng = StableRNG(seed))
    @test jsol[rn.E] == jsol_merged[rn_merged.E]
    @test jsol[rn.A] == jsol[rn.B] == jsol_merged[rn_merged.B]

    # Hybrid dynamics (reaction 3 as jumps, the others as ODEs).
    hprob = HybridProblem(rn, u0, (0.0, 10.0), ps;
        default_scale = PhysicalScale.ODE, physical_scales = [3 => PhysicalScale.Jump])
    hprob_merged = HybridProblem(rn_merged, u0_merged, (0.0, 10.0), ps_merged;
        default_scale = PhysicalScale.ODE, physical_scales = [3 => PhysicalScale.Jump])
    hsol = solve(hprob, Tsit5(); rng = StableRNG(seed), saveat = 1.0)
    hsol_merged = solve(hprob_merged, Tsit5(); rng = StableRNG(seed), saveat = 1.0)
    @test hsol[rn.E] ≈ hsol_merged[rn_merged.E]

    # The eliminated parameter is bound: not tunable, and read as the one replacing it.
    osys = oprob.f.sys
    @test issetequal(MT.tunable_parameters(osys, parameters(osys)),
        MT.tunable_parameters(oprob_merged.f.sys, parameters(oprob_merged.f.sys)))
    @test oprob.ps[rn.k1] == oprob.ps[rn.k2] == 0.5
end

# The documented composition example: an alias between the composing system's own species and a
# subsystem's species.
let
    @species A(t) B(t) C(t) X(t)
    @parameters k1 k2
    @named sub1 = ReactionSystem([Reaction(k1, [A], [B])], t)
    @named sub2 = ReactionSystem([Reaction(k2, [C], [X])], t)
    sys = complete(compose(sub1, [sub2]; aliases = [sub1.A ~ sub2.C]))
    @test issetequal(unknowns(sys), [sys.B, sys.sub2.C, sys.sub2.X])

    prob = ODEProblem(sys, [sys.A => 1.0, sys.B => 0.0, sys.sub2.X => 0.0], (0.0, 1.0),
        [sys.k1 => 1.0, sys.sub2.k2 => 2.0])
    sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    @test sol[sys.sub2.C][end] ≈ exp(-3.0)
    @test sol[sys.A] == sol[sys.sub2.C]
    @test sol[sys.B][end] ≈ (1 - exp(-3.0)) / 3

    # The same aliasing, composing both systems into a parent.
    @named parent = ReactionSystem(Reaction[], t)
    sys2 = complete(compose(parent, [sub1, sub2]; aliases = [sub1.A ~ sub2.C]))
    @test issetequal(unknowns(sys2), [sys2.sub1.B, sys2.sub2.C, sys2.sub2.X])
end

# Alias chains, and querying aliases before and after completing (nothing is eliminated before).
let
    @species A(t) B(t) C(t) X(t)
    @parameters k1 k2 k3
    base = ReactionSystem([Reaction(k1, [A], [X]), Reaction(k2, [B], [X])], t; aliases = [A ~ B],
        name = :base)
    @named sub = ReactionSystem([Reaction(k3, [C], nothing)], t)
    @test isempty(alias_map(base))
    @test issetequal(unknowns(base), [A, B, X])
    @test isequal(alias_map(complete(base)), Dict(unwrap(A) => unwrap(B)))

    composed = compose(base, [sub]; aliases = [base.B ~ sub.C])
    sub_C = unwrap(sub.C)
    @test isempty(alias_map(composed))
    @test issetequal(species(composed), [unwrap(A), unwrap(B), unwrap(X), sub_C])
    @test length(aliases(composed)) == 2
    rs = complete(composed)
    @test isequal(alias_map(rs), Dict(unwrap(A) => sub_C, unwrap(B) => sub_C))
    @test issetequal(unknowns(rs), [unwrap(X), sub_C])
    @test isempty(aliases(rs))
end

# `flatten` keeps aliases (namespaced); completing the flattened system resolves them as
# completing the hierarchical one does.
let
    @species A(t) B(t) C(t) X(t)
    @parameters k1 k2 p
    @named sub = ReactionSystem([Reaction(k1, [A], [C]), Reaction(k2, [B], [C])], t;
        aliases = [A ~ B, k1 ~ k2])
    @named parent = ReactionSystem([Reaction(p, nothing, [X])], t; systems = [sub])
    flatparent = Catalyst.flatten(parent)
    @test isempty(MT.get_systems(flatparent))
    @test issetequal(Catalyst.get_aliases(flatparent), [sub.A ~ sub.B, sub.k1 ~ sub.k2])
    @test isempty(alias_map(flatparent))

    rs = complete(parent)
    @test isequal(alias_map(rs), Dict(unwrap(rs.sub.A) => unwrap(rs.sub.B), unwrap(rs.sub.k1) => unwrap(rs.sub.k2)))
    @test isequal(alias_map(complete(flatparent)), alias_map(rs))
    prob = ODEProblem(rs, [rs.sub.A => 2.0, rs.sub.C => 0.0, rs.X => 0.0], 1.0, [rs.sub.k1 => 1.0, rs.p => 1.0])
    @test prob[rs.sub.B] == 2.0
    @test prob.ps[rs.sub.k2] == 1.0
end

# A flattened system keeps its aliases when extended, and they are resolved when completing (in
# either order).
let
    @species C(t) D(t)
    @parameters p q r
    @named base = ReactionSystem([Reaction(p, [C], nothing), Reaction(q, [D], nothing)], t;
        aliases = [C ~ D, p ~ q])
    @named extra = ReactionSystem([Reaction(r, nothing, [D])], t)
    flatbase = Catalyst.flatten(base)
    for rs in (complete(extend(flatbase, extra)), complete(extend(extra, flatbase)))
        @test isequal(alias_map(rs), Dict(unwrap(C) => unwrap(D), unwrap(p) => unwrap(q)))
        prob = ODEProblem(rs, [C => 1.0], 1.0, [p => 1.0, r => 1.0])
        @test prob[D] == 1.0
        @test prob.ps[q] == prob.ps[p] == 1.0
    end
end

# Aliases given to `compose` may name the symbols of the composing system's existing subsystems
# through its name.
let
    @species E(t) F(t) G(t)
    @parameters p q r
    @named a = ReactionSystem([Reaction(p, [E], nothing)], t)
    @named b = ReactionSystem([Reaction(q, [F], nothing)], t)
    @named root = ReactionSystem([Reaction(r, nothing, [G])], t; systems = [a])
    for alias in (root.a.E ~ b.F, a.E ~ b.F)
        rs = complete(compose(root, [b]; aliases = [alias]))
        @test isequal(alias_map(rs), Dict(unwrap(rs.a.E) => unwrap(rs.b.F)))
    end
end


### Kinds of Aliased Symbols ###

# Constant species, boundary condition species and parameters.
let
    @species X(t) Y(t) B1(t) [isbcspecies = true] B2(t) [isbcspecies = true]
    @parameters k K1 [isconstantspecies = true] K2 [isconstantspecies = true]

    # Constant species are merged in reactions: `K1 + K2 + X --> Y` is `2K2 + X --> Y`.
    rs = complete(ReactionSystem([Reaction(k, [K1, K2, X], [Y])], t, [X, Y], [k, K1, K2];
        aliases = [K1 ~ K2], name = :rs))
    @test isequal(reactions(rs), [Reaction(k, [K2, X], [Y], [2, 1], [1])])
    @test same_value(oderatelaw(only(reactions(rs))), k * K2^2 * X / 2, Dict(k => 1.3, K2 => 2.0, X => 3.0))

    # Boundary condition species.
    rs = complete(ReactionSystem([Reaction(k, [B1], [B1, Y]), Reaction(k, [B2, Y], [B2])], t;
        aliases = [B1 ~ B2], name = :rs))
    @test issetequal(unknowns(rs), [B2, Y])

    # Parameters: both names stay readable, remake with the surviving name updates both, and
    # the eliminated name can't be given a value through remake.
    @parameters k1 k2
    rs = complete(ReactionSystem([Reaction(k1, [X], nothing), Reaction(k2, nothing, [X])], t;
        aliases = [k1 ~ k2], name = :rs))
    @test issetequal(parameters(rs), [k1, k2])
    @test issetequal(uneliminated_parameters(rs), [k2])
    prob = ODEProblem(rs, [X => 1.0], (0.0, 1.0), [k1 => 0.4])
    @test prob.ps[k1] == prob.ps[k2] == 0.4
    prob2 = remake(prob; p = [k2 => 5.0])
    @test prob2.ps[k1] == prob2.ps[k2] == 5.0
    @test_throws Exception remake(prob; p = [k1 => 5.0])
end

# Non-species variables: the eliminated variable's equation is rewritten; two defining
# equations for the same variable are an error.
let
    @species X(t)
    @variables V(t) W(t)
    @parameters k
    rs = complete(ReactionSystem([Reaction(k, nothing, [X]), D(V) ~ -V], t, [X, V, W], [k];
        aliases = [V ~ W], name = :rs))
    @test issetequal(unknowns(rs), [X, W])
    @test any(isequal(D(W) ~ -W), equations(rs))
    @test_throws ArgumentError complete(ReactionSystem(
        [Reaction(k, nothing, [X]), D(V) ~ -V, D(W) ~ 1 - W], t; aliases = [V ~ W], name = :rs))
end

# Parametric stoichiometry: aliased stoichiometry parameters are substituted.
let
    @species A(t) B(t) C(t)
    @parameters k n1::Int n2::Int
    rs = complete(ReactionSystem([Reaction(k, [A], [B], [n1], [1]), Reaction(k, [B], [C], [n2], [1])],
        t, [A, B, C], [k, n1, n2]; aliases = [n1 ~ n2], name = :rs))
    @test all(rx -> isequal(only(rx.substoich), unwrap(n2)), reactions(rs))
    @test !any(eq -> occursin("n1", string(eq)), equations(ode_model(rs)))
end


### Reactions That Become Self-Loops ###

# An alias that turns a reaction into a self-loop is an error naming the reaction; a self-loop
# that existed before is kept. Reactions are never dropped, so reaction indices (e.g. in
# `physical_scales`) stay valid.
let
    @species A(t) B(t) C(t)
    @parameters k1 k2 k3
    rs = ReactionSystem([Reaction(k1, [B], [C]), Reaction(k2, [A], [B])], t; aliases = [A ~ B], name = :rs)
    err = try
        complete(rs)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("reaction 2", err.msg)
    @test occursin("A(t) ~ B(t)", err.msg)

    rs = complete(ReactionSystem([Reaction(k1, [B], [B]), Reaction(k2, [B], [C])], t, [B, C], [k1, k2, k3];
        aliases = [k3 ~ k2], name = :rs))
    @test length(reactions(rs)) == 2

    rs = complete(ReactionSystem([Reaction(k1, [A], [C]), Reaction(k2, [B], [C])], t; aliases = [A ~ B], name = :rs))
    @test length(reactions(rs)) == 2
    sys = hybrid_model(rs; default_scale = PhysicalScale.ODE, physical_scales = [2 => PhysicalScale.Jump])
    @test isequal(only(MT.jumps(sys)).rate, unwrap(k2 * B))
end


### Other System Components ###

# Events (one that fires, and the settings of a rebuilt one), observables, tstops, noise scaling
# and poissonian rates use the surviving symbols.
let
    rn = @reaction_network begin
        @aliases A ~ B
        @observables Btot ~ A + C
        @discrete_events [1.0] => [A => A + 5.0]
        @continuous_events [A ~ 100.0] => [A => 0.0]
        (k, d), 0 <--> A
        k, B --> C
    end
    sol = solve(ODEProblem(rn, [rn.B => 1.0, rn.C => 0.0], (0.0, 2.0), [rn.k => 0.0, rn.d => 0.0]), Tsit5())
    @test sol(0.5; idxs = rn.B) ≈ 1.0
    @test sol(1.5; idxs = rn.B) ≈ 6.0
    @test sol(1.5; idxs = rn.A) ≈ 6.0
    @test isequal(MT.observed(rn)[1].rhs, unwrap(rn.B + rn.C))

    # A rebuilt continuous event keeps its settings.
    @species X(t) Y(t)
    @parameters p
    cevent = MT.SymbolicContinuousCallback([X ~ 2.0], [X ~ 0.0]; affect_neg = [X ~ 1.0],
        rootfind = SciMLBase.RightRootFind)
    rs = complete(ReactionSystem([Reaction(p, nothing, [Y])], t, [X, Y], [p];
        aliases = [X ~ Y], continuous_events = [cevent], name = :rs))
    evt = only(MT.get_continuous_events(rs))
    @test evt.rootfind == SciMLBase.RightRootFind
    @test isequal(evt.conditions, [Y ~ 2.0])
    @test isequal(evt.affect_neg.affect, [Y ~ 1.0])

    # Noise scaling metadata and tstops.
    @parameters η1 η2 τ1 τ2
    rx = Reaction(p, [X], nothing; metadata = [:noise_scaling => η1])
    rs = complete(ReactionSystem([rx], t, [X], [p, η1, η2, τ1, τ2]; aliases = [η1 ~ η2, τ1 ~ τ2],
        tstops = [τ1], name = :rs))
    @test isequal(Catalyst.getnoisescaling(only(reactions(rs))), unwrap(η2))
    @test isequal(only(MT.get_tstops(rs)), unwrap(τ2))

    # Poissonian rates.
    rn = @reaction_network begin
        @species X(t)
        @aliases X ~ Y
        @poissonians dN(1 + X)
        @variables V(t)
        @equations D(V) ~ dN
        k, Y --> 0
    end
    @test isequal(MT.getpoissonianrate(only(MT.get_poissonians(rn))), unwrap(1 + rn.Y))
end

# Initial conditions: values are resolved before being compared.
let
    @species X(t) Y(t)
    @parameters k1 k2 d
    rs = complete(ReactionSystem([Reaction(d, [Y], nothing)], t, [X, Y], [k1, k2, d];
        aliases = [X ~ Y, k1 ~ k2], initial_conditions = Dict(X => k1, Y => k2), name = :rs))
    @test isequal(MT.initial_conditions(rs)[unwrap(Y)], unwrap(k2))
    @test_throws ArgumentError complete(ReactionSystem([Reaction(d, [Y], nothing)], t, [X, Y], [d];
        aliases = [X ~ Y], initial_conditions = Dict(X => 1.0, Y => 2.0), name = :rs))
end

# A default value of an eliminated parameter moves to the parameter replacing it (the eliminated
# parameter, now bound, keeps none); different defaults for the two are an error.
let
    @species X(t)
    @parameters k1 = 0.5 k2 k3 = 0.7
    rxs = [Reaction(k1, [X], nothing), Reaction(k2, nothing, [X])]
    rs = complete(ReactionSystem(rxs, t; aliases = [k1 ~ k2], name = :rs))
    prob = ODEProblem(rs, [X => 1.0], 1.0)
    @test prob.ps[k1] == prob.ps[k2] == 0.5
    @test remake(prob; p = [k2 => 2.0]).ps[k1] == 2.0
    @test_throws ArgumentError complete(ReactionSystem(rxs, t; aliases = [k1 ~ k3], name = :rs))
end


### Problem Inputs ###

# Eliminated names are accepted in every input form, conflicting values are an error, and
# remake uses the surviving names.
let
    rn = @reaction_network begin
        @species B(t)
        @aliases A ~ B
        (k1, k2), A <--> C
    end
    tspan = (0.0, 1.0)
    ps = [:k1 => 1.0, :k2 => 2.0]
    for u0 in ([rn.A => 3.0, rn.C => 1.0], SA[rn.A => 3.0, rn.C => 1.0], Dict([rn.A => 3.0, rn.C => 1.0]),
            (rn.A => 3.0, rn.C => 1.0), [:A => 3.0, :C => 1.0], (:A => 3.0, :C => 1.0))
        prob = ODEProblem(rn, u0, tspan, ps)
        @test prob[rn.B] == 3.0
        @test prob[rn.C] == 1.0
    end
    for u0 in ([rn.A => 3.0, rn.B => 2.0, rn.C => 1.0], [:A => 3.0, :B => 2.0, :C => 1.0],
            [rn.A => 3.0, :B => 2.0, :C => 1.0], (:B => 2.0, rn.A => 3.0, :C => 1.0),
            Dict(:A => 3.0, :B => 2.0, :C => 1.0))
        @test_throws ArgumentError ODEProblem(rn, u0, tspan, ps)
    end
    @test ODEProblem(rn, [rn.A => 3.0, rn.B => 3.0, rn.C => 1.0], tspan, ps)[rn.B] == 3.0
    @test ODEProblem(rn, [:A => 3.0, :B => 3.0, :C => 1.0], tspan, ps)[rn.B] == 3.0

    # Tuples keep their value types (integers for the species here).
    jprob = JumpProblem(rn, (rn.A => 3, rn.C => 1), tspan, (rn.k1 => 1.0, rn.k2 => 2.0))
    @test eltype(jprob.prob.u0) == Int

    # remake with the surviving symbols.
    prob = ODEProblem(rn, [rn.A => 3.0, rn.C => 1.0], tspan, ps)
    @test issetequal(unknowns(rn), [rn.B, rn.C])
    @test remake(prob; u0 = [rn.B => 5.0])[rn.B] == 5.0

    # The two-stage workflow (an MTK problem constructor) takes the surviving symbols.
    osys = complete(ode_model(rn))
    prob = ODEProblem(osys, [rn.B => 3.0, rn.C => 1.0, rn.k1 => 1.0, rn.k2 => 2.0], tspan)
    @test prob[rn.A] == 3.0
end


### Bindings ###

# Eliminated discrete parameters track updates; bindings survive flattening and conversions.
let
    rn = @reaction_network begin
        @discretes d1(t) d2(t)
        @aliases d1 ~ d2
        @discrete_events [1.0] => [d2 => 3.0]
        d1, X --> 0
    end
    prob = ODEProblem(rn, [:X => 1.0], (0.0, 2.0), [:d2 => 1.0])
    @test prob.ps[rn.d1] == prob.ps[rn.d2] == 1.0
    sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    @test sol[rn.X][end] ≈ exp(-4.0)
    @test sol.ps[rn.d2] == [1.0, 3.0]
    @test sol.ps[rn.d1] == 3.0

    @species X(t)
    @parameters a b c
    @named sub = ReactionSystem([Reaction(a, [X], nothing), Reaction(b + c, nothing, [X])], t, [X], [a, b, c];
        bindings = [b => 2a, c => a + 1])
    @named parent = ReactionSystem(Reaction[], t; systems = [sub])
    @test length(MT.get_bindings(Catalyst.flatten(parent))) == 2
    rs = complete(parent)
    for sys in (ode_model(rs), sde_model(rs), jump_model(rs), ss_ode_model(rs))
        @test length(MT.bindings(sys)) == 2
    end
end


### Declared (Unresolved) Systems and Errors ###

# Unresolved systems list both names; conservation laws need resolved aliases.
let
    nc = @network_component begin
        @aliases A ~ B
        (k1, k2), A <--> C
        k3, B --> C
    end
    @test issetequal(MT.getname.(species(nc)), [:A, :B, :C])
    @test Catalyst.aliases_present(nc)
    @test_throws ArgumentError conservationlaws(nc)
    @test size(conservationlaws(complete(nc)), 1) == 1
end

# A system completed without flattening keeps its aliases; simulating it resolves them, also
# accepting eliminated names.
let
    @species A(t) B(t) C(t)
    @parameters k1 k2
    @named sub = ReactionSystem([Reaction(k1, [A], [C])], t)
    @named rs = ReactionSystem([Reaction(k2, [B], [C])], t; systems = [sub], aliases = [sub.A ~ B])
    crs = complete(rs; flatten = false)
    @test Catalyst.aliases_present(crs)
    @test isempty(alias_map(crs))
    u0 = [crs.sub.A => 2.0, crs.sub.C => 0.0, crs.C => 0.0]
    ps = [crs.sub.k1 => 1.0, crs.k2 => 3.0]
    prob = ODEProblem(crs, u0, 1.0, ps)
    @test prob[crs.B] == 2.0
    sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    @test sol[crs.B][end] ≈ 2.0 * exp(-4.0)
end

# Invalid aliases are reported when the system is completed.
let
    @species A(t) B(t) [unit = us"m"] C(t)
    @variables V(t)
    @parameters k k2
    rx = Reaction(k, [A], [C])
    for aliases in ([A ~ k], [A ~ V], [A ~ B], [A ~ C, C ~ A])
        @test_throws ArgumentError complete(ReactionSystem([rx, Reaction(k, [B], nothing), D(V) ~ -V], t;
            aliases, name = :rs))
    end
    @test_throws "not a species, variable or parameter" complete(ReactionSystem([rx], t, [A, C], [k];
        aliases = [A ~ B], name = :rs))
    @test_throws ArgumentError ReactionSystem([rx], t; aliases = [A ~ A], name = :rs)
    @test_throws ArgumentError ReactionSystem([rx], t; aliases = [A ~ C, A ~ B], name = :rs)
    @test_throws ArgumentError complete(ReactionSystem([rx], t, [A, C], [k, k2];
        jumps = [MassActionJump([k], [[A => 1]], [[A => -1]])], aliases = [k ~ k2], name = :rs))
end

# Serialization gives a clear error (spatial systems are tested in
# spatial_modelling/dspace_reaction_systems.jl).
let
    rn = @reaction_network begin
        @aliases A ~ B
        (k1, k2), A <--> B2
        k3, B --> 0
    end
    @test_throws Exception save_reactionsystem(joinpath(mktempdir(), "rn.jl"), rn)
end

# Symbols in `@aliases` must be part of the network.
let
    @test_throws Exception @eval @reaction_network begin
        @aliases A ~ B
        k, A --> 0
    end
end

# `@aliases` in a `@network_component`, then composed.
let
    nc = @network_component begin
        @aliases A ~ B
        k1, A --> C
        k2, B --> 0
    end
    @species Z(t)
    @parameters kz
    @named sub = ReactionSystem([Reaction(kz, nothing, [Z])], t)
    sys = complete(compose(nc, [sub]))
    @test issetequal(MT.getname.(unknowns(sys)), [:B, :C, Symbol("sub₊Z")])
    @test Dict(MT.getname(k) => MT.getname(v) for (k, v) in alias_map(sys)) == Dict(:A => :B)
end
