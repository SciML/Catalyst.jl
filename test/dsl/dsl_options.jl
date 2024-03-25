### Prepares Tests ###

### Prepares Tests ###

# Fetch packages.
using Catalyst, ModelingToolkit, OrdinaryDiffEq, StochasticDiffEq, Plots, Test
using Symbolics: unwrap

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# Sets the default `t` to use.
t = default_t()

### Tests `@parameters`, `@species`, and `@variables` Options ###

# Sets the default `t` to use.
t = default_t()

# Fetch test networks and functions.
include("../test_networks.jl")
include("../test_functions.jl")

### Declares Testing Functions ###

function unpacksys(sys)
    get_eqs(sys), get_iv(sys), get_unknowns(sys), get_ps(sys), nameof(sys), get_systems(sys)
end

opname(x) = istree(x) ? nameof(operation(x)) : nameof(x)
alleq(xs, ys) = all(isequal(x, y) for (x, y) in zip(xs, ys))

# Gets all the reactants in a set of equations.
function all_reactants(eqs)
    all_reactants = []
    for eq in eqs
        append!(all_reactants, opname.(eq.substrates))
        append!(all_reactants, opname.(eq.products))
    end
    return Set{Symbol}(unique(all_reactants))
end

# Gets all parameters (where every reaction rate is constant)
function all_parameters(eqs)
    return Set(unique(map(eq -> opname(eq.rate), eqs)))
end

# Perform basic tests.
function basic_test(rn, N, unknowns_syms, p_syms)
    eqs, iv, unknowns, ps, name, systems = unpacksys(rn)
    @test length(eqs) == N
    @test opname(iv) == :t
    @test length(unknowns) == length(unknowns_syms)
    @test issetequal(map(opname, unknowns), unknowns_syms)
    @test all_reactants(eqs) == Set(unknowns_syms)
    @test length(ps) == length(p_syms)
    @test issetequal(map(opname, ps), p_syms)
end

### Basic Tests ###

# Test basic properties of networks.
let
    basic_test(reaction_networks_standard[1], 10, [:X1, :X2, :X3],
               [:p1, :p2, :p3, :k1, :k2, :k3, :k4, :d1, :d2, :d3])
    @test all_parameters(get_eqs(reaction_networks_standard[1])) ==
          Set([:p1, :p2, :p3, :k1, :k2, :k3, :k4, :d1, :d2, :d3])
    basic_test(reaction_networks_standard[2], 3, [:X1, :X2], [:v1, :K1, :v2, :K2, :d])
    basic_test(reaction_networks_standard[3], 10, [:X1, :X2, :X3, :X4],
               [:v1, :K1, :v2, :K2, :k1, :k2, :k3, :k4, :d])
    basic_test(reaction_networks_standard[4], 8, [:X1, :X2, :X3, :X4],
               [:v1, :K1, :v2, :K2, :v3, :K3, :v4, :K4, :d1, :d2, :d3, :d4])
    basic_test(reaction_networks_standard[5], 8, [:X1, :X2, :X3, :X4],
               [:p, :k1, :k2, :k3, :k4, :k5, :k6, :d])
    @test all_parameters(get_eqs(reaction_networks_standard[5])) ==
          Set([:p, :k1, :k2, :k3, :k4, :k5, :k6, :d])
    basic_test(reaction_networks_hill[1], 4, [:X1, :X2],
               [:v1, :v2, :K1, :K2, :n1, :n2, :d1, :d2])
    basic_test(reaction_networks_constraint[1], 6, [:X1, :X2, :X3],
               [:k1, :k2, :k3, :k4, :k5, :k6])
    basic_test(reaction_networks_real[1], 4, [:X, :Y], [:A, :B])
    basic_test(reaction_networks_weird[1], 2, [:X], [:p, :d])
    basic_test(reaction_networks_weird[2], 4, [:X, :Y, :Z], [:k1, :k2, :k3, :k4])
end

# Compares networks to networks created using different arrow types.
let
    networks_1 = []
    networks_2 = []

    different_arrow_1 = @reaction_network begin
        (p1, p2, p3), ∅ > (X1, X2, X3)
        (k1, k2), X2 ↔ X1 + 2X3
        (k3, k4), X1 ⟷ X3
        (d1, d2, d3), (X1, X2, X3) → ∅
    end
    push!(networks_1, reaction_networks_standard[1])
    push!(networks_2, different_arrow_1)

    different_arrow_2 = @reaction_network begin
        mmr(X2, v1, K1), ∅ → X1
        mm(X1, v2, K2), ∅ ↣ X2
        d, X1 + X2 ↦ ∅
    end
    push!(networks_1, reaction_networks_standard[2])
    push!(networks_2, different_arrow_2)

    different_arrow_3 = @reaction_network begin
        mm(X2, v1, K1), ∅ ⇾ X1
        mm(X3, v2, K2), ∅ ⟶ X2
        (k1, k2), X1 ⇄ X3
        (k3, k4), X3 + X2 ⇆ X4 + X1
        d, (X1, X2, X3, X4) ⟼ ∅
    end
    push!(networks_1, reaction_networks_standard[3])
    push!(networks_2, different_arrow_3)

    different_arrow_4 = @reaction_network begin
        mmr(X4, v1, K1), ∅ ⥟ X1
        mmr(X1, v2, K2), ∅ ⥟ X2
        mmr(X2, v3, K3), ∅ ⇀ X3
        mmr(X3, v4, K4), ∅ ⇁ X4
        (d1, d2, d3, d4), (X1, X2, X3, X4) --> ∅
    end
    push!(networks_1, reaction_networks_standard[4])
    push!(networks_2, different_arrow_4)

    # Yes the name is different, I wanted one with several single direction arrows.
    different_arrow_8 = @reaction_network begin
        p, 2X1 < ∅
        k1, X2 ← X1
        (k2, k3), X3 ⟻ X2
        d, ∅ ↼ X3
    end
    push!(networks_1, reaction_networks_standard[8])
    push!(networks_2, different_arrow_8)

    for (rn_1, rn_2) in zip(networks_1, networks_2)
        for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
            u0 = rnd_u0(rn_1, rng; factor)
            p = rnd_ps(rn_1, rng; factor)
            t = rand(rng)
            
            @test f_eval(rn_1, u0, p, t) ≈ f_eval(rn_2, u0, p, t)
            @test jac_eval(rn_1, u0, p, t) ≈ jac_eval(rn_2, u0, p, t)
            @test g_eval(rn_1, u0, p, t) ≈ g_eval(rn_2, u0, p, t)
        end
    end
end

# Checks that some created networks are identical.
let
    # Declares network.
    differently_written_5 = @reaction_network begin
        q, ∅ → Y1
        (l1, l2), Y1 ⟷ Y2
        (l3, l4), Y2 ⟷ Y3
        (l5, l6), Y3 ⟷ Y4
        c, Y4 → ∅
    end

    # Computes initial conditions/parameter values.
    u0_vals = rand(rng, length(species(differently_written_5)))
    ps_vals = rand(rng, length(parameters(differently_written_5)))
    u0_1 = [:X1 => u0_vals[1], :X2 => u0_vals[2], :X3 => u0_vals[3], :X4 => u0_vals[4]]
    u0_2 = [:Y1 => u0_vals[1], :Y2 => u0_vals[2], :Y3 => u0_vals[3], :Y4 => u0_vals[4]]
    ps_1 = [:p  => ps_vals[1], :k1  => ps_vals[2], :k2  => ps_vals[3], :k3  => ps_vals[4], 
            :k4  => ps_vals[5], :k5  => ps_vals[6], :k6  => ps_vals[7], :d => ps_vals[8]]
    ps_2 = [:q  => ps_vals[1], :l1  => ps_vals[2], :l2  => ps_vals[3], :l3  => ps_vals[4], 
            :l4  => ps_vals[5], :l5  => ps_vals[6], :l6  => ps_vals[7], :c => ps_vals[8]]
    t = rand(rng)

    # Checks equivalence.
    rn_1 = reaction_networks_standard[5]
    rn_2 = differently_written_5
    @test f_eval(rn_1, u0_1, ps_1, t) ≈ f_eval(rn_2, u0_2, ps_2, t)
    @test jac_eval(rn_1, u0_1, ps_1, t) ≈ jac_eval(rn_2, u0_2, ps_2, t)
    @test g_eval(rn_1, u0_1, ps_1, t) ≈ g_eval(rn_2, u0_2, ps_2, t)
end

# Compares networks to networks written in different ways.
let
    networks_1 = []
    networks_2 = []

    # Unfold reactions.
    differently_written_6 = @reaction_network begin
        p1, ∅ → X1
        p2, ∅ → X2
        k1, 2X1 → X3
        k2, X3 → 2X1
        k3, X2 + X3 → 4X4
        k4, 4X4 → X2 + X3
        k5, X4 + X1 → 2X3
        k6, 2X3 → X4 + X1
        d, X1 → ∅
        d, X2 → ∅
        d, X3 → ∅
        d, X4 → ∅
    end
    push!(networks_1, reaction_networks_standard[6])
    push!(networks_2, differently_written_6)

    # Ignore mass action.
    differently_written_7 = @reaction_network begin
        @parameters p1 p2 p3 k1 k2 k3 v1 K1 d1 d2 d3 d4 d5
        (p1, p2, p3), ∅ ⇒ (X1, X2, X3)
        (k1 * X1 * X2^2 / 2, k2 * X4), X1 + 2X2 ⟺ X4
        (mm(X3, v1, K1) * X4, k3 * X5), X4 ⇔ X5
        (d1 * X1, d2 * X2, d3 * X3, d4 * X4, d5 * X5), ∅ ⟽ (X1, X2, X3, X4, X5)
    end
    push!(networks_1, reaction_networks_standard[7])
    push!(networks_2, differently_written_7)

    # Ignore mass action new arrows.
    differently_written_8 = @reaction_network begin
        @parameters p1 p2 p3 k1 k2 k3 v1 K1 d1 d2 d3 d4 d5
        (p1, p2, p3), ∅ => (X1, X2, X3)
        (k1 * X1 * X2^2 / 2, k2 * X4), X1 + 2X2 ⟺ X4
        (mm(X3, v1, K1) * X4, k3 * X5), X4 ⇔ X5
        (d1 * X1, d2 * X2, d3 * X3, d4 * X4, d5 * X5), ∅ <= (X1, X2, X3, X4, X5)
    end
    push!(networks_1, reaction_networks_standard[7])
    push!(networks_2, differently_written_8)

    @test ModelingToolkit.defaults(rn27) == defs29
    @test merge(ModelingToolkit.defaults(rn28), defs28) == ModelingToolkit.defaults(rn27)
end

# Tests that parameter type designation works.
let
    # Creates a model.
    rn = @reaction_network begin
        @parameters begin
            k1
            l1
            k2::Float64 = 2.0
            l2::Float64
            k3::Int64 = 2, [description="A parameter"]
            l3::Int64
            k4::Float32, [description="Another parameter"]
            l4::Float32
            k5::Rational{Int64}
            l5::Rational{Int64}
        end
        (k1,l1), X1 <--> Y1
        (k2,l2), X2 <--> Y2
        (k3,l3), X3 <--> Y3
        (k4,l4), X4 <--> Y4
        (k5,l5), X5 <--> Y5
    end

    # Checks parameter types.
    @test unwrap(rn.k1) isa SymbolicUtils.BasicSymbolic{Real}
    @test unwrap(rn.l1) isa SymbolicUtils.BasicSymbolic{Real}
    @test unwrap(rn.k2) isa SymbolicUtils.BasicSymbolic{Float64}
    @test unwrap(rn.l2) isa SymbolicUtils.BasicSymbolic{Float64}
    @test unwrap(rn.k3) isa SymbolicUtils.BasicSymbolic{Int64}
    @test unwrap(rn.l3) isa SymbolicUtils.BasicSymbolic{Int64}
    @test unwrap(rn.k4) isa SymbolicUtils.BasicSymbolic{Float32}
    @test unwrap(rn.l4) isa SymbolicUtils.BasicSymbolic{Float32}
    @test unwrap(rn.k5) isa SymbolicUtils.BasicSymbolic{Rational{Int64}}
    @test unwrap(rn.l5) isa SymbolicUtils.BasicSymbolic{Rational{Int64}}

    # Checks that other parameter properties are assigned properly.
    @test !ModelingToolkit.hasdefault(unwrap(rn.k1))
    @test ModelingToolkit.getdefault(unwrap(rn.k2)) == 2.0
    @test ModelingToolkit.getdefault(unwrap(rn.k3)) == 2
    @test ModelingToolkit.getdescription(unwrap(rn.k3)) == "A parameter"
    @test ModelingToolkit.getdescription(unwrap(rn.k4)) == "Another parameter"
    @test !ModelingToolkit.hasdescription(unwrap(rn.k5))
end

# Test @variables in DSL.
let
    rn = @reaction_network tester begin
        @parameters k1
        @variables V1(t) V2(t) V3(t)
        @species B1(t) B2(t)
        (k1*k2 + V3), V1*A + 2*B1 --> V2*C + B2
    end

    @parameters k1 k2
    @variables V1(t) V2(t) V3(t)
    @species A(t) B1(t) B2(t) C(t)
    rx = Reaction(k1*k2 + V3, [A, B1], [C, B2], [V1, 2], [V2, 1])
    @named tester = ReactionSystem([rx], t)
    @test tester == rn

    sts = (A, B1, B2, C, V1, V2, V3)
    spcs = (A, B1, B2, C)
    @test issetequal(unknowns(rn), sts)
    @test issetequal(species(rn), spcs)

    @test_throws ArgumentError begin
        rn = @reaction_network begin
            @variables K
            k, K*A --> B
        end
    end
end

### Test Independent Variable Designations ###

# Test ivs in DSL.
let
    rn = @reaction_network ivstest begin
        @ivs s x
        @parameters k2
        @variables D(x) E(s) F(s,x)
        @species A(s,x) B(s) C(x)
        k*k2*D, E*A +B --> F*C + C2
    end

    @parameters k k2
    @variables s x D(x) E(s) F(s,x)
    @species A(s,x) B(s) C(x) C2(s,x)
    rx = Reaction(k*k2*D, [A, B], [C, C2], [E, 1], [F, 1])
    @named ivstest = ReactionSystem([rx], s; spatial_ivs = [x])

    @test ivstest == rn
    @test issetequal(unknowns(rn), [D, E, F, A, B, C, C2])
    @test issetequal(species(rn), [A, B, C, C2])
    @test isequal(ModelingToolkit.get_iv(rn), s)
    @test issetequal(Catalyst.get_sivs(rn), [x])
end


### Observables ###

# Test basic functionality.
# Tests various types of indexing.
let
    rn = @reaction_network begin
        @observables begin
            X ~ Xi + Xa
            Y ~ Y1 + Y2
        end
    end
    @unpack X, Xi, Xa, Y, Y1, Y2, p, d, k1, k2, k3, k4 = rn

    # Test that ReactionSystem have the correct properties.
    @test length(species(rn)) == 4
    @test length(unknowns(rn)) == 4
    @test length(observed(rn)) == 2
    @test length(equations(rn)) == 6

    @test isequal(observed(rn)[1], X ~ Xi + Xa)
    @test isequal(observed(rn)[2], Y ~ Y1 + Y2)

    # Tests correct indexing of solution.
    u0 = [Xi => 0.0, Xa => 0.0, Y1 => 1.0, Y2 => 2.0]
    ps = [p => 1.0, d => 0.2, k1 => 1.5, k2 => 1.5, k3 => 5.0, k4 => 5.0]

    oprob = ODEProblem(rn, u0, (0.0, 1000.0), ps)
    sol = solve(oprob, Tsit5())
    @test sol[X][end] ≈ 10.0
    @test sol[Y][end] ≈ 3.0
    @test sol[rn.X][end] ≈ 10.0
    @test sol[rn.Y][end] ≈ 3.0
    @test sol[:X][end] ≈ 10.0
    @test sol[:Y][end] ≈ 3.0

    # Tests that observables can be used for plot indexing.
    @test_broken false # plot(sol; idxs=X).series_list[1].plotattributes[:y][end] ≈ 10.0
    @test plot(sol; idxs=rn.X).series_list[1].plotattributes[:y][end] ≈ 10.0
    @test plot(sol; idxs=:X).series_list[1].plotattributes[:y][end] ≈ 10.0
    @test plot(sol; idxs=[X, Y]).series_list[2].plotattributes[:y][end] ≈ 3.0
    @test plot(sol; idxs=[rn.X, rn.Y]).series_list[2].plotattributes[:y][end] ≈ 3.0
    @test_broken false # plot(sol; idxs=[:X, :Y]).series_list[2].plotattributes[:y][end] ≈ 3.0
end

# Compares networks to networks written without parameters,
let
    networks_1 = []
    networks_2 = []
    parameter_sets = []

    # Different parameter and variable names.
    no_parameters_9 = @reaction_network begin
        (1.5, 1, 2), ∅ ⟶ (X1, X2, X3)
        (0.01, 2.3, 1001), (X1, X2, X3) ⟶ ∅
        (π, 42), X1 + X2 ⟷ X3
        (19.9, 999.99), X3 ⟷ X4
        (sqrt(3.7), exp(1.9)), X4 ⟷ X1 + X2
    end
    push!(networks_1, reaction_networks_standard[9])
    push!(networks_2, no_parameters_9)
    push!(parameter_sets, [:p1 => 1.5, :p2 => 1, :p3 => 2, :d1 => 0.01, :d2 => 2.3, :d3 => 1001, 
                           :k1 => π, :k2 => 42, :k3 => 19.9, :k4 => 999.99, :k5 => sqrt(3.7), :k6 => exp(1.9)])

    @variables X(t) Y(t)
    @species x(t), y(t), x2y(t)
    @parameters k kB kD d
    r1 = Reaction(k, nothing, [x], nothing, [1])
    r2 = Reaction(k, nothing, [y], nothing, [1])
    r3 = Reaction(kB, [x, y], [x2y], [2, 1], [1])
    r4 = Reaction(kD, [x2y], [x, y], [1], [2, 1])
    r5 = Reaction(d, [x], nothing, [1], nothing)
    r6 = Reaction(d, [y], nothing, [1], nothing)
    r7 = Reaction(d, [x2y], nothing, [1], nothing)
    obs_eqs = [X ~ x + 2x2y, Y ~ y + x2y]
    @named rn_prog = ReactionSystem([r1, r2, r3, r4, r5, r6, r7], t, [x, y, x2y], [k, kB, kD, d]; observed = obs_eqs)
    rn_prog = complete(rn_prog)


    for (rn_1, rn_2, p_1) in zip(networks_1, networks_2, parameter_sets)
        for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
            u0 = rnd_u0(rn_1, rng; factor)
            t = rand(rng)
            
            @test f_eval(rn_1, u0, p_1, t) ≈ f_eval(rn_2, u0, [], t)
            @test jac_eval(rn_1, u0, p_1, t) ≈ jac_eval(rn_2, u0, [], t)
            @test g_eval(rn_1, u0, p_1, t) ≈ g_eval(rn_2, u0, [], t)
        end
    end
end

# Tests for complicated observable formula.
# Tests using a single observable (without begin/end statement).
# Tests using observable component not part of reaction.
# Tests using parameters in observables formula.
let
    rn = @reaction_network begin
        @parameters op_1 op_2
        @species X4(t)
        @observables X ~ X1^2 + op_1*(X2 + 2X3) + X1*X4/op_2 + p
        (p,d), 0 <--> X1
        (k1,k2), X1 <--> X2
        (k3,k4), X2 <--> X3
    end

    u0 = Dict([:X1 => 1.0, :X2 => 2.0, :X3 => 3.0, :X4 => 4.0])
    ps = Dict([:p => 1.0, :d => 0.2, :k1 => 1.5, :k2 => 1.5, :k3 => 5.0, :k4 => 5.0, :op_1 => 1.5, :op_2 => 1.5])

    rxs_1 = [Reaction(p, nothing, [X1], nothing, [2]),
        Reaction(k1, [X1], [X2], [1], [1]),
        Reaction(k2, [X2], [X3], [1], [1]),
        Reaction(k3, [X2], [X3], [1], [1]),
        Reaction(d, [X3], nothing, [1], nothing)]
    @named rs_1 = ReactionSystem(rxs_1, t, [X1, X2, X3], [p, k1, k2, k3, d])
    push!(identical_networks_4, reaction_networks_standard[8] => rs_1)

    @test sol[:X][1] == u0[:X1]^2 + ps[:op_1]*(u0[:X2] + 2*u0[:X3]) + u0[:X1]*u0[:X4]/ps[:op_2] + ps[:p]
end

    rxs_3 = [Reaction(k1, [X1], [X2], [1], [1]),
        Reaction(0, [X2], [X3], [1], [1]),
        Reaction(k2, [X3], [X4], [1], [1]),
        Reaction(k3, [X4], [X5], [1], [1])]
    @named rs_3 = ReactionSystem(rxs_3, t, [X1, X2, X3, X4, X5], [k1, k2, k3])
    push!(identical_networks_4, reaction_networks_weird[7] => rs_3)

    for networks in identical_networks_4
        @test isequal(get_iv(networks[1]), get_iv(networks[2]))
        @test alleq(get_unknowns(networks[1]), get_unknowns(networks[2]))
        @test alleq(get_ps(networks[1]), get_ps(networks[2]))
        @test ModelingToolkit.get_systems(networks[1]) ==
              ModelingToolkit.get_systems(networks[2])
        @test length(get_eqs(networks[1])) == length(get_eqs(networks[2]))
        for (e1, e2) in zip(get_eqs(networks[1]), get_eqs(networks[2]))
            @test isequal(e1.rate, e2.rate)
            @test isequal(e1.substrates, e2.substrates)
            @test isequal(e1.products, e2.products)
            @test isequal(e1.substoich, e2.substoich)
            @test isequal(e1.prodstoich, e2.prodstoich)
            @test isequal(e1.netstoich, e2.netstoich)
            @test isequal(e1.only_use_rate, e2.only_use_rate)
        end
    end
end

### Tests Usage of Various Symbols ###

# Tests that time is handled properly.
let
    time_network = @reaction_network begin
        (t, k2), X1 ↔ X2
        (k3, t), X2 ↔ X3
        (t, k6), X3 ↔ X1
    end
    @test ModelingToolkit.getdescription(observed(rn)[1].lhs) == "my_description"
end

# Declares observables implicitly/explicitly.
# Cannot test `isequal(rn1, rn2)` because the two sets of observables have some obscure Symbolics
# substructure that is different.
let
    # Basic case.
    rn1 = @reaction_network rn_observed begin
        @observables X ~ X1 + X2
        k, 0 --> X1 + X2
    end
    rn2 = @reaction_network rn_observed begin
        @variables X(t)
        @observables X ~ X1 + X2
        k, 0 --> X1 + X2
    end
    @test isequal(observed(rn1)[1].rhs, observed(rn2)[1].rhs)
    @test isequal(observed(rn1)[1].lhs.metadata, observed(rn2)[1].lhs.metadata)
    @test isequal(unknowns(rn1), unknowns(rn2))

        @test f_eval(reaction_networks_constraint[1], u, p_1, τ) ≈ f_eval(time_network, u, p_2, τ)
        @test jac_eval(reaction_networks_constraint[1], u, p_1, τ) ≈ jac_eval(time_network, u, p_2, τ)
        @test g_eval(reaction_networks_constraint[1], u, p_1, τ) ≈ g_eval(time_network, u, p_2, τ)
    end
end

# Check that various symbols can be used as species/parameter names.
let
    @reaction_network begin
        (a, A), n ⟷ N
        (b, B), o ⟷ O
        (c, C), p ⟷ P
        (d, D), q ⟷ Q
        (e, E), r ⟷ R
        (f, F), s ⟷ S
        (g, G), u ⟷ U
        (h, H), v ⟷ V
        (j, J), w ⟷ W
        (k, K), x ⟷ X
        (l, L), y ⟷ Y
        (m, M), z ⟷ Z
    end
    @test isequal(observed(rn1)[1].lhs, X)
    @test ModelingToolkit.getdescription(rn1.X) == "An observable"
    @test isspecies(rn1.X)
    @test length(unknowns(rn1)) == 2

    @reaction_network begin (1.0, 1.0), i ⟷ T end

    @reaction_network begin
        (å, Å), ü ⟷ Ü
        (ä, Ä), ñ ⟷ Ñ
        (ö, Ö), æ ⟷ Æ
    end

    @reaction_network begin
        (α, Α), ν ⟷ Ν
        (β, Β), ξ ⟷ Ξ
        (γ, γ), ο ⟷ Ο
        (δ, Δ), Π ⟷ Π
        (ϵ, Ε), ρ ⟷ Ρ
        (ζ, Ζ), σ ⟷ Σ
        (η, Η), τ ⟷ Τ
        (θ, Θ), υ ⟷ Υ
        (ι, Ι), ϕ ⟷ Φ
        (κ, κ), χ ⟷ Χ
        (λ, Λ), ψ ↔ Ψ
        (μ, Μ), ω ⟷ Ω
    end
end

# Tests specific declaration of observables as species/variables.
let
    rn = @reaction_network begin
        k1, S + I --> 2I
        k2, I --> R
    end

    @test isspecies(rn.X)
    @test !isspecies(rn.Y)
    @test !isspecies(rn.Z)
end

# Tests various erroneous declarations throw errors.
let
    # Independent variable in @observables.
    @test_throws Exception @eval @reaction_network begin
        @observables X(t) ~ X1 + X2
        k, 0 --> X1 + X2
    end

    rn2 = @reaction_network arrowtest begin
        a1, C --> 0
        a2, 0 --> C
        k1, A + B --> C
        k2, C --> A + B
        b1, B --> 0
    end

    @test rn1 == rn2
end

# Tests arrow variants in "@reaction" macro .
let
    @test isequal((@reaction k, 0 --> X), (@reaction k, X <-- 0))
    @test isequal((@reaction k, 0 --> X), (@reaction k, X ⟻ 0))
    @test isequal((@reaction k, 0 --> X), (@reaction k, 0 → X))
    @test isequal((@reaction k, 0 --> X), (@reaction k, 0 ⥟ X))
end

# Test that symbols with special mean, or that are forbidden, are handled properly.
let
    test_network = @reaction_network begin t * k, X --> ∅ end
    @test length(species(test_network)) == 1
    @test length(parameters(test_network)) == 1

    test_network = @reaction_network begin π, X --> ∅ end
    @test length(species(test_network)) == 1
    @test length(parameters(test_network)) == 0
    @test reactions(test_network)[1].rate == π

    test_network = @reaction_network begin pi, X --> ∅ end
    @test length(species(test_network)) == 1
    @test length(parameters(test_network)) == 0
    @test reactions(test_network)[1].rate == pi

    test_network = @reaction_network begin ℯ, X --> ∅ end
    @test length(species(test_network)) == 1
    @test length(parameters(test_network)) == 0
    @test reactions(test_network)[1].rate == ℯ


### Test `@equations` Option for Coupled CRN/Equations Models ###

# Checks creation of basic network.
# Check indexing of output solution.
# Check that DAE is solved correctly.
let
    rn = @reaction_network rn begin
        @parameters k
        @variables X(t) Y(t)
        @equations begin
            X + 5 ~ k*S
            3Y + X  ~ S + X*d
        end
        (p,d), 0 <--> S
    end

    @unpack X, Y, S, k, p, d = rn

    # Checks that the internal structures have the correct lengths.
    @test length(species(rn)) == 1
    @test length(unknowns(rn)) == 3
    @test length(reactions(rn)) == 2
    @test length(equations(rn)) == 4
    @test !has_diff_equations(rn)
    @test isequal(diff_equations(rn), [])
    @test has_alg_equations(rn)
    @test isequal(alg_equations(rn), [X + 5 ~ k*S, 3Y + X  ~ S + X*d])

    # Checks that the internal structures contain the correct stuff, and are correctly sorted.
    @test isspecies(unknowns(rn)[1])
    @test !isspecies(unknowns(rn)[2])
    @test !isspecies(unknowns(rn)[3])
    @test equations(rn)[1] isa Reaction
    @test equations(rn)[2] isa Reaction
    @test equations(rn)[3] isa Equation
    @test equations(rn)[3] isa Equation
    @test isequal(equations(rn)[3], X + 5 ~ k*S)
    @test isequal(equations(rn)[4], 3Y + X  ~ S + X*d)

    # Checks that simulations has the correct output
    u0 = Dict([S => 1 + rand(rng), X => 1 + rand(rng), Y => 1 + rand(rng)])
    ps = Dict([p => 1 + rand(rng), d => 1 + rand(rng), k => 1 + rand(rng)])
    oprob = ODEProblem(rn, u0, (0.0, 10000.0), ps; structural_simplify=true)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[S][end] ≈ sol.ps[p]/sol.ps[d]
    @test sol[X] .+ 5 ≈ sol.ps[k] .*sol[S]
    @test 3*sol[Y] .+ sol[X] ≈ sol[S] .+ sol[X].*sol.ps[d]
end

# Checks that block form is not required when only a single equation is used.
let
    rn1 = @reaction_network rn begin
        @parameters k
        @variables X(t)
        @equations X + 2 ~ k*S
        (p,d), 0 <--> S
    end
    rn2 = @reaction_network rn begin
        @parameters k
        @variables X(t)
        @equations begin
            X + 2 ~ k*S
        end
        (p,d), 0 <--> S
    end
    @test rn1 == rn2
end

# Tries for reaction system without any reactions (just an equation).
# Tries with interpolating a value into an equation.
# Tries using rn.X notation for designating variables.
# Tries for empty parameter vector.
let
    c = 6.0
    rn = complete(@reaction_network begin
        @variables X(t)
        @equations 2X ~ $c - X
    end)

    u0 = [rn.X => 0.0]
    ps = []
    oprob = ODEProblem(rn, u0, (0.0, 100.0), ps; structural_simplify=true)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[rn.X][end] ≈ 2.0
end

# Checks hierarchical model.
let
    base_rn = @network_component begin
        @variables V1(t)
        @equations begin
            X*3V1 ~ X - 2
        end
        (p,d), 0 <--> X
    end
    @unpack X, V1, p, d = base_rn

    internal_rn = @network_component begin
        @variables V2(t)
        @equations begin
            X*4V2 ~ X - 3
        end
        (p,d), 0 <--> X
    end

    rn = complete(compose(base_rn, [internal_rn]))

    u0 = [V1 => 1.0, X => 3.0, internal_rn.V2 => 2.0, internal_rn.X => 4.0]
    ps = [p => 1.0, d => 0.2, internal_rn.p => 2.0, internal_rn.d => 0.5]
    oprob = ODEProblem(rn, u0, (0.0, 1000.0), ps; structural_simplify=true)
    sol = solve(oprob, Rosenbrock23(); abstol=1e-9, reltol=1e-9)

    @test sol[X][end] ≈ 5.0
    @test sol[X][end]*3*sol[V1][end] ≈ sol[X][end] - 2
    @test sol[internal_rn.X][end] ≈ 4.0
end

# Check for combined differential and algebraic equation.
# Check indexing of output solution using Symbols.
let
    rn = @reaction_network rn begin
        @parameters k
        @variables X(t) Y(t)
        @equations begin
            X + 5 ~ k*S
            D(Y) ~ X + S - 5*Y
        end
        (p,d), 0 <--> S
    end
    @unpack X, Y, S, p, d, k = rn

    # Checks that the internal structures have the correct lengths.
    @test length(species(rn)) == 1
    @test length(unknowns(rn)) == 3
    @test length(reactions(rn)) == 2
    @test length(equations(rn)) == 4
    @test has_diff_equations(rn)
    @test length(diff_equations(rn)) == 1
    @test has_alg_equations(rn)
    @test length(alg_equations(rn)) == 1

    # Checks that the internal structures contain the correct stuff, and are correctly sorted.
    @test isspecies(unknowns(rn)[1])
    @test !isspecies(unknowns(rn)[2])
    @test !isspecies(unknowns(rn)[3])
    @test equations(rn)[1] isa Reaction
    @test equations(rn)[2] isa Reaction
    @test equations(rn)[3] isa Equation
    @test equations(rn)[3] isa Equation

    # Checks that simulations has the correct output
    u0 = Dict([S => 1 + rand(rng), X => 1 + rand(rng), Y => 1 + rand(rng)])
    ps = Dict([p => 1 + rand(rng), d => 1 + rand(rng), k => 1 + rand(rng)])
    oprob = ODEProblem(rn, u0, (0.0, 10000.0), ps; structural_simplify=true)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[:S][end] ≈ sol.ps[:p]/sol.ps[:d]
    @test sol[:X] .+ 5 ≈ sol.ps[:k] .*sol[:S]
    @test 5*sol[:Y][end] ≈ sol[:S][end] + sol[:X][end]
end

# Tests that various erroneous declarations throw errors.
let
    # Using = instead of ~ (for equation).
    @test_throws Exception @eval @reaction_network begin
        @variables X(t)
        @equations X = 1 - S
        (p,d), 0 <--> S
    end

    # Equation with component undeclared elsewhere.
    @test_throws Exception @eval @reaction_network begin
        @equations X ~ p - S
        (P,D), 0 <--> S
    end
end

# test combinatoric_ratelaws DSL option
let
    rn = @reaction_network begin
        @combinatoric_ratelaws false
        (k1,k2), 2A <--> B
        end
    combinatoric_ratelaw = Catalyst.get_combinatoric_ratelaws(rn)
    @test combinatoric_ratelaw == false
    rl = oderatelaw(reactions(rn)[1]; combinatoric_ratelaw)
    @unpack k1, A = rn
    @test isequal(rl, k1*A^2)

    rn2 = @reaction_network begin
        @combinatoric_ratelaws true
        (k1,k2), 2A <--> B
        end
    combinatoric_ratelaw = Catalyst.get_combinatoric_ratelaws(rn2)
    @test combinatoric_ratelaw == true
    rl = oderatelaw(reactions(rn2)[1]; combinatoric_ratelaw)
    @unpack k1, A = rn2
    @test isequal(rl, k1*A^2/2)

    crl = false
    rn3 = @reaction_network begin
        @combinatoric_ratelaws $crl
        (k1,k2), 2A <--> B
        end
    combinatoric_ratelaw = Catalyst.get_combinatoric_ratelaws(rn3)
    @test combinatoric_ratelaw == crl
    rl = oderatelaw(reactions(rn3)[1]; combinatoric_ratelaw)
    @unpack k1, A = rn3
    @test isequal(rl, k1*A^2)
end