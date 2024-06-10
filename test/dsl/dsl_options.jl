#! format: off

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

# Test creating networks with/without options.
let
    @reaction_network begin (k1, k2), A <--> B end
    @reaction_network begin
        @parameters k1 k2
        (k1, k2), A <--> B
    end
    @reaction_network begin
        @parameters k1 k2
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    @reaction_network begin
        @species A(t) B(t)
        (k1, k2), A <--> B
    end

    @reaction_network begin
        @parameters begin
            k1
            k2
        end
        (k1, k2), A <--> B
    end
    @reaction_network begin
        @species begin
            A(t)
            B(t)
        end
        (k1, k2), A <--> B
    end
    @reaction_network begin
        @parameters begin
            k1
            k2
        end
        @species begin
            A(t)
            B(t)
        end
        (k1, k2), A <--> B
    end

    n1 = @reaction_network rnname begin (k1, k2), A <--> B end
    n2 = @reaction_network rnname begin
        @parameters k1 k2
        (k1, k2), A <--> B
    end
    n3 = @reaction_network rnname begin
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    n4 = @reaction_network rnname begin
        @parameters k1 k2
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    n5 = @reaction_network rnname begin
        (k1, k2), A <--> B
        @parameters k1 k2
    end
    n6 = @reaction_network rnname begin
        (k1, k2), A <--> B
        @species A(t) B(t)
    end
    n7 = @reaction_network rnname begin
        (k1, k2), A <--> B
        @parameters k1 k2
        @species A(t) B(t)
    end
    n8 = @reaction_network rnname begin
        @parameters begin
            k1
            k2
        end
        (k1, k2), A <--> B
    end
    n9 = @reaction_network rnname begin
        @species begin
            A(t)
            B(t)
        end
        (k1, k2), A <--> B
    end
    n10 = @reaction_network rnname begin
        @parameters begin
            k1
            k2
        end
        @species begin
            A(t)
            B(t)
        end
        (k1, k2), A <--> B
    end
    @test all(==(n1), (n2, n3, n4, n5, n6, n7, n8, n9, n10))
end

# Tests that when either @species or @parameters is given, the other is inferred properly.
let
    rn1 = @reaction_network begin
        k*X, A + B --> 0
    end
    @test issetequal(species(rn1), @species A(t) B(t))
    @test issetequal(parameters(rn1), @parameters k X)

    rn2 = @reaction_network begin
        @species A(t) B(t) X(t)
        k*X, A + B --> 0
    end
    @test issetequal(species(rn2), @species A(t) B(t) X(t))
    @test issetequal(parameters(rn2), @parameters k)

    rn3 = @reaction_network begin
        @parameters k
        k*X, A + B --> 0
    end
    @test issetequal(species(rn3), @species A(t) B(t))
    @test issetequal(parameters(rn3), @parameters k X)

    rn4 = @reaction_network begin
        @species A(t) B(t) X(t)
        @parameters k
        k*X, A + B --> 0
    end
    @test issetequal(species(rn4), @species A(t) B(t) X(t))
    @test issetequal(parameters(rn4), @parameters k)

    rn5 = @reaction_network begin
        @parameters k B [isconstantspecies=true]
        k*X, A + B --> 0
    end
    @test issetequal(species(rn5), @species A(t))
    @test issetequal(parameters(rn5), @parameters k B X)
end

# Test inferring with stoichiometry symbols and interpolation.
let
    @parameters k g h gg X y [isconstantspecies = true]
    t = Catalyst.DEFAULT_IV
    @species A(t) B(t) BB(t) C(t)

    rni = @reaction_network inferred begin
        $k*X, $y + g*A + h*($gg)*B + $BB * C --> k*C
    end
    @test issetequal(species(rni), [A, B, BB, C])
    @test issetequal(parameters(rni), [k, g, h, gg, X, y])

    rnii = @reaction_network inferred begin
        @species BB(t)
        @parameters y [isconstantspecies = true]
        k*X, y + g*A + h*($gg)*B + BB * C --> k*C
    end
    @test rnii == rni
end

# Tests that when some species or parameters are left out, the others are set properly.
let
    rn6 = @reaction_network begin
        @species A(t)
        k*X, A + B --> 0
    end
    @test issetequal(species(rn6), @species A(t) B(t))
    @test issetequal(parameters(rn6), @parameters k X)

    rn7 = @reaction_network begin
        @species A(t) X(t)
        k*X, A + B --> 0
    end
    @test issetequal(species(rn7), @species A(t) X(t) B(t))
    @test issetequal(parameters(rn7), @parameters k)

    rn7 = @reaction_network begin
        @parameters B [isconstantspecies=true]
        k*X, A + B --> 0
    end
    @test issetequal(species(rn7), @species A(t))
    @test issetequal(parameters(rn7), @parameters B k X)

    rn8 = @reaction_network begin
        @parameters B [isconstantspecies=true] k
        k*X, A + B --> 0
    end
    @test issetequal(species(rn8), @species A(t))
    @test issetequal(parameters(rn8), @parameters B k X)

    rn9 = @reaction_network begin
        @parameters k1 X1
        @species A1(t) B1(t)
        k1*X1, A1 + B1 --> 0
        k2*X2, A2 + B2 --> 0
    end
    @test issetequal(species(rn9), @species A1(t) B1(t) A2(t) B2(t))
    @test issetequal(parameters(rn9), @parameters k1 X1 k2 X2)

    rn10 = @reaction_network begin
        @parameters k1 X2 B2 [isconstantspecies=true]
        @species A1(t) X1(t)
        k1*X1, A1 + B1 --> 0
        k2*X2, A2 + B2 --> 0
    end
    @test issetequal(species(rn10), @species A1(t) X1(t) B1(t) A2(t))
    @test issetequal(parameters(rn10), @parameters k1 X2 B2 k2)

    rn11 = @reaction_network begin
        @parameters k1 k2
        @species X1(t)
        k1*X1, A1 + B1 --> 0
        k2*X2, A2 + B2 --> 0
    end
    @test issetequal(species(rn11), @species X1(t) A1(t) A2(t) B1(t) B2(t))
    @test issetequal(parameters(rn11), @parameters k1 k2 X2)
end

# Checks that some created networks are identical.
let
    rn12 = @reaction_network rnname begin (k1, k2), A <--> B end
    rn13 = @reaction_network rnname begin
        @parameters k1 k2
        (k1, k2), A <--> B
    end
    rn14 = @reaction_network rnname begin
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    rn15 = @reaction_network rnname begin
        @parameters k1 k2
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    @test all(==(rn12), (rn13, rn14, rn15))
end

# Checks that the rights things are put in vectors.
let
    rn18 = @reaction_network rnname begin
        @parameters p d1 d2
        @species A(t) B(t)
        p, 0 --> A
        1, A --> B
        (d1, d2), (A, B) --> 0
    end
    rn19 = @reaction_network rnname begin
        p, 0 --> A
        1, A --> B
        (d1, d2), (A, B) --> 0
    end
    @test rn18 == rn19

    @parameters p d1 d2
    @species A(t) B(t)
    @test isequal(parameters(rn18)[1], p)
    @test isequal(parameters(rn18)[2], d1)
    @test isequal(parameters(rn18)[3], d2)
    @test isequal(species(rn18)[1], A)
    @test isequal(species(rn18)[2], B)

    rn20 = @reaction_network rnname begin
        @species X(t)
        @parameters S
        mm(X,v,K), 0 --> Y
        (k1,k2), 2Y <--> Y2
        d*Y, S*(Y2+Y) --> 0
    end
    rn21 = @reaction_network rnname begin
        @species X(t) Y(t) Y2(t)
        @parameters v K k1 k2 d S
        mm(X,v,K), 0 --> Y
        (k1,k2), 2Y <--> Y2
        d*Y, S*(Y2+Y) --> 0
    end
    rn22 = @reaction_network rnname begin
        @species X(t) Y2(t)
        @parameters d k1
        mm(X,v,K), 0 --> Y
        (k1,k2), 2Y <--> Y2
        d*Y, S*(Y2+Y) --> 0
    end
    @test all(==(rn20), (rn21, rn22))
    @parameters v K k1 k2 d S
    @species X(t) Y(t) Y2(t)
    @test issetequal(parameters(rn22),[v K k1 k2 d S])
    @test issetequal(species(rn22), [X Y Y2])
end

# Tests that defaults work.
let
    rn26 = @reaction_network rnname begin
        @parameters p=1.0 d1 d2=5
        @species A(t) B(t)=4
        p, 0 --> A
        1, A --> B
        (d1, d2), (A, B) --> 0
    end

    rn27 = @reaction_network rnname begin
    @parameters p1=1.0 p2=2.0 k1=4.0 k2=5.0 v=8.0 K=9.0 n=3 d=10.0
    @species X(t)=4.0 Y(t)=3.0 X2Y(t)=2.0 Z(t)=1.0
        (p1,p2), 0 --> (X,Y)
        (k1,k2), 2X + Y --> X2Y
        hill(X2Y,v,K,n), 0 --> Z
        d, (X,Y,X2Y,Z) --> 0
    end
    u0_27 = []
    p_27 = []

    rn28 = @reaction_network rnname begin
    @parameters p1=1.0 p2 k1=4.0 k2 v=8.0 K n=3 d
    @species X(t)=4.0 Y(t) X2Y(t) Z(t)=1.0
        (p1,p2), 0 --> (X,Y)
        (k1,k2), 2X + Y --> X2Y
        hill(X2Y,v,K,n), 0 --> Z
        d, (X,Y,X2Y,Z) --> 0
    end
    u0_28 = symmap_to_varmap(rn28, [:p2=>2.0, :k2=>5.0, :K=>9.0, :d=>10.0])
    p_28 = symmap_to_varmap(rn28, [:Y=>3.0, :X2Y=>2.0])
    defs28 = Dict(Iterators.flatten((u0_28, p_28)))

    rn29 = @reaction_network rnname begin
    @parameters p1 p2 k1 k2 v K n d
    @species X(t) Y(t) X2Y(t) Z(t)
        (p1,p2), 0 --> (X,Y)
        (k1,k2), 2X + Y --> X2Y
        hill(X2Y,v,K,n), 0 --> Z
        d, (X,Y,X2Y,Z) --> 0
    end
    u0_29 = symmap_to_varmap(rn29, [:p1=>1.0, :p2=>2.0, :k1=>4.0, :k2=>5.0, :v=>8.0, :K=>9.0, :n=>3, :d=>10.0])
    p_29 = symmap_to_varmap(rn29, [:X=>4.0, :Y=>3.0, :X2Y=>2.0, :Z=>1.0])
    defs29 = Dict(Iterators.flatten((u0_29, p_29)))

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
        (p,d), 0 <--> Xi
        (k1,k2), Xi <--> Xa
        (k3,k4), Y1 <--> Y2
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

# Compares programmatic and DSL system with observables.
let
    # Model declarations.
    rn_dsl = @reaction_network begin
        @observables begin
            X ~ x + 2x2y
            Y ~ y + x2y
        end
        k, 0 --> (x, y)
        (kB, kD), 2x + y <--> x2y
        d, (x,y,x2y) --> 0
    end

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

    # Make simulations.
    u0 = [x => 1.0, y => 0.5, x2y => 0.0]
    tspan = (0.0, 15.0)
    ps = [k => 1.0, kD => 0.1, kB => 0.5, d => 5.0]
    oprob_dsl = ODEProblem(rn_dsl, u0, tspan, ps)
    oprob_prog = ODEProblem(rn_prog, u0, tspan, ps)

    sol_dsl = solve(oprob_dsl, Tsit5(); saveat=0.1)
    sol_prog = solve(oprob_prog, Tsit5(); saveat=0.1)

    # Tests observables equal in both cases.
    @test oprob_dsl[:X] == oprob_prog[:X]
    @test oprob_dsl[:Y] == oprob_prog[:Y]
    @test sol_dsl[:X] == sol_prog[:X]
    @test sol_dsl[:Y] == sol_prog[:Y]
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

    oprob = ODEProblem(rn, u0, (0.0, 1000.0), ps)
    sol = solve(oprob, Tsit5())

    @test sol[:X][1] == u0[:X1]^2 + ps[:op_1]*(u0[:X2] + 2*u0[:X3]) + u0[:X1]*u0[:X4]/ps[:op_2] + ps[:p]
end

# Checks that ivs are correctly found.
let
    rn = @reaction_network begin
        @ivs t x y
        @species V1(t) V2(t,x) V3(t, y) W1(t) W2(t, y)
        @observables begin
            V ~ V1 + 2V2 + 3V3
            W ~ W1 + W2
        end
    end
    V,W = getfield.(observed(rn), :lhs)
    @test isequal(arguments(ModelingToolkit.unwrap(V)), Any[Catalyst.get_iv(rn), Catalyst.get_sivs(rn)[1], Catalyst.get_sivs(rn)[2]])
    @test isequal(arguments(ModelingToolkit.unwrap(W)), Any[Catalyst.get_iv(rn), Catalyst.get_sivs(rn)[2]])
end

# Checks that metadata is written properly.
let
    rn = @reaction_network rn_observed begin
        @observables (X, [description="my_description"]) ~ X1 + X2
        k, 0 --> X1 + X2
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

    # Case with metadata.
    rn3 = @reaction_network rn_observed begin
        @observables (X,  [description="description"]) ~ X1 + X2
        k, 0 --> X1 + X2
    end
    rn4 = @reaction_network rn_observed begin
        @variables X(t) [description="description"]
        @observables X ~ X1 + X2
        k, 0 --> X1 + X2
    end
    @test isequal(observed(rn3)[1].rhs, observed(rn4)[1].rhs)
    @test isequal(observed(rn3)[1].lhs.metadata, observed(rn4)[1].lhs.metadata)
    @test isequal(unknowns(rn3), unknowns(rn4))
end

# Tests for interpolation into the observables option.
let
    # Interpolation into lhs.
    @species X [description="An observable"]
    rn1 = @reaction_network begin
        @observables $X ~ X1 + X2
        (k1, k2), X1 <--> X2
    end
    @test isequal(observed(rn1)[1].lhs, X)
    @test ModelingToolkit.getdescription(rn1.X) == "An observable"
    @test isspecies(rn1.X)
    @test length(unknowns(rn1)) == 2

    # Interpolation into rhs.
    @parameters n [description="A parameter"]
    @species S(t)
    rn2 = @reaction_network begin
        @observables Stot ~ $S + $n*Sn
        (kB, kD), $n*S <--> Sn
    end
    @unpack Stot, Sn, kD, kB = rn2

    u0 = Dict([S => 5.0, Sn => 1.0])
    ps = Dict([n => 2, kB => 1.0, kD => 1.0])
    oprob = ODEProblem(rn2, u0, (0.0, 1.0), ps)

    @test issetequal(Symbolics.get_variables(observed(rn2)[1].rhs), [S, n, Sn])
    @test oprob[Stot] == u0[S] + ps[n]*u0[Sn]
    @test length(unknowns(rn2)) == 2
end

# Tests specific declaration of observables as species/variables.
let
    rn = @reaction_network begin
        @species X(t)
        @variables Y(t)
        @observables begin
            X ~ X + 2X2
            Y ~ Y1 + Y2
            Z ~ X + Y
        end
        (kB,kD), 2X <--> X2
        (k1,k2), Y1 <--> Y2
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

    # System with observable in observable formula.
    @test_throws Exception @eval @reaction_network begin
        @observables begin
            X ~ X1 + X2
            X2 ~ 2X
        end
        (p,d), 0 <--> X1 + X2
    end

    # Multiple @observables options
    @test_throws Exception @eval @reaction_network begin
        @observables X ~ X1 + X2
        @observables Y ~ Y1 + Y2
        k, 0 --> X1 + X2
        k, 0 --> Y1 + Y2
    end
    @test_throws Exception @eval @reaction_network begin
        @observables begin
            X ~ X1 + X2
        end
        @observables begin
            X ~ 2(X1 + X2)
        end
        (p,d), 0 <--> X1 + X2
    end

    # Default value for compound.
    @test_throws Exception @eval @reaction_network begin
        @observables (X = 1.0) ~ X1 + X2
        k, 0 --> X1 + X2
    end

    # Forbidden symbols as observable names.
    @test_throws Exception @eval @reaction_network begin
        @observables t ~ t1 + t2
        k, 0 --> t1 + t2
    end
    @test_throws Exception @eval @reaction_network begin
        @observables im ~ i + m
        k, 0 --> i + m
    end

    # Non-trivial observables expression.
    @test_throws Exception @eval @reaction_network begin
        @observables X - X1 ~ X2
        k, 0 --> X1 + X2
    end

    # Occurrence of undeclared dependants.
    @test_throws Exception @eval @reaction_network begin
        @observables X ~ X1 + X2
        k, 0 --> X1
    end

    # Interpolation and explicit declaration of an observable.
    @variables X(t)
    @test_throws Exception @eval @reaction_network begin
        @variables X(t)
        @observables $X ~ X1 + X2
        (k1,k2), X1 <--> X2
    end
end


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