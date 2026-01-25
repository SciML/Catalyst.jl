#! format: off

### Prepares Tests ###

# Fetch packages.
using Catalyst, NonlinearSolve, OrdinaryDiffEqTsit5, SparseArrays, StochasticDiffEq, Test
using LinearAlgebra: norm
using ModelingToolkitBase: value

# Sets the default `t` to use.
t = default_t()

# Fetch test networks.
include("../test_networks.jl")

### Tests Basic Getters ###

# Checks various getter functions.
# Uses several system-modifying functions, and should probably be rewritten not to use these.
let
    @parameters k1 k2
    @species S(t) I(t) R(t)
    rxs = [Reaction(k1, [S, I], [I], [1, 1], [2]),
        Reaction(k2, [I], [R])]
    @named rs = ReactionSystem(rxs, t, [S, I, R], [k1, k2])

    specset = Set([value(S) => 1, value(I) => 2, value(R) => 3])
    @test issetequal(specset, speciesmap(rs))

    pset = Set([value(k1) => 1, value(k2) => 2])
    @test issetequal(pset, paramsmap(rs))

    rxs = [Reaction(k1, [S, I], [I], [1, 1], [2]),
        Reaction(k2, [I], [R])]
    @named rs = ReactionSystem(rxs, t, [S, I, R], [k1, k2])
    rxs2 = [Reaction(k2, [I], [R], [1], [1]),
        Reaction(k1, [S, I], [I], [1, 1], [2])]
    rs2 = ReactionSystem(rxs2, t, [R, I, S], [k2, k1], name = :rs)
    @test rs == rs2

    @parameters k3 k4
    @species D(t)
    rxs3 = [Reaction(k3, [S], [D]), Reaction(k4, [S, I], [D])]
    @named rs3 = ReactionSystem(rxs3, t)
    @test issetequal(species(rs3), [S, D, I])
    @test issetequal(parameters(rs3), [k3, k4])
    @named rs = extend(rs3, rs)
    rxs2b = [Reaction(k3, [S], [D]), Reaction(k4, [S, I], [D])]
    @named rs2b = ReactionSystem(rxs2b, t)
    @named rs2 = extend(rs2b, rs2; name = :rs)
    @test rs2 == rs

    rxs = [Reaction(k1, [S, I], [I], [1, 1], [2]),
        Reaction(k2, [I], [R])]
    @named rs = ReactionSystem(rxs, t)
    rxs3 = [Reaction(k3, [S], [D]), Reaction(k4, [S, I], [D])]
    @named rs3 = ReactionSystem(rxs3, t)
    rs4 = extend(rs, rs3; name = :rs)
    @test rs2 == rs4

    rxs = [Reaction(k1 * S, [S, I], [I], [2, 3], [2]),
        Reaction(k2 * R, [I], [R])]
    @named rs = ReactionSystem(rxs, t)
    deps = dependents(rxs[2], rs)
    @test issetequal(deps, [R, I])
    @test issetequal(dependents(rxs[1], rs), dependants(rxs[1], rs))
 end

# Tests `substoichmat` and `prodstoichmat` getters.
let
    rnmat = @reaction_network begin
        α, S + 2I --> 2I
        β, 3I --> 2R + S
    end

    smat = [1 0;
            2 3;
            0 0]
    pmat = [0 1;
            2 0;
            0 2]
    @test smat == substoichmat(rnmat) == Matrix(substoichmat(rnmat, sparse = true))
    @test pmat == prodstoichmat(rnmat) == Matrix(prodstoichmat(rnmat, sparse = true))
end

# Tests `reactionrates`, and `symmap_to_varmap` getters.
let
    rn = @reaction_network begin
        (p,d), 0 <--> X
        (kB,kD), 2X <--> X
    end
    @unpack p, d, kB, kD = rn
    issetequal(reactionrates(rn), [p, d, kB, kD])
    isequal(symmap_to_varmap(rn, [:p => 1.0, :kB => 3.0]), [p => 1.0, kB => 3.0])
end

### Test Intermediate Complexes Reaction Networks ###

# Test function.
function testnetwork(rn, B, Z, Δ, lcs, d, subrn, lcd; skiprxtest = false)
    B2 = reactioncomplexes(rn)[2]
    @test B == B2 == Matrix(reactioncomplexes(rn, sparse = true)[2])
    @test B == incidencemat(rn)
    @test Z == complexstoichmat(rn) == Matrix(complexstoichmat(rn, sparse = true))
    @test Δ == complexoutgoingmat(rn) == Matrix(complexoutgoingmat(rn, sparse = true))
    ig = incidencematgraph(rn)
    lcs2 = linkageclasses(rn)
    @test lcs2 == linkageclasses(incidencematgraph(sparse(B))) == lcs
    @test deficiency(rn) == d
    if !skiprxtest
        @test all(issetequal.(subrn, reactions.(subnetworks(rn))))
    end
    @test linkagedeficiencies(rn) == lcd
    @test sum(linkagedeficiencies(rn)) <= deficiency(rn)
end

# Mass-action non-catalytic.
let
    rn = @reaction_network begin
        k₁, 2A --> B
        k₂, A --> C
        k₃, C --> D
        k₄, B + D --> E
    end
    Z = [2 0 1 0 0 0 0;
        0 1 0 0 0 1 0;
        0 0 0 1 0 0 0;
        0 0 0 0 1 1 0;
        0 0 0 0 0 0 1]
    B = [-1 0 0 0;
        1 0 0 0;
        0 -1 0 0;
        0 1 -1 0;
        0 0 1 0;
        0 0 0 -1;
        0 0 0 1]
    Δ = [-1 0 0 0; 0 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 0; 0 0 0 -1; 0 0 0 0]
    lcs = [[1, 2], [3, 4, 5], [6, 7]]
    r = reactions(rn)
    subrn = [[r[1]], [r[2], r[3]], [r[4]]]
    lcd = [0, 0, 0]
    testnetwork(rn, B, Z, Δ, lcs, 0, subrn, lcd)

    # Constant and BC species test.
    @parameters F [isconstantspecies = true]
    crn = @reaction_network begin
        k₁, 2A --> B
        k₂, A --> C + $F
        k₃, C --> D
        k₄, B + D + $F --> E
    end
    testnetwork(crn, B, Z, Δ, lcs, 0, subrn, lcd; skiprxtest = true)
end

# Mass-action rober.
let
    rn = @reaction_network begin
        k₁, A --> B
        k₂, B + B --> C + B
        k₃, B + C --> A + C
    end
    Z = [1 0 0 0 1;
        0 1 2 1 0;
        0 0 0 1 1]
    B = [-1 0 0;
        1 0 0;
        0 -1 0;
        0 1 -1;
        0 0 1]
    Δ = [-1 0 0; 0 0 0; 0 -1 0; 0 0 -1; 0 0 0]
    lcs = [[1, 2], [3, 4, 5]]
    r = reactions(rn)
    subrn = [[r[1]], [r[2], r[3]]]
    lcd = [0, 0]
    testnetwork(rn, B, Z, Δ, lcs, 1, subrn, lcd)
end

# Some rational functions as rates.
let
    rn = @reaction_network begin
        k₁, ∅ --> X₁
        (k₂ / (1 + X₁ * X₂ + X₃ * X₄), k₃ / (1 + X₁ * X₂ + X₃ * X₄)), 2X₁ + X₂ ↔ 3X₃ + X₄
        k₄, X₄ --> ∅
    end
    Z = [0 1 2 0 0;
        0 0 1 0 0;
        0 0 0 3 0;
        0 0 0 1 1]
    B = [-1 0 0 1;
        1 0 0 0;
        0 -1 1 0;
        0 1 -1 0;
        0 0 0 -1]
    Δ = [-1 0 0 0; 0 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]
    lcs = [[1, 2, 5], [3, 4]]
    r = reactions(rn)
    subrn = [[r[1], r[4]], [r[2], r[3]]]
    lcd = [0, 0]
    testnetwork(rn, B, Z, Δ, lcs, 0, subrn, lcd)
end

# Repressilator.
let
    rn = @reaction_network begin
        hillr(P₃, α, K, n), ∅ --> m₁
        hillr(P₁, α, K, n), ∅ --> m₂
        hillr(P₂, α, K, n), ∅ --> m₃
        (δ, γ), m₁ ↔ ∅
        (δ, γ), m₂ ↔ ∅
        (δ, γ), m₃ ↔ ∅
        β, m₁ --> m₁ + P₁
        β, m₂ --> m₂ + P₂
        β, m₃ --> m₃ + P₃
        μ, P₁ --> ∅
        μ, P₂ --> ∅
        μ, P₃ --> ∅
    end
    Z = [0 1 0 0 1 0 0 0 0 0;
        0 0 1 0 0 1 0 0 0 0;
        0 0 0 1 0 0 1 0 0 0;
        0 0 0 0 1 0 0 1 0 0;
        0 0 0 0 0 1 0 0 1 0;
        0 0 0 0 0 0 1 0 0 1]
    B = [-1 -1 -1 1 -1 1 -1 1 -1 0 0 0 1 1 1;
        1 0 0 -1 1 0 0 0 0 -1 0 0 0 0 0;
        0 1 0 0 0 -1 1 0 0 0 -1 0 0 0 0;
        0 0 1 0 0 0 0 -1 1 0 0 -1 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]
    Δ = [-1 -1 -1 0 -1 0 -1 0 -1 0 0 0 0 0 0; 0 0 0 -1 0 0 0 0 0 -1 0 0 0 0 0;
        0 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0;
        0 0 0 0 0 0 0 -1 0 0 0 -1 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]
    lcs = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
    r = reactions(rn)
    subrn = [[r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[13], r[14], r[15],
        r[10], r[11], r[12]]]
    lcd = [3]
    testnetwork(rn, B, Z, Δ, lcs, 3, subrn, lcd)
end

# Brusselator.
let
    rn = @reaction_network begin
        A, ∅ → X
        1, 2X + Y → 3X
        B, X → Y
        1, X → ∅
    end
    Z = [0 1 2 3 0;
        0 0 1 0 1]
    B = [-1 0 0 1;
        1 0 -1 -1;
        0 -1 0 0;
        0 1 0 0;
        0 0 1 0]
    Δ = [-1 0 0 0; 0 0 -1 -1; 0 -1 0 0; 0 0 0 0; 0 0 0 0]
    lcs = [[1, 2, 5], [3, 4]]
    r = reactions(rn)
    subrn = [[r[1], r[4], r[3]], [r[2]]]
    lcd = [0, 0]
    testnetwork(rn, B, Z, Δ, lcs, 1, subrn, lcd)
end

# Some rational functions as rates.
let
    rn = @reaction_network begin
        (k₁, k₋₁), X₁ + X₂ <--> X₃ + 2X₄
        (k₂ / (1 + X₄ * X₅ + X₆ * X₇), k₋₂ / (1 + X₄ * X₅ + X₆ * X₇)), 3X₄ + X₅ <--> X₆ + X₇
        (k₃ / (1 + X₇ + X₈ + X₉ + X₁₀), k₋₃ / (1 + X₇ + X₈ + X₉ + X₁₀)), 5X₇ + X₈ <--> X₉ + X₁₀
    end
    Z = [1 0 0 0 0 0;
        1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 2 3 0 0 0;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 1 5 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        0 0 0 0 0 1]
    B = [-1 1 0 0 0 0;
        1 -1 0 0 0 0;
        0 0 -1 1 0 0;
        0 0 1 -1 0 0;
        0 0 0 0 -1 1;
        0 0 0 0 1 -1]
    Δ = [-1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0; 0 0 0 0 0 -1]
    lcs = [[1, 2], [3, 4], [5, 6]]
    r = reactions(rn)
    subrn = [[r[1], r[2]], [r[3], r[4]], [r[5], r[6]]]
    lcd = [0, 0, 0]
    testnetwork(rn, B, Z, Δ, lcs, 0, subrn, lcd)
end

### Other Tests ###

let
    myrn = [reaction_networks_standard; reaction_networks_hill; reaction_networks_real]
    for i in 1:length(myrn)
        local rcs, B = reactioncomplexes(myrn[i])
        @test B == Matrix(reactioncomplexes(myrn[i], sparse = true)[2])
        local Z = complexstoichmat(myrn[i])
        @test Z == Matrix(complexstoichmat(myrn[i], sparse = true))
        @test Z * B == netstoichmat(myrn[i]) == Matrix(netstoichmat(myrn[i], sparse = true))
    end
end

# Test @unpack with symbolic variables.
let
    tspan = (0.0, 250.0)

    function unpacktest(rn)
        @unpack S1, I1, R1, α1, β1 = rn
        u₀ = [S1 => 999.0, I1 => 1.0, R1 => 0.0]
        p = [α1 => 1e-4, β1 => 0.01]
        op = ODEProblem(rn, u₀, (0.0, 250.0), p)
        solve(op, Tsit5())
    end
    rn = @reaction_network begin
        α1, S1 + I1 --> 2I1
        β1, I1 --> R1
    end
    sol = unpacktest(rn)

    # Test symmap_to_varmap.
    sir = @network_component sir begin
        β, S + I --> 2I
        ν, I --> R
    end
    subsys = @network_component subsys begin
        k, A --> B
    end
    @named sys = compose(sir, [subsys])
    sir = complete(sir)
    sys = complete(sys)

    symmap = [:S => 1.0, :I => 1.0, :R => 1.0, :subsys₊A => 1.0, :subsys₊B => 1.0]
    u0map = symmap_to_varmap(sys, symmap)
    pmap = symmap_to_varmap(sys, [:β => 1.0, :ν => 1.0, :subsys₊k => 1.0])
    @test isequal(u0map[4][1], subsys.A)
    @test isequal(u0map[1][1], @nonamespace sir.S)

    u0map = symmap_to_varmap(sir, [:S => 999.0, :I => 1.0, :R => 0.0])
    pmap = symmap_to_varmap(sir, [:β => 1e-4, :ν => 0.01])
    op = ODEProblem(sir, u0map, tspan, pmap)
    sol2 = solve(op, Tsit5())
    @test norm(sol.u - sol2.u) ≈ 0 atol = 1e-10

    u0map = [:S => 999.0, :I => 1.0, :R => 0.0]
    pmap = (:β => 1e-4, :ν => 0.01)
    op = ODEProblem(sir, u0map, tspan, pmap)
    sol3 = solve(op, Tsit5())
    @test norm(sol.u - sol3.u) ≈ 0 atol = 1e-10
end


# Tests non-integer stoichiometry.
let
    function test_stoich(T, rn)
        @test eltype(substoichmat(rn)) == T
        @test eltype(prodstoichmat(rn)) == T
        @test eltype(netstoichmat(rn)) == T
        @test eltype(substoichmat(rn; sparse = true)) == T
        @test eltype(prodstoichmat(rn; sparse = true)) == T
        @test eltype(netstoichmat(rn; sparse = true)) == T
        nothing
    end

    rn = @reaction_network ABtoC begin
        (k₊,k₋), 3.4*A + 2B <--> 2.5*C
    end
    test_stoich(Float64, rn)

    rn2 = @reaction_network ABtoC begin
        (k₊,k₋), 3*A + 2B <--> 2*C
    end
    test_stoich(Int, rn2)
end

### Test Polynomial Transformation Functionality ###

# Tests normal network.
let
    rn = @reaction_network begin
        (p,d), 0 <--> X
        (kB,kD), 2X <--> X2
    end
    ns = make_rre_algeqs(rn)
    neweqs = getfield.(equations(ns), :rhs)
    poly = Catalyst.to_multivariate_poly(neweqs)
    @test length(poly) == 2
end

# Tests network with a fraction.
let
    rn = @reaction_network begin
        (p/X,d), 0 <--> X
    end
    ns = make_rre_algeqs(rn)
    neweqs = getfield.(equations(ns), :rhs)
    poly = Catalyst.to_multivariate_poly(neweqs)
    @test length(poly) == 1
end

# Test empty network.
let
    rn = @reaction_network
    ns = make_rre_algeqs(rn)
    neweqs = getfield.(equations(ns), :rhs)
    @test_throws AssertionError Catalyst.to_multivariate_poly(neweqs)
end

# Tests `isautonomous` function.
let
    # Using default iv.
    rn1 = @reaction_network begin
        (p + X*(p1/(t+p3)),d), 0 <--> X
        (kB,kD), 2X <--> X
    end
    rn2 = @reaction_network begin
        (hill(X, v/t, K, n),d), 0 <--> X
        (kB,kD), 2X <--> X
    end
    rn3 = @reaction_network begin
        (p + X*(p1+p2),d), 0 <--> X
        (kB,kD), 2X <--> X
    end
    @test !isautonomous(rn1)
    @test !isautonomous(rn2)
    @test isautonomous(rn3)

    # Using non-default iv.
    rn4 = @reaction_network begin
        @ivs i1 i2
        (p + X*(p1/(1+i1)),d), 0 <--> X
        (kB,kD), 2X <--> X
    end
    rn5 = @reaction_network begin
        @ivs i1 i2
        (p + X*(i2+p2),d), 0 <--> X
        (kB,kD), 2X <--> X
    end
    rn6 = @reaction_network begin
        @ivs i1 i2
        (hill(X, v/i1, i2, n),d), 0 <--> X
        (kB,kD), 2X <--> X
    end
    rn7 = @reaction_network begin
        @ivs i1 i2
        (p + X*(p1+p2),d), 0 <--> X
        (kB,kD), 2X <--> X
    end
    @test !isautonomous(rn4)
    @test !isautonomous(rn5)
    @test !isautonomous(rn6)
    @test isautonomous(rn7)

    # Using a coupled CRN/equation model.
    rn8 = @reaction_network begin
        @equations D(V) ~ X/(1+t) - V
        (p,d), 0 <--> X
    end
    @test !isautonomous(rn8)

    # Using a registered function.
    f(d,t) = d/(1 + t)
    Symbolics.@register_symbolic f(d,t)
    rn9 = @reaction_network begin
        f(d,t), X --> 0
    end
    @test !isautonomous(rn9)
end
