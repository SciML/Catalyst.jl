### Prepares Tests ###
# Tests for network structural information: associated matrices, graphs, linkage classes, etc.

# Fetch packages.
using Catalyst, LinearAlgebra, Test, SparseArrays
using SymbolicUtils: _iszero

# Sets stable rng number.
using StableRNGs
rng = StableRNG(514)

### Basic Tests ###

# Tests basic `ReactionComplex` properties.
let
    rcs1 = Catalyst.ReactionComplex([1, 2], [1, 3])
    rcs2 = Catalyst.ReactionComplex([1, 2], [1, 3])
    rcs3 = Catalyst.ReactionComplex([3], [2])
    @test rcs1 == rcs2
    @test rcs2 != rcs3
    @test length(rcs1) == 2
    @test length(rcs3) == 1
end

# Tests network analysis functions on MAPK network (by comparing to manually computed outputs).
let
    MAPK = @reaction_network MAPK begin
        (k₁, k₂),KKK + E1 <--> KKKE1
        k₃, KKKE1 --> KKK_ + E1
        (k₄, k₅), KKK_ + E2 <--> KKKE2
        k₆, KKKE2 --> KKK + E2
        (k₇, k₈), KK + KKK_ <--> KK_KKK_
        k₉, KK_KKK_ --> KKP + KKK_
        (k₁₀, k₁₁), KKP + KKK_ <--> KKPKKK_
        k₁₂, KKPKKK_ --> KKPP + KKK_
        (k₁₃, k₁₄), KKP + KKPase <--> KKPKKPase
        k₁₅, KKPPKKPase --> KKP + KKPase
        k₁₆,KKPKKPase --> KK + KKPase
        (k₁₇, k₁₈), KKPP + KKPase <--> KKPPKKPase
        (k₁₉, k₂₀), KKPP + K <--> KKPPK
        k₂₁, KKPPK --> KKPP + KP
        (k₂₂, k₂₃), KKPP + KP <--> KPKKPP
        k₂₄, KPKKPP --> KPP + KKPP
        (k₂₅, k₂₆), KP + KPase <--> KPKPase
        k₂₇, KKPPKPase --> KP + KPase
        k₂₈, KPKPase --> K + KPase
        (k₂₉, k₃₀), KPP + KPase <--> KKPPKPase
    end
    rcs, B = reactioncomplexes(MAPK)
    @test length(rcs) == 26
    num_lcs = length(linkageclasses(MAPK))
    @test num_lcs == 6
    δ = deficiency(MAPK)
    @test δ == 5
    @test all(==(0), linkagedeficiencies(MAPK))
    @test isreversible(MAPK) == false
    @test isweaklyreversible(MAPK, subnetworks(MAPK)) == false
    cls = conservationlaws(MAPK)
    @test Catalyst.get_networkproperties(MAPK).rank == 15

    k = rand(rng, numparams(MAPK))
    rates = Dict(zip(parameters(MAPK), k))
    @test Catalyst.iscomplexbalanced(MAPK, rates) == false
    cyclemat = Catalyst.cycles(MAPK)
    S = netstoichmat(MAPK)
    for i in 1:size(S, 2)-1
        if S[:,i] == -S[:,i+1]
           cycle = [(j == i) || (j == i+1) ? 1 : 0 for j in 1:size(S,2)]
           @test rank(cyclemat) == rank(hcat(cyclemat, cycle))
        end
    end
end

# Tests network analysis functions on a second network (by comparing to manually computed outputs).
let
    rn2 = @reaction_network begin
        (k₁, k₂), E + S1 <--> ES1
        (k₃, k₄), E + S2 <--> ES2
        (k₅, k₆),  S2 + ES1 <--> ES1S2
        (k₆, k₇), ES1S2 <--> S1 + ES2
        k₈, ES1S2 --> E+P
        (k₉, k₁₀), S1 <--> 0
        (k₁₀, k₁₁), 0 <--> S2
        k₁₂, P --> 0
    end

    rcs, B = reactioncomplexes(rn2)
    @test length(rcs) == 12
    @test length(linkageclasses(rn2)) == 4
    @test deficiency(rn2) == 2
    @test all(==(0), linkagedeficiencies(rn2))
    @test isreversible(rn2) == false
    @test isweaklyreversible(rn2, subnetworks(rn2)) == false
    cls = conservationlaws(rn2)
    @test Catalyst.get_networkproperties(rn2).rank == 6

    k = rand(rng, numparams(rn2))
    rates = Dict(zip(parameters(rn2), k))
    @test Catalyst.iscomplexbalanced(rn2, rates) == false
end

# Tests network analysis functions on third network (by comparing to manually computed outputs).
let
    rn3 = @reaction_network begin
        (k₁, k₂), A11 <--> 0
        (k₃, k₄), A11 <--> A13
        (k₅, k₆),  0 <--> A12
        (k₆, k₇), 0 <--> A2
        k₈, A10 --> 0
        (k₉, k₁₀), A12 <--> A6
        (k₁₁, k₁₂), A6<--> A4
        (k₁₃, k₁₄), A4 <--> A3
        k₁₅, A8 --> A9
        (k₁₆,k₁₇), A8 <--> A3 + A11
        k₁₈, A9 --> A3 + A10
        k₁₉, A2+A4 --> A2 + A6
    end
    rcs, B = reactioncomplexes(rn3)
    @test length(rcs) == 15
    @test length(linkageclasses(rn3)) == 3
    @test deficiency(rn3) == 2
    @test all(==(0), linkagedeficiencies(rn3))
    @test isreversible(rn3) == false
    @test isweaklyreversible(rn3, subnetworks(rn3)) == false
    cls = conservationlaws(rn3)
    @test Catalyst.get_networkproperties(rn3).rank == 10

    k = rand(rng, numparams(rn3))
    rates = Dict(zip(parameters(rn3), k))
    @test Catalyst.iscomplexbalanced(rn3, rates) == false
end

### STRONG LINKAGE CLASS TESTS

# a) Checks that strong/terminal linkage classes are correctly found. Should identify the (A, B+C) linkage class as non-terminal, since B + C produces D
let
    rn = @reaction_network begin
        (k1, k2), A <--> B + C
        k3, B + C --> D
        k4, D --> E
        (k5, k6), E <--> 2F
        k7, 2F --> D
        (k8, k9), D + E <--> G
    end

    rcs, D = reactioncomplexes(rn)
    slcs = stronglinkageclasses(rn)
    tslcs = terminallinkageclasses(rn)
    @test length(slcs) == 3
    @test length(tslcs) == 2
    @test issubset([[1,2], [3,4,5], [6,7]], slcs)
    @test issubset([[3,4,5], [6,7]], tslcs)
end

# b) Makes the D + E --> G reaction irreversible. Thus, (D+E) becomes a non-terminal linkage class. Checks whether correctly identifies both (A, B+C) and (D+E) as non-terminal
let
    rn = @reaction_network begin
        (k1, k2), A <--> B + C
        k3, B + C --> D
        k4, D --> E
        (k5, k6), E <--> 2F
        k7, 2F --> D
        k8, D + E --> G
    end

    rcs, D = reactioncomplexes(rn)
    slcs = stronglinkageclasses(rn)
    tslcs = terminallinkageclasses(rn)
    @test length(slcs) == 4
    @test length(tslcs) == 2
    @test issubset([[1,2], [3,4,5], [6], [7]], slcs)
    @test issubset([[3,4,5], [7]], tslcs)
end

# From a), makes the B + C <--> D reaction reversible. Thus, the non-terminal (A, B+C) linkage class gets absorbed into the terminal (A, B+C, D, E, 2F) linkage class, and the terminal linkage classes and strong linkage classes coincide.
let
    rn = @reaction_network begin
        (k1, k2), A <--> B + C
        (k3, k4), B + C <--> D
        k5, D --> E
        (k6, k7), E <--> 2F
        k8, 2F --> D
        (k9, k10), D + E <--> G
    end

    rcs, D = reactioncomplexes(rn)
    slcs = stronglinkageclasses(rn)
    tslcs = terminallinkageclasses(rn)
    @test length(slcs) == 2
    @test length(tslcs) == 2
    @test issubset([[1,2,3,4,5], [6,7]], slcs)
    @test issubset([[1,2,3,4,5], [6,7]], tslcs)
end

# Simple test for strong and terminal linkage classes
let
    rn = @reaction_network begin
        (k1, k2), A <--> 2B
        k3, A --> C + D
        (k4, k5), C + D <--> E
        k6, 2B --> F
        (k7, k8), F <--> 2G
        (k9, k10), 2G <--> H
        k11, H --> F
    end

    rcs, D = reactioncomplexes(rn)
    slcs = stronglinkageclasses(rn)
    tslcs = terminallinkageclasses(rn)
    @test length(slcs) == 3
    @test length(tslcs) == 2
    @test issubset([[1,2], [3,4], [5,6,7]], slcs)
    @test issubset([[3,4], [5,6,7]], tslcs)
end

# Cycle Test: Open Reaction Network
let
    rn = @reaction_network begin
        k1, 0 --> X1
        k2, X1 --> 0
        k3, X1 --> X2
        (k4, k5), X2 <--> X3
        (k6, k7), X3 <--> 0
    end

    # 0 --> X1 --> X2 --> X3 --> 0
    cycle = [1, 0, 1, 1, 0, 1, 0]
    cyclemat = Catalyst.cycles(rn)
    @test rank(cyclemat) == rank(hcat(cyclemat, cycle))
end

# From stoichiometric matrix. Reference: Trinh, 2008, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2909134/
let
   S = [1 -1 0 0 -1 0 0 0 0;
        0 0 0 0 1 -1 -1 -1 0;
        0 1 -1 0 0 1 0 0 0;
        0 0 1 0 0 0 0 0 -1;
        0 0 1 -1 0 0 2 0 0]

   EFMs = [1 0 1 1 0 1 1 1;
           1 0 0 1 0 0 1 0;
           0 1 0 1 0 0 0 1;
           0 1 0 1 2 2 2 1;
           0 0 1 0 0 1 0 1;
           -1 1 0 0 0 0 -1 1;
           0 0 0 0 1 1 1 0;
           1 -1 1 0 -1 0 0 0;
           0 1 0 1 0 0 0 1]

   cyclemat = Catalyst.cycles(S)
   for i in 1:size(EFMs, 2)
       EFM = EFMs[:, i]
       @test rank(cyclemat) == rank(hcat(cyclemat, EFM))
   end
end

# No cycles should exist in the following network (the graph is treelike and irreversible)

let
    rn = @reaction_network begin
        k1, A + B --> C + D
        k2, C + D --> E + F
        k3, C + D --> 2G + H
        k4, 2G + H --> 3I
        k5, E + F --> J
        k6, 3I --> K
    end

    S = netstoichmat(rn)
    cyclemat = Catalyst.cycles(S)
    @test isempty(cyclemat)
end


### Other Network Properties Tests ###

# Tests outgoing complexes matrices (1).
# Checks using dense and sparse representation.
let
    # Declares network.
    rs = @reaction_network begin
        k1, X1 + X2 --> X3 + X4
        (k2,k2), X3 + X4 <--> X1
        k3, X1 --> X2
        k4, X1 + X2 --> X2
    end

    # Compares to manually computed matrix.
    cmplx_out_mat = [
        -1 0 0 0 -1;
        0 -1 0 0 0;
        0 0 -1 -1 0;
        0 0 0 0 0;
    ]
    complexoutgoingmat(rs) == cmplx_out_mat
    complexoutgoingmat(rs; sparse = true) == sparse(cmplx_out_mat)
end

# Tests outgoing complexes matrices (2).
# Checks using dense and sparse representation.
let
    # Declares network.
    rs = @reaction_network begin
        k1, X1 --> X2
        k2, X2 --> X3
        k3, X3 --> X4
        k4, X3 --> X5
        k5, X2 --> X1
        k6, X1 --> X2
    end

    # Compares to manually computed matrix.
    cmplx_out_mat = [
        -1 0 0 0 0 -1;
        0 -1 0 0 -1 0;
        0 0 -1 -1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
    ]
    complexoutgoingmat(rs)
    complexoutgoingmat(rs; sparse = true) == sparse(cmplx_out_mat)
end

### Tests for the matrices and vectors that are in the species-formation rate function
# ẋ = S * K * Φ(t)
# ẋ = Y * A_K * Φ(t)

@test_broken let
    return false # Sometimes some fail due to https://github.com/JuliaSymbolics/Symbolics.jl/issues/1739
    MAPK = @reaction_network MAPK begin
        (k₁, k₂),KKK + E1 <--> KKKE1
        k₃, KKKE1 --> KKK_ + E1
        (k₄, k₅), KKK_ + E2 <--> KKKE2
        k₆, KKKE2 --> KKK + E2
        (k₇, k₈), KK + KKK_ <--> KK_KKK_
        k₉, KK_KKK_ --> KKP + KKK_
        (k₁₀, k₁₁), KKP + KKK_ <--> KKPKKK_
        k₁₂, KKPKKK_ --> KKPP + KKK_
        (k₁₃, k₁₄), KKP + KKPase <--> KKPKKPase
        k₁₅, KKPPKKPase --> KKP + KKPase
        k₁₆,KKPKKPase --> KK + KKPase
        (k₁₇, k₁₈), KKPP + KKPase <--> KKPPKKPase
        (k₁₉, k₂₀), KKPP + K <--> KKPPK
        k₂₁, KKPPK --> KKPP + KP
        (k₂₂, k₂₃), KKPP + KP <--> KPKKPP
        k₂₄, KPKKPP --> KPP + KKPP
        (k₂₅, k₂₆), KP + KPase <--> KPKPase
        k₂₇, KKPPKPase --> KP + KPase
        k₂₈, KPKPase --> K + KPase
        (k₂₉, k₃₀), KPP + KPase <--> KKPPKPase
    end

    Φ = Catalyst.massactionvector(MAPK)
    specs = species(MAPK)
    truevec = [MAPK.KKK * MAPK.E1,
               MAPK.KKKE1,
               MAPK.KKK_ * MAPK.E1,
               MAPK.KKK_ * MAPK.E2,
               MAPK.KKKE2,
               MAPK.KKK * MAPK.E2,
               MAPK.KK * MAPK.KKK_,
               MAPK.KK_KKK_,
               MAPK.KKP * MAPK.KKK_,
               MAPK.KKPKKK_,
               MAPK.KKPP * MAPK.KKK_,
               MAPK.KKP * MAPK.KKPase,
               MAPK.KKPKKPase,
               MAPK.KKPPKKPase,
               MAPK.KK * MAPK.KKPase,
               MAPK.KKPP * MAPK.KKPase,
               MAPK.KKPP * MAPK.K,
               MAPK.KKPPK,
               MAPK.KKPP * MAPK.KP,
               MAPK.KPKKPP,
               MAPK.KPP * MAPK.KKPP,
               MAPK.KP * MAPK.KPase,
               MAPK.KPKPase,
               MAPK.KKPPKPase,
               MAPK.K * MAPK.KPase,
               MAPK.KPP * MAPK.KPase,
              ]
    @test isequal(Φ, truevec)

    K = Catalyst.fluxmat(MAPK)
    # Construct flux matrix from incidence matrix
    mat = Matrix{Any}(zeros(30, 26))
    D = incidencemat(MAPK)
    rates = reactionrates(MAPK)
    for (i, col) in enumerate(eachcol(D))
        sub = findfirst(==(-1), col)
        mat[i, sub] = rates[i]
    end
    @test isequal(K, mat)
    @test isequal(K[1, 1], MAPK.k₁)
    @test all(==(0), K[1, 2:end])
    @test isequal(K[2, 2], MAPK.k₂)
    @test all(==(0), vcat(K[2,1], K[2,3:end]))
    @test isequal(K[3, 2], MAPK.k₃)
    @test all(==(0), vcat(K[3,1], K[3,3:end]))
    @test count(k -> !isequal(k, 0), K) == length(reactions(MAPK))

    A_k = Catalyst.laplacianmat(MAPK)
    @test all(col -> sum(col) == 0, eachcol(A_k))

    S = netstoichmat(MAPK)
    Y = complexstoichmat(MAPK)
    @test isequal(S*K, Y*A_k)

    eqs = Catalyst.assemble_oderhs(MAPK, specs)
    @test all(_iszero, simplify(eqs - S*K*Φ))
    @test all(_iszero, simplify(eqs - Y*A_k*Φ))

    # Test using numbers
    k = rand(rng, numparams(MAPK))
    ratevec = collect(zip(parameters(MAPK), k))
    ratemap = Dict(ratevec)
    ratetup = Tuple(ratevec)

    @test Catalyst.fluxmat(MAPK, ratemap) == Catalyst.fluxmat(MAPK, ratevec) == Catalyst.fluxmat(MAPK, ratetup)

    K = Catalyst.fluxmat(MAPK, ratemap)
    A_k = Catalyst.laplacianmat(MAPK, ratemap)
    @test all(col -> sum(col) == 0, eachcol(A_k))

    numeqs = similar(eqs)
    for i in 1:length(eqs)
        numeqs[i] = substitute(eqs[i], ratemap)
    end
    @test all(_iszero, simplify(numeqs - S*K*Φ))
    @test all(_iszero, simplify(numeqs - Y*A_k*Φ))
end

# Test handling for weird complexes and combinatoric rate laws.
@test_broken let
    return false # Sometimes some fail due to https://github.com/JuliaSymbolics/Symbolics.jl/issues/1739
    rn = @reaction_network begin
        k1, 2X + Y + 3Z --> ∅
        (k2, k3), 2Y + 2Z <--> 3X
    end

    Φ = Catalyst.massactionvector(rn)
    specs = species(rn)
    crvec = [rn.X^2/2 * rn.Y * rn.Z^3/6,
             1.,
             rn.Y^2/2 * rn.Z^2/2,
             rn.X^3/6]
    @test isequal(Φ, crvec)
    ncrvec = [rn.X^2 * rn.Y * rn.Z^3,
              1.,
              rn.Y^2 * rn.Z^2,
              rn.X^3]
    Φ_2 = Catalyst.massactionvector(rn; combinatoric_ratelaws = false)
    @test isequal(Φ_2, ncrvec)

    # Test that the ODEs generated are the same.
    eqs = Catalyst.assemble_oderhs(rn, specs)
    S = netstoichmat(rn)
    Y = complexstoichmat(rn)
    K = fluxmat(rn)
    A_k = laplacianmat(rn)
    @test all(_iszero, simplify(eqs - S*K*Φ))
    @test all(_iszero, simplify(eqs - Y*A_k*Φ))

    eq_ncr = Catalyst.assemble_oderhs(rn, specs; combinatoric_ratelaws = false)
    @test all(_iszero, simplify(eq_ncr - S*K*Φ_2))
    @test all(_iszero, simplify(eq_ncr - Y*A_k*Φ_2))

    # Test that the ODEs with rate constants are the same.
    k = rand(rng, numparams(rn))
    ratevec = collect(zip(parameters(rn), k))
    ratemap = Dict(ratevec)
    K = fluxmat(rn, ratemap)
    A_k = laplacianmat(rn, ratemap)

    numeqs = similar(eqs)
    for i in 1:length(eqs)
        numeqs[i] = substitute(eqs[i], ratemap)
    end
    # Broken but the difference is just numerical, something on the order of 1e-17 times a term
    @test all(_iszero, simplify(numeqs - S*K*Φ))
    @test all(_iszero, simplify(numeqs - Y*A_k*Φ))

    numeqs_ncr = similar(eq_ncr)
    for i in 1:length(eq_ncr)
        numeqs_ncr[i] = substitute(eq_ncr[i], ratemap)
    end
    @test all(_iszero, simplify(numeqs_ncr - S*K*Φ_2))
    @test all(_iszero, simplify(numeqs_ncr - Y*A_k*Φ_2))

    # Test that handling of species concentrations is correct.
    u0vec = [:X => 3., :Y => .5, :Z => 2.]
    u0map = Dict(u0vec)
    u0tup = Tuple(u0vec)

    Φ = Catalyst.massactionvector(rn, u0vec)
    @test isequal(Φ[1], 3.)
    Φ_2 = Catalyst.massactionvector(rn, u0tup; combinatoric_ratelaws = false)
    @test isequal(Φ_2[1], 36.)
    Φ = Catalyst.massactionvector(rn, u0map)
    @test isequal(Φ[1], 3.)

    # Test full simplification.
    u0map = symmap_to_varmap(rn, u0map)
    numeqs = [substitute(eq, u0map) for eq in numeqs]
    @test isapprox(numeqs, S*K*Φ)
    @test isapprox(numeqs, Y*A_k*Φ)

    numeqs_ncr = [substitute(eq, u0map) for eq in numeqs_ncr]
    @test isapprox(numeqs_ncr, S*K*Φ_2)
    @test isapprox(numeqs_ncr, Y*A_k*Φ_2)
end
