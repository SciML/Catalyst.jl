### Prepares Tests ###

# Fetch packages.
using Catalyst, LinearAlgebra, Test, SparseArrays

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
    # i=0;
    # for lcs in linkageclasses(MAPK)
    #     i=i+1
    #     println("Linkage no ",i)
    #     for comps in rcs[lcs]
    #         if comps.speciesids ≠ Int64[]
    #             println(sum(species(rn2)[comps.speciesids]))
    #         else
    #             println("0")
    #         end
    #     end
    #     println("-----------")
    # end

    # Testing if cycles identifies reversible reactions as cycles (one forward, one reverse) 
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
        (k₆, k₇), ES1S2 --> S1 + ES2
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
    # i=0;
    # for lcs in linkageclasses(rn2)
    #     i=i+1
    #     println("Linkage no ",i)
    #     for comps in rcs[lcs]
    #         if comps.speciesids ≠ Int64[]
    #             println(sum(species(rn2)[comps.speciesids]))
    #         else
    #             println("0")
    #         end
    #     end
    #     println("-----------")
    # end
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
    # i=0;
    # for lcs in linkageclasses(rn3)
    #     i=i+1
    #     println("Linkage no ",i)
    #     for comps in rcs[lcs]
    #         if comps.speciesids ≠ Int64[]
    #             println(sum(species(rn3)[comps.speciesids]))
    #         else
    #             println("0")
    #         end
    #     end
    #     println("-----------")
    # end
end

let
    rn4 = @reaction_network begin
        (k1, k2), C1 <--> C2
        (k3, k4), C2 <--> C3
        (k5, k6), C3 <--> C1
    end

    k = rand(rng, numparams(rn4))
    rates = Dict(zip(parameters(rn4), k))
    @test Catalyst.iscomplexbalanced(rn4, rates) == true
end
    
### Tests Reversibility ###

# Test function.
function testreversibility(rn, B, rev, weak_rev)
    @test isreversible(rn) == rev
    subrn = subnetworks(rn)
    @test isweaklyreversible(rn, subrn) == weak_rev
end

# Tests reversibility for networks with known reversibility.
let
    rn = @reaction_network begin
        (k2, k1), A1 <--> A2 + A3
        k3, A2 + A3 --> A4
        k4, A4 --> A5
        (k6, k5), A5 <--> 2A6
        k7, 2A6 --> A4
        k8, A4 + A5 --> A7
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)

    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == false 
end

let
    rn = @reaction_network begin
        (k2, k1), A1 <--> A2 + A3
        k3, A2 + A3 --> A4
        k4, A4 --> A5
        (k6, k5), A5 <--> 2A6
        k7, A4 --> 2A6
        (k9, k8), A4 + A5 <--> A7
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)

    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == false 
end

let
    rn = @reaction_network begin
        k1, A --> B
        k2, A --> C
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)
    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == false 
end

let
    rn = @reaction_network begin
        k1, A --> B
        k2, A --> C
        k3, B + C --> 2A
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)

    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == false 
end

let
    rn = @reaction_network begin
        (k2, k1), A <--> 2B
        (k4, k3), A + C <--> D
        k5, D --> B + E
        k6, B + E --> A + C
    end
    rev = false
    weak_rev = true
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)

    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == true 
end

let
    rn = @reaction_network begin
        (k2, k1), A + E <--> AE
        k3, AE --> B + E
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)

    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == false 
end

let
    rn = @reaction_network begin
        (k2, k1), A + E <--> AE
        (k4, k3), AE <--> B + E
    end
    rev = true
    weak_rev = true
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)

    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == true  
end

let
    rn = @reaction_network begin (k2, k1), A + B <--> 2A end
    rev = true
    weak_rev = true
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)

    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == true 
end

let
    rn = @reaction_network begin
        k1, A + B --> 3A
        k2, 3A --> 2A + C
        k3, 2A + C --> 2B
        k4, 2B --> A + B
    end
    rev = false
    weak_rev = true
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)

    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == true 
end

let
    rn = @reaction_network begin
        (k2, k1), A + E <--> AE
        (k4, k3), AE <--> B + E
        k5, B --> 0
        k6, 0 --> A
    end
    rev = false
    weak_rev = false
    testreversibility(rn, reactioncomplexes(rn)[2], rev, weak_rev)

    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == false 
end

let
    rn = @reaction_network begin
        k1, 3A + 2B --> 3C 
        k2, B + 4D --> 2E
        k3, 2E --> 3C
        (k4, k5), B + 4D <--> 3A + 2B
        k6, F --> B + 4D
        k7, 3C --> F
    end

    k = rand(rng, numparams(rn))
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.iscomplexbalanced(rn, rates) == true 
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
        (k8, k9), D + E --> G
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

# Tests that `iscomplexbalanced` works for different rate inputs.
# Tests that non-valid rate input yields and error
let
    # Declares network.
    rn = @reaction_network begin
        k1, 3A + 2B --> 3C 
        k2, B + 4D --> 2E
        k3, 2E --> 3C
        (k4, k5), B + 4D <--> 3A + 2B
        k6, F --> B + 4D
        k7, 3C --> F
    end

    # Declares rate alternatives.
    k = rand(rng, numparams(rn))
    rates_vec = Pair.(parameters(rn), k)
    rates_tup = Tuple(rates_vec)
    rates_dict = Dict(rates_vec)
    rates_invalid = k

    # Tests that inputs are handled correctly.
    @test Catalyst.iscomplexbalanced(rn, rates_vec) == Catalyst.iscomplexbalanced(rn, rates_tup)
    @test Catalyst.iscomplexbalanced(rn, rates_tup) == Catalyst.iscomplexbalanced(rn, rates_dict)
    @test_throws Exception Catalyst.iscomplexbalanced(rn, k)
end

# Tests rate matrix computation for various input types.
let
    # Declares network and its known rate matrix.
    rn = @reaction_network begin
        (k2, k1), A1 <--> A2 + A3
        k3, A2 + A3 --> A4
        k4, A4 --> A5
        (k6, k5), A5 <--> 2A6
        k7, 2A6 --> A4
        k8, A4 + A5 --> A7
    end
    rate_mat = [
        0.0  1.0  0.0  0.0  0.0  0.0  0.0;
        2.0  0.0  3.0  0.0  0.0  0.0  0.0;
        0.0  0.0  0.0  4.0  0.0  0.0  0.0;
        0.0  0.0  0.0  0.0  5.0  0.0  0.0;
        0.0  0.0  7.0  6.0  0.0  0.0  0.0;
        0.0  0.0  0.0  0.0  0.0  0.0  8.0;
        0.0  0.0  0.0  0.0  0.0  0.0  0.0;
    ]

    # Declares rate alternatives.
    rate_vals = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    rates_vec = Pair.(parameters(rn), rate_vals)
    rates_tup = Tuple(rates_vec)
    rates_dict = Dict(rates_vec)
    rates_invalid = reshape(rate_vals, 1, 8)

    # Tests that all input types generates the correct rate matrix.
    Catalyst.ratematrix(rn, rate_vals) == rate_mat
    Catalyst.ratematrix(rn, rates_vec) == rate_mat
    Catalyst.ratematrix(rn, rates_tup) == rate_mat
    Catalyst.ratematrix(rn, rates_dict) == rate_mat
    @test_throws Exception Catalyst.iscomplexbalanced(rn, rates_invalid)
end

### CONCENTRATION ROBUSTNESS TESTS

# Check whether concentration-robust species are correctly identified for two well-known reaction networks: the glyoxylate IDHKP-IDH system, and the EnvZ_OmpR signaling pathway. 

let
    IDHKP_IDH = @reaction_network begin
        (k1, k2), EIp + I <--> EIpI
        k3, EIpI --> EIp + Ip
        (k4, k5), E + Ip <--> EIp
        k6, EIp --> E + I
    end

    @test Catalyst.robustspecies(IDHKP_IDH) == [2]
end

let
    EnvZ_OmpR = @reaction_network begin
        (k1, k2), X <--> XT
        k3, XT --> Xp
        (k4, k5), Xp + Y <--> XpY
        k6, XpY --> X + Yp
        (k7, k8), XT + Yp <--> XTYp
        k9, XTYp --> XT + Y
    end

    @test Catalyst.robustspecies(EnvZ_OmpR) == [6]
end

### DEFICIENCY ONE TESTS

# Fails because there are two terminal linkage classes in the linkage class
let 
    rn = @reaction_network begin
        k1, A + B --> 2B
        k2, A + B --> 2A
    end

    @test Catalyst.satisfiesdeficiencyone(rn) == false
end

# Fails because linkage deficiencies do not sum to total deficiency
let 
    rn = @reaction_network begin
        (k1, k2), A <--> 2A
        (k3, k4), A + B <--> C
        (k5, k6), C <--> B 
    end

    @test Catalyst.satisfiesdeficiencyone(rn) == false
end

# Fails because a linkage class has deficiency two
let 
    rn = @reaction_network begin
        k1, 3A --> A + 2B
        k2, A + 2B --> 3B
        k3, 3B --> 2A + B
        k4, 2A + B --> 3A
    end

    @test Catalyst.satisfiesdeficiencyone(rn) == false
end

let
    rn = @reaction_network begin
        (k1, k2), 2A <--> D
        (k3, k4), D <--> A + B
        (k5, k6), A + B <--> C
        (k7, k8), C <--> 2B
        (k9, k10), C + D <--> E + F
        (k11, k12), E + F <--> H
        (k13, k14), H <--> C + E
        (k15, k16), C + E <--> D + F
        (k17, k18), A + D <--> G
        (k19, k20), G <--> B + H
    end

    @test Catalyst.satisfiesdeficiencyone(rn) == true
end

### Some tests for deficiency zero networks. 

let
    rn = @reaction_network begin
        (k1, k2), A <--> 2B
        (k3, k4), A + C <--> D
        k5, D --> B + E
        k6, B + E --> A + C
    end

    # No longer weakly reversible
    rn2 = @reaction_network begin
        (k1, k2), A <--> 2B
        (k3, k4), A + C <--> D
        k5, B + E --> D
        k6, B + E --> A + C
    end

    # No longer weakly reversible
    rn3 = @reaction_network begin
        k1, A --> 2B
        (k3, k4), A + C <--> D
        k5, D --> B + E 
        k6, B + E --> A + C
    end

    # Weakly reversible but deficiency one
    rn4 = @reaction_network begin
        (k1, k2), A <--> 2A
        (k3, k4), A + B <--> C
        (k5, k6), C <--> B
    end

    @test satisfiesdeficiencyzero(rn) == true
    @test satisfiesdeficiencyzero(rn2) == false 
    @test satisfiesdeficiencyzero(rn3) == false 
    @test satisfiesdeficiencyzero(rn4) == false
end

