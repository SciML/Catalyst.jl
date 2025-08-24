# Tests for properties from chemical reaction network theory: deficiency theorems, complex/detailed balance, etc.
using Catalyst, StructuralIdentifiability, LinearAlgebra, Test
rng = StableRNG(514)
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

# Test that incomplete rate maps error.
let
    rn = @reaction_network begin
        (k1, k2), C1 <--> C2
        (k3, k4), C2 <--> C3
        (k5, k6), C3 <--> C1
    end

    incorrect_params = Dict(:k1 => 0.5)
    @test_throws ErrorException Catalyst.iscomplexbalanced(rn, incorrect_params)
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
    @test Catalyst.adjacencymat(rn, rates_vec) == rate_mat
    @test Catalyst.adjacencymat(rn, rates_tup) == rate_mat
    @test Catalyst.adjacencymat(rn, rates_dict) == rate_mat
    @test_throws Exception Catalyst.adjacencymat(rn, rate_vals)

    # Tests that throws error in rate matrix.
    incorrect_param_dict = Dict(:k1 => 1.0)

    @test_throws ErrorException Catalyst.adjacencymat(rn, 123)
    @test_throws ErrorException Catalyst.adjacencymat(rn, incorrect_param_dict)

    @test_throws Exception Catalyst.iscomplexbalanced(rn, rates_invalid)

    # Test sparse matrix
    @test Catalyst.adjacencymat(rn, rates_vec; sparse = true) == rate_mat
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

let
    # Define a reaction network with bi-directional reactions
    non_deficient_network = @reaction_network begin
        (k1, k2), A <--> B
        (k3, k4), B <--> C
    end

    # Test: Check that the error is raised for networks with deficiency != 1
    @test_throws ErrorException Catalyst.robustspecies(non_deficient_network)
end

### Complex balance and reversibility tests ###

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
    rn4 = @reaction_network begin
        (k1, k2), C1 <--> C2
        (k3, k4), C2 <--> C3
        (k5, k6), C3 <--> C1
    end

    k = rand(rng, numparams(rn4))
    rates = Dict(zip(parameters(rn4), k))
    @test Catalyst.iscomplexbalanced(rn4, rates) == true
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
    @test Catalyst.isdetailedbalanced(rn, rates) == false
end


### DEFICIENCY THEOREMS TESTS

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

    @test Catalyst.satisfiesdeficiencyzero(rn) == true
    @test Catalyst.satisfiesdeficiencyzero(rn2) == false
    @test Catalyst.satisfiesdeficiencyzero(rn3) == false
    @test Catalyst.satisfiesdeficiencyzero(rn4) == false
end

### Detailed balance tests

# The following network is conditionally complex balanced - it only

# Reversible, forest-like deficiency zero network - should be detailed balance for any choice of rate constants.
let
    rn = @reaction_network begin
        (k1, k2), A <--> B + C
        (k3, k4), A <--> D
        (k5, k6), A + D <--> E
        (k7, k8), A + D <--> G
        (k9, k10), G <--> 2F
        (k11, k12), A + E <--> H
    end

    k1 = rand(rng, numparams(rn))
    rates1 = Dict(zip(parameters(rn), k1))
    k2 = rand(StableRNG(232), numparams(rn))
    rates2 = Dict(zip(parameters(rn), k2))

    rcs, D = reactioncomplexes(rn)
    @test Catalyst.isforestlike(rn) == true
    @test Catalyst.isdetailedbalanced(rn, rates1) == true
    @test Catalyst.isdetailedbalanced(rn, rates2) == true
end

# Simple connected reversible network
let
    rn = @reaction_network begin
        (k1, k2), A <--> B
        (k3, k4), B <--> C
        (k5, k6), C <--> A
    end

    rcs, D = reactioncomplexes(rn)
    rates1 = [:k1=>1.0, :k2=>1.0, :k3=>1.0, :k4=>1.0, :k5=>1.0, :k6=>1.0]
    @test Catalyst.isdetailedbalanced(rn, rates1) == true
    rates2 = [:k1=>2.0, :k2=>1.0, :k3=>1.0, :k4=>1.0, :k5=>1.0, :k6=>1.0]
    @test Catalyst.isdetailedbalanced(rn, rates2) == false
end

# Independent cycle tests: the following reaction entwork has 3 out-of-forest reactions.
let
    rn = @reaction_network begin
        (k1, k2), A <--> B + C
        (k3, k4), A <--> D
        (k5, k6), B + C <--> D
        (k7, k8), A + D <--> E
        (k9, k10), G <--> 2F
        (k11, k12), A + D <--> G
        (k13, k14), G <--> E
        (k15, k16), 2F <--> E
        (k17, k18), A + E <--> H
    end

    rcs, D = reactioncomplexes(rn)
    k = rand(rng, numparams(rn))
    p = parameters(rn)
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.isdetailedbalanced(rn, rates) == false

    # Adjust rate constants to obey the independent cycle conditions.
    rates[p[6]] = rates[p[1]]*rates[p[4]]*rates[p[5]] / (rates[p[2]]*rates[p[3]])
    rates[p[14]] = rates[p[13]]*rates[p[11]]*rates[p[8]] / (rates[p[12]]*rates[p[7]])
    rates[p[16]] = rates[p[8]]*rates[p[15]]*rates[p[9]]*rates[p[11]] / (rates[p[7]]*rates[p[12]]*rates[p[10]])
    @test Catalyst.isdetailedbalanced(rn, rates) == true
end

# Deficiency two network: the following reaction network must satisfy both the independent cycle conditions and the spanning forest conditions
let
    rn = @reaction_network begin
        (k1, k2), 3A <--> A + 2B
        (k3, k4), A + 2B <--> 3B
        (k5, k6), 3B <--> 2A + B
        (k7, k8), 2A + B <--> 3A
        (k9, k10), 3A <--> 3B
    end

    rcs, D = reactioncomplexes(rn)
    @test Catalyst.edgeindex(D, 1, 2) == 1
    @test Catalyst.edgeindex(D, 4, 3) == 6
    k = rand(rng, numparams(rn))
    p = parameters(rn)
    rates = Dict(zip(parameters(rn), k))
    @test Catalyst.isdetailedbalanced(rn, rates) == false

    # Adjust rate constants to fulfill independent cycle conditions.
    rates[p[8]] = rates[p[7]]*rates[p[5]]*rates[p[9]] / (rates[p[6]]*rates[p[10]])
    rates[p[3]] = rates[p[2]]*rates[p[4]]*rates[p[9]] / (rates[p[1]]*rates[p[10]])
    @test Catalyst.isdetailedbalanced(rn, rates) == false
    # Should still fail - doesn't satisfy spanning forest conditions.

    # Adjust rate constants to fulfill spanning forest conditions.
    cons = rates[p[6]] / rates[p[5]]
    rates[p[1]] = rates[p[2]] * cons
    rates[p[9]] = rates[p[10]] * cons^(3/2)
    rates[p[8]] = rates[p[7]]*rates[p[5]]*rates[p[9]] / (rates[p[6]]*rates[p[10]])
    rates[p[3]] = rates[p[2]]*rates[p[4]]*rates[p[9]] / (rates[p[1]]*rates[p[10]])
    @test Catalyst.isdetailedbalanced(rn, rates) == true
end
