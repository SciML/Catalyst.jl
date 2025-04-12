### File declaring various reaction networks for the tests to be run on ###

# Declares the vectors which contains the various test sets.
reaction_networks_standard = Vector{ReactionSystem}(undef, 10)
reaction_networks_hill = Vector{ReactionSystem}(undef, 10)
reaction_networks_conserved = Vector{ReactionSystem}(undef, 10)
reaction_network_conslaws = Vector{Matrix{Int}}(undef, 10)
reaction_networks_real = Vector{ReactionSystem}(undef, 4)
reaction_networks_weird = Vector{ReactionSystem}(undef, 10)

### Standard reaction networks. ###
reaction_networks_standard[1] = @reaction_network rns1 begin
    (p1, p2, p3), ∅ → (X1, X2, X3)
    (k1, k2), X2 ⟷ X1 + 2X3
    (k3, k4), X1 ⟷ X3
    (d1, d2, d3), (X1, X2, X3) → ∅
end

reaction_networks_standard[2] = @reaction_network rns2 begin
    mmr(X2, v1, K1), ∅ → X1
    mm(X1, v2, K2), ∅ → X2
    d, X1 + X2 → ∅
end

reaction_networks_standard[3] = @reaction_network rns3 begin
    mm(X2, v1, K1), ∅ → X1
    mm(X3, v2, K2), ∅ → X2
    (k1, k2), X1 ⟷ X3
    (k3, k4), X3 + X2 ⟷ X4 + X1
    d, (X1, X2, X3, X4) → ∅
end

reaction_networks_standard[4] = @reaction_network rns4 begin
    mmr(X4, v1, K1), ∅ → X1
    mmr(X1, v2, K2), ∅ → X2
    mmr(X2, v3, K3), ∅ → X3
    mmr(X3, v4, K4), ∅ → X4
    (d1, d2, d3, d4), (X1, X2, X3, X4) → ∅
end

reaction_networks_standard[5] = @reaction_network rns5 begin
    p, ∅ → X1
    (k1, k2), X1 ⟷ X2
    (k3, k4), X2 ⟷ X3
    (k5, k6), X3 ⟷ X4
    d, X4 → ∅
end

reaction_networks_standard[6] = @reaction_network rns6 begin
    (p1, p2), ∅ → (X1, X2)
    (k1, k2), 2X1 ⟷ X3
    (k3, k4), X2 + X3 ⟷ 4X4
    (k5, k6), X4 + X1 ⟷ 2X3
    d, (X1, X2, X3, X4) → ∅
end

reaction_networks_standard[7] = @reaction_network rns7 begin
    (p1, p2, p3), ∅ → (X1, X2, X3)
    (k1, k2), X1 + 2X2 ⟷ X4
    (mm(X3, v1, K1), k3), X4 ⟷ X5
    (d1, d2, d3, d4, d5), (X1, X2, X3, X4, X5) → ∅
end

reaction_networks_standard[8] = @reaction_network rns8 begin
    p, ∅ → 2X1
    k1, X1 → X2
    k2, X2 → X3
    k3, X2 → X3
    d, X3 → ∅
end

reaction_networks_standard[9] = @reaction_network rns9 begin
    (p1, p2, p3), ∅ ⟶ (X1, X2, X3)
    (d1, d2, d3), (X1, X2, X3) ⟶ ∅
    (k1, k2), X1 + X2 ⟷ X3
    (k3, k4), X3 ⟷ X4
    (k5, k6), X4 ⟷ X1 + X2
end

reaction_networks_standard[10] = @reaction_network rns10 begin
    p, ∅ ⟶ X1
    (k1, k2), (X1,X1) → X2
    (k3, k4), (X2,X2) → X3
    (k5, k6), (X3,X3) → X4
    (k7, k8), (X4,X4) → X5
    d, X5 ⟶ ∅
end

### Network with Hill functions ###.

reaction_networks_hill[1] = @reaction_network rnh1 begin
    hillr(X2, v1, K1, n1), ∅ → X1
    hillr(X1, v2, K2, n2), ∅ → X2
    (d1, d2), (X1, X2) → ∅
end

reaction_networks_hill[2] = @reaction_network rnh2 begin
    hillr(X3, v1, K1, n1), ∅ → X1
    hillr(X1, v2, K2, n2), ∅ → X2
    hillr(X2, v3, K3, n3), ∅ → X3
    (d1, d2, d3), (X1, X2, X3) → ∅
end

reaction_networks_hill[3] = @reaction_network rnh3 begin
    hillr(X2, v1, K1, n1), ∅ → X1
    hill(X1, v2, K2, n2), ∅ → X2
    d, X1 + X2 → ∅
end

reaction_networks_hill[4] = @reaction_network rnh4 begin
    hillr(X2, v1, K1, n1) * hillr(X3, v1, K1, n1), ∅ → X1
    hillr(X1, v2, K2, n2) * hillr(X3, v2, K2, n2), ∅ → X2
    hillr(X1, v3, K3, n3) * hillr(X2, v3, K3, n3), ∅ → X3
    (d1, d2, d3), (X1, X2, X3) ⟶ ∅
end

reaction_networks_hill[5] = @reaction_network rnh5 begin
    hillr(X2, v1, K1, n1) * hill(X4, v1, K1, n1), ∅ → X1
    hill(X5, v2, K2, n2), ∅ → X2
    hill(X3, v3, K3, n3), ∅ → X3
    hillr(X1, v4, K4, n4), ∅ → X4
    hill(X2, v5, K5, n5), ∅ → X5
    (k1, k2), X2 ⟷ X1 + 2X4
    (k3, k4), X4 ⟷ X3
    (k5, k6), 3X5 + X1 ⟷ X2
    (d1, d2, d3, d4, d5), (X1, X2, X3, X4, X5) ⟶ ∅
end

reaction_networks_hill[6] = @reaction_network rnh6 begin
    v / 10 + hill(X1, v, K, n), ∅ → X1
    d, X1 → ∅
end

reaction_networks_hill[7] = @reaction_network rnh7 begin
    v / 10 + hill(X1, v, K, n), ∅ → X1 + X2
    (k1, k2), X1 + X2 ↔ X3
    k3, X3 → X1
    d, (X1, X2, X3) → ∅
end

reaction_networks_hill[8] = @reaction_network rnh8 begin
    hill(X2, v1, K1, n1), ∅ → X1
    hillr(X1, v2, K2, n2) * hill(X3, v3, K3, n3), ∅ → X2
    hill(X2, v4, K4, n4), ∅ → X3
    (d1, d2, d3), (X1, X2, X3) → ∅
end

reaction_networks_hill[9] = @reaction_network rnh9 begin
    hill(X1, v1, K1, n1) * hillr(X1, v2, K2, n2), ∅ → X1
    d, X1 → ∅
end

reaction_networks_hill[10] = @reaction_network rnh10 begin
    hill(X2, v1, K1, n1), ∅ → X1
    hillr(X4, v2, K2, n2), ∅ → X2
    (k1, k2), 2X1 + X2 ⟷ X3
    (k3, k4), 2X2 + X3 ⟷ X4
    (k5, k6), X1 + X2 + X3 + X4 ⟷ X5 + X6
    (d1, d2), (X5, X6) → ∅
end

### Reaction networks were some linear combination concentrations remain fixed (steady state values depends on initial conditions). ###

reaction_networks_conserved[1] = @reaction_network rnc1 begin
    (k1, k2), X1 ↔ X2
    (k3, k4), X2 ↔ X3
    (k5, k6), X3 ↔ X1
end
reaction_network_conslaws[1] = [1 1 1]

reaction_networks_conserved[2] = @reaction_network rnc2 begin
    (k1, k2), X1 ↔ 2X1
    (k3, k4), X1 + X2 ↔ X3
    (k5, k6), X3 ↔ X2
end
reaction_network_conslaws[2] = [0 1 1]

reaction_networks_conserved[3] = @reaction_network rnc3 begin
    (k1, k2 * X5), X1 ↔ X2
    (k3 * X5, k4), X3 ↔ X4
    (p + k5 * X2 * X3, d), ∅ ↔ X5
end
reaction_network_conslaws[3] = [0 0 1 1 0; 1 1 0 0 0]

reaction_networks_conserved[4] = @reaction_network rnc4 begin
    (k1, k2), X1 + X2 ↔ X3
    (mm(X3, v, K), d), ∅ ↔ X4
end
reaction_network_conslaws[4] = [0 1 1 0; -1 1 0 0]

reaction_networks_conserved[5] = @reaction_network rnc5 begin
    (k1, k2), X1 ↔ 2X2
    (k3, k4), 2X2 ↔ 3X3
    (k5, k6), 3X3 ↔ 4X4
end
reaction_network_conslaws[5] = [12 6 4 3]

reaction_networks_conserved[6] = @reaction_network rnc6 begin
    mmr(X1, v1, K1), X1 → X2
    mmr(X2, v2, K2), X2 → X3
    mmr(X3, v3, K3), X3 → X1
end
reaction_network_conslaws[6] = [1 1 1]

reaction_networks_conserved[7] = @reaction_network rnc7 begin
    (k1, k2), X1 + X2 ↔ X3
    (mm(X3, v, K), d), ∅ ↔ X2
    (k3, k4), X2 ↔ X4
end
reaction_network_conslaws[7] = [1 0 1 0]

reaction_networks_conserved[8] = @reaction_network rnc8 begin
    (k1, k2), X1 + X2 ↔ X3
    (mm(X3, v1, K1), mm(X4, v2, K2)), X3 ↔ X4
end
reaction_network_conslaws[8] = [-1 1 0 0; 0 1 1 1]

reaction_networks_conserved[9] = @reaction_network rnc9 begin
    (k1, k2), X1 + X2 ↔ X3
    (k3, k4), X3 + X4 ↔ X5
    (k5, k6), X5 + X6 ↔ X7
end
reaction_network_conslaws[9] = [1 0 1 0 1 0 1; -1 1 0 0 0 0 0; 0 0 0 1 1 0 1;
                                   0 0 0 0 0 1 1]

reaction_networks_conserved[10] = @reaction_network rnc10 begin
    kDeg, (w, w2, w2v, v, w2v2, vP, σB, w2σB) ⟶ ∅
    kDeg, vPp ⟶ phos
    (kBw, kDw), 2w ⟷ w2
    (kB1, kD1), w2 + v ⟷ w2v
    (kB2, kD2), w2v + v ⟷ w2v2
    kK1, w2v ⟶ w2 + vP
    kK2, w2v2 ⟶ w2v + vP
    (kB3, kD3), w2 + σB ⟷ w2σB
    (kB4, kD4), w2σB + v ⟷ w2v + σB
    (kB5, kD5), vP + phos ⟷ vPp
    kP, vPp ⟶ v + phos
    v0 * ((1 + F * σB) / (K + σB)), ∅ ⟶ σB
    λW * v0 * ((1 + F * σB) / (K + σB)), ∅ ⟶ w
    λV * v0 * ((1 + F * σB) / (K + σB)), ∅ ⟶ v
end;
reaction_network_conslaws[10] = [0 0 0 0 0 0 0 0 1 1]

### Reaction networks that are actual models that have been used ###

# Brusselator.
reaction_networks_real[1] = @reaction_network rnr1 begin
    A, ∅ → X
    1, 2X + Y → 3X
    B, X → Y
    1, X → ∅
end

# The B.subtilis σV Lysozyme stress response.
reaction_networks_real[2] = @reaction_network rnr2 begin
    v0 + hill(σ, v, K, n), ∅ → (σ + A)
    deg, (σ, A, Aσ) → ∅
    (kB, kD), A + σ ↔ Aσ
    S * kC, Aσ → σ
end

# A cell cycle model
reaction_networks_real[3] = @reaction_network rnr3 begin
    k1, 0 --> Y
    k2p, Y --> 0
    k2pp * P, Y --> 0
    (k3p + k3pp * A) / (J3 + Po), Po --> P
    (k4 * m) / (J4 + P), Y + P --> Y + Po
end
#@add_constraint cc_network P+Po=1
#reaction_networks_real[3] = cc_network

# A bistable switch
reaction_networks_real[4] = @reaction_network rnr4 begin
    d, (X, Y) → ∅
    hillr(Y, v1, K1, n1), ∅ → X
    hillr(X, v2, K2, n2), ∅ → Y
end

### Reaction networks that contain weird functions, stuff, and other oddities ###

reaction_networks_weird[1] = @reaction_network rnw1 begin
    exp(p), ∅ → X
    d * exp(X), X → ∅
end

reaction_networks_weird[2] = @reaction_network rnw2 begin
    k1, ∅ → X
    k2 * log(12 + X), X → Y
    k3 * log(3 + Y), Y → Z
    log(5, 6 + k4), Z → ∅
end

reaction_networks_weird[3] = @reaction_network rnw3 begin
    d, (X, Y) → ∅
    hillr(hill(Y, v2, K2, n2), v1, K1, n1), ∅ → X
    hillr(hill(X, v1, K1, n1), v2, K2, n2), ∅ → Y
end

reaction_networks_weird[4] = @reaction_network rnw4 begin
    p, ∅ → X1
    (k1, X2^X3), X1 ⟷ X2
    (X2 / X4, k4), X2 ⟷ X3
    (k5, X1 * X2 * X4), X3 ⟷ X4
    d * X2, X4 → ∅
end

reaction_networks_weird[5] = @reaction_network rnw5 begin
    (k1 * tanh(X1), k2), X1 ↔ 2X1
    (2 + sin(X1), 3 + sin(k3)), X1 + X2 ↔ X3
    (k5, 4 + cos(X3 + X1) + sin(X2)), X3 ↔ X2
end

reaction_networks_weird[6] = @reaction_network rnw6 begin
    (p, d), ∅ ↔ X1
    (p, d), ∅ ↔ X1
    (p, d), ∅ ↔ X1
    (p, d), ∅ ↔ X1
    (p, d), ∅ ↔ X1
    (p, d), ∅ ↔ X1
end

reaction_networks_weird[7] = @reaction_network rnw7 begin
    k1, X1 → X2
    0, X2 → X3
    k2, X3 → X4
    k3, X4 → X5
end

reaction_networks_weird[8] = @reaction_network rnw8 begin
    k1 + X3, X1 → X2
    k2 * X4, X2 → X3
    k3 / (X5 + 0.01), X3 → X4
    X6 / k4, X4 → X5
    X7^k5, X5 → X6
    k6^X8, X6 → X7
    sqrt(abs(k7 * X9)), X7 → X8
    cbrt(abs(k8 + X1)), X8 → X9
    X2^3 + 2X2^2 + k9, X9 → X1
end

reaction_networks_weird[9] = @reaction_network rnw9 begin
    p, ∅ → X1
    hill(hill(hill(hill(X1, v1, K1, n1), v2, K2, n2), v3, K3, n3), v4, K4, n4), X1 → ∅
end

reaction_networks_weird[10] = @reaction_network rnw10 begin
    d, 5X1 → 4X1
end

### Gathers all networks in a single array ###
reaction_networks_all = [reaction_networks_standard...,
    reaction_networks_hill...,
    reaction_networks_conserved...,
    reaction_networks_real...,
    reaction_networks_weird...]
