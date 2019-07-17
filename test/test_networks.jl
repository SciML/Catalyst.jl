### File declaring all the various reaction networks the package is tested on ###
#A large number of reaction networks are added to test as many different cases as possible.
#More test networks can be added favourably.

#Decalares 3 vectors to contain the various test sets (and a fourth with fixed concentrations coupled with the last test set).
reaction_networks_standard = Vector{DiffEqBase.AbstractReactionNetwork}(undef,10)
reaction_networks_hill = Vector{DiffEqBase.AbstractReactionNetwork}(undef,10)
reaction_networks_fixed_conc = Vector{DiffEqBase.AbstractReactionNetwork}(undef,10)
min_reaction_networks_standard = Vector{DiffEqBase.AbstractReactionNetwork}(undef,10)

#Standard reaction networks.
reaction_networks_standard[1] = @reaction_network begin
    (p1,p2,p3), ∅ → (X1,X2,X3)
    (k1,k2), X2 ⟷ X1 + 2X3
    (k3,k4), X1 ⟷ X3
    (d1,d2,d3), (X1,X2,X3) → ∅
end p1 p2 p3 k1 k2 k3 k4 d1 d2 d3

reaction_networks_standard[2] = @reaction_network begin
    mmR(X2,v1,K1), ∅ → X1
    mm(X1,v2,K2), ∅ → X2
    d, X1+X2 → ∅
end v1 K1 v2 K2 d

reaction_networks_standard[3] = @reaction_network begin
    mm(X2,v1,K1), ∅ → X1
    mm(X3,v2,K2), ∅ → X2
    (k1,k2), X1 ⟷ X4
    (k3,k4), X4 + X2 ⟷ X3 +X1
    d, (X1,X2,X3,X4) → ∅
end v1 K1 v2 K2 k1 k2 k3 k4 d

reaction_networks_standard[4] = @reaction_network begin
    mmR(X4,v1,K1), ∅ → X1
    mmR(X1,v2,K2), ∅ → X2
    mmR(X2,v3,K3), ∅ → X3
    mmR(X3,v4,K4), ∅ → X4
    (d1,d2,d3,d4), (X1,X2,X3,X4) → ∅
end v1 K1 v2 K2 v3 K3 v4 K4 d1 d2 d3 d4

reaction_networks_standard[5] = @reaction_network begin
    p, ∅ → X1
    (k1,k2), X1 ⟷ X2
    (k3,k4), X2 ⟷ X3
    (k5,k6), X3 ⟷ X4
    d, X4 → ∅
end p k1 k2 k3 k4 k5 k6 d

reaction_networks_standard[6] = @reaction_network begin
    (p1,p2), ∅ → (X1,X2)
    (k1,k2), 2X1 ⟷ X3
    (k3,k4), X2 + X3 ⟷ 4X4
    (k5,k6), X4 + X1 ⟷ 2X3
    d, (X1,X2,X3,X4) → ∅
end p1 p2 k1 k2 k3 k4 k5 k6 d

reaction_networks_standard[7] = @reaction_network begin
    (p1,p2,p3), ∅ → (X1,X2,X3)
    (k1,k2), X1 + 2X2 ⟷ X4
    (mm(X3,v1,K1),k3), X4 ⟷ X5
    (d1,d2,d3,d4,d5), (X1,X2,X3,X4,X5) → ∅
end p1 p2 p3 k1 k2 k3 v1 K1 d1 d2 d3 d4 d5

reaction_networks_standard[8] = @reaction_network begin
    p, ∅ → 2X1
    k1, X1 → X2
    (k2, k3), X2 → X3
    d, X3 → ∅
end p k1 k2 k3 d

reaction_networks_standard[9] = @reaction_network begin
    (p1,p2,p3), ∅ ⟶ (X1,X2,X3)
    (d1,d2,d3), (X1,X2,X3) ⟶ ∅
    (k1,k2), X1 + X2 ⟷ X3
    (k3,k4), X3 ⟷ X4
    (k5,k6), X4 ⟷ X1 + X2
end p1 p2 p3 d1 d2 d3 k1 k2 k3 k4 k5 k6

reaction_networks_standard[10] = @reaction_network begin
    p, ∅ ⟶ X1
    (k1,k2), X1 → X2
    (k3,k4), X2 → X3
    (k5,k6), X3 → X4
    (k7,k8), X4 → X5
    d, X5  ⟶ ∅
end p k1 k2 k3 k4 k5 k6 k7 k8 d


#Hill type reaction networks (these contaisn a parameter in the exponent and should have integer parameters).
reaction_networks_hill[1] = @reaction_network begin
    hillR(X2,v1,K1,n1),   ∅ → X1
    hillR(X1,v2,K2,n2),   ∅ → X2
    (d1,d2),          (X1,X2) → ∅
end v1 v2 K1 K2 n1 n2 d1 d2

reaction_networks_hill[2] = @reaction_network begin
    hillR(X3,v1,K1,n1),   ∅ → X1
    hillR(X1,v2,K2,n2),   ∅ → X2
    hillR(X2,v3,K3,n3),   ∅ → X3
    (d1,d2,d3),          (X1,X2,X3) → ∅
end v1 v2 v3 K1 K2 K3 n1 n2 n3 d1 d2 d3

reaction_networks_hill[3] = @reaction_network begin
    hillR(X2,v1,K1,n1), ∅ → X1
    hill(X1,v2,K2,n2), ∅ → X2
    d, X1+X2 → ∅
end v1 K1 n1 v2 K2 n2 d

reaction_networks_hill[4] = @reaction_network begin
    hillR(X2,v1,K1,n1)*hillR(X3,v1,K1,n1), ∅ → X1
    hillR(X1,v2,K2,n2)*hillR(X3,v2,K2,n2), ∅ → X2
    hillR(X1,v3,K3,n3)*hillR(X2,v3,K3,n3), ∅ → X3
    (d1,d2,d3), (X1,X2,X3)  ⟶ ∅
end v1 K1 n1 v2 K2 n2 v3 K3 n3 d1 d2 d3

reaction_networks_hill[5] = @reaction_network begin
    hillR(X2,v1,K1,n1)*hill(X4,v1,K1,n1), ∅ → X1
    hill(X5,v2,K2,n2), ∅ → X2
    hill(X3,v3,K3,n3), ∅ → X3
    hillR(X1,v4,K4,n4), ∅ → X4
    hill(X2,v5,K5,n5), ∅ → X5
    (k1,k2), X2 ⟷ X1 + 2X4
    (k3,k4), X4 ⟷ X3
    (k5,k6), 3X5 + X1 ⟷ X2
    (d1,d2,d3,d4,d5), (X1,X2,X3,X4,X5)  ⟶ ∅
end v1 K1 n1 v2 K2 n2 v3 K3 n3 v4 K4 n4 v5 K5 n5 k1 k2 k3 k4 k5 k6 d1 d2 d3 d4 d5

reaction_networks_hill[6] = @reaction_network begin
    v/10+hill(X1,v,K,n), ∅ → X1
    d, X1 → ∅
end v K n d

reaction_networks_hill[7] = @reaction_network begin
    v/10 + hill(X1,v,K,n), ∅ → X1+X2
    (k1,k2), X1 + X2 ↔ X3
    k3, X3 → X1
    d, (X1,X2,X3) → ∅
end v K n k1 k2 k3 d

reaction_networks_hill[8] = @reaction_network begin
    hill(X2,v1,K1,n1), ∅ → X1
    hillR(X1,v2,K2,n2)*hill(X3,v3,K3,n3), ∅ → X2
    hill(X2,v4,K4,n4), ∅ → X3
    (d1,d2,d3), (X1,X2,X3) → ∅
end v1 K1 n1 v2 K2 n2 v3 K3 n3 v4 K4 n4 d1 d2 d3

reaction_networks_hill[9] = @reaction_network begin
    hill(X1,v1,K1,n1)*hillR(X1,v2,K2,n2), ∅ → X1
    d, X1 → ∅
end v1 K1 n1 v2 K2 n2 d

reaction_networks_hill[10] = @reaction_network begin
    hill(X2,v1,K1,n1), ∅ → X1
    hillR(X4,v2,K2,n2), ∅ → X2
    (k1,k2), 2X1 + X2 ⟷ X3
    (k3,k4), 2X2 + X3 ⟷ X4
    (k5,k6), X1 + X2 + X3 + X4 ⟷ X5 + X6
    (d1,d2), (X5,X6) → ∅
end v1 K1 n1 v2 K2 n2 k1 k2 k3 k4 k5 k6 d1 d2


#Reaction networks with fixed concentrations.
#The fixed concentration(s) is declared underneth the reaction network.
reaction_networks_fixed_conc_1 = @reaction_network begin
    (k1,k2), X1 ↔ X2
    (k3,k4), X2 ↔ X3
    (k5,k6), X3 ↔ X1
end k1 k2 k3 k4 k5 k6 FC
@add_constraint reaction_networks_fixed_conc_1 X1+X2+X3=FC
reaction_networks_fixed_conc[1] = reaction_networks_fixed_conc_1


reaction_networks_fixed_conc_2 = @reaction_network begin
    (k1,k2), X1 ↔ 2X1
    (k3,k4), X1 + X2 ↔ X3
    (k5,k6), X3 ↔ X2
end k1 k2 k3 k4 k5 k6 FC
@add_constraint reaction_networks_fixed_conc_2 X2+X3=FC
reaction_networks_fixed_conc[2] = reaction_networks_fixed_conc_2

reaction_networks_fixed_conc_3 = @reaction_network begin
    (k1,k2*X5), X1 ↔ X2
    (k3*X5,k4), X3 ↔ X4
    (p+k5*X2*X3,d), ∅ ↔ X5
end k1 k2 k3 k4 p k5 d FC1 FC2
@add_constraints reaction_networks_fixed_conc_3 begin
    X1+X2=FC1
    X3+X4=FC2
end
reaction_networks_fixed_conc[3] = reaction_networks_fixed_conc_3

reaction_networks_fixed_conc_4 = @reaction_network begin
    (k1,k2), X1 + X2 ↔ X3
    (mm(X3,v,K),d), ∅ ↔ X4
end k1 k2 v K d FC1 FC2
@add_constraints reaction_networks_fixed_conc_4 begin
    X1+X3=FC1
    X2+X3=FC2
end
reaction_networks_fixed_conc[4] = reaction_networks_fixed_conc_4

reaction_networks_fixed_conc_5 = @reaction_network begin
    (k1,k2), X1 ↔ 2X2
    (k3,k4), 2X2 ↔ 3X3
    (k5,k6), 3X3 ↔ 4X4
end k1 k2 k3 k4 k5 k6 FC
@add_constraint reaction_networks_fixed_conc_5 24X1+12X2+8X3+6X4=FC
reaction_networks_fixed_conc[5] = reaction_networks_fixed_conc_5

reaction_networks_fixed_conc_6 = @reaction_network begin
    mmR(X1,v1,K1), X1 → X2
    mmR(X2,v2,K2), X2 → X3
    mmR(X3,v3,K3), X3 → X1
end v1 K1 v2 K2 v3 K3 FC
@add_constraint reaction_networks_fixed_conc_6 X1+X2+X3=FC
reaction_networks_fixed_conc[6] = reaction_networks_fixed_conc_6

reaction_networks_fixed_conc_7 = @reaction_network begin
    (k1,k2), X1 + X2 ↔ X3
    (mm(X3,v,K),d), ∅ ↔ X2
    (k3,k4), X2 ↔ X4
end k1 k2 k3 k4 v K d FC
@add_constraint reaction_networks_fixed_conc_7 X1+X3=FC
reaction_networks_fixed_conc[7] = reaction_networks_fixed_conc_7

reaction_networks_fixed_conc_8 = @reaction_network begin
    (k1,k2), X1 + X2 ↔ X3
    (mm(X3,v1,K1),mm(X4,v2,K2)), X3 ↔ X4
end k1 k2 v1 K1 v2 K2 FC1 FC2
@add_constraints reaction_networks_fixed_conc_8 begin
    X1+X2+X4=FC1
    X2+X3+X4=FC2
end
reaction_networks_fixed_conc[8] = reaction_networks_fixed_conc_8

reaction_networks_fixed_conc_9 = @reaction_network begin
    (k1,k2), X1 + X2 ↔ X3
    (k3,k4), X3 + X4 ↔ X5
    (k5,k6), X5 + X6 ↔ X7
end k1 k2 k3 k4 k5 k6 FC1 FC2 FC3 FC4
@add_constraints reaction_networks_fixed_conc_9 begin
    X1+X3+X5+X7=FC1
    X2+X3+X5+X7=FC2
    X4+X5+X7=FC3
    X6+X7=FC4
end
reaction_networks_fixed_conc[9] = reaction_networks_fixed_conc_9

reaction_networks_fixed_conc_10 = @reaction_network rnType begin
    kDeg,       (w,w2,w2v,v,w2v2,vP,σB,w2σB) ⟶ ∅
    kDeg,       vPp ⟶ phos
    (kBw,kDw),  2w ⟷ w2
    (kB1,kD1),  w2 + v ⟷ w2v
    (kB2,kD2),  w2v + v ⟷ w2v2
    kK1,        w2v ⟶ w2 + vP
    kK2,        w2v2 ⟶ w2v + vP
    (kB3,kD3),  w2 + σB ⟷ w2σB
    (kB4,kD4),  w2σB + v ⟷ w2v + σB
    (kB5,kD5),  vP + phos ⟷ vPp
    kP,         vPp ⟶ v + phos
    v0*((1+F*σB)/(K+σB)),     ∅ ⟶ σB
    λW*v0*((1+F*σB)/(K+σB)),  ∅ ⟶ w
    λV*v0*((1+F*σB)/(K+σB)),  ∅ ⟶ v
end kBw kDw kD kB1 kB2 kB3 kB4 kB5 kD1 kD2 kD3 kD4 kD5 kK1 kK2 kP kDeg v0 F K λW λV pTot;
@add_constraint reaction_networks_fixed_conc_10 vPp+phos=pTot
reaction_networks_fixed_conc[10] = reaction_networks_fixed_conc_10

#Minimal reaction network version of the standard reaction networks.
min_reaction_networks_standard[1] = @min_reaction_network begin
    (p1,p2,p3), ∅ → (X1,X2,X3)
    (k1,k2), X2 ⟷ X1 + 2X3
    (k3,k4), X1 ⟷ X3
    (d1,d2,d3), (X1,X2,X3) → ∅
end p1 p2 p3 k1 k2 k3 k4 d1 d2 d3

min_reaction_networks_standard[2] = @min_reaction_network begin
    mmR(X2,v1,K1), ∅ → X1
    mm(X1,v2,K2), ∅ → X2
    d, X1+X2 → ∅
end v1 K1 v2 K2 d

min_reaction_networks_standard[3] = @min_reaction_network begin
    mm(X2,v1,K1), ∅ → X1
    mm(X3,v2,K2), ∅ → X2
    (k1,k2), X1 ⟷ X4
    (k3,k4), X4 + X2 ⟷ X3 +X1
    d, (X1,X2,X3,X4) → ∅
end v1 K1 v2 K2 k1 k2 k3 k4 d

min_reaction_networks_standard[4] = @min_reaction_network begin
    mmR(X4,v1,K1), ∅ → X1
    mmR(X1,v2,K2), ∅ → X2
    mmR(X2,v3,K3), ∅ → X3
    mmR(X3,v4,K4), ∅ → X4
    (d1,d2,d3,d4), (X1,X2,X3,X4) → ∅
end v1 K1 v2 K2 v3 K3 v4 K4 d1 d2 d3 d4

min_reaction_networks_standard[5] = @min_reaction_network begin
    p, ∅ → X1
    (k1,k2), X1 ⟷ X2
    (k3,k4), X2 ⟷ X3
    (k5,k6), X3 ⟷ X4
    d, X4 → ∅
end p k1 k2 k3 k4 k5 k6 d

min_reaction_networks_standard[6] = @min_reaction_network begin
    (p1,p2), ∅ → (X1,X2)
    (k1,k2), 2X1 ⟷ X3
    (k3,k4), X2 + X3 ⟷ 4X4
    (k5,k6), X4 + X1 ⟷ 2X3
    d, (X1,X2,X3,X4) → ∅
end p1 p2 k1 k2 k3 k4 k5 k6 d

min_reaction_networks_standard[7] = @min_reaction_network begin
    (p1,p2,p3,p4,p5), ∅ → (X1,X2,X3,X4,X5)
    (k1,k2), X1 + 3X2 ⟷ 2X3 + 2X4
    (k3,k4), X2 + X3 ⟷ X5
    (mm(X1,v1,K1),k5), X5 ⟷ X1
    (mm(X2,v2,K2),k6), X1 ⟷ X2
    (d1,d2,d3,d4,d5), (X1,X2,X3,X4,X5) → ∅
end p1 p2 p3 p4 p5 k1 k2 k3 k4 k5 k6 v1 K1 v2 K2 d1 d2 d3 d4 d5

min_reaction_networks_standard[8] = @min_reaction_network begin
    p, ∅ → 2X1
    k1, X1 → X2
    (k2, k3), X2 → X3
    d, X3 → ∅
end p k1 k2 k3 d

min_reaction_networks_standard[9] = @min_reaction_network begin
    (p1,p2,p3), ∅ ⟶ (X1,X2,X3)
    (d1,d2,d3), (X1,X2,X3) ⟶ ∅
    (k1,k2), X1 + X2 ⟷ X3
    (k3,k4), X3 ⟷ X4
    (k5,k6), X4 ⟷ X1 + X2
end p1 p2 p3 d1 d2 d3 k1 k2 k3 k4 k5 k6

min_reaction_networks_standard[10] = @min_reaction_network begin
    p, ∅ ⟶ X1
    (k1,k2), X1 ⟷ X2
    (k3,k4), X2 ⟷ X3
    (k5,k6), X3 ⟷ X4
    (k7,k8), X4 ⟷ X5
    d, X5  ⟶ ∅
end p k1 k2 k3 k4 k5 k6 k7 k8 d
