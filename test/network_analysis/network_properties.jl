#! format: off

using Catalyst, Test

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


#########################
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


######################

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
