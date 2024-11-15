using Catalyst, GraphMakie, GLMakie
include("../test_networks.jl")
# Test that speciesreactiongraph is generated correctly
let
    brusselator = @reaction_network begin
        A, ∅ --> X
        1, 2X + Y --> 3X
        B, X --> Y
        1, X --> ∅
    end

    srg = Catalyst.species_reaction_graph(brusselator)
    s = length(species(brusselator))
    edgel = Edge.([(s+1, 1),
                   (1, s+2),
                   (2, s+2),
                   (s+2, 1),
                   (s+3, 2),
                   (1, s+3),
                   (1, s+4)])
    @test all(∈(collect(Graphs.edges(srg))), edgel)

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
    srg = Catalyst.species_reaction_graph(MAPK)
    @test nv(srg) == length(species(MAPK)) + length(reactions(MAPK))
    @test ne(srg) == 90
end

# Test that rate edges are inferred correctly
let
    rn = @reaction_network begin
        k, A --> B
        k * C, A --> C
        k * B, B --> C
    end
    srg = SRGraphWrap(rn)
    s = length(species(rn))
    @test Edge(3, s+2) ∈ srg.rateedges 
    @test Edge(2, s+3) ∈ srg.rateedges 
    # Since B is both a dep and a reactant
    @test count(==(Edge(2, s+3)), edges(srg)) == 2
end
