using Catalyst, GraphMakie, CairoMakie, Graphs, SparseArrays
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
    edgel = Graphs.Edge.([(s+1, 1),
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

    # Test that figures are generated properly.
    f = plot_network(MAPK)
    save("fig.png", f)
    @test isfile("fig.png")
    rm("fig.png")
    f = plot_network(brusselator)
    save("fig.png", f)
    @test isfile("fig.png")
    rm("fig.png")

    f = plot_complexes(MAPK); save("fig.png", f)
    @test isfile("fig.png")
    rm("fig.png")
    f = plot_complexes(brusselator); save("fig.png", f)
    @test isfile("fig.png")
    rm("fig.png")
end

CGME = Base.get_extension(parentmodule(ReactionSystem), :CatalystGraphMakieExtension)
# Test that rate edges are inferred correctly. We should see two for the following reaction network.
let
    # Two rate edges, one to species and one to product
    rn = @reaction_network begin
        k, A --> B
        k * C, A --> C
        k * B, B --> C
    end
    srg = CGME.SRGraphWrap(rn)
    s = length(species(rn))
    @test ne(srg) == 8
    @test Graphs.Edge(2, s+3) ∈ srg.multiedges
    # Since B is both a dep and a reactant
    @test count(==(Graphs.Edge(2, s+3)), edges(srg)) == 2

    f = plot_network(rn)
    save("fig.png", f)
    @test isfile("fig.png")
    rm("fig.png")
    f = plot_complexes(rn); save("fig.png", f)
    @test isfile("fig.png")
    rm("fig.png")

    # Two rate edges, both to reactants
    rn = @reaction_network begin
        k, A --> B
        k * A, A --> C
        k * B, B --> C
    end
    srg = CGME.SRGraphWrap(rn)
    s = length(species(rn))
    @test ne(srg) == 8
    # Since A, B is both a dep and a reactant
    @test count(==(Graphs.Edge(1, s+2)), edges(srg)) == 2
    @test count(==(Graphs.Edge(2, s+3)), edges(srg)) == 2
end

function test_edgeorder(rn)
    # The initial edgelabels in `plot_complexes` is given by the order of reactions in reactions(rn).
    D = incidencemat(rn; sparse=true)
    rxs = reactions(rn)
    edgelist = Vector{Graphs.SimpleEdge{Int}}()
    rows = rowvals(D)
    vals = nonzeros(D)

    for (i, rx) in enumerate(rxs)
        inds = nzrange(D, i)
        val = vals[inds]
        row = rows[inds]
        (sub, prod) = val[1] == -1 ? (row[1], row[2]) : (row[2], row[1])
        push!(edgelist, Graphs.SimpleEdge(sub, prod))
    end

    img, rxorder = CGME.ComplexGraphWrap(rn)

    # Label iteration order is given by edgelist[rxorder]. Actual edge drawing iteration order is given by edges(g)
    @test edgelist[rxorder] == Graphs.edges(img)
    return rxorder
end

# Test edge order for complexes.
let
    # Multiple edges
    rn = @reaction_network begin
        k1, A --> B
        (k2, k3), C <--> D
        k4, A --> B
    end
    rxorder = test_edgeorder(rn)
    edgelabels = [repr(rx.rate) for rx in reactions(rn)]
    # Test internal order of labels is preserved
    @test edgelabels[rxorder][1] == "k1"
    @test edgelabels[rxorder][2] == "k4"

    # Multiple edges with species dependencies
    rn = @reaction_network begin
        k1, A --> B
        (k2, k3), C <--> D
        k4, A --> B
        hillr(D, α, K, n), C --> D
        k5*B, A --> B
    end
    rxorder = test_edgeorder(rn)
    edgelabels = [repr(rx.rate) for rx in reactions(rn)]
    labels = ["k1", "k4", "k5*B(t)", "k2", "Catalyst.hillr(D(t), α, K, n)", "k3"]
    @test edgelabels[rxorder] == labels

    rs = @reaction_network begin
        ka, Depot --> Central
        (k12, k21), Central <--> Peripheral
        ke, Central --> 0
    end
    test_edgeorder(rs)

    rn = @reaction_network begin
        (k1, k2), A <--> B
        k3, C --> B
        (α, β), (A, B) --> C
        k4, B --> A
        (k5, k6), B <--> A
        k7, B --> C
        (k8, k9), C <--> A
        (k10, k11), (A, C) --> B
        (k12, k13), (C, B) --> A
    end
    rxorder = test_edgeorder(rn)
    edgelabels = [repr(rx.rate) for rx in reactions(rn)]
    @test edgelabels[rxorder][1:3] == ["k1", "k6", "k10"]
end

# Test that array species are labeled correctly in plots.
let
    rn = @reaction_network begin
        @species (S(t))[1:3]
        k1, S[1] --> S[2]
        k2, S[2] --> S[3]
    end
    f = plot_network(rn)
    save("fig.png", f)
    @test isfile("fig.png")
    rm("fig.png")
    f = plot_complexes(rn)
    save("fig.png", f)
    @test isfile("fig.png")
    rm("fig.png")
end
