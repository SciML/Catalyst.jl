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

# Tests that the species-reaction graph has a vertex for each species and reaction, and that its
# edges are the sparsity pattern of the substrate and product stoichiometry matrices (which
# exclude constant species).
let
    const_rns = [
        @reaction_network(begin
            @parameters X [isconstantspecies = true]
            k, X + Y --> XY
        end),
        @reaction_network(begin
            @parameters X [isconstantspecies = true]
            k1, 2X + Y --> X + Z
            k2, Z --> 3X + Y
        end),
    ]
    for rn in [reaction_networks_standard; reaction_networks_hill;
            reaction_networks_conserved; reaction_networks_real; const_rns]
        srg = Catalyst.species_reaction_graph(rn)
        ns, nr = numspecies(rn), numreactions(rn)
        S, P = substoichmat(rn), prodstoichmat(rn)
        @test nv(srg) == ns + nr
        @test Set(edges(srg)) ==
            Set([[Graphs.Edge(i, ns + j) for i in 1:ns, j in 1:nr if !iszero(S[i, j])];
                [Graphs.Edge(ns + j, i) for i in 1:ns, j in 1:nr if !iszero(P[i, j])]])
    end

    srg = Catalyst.species_reaction_graph(const_rns[1])
    @test nv(srg) == 3
    @test Set(edges(srg)) == Set([Graphs.Edge(1, 3), Graphs.Edge(3, 2)])
end

# Returns the edges drawn by `plot_network` as (source label, destination label, colour, edge
# label) tuples. These are sorted, so that comparisons also check how often each edge is drawn.
function sr_edge_table(rn)
    f, ax, p = plot_network(rn)
    g, lbl = p[1][], p.ilabels[]
    sort!([(lbl[src(e)], lbl[dst(e)], p.edge_color[][i], p.elabels[][i])
        for (i, e) in enumerate(edges(g))])
end

# Tests that constant species are drawn by `plot_network` as grey nodes (after the reaction
# nodes), and that a constant substrate is drawn the same as a constant species in the rate.
let
    rn1 = @reaction_network begin
        @parameters X [isconstantspecies = true]
        k, X + Y --> XY
    end
    rn2 = @reaction_network begin
        @parameters X [isconstantspecies = true]
        k * X, Y --> XY
    end
    for rn in (rn1, rn2)
        @test sr_edge_table(rn) == sort!([("X", "R1", :red, ""), ("Y", "R1", :black, "1"),
            ("R1", "XY", :black, "1")])
        f, ax, p = plot_network(rn)
        @test p.ilabels[] == ["Y", "XY", "R1", "X"]
        @test p.node_color[] == [:skyblue3, :skyblue3, :green, :grey]
    end
end

# Tests that constant substrates are labelled with non-unit stoichiometries, and that a constant
# species that is both a substrate and in the rate is drawn with a single edge.
let
    rn = @reaction_network begin
        @parameters X [isconstantspecies = true]
        k, 2X + Y --> Z
    end
    @test sr_edge_table(rn) == sort!([("X", "R1", :red, "2"), ("Y", "R1", :black, "1"),
        ("R1", "Z", :black, "1")])

    rn = @reaction_network begin
        @parameters X [isconstantspecies = true]
        k * X, X + Y --> Z
    end
    @test sr_edge_table(rn) == sort!([("X", "R1", :red, ""), ("Y", "R1", :black, "1"),
        ("R1", "Z", :black, "1")])
end

# Tests that constant products, and constant substrates of `only_use_rate` reactions that do not
# appear in the rate, are drawn with grey edges (as they do not affect the dynamics).
let
    rn = @reaction_network begin
        @parameters X [isconstantspecies = true]
        k, Y --> 2X + Z
    end
    @test sr_edge_table(rn) == sort!([("Y", "R1", :black, "1"), ("R1", "X", :grey, "2"),
        ("R1", "Z", :black, "1")])

    rn = @reaction_network begin
        @parameters X [isconstantspecies = true] W [isconstantspecies = true]
        k * W, X + W + Y => Z
    end
    @test sr_edge_table(rn) == sort!([("X", "R1", :grey, ""), ("W", "R1", :red, ""),
        ("Y", "R1", :black, "1"), ("R1", "Z", :black, "1")])
end

# Tests that a constant species appearing in several reactions is drawn as a single node, and
# that constant array species are labelled correctly.
let
    rn = @reaction_network begin
        @parameters E [isconstantspecies = true]
        kB, S + E --> SE
        kP, SE --> P + E
    end
    @test sr_edge_table(rn) == sort!([("S", "R1", :black, "1"), ("E", "R1", :red, ""),
        ("R1", "SE", :black, "1"), ("SE", "R2", :black, "1"), ("R2", "P", :black, "1"),
        ("R2", "E", :grey, "")])
    f, ax, p = plot_network(rn)
    @test p.ilabels[] == ["S", "SE", "P", "R1", "R2", "E"]
    @test p.node_color[] == [:skyblue3, :skyblue3, :skyblue3, :green, :green, :grey]

    rn = @reaction_network begin
        @parameters x[1:2] [isconstantspecies = true]
        @species (X(t))[1:2]
        k, X[1] + x[1] --> X[2]
    end
    @test sr_edge_table(rn) == sort!([("X[1]", "R1", :black, "1"), ("x[1]", "R1", :red, ""),
        ("R1", "X[2]", :black, "1")])
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
    srg = first(CGME.SRGraphWrap(rn))
    s = length(species(rn))
    @test ne(srg) == 8
    @test Graphs.Edge(2, s+3) ∈ srg.multiedges
    # Since B is both a dep and a reactant
    @test count(==(Graphs.Edge(2, s+3)), edges(srg)) == 2
    @test sr_edge_table(rn) == sort!([("A", "R1", :black, "1"), ("R1", "B", :black, "1"),
        ("A", "R2", :black, "1"), ("C", "R2", :red, ""), ("R2", "C", :black, "1"),
        ("B", "R3", :black, "1"), ("B", "R3", :red, ""), ("R3", "C", :black, "1")])

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
    srg = first(CGME.SRGraphWrap(rn))
    s = length(species(rn))
    @test ne(srg) == 8
    # Since A, B is both a dep and a reactant
    @test count(==(Graphs.Edge(1, s+2)), edges(srg)) == 2
    @test count(==(Graphs.Edge(2, s+3)), edges(srg)) == 2
    @test sr_edge_table(rn) == sort!([("A", "R1", :black, "1"), ("R1", "B", :black, "1"),
        ("A", "R2", :black, "1"), ("A", "R2", :red, ""), ("R2", "C", :black, "1"),
        ("B", "R3", :black, "1"), ("B", "R3", :red, ""), ("R3", "C", :black, "1")])
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
    @test rxorder == [1, 4, 6, 2, 5, 3]

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
