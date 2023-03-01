using Catalyst, Graphviz_jll

### Basic Tests ###
#let
    rn = @reaction_network begin
        α, S + I --> 2I
        β, I --> R
        S^2, R --> 0
    end

    # Check can make a graph.
    gr = Graph(rn)

    # Check can save a graph.
    fname = Base.Filesystem.tempname()
    savegraph(gr, "$fname.svg", "svg")

    rcgr = complexgraph(rn)
    fname = Base.Filesystem.tempname()
    savegraph(rcgr, "$fname.svg", "svg")
#end

# these are broken in the jll, see
# https://github.com/JuliaPackaging/Yggdrasil/issues/1428
# savegraph(gr, "$fname.pdf", "pdf")
# savegraph(gr, "$fname.png", "png")
