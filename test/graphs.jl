using Catalyst, Graphviz_jll 

rn = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
    S^2, R --> 0
end α β

# check can make a graph
gr = Graph(rn)

# check can save a graph
fname = Base.Filesystem.tempname()
savegraph(gr, fname, "png")

