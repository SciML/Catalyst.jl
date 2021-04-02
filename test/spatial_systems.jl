using MacroTools, Test
include("../src/Catalyst.jl")

G = Catalyst.PeriodicGrid(3)
H = Catalyst.PeriodicGrid(2)

systems = [ 
(quote
    sigma[i], 0 --> X[i]
    D, X[i] --> X[j]
    delta, 2X[i]-->0
end, 
[Meta.parse("PeriodicGrid(3)[i,j,k]")],  
quote
    var"sigma_(1)", 0 → var"X_(1)"
    var"sigma_(2)", 0 → var"X_(2)"
    var"sigma_(3)", 0 → var"X_(3)"
    D, var"X_(3)" → var"X_(1)"
    D, var"X_(1)" → var"X_(2)"
    D, var"X_(2)" → var"X_(3)"
    D, var"X_(2)" → var"X_(1)"
    D, var"X_(3)" → var"X_(2)"
    D, var"X_(1)" → var"X_(3)"
    delta, 2var"X_(1)" → 0
    delta, 2var"X_(2)" → 0
    delta, 2var"X_(3)" → 0
end       
),
(quote
    sigma, 0 --> X[i]
    D, X[i] --> X[j]
    delta, 2X[i]-->0
    kappa[k], 0 --> Y[k]
    iota[i], X[i] + Y[k] --> 0
end, 
[Meta.parse("PeriodicGrid(3)[i,j]"), Meta.parse("PeriodicGrid(2)[k]")],  
quote
    sigma, 0 → var"X_(1)"
    sigma, 0 → var"X_(2)"
    sigma, 0 → var"X_(3)"
    D, var"X_(3)" → var"X_(1)"
    D, var"X_(1)" → var"X_(2)"
    D, var"X_(2)" → var"X_(3)"
    D, var"X_(2)" → var"X_(1)"
    D, var"X_(3)" → var"X_(2)"
    D, var"X_(1)" → var"X_(3)"
    delta, 2var"X_(1)" → 0
    delta, 2var"X_(2)" → 0
    delta, 2var"X_(3)" → 0
    var"kappa_(1)", 0 → var"Y_(1)"
    var"kappa_(2)", 0 → var"Y_(2)"
    var"iota_(1)", var"X_(1)" + var"Y_(1)" → 0
    var"iota_(2)", var"X_(2)" + var"Y_(1)" → 0
    var"iota_(3)", var"X_(3)" + var"Y_(1)" → 0
    var"iota_(1)", var"X_(1)" + var"Y_(2)" → 0
    var"iota_(2)", var"X_(2)" + var"Y_(2)" → 0
    var"iota_(3)", var"X_(3)" + var"Y_(2)" → 0
end   
)
]

for (expr, graphs, target) in systems
    expr = MacroTools.striplines(expr)
    target = MacroTools.striplines(target)
    graphs = MacroTools.striplines.(graphs)
    expr_spatial = Catalyst.make_spatial_reaction_network(expr, graphs)
                    
    @test target == expr_spatial
end