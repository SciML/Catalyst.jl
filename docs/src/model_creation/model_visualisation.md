# [Model Visualisation](@id visualisation)
Catalyst-created `ReactionSystem` models can be visualised either as LaTeX code (of either the model reactions or its equations) or as a network graph. This section describes both functionalities.

## [Displaying models using LaTeX](@id visualisation_latex)
Once a model has been created, the [Latexify.jl](https://github.com/korsbo/Latexify.jl) package can be used to generate LaTeX code of the model. This can either be used for easy model inspection (e.g. to check which equations are being simulated), or to generate code which can be directly pasted into a LaTeX document.

Let us consider a simple [Brusselator model](@ref basic_CRN_library_brusselator):
```@example visualisation_latex
using Catalyst
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
```
To display its reaction (using LaTeX formatting) we run `latexify` with our model as input:
```@example visualisation_latex
using Latexify
latexify(brusselator)
```
Here, we note that the output of `latexify(brusselator)` is identical to how a model is displayed by default. Indeed, the reason is that Catalyst internally uses Latexify's `latexify` function to display its models. It is also possible to display the ODE equations a model would generate by adding the `form = :ode` argument:
```@example visualisation_latex
latexify(brusselator; form = :ode)
```
!!! note
    Internally, `latexify(brusselator; form = :ode)` calls `latexify(convert(ODESystem, brusselator))`. Hence, if you have already [generated the `ODESystem` corresponding to your model](@ref ref), it can be used directly as input to `latexify`.

!!! note 
    It should be possible to also generate SDEs through the `form = :sde` input. This feature is, however, currently broken.

If you wish to copy the output to your [clipboard](https://en.wikipedia.org/wiki/Clipboard_(computing)) (e.g. so that you can paste it into a LaTeX document), run `copy_to_clipboard(true)` before you run `latexify`. A more throughout description of Latexify's features can be found in [its documentation](https://korsbo.github.io/Latexify.jl/stable/).

!!! note
    For a model to be nicely displayed you have to use an IDE that actually supports this (such as a [notebook](https://jupyter.org/)). Other environments (such as [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/)) will simply return the full LaTeX code which would generate the desired expression. 

## [Displaying model networks](@id visualisation_graphs)
A network graph showing a Catalyst model's species and reactions can be displayed using the `Graph` function. This first requires [Graphviz](https://graphviz.org/) to be installed and command line accessible. Here, we first declare a [Brusselator model](@ref basic_CRN_library_brusselator) and then displays its network topology:
```@example visualisation_graphs
using Catalyst
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
Graph(brusselator)
nothing # hide
```
!["Brusselator Graph"](../assets/network_graphs/brusselator_graph.png)

The network graph represents species as blue nodes and reactions as orange dots. Black arrows from species to reactions indicate substrates, and are labelled with their respective stoichiometries. Similarly, black arrows from reactions to species indicate products (also labelled with their respective stoichiometries). If there are any reactions where a species affect the rate, but does not participate as a reactant, this is displayed with a dashed red arrow. This can be seen in the following [Repressilator model](@ref basic_CRN_library_repressilator):
```@example visualisation_graphs
repressilator = @reaction_network begin
    hillr(Z,v,K,n), ∅ --> X
    hillr(X,v,K,n), ∅ --> Y
    hillr(Y,v,K,n), ∅ --> Z
    d, (X, Y, Z) --> ∅
end
Graph(repressilator)
nothing # hide
```
!["Repressilator Graph"](../assets/network_graphs/repressilator_graph.png)

A generated graph can be saved using the `savegraph` function:
```julia
repressilator_graph = Graph(repressilator)
savegraph(repressilator_graph, "repressilator_graph.png")
```

Finally, a [network's reaction complexes](@ref network_analysis_reaction_complexes) (and the reactions in between these) can be displayed using the `complexgraph(brusselator)` function:
```@example visualisation_graphs
complexgraph(brusselator)
nothing # hide
```
!["Repressilator Complex Graph"](../assets/network_graphs/repressilator_complex_graph.png)
Here, reaction complexes are displayed as blue nodes, and reactions in between these as black arrows.
