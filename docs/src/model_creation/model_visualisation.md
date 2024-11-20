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
brusselator # hide
```
Here, we note that the output of `latexify(brusselator)` is identical to how a model is displayed by default. Indeed, the reason is that Catalyst internally uses Latexify's `latexify` function to display its models. It is also possible to display the ODE equations a model would generate by adding the `form = :ode` argument:
```@example visualisation_latex
latexify(brusselator; form = :ode)
```
!!! note
    Internally, `latexify(brusselator; form = :ode)` calls `latexify(convert(ODESystem, brusselator))`. Hence, if you have already generated the `ODESystem` corresponding to your model, it can be used directly as input to `latexify`.

!!! note 
    It should be possible to also generate SDEs through the `form = :sde` input. This feature is, however, currently broken.

If you wish to copy the output to your [clipboard](https://en.wikipedia.org/wiki/Clipboard_(computing)) (e.g. so that you can paste it into a LaTeX document), run `copy_to_clipboard(true)` before you run `latexify`. A more throughout description of Latexify's features can be found in [its documentation](https://korsbo.github.io/Latexify.jl/stable/).

!!! note
    For a model to be nicely displayed you have to use an IDE that actually supports this (such as a [notebook](https://jupyter.org/)). Other environments (such as [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/)) will simply return the full LaTeX code which would generate the desired expression. 

## [Displaying model networks](@id visualisation_graphs)
Catalyst uses `GraphMakie` to display representations of chemical reaction networks, including the complex graph and the species-reaction graph (which is similar to the [Petri net](https://en.wikipedia.org/wiki/Petri_net) representation). To get started, import Catalyst and GraphMakie to load the `CatalystGraphMakieExtension` extension, and then load a Makie backend (`GLMakie` is recommended).

```@example visualisation_graphs
using Catalyst, GraphMakie
using GLMakie
```

A network graph showing a Catalyst model's species and reactions can be displayed using the `Graph` function. This first requires [Graphviz](https://graphviz.org/) to be installed and command line accessible. 

Let's declare a [Brusselator model](@ref basic_CRN_library_brusselator) to see this plotting functionality. The functions `plot_network` and `plot_complexes` are used to create the species-reaction and complex graphs, respectively. For a more thorough description of these two representations, please see the [network visualization](@network_visualization) section of the API, but the gist is that the species-reaction graph has species and reactions as nodes, and the complex graph has reaction complexes as nodes. Below we will plot the species-reaction graph using `plot_network`. 
```@example visualisation_graphs
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
plot_network(brusselator)
nothing # hide
```
!["Brusselator Graph"](../assets/network_graphs/brusselator_graph_makie.png)

The species-reaction graph (or network graph) represents species as blue nodes and reactions as green dots. Black arrows from species to reactions indicate substrates, and are labelled with their respective stoichiometries. Similarly, black arrows from reactions to species indicate products (also labelled with their respective stoichiometries). If there are any reactions where a species affect the rate, but does not participate as a reactant, this is displayed with a dashed red arrow. This can be seen in the following [Repressilator model](@ref basic_CRN_library_repressilator):
```@example visualisation_graphs
repressilator = @reaction_network begin
    hillr(Z,v,K,n), ∅ --> X
    hillr(X,v,K,n), ∅ --> Y
    hillr(Y,v,K,n), ∅ --> Z
    d, (X, Y, Z) --> ∅
end
plot_network(repressilator)
nothing # hide
```
!["Repressilator Graph"](../assets/network_graphs/repressilator_graph_makie.png)

A generated graph can be saved using Makie's `save` function. 
```julia
repressilator_graph = plot_network(repressilator)
save("repressilator_graph.png", repressilator_graph)
```

Finally, a [network's reaction complexes](@ref network_analysis_reaction_complexes) (and the reactions in between these) can be displayed using the `plot_complexes(brusselator)` function:
```@example visualisation_graphs
plot_complexes(brusselator)
nothing # hide
```
!["Brusselator Complex Graph"](../assets/network_graphs/brusselator_complex_graph_makie.png)
Here, reaction complexes are displayed as blue nodes, and reactions between complexes are displayed as black arrows. Edges are labeled with their rate expressions.

Makie graphs can be made to be interactive, allowing one to drag nodes and edges. To do this, we retrieve the axis from the GraphMakie plot, and then register the interactions. Note that this can only be done if `GLMakie` is the installed Makie backend. See the [GraphMakie docs](https://graph.makie.org/stable/#Predefined-Interactions) for more information.  

```@example visualization_graphs
using GLMakie
f, ax, p = plot_network(brusselator)
deregister_interaction!(ax, :rectanglezoom)
register_interaction!(ax, :ndrag, NodeDrag(p))
register_interaction!(ax, :edrag, EdgeDrag(p))
nothing
```

The equivalent of `show` for Makie objects is the `display` function. 
```@example visualization_graphs
f = plot_network(brusselator)
display(f)
nothing #hide
```
