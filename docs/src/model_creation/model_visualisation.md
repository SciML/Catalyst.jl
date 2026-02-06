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
    Internally, `latexify(brusselator; form = :ode)` calls `latexify(ode_model(brusselator))`. Hence, if you have already generated the `ODESystem` corresponding to your model, it can be used directly as input to `latexify`.

!!! note 
    It should be possible to also generate SDEs through the `form = :sde` input. This feature is, however, currently broken.

If you wish to copy the output to your [clipboard](https://en.wikipedia.org/wiki/Clipboard_(computing)) (e.g. so that you can paste it into a LaTeX document), run `copy_to_clipboard(true)` before you run `latexify`. A more throughout description of Latexify's features can be found in [its documentation](https://korsbo.github.io/Latexify.jl/stable/).

!!! note
    For a model to be nicely displayed you have to use an IDE that actually supports this (such as a [notebook](https://jupyter.org/)). Other environments (such as [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/)) will simply return the full LaTeX code which would generate the desired expression. 

## [Displaying model networks](@id visualisation_graphs)
Catalyst uses [GraphMakie.jl](https://github.com/MakieOrg/GraphMakie.jl) to display representations of chemical reaction networks, including the complex graph and the species-reaction graph (which is similar to the [Petri net](https://en.wikipedia.org/wiki/Petri_net) representation). To get started, import Catalyst, GraphMakie, and NetworkLayout to load the `CatalystGraphMakieExtension` extension, and then load a Makie backend ([`CairoMakie`](https://github.com/MakieOrg/Makie.jl) is a good lightweight choice). Here, while Catalyst primarily uses [Plots.jl](https://github.com/JuliaPlots/Plots.jl) for plotting, [Makie](https://github.com/MakieOrg/Makie.jl) is used for displaying network graphs. Makie can also be used for plotting more generally (and is also a preferred option for some).

```@example visualisation_graphs
using Catalyst, GraphMakie, NetworkLayout
using CairoMakie
nothing # hide
```

Let's declare a [Brusselator model](@ref basic_CRN_library_brusselator) to see this plotting functionality. The functions `plot_network` and `plot_complexes` are used to create the species-reaction and complex graphs, respectively. For a more thorough description of these two representations, please see the [network visualization](@ref network_visualization) section of the API, but the gist is that the species-reaction graph has species and reactions as nodes, and the complex graph has reaction complexes as nodes. Below we will plot the species-reaction graph using `plot_network`. 
```@example visualisation_graphs
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
plot_network(brusselator)
```

The species-reaction graph (or network graph) represents species as blue nodes and reactions as green dots. Black arrows from species to reactions indicate substrates, and are labelled with their respective stoichiometries. Similarly, black arrows from reactions to species indicate products (also labelled with their respective stoichiometries). If there are any reactions where a species affect the rate, but does not participate as a reactant, this is displayed with a dashed red arrow. This can be seen in the following [Repressilator model](@ref basic_CRN_library_repressilator):
```@example visualisation_graphs
repressilator = @reaction_network begin
    hillr(Z,v,K,n), ∅ --> X
    hillr(X,v,K,n), ∅ --> Y
    hillr(Y,v,K,n), ∅ --> Z
    d, (X, Y, Z) --> ∅
end
plot_network(repressilator)
```

A generated graph can be saved using Makie's `save` function. 
```julia
repressilator_graph = plot_network(repressilator)
save("repressilator_graph.png", repressilator_graph)
```

Finally, a [network's reaction complexes](@ref network_analysis_reaction_complexes) (and the reactions in between these) can be displayed using the `plot_complexes(brusselator)` function:
```@example visualisation_graphs
plot_complexes(brusselator)
```
Here, reaction complexes are displayed as blue nodes, and reactions between complexes are displayed as black arrows. Red arrows indicate that the rate constantof a reaction has a species-dependence. Edges can be optionally labeled with their rate expressions by calling with the option `show_rate_labels`.
```@example visualisation_graphs
plot_complexes(brusselator, show_rate_labels = true)
```

## Customizing Plots
In this section we demonstrate some of the ways that plot objects can be manipulated to give nicer images. Let's start with our brusselator plot once again. Note that the `plot` function returns three objects: the `Figure`, the `Axis`, and the `Plot`, which can each be customized independently. See the general [Makie documentation](https://docs.makie.org/stable/) for more information.

```@example visualisation_graphs
f, ax, p = plot_complexes(brusselator, show_rate_labels = true)
```

It seems like a bit of the top node is cut off. Let's increase the top and bottom margins by increasing `yautolimitmargin`.
```@example visualisation_graphs
ax.yautolimitmargin = (0.3, 0.3) # defaults to (0.15, 0.15)
ax.aspect = DataAspect()
f
```

There are many keyword arguments that can be passed to `plot_network` or `plot_complexes` to change the look of the graph (which get passed to the `graphplot` Makie recipe). Let's change the color of the nodes and make the inner labels a bit smaller. Let's also give the plot a title. 
```@example visualisation_graphs
f, ax, p = plot_complexes(brusselator, show_rate_labels = true, node_color = :yellow, ilabels_fontsize = 10)
ax.title = "Brusselator"
f
```

Most of the kwargs that modify the nodes or edges will also accept a vector with the same length as the number of nodes or edges, respectively. See [here](https://graph.makie.org/stable/#The-graphplot-Recipe) for a full list of keyword arguments to `graph_plot`. Note that `plot_complexes` and `plot_network` default to `layout = Stress()` rather than `layout = Spring()`, since `Stress()` is better at generating plots with fewer edge crossings. More layout options and customizations (such as pinning nodes to certain positions) can be found in the [`NetworkLayout` documentation](https://juliagraphs.org/NetworkLayout.jl/stable/).

Once a graph is already created we can also change the keyword arguments by modifying the fields of the `Plot` object `p`.
```@example visualisation_graphs
p.node_color = :orange
f
```

Custom node positions can also be given, if the automatic layout is unsatisfactory.
```@example visualisation_graphs
fixedlayout = [(0,0), (1,0), (0,1), (1,1), (2,0)]
p.layout = fixedlayout
autolimits!(ax)
f
```

Makie graph plots can be made to be interactive, allowing one to drag nodes and edges. To do this, we retrieve the axis from the GraphMakie plot, and then register the interactions. **Note that this can only be done if `GLMakie` is the installed Makie backend. See the [GraphMakie docs](https://graph.makie.org/stable/#Predefined-Interactions) for more information about the types of interactions one can register.** Below is a non-interactive code example that shows how to do this:

```julia
using GLMakie
f, ax, p = plot_network(brusselator)
deregister_interaction!(ax, :rectanglezoom)
register_interaction!(ax, :ndrag, NodeDrag(p))
register_interaction!(ax, :edrag, EdgeDrag(p))
```

The equivalent of `show` for Makie plots is the `display` function. 
```julia
f = plot_network(brusselator)
display(f)
```

Once you are happy with the graph plot, you can save it using the `save` function. 
```julia
save("fig.png", f)
```
