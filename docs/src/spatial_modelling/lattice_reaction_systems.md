# [Introduction to Spatial Modelling with Catalyst](@id spatial_lattice_modelling_intro)
Catalyst supports the expansion of non-spatial [`ReactionSystem`](@ref)s (created using e.g. the `@reaction_network` DSL) to spatial domains. Spatial simulation of Catalyst models is a work in progress. Currently, the following is supported:
- Spatial ODE and Jump simulations.
- Discrete spatial domains.
- Constant-rate transportation reactions (species moving spatially at constant rates). 

Features for which support is planned in future updates include:
- Models on continuous domains with automatic discretisation (these models can already be simulated if the user provides a discretisation).
- SDE simulations.
- Transport reactions with non-constant rates as well as more general spatial reactions.

This tutorial introduces spatial modelling on discrete domains, here called lattices. It describes the basics of creating and simulating such models. To do so, it uses ODE simulations as examples. Additional tutorials provide further details on how to interact with [spatial simulation structures](@ref lattice_simulation_plotting) and [plot spatial simulations](@ref lattice_simulation_plotting), and also provide further details on [ODE](@ref spatial_lattice_ode_simulations) and [jump](@ref spatial_lattice_jump_simulations) simulations, respectively.


## [Basic example of a spatial simulation on a discrete domain](@id spatial_lattice_modelling_intro_example)
To perform discrete-space spatial simulations, the user must first define a [`LatticeReactionSystem`](@ref). These combine:
- A (non-spatial) `ReactionSystem`(@ref) model (created using standard Catalyst syntax).
- A vector of spatial reactions, describing how species can move spatially across the domain.
- A lattice defining the spatial domain's compartments and how they are connected.

Here, as an example, we will simulate a spatial [two-state model](@ref basic_CRN_library_two_states). To do so, we first define our (non-spatial) model, the spatial reactions, and the lattice. These are then bundled into a [`LatticeReactionSystem`](@ref).
```@example spatial_intro_basics
using Catalyst
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
diffusion_rx = @transport_reaction D X1
lattice = CartesianGrid((5,5))
lrs = LatticeReactionSystem(two_state_model, [diffusion_rx], lattice)
nothing # hide
```
This model contains:
- A single spatial reaction, a transport reaction where $X1$ moves at constant rate $D$ between adjacent compartments.
- A 2d Cartesian grid of 5x5 compartments to simulate our model on.

More details on spatial reactions are available [here](@ref spatial_lattice_modelling_intro_spatial_reactions). In addition to Cartesian grid lattices (in 1, 2, and 3 dimensions), masked and unstructured (graph) lattices are also supported. The different lattice types described in more detail [here](@ref spatial_lattice_modelling_intro_lattices).

Once created, [`LatticeReactionSystem`](@ref)s can be used as input to various problem types, which then can be simulated using the same syntax as non-spatial models. Here, we prepare an ODE simulation by creating an `ODEProblem`:
```@example spatial_intro_basics
u0 = [:X1 => rand(5, 5), :X2 => 2.0]
tspan = (0.0, 10.0)
ps = [:k1 => 1.0, :k2 => 2.0, :D => 0.2]
oprob = ODEProblem(lrs, u0, tspan, ps)
nothing # hide
```
In this example we used non-uniform values for $X1$'s initial condition, but uniform values for the remaining initial condition and parameter values. More details of uniform and non-uniform initial conditions and parameter values are provided [here](@ref spatial_lattice_modelling_intro_simulation_inputs). We also note that the diffusion reaction introduces a new parameter, $D$ (determining $X1$'s diffusion rate), whose value must be designated in the parameter vector. 

We can now simulate our model:
```@example spatial_intro_basics
using OrdinaryDiffEq
sol = solve(oprob)
nothing # hide
```
We note that simulations of spatial models are often computationally expensive. Advice on the performance of spatial ODE simulations is provided [here](@ref spatial_lattice_ode_simulations_solvers). 

Finally, we can access "$X1$'s value across the simulation using
```@example spatial_intro_basics
lat_getu(sol, :X1, lrs)
```
and plot the simulation using
```@example spatial_intro_basics
import CairoMakie
lattice_animation(sol, :X1, lrs, "lattice_simulation_2d.mp4")
```
![](./lattice_simulation_2d.mp4)
More information on how to retrieve values from spatial simulations can be found [here](@ref lattice_simulation_structure_interaction_simulation_species), and for plotting them, [here](@ref lattice_simulation_plotting). Finally, a list of functions for querying `LatticeReactionSystems` for various properties can be found [here](@ref api_lattice_simulations).

## [Spatial reactions](@id spatial_lattice_modelling_intro_spatial_reactions)
Spatial reactions describe reaction events which involve species across two connected compartments. Currently, only so-called *transportation reactions* are supported. These consist of:
- A rate at which the reaction occurs. As for non-spatial reactions, this can be any expression. However, currently, it may only consist of parameters and other constants. 
- A single species which is transported from one compartment to an adjacent one.

At the occurrence of a transport reaction, the specific species moves to the adjacent compartment. Many common spatial models can be represented using transport reactions only. These can model phenomena such as diffusion or constant flux. A transportation reaction can be created using the `@transportation_reaction` macro. E.g. above we used
```@example spatial_intro_spat_rxs
using Catalyst # hide
diffusion_rx = @transport_reaction D X1
nothing # hide
```
to create a reaction where species $X$ moves at a constant rate $D$ between adjacent compartments (in the ODE this creates terms $D\cdot X1_i$, where $X1_i$ is the concentration of $X1$ in compartment $i$). Transport reactions may have rates depending on several parameters. E.g. to model a system with two species $X1$ and $X2$, where both species are transported at a rate which depends both on the species, but also on some non-uniform parameter which is unique to each connection (e.g. representing the area connecting two cells in a tissue) we could do:
```@example spatial_intro_spat_rxs
dr_X1 = @transport_reaction D1*a X1
dr_X2 = @transport_reaction D2*a X2
nothing # hide
```

!!! note
    Any species which occurs is occurs in a transport reaction that is used to construct a [`LatticeReactionSystem`](@ref) must also occur in the corresponding non-spatial [`ReactionSystem`](@ref).

### [Creating transport reactions programmatically](@id spatial_lattice_modelling_intro_spatial_reactions_programmatic)
If models are created [programmatically](@ref programmatic_CRN_construction) it is also possible to create transportation reactions programmatically. To do so, use the `TransportReaction` constructor, providing first the rate and then the transported species:
```@example spatial_intro_spat_rxs
@variables t
@species X1(t) X2(t)
@parameters k1 k2 D [edgeparameter=true]
tr_X = TransportReaction(D, X1)
nothing # hide
```
Note that in this example, we specifically designate $D$ as an [edge parameter](@ref spatial_lattice_modelling_intro_simulation_edge_parameters).

## [Defining discrete spatial domains (lattices)](@id spatial_lattice_modelling_intro_lattices)
Discrete spatial domains can represent:
1. Systems which are composed of a (finite number of) compartments, where each compartment can be considered well-mixed (e.g. can be modelled non-spatially) and where (potentially) species can move between adjacent compartments. Tissues, where each compartment corresponds to a biological cell, are examples of such systems.
2. Systems that are continuous in nature, but have been approximated as a discrete domain. Future Catalyst updates will include the ability for the definition, and automatic discretisation, of continuous domains. Currently, however, the user has to perform this discretisation themselves.

Catalyst supports three distinct types of lattices:
- [Cartesian lattices](@ref spatial_lattice_modelling_intro_lattices_cartesian). These are grids where each grid point corresponds to a compartment. Spatial transportation is permitted between adjacent compartments.
- [Masked lattices](@ref spatial_lattice_modelling_intro_lattices_masked). In these grids, only a subset of the grid point actually corresponds to compartments. Spatial transportation is permitted between adjacent compartments.
- [Unstructured (or graph) lattices](@ref spatial_lattice_modelling_intro_lattices_graph). These are defined by graphs, where vertices correspond to compartments and edges connect adjacent compartments.

Here, Cartesian lattices are a subset of the masked lattices, which are a subset of the unstructured lattices. If possible, it is advantageous to use as narrow a lattice definition as possible (this may both improve simulation performance and simplify syntax). Cartesian and masked lattices can be defined as one, two, and three-dimensional. By default, these lattices assume that diagonally neighbouring compartments are non-adjacent (do not permit direct movement of species in between themselves). To change this, provide the `diagonally_adjacent = true` argument to your [`LatticeReactionSystem`](@ref) when it is created.

### [Defining Cartesian lattices](@id spatial_lattice_modelling_intro_lattices_cartesian)
A Cartesian lattice is defined using the `CartesianGrid` function, which takes a single argument. For a 1d grid, simply provide the length of the grid as a single argument:
```@example spatial_intro_lattices
using Catalyst # hide
cgrid_1d = CartesianGrid(5)
nothing # hide
```
For 2d and 3d grids, we instead provide a Tuple with the length of the grid in each dimension:
```@example spatial_intro_lattices
cgrid_2d = CartesianGrid((3, 9))
cgrid_3d = CartesianGrid((2, 4, 8))
nothing # hide
```

### [Defining masked lattices](@id spatial_lattice_modelling_intro_lattices_masked)
Masked lattices are defined through 1d, 2d, or 3d Boolean arrays. Each position in the array is `true` if it corresponds to a compartment, and `false` if it does not. E.g. to define a 1d grid corresponding to two disjoint sets, each with 3 compartments, use:
```@example spatial_intro_lattices
rgrid_1d = [true, true, true, false, true, true, true]
nothing # hide
```
To define a 2d grid corresponding to the shape of an (laying) "8", we can use:
```@example spatial_intro_lattices
rgrid_2d = [
    true  true  true  true  true;
    true  false true  false true;
    true  true  true  true  true
]
nothing # hide
```
Finally, a 4x5x6 3d grid of randomly distributed compartments can be created using:
```@example spatial_intro_lattices
rgrid_3d = rand([true, false], 4, 5, 6)
nothing # hide
```

### [Defining unstructured lattices](@id spatial_lattice_modelling_intro_lattices_graph)
To define unstructured lattices, we must first import the [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) package. Next, we can either use some [pre-defined formula for building graphs](https://juliagraphs.org/Graphs.jl/stable/core_functions/simplegraphs_generators/#Generators-for-common-graphs), or [build a graph from scratch](https://juliagraphs.org/Graphs.jl/stable/first_steps/construction/). Here we create a cyclic graph (where each compartment is connected to exactly two other compartments):
```@example spatial_intro_lattices
using Graphs
cycle_graph(7)
nothing # hide
```

Since graphs can represent any network of connected compartments, they do not have dimensions (like Cartesian or masked lattices). Another feature of graph lattices is that they can have non-symmetric connections, i.e. pairs of compartments where spatial movement of species is only permitted in one direction (in practice, this can be done for Cartesian and masked lattices as well, by [defining non-uniform spatial rates](@ref spatial_lattice_modelling_intro_simulation_inputs) and setting them to zero in one direction). This can be done by using a [*directed graph*](https://juliagraphs.org/Graphs.jl/dev/algorithms/digraph/) as input. E.g. here we define a directed cyclic graph, where movement is only allowed in one direction of the cycle:
```@example spatial_intro_lattices
cycle_digraph(7)
nothing # hide
```

## [Non-uniform initial conditions and parameter values](@id spatial_lattice_modelling_intro_simulation_inputs)
For spatial models, initial conditions and parameter values are provided similarly as for non-spatial models. Wherever a single value is provided, it is used *uniformly* across the lattice. E.g. if we, for our [previous two-state model](@ref spatial_lattice_modelling_intro_example), set
```@example spatial_intro_nonuniform_vals
u0 = [:X1 => 1.0, :X2 => 2.0]
ps = [:k1 => 1.0, :k2 => 2.0, :D => 0.2]
nothing # hide
```
The initial condition will be $1.0$ for $X1$ across compartments, and $2.0$ for $X2$. Furthermore, for each simulation, in each compartment, the value of $1.0$ will be used for $k1$ and $2.0$ for $k2$. Finally, the transportation rate of $X1$ (set by the parameter $D$) will be $0.2$ across all connections. To set non-uniform values, non-scalar values must be provided in the map. How to do this depends on how the [lattice was defined](@ref spatial_lattice_modelling_intro_lattices). Furthermore, some parameters that are part of spatial reactions may have their value tied to *connections between compartments*, rather than compartments (we call these *edge parameters*). These are handled slightly differently. How to designate parameters as either *edge parameters* or *compartment parameters* is described [here](@ref spatial_lattice_modelling_intro_simulation_edge_parameters).

Below we describe how to set non-uniform values in the various cases.

### [Non-uniform compartment values for Cartesian lattices](@id spatial_lattice_modelling_intro_simulation_inputs_cartesian)
To provide non-uniform values across a Cartesian lattice, simply provide the values in an array of the same dimension and size as the Cartesian lattice. E.g. for a 5x10 Cartesian lattice:
```@example spatial_intro_nonuniform_vals
using Catalyst # hide
ccart_lattices = CartesianGrid((5, 10))
nothing # hide
```
random values (uniformly distributed between $0$ and $1$) can be provided using
```@example spatial_intro_nonuniform_vals
[:X1 => rand(5, 10), :X2 => 10.0]
nothing # hide
```
Non-uniform values for parameters (which values are tied to compartments) are provided similarly.

### [Non-uniform compartment values for masked lattices](@id spatial_lattice_modelling_intro_simulation_inputs_masked)
Non-uniform values for masked lattices are provided in the same manner as for Cartesian lattices (however, values at coordinates that do not hold compartments are ignored). E.g. To provide random values for a masked lattice contained within a 5x10 Cartesian lattices we can again set:
```@example spatial_intro_nonuniform_vals
[:X1 => rand(5, 10), :X2 => 10.0]
nothing # hide
```
If we want, it is also possible to provide the values as a [*sparse array*](https://github.com/JuliaSparse/SparseArrays.jl) with values only in the coordinates that corresponds to compartments.

### [Non-uniform compartment values for unstructured lattices](@id spatial_lattice_modelling_intro_simulation_inputs_graphs)
In graphs (which are used to represent unstructured lattices) each vertex (i.e. compartment) has a specific index. To set non-uniform values for unstructured lattices, provide a vector where the $i$'th value corresponds to the value in the compartment with index $i$ in the graph. E.g. for a graph with 5 vertexes, where we want $X$ to be zero in all compartments bar one (where it is $1.0$) we use:
```@example spatial_intro_nonuniform_vals
[:X1 => [0.0, 0.0, 0.0, 0.0, 1.0], :X2 => 10.0]
nothing # hide
```

### [Non-uniform values for edge-parameters](@id spatial_lattice_modelling_intro_simulation_inputs_edge_parameters)
Adjacent compartments are connected by edges (with which compartments are connected by edges being defined by the lattice). For unstructured lattices, it is possible (if a directed graph was used) to have edges from one compartment to another, but not in the opposite direction. For a lattice with $N$ compartments, edge values are set by a $NxN$ matrix, where value $(i,j)$ corresponds to the parameter's values in the edge going *from* compartment $i$ *to* compartment $j$. This matrix can be either [sparse or non-sparse](https://docs.julialang.org/en/v1/stdlib/SparseArrays/). In the latter cases, values corresponding to non-existing edges are ignored. 

Let's consider a 1d Cartesian lattice with 4 compartments. Here, an edge parameter's values are provided in a 4x4 matrix. For [the Brusselator model described previously](@ref spatial_lattice_modelling_intro_example), $D$'s value was tied to edges. If we wish to set the value of $D$ to various values between $0.1$ and $0.4$ we can do:
```@example spatial_intro_nonuniform_vals
ps = [:k1 => 1.0, :k2 => 2.0, 
      :D => [
        0.0 0.1 0.0 0.0;
        0.1 0.0 0.2 0.0;
        0.0 0.2 0.0 0.3;
        0.0 0.0 0.3 0.0]
]
nothing # hide
```
Here, the value at index $i,j$ corresponds to $D$'s value in the edge from compartment $i$ to compartment $j$. `0.0` is used for elements that do not correspond to an edge. The [`make_edge_p_values`](@ref) and [`make_directed_edge_values`](@ref) provide convenient interfaces for generating non-uniform edge parameter values.

## [Edge parameters and compartment parameters](@id spatial_lattice_modelling_intro_simulation_edge_parameters)
Parameters can be divided into *edge parameters* and *compartment parameters* (initial condition values are always tied to compartments). Here, edge parameters have their values tied to edges, while compartment parameters have their values tied to compartments. All parameters that are part of the rates (or stoichiometries) of non-spatial reactions must be compartment parameters. Parameters that are part of spatial reactions can be either compartment parameters or edge parameters. When a spatial reaction's rate is computed, edge parameters fetch their values for from the edge of the transition, and compartment parameters from the compartment *from which the edge originates*.

When a [`LatticeReactionSystem`](@ref) is created, its parameters is the union of all parameters occurring in the (non-spatial) [`ReactionSystem`](@ref) and in all spatial reactions. By default, parameters occurring only in spatial reactions are considered edge parameters (and if they occur in the non-spatial [`ReactionSystem`](@ref) they are considered compartment parameters). It is, however, possible to designate a parameter specifically as an edge parameter (or not), by using the `edgeparameter` [metadata](@ref dsl_advanced_options_species_and_parameters_metadata). E.g. to designate that `D` (when declared in a non-spatial [`ReactionSystem`](@ref) using the DSL) is an edge parameter, not a compartment parameter, we use:
```@example spatial_intro_edge_ps
using Catalyst # hide
two_state_model = @reaction_network begin
    @parameters D [edgeparameter=true]
    (k1,k2), X1 <--> X2
end
nothing # hide
```

To learn the compartment and edge parameters of a `LatticeReaction`, the [`vertex_parameters`](@ref) and [`edge_parameters`](@ref) functions can be used:
```@example spatial_intro_edge_ps
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
diffusion_rx = @transport_reaction D X1
lattice = CartesianGrid((20,20))
lrs = LatticeReactionSystem(two_state_model, [diffusion_rx], lattice)
edge_parameters(lrs)
```

## [Spatial modelling limitations](@id spatial_lattice_modelling_intro_limitations)
Many features which are supported for non-spatial `ReactionSystem`s are currently unsupported for [`LatticeReactionSystem`](@ref)s. This includes [observables](@ref dsl_advanced_options_observables), [algebraic and differential equations](@ref constraint_equations), [hierarchical models](@ref compositional_modeling), and [events](@ref constraint_equations_events). It is possible that these features will be supported in the future. Furthermore, [removal of conserved quantities](@ref network_analysis_conservation_laws) is not supported when creating spatial `ODEProblem`s.

If you are using Catalyst's features for spatial modelling, please give us feedback on how we can improve these features. Additionally, just letting us know that you use these features is useful, as it helps inform us whether continued development of spatial modelling features is worthwhile. 
