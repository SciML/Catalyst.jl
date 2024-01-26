# [Introduction to spatial modelling with Catalyst](@id spatial_lattice_modelling_intro)
Catalyst supports the expansion of non-spatial `ReactionSystem` (created with e.g. the `@reaction_network` DSL) to spatial domains. Spatial simulation of Catalyst models is a work in progress. Currently, the following is supported:
- Spatial ODE and Jump simulations.
- Discrete spatial domains.
- Constant-rate transportation reactions (a species moving at a constant rate across the domain). 

Features which support is planned in future updates include:
- Models on continuous domains (that are automatically discretised).
- SDE simulations.
- Transportation reactions with non-constant rates and general spatial reactions.

This tutorial introduces spatial modelling on discrete domains. It describes the basics of creating such models, simulating them, and interacting with the result. To do so, it uses ODE simulations as examples. Additional tutorials provide [in-depth description of spatial ODE simulations](@ref spatial_lattice_ode_simulations), and [introduces spatial jump simulations](@ref spatial_lattice_jump_simulations).

## [Basic example of spatial simulations on a discrete domain](@id spatial_lattice_modelling_intro_example)
To perform discrete-space spatial simulations, the user must first define a `LatticeReactionSystem`. These combine:
- A (non-spatial) `ReactionSystem` model (created using standard Catalyst Syntax).
- A vector of spatial reactions, describing how species can move spatially across the domain.
- A lattice defining the spatial domain's compartments and how they are connected.

Here, as an example, we will simulate a spatial [Brusselator](https://en.wikipedia.org/wiki/Brusselator). To do so, we first define our (non-spatial) model, the spatial reactions, and the lattice. These are then bundled into a `LatticeReactionSystem`.
```@example spatial_intro_1
using Catalyst
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
diffusion_reaction = @transport_reaction D X
lattice = CartesianGrid((20,20))
lrs = LatticeReactionSystem(brusselator, [diffusion_reaction], lattice)
nothing # hide
```
Here we define
- A single spatial reaction, a transport reaction where $X$ moves at constant rate $D$ between adjacent compartments.
- A 2d Cartesian grid of $20x20$ compartments to simulate our model on.

More details on spatial reactions are available [here](@ref spatial_lattice_modelling_intro_spatial_reactions). In addition to Cartesian grid lattices (in 1, 2, or 3 dimensions), masked and unstructured (graph) grids are also supported. These are described in more detail [here](@ref spatial_lattice_modelling_intro_lattices).

Once created, `LatticeReactionSystem`s can be used as input to various problem types, which then can be simulated in the same manner as non-spatial models. Here, we prepare an ODE simulation by creating an `ODEProblem`:
```@example spatial_intro_1
u0 = [:X => rand(20,20), :Y => 10.0]
tspan = (0.0, 200.0)
ps = [:A => 1.0, :B => 4.0, :D => 0.2]
oprob = ODEProblem(lrs, u0, tspan, ps)
nothing # hide
```
Here we used non-uniform values for $X$'s initial condition, but uniform values for the remaining initial condition and parameter values. More details of uniform and non-uniform initial conditions and parameter values are provided [here](@ref spatial_lattice_modelling_intro_simulation_inputs). Furthermore, the diffusion reaction introduces a parameter $D$ (determining $X$'s diffusion rate) for which we must provide a value. 

We can now simulate our model (using the same syntax as for non-spatial ones):
```@example spatial_intro_1
sol = solve(oprob, QNDF())
nothing # hide
```
We note that simulations of spatial models are often computationally expensive. Here we use the `QNDF` solver, which typically performs well for large systems. More advice on the performance of spatial simulations is provided [here](@ref ref). 

Finally, we can access the value of $X$ in comaprtment $(3,4)$ using
```@example spatial_intro_1
sol[:X][3,4]
```
Alternatively, we can create an animation of our simulation using:
```@example spatial_intro_1
using Makie
lattice_animation(sol, lrs)
```
More information on how to retrieve values from spatial simulations can be found [here](@ref spatial_lattice_modelling_intro_solutions), and for plotting them, [here](@ref spatial_lattice_modelling_intro_plotting). Finally, a list of functions for querying `LatticeReactionSystems` for various properties can be found [here](@ref ref).

## [Spatial reactions](@id spatial_lattice_modelling_intro_spatial_reactions)
Spatial reactions describe reaction events which involve species across two connected compartments. Currently, only so-called *transportation reactions* are supported. These consist of:
- A rate at which it occurs. As for non-spatial reactions, this can be any expression. However, currently, it may only consist of parameters. 
- A single species which is transported from one compartment to an adjacent one.

At the occurrence of a transport reaction, the specific species moves to the adjacent compartment. Most common spatial models can be represented using transport reactions only. These can model phenomena such as diffusion or constant flux. A transportation reaction can be created using the `@transportation_reaction` macro. E.g. above we used
```@example spatial_intro_3
using Catalyst # hide
diffusion_reaction = @transport_reaction D X
nothing # hide
```
to create a reaction where species $X$ moves at a constant rate $D$ between adjacent compartments (this rate also scales linearly with the concentration of $X$ in the source compartment). Transport reactions may have rates depending on several parameters. E.g. to model a Brusselator where both species ($X$ and $Y$) are transported at a rate which depends both on the species, but also on some non-uniform parameter which is unique to each connection (e.g. representing the area connecting two cells in a tissue) we could do:
```@example spatial_intro_3
dr_X = @transport_reaction Dx*A X
dr_Y = @transport_reaction Dy*A Y
nothing # hide
```

### Creating transport reactions programmatically
If models are created [programmatically](@ref ref) it is also possible to create transportation reactions programmatically. To do so, use the `TransportatioNReaction` constructor, providing first the rate, and then the transported species:
```@example spatial_intro_3
@variables t
@species X(t) Y(t)
@parameters A B D
tr_X = TransportationReaction(D, X)
nothing # hide
```

## [Defining discrete spatial domains (lattices)](@id spatial_lattice_modelling_intro_lattices)
Discrete spatial domains can represent:
1. Systems which are composed of a (finite number of) compartments, where each compartment can be considered well-mixed (e.g. can be modelled non-spatially) and (potentially) species can move between different compartments. Tissues, where each compartment corresponds oa biological cell, can be examples of such systems.
2. Systems that are continuous in nature, but have been approximated as a discrete domain. Future updates will include the ability for definition, and automatic discretization, of continuous domains. Currently, however, the user has to perform this discretisation themselves.

Catalyst supports three distinct types of lattices:
- Cartesian grids. These are grids, where each grid point corresponds to a compartment. Spatial transportation is permitted between adjacent compartments.
- Masked grids. These are defined grids, however, only a subset of the grid point actually corresponds to compartments. Spatial transportation is permitted between adjacent compartments.
- Unstructured (or graph) grids. These are defined by graphs, where vertices correspond to compartments and edges connect adjacent compartments.

Here, Cartesian grids are a subset of the masked grids, which are a subset of the unstructured grids. If possible, it is advantageous to use as narrow a lattice definition as possible (this may both improve simulation performance, but also simplifies syntax). Cartesian and masked grids can be defined both in one, two, and three dimensions. By default, these grids assume that diagonally neighbouring compartments are non-adjacent (do not permit direct movement of species in between themselves). To change this, provide the `diagonally_adjacent=true` argument to your `LatticeReactionSystem` when it is created.

### Defining Cartesian grids.
A Cartesian grid is defined using the `CartesianGrid` function, which takes a single argument. For a 1d grid, simply provide the length of the grid as a single argument:
```@example spatial_intro_3
using Catalyst # hide
cgrid_1d = CartesianGrid(5)
nothing # hide
```
For 2d and 3d grids, instead provide a Tuple with the length of the grid in each dimension:
```@example spatial_intro_3
cgrid_2d = CartesianGrid((3,9))
cgrid_3d = CartesianGrid((2,4,8))
nothing # hide
```

### Defining masked grids.
Masked grids are defined through 1d, 2d, or 3d Boolean arrays. Each position in the array is `true` if it corresponds to a compartment, and `false` if it does not. E.g. to define a 1d grid corresponding to two disjoint sets, each with 3 compartments, use:
```@example spatial_intro_3
rgrid_1d = [true, true, true, false, true, true, true]
nothing # hide
```
To define a 2d grid corresponding to the shape of an "8", we can use:
```@example spatial_intro_3
rgrid_2d = [
    true  true  true  true  true;
    true  false true  false true;
    true  true  true  true  true
]
nothing # hide
```
Finally, a 4x5x6 3d grid of randomly distributed compartments can be created using:
```@example spatial_intro_3
rgrid_3d = rand(4,5,6, [true, false])
nothing # hide
```

### Defining unstructured grids.
To define unstructured grids, we must first import the [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) package. Next, we can either use some [pre-defined formula for building graphs](https://juliagraphs.org/Graphs.jl/dev/core_functions/simplegraphs_generators/#Generators-for-common-graphs), or [build a graph from scratch](https://juliagraphs.org/Graphs.jl/dev/first_steps/construction/). Here we create a cyclic graph (where each compartment is connected to exactly two others):
```@example spatial_intro_3
cycle_graph(7)
nothing # hide
```

Since graphs can represent any network of connected compartments, they do not have dimensions (like Cartesian or masked grids). Another feature of graph lattices is that they can have non-symmetric connections, i.e. pairs of compartments where spatial movement of species is only permitted in one direction (in practise, this can be done for Cartesian and masked grids as well, by [defining non-uniform spatial rates](@ref spatial_lattice_modelling_intro_simulation_inputs) and setting them to zero in one direction). This can be done by using a [*directed graph*](https://juliagraphs.org/Graphs.jl/dev/algorithms/digraph/) as input. E.g. here we define a directed cyclic graph, where movement is only allowed in one direction of the cycle:
```@example spatial_intro_3
cycle_digraph(7)
nothing # hide
```

## [Non-uniform initial conditions and parameter values](@id spatial_lattice_modelling_intro_simulation_inputs)
For spatial models, initial conditions and parameter values are provided similarly as for non-spatial models. Wherever a single value is provided, it is used *uniformly* across the lattice. E.g. if we, for our previous Brusselator model, set
```@example spatial_intro_4
u0 = [:X => 1.0, :Y => 10.0]
ps = [:A => 1.0, :B => 4.0, :D => 0.2]
nothing # hide
```
The initial condition will be $1.0$ for $X$ across compartments, and of $2.0$ for $Y$. Furthermore, for each Brusselator simulation, in each compartment, the value of $1.0$ will be used for $A$ and $4.0$ for $B$. Finally, the transportation rate of $X$ (set by the parameter $D$) will be $0.2$ across all connections. To set non-uniform values, non-scalar values must be provided in the map. How to do this depends on how the [lattice was defined](@ref spatial_lattice_modelling_intro_lattices). Furthermore, some parameters that are part of spatial reactions may have their value tied to *connections between compartments*, rather than compartments (we call these *edge parameters*). These are handled slightly differently. How to designate parameters as either *edge parameters* or *compartment parameters* is described [here](@ref spatial_lattice_modelling_intro_simulation_edge_parameters).

Below we describe how to set non-uniform values in the various cases.

### Non-uniform compartment values for Cartesian grids
To provide non-uniform values across a Cartesian grid, simply provide the values in an array of the same dimension and size as the Cartesian grid. E.g. for a $5x10$ Cartesian grid:
```@example spatial_intro_4
cgrid = CartesianGrid((5,10))
nothing # hide
```
random values (uniformly distributed between $0$ and $1$) can be provided using
```@example spatial_intro_4
[:X => rand(5,10), :Y => 10.0]
nothing # hide
```
Non-uniform values for parameters (which values are tied to compartments) are provided similarly.

### Non-uniform compartment values for masked grids
Non-uniform values for masked grids are provided in the same manner as for Cartesian grids (however, values at coordinates that do not hold compartments are ignored). E.g. To provide random values for a masked grid contained within a $5x10$ Cartesian grid we can again set:
```@example spatial_intro_4
[:X => rand(5,10), :Y => 10.0]
nothing # hide
```
If we want, it is also possible to provide the values as a [*sparse array*](https://github.com/JuliaSparse/SparseArrays.jl) with values only in the coordinates corresponding to compartments.

### [Non-uniform compartment values for unstructured grids](@id spatial_lattice_modelling_intro_simulation_inputs_unstructured)
In graphs (which are used to represent unstructured grids) each vertex has a specific index. To set non-uniform values for unstructured grids, provide a vector, where the i'th value corresponds to the value in the compartment with index i in the graph. E.g. for a graph with 5 vertexes, where we want $X$ to be zero in all compartments bar one (where it is $1.0$) we use:
```@example spatial_intro_4
[:X => [0.0, 0.0, 0.0, 0.0, 1.0], :Y => 10.0]
nothing # hide
```

### Non-uniform values for edge-parameters
Adjacent compartments are connected by edges (with which compartments are connected by edges being defined by the lattice). For unstructured lattices, it is possible (if a directed graph was used) to have edges from one compartment to another, but not in the opposite direction. For a lattice with $n$ compartments, edge values are set by a $nxn$ matrix, where value $(i,j)$ corresponds to the parameter's values in the edge going *from* compartment i *to* compartment j. This matrix can be either [sparse or non-sparse](https://docs.julialang.org/en/v1/stdlib/SparseArrays/). In the latter cases, values corresponding to non-existing edges are ignored. 

E.g. let's consider a 1d Cartesian grid with 4 compartments. Here, a spatial parameter's values are provided in a $4x4$ matrix. For [the Brusselator model described previously](@ref spatial_lattice_modelling_intro_example), $D$'s value was tied to edges. If we wish to set the value of $D$ to various values between $0.1$ and $0.4$ we can do:
```@example spatial_intro_4
ps = [:A => 1.0, :B => 4.0, 
      :D => [
        0.0 0.1 0.0 0.0;
        0.1 0.0 0.2 0.0;
        0.0 0.2 0.0 0.3;
        0.0 0.0 0.3 0.0]
]
nothing # hide
```
Here, `0.0` is used for indexes that do not correspond to an edge.

## [Edge parameters and compartment parameters](@id spatial_lattice_modelling_intro_simulation_edge_parameters)
Parameters can be divided into *edge parameters* and *compartment parameters* (initial condition values are always tied to compartments). Here, edge parameters have their values tied to edges, while compartment parameters have their values tied to compartments. All parameters that are part of the rates (or stoichiometries) of non-spatial reactions must be compartment parameters. Parameters that are part of spatial reactions can be either compartment parameters or edge parameters. When a spatial reaction's rate is computed, edge parameters fetch their values for from the edge of the transition, and compartment parameters from the compartment in which *the edge originates*.

When a `LatticeReactionSystem` is created, its parameter is the union of all parameters occurring in the (non-spatial) `ReactionSystem` and in all spatial reactions. By default, parameters occurring only in spatial reactions are considered edge parameters (and if they occur in the non-spatial `ReactionSystem` they are considered compartment parameters). It is however possible to designate a parameter specifically as an edge parameter or not, by using the `edgeparameter` [metadata](@ref ref). E.g. to designate that `D` (when declared in a non-spatial `ReactionSystem` using the DSL) is an edge parameter, not a compartment parameter, we use:
```@example spatial_intro_4
using Catalyst # hide
brusselator = @reaction_network begin
    @parameters D [edgeparameter=true]
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
```

To learn the compartment and edge parameters of a `LatticeReaction`, the `compartment_parameters` and `edge_parameters` functions can be used:
```@example spatial_intro_4
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
diffusion_reaction = @transport_reaction D X
lattice = CartesianGrid((20,20))
lrs = LatticeReactionSystem(brusselator, [diffusion_reaction], lattice)
edge_parameters(lrs)
```

## [Interacting with spatial simulation solutions and integrators](@id spatial_lattice_modelling_intro_solutions)

We [have previously described](@ref simulation_structure_interfacing) how to interface with problems, integrators, and solutions produced by non-spatial models. Here, for a simulation solution `sol`, to retrieve a vector with all values of $X$ we simply use:
```julia
sol[:X]
```
A similar syntax can be used for spatial models. However, the output will be value of $X$ across the entire spatial grid. How to retrieve (and for integrators and problems, set) values at specific compartments, a specific notation will have to be used. This also depends on the type of lattice used.

### Retrieving solution values from Cartesian grids
Let us consider the solution to a `LatticeReactionSystem` where a Cartesian grid was used:
```@example spatial_intro_5
using Catalyst, OrdinaryDiffEq
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
diffusion_reaction = @transport_reaction D X
lattice = CartesianGrid((20,20))
lrs = LatticeReactionSystem(brusselator, [diffusion_reaction], lattice)

u0 = [:X => rand(20,20), :Y => 10.0]
tspan = (0.0, 200.0)
ps = [:A => 1.0, :B => 4.0, :D => 0.2]
oprob = ODEProblem(lrs, u0, tspan, ps)
sol = solve(oprob, QNDF())
nothing # hide
```
To retrieve the value of $X$ at all time points we use
```@example spatial_intro_5
sol[:X]
```
Here, `sol[:X]` is a `Vector{Matrix{Float64}}`. Each element in the vector corresponds to the values of $X$ at that time point. In this case, a 2d Cartesian grid was used, so the values of $X$ is a matrix index $i,j$ is the value of $X$ in compartment $i,j$. For 1d and 3d Cartesian grids, The values of $X$ in each time point are instead a `Vector{Float64}` or a `Array{Float64,3}`, respectively.

### Retrieving solution values from masked grids
If a masked grid had been used instead, e.g for 
```@example spatial_intro_5
lattice = [true true true; true false true; true true true]
lrs = LatticeReactionSystem(brusselator, [diffusion_reaction], lattice)
nothing # hide
```
The value of `sol[:X]` returns a vector of sparse arrays, containing the values of $X$ in all the compartments. The dimension of the sparse array is equal to the dimension of the masked grid which was provided as the lattice.

### Retrieving solution values from unstructured grids
If an unstructured grid (graph) was used as a lattice
```@example spatial_intro_5
lattice = cycle_graph(7)
lrs = LatticeReactionSystem(brusselator, [diffusion_reaction], lattice)
nothing # hide
```
The value of `sol[:X]` returns a `Vector{Vector{Float64}}`. Here, the top-layer vector corresponds to the time points, and the bottom-layer to the value of $X$ in each compartment. Previously, [we described how to give initial conditions for unstructured grids by providing a vector where each element corresponds to a vertex in the grid](@ref spatial_lattice_modelling_intro_simulation_inputs_unstructured). Similarly, here, the i'th element in the bottom-layer vector is the value of $X$ in the i'th vertex of the graph that was provided as the lattice.

### Retrieving and setting values for problems and integrators

## [Plotting Spatial simulations](@id spatial_lattice_modelling_intro_plotting)

To aid the visualisation of spatial simulations, we provide 6 different plot recipes, covering the following situations:
- For a 1d structured (Cartesian or masked) grid, plot its value at a single time point.
- For a 1d structured grid, plot a kymograph of a single species's value over time.
- For a 2d structured grid, plot a single species's value at a single time point.
- For a 2d structured grid, create an animation of a single species's value over time.
- For an unstructured grid, where the user provides the coordinates of each vertex, plot a species's value at a single time point. 
- For an unstructured grid, where the user provides the coordinates of each vertex,  create an animation of a single species's value over time.

### Plotting 1d structured lattice simulation at single time points

### Plotting 1d structured lattice simulation over time

### Plotting 2d structured lattice simulation at single time points

### Animating 2d structured lattice simulation over time

### Plotting unstructured lattice simulation at single time points

### Animating unstructured lattice simulation over time
