# [Introduction to spatial modelling with Catalyst](@id spatial_lattice_modelling_intro)
Catalyst supports the expansion of non-spatial `ReactionSystem` (created with e.g. the `@reaction_network` DSL) to spatial domains. Spatial simulation of Catalyst models is a work in progress. Currently, the following is supported:
- Spatial ODE and Jump simulations.
- Discrete spatial domains.
- Constant-rate transportation reactions (a species moving at constant rate across the domain). 

Features which support is planned in future updates include:
- Models on continuous domains (that are automatically discretised).
- SDE simulations.
- Transportation reactions with non-constant rates and general spatial reactions.

This tutorial introduces how to create, simulate, and interact with, spatial models in Catalyst. To do so, it uses ODE simulations as examples. Additional tutorials provides [in-depth description of spatial ODE simulations](@ref spatial_lattice_ode_simulations), and [introduces spatial jump simulations](@ref spatial_lattice_jump_simulations).

## Spatial simulations on discrete domains
Discrete spatial domains can represent:
1. Systems which are composed of a (finite number of) compartments, where each compartment can be considered well-mixed (e.g. can be modelled non-spatially) and (potentially) species can move between different compartments. Tissues, where each compartment correspond oa biological cell, can be examples of such systems.
2. Systems that are continuous in nature, but has been approximated as a discrete domain. Future updates will include teh ability for definition, and automatic discretization, of continuous domains. Currently, however, the user have to perform this discretisation themselves.

To perform discrete-space spatial simulations, the user must first define a `LatticeReactionSystem`. These combine:
- A (non-spatial) `ReactionSystem` model (created using standard Catalyst Syntax).
- A vector of spatial reactions, describing how species can move spatially across the domain.
- A lattice defining all compartments and how they are connected.

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
Here we:
- Use a designate a single spatial reaction, a transport reaction where $X$ moves at constant rate $D$ between adjacent compartments.
- Define a 2d Cartesian grid of $20x20$ compartments to simulate our model on.

More details on spatial reactions are available [here](@ref spatial_lattice_modelling_intro_spatial_reactions). In addition to Cartesian grid lattices (1d, 2d, or 3d), regular (but non-Cartesian) and irregular grids are also supported. These are described in more details [here](@ref spatial_lattice_modelling_intro_lattices).

Once created, `LatticeReactionSystem`s can be used as input to various problem types, which then can be simulated in the same manner as non-spatial models. Here, we prepare an ODE simulation by creating and `ODEProblem`:
```@example spatial_intro_1
u0 = [:X => rand(20,20), :Y => 10.0]
tspan = (0.0, 200.0)
ps = [:A => 1.0, :B => 4.0, :D => 0.2]
oprob = ODEProblem(lrs, u0, tspan, ps)
nothing # hide
```
Here we note that we uses non-uniform values for $X$'s initial condition, but uniform values for the remaining initial condition and parameter values. More details of uniform and non-uniform initial conditions and parameter values are provided [here](@ref spatial_lattice_modelling_intro_simulation_inputs). Furthermore, the diffusion reaction introduces a parameter $D$ (determining $X$'s diffusion rate) for which we must provide a value. We can now simulate our model (using the same syntax as for non-spatial ones):
```@example spatial_intro_1
sol = solve(oprob, QNDF())
nothing # hide
```
We note that simulation of spatial models is often computationally expensive. Here we use the `QNDF` solver, which typically performs well for large system. More advice on performance fo spatial simulation is provided [here](@ref ref). Finally, we can create an animation of our simulation using:
```@example spatial_intro_1
using Makie
lattice_animation(sol, lrs)
```
More information on plotting spatial simulations can be found [here](@ref spatial_lattice_modelling_intro_plotting), and on how to retrieve values from them can be found [here](@ref spatial_lattice_modelling_intro_solutions). Finally, a list of functions for querying `LatticeReactioSystems` for various properties can be found [here]().

## [Spatial reactions](@id spatial_lattice_modelling_intro_spatial_reactions)
Spatial reactions describe reaction events which involves species across two connected compartments. Currently, only so-called *transportation reactions* are supported. These consists of:
- A rate at which it occurs. Like for non-spatial reactions, this can be any expression. However, currently, it may only consist of parameters. 
- A single species which is transported from one compartment to an adjacent one.

At the occurrence of a transport reaction, the specific species moves to the adjacent compartment. Most common spatial models can be represented using transport reactions only. These can model phenomenon such as diffusion or constant flux. A transportation reaction can be created using the @transportation_reaction macro. E.g. above we used
```@example spatial_intro_3
using Catalyst # hide
diffusion_reaction = @transport_reaction D X
nothing # hide
```
to create a reaction where species $X$ moves at constant rate $D$ between adjacent compartment (this rate also scales linearly with the concentration of $X$ in the source compartment). Transport reactions may have rates depending on several parameters. E.g. to model a Brusselator where both species ($X$ and $Y$) are transported at a rate which depend both on the species, but also one some non-uniform parameter which is unique to each connection (e.g. representing the area connecting two cells in a tissue) we could do:
```@example spatial_intro_3
dr_X = @transport_reaction Dx*A X
dr_Y = @transport_reaction Dy*A Y
nothing # hide
```

If models are created [programmatically](@ref ref) it is also possible to create transportation reactions programmatically. To do so, use the `TransportatioNReaction` constructor, providing first the rate, and then the tranpsorted species:
```@example spatial_intro_3
@variables t
@species X(t) Y(t)
@parameters A B D
dr_X = TransportationReaction(D, X)
nothing # hide
```

## [Lattices](@id spatial_lattice_modelling_intro_lattices)
Catalyst supports three distinct types of lattices:
- Cartesian grids. These are $lxmxn$ grids, where each grid point correspond to a compartment. Spatial transportation is permitted between adjacent compartments.
- Regular grids. These are defined on $lxmxn$ grids, however, only a subset of the grid point actually correspond to compartments. Spatial transportation is permitted between adjacent compartments.
- Irregular grids. These are defined by graphs, where vertices corresponds to compartments and edges connect adjacent compartments.

Here, Cartesian grids are a subset of the regular grids, which are a subset of the irregular grids. If possible, it is advantageous to use as narrow a lattice definition as possible (this may both improve simulation performance, but also simplifies syntax). Cartesian and regular grids can be defined both in one, two, and three dimensions. By default, these grids assume that diagonally neighbouring compartments are non-adjacent (does not permit direct movement of species in between themselves). To change this, provide the `diagonally_adjacent=true` argument to your `LatticeReactionSystem` when it is created.

Whether these should consider diagonally neighbouring compartments as adjacent is optional.

### Defining Cartesian grids.
A Cartesian grid is defined using the `CartesianGrid` function, which takes a single argument. For a 1d grid, simply provide the length of the grid as a single argument:
```@example spatial_intro_3
using Catalyst # hide
cgrid_1d = CartesianGrid(5)
nothing # hide
```
For 2d and 3d grids, instead provide a Tuple with the grids length in each dimension:
```@example spatial_intro_3
cgrid_2d = CartesianGrid((3,9))
cgrid_3d = CartesianGrid((2,4,8))
nothing # hide
```

### Defining regular grids.
Regular grids are defined through 1d, 2d, or 3d Boolean arrays. Each position in the array is `true` if it correspond to a compartment, and `false` if it does no. E.g. to define a 1d grid corresponding to two disjoint sets, each with 3 compartments, use:
```@example spatial_intro_3
rgrid_1d = [true, true, true, false, true, true]
nothing # hide
```
To define a 2d grid corresponding to the shape of a "8", we can use:
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

### Defining irregular grids.
To define irregular grids, we must first fetch the [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) package. Next, we can either use some [pre-defined formula for building graphs](https://juliagraphs.org/Graphs.jl/dev/core_functions/simplegraphs_generators/#Generators-for-common-graphs), or [build a graph from scratch](https://juliagraphs.org/Graphs.jl/dev/first_steps/construction/). Here we create a cyclic graph (where each compartment is connected to exactly two other):
```@example spatial_intro_3
cycle_graph(7)
nothing # hide
```
Since graphs can represent any network of connected compartments, they do not have dimensions (like Cartesian or regular grids). A feature of graph lattices are they they can have non-symmetric connection, i.e. pairs of compartments where spatial movement of species is only permitted in one direction (in practise, this can be done for Cartesian and regular grids as well, by [defining non-uniform spatial rates](@ref spatial_lattice_modelling_intro_simulation_inputs) and setting them to zero in either direction). This can be done by using a [*directed graph*](https://juliagraphs.org/Graphs.jl/dev/algorithms/digraph/) as input. E.g. here we define a directed cyclic graph, where movement is only allowed in one direction of teh cycle:
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
The initial condition will be $1.0$ for $X$ across compartments, and of $2.0$ for $Y$. Furthermore, for each Brusselator simulation, in each compartment, the value of $1.0$ will be used for $A$ and $4.0$ for $B$. Finally, the transportation rate of $X$ (set by the parameter $D$) will be $0.2$ across all connections. To set non-uniform values, non-scalar values must be provided in the map. How to do this depends on how the [lattice was defined](@ref spatial_lattice_modelling_intro_lattices). Furthermore some parameters that are part of spatial reactions may have their value tied to *connections between compartments*, rather than compartments (so called edge-parameters). These are handled slightly different. How to designate parameters as either *edge parameters* or *compartment parameters* is described [here](@ref spatial_lattice_modelling_intro_simulation_edge_parameters).

Below we describe how to set non-uniform values in the various cases.

### Non-uniform compartment values for cartesian grids
To provide non-uniform values across a Cartesian grid, simply provides the values in an array of the same dimension and size as the cartesian grid. E.g. for a $5x10$ Cartesian grid:
```@example spatial_intro_4
cgrid = CartesianGrid((5,10))
nothing # hide
```
random values (uniformly distributed between $0$ and $1$) can be provided using
```@example spatial_intro_4
[:X => rand(5,10), :Y => 10.0]
nothing # hide
```

Alternatively, the fact that arrays can be indexed with scalars can be used to provide values as a vector:
```@example spatial_intro_4
[:X => rand(50), :Y => 10.0]
nothing # hide
```
Here, the value of compartment $(i,j)$ is given by the corresponding coordinate in the input vector. 

### Non-uniform compartment values for regular grids
Non-uniform values for regular grids are provided in the same manner as for Cartesian grids (however, values at coordinates that does not hold compartments are ignored). E.g. To provide random values for a regular grid contained within a $5x10$ Cartesian grid we can again set:
```@example spatial_intro_4
[:X => rand(5,10), :Y => 10.0]
nothing # hide
```
If we want, it is also possible to provide teh values as a [*sparse array*](https://github.com/JuliaSparse/SparseArrays.jl) with values only in the coordinates corresponding to compartments.

Again, values can also be provided in vector form. However, here it is important to keep in mind that empty compartments are skipped when the vector indexes are determined.

### Non-uniform compartment values for irregular grids
In graphs (which are used to represent irregular grids) each vertex have a specific index. To set non-uniform values for irregular grids, provide a vector, where the i'th value correspond to teh value in the compartment with index i in the graph. E.g. for a graph with 5 vertexes, where we want $X$ to be zero in all compartments bar one (where it is $1$) we use:
```@example spatial_intro_4
[:X => [0.0, 0.0, 0.0, 0.0, 1.0], :Y => 10.0]
nothing # hide
```

### Non-uniform values for edge-parameters
Adjacent compartments are connected by edges. Which compartments are connected bye edges depend on how the lattice was created. For irregular lattices, it is possible (if a directed graph was used) to have edges from one compartment to another, but not in the opposite direction. For a lattice with $n$ compartments, edges values are set by a $nxn$ matrix, where value $(i,j)$ correspond to the parameters values in the edge going *from* compartment i *to* compartment j. This matrix can be either sparse or non-sparse. In the former cases, values corresponding to non-existing edges are ignored. 

E.g. let's consider a 1d Cartesian grid with 4 values. Here, a spatial parameters values i provided in a $4x4$ matrix. E.g. to set the value of $D$ to variosu values between $0.1$ and $0.4$ we can do:
```@example spatial_intro_4
ps = [:A => 1.0, :B => 4.0, 
      :D => [
        0.0 0.1 0.0 0.0;
        0.1 0.0 0.2 0.0;
        0.0 0.2 0.0 0.3;
        0.0 0.0 0.3 0.0
]]
nothing # hide
```
Here, `0.0` is used for indexes that do not correspond to an edge.

To help create such matrices, we provide a a function, `to_sparse_matrix` which converts a `Dict{Int64,Dict{Int64,T}}` to a `SparseMatrix{T}`. Here, if the dictionary is called `val_dict`, if `val_dict[i][j]` is a value, the matrix will have that value in index $(i,j)$, else, that position will be empty.

## [Edge parameters and compartment parameters](@id spatial_lattice_modelling_intro_simulation_edge_parameters)
Parameters can be divided into *edge parameters* and *compartment parameters* (initial condition values are always tied to compartments). Here, edge parameters have their values tied to edges, while compartment parameters have their values tied to compartments. All parameters that are part of the rates (or stoichiometries) of non-spatial reaction systems must be compartment parameters. Parameters that are part of spatial reactions can be either compartment parameters or edge parameters. When a spatial reaction's rate is computed, edge parameters fetch their values for from the edge of the transition, and compartment parameters from the compartment in which *the edge originates*.

When a `LatticeReactionSystem` is created, its parameter is the union of all parameters occurring in the (non-spatial) `ReactionSystem` and in all spatial reactions. By default, parameters occurring only in spatial reaction are considered edge parameters, and if they occur in the non-spatial `ReactionSystem` they are considered compartment parameters. It is however possible to designate a parameter specifically as an edge parameter or not, by using the `edgeparameter` metadata. E.g. to designate that `D` (when declared in a non-spatial `ReactionSystem` using the DSL) is a edge parameter, not a compartment parameter, we use:
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

## [Interacting with spatial simulation solutions and integrators](@id spatial_lattice_modelling_intro_solutions)

We [have previously described](@ref simulation_structure_interfacing) how to interface with problems, integrators, and solutions produced by non-spatial models. Here, for a simulation solution `sol`, to retrive a vector with all values of $X$ we simply use:
```@example spatial_intro_5
sol[:X]
nothing # hide
```
A similar syntax can be used for spatial models. However, the output will be value of $X$ across the entire spatial grid. How to retrieve (and for integrators and problems, set) values at specific compartments, specific notation will have to be used. This also depend on the type of lattice used.

### Retrieving solution values from Cartesian grids

### Retrieving solution values from regular grids

### Retrieving solution values from irregular grids

### Retrieving and setting values for problems and integrators

## [Plotting Spatial simulations](@id spatial_lattice_modelling_intro_plotting)


