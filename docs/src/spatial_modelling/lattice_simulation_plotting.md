# [Plotting Spatial Simulations](@id lattice_simulation_plotting)

To aid the investigation of spatial simulations we have implemented several helper functions for creating plots and animation of these simulations. The development of these functions are currently under development. Currently, Catalyst exports three different such functions:
- [`lattice_plot`](@ref): Which plots a lattice simulation at a specified time point.
- [`lattice_animation`](@ref): Which creates an animation of a lattice simulation across the entire time series.
- [`lattice_kymograph`](@ref): Which, for a 1d Cartesian or masked lattice based simulation, creates a kymograph of the simulation across time and space.

The first two functions can be applied to [graph](@ref spatial_lattice_modelling_intro_lattices_graph) and 1d and 2d [Cartesian](@ref spatial_lattice_modelling_intro_lattices_cartesian) and [masked](@ref spatial_lattice_modelling_intro_lattices_masked) lattice based simulations, while `lattice_kymograph` can only be applied to 1d Cartesian and masked lattice based simulations. Currently, there is no functionality for plotting simulations based on 3d Cartesian and masked lattice. Here we will demonstrate all plotting functions using ODE simulations, but they work equally well for jump simulations.

!!! note
    The plotting interfaces used for non-spatial Catalyst simulations have seen lots of work to ensure high quality plots. However, the corresponding functions for spatial simulations are primarily intended to aid the user to investigate their simulation results. Thus, they might not be fully suitable for e.g. creating publication-quality graphics. If you are using these functions, please let us now. This helps inform us whether continued development of spatial modelling features is worthwhile. 

!!! note
    To create animations we use [Makie.jl](https://docs.makie.org/stable/), which is an alternative plotting package to [Plots.jl](https://github.com/JuliaPlots/Plots.jl) (which is typically the preferred plotting package within the context of Catalyst). Generally, Makie is good at creating animations, hence we use it here (however, it is also a [popular competitor to Plots.jl for general-purpose plotting](https://juliapackagecomparisons.github.io/pages/plotting/)). 

!!! warning
    These plotting interfaces are a work in progress. Hence, they and their interfaces may see more change that Catalyst features typically do. This include *the possibility of breaking changes without breaking releases to Catalyst*.

## [Animation and plotting of 1d Cartesian or masked lattice simulations](@id lattice_simulation_plotting_1d_grids)
Let us consider a spatial simulation on a 1d Cartesian grid lattice:
```@example lattice_plotting_1d
using Catalyst, OrdinaryDiffEqDefault
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
diffusion_rx = @transport_reaction D X1
lattice = CartesianGrid(10)
lrs = LatticeReactionSystem(two_state_model, [diffusion_rx], lattice)

X1_0 = ones(10)
X1_0[5] = 10.0
u0 = [:X1 => X1_0, :X2 => 1.0]
tspan = (0.0, 5.0)
ps = [:k1 => 1.0, :k2 => 2.0, :D => 0.5]
oprob = ODEProblem(lrs, u0, tspan, ps)
sol = solve(oprob)
nothing # hide
```
To plot the simulation at a specific time point we use the `lattice_plot` function. In addition to the simulation we wish to plot, it takes the species we wish to plot and the [`LatticeReactionSystem`](@ref) which generated the simulation as arguments. It also takes the time point at which we wish to plot the simulation as an additional argument (if not provided, the simulation will be plotted at the final time point). To use the [`lattice_plot`](@ref) function (or any other of Catalyst's spatial plotting functions) we also need to load the [CairoMakie](https://github.com/MakieOrg/Makie.jl) package (here, `import CairoMakie` is enough, and `using CairoMakie` is not required).
```@example lattice_plotting_1d
import CairoMakie
lattice_plot(sol, :X1, lrs; t = 2.0)
```

If we instead wish to create an animation of our solution across the entire simulation, we can use the [`lattice_animation`](@ref) function. This takes a fourth required argument, a file to which the animation is saved.
```@example lattice_plotting_1d
lattice_animation(sol, :X1, lrs, "lattice_simulation_1d.mp4")
```
![](./lattice_simulation_1d.mp4)
Since we animate the solution across the entire simulation, we do not need to provide a `t` value. However, there are some additional (optional) arguments we might wish to provide:
- `nframes = 200`: The number of frames in the animation (these are evenly samples across the simulation).
- `framerate = 20`: The frame rate of the animation.

Finally, we can display a kymograph of our simulation across the full time span using [`lattice_kymograph`](@ref).
```@example lattice_plotting_1d
lattice_kymograph(sol, :X1, lrs)
```
Here, we require neither a filename nor a `t` to be provided. However, the `nframes` argument can still be used to determine how frequently (in time) we wish to sample our simulation.

For more information of either function, and additional optional arguments, please read their corresponding api sections ([`lattice_plot`](@ref), [`lattice_animation`](@ref), and [`lattice_kymograph`](@ref)).

## [Animation and plotting of 2d Cartesian or masked lattice simulations](@id lattice_simulation_plotting_2d_grids)
Two-dimensional lattice simulations can be plotted in the same manner as one-dimensional ones. However, instead of displaying a species's value as a line plot, it is displayed as a heatmap. E.g. here we simulate a spatial [Brusselator](@ref basic_CRN_library_brusselator) model and display the value of $X$ at a designated time point.
```@example lattice_plotting_2d
using Catalyst, OrdinaryDiffEqBDF
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
diffusion_rx = @transport_reaction D X
lattice = CartesianGrid((20,20))
lrs = LatticeReactionSystem(brusselator, [diffusion_rx], lattice)

u0 = [:X => rand(20, 20), :Y => 10.0]
tspan = (0.0, 20.0)
ps = [:A => 1.0, :B => 4.0, :D => 0.2]
oprob = ODEProblem(lrs, u0, tspan, ps; jac = true, sparse = true)
sol = solve(oprob, FBDF())

import CairoMakie
lattice_plot(sol, :X, lrs; t = 18.0)
```
Here, the intensity of the colour scale describes $X$'s value in each compartment. The color scale is normalised across the entire simulation (not the plotted time step).

An animation of the solution can be created in a similar manner as for [the one-dimensional case](@ref lattice_simulation_plotting_1d_grids):
```@example lattice_plotting_2d
lattice_animation(sol, :X, lrs, "lattice_simulation_2d.mp4")
```
![](./lattice_simulation_2d.mp4)

Again, please check the API pages for the [`lattice_plot`](@ref) and [`lattice_animation`](@ref) functions to see more details of their various options.

## [Animation and plotting of graph lattice simulations](@id lattice_simulation_plotting_graphs)
Finally, we consider lattice simulations on graph lattices. We first simulate a simple [birth-death process](@ref basic_CRN_library_bd) on a (6-node cyclic) graph lattice.
```@example lattice_plotting_graphs
using Catalyst, Graphs, OrdinaryDiffEqDefault
rs = @reaction_network begin
    (p,d), 0 <--> X
end
diffusion_rx = @transport_reaction D X
lattice = cycle_graph(6)
lrs = LatticeReactionSystem(rs, [diffusion_rx], lattice)

u0 = [:X => rand(6)]
tspan = (0.0, 1.0)
ps = [:p => 1.0, :d => 0.5, :D => 0.1]
oprob = ODEProblem(lrs, u0, tspan, ps)
osol = solve(oprob)
nothing # hide
```
As for Cartesian and masked lattice-based simulations we can plot this simulation using the `lattice_plot` function. In this case, however, we need to first load *both* [CairoMakie](https://github.com/MakieOrg/Makie.jl) and the [GraphMakie](https://github.com/MakieOrg/GraphMakie.jl) packages.
```@example lattice_plotting_graphs
import CairoMakie, GraphMakie
lattice_plot(osol, :X, lrs)
```
Here we provide no `t` argument, and hence plot the solution at the final time point. As for 2d simulations, a colour scale described $X$'s concentration in each compartment. Some arguments, are not relevant when plotting simulations on Cartesian and masked lattices, but which are relevant to graph lattices are:
- `node_size`: Sets the size of each compartment in the plot.
- `layout`: A vector of x and y values which can be used to determine the position of each compartment in the plot.

E.g. here we manually set where each compartment should be plotted:
```@example lattice_plotting_graphs
isq2 = 1/sqrt(2)
layout = [(0, 0), (isq2, isq2), (isq2, 1 + isq2), (0, 1 + 2isq2), (-isq2, 1 + isq2), (-isq2, isq2)]
lattice_plot(osol, :X, lrs; layout)
```

Finally, animations of graph lattice simulation work similarly to 2d ones, but accept the additional arguments relevant to plotting graphs.

## [Final notes](@id lattice_simulation_plotting_notes)
If you are using these interfaces, but there is some feature that is missing, you might wish to consider modifying the original code. This can be found [here](https://github.com/SciML/Catalyst.jl/blob/master/ext/CatalystCairoMakieExtension/cairo_makie_extension_spatial_modelling.jl), from which you can copy any code you need to make your own plotting interfaces. If you do so, please provide any feedback by raising [an issue](https://github.com/SciML/Catalyst.jl/issues) on the Catalyst GitHub page. As mentioned, these plotting interfaces are a work in progress, and input from users is valuable to improve them further. Many of Makie's plotting arguments, even those not described here, are handled by these functions, and someone familiar with the package should be able to use these to customise the plots further.

It should also be noted that these interfaces has note been optimised for performance, and the generation of an animation can often surpass 1 second for larger models. Again, this can likely be improved, and if performance is an problem do raise an issue, in which case an additional effort can be made to improve performance.

