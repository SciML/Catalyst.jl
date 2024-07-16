# [Plotting Spatial Simulations](@id lattice_simulation_plotting)

To aid the investigation of spatial simulations we have implemented a helper functions for creating plots and animation of the simulations. The development of these functions are currently under development, with entries being added here as they are created. A list of currently implemented spatial plotting functions can be found in the [api](@ref ).

!!! note
    The plotting interfaces used for non-spatial Catalyst simulations have seen lots of work to ensure high quality plots. However, the corresponding functions for spatial simulations are primarily intended to aid the user to investigate their simulation results. Thus, they might not be fully suitable for e.g. creating publication-quality graphics. If you are using these functions, please let us now. This helps inform us whether continued development of spatial modelling features is worthwhile. 

!!! note
    To create animations we use [Makie.jl](https://docs.makie.org/stable/), which is an alternative plotting package to [Plots.jl](https://github.com/JuliaPlots/Plots.jl) (which is typically the preferred plotting package within the context of Catalyst). Generally, Makie is good at creating animations, hence we use it here (however, it is also a [popular competitor to Plots.jl for general-purpose plotting](https://juliapackagecomparisons.github.io/pages/plotting/)). 

!!! warn
    These plotting interfaces are very much a work in progress. Hence, they and their interfaces may see more change that Catalyst features typically do. This include *the possibility of breaking changes without breaking releases to Catalyst*.

## [Animations of 2d Cartesian or masked lattices](@id lattice_simulation_plotting_2d_grid_animations)
Spatial simulations based on two dimensional [Cartesian] or [masked] lattices can both be animated using the same `lattice_animation` function. It four three arguments, the simulation you which animate, the species which values to display, the name of the file to which you wish to save the simulation, and the `LatticeReactionSystem` which was simulated. Here we first create and simulate a spatial [Brusselator](@ref basic_CRN_library_brusselator) model.
```@example lattice_plotting
using Catalyst, OrdinaryDiffEq
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
tspan = (0.0, 40.0)
ps = [:A => 1.0, :B => 4.0, :D => 0.2]
oprob = ODEProblem(lrs, u0, tspan, ps)
sol = solve(oprob, FBDF())
nothing # hide
```
Next, we can animate it using the `lattice_animation` function, which requires the [CairoMakie.jl](https://github.com/JuliaPlots/CairoMakie.jl) package:
```@example lattice_plotting
using CairoMakie
lattice_animation(sol, :X, lrs, "brusselator.mp4")
nothing # hide
```
```@raw html
<video autoplay loop muted playsinline controls src="./brusselator.mp4" />
```
```@example lattice_plotting
rm("brusselator.mp4")
```

Here, `lattice_animation` takes various additional optional arguments that can be useful:
- `colormap`: The [colour map](https://docs.makie.org/v0.21/explanations/colors#misc) to used in the heatmap to display the change in species amounts. The default value is `:BuGn_7`.
- `nframes`: The number of frames in the animation. These are sampled evenly across the simulation timespan. The default value is `200`.
- `framerate`: The animations framerate. The default value is `20`.
- Many additional arguments which can be provided to normal Makie heat maps can be provided to `lattice_animation` as well.