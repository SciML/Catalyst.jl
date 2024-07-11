# [Plotting Spatial Simulations](@id lattice_simulation_plotting)

To aid the investigation of spatial simulations we have implemented a helper functions for creating plots and animation of the simulations. The development of these functions are currently under development, with entries being added here as they are created. A list of currently implemented spatial plotting functions can be found in the [api](@ref ).

!!! note
    The plotting interfaces used for non-spatial Catalyst simulations have seen lots of work to ensure high quality plots. However, the corresponding functions for spatial simulations are primarily intended to aid the user to investigate their simulation results. Thus, they might not be fully suitable for e.g. creating publication-quality graphics. If you are using these functions, please let us now. This helps inform us whether continued development of spatial modelling features is worthwhile. 

!!! note
    To create animations we use [Makie.jl](https://docs.makie.org/stable/), which is an alternative plotting package to [Plots.jl](https://github.com/JuliaPlots/Plots.jl) (which is typically the preferred plotting package within the context of Catalyst). Generally, Makie is good at creating animations, hence we use it here (however, it is also a [popular competitor to Plots.jl for general-purpose plotting](https://juliapackagecomparisons.github.io/pages/plotting/)). 

!!! warn
    These plotting interfaces are very much a work in progress. Hence, they and their interfaces may see more change that Catalyst features typically do. This include *the possibility of breaking changes without breaking releases to Catalyst*.

## [Animations of 2d Cartesian or masked lattices](@id lattice_simulation_plotting_2d_grid_animations)
Spatial simulations based on two dimensional [Cartesian] or [masked] lattices can both be animated using the same `lattice_animation` function. It four three arguments, the simulation you which animate, the species which values to display, the name of the file to which you wish to save the simulation, and the `LatticeReactionSystem` which was simulated. 
