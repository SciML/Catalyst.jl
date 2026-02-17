### Exported Plotting Functions ###

# Only contain the API strings. Actual functions are extended by the `CairoMakie` and `GraphMakie`
# extensions which contain the actual code.

"""
    dspace_plot(sol, sp, dsrs::DiscreteSpaceReactionSystem, filename::String; t = sol.tspan[2], kwargs...)

Creates a plot of a `DiscreteSpaceReactionSystem` simulation. The plot is created at the time point
specified by `t` (defaults to the simulation's final time point).

Arguments (all spaces):
- `sol`: The simulation we wish to plot.
- `sp`: The species whose values we wish to plot. Can be provided either in its symbolic form or as a symbol.
- `dsrs`: The `DiscreteSpaceReactionSystem` which was simulated.
- `t = sol.t[end]`: The time point at which we wish to plot the solution

In addition, depending on the type of space used, the following optional arguments might be relevant.

Arguments (1d spaces):
- `markersize = 20`: The size of the markers marking each compartment's value.

Arguments (Graph & 2d spaces):
- `colormap = :BuGn_7`: The colour map with which we display the species amounts in the animation.
- `plot_min = nothing`: The minimum value for the colour scale (values less than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value (across the entire simulation, not just at the plotted time value).
- `plot_max = nothing`: The maximum value for the colour scale (values more than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value (across the entire simulation, not just at the plotted time value).

Arguments (Graph spaces):
- `node_size = 50`: The size of the compartments in the plot.
- `layout = Spring()`: The layout for the graph nodes in the plot. Can be provided as a vector, where the i'th element is a 2-valued tuple (determining the i'th compartment's y and x positions, respectively).

Notes: 
- For masked spaces, there are no value displayed for grid points which do not correspond to a compartments.
- The current plotting interface is a work in progress, and modifications are expected. if you have any feedback, please contact the package authors.
- Additional arguments can be passed to `dspace_plot`, which then will be passed to Makie's `lines` plotting command.
"""
function dspace_plot end

"""
    dspace_animation(sol, sp, dsrs::DiscreteSpaceReactionSystem, filename::String; kwargs...)

Creates an animation of a `DiscreteSpaceReactionSystem` simulation. The animation is saved to a file, 
whose name is provided in the `filename` argument.

Arguments (all spaces):
- `sol`: The simulation we wish to animate.
- `sp`: The species which values we wish to animate. Can be provided either in its symbolic form or as a symbol.
- `dsrs`: The `DiscreteSpaceReactionSystem` which was simulated.
- `filename`: The name of the file to which we wish to save the animation.
- `nframes = 200`: The number of frames in the animation (these are evenly samples across the simulation).
- `framerate = 20`: The frame rate of the animation.
- `ttitle = true`: Whether to add a title showing the simulation's time throughout the animation.

In addition, depending on the type of space used, the following optional arguments might be relevant.

Arguments (1d spaces):
- `markersize = 20`: The size of the markers marking each compartment's value.
- `plot_min = nothing`: The y-scale's minimum. If `nothing`, use the simulation's minimum value.
- `plot_max = nothing`: The y-scale's maximum. If `nothing`, use the simulation's maximum value.

Arguments (Graph & 2d spaces):
- `colormap = :BuGn_7`: The colour map with which we display the species amounts in the animation.
- `plot_min = nothing`: The minimum value for the colour scale (values less than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value.
- `plot_max = nothing`: The maximum value for the colour scale (values more than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value.

Arguments (Graph spaces):
- `node_size = 50`: The size of the compartments in the plot.
- `layout = Spring()`: The layout for the graph nodes in the plot. Can be provided as a vector, where the i'th element is a 2-valued tuple (determining the i'th compartment's y and x positions, respectively).

Notes: 
- For masked spaces, there are no value displayed for grid points which do not correspond to a compartments.
- The current animation interface if a work in progress, and modifications are expected. if you have any feedback, please contact the package authors.
- Additional arguments can be passed to `dspace_animation`, which then will be passed to Makie's `heatmap` plotting command.
"""
function dspace_animation end

"""
    dspace_kymograph(sol, sp, dsrs::DiscreteSpaceReactionSystem, kwargs...)

Creates a kymograph of a `DiscreteSpaceReactionSystem` simulation based on a Cartesian or masked space. 
The plot shows the compartments on the y-axis, and the time development of the system's state along
the x-axis. Species amounts are shown as a heatmap.

Arguments (all spaces):
- `sol`: The simulation we wish to plot.
- `sp`: The species whose values we wish to plot. Can be provided either in its symbolic form or as a symbol.
- `dsrs`: The `DiscreteSpaceReactionSystem` which was simulated.
- `colormap = :BuGn_7`: The colour map with which we display the species amounts in the kymograph.
- `nframes = 200`: The number of time samples which the time series is sampled with.
- `plot_min = nothing`: The minimum value for the colour scale (values less than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value.
- `plot_max = nothing`: The maximum value for the colour scale (values more than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value.

Notes: 
- For masked spaces, there are no value displayed for grid points which do not correspond to a compartments.
- The current plotting interface is a work in progress, and modifications are expected. if you have any feedback, please contact the package authors.
- Additional arguments can be passed to `dspace_plot`, which then will be passed to Makie's `heatmap` plotting command.
"""
function dspace_kymograph end

### Utility Functions ###
# Some of this is used by both extension

# Handles masked spaces by converting sparse arrays to dense arrays where missing values are
# replaced by `NaN` (which Makie can handle well in plotting). The `deepcopy` bit is due to a weird
# bug that I really cannot figure out. The conversion to Float64 is for jump simulations, where
# `NaN`s cannot be presented in Int64 vectors.
function demask_vals(vals, dsrs::DiscreteSpaceReactionSystem)
    has_masked_dspace(dsrs) || return vals
    idxs = findall(!b for b in Catalyst.dspace(dsrs))
    return [demask_vals(vs, idxs) for vs in vals]
end
function demask_vals(v::AbstractSparseArray{T, Int64, 1}, idxs) where {T}
    v = deepcopy(v)
    (eltype(v) <: Int) && (v = Float64.(v))
    foreach(idx -> v[idx] = NaN, idxs)
    return Vector{Float64}(v)
end
function demask_vals(v::AbstractSparseArray{T, Int64, 2}, idxs) where {T}
    v = deepcopy(v)
    (eltype(v) <: Int) && (v = Float64.(v))
    foreach(idx -> v[idx] = NaN, idxs)
    return Matrix(v)
end

# Extract the values from the solution into a value vector (sampled at designated time points).
# Also returns the `plot_min` and `plot_max` values.
function extract_vals(sol, sp, dsrs::DiscreteSpaceReactionSystem, plot_min, plot_max, t)
    vals = spat_getu(sol, sp, dsrs; t)
    isnothing(plot_min) && (plot_min = minimum(minimum(val) for val in vals))
    isnothing(plot_max) && (plot_max = maximum(maximum(val) for val in vals))
    vals = demask_vals(vals, dsrs)
    return vals, plot_min, plot_max
end

# For a Cartesian or masked space, return the grid axis vectors.
function extract_grid_axes(dsrs::DiscreteSpaceReactionSystem)
    x_vals = LinRange(1, grid_size(dsrs)[1], grid_size(dsrs)[1])
    y_vals = LinRange(1, grid_size(dsrs)[2], grid_size(dsrs)[2])
    return x_vals, y_vals
end
