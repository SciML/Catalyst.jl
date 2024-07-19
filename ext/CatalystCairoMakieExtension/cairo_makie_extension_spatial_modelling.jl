### Lattice Master Functions ###

"""
    lattice_plot(sol, sp, lrs::LatticeReactionSystem, filename::String; t = sol.tspan[2], kwargs...)

Creates a plot of a `LatticeReactionSystem` simulation based on a Cartesian or masked lattice. 
The plot is created at the time point specified by `t` (defaults to the simulation's final time point).

Arguments (all lattices):
- `sol`: The simulation we wish to plot.
- `sp`: The species whose values we wish to plot. Can be provided either in its symbolic form or as a symbol.
- `lrs`: The `LatticeReactionSystem` which was simulated.
- `t = sol.t[end]`: The time point at which we wish to plot the solution

In addition, depending on whether a 1d or 2d lattice was used, the following optional arguments
might be relevant.

Arguments (1d lattices):
- `markersize = 20`: The size of the markers marking each compartment's value.

Arguments (2d lattices):
- `colormap = :BuGn_7`: The colour map with which we display the species amounts in the animation.
- `plot_min = nothing`: The minimum value for the colour scale (values less than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value (across the entire simulation, not just at the plotted time value).
- `plot_max = nothing`: The maximum value for the colour scale (values more than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value (across the entire simulation, not just at the plotted time value).

Notes: 
- For masked lattices, there are no value displayed for grid points which do not correspond to a compartments.
- The current plotting interface is a work in progress, and modifications are expected. if you have any feedback, please contact the package authors.
- Additional arguments can be passed to `lattice_plot`, which then will be passed to Makie's `lines` plotting command.
"""
function lattice_plot(sol, sp, lrs::LatticeReactionSystem; t = sol.t[end],
        kwargs...)
    if has_graph_lattice(lrs)
        return lattice_plot_graph(sol, sp, lrs; t, kwargs...)
    elseif grid_dims(lrs) == 1
        return lattice_plot_1d(sol, sp, lrs; t, kwargs...)
    elseif grid_dims(lrs) == 2
        return lattice_plot_2d(sol, sp, lrs; t, kwargs...)
    else
        throw(ArgumentError("The lattice which the input lattice simulation is based on is currently not supported by `lattice_plot`. 1d, 2d (masked and Cartesian) and graph lattices are currently the only supported lattice types."))
    end
end

"""
    lattice_animation(sol, sp, lrs::LatticeReactionSystem, filename::String; kwargs...)

Creates an animation of a `LatticeReactionSystem` simulation based on a Cartesian or masked lattice. 
The animation is saved to a file, whose name is provided in the `filename` argument.

Arguments (all lattices):
- `sol`: The simulation we wish to animate.
- `sp`: The species which values we wish to animate. Can be provided either in its symbolic form or as a symbol.
- `lrs`: The `LatticeReactionSystem` which was simulated.
- `filename`: The name of the file to which we wish to save the animation.
- `nframes = 200`: The number of frames in the animation (these are evenly samples across the simulation).
- `framerate = 20`: The frame rate of the animation.
- `ttitle = true`: Whether to add a title showing the simulation's time throughout the animation.

In addition, depending on whether a 1d or 2d lattice was used, the following optional arguments
might be relevant.

Arguments (1d lattices):
- `markersize = 20`: The size of the markers marking each compartment's value.
- `plot_min = nothing`: The y-scale's minimum. If `nothing`, use the simulation's minimum value.
- `plot_max = nothing`: The y-scale's maximum. If `nothing`, use the simulation's maximum value.

Arguments (2d lattices):
- `colormap = :BuGn_7`: The colour map with which we display the species amounts in the animation.
- `plot_min = nothing`: The minimum value for the colour scale (values less than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value.
- `plot_max = nothing`: The maximum value for the colour scale (values more than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value.

Notes: 
- For masked lattices, there are no value displayed for grid points which do not correspond to a compartments.
- The current animation interface if a work in progress, and modifications are expected. if you have any feedback, please contact the package authors.
- Additional arguments can be passed to `lattice_animation`, which then will be passed to Makie's `heatmap` plotting command.
"""
function lattice_animation(sol, sp, lrs::LatticeReactionSystem, filename::String;
        kwargs...)
    if has_graph_lattice(lrs)
        return lattice_animation_graph(sol, sp, lrs, filename; kwargs...)
    elseif grid_dims(lrs) == 1
        return lattice_animation_1d(sol, sp, lrs, filename; kwargs...)
    elseif grid_dims(lrs) == 2
        return lattice_animation_2d(sol, sp, lrs, filename; kwargs...)
    else
        throw(ArgumentError("The lattice which the input lattice simulation is based on is currently not supported by `lattice_plot`. 1d, 2d (masked and Cartesian) and graph lattices are currently the only supported lattice types."))
    end
end

### 1d Lattice Simulation Plots/Animations ###

# Internal function which handles the plotting of a lattice simulation on a 1d lattice (Cartesian
# or graph). Here, unlike for 2d and graph lattices, `sp` can be a vector of species (in which case
# each on is displayed).
function lattice_plot_1d(sol, sp, lrs::LatticeReactionSystem; t = nothing, markersize = 20, kwargs...)
    vals = lat_getu(sol, sp, lrs; t = [t])
    vals = demask_vals(vals, lrs)[1]
    scatterlines(vals; axis = (xlabel = "Compartment", ylabel = "$(sp)"),
        markersize = markersize, kwargs...)
end

# Internal function which handles the animation of a lattice simulation on a 1d lattice (Cartesian
# or graph). Here, unlike for 2d and graph lattices, `sp` can be a vector of species (in which case
# each on is displayed).
function lattice_animation_1d(sol, sp, lrs::LatticeReactionSystem, filename::String;
        markersize = 20, plot_min = nothing, plot_max = nothing, nframes = 200, framerate = 20,
        ttitle = true, kwargs...)

    # Prepares the inputs to the figure.
    t = LinRange(sol.prob.tspan[1], sol.prob.tspan[2], nframes)
    vals, plot_min, plot_max = extract_vals(sol, sp, lrs, plot_min, plot_max, t)

    # Creates the base figure (which is modified in the animation).
    fig, ax, plt = scatterlines(vals[1]; 
        axis = (xlabel = "Compartment", ylabel = "$(sp)", 
        limits = (nothing, nothing, plot_min, plot_max)),
        markersize = markersize, kwargs...)
    ttitle && (ax.title = "Time: $(round(t[1]; sigdigits = 3))")

    # Creates the animation.
    record(fig, filename, 1:1:nframes; framerate) do i
        for vertex in 1:grid_size(lrs)[1]
            plt[1].val[vertex] = [vertex, vals[i][vertex]]
        end
        ttitle && (ax.title = "Time: $(round(t[i]; sigdigits = 3))")
    end
end

"""
    lattice_kymograph(sol, sp, lrs::LatticeReactionSystem, kwargs...)

Creates a kymograph of a `LatticeReactionSystem` simulation based on a Cartesian or masked lattice. 
The plot shows the compartments on the y-axis, and the time development of the system's state along
the x-axis. Species amounts are shown as a heatmap.

Arguments (all lattices):
- `sol`: The simulation we wish to plot.
- `sp`: The species whose values we wish to plot. Can be provided either in its symbolic form or as a symbol.
- `lrs`: The `LatticeReactionSystem` which was simulated.
- `colormap = :BuGn_7`: The colour map with which we display the species amounts in the kymograph.
- `nframes = 200`: The number of time samples which the time series is sampled with.
- `plot_min = nothing`: The minimum value for the colour scale (values less than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value.
- `plot_max = nothing`: The maximum value for the colour scale (values more than this will be set at this value when the colour scale is computed). If `nothing`, use the simulation's minimum value.

Notes: 
- For masked lattices, there are no value displayed for grid points which do not correspond to a compartments.
- The current plotting interface is a work in progress, and modifications are expected. if you have any feedback, please contact the package authors.
- Additional arguments can be passed to `lattice_plot`, which then will be passed to Makie's `heatmap` plotting command.
"""
function lattice_kymograph(sol, sp, lrs::LatticeReactionSystem; colormap = :BuGn_7,
        plot_min = nothing, plot_max = nothing, nframes = 200, kwargs...)
    # Prepares the inputs to the figure.
    t = LinRange(sol.prob.tspan[1], sol.prob.tspan[2], nframes)
    vals, plot_min, plot_max = extract_vals(sol, sp, lrs, plot_min, plot_max, t)
    vals = hcat(vals...)'

    # Creates the figure.
    y_vals = LinRange(1, grid_size(lrs)[1], grid_size(lrs)[1])
    fig, ax, hm = heatmap(t, y_vals, vals;
        axis = (xlabel = "Time", ylabel = "Compartment", xgridvisible = false,
        ygridvisible = false), colormap, colorrange = (plot_min, plot_max))
end

### 2d Lattice Simulation Plots/Animations ###

# Internal function which handles the plotting of a lattice simulation on a 2d lattice (Cartesian
# or graph).
function lattice_plot_2d(sol, sp, lrs::LatticeReactionSystem; t = nothing, 
        colormap = :BuGn_7, plot_min = nothing, plot_max = nothing, kwargs...)
    # Prepares the inputs to the figure (the `extract_vals` call only finds limits).
    _, plot_min, plot_max = extract_vals(sol, sp, lrs, plot_min, plot_max, nothing)
    vals = lat_getu(sol, sp, lrs; t = [t])
    vals = demask_vals(vals, lrs)[1]
    x_vals, y_vals = extract_grid_axes(lrs)

    # Creates the figure.
    fig, ax, hm = heatmap(x_vals, y_vals, vals; 
        axis = (xlabel = "Time", ylabel = "Compartment", xgridvisible = false,
        ygridvisible = false), colormap, colorrange = (plot_min, plot_max), kwargs...)
end

# Internal function which handles the animation of a lattice simulation on a 2d lattice (Cartesian
# or graph).
function lattice_animation_2d(sol, sp, lrs::LatticeReactionSystem, filename::String; 
        colormap = :BuGn_7,  nframes = 200, framerate = 20, plot_min = nothing, 
        plot_max = nothing, ttitle = true, kwargs...)

    # Prepares the inputs to the figure.
    t = LinRange(sol.prob.tspan[1], sol.prob.tspan[2], nframes)
    vals, plot_min, plot_max = extract_vals(sol, sp, lrs, plot_min, plot_max, t)
    x_vals, y_vals = extract_grid_axes(lrs)

    # Creates the base figure (which is modified in the animation).
    fig, ax, hm = heatmap(x_vals, y_vals, vals[1]; 
        axis = (xgridvisible = false, ygridvisible = false) ,
        colormap, colorrange = (plot_min, plot_max),
        kwargs...)
    ttitle && (ax.title = "Time: $(round(t[1]; sigdigits = 3))")

    # Creates the animation.
    record(fig, filename, 1:1:nframes; framerate) do i
        ttitle && (ax.title = "Time: $(round(t[i]; sigdigits = 3))")
        hm[3] = vals[i]
    end
end

### Graph Lattice Simulation Plots/Animations ###

# Internal function which handles the plotting of a lattice simulation on a graph lattice.
# Takes the additional, required, kwarg: `vert_positions`, which determines the positions of the
# vertices in the graph plot.
# TODO: Will probably require a separate GraphMakie extension.
function lattice_plot_graph(sol, sp, lrs::LatticeReactionSystem; kwargs...)
    error("Graph lattice plotting not supported yet")
end

# Internal function which handles the animation of a lattice simulation on a graph lattice.
# Takes the additional, required, kwarg: `vert_positions`, which determines the positions of the
# vertices in the graph plot.
# TODO: Will probably require a separate GraphMakie extension.
function lattice_animation_graph(sol, sp, lrs::LatticeReactionSystem, filename::String;
        kwargs...)
    error("Graph lattice animations not supported yet")
end

### Utility Functions ###

# Handles masked lattices by converting sparse arrays to dense arrays where missing values are
# replaced by `NaN` (which Makie can handle well in plotting). The `deepcopy` bit is due to a weird
# bug that I really cannot figure out. The conversion to Float64 is for jump simulations, where
# `NaN`s cannot be presented in Int64 vectors.
function demask_vals(vals, lrs::LatticeReactionSystem)
    has_masked_lattice(lrs) || return vals
    idxs = findall(!b for b in Catalyst.lattice(lrs))
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
function extract_vals(sol, sp, lrs::LatticeReactionSystem, plot_min, plot_max, t)
    vals = lat_getu(sol, sp, lrs; t)
    isnothing(plot_min) && (plot_min = minimum(minimum(val) for val in vals))
    isnothing(plot_max) && (plot_max = maximum(maximum(val) for val in vals))
    vals = demask_vals(vals, lrs)
    return vals, plot_min, plot_max
end

# For a Cartesian or masked lattice, return the grid axis vectors.
function extract_grid_axes(lrs::LatticeReactionSystem)
    x_vals = LinRange(1, grid_size(lrs)[1], grid_size(lrs)[1])
    y_vals = LinRange(1, grid_size(lrs)[2], grid_size(lrs)[2])
    return x_vals, y_vals
end
