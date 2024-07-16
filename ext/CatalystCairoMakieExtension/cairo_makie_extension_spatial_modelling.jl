### Lattice Simulation Plotting ###

### Lattice Simulation Animations ###

"""
    lattice_animation(sol::ODESolution, sp, lrs::LatticeReactionSystem, filename::String;
        colormap = :BuGn_7, nframes = 200, framerate = 20, plot_max = nothing, kwargs...)

Creates an animation of a `LatticeReactionSystem` simulation based on a 2d Cartesian lattice. The
animation is saved to a file, which name is provided in the `filename` argument.

Arguments:
- `sol`: The simulation we wish to animate.
- `sp`: The species which values we wish to animate.
- `filename`: The name f the file to which we wish to save the animation.
- `lrs`: The `LatticeReactionSystem` which was simulated.
- `colormap = :BuGn_7`: The color map with which we display the animation.
- `nframes = 200`: The number of frames in the animation.
- `framerate = 20`: The frame rate of the animation.
- `plot_min = nothing`: The minimum value for the color scale. Each value in the simulation will be 
rescaled according to ` max(0.0, min(plot_max - plot_min, val - plot_min))/(plot_max-plot_min)`, 
where `plot_min` is the minimum value of the scale. If `plot_min = nothing`, it is set to the minimal 
value that appears in the simulation (for our species of interest).
- `plot_max = nothing`: The maximum value for the color scale. Each value in the simulation will be 
rescaled according to ` max(0.0, min(plot_max - plot_min, val - plot_min))/(plot_max-plot_min)`, 
where `plot_max` is the maximum value of the scale. If `plot_max = nothing`, it is set to the maximum 
value that appears in the simulation (for our species of interest).
"""
function Catalyst.lattice_animation(sol, sp, lrs::LatticeReactionSystem, filename::String;
        colormap = :BuGn_7, nframes = 200, framerate = 20, plot_min = nothing, 
        plot_max = nothing, kwargs...)
    # Error checks.
    if has_graph_lattice(lrs)
        error("The `lattice_animation` function does not currently support animations of simulations based on unstructured (graph) lattices.")
    end
    if grid_dims(lrs) != 2
        error("The `lattice_animation` function does not currently support animations of simulations based on lattices of dimensions other than 2.")
    end

    # Gets the values to plot, and the axes to plot them along.
    vals = lat_getu(sol, sp, lrs;
        t = LinRange(sol.prob.tspan[1], sol.prob.tspan[2], nframes))
    x_vals = LinRange(1, grid_size(lrs)[1], grid_size(lrs)[1])
    y_vals = LinRange(1, grid_size(lrs)[2], grid_size(lrs)[2])

    # Rescales all values by the `plot_max` value.
    isnothing(plot_min) && (plot_min = minimum(minimum(val) for val in vals))
    isnothing(plot_max) && (plot_max = maximum(maximum(val) for val in vals))
    has_masked_lattice(lrs) && (vals = [densify_vals(smat, plot_min) for smat in vals])
    vals = [[scale_val(v, plot_min, plot_max) for v in val] for val in vals]

    # Creates the base figure.
    fig, ax, hm = heatmap(x_vals, y_vals, vals[1]; colormap, kwargs...)
    
    # Creates the animation.
    record(fig, filename, 1:1:nframes; framerate) do i
        hm[3] = vals[i]
    end
end


### Lattice Simulation Kymographs ###



### Graph Lattice Simulation Plotting & Animation ###

### Utility Functions ###

# If we have a masked lattice, sets empty lattice points to the minimum value of the heatmap.
# TODO: Change so that we instead can support an entirely different value outside of the heatmap
# colour range.
function densify_vals(smat, plot_min)
    mat = Matrix(smat)
    
    # Utilises the fact that if any matrix values are zero not due to sparsity, then the minimum
    # value of `plot_min` is 0 anyway, so cannot introduce an error by changing these (if `plot_min`)
    # is larger, zero values will be set to `plot_min` later anyway.
    foreach(i -> (mat[i] == 0) && (mat[i] = plot_min), eachindex(mat))
    return mat
end

# Rescales a value between a given maximum and minimum value.
function scale_val(val, min_val, max_val)
    return max(0.0, min(max_val - min_val, val - min_val))/(max_val - min_val)
end