### Lattice Simulation Animations ###

"""
    lattice_animation(sol::ODESolution, sp, filename::String, lrs::LatticeReactionSystem;
        colormap = :BuGn_7, nframes = 200, framerate = 20, plot_max = nothing, kwargs...)

Creates an animation of a `LatticeReactionSystem` simulation based on a 2d Cartesian lattice. The
animation is saved to a file, which name is provided in the `filename` argument.

Arguments:
- `sol`: The simulation we wish to animate.
- `sp`: The species which values we wish to animate.
- `filename`: The name f the file to which we wish to save the animation.
- `lrs`: The `LatticeReactionSystem` which was simulated.
- `colormap = :BuGn_7`: The colormap with which we display the animation.
- `nframes = 200`: The number of frames in the animation.
- `framerate = 20`: The frame rate of the animation.
- `plot_max = nothing`: A selection of the maximum value of the solution. If set to `nothing`, the 
maximum value of the simulation is used. All values in the animation is scaled by this value.
"""
function Catalyst.lattice_animation(sol, sp, filename::String,
        lrs::LatticeReactionSystem; colormap = :BuGn_7, nframes = 200, framerate = 20,
        plot_max = nothing, kwargs...)
    # Error checks.
    if has_graph_lattice(lrs)
        error("The `lattice_animation` function does not currently support animations of simulations based on unstructured (graph) lattices.")
    end
    if has_masked_lattice(lrs)
        error("The `lattice_animation` function does not currently support animations of simulations based on masked (graph) lattices.")
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
    isnothing(plot_max) && (plot_max = maximum(maximum(val) for val in vals))
    vals = [[v / plot_max for v in val] for val in vals]

    # Creates the base figure.
    fig, ax, hm = heatmap(x_vals, y_vals, vals[1]; colormap, kwargs...)

    # Creates the animation.
    record(fig, filename, 1:1:nframes; framerate) do i
        hm[3] = vals[i]
    end
end
