### 1d Lattice Simulation Plots/Animations ###

# Internal dispatch for the plotting of a lattice simulation on a 1d lattice (Cartesian or masked). 
function lattice_plot(
        sol, sp, lrs::LatticeReactionSystem{Q, R, <:Catalyst.GridLattice{1, S}, T};
        t = sol.t[end], markersize = 20, kwargs...) where {Q, R, S, T}

    # Prepares and creates the plot.
    vals = lat_getu(sol, sp, lrs; t = [t])
    vals = Catalyst.demask_vals(vals, lrs)[1]
    return scatterlines(vals; axis = (xlabel = "Compartment", ylabel = "$(sp)"),
        markersize = markersize, kwargs...)
end

# Internal dispatch for the animation of a lattice simulation on a 1d lattice (Cartesian or masked). 
function lattice_animation(
        sol, sp, lrs::LatticeReactionSystem{Q, R, <:Catalyst.GridLattice{1, S}, T},
        filename::String;
        markersize = 20, plot_min = nothing, plot_max = nothing, nframes = 200, framerate = 20,
        ttitle = true, kwargs...) where {Q, R, S, T}

    # Prepares the inputs to the figure.
    t = LinRange(sol.prob.tspan[1], sol.prob.tspan[2], nframes)
    vals, plot_min, plot_max = Catalyst.extract_vals(sol, sp, lrs, plot_min, plot_max, t)

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
    nothing
end

# Internal dispatch for the kymographs of a lattice simulation on a 1d lattice (Cartesian or masked). 
function lattice_kymograph(
        sol, sp, lrs::LatticeReactionSystem{Q, R, <:Catalyst.GridLattice{1, S}, T};
        colormap = :BuGn_7,
        plot_min = nothing, plot_max = nothing, nframes = 200, kwargs...) where {Q, R, S, T}

    # Prepares the inputs to the figure.
    t = LinRange(sol.prob.tspan[1], sol.prob.tspan[2], nframes)
    vals, plot_min, plot_max = Catalyst.extract_vals(sol, sp, lrs, plot_min, plot_max, t)
    vals = hcat(vals...)'

    # Creates the figure.
    y_vals = LinRange(1, grid_size(lrs)[1], grid_size(lrs)[1])
    return heatmap(t,
        y_vals,
        vals;
        axis = (xlabel = "Time", ylabel = "Compartment",
            xgridvisible = false, ygridvisible = false),
        colormap,
        colorrange = (plot_min, plot_max))
end

### 2d Lattice Simulation Plots/Animations ###

# Internal dispatch for the plotting of a lattice simulation on a 2d lattice (Cartesian or masked). 
function lattice_plot(sol, sp,
        lrs::LatticeReactionSystem{Q, R, <:Catalyst.GridLattice{2, S}, T}; t = sol.t[end],
        colormap = :BuGn_7, plot_min = nothing, plot_max = nothing,
        kwargs...) where {Q, R, S, T}

    # Prepares the inputs to the figure (the `extract_vals` call only finds limits).
    _, plot_min, plot_max = Catalyst.extract_vals(sol, sp, lrs, plot_min, plot_max, nothing)
    vals = lat_getu(sol, sp, lrs; t = [t])
    vals = Catalyst.demask_vals(vals, lrs)[1]
    x_vals, y_vals = Catalyst.extract_grid_axes(lrs)

    # Creates the figure.
    return heatmap(x_vals,
        y_vals,
        vals;
        axis = (xlabel = "Time", ylabel = "Compartment",
            xgridvisible = false, ygridvisible = false),
        colormap,
        colorrange = (plot_min, plot_max),
        kwargs...)
end

# Internal dispatch for the animation of a lattice simulation on a 2d lattice (Cartesian or masked). 
function lattice_animation(
        sol, sp, lrs::LatticeReactionSystem{Q, R, <:Catalyst.GridLattice{2, S}, T},
        filename::String;
        colormap = :BuGn_7, nframes = 200, framerate = 20, plot_min = nothing,
        plot_max = nothing, ttitle = true, kwargs...) where {Q, R, S, T}

    # Prepares the inputs to the figure.
    t = LinRange(sol.prob.tspan[1], sol.prob.tspan[2], nframes)
    vals, plot_min, plot_max = Catalyst.extract_vals(sol, sp, lrs, plot_min, plot_max, t)
    x_vals, y_vals = Catalyst.extract_grid_axes(lrs)

    # Creates the base figure (which is modified in the animation).
    fig, ax, hm = heatmap(x_vals, y_vals, vals[1];
        axis = (xgridvisible = false, ygridvisible = false),
        colormap, colorrange = (plot_min, plot_max),
        kwargs...)
    ttitle && (ax.title = "Time: $(round(t[1]; sigdigits = 3))")

    # Creates the animation.
    record(fig, filename, 1:1:nframes; framerate) do i
        ttitle && (ax.title = "Time: $(round(t[i]; sigdigits = 3))")
        hm[3] = vals[i]
    end
    nothing
end

### 3d Lattice Simulation Plots/Animations (Errors Only) ###

# Internal dispatch for the plotting of a lattice simulation on a 3d lattice (Cartesian or masked). 
function lattice_plot(
        sol, sp, lrs::LatticeReactionSystem{Q, R, <:Catalyst.GridLattice{3, S}, T};
        kwargs...) where {Q, R, S, T}
    throw(ArgumentError("The `lattice_plot` function does not support 3d Cartesian/masked lattices."))
end

# Internal dispatch for the animation of a lattice simulation on a 3d lattice (Cartesian or masked). 
function lattice_animation(
        sol, sp, lrs::LatticeReactionSystem{Q, R, <:Catalyst.GridLattice{3, S}, T},
        filename::String; kwargs...) where {Q, R, S, T}
    throw(ArgumentError("The `lattice_animation` function does not support 3d Cartesian/masked lattices."))
end
