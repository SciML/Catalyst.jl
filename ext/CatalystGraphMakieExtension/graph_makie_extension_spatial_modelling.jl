### Graph Lattice Simulation Plots/Animations ###

# Internal dispatch for the plotting of a lattice simulation on a unstructured (graph) lattice. 
function lattice_plot(sol, sp, lrs::LatticeReactionSystem{Q, R, <:AbstractGraph, T};
        t = sol.t[end], plot_min = nothing, plot_max = nothing, colormap = :BuGn_7,
        node_size = 50, kwargs...) where {Q, R, T}
    # Prepares the inputs to the figure.
    plot_graph = SimpleGraph(Catalyst.lattice(lrs))
    _, plot_min, plot_max = Catalyst.extract_vals(sol, sp, lrs, plot_min, plot_max, nothing)
    vals = lat_getu(sol, sp, lrs; t = [t])
    vals = Catalyst.demask_vals(vals, lrs)[1]
    println(kwargs)
    # Creates the figure.
    return graphplot(plot_graph; node_color = vals,
        node_attr = (colorrange = (plot_min, plot_max), colormap), node_size, kwargs...
    )
end

# Internal dispatch for the animation of a lattice simulation on a unstructured (graph) lattice. 
function lattice_animation(
        sol, sp, lrs::LatticeReactionSystem{Q, R, <:AbstractGraph, T}, filename::String;
        t = sol.t[end], nframes = 200, framerate = 20, plot_min = nothing, plot_max = nothing,
        colormap = :BuGn_7, node_size = 50, ttitle::Bool = true, kwargs...) where {Q, R, T}
    # Prepares the inputs to the figure.
    plot_graph = SimpleGraph(Catalyst.lattice(lrs))
    t = LinRange(sol.prob.tspan[1], sol.prob.tspan[2], nframes)
    vals, plot_min, plot_max = Catalyst.extract_vals(sol, sp, lrs, plot_min, plot_max, t)

    # Creates the base figure (which is modified in the animation).
    frame = Makie.Observable(1)
    axis_kwarg = if ttitle
        (; axis = (; title = Makie.@lift(string("Time: ", round(t[$frame]; sigdigits = 3)))))
    else
        (;)
    end
    fig, ax, plt = graphplot(plot_graph; node_color = Makie.@lift(vals[$frame]),
        node_attr = (colorrange = (plot_min, plot_max), colormap), node_size, axis_kwarg..., kwargs...)

    # Creates the animation.
    GraphMakie.record(fig, filename, 1:1:nframes; framerate) do i
        frame[] = i
    end
    return nothing
end
