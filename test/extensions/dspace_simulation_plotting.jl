### Preparations ###

# Fetch packages.
using Catalyst, CairoMakie, GraphMakie, Graphs
using JumpProcesses, OrdinaryDiffEqTsit5, Test


### Checks Basic Plot Cases ###

# Basic test for animations and plots of 1d Cartesian and masked space simulations.
# Checks both animations and plots.
# Checks for both ODE and jump simulations.
# TODO: Should be expanded and made more comprehensive.
let
    # Creates the `DiscreteSpaceReactionSystem` model.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    diffusion_rx = @transport_reaction D X
    for space in [CartesianGrid(3), [true, true, false]]
        dsrs = DiscreteSpaceReactionSystem(rs, [diffusion_rx], space)

        # Simulates the model (using ODE and jumps).
        u0 = [:X => [1, 2, 3]]
        tspan = (0.0, 1.0)
        ps = [:p => 1.0, :d => 1.0, :D => 0.2]
        oprob = ODEProblem(dsrs, u0, tspan, ps)
        osol = solve(oprob, Tsit5())
        jprob = JumpProblem(dsrs, u0, tspan, ps)
        jsol = solve(jprob, SSAStepper())

        for sol in [osol, jsol]
            # Plots the simulation and checks that a stored value is correct.
            fig, ax, plt = dspace_plot(sol, :X, dsrs; t = 1.0)
            @test_broken plt[1].val[1][2] ≈ sol.u[end][1] # Interface for accessing internals in makie plots have changed. Need to update test.

            # Attempts to animate the simulation (using various arguments). Deletes the saved file.
            dspace_animation(sol, :X, dsrs, "animation_tmp.mp4"; nframes = 10, framerate = 10, colormap = :BuGn_6)
            @test isfile("animation_tmp.mp4")
            rm("animation_tmp.mp4")

            # Plots the kymograph and checks that a stored value is correct.
            fig, ax, hm = dspace_kymograph(sol, :X, dsrs)
            @test_broken hm[3].val[end,1] ≈ sol.u[end][1] # Interface for accessing internals in makie plots have changed. Need to update test.
        end
    end
end

# Basic test for animations and plots of 2d Cartesian and masked space simulations.
# Checks both animations and plots.
# Checks for both ODE and jump simulations.
# TODO: Should be expanded and made more comprehensive.
let
    # Creates the `DiscreteSpaceReactionSystem` model.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    diffusion_rx = @transport_reaction D X
    for space in [CartesianGrid((2,2)), [true true; false true]]
        dsrs = DiscreteSpaceReactionSystem(rs, [diffusion_rx], space)

        # Simulates the model (using ODE and jumps).
        u0 = [:X => [1 2; 3 4]]
        tspan = (0.0, 1.0)
        ps = [:p => 1.0, :d => 1.0, :D => 0.2]
        oprob = ODEProblem(dsrs, u0, tspan, ps)
        osol = solve(oprob, Tsit5())
        jprob = JumpProblem(dsrs, u0, tspan, ps)
        jsol = solve(jprob, SSAStepper())

        for sol in [osol, jsol]
            # Plots the simulation and checks that a stored value is correct.
            fig, ax, hm = dspace_plot(sol, :X, dsrs; t = 1.0)
            @test_broken hm[3].val[1] ≈ sol.u[end][1] # Interface for accessing internals in makie plots have changed. Need to update test.

            # Attempts to animate the simulation (using various arguments). Deletes the saved file.
            dspace_animation(sol, :X, dsrs, "animation_tmp.mp4"; nframes = 10, framerate = 10, colormap = :BuGn_6)
            @test isfile("animation_tmp.mp4")
            rm("animation_tmp.mp4")
        end
    end
end

# Basic tests that plotting functions yield correct errors for 3d Cartesian and masked space simulations.
let
    rs = @reaction_network begin
        d, X --> 0
    end
    diffusion_rx = @transport_reaction D X
    space = CartesianGrid((2,2,2))
    dsrs = DiscreteSpaceReactionSystem(rs, [diffusion_rx], space)
    oprob = ODEProblem(dsrs, [:X => 1.0], 1.0, [:d => 1.0, :D => 0.2])
    osol = solve(oprob, Tsit5())

    @test_throws ArgumentError dspace_plot(osol, :X, dsrs)
    @test_throws ArgumentError dspace_animation(osol, :X, dsrs, "animation_tmp.mp4")
end

# Basic test for animations and plots of graph space simulations.
# Checks both animations and plots.
# Checks for both ODE and jump simulations.
# TODO: Should be expanded and made more comprehensive.
let
    # Creates the `DiscreteSpaceReactionSystem` model.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    diffusion_rx = @transport_reaction D X
    space = Graphs.SimpleGraphs.cycle_graph(4)
    dsrs = DiscreteSpaceReactionSystem(rs, [diffusion_rx], space)

    # Simulates the model (using ODE and jumps).
    u0 = [:X => [1, 2, 3, 4]]
    tspan = (0.0, 1.0)
    ps = [:p => 1.0, :d => 1.0, :D => 0.2]
    oprob = ODEProblem(dsrs, u0, tspan, ps)
    osol = solve(oprob, Tsit5())
    jprob = JumpProblem(dsrs, u0, tspan, ps)
    jsol = solve(jprob, SSAStepper())

    for sol in [osol, jsol]
        # Plots the simulation and checks that a stored value is correct.
        fig, ax, plt = dspace_plot(sol, :X, dsrs; t = 0.0)
        @test plt.node_color[] == osol.u[1]

        # Attempts to animate the simulation (using various arguments). Deletes the saved file.
        dspace_animation(sol, :X, dsrs, "animation_tmp.mp4"; nframes = 10, framerate = 10, colormap = :BuGn_6)
        @test isfile("animation_tmp.mp4")
        rm("animation_tmp.mp4")
    end
end
