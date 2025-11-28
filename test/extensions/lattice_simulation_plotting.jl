### Preparations ###

# Fetch packages.
using Catalyst, CairoMakie, GraphMakie, Graphs
using JumpProcesses, OrdinaryDiffEqTsit5, Test


### Checks Basic Plot Cases ###

# Basic test for animations and plots of 1d Cartesian and masked lattice simulations.
# Checks both animations and plots.
# Checks for both ODE and jump simulations.
# TODO: Should be expanded and made more comprehensive.
let
    # Creates the `LatticeReactionSystem` model.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    diffusion_rx = @transport_reaction D X
    for lattice in [CartesianGrid(3), [true, true, false]]
        lrs = LatticeReactionSystem(rs, [diffusion_rx], lattice)

        # Simulates the model (using ODE and jumps).
        u0 = [:X => [1, 2, 3]]
        tspan = (0.0, 1.0)
        ps = [:p => 1.0, :d => 1.0, :D => 0.2]
        oprob = ODEProblem(lrs, u0, tspan, ps)
        osol = solve(oprob, Tsit5())
        dprob = DiscreteProblem(lrs, u0, tspan, ps)
        jprob = JumpProblem(lrs, dprob, NSM())
        jsol = solve(jprob, SSAStepper())

        for sol in [osol, jsol]
            # Plots the simulation and checks that a stored value is correct.
            fig, ax, plt = lattice_plot(sol, :X, lrs; t = 1.0)

            # Attempts to animate the simulation (using various arguments). Deletes the saved file.
            lattice_animation(sol, :X, lrs, "animation_tmp.mp4"; nframes = 10, framerate = 10, colormap = :BuGn_6)
            @test isfile("animation_tmp.mp4")
            rm("animation_tmp.mp4")

            # Plots the kymograph and checks that a stored value is correct.
            fig, ax, hm = lattice_kymograph(sol, :X, lrs)
        end
    end
end

# Basic test for animations and plots of 2d Cartesian and masked lattice simulations.
# Checks both animations and plots.
# Checks for both ODE and jump simulations.
# TODO: Should be expanded and made more comprehensive.
let
    # Creates the `LatticeReactionSystem` model.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    diffusion_rx = @transport_reaction D X
    for lattice in [CartesianGrid((2,2)), [true true; false true]]
        lrs = LatticeReactionSystem(rs, [diffusion_rx], lattice)

        # Simulates the model (using ODE and jumps).
        u0 = [:X => [1 2; 3 4]]
        tspan = (0.0, 1.0)
        ps = [:p => 1.0, :d => 1.0, :D => 0.2]
        oprob = ODEProblem(lrs, u0, tspan, ps)
        osol = solve(oprob, Tsit5())
        dprob = DiscreteProblem(lrs, u0, tspan, ps)
        jprob = JumpProblem(lrs, dprob, NSM())
        jsol = solve(jprob, SSAStepper())

        for sol in [osol, jsol]
            # Plots the simulation and checks that a stored value is correct.
            fig, ax, hm = lattice_plot(sol, :X, lrs; t = 1.0)

            # Attempts to animate the simulation (using various arguments). Deletes the saved file.
            lattice_animation(sol, :X, lrs, "animation_tmp.mp4"; nframes = 10, framerate = 10, colormap = :BuGn_6)
            @test isfile("animation_tmp.mp4")
            rm("animation_tmp.mp4")
        end
    end
end

# Basic tests that plotting functions yield correct errors for 3d Cartesian and masked lattice simulations.
let
    rs = @reaction_network begin
        d, X --> 0
    end
    diffusion_rx = @transport_reaction D X
    lattice = CartesianGrid((2,2,2))
    lrs = LatticeReactionSystem(rs, [diffusion_rx], lattice)
    oprob = ODEProblem(lrs, [:X => 1.0], 1.0, [:d => 1.0, :D => 0.2])
    osol = solve(oprob, Tsit5())

    @test_throws ArgumentError lattice_plot(osol, :X, lrs)
    @test_throws ArgumentError lattice_animation(osol, :X, lrs, "animation_tmp.mp4")
end

# Basic test for animations and plots of graph lattices simulations.
# Checks both animations and plots.
# Checks for both ODE and jump simulations.
# TODO: Should be expanded and made more comprehensive.
let
    # Creates the `LatticeReactionSystem` model.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    diffusion_rx = @transport_reaction D X
    lattice = Graphs.SimpleGraphs.cycle_graph(4)
    lrs = LatticeReactionSystem(rs, [diffusion_rx], lattice)

    # Simulates the model (using ODE and jumps).
    u0 = [:X => [1, 2, 3, 4]]
    tspan = (0.0, 1.0)
    ps = [:p => 1.0, :d => 1.0, :D => 0.2]
    oprob = ODEProblem(lrs, u0, tspan, ps)
    osol = solve(oprob, Tsit5())
    dprob = DiscreteProblem(lrs, u0, tspan, ps)
    jprob = JumpProblem(lrs, dprob, NSM())
    jsol = solve(jprob, SSAStepper())

    for sol in [osol, jsol]
        # Plots the simulation and checks that a stored value is correct.
        fig, ax, plt = lattice_plot(sol, :X, lrs; t = 0.0)

        # Attempts to animate the simulation (using various arguments). Deletes the saved file.
        lattice_animation(sol, :X, lrs, "animation_tmp.mp4"; nframes = 10, framerate = 10, colormap = :BuGn_6)
        @test isfile("animation_tmp.mp4")
        rm("animation_tmp.mp4")
    end
end
