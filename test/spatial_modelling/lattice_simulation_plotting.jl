### Preparations ###

# Fetch packages.
using Catalyst, CairoMakie, OrdinaryDiffEq, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)


### Checks Basic Plot Cases ###

# Basic test for animations and plots of 1d Cartesian and masked lattices.
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
            @test plt[1].val[1][2] ≈ sol.u[end][1]
            
            # Attempts to animate the simulation (using various arguments). Deletes the saved file.
            lattice_animation(sol, :X, lrs, "animation_tmp.mp4"; nframes = 10, framerate = 10, colormap = :BuGn_6)
            @test isfile("animation_tmp.mp4")
            rm("animation_tmp.mp4")

            # Plots the kymograph and checks that a stored value is correct.
            fig, ax, hm = lattice_kymograph(sol, :X, lrs)
            hm[3].val[end,1] ≈ sol.u[end][1]
        end
    end
end

# Basic test for animations and plots of 2d Cartesian and masked lattices.
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
            @test hm[3].val[1] ≈ sol.u[end][1]
                    
            # Attempts to animate the simulation (using various arguments). Deletes the saved file.
            lattice_animation(sol, :X, lrs, "animation_tmp.mp4"; nframes = 10, framerate = 10, colormap = :BuGn_6)
            @test isfile("animation_tmp.mp4")
            rm("animation_tmp.mp4")
        end
    end
end
