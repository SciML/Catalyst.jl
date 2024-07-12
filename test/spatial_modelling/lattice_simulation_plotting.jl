### Preparations ###

# Fetch packages.
using Catalyst, CairoMakie, OrdinaryDiffEq, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)


### 2d Grid Animations ###

# Basic test for animation of 2d Cartesian grid.
# Checks that animation is successfully created. 
let
    # Creates the `LatticeReactionSystem` model.
    rs = @reaction_network begin
        (p,d), 0 <--> X
    end
    diffusion_rx = @transport_reaction D X
    lattice = CartesianGrid((2,2))
    lrs = LatticeReactionSystem(rs, [diffusion_rx], lattice)
    
    # Simulates the model.
    u0 = [:X => rand(rng, 2, 2)]
    tspan = (0.0, 1.0)
    ps = [:p => 1.0, :d => 1.0, :D => 0.2]
    oprob = ODEProblem(lrs, u0, tspan, ps)
    sol = solve(oprob, Tsit5())
    
    # Attempts to animate the simulation (using various arguments). Deletes the saved file.
    lattice_animation(sol, :X, "tmp.mp4", lrs; nframes = 10, framerate = 10, colormap = :BuGn_6)
    @test isfile("tmp.mp4")
    rm("../tmp.jl")
end
