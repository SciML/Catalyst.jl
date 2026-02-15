### Extensions Tests ###
# This file can be run directly for the Extensions workflow,
# or included from runtests.jl for GROUP=="Extensions" or GROUP=="All".

using SafeTestsets, Test

# Note: We use @__DIR__ directly in each include because @safetestset
# runs code in a separate module where outer constants aren't visible.

@time begin
    @time @safetestset "Graph visualization" begin include(joinpath(@__DIR__, "extensions", "graphmakie.jl")) end
    @time @safetestset "BifurcationKit Extension" begin include(joinpath(@__DIR__, "extensions", "bifurcation_kit.jl")) end
    @time @safetestset "HomotopyContinuation Extension" begin include(joinpath(@__DIR__, "extensions", "homotopy_continuation.jl")) end
    @time @safetestset "Structural Identifiability Extension" begin include(joinpath(@__DIR__, "extensions", "structural_identifiability.jl")) end
    @time @safetestset "Steady State Stability Computations" begin include(joinpath(@__DIR__, "extensions", "stability_computation.jl")) end

    # Test spatial plotting, using CairoMakie and GraphMakie
    @time @safetestset "Discrete Space Simulation Plotting" begin include(joinpath(@__DIR__, "extensions", "dspace_simulation_plotting.jl")) end
end
