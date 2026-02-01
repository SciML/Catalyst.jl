### Extensions Tests ###
# This file can be run directly for the Extensions workflow,
# or included from runtests.jl for GROUP=="Extensions" or GROUP=="All".

using SafeTestsets, Test

# Use @__DIR__ so paths work whether this file is run directly or included
const EXTENSIONS_TEST_DIR = @__DIR__

@time begin
    @time @safetestset "Graph visualization" begin include(joinpath(EXTENSIONS_TEST_DIR, "extensions", "graphmakie.jl")) end
    @time @safetestset "BifurcationKit Extension" begin include(joinpath(EXTENSIONS_TEST_DIR, "extensions", "bifurcation_kit.jl")) end
    @time @safetestset "HomotopyContinuation Extension" begin include(joinpath(EXTENSIONS_TEST_DIR, "extensions", "homotopy_continuation.jl")) end
    @time @safetestset "Structural Identifiability Extension" begin include(joinpath(EXTENSIONS_TEST_DIR, "extensions", "structural_identifiability.jl")) end
    @time @safetestset "Steady State Stability Computations" begin include(joinpath(EXTENSIONS_TEST_DIR, "extensions", "stability_computation.jl")) end

    # Test spatial plotting, using CairoMakie and GraphMakie
    @time @safetestset "Lattice Simulation Plotting" begin include(joinpath(EXTENSIONS_TEST_DIR, "extensions", "lattice_simulation_plotting.jl")) end
end
