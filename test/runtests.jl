### Preparations ###

# Required for `@safetestset` and `@testset`, respectively.
using SafeTestsets, Test, Pkg

# Required for running parallel test groups (copied from ModelingToolkit).
const GROUP = get(ENV, "GROUP", "All")

function activate_extensions_env()
    Pkg.activate("extensions")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
end

### Run Tests ###
@time begin


    # Tests extensions.
    if GROUP == "All" || GROUP == "Extensions"
        activate_extensions_env()

        @time @safetestset "Graph visualization" begin include("extensions/graphmakie.jl") end
        @time @safetestset "BifurcationKit Extension" begin include("extensions/bifurcation_kit.jl") end
        @time @safetestset "HomotopyContinuation Extension" begin include("extensions/homotopy_continuation.jl") end

        # BROKEN
        # @time @safetestset "Structural Identifiability Extension" begin include("extensions/structural_identifiability.jl") end

        # Tests stability computation (but requires the HomotopyContinuation extension).
        #@time @safetestset "Steady State Stability Computations" begin include("extensions/stability_computation.jl") end

        # Test spatial plotting, using CairoMakie and GraphMakie
        @time @safetestset "Lattice Simulation Plotting" begin include("extensions/lattice_simulation_plotting.jl") end
    end

end # @time
