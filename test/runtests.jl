### Preparations ###

# Required for `@safetestset` and `@testset`, respectively.
using SafeTestsets, Test

# Required for running parallel test groups (copied from ModelingToolkit).
#const GROUP = get(ENV, "GROUP", "All")


### Run Tests ###
@time begin


        @time @safetestset "Lattice Reaction Systems" begin include("spatial_reaction_systems/lattice_reaction_systems.jl") end
        @time @safetestset "Lattice Reaction Systems Lattice Types" begin include("spatial_reaction_systems/lattice_reaction_systems_lattice_types.jl") end
        @time @safetestset "ODE Lattice Systems Simulations" begin include("spatial_reaction_systems/lattice_reaction_systems_ODEs.jl") end
        @time @safetestset "Jump Lattice Systems Simulations" begin include("spatial_reaction_systems/lattice_reaction_systems_jumps.jl") end

    #if GROUP == "All" || GROUP == "Visualisation-Extensions"
        # Tests network visualisation.
        @time @safetestset "Latexify" begin include("visualisation/latexify.jl") end
        # Disable on Macs as can't install GraphViz via jll
        if !Sys.isapple()
            @time @safetestset "Graphs Visualisations" begin include("visualisation/graphs.jl") end
        end

        # Tests extensions.
        # @time @safetestset "BifurcationKit Extension" begin include("extensions/bifurcation_kit.jl") end
        # @time @safetestset "HomotopyContinuation Extension" begin include("extensions/homotopy_continuation.jl") end
        # @time @safetestset "Structural Identifiability Extension" begin include("extensions/structural_identifiability.jl") end
    #end

end # @time
