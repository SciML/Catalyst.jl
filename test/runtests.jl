### Fetch the require packages ###
using SafeTestsets

### Run the tests ###
@time begin

    @time @safetestset "Symbolic Stoichiometry" begin include("miscellaneous_tests/symbolic_stoichiometry.jl") end
    @time @safetestset "NonlinearProblems and Steady State Solving" begin include("miscellaneous_tests/nonlinear_solve.jl") end
    @time @safetestset "Events" begin include("miscellaneous_tests/events.jl") end
    @time @safetestset "Compound species" begin include("miscellaneous_tests/compound_macro.jl") end
    @time @safetestset "Reaction balancing" begin include("miscellaneous_tests/reaction_balancing.jl") end
    @time @safetestset "Units" begin include("miscellaneous_tests/units.jl") end

    ### Reaction network analysis. ###
    @time @safetestset "Conservation Laws" begin include("network_analysis/conservation_laws.jl") end
    @time @safetestset "Network Properties" begin include("network_analysis/network_properties.jl") end

    ### Tests ODE, SDE, PDE, and Gillespie Simulations. ###
    @time @safetestset "ODE System Simulations" begin include("model_simulation/simulate_ODEs.jl") end
    @time @safetestset "Automatic Jacobian Construction" begin include("model_simulation/make_jacobian.jl") end
    @time @safetestset "U0 and Parameters Input Variants" begin include("model_simulation/u0_n_parameter_inputs.jl") end
    @time @safetestset "SDE System Simulations" begin include("model_simulation/simulate_SDEs.jl") end
    @time @safetestset "Jump System Simulations" begin include("model_simulation/simulate_jumps.jl") end

    ### Tests Spatial Network Simulations. ###
    @time @safetestset "PDE Systems Simulations" begin include("spatial_reaction_systems/simulate_PDEs.jl") end
    @time @safetestset "Lattice Reaction Systems" begin include("spatial_reaction_systems/lattice_reaction_systems.jl") end
    # Requires internal rework due to permitting non-default parameter types. However, rework have already been carried out in a separate branch.
    # @time @safetestset "ODE Lattice Systems Simulations" begin include("spatial_reaction_systems/lattice_reaction_systems_ODEs.jl") end

    ### Tests network visualization. ###
    @time @safetestset "Latexify" begin include("visualization/latexify.jl") end
    # @time @safetestset "Basic Plotting" begin include("visualization/plotting.jl") end
    # Disable on Macs as can't install GraphViz via jll
    if !Sys.isapple()
        @time @safetestset "Graphs" begin include("visualization/graphs.jl") end
    end
    
    ### Tests extensions. ###
    @time @safetestset "BifurcationKit Extension" begin include("extensions/bifurcation_kit.jl") end
    @time @safetestset "HomotopyContinuation Extension" begin include("extensions/homotopy_continuation.jl") end
    @time @safetestset "Structural Identifiability Extension" begin include("extensions/structural_identifiability.jl") end
    
end # @time
