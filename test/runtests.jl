### Preparations ###

# Fetch packages.
using SafeTestsets

# Designate test groups.
const GROUP = get(ENV, "GROUP", "All")


### Run Tests ###
@time begin

    # Tests the `ReactionSystem` structure and its properties.
    if GROUP == "All" || GROUP == "ReactionSystem"
        @testset "ReactionSystem" begin
            @time @safetestset "Reaction Structure" begin include("reactionsystem_core/reaction_structure.jl") end
            @time @safetestset "ReactionSystem Structure" begin include("reactionsystem_core/reactionsystem_structure.jl") end
            @time @safetestset "Higher Order Reactions" begin include("reactionsystem_core/higher_order_reactions.jl") end
            @time @safetestset "Symbolic Stoichiometry" begin include("reactionsystem_core/symbolic_stoichiometry.jl") end
            @time @safetestset "Parameter Type Designation" begin include("reactionsystem_core/designating_parameter_types.jl") end
            @time @safetestset "Custom CRN Functions" begin include("reactionsystem_core/custom_functions.jl") end
            # @time @safetestset "Coupled CRN/Equation Systems" begin include("reactionsystem_structure/coupled_equation_reaction_systems.jl") end
            @time @safetestset "Events" begin include("reactionsystem_core/events.jl") end
        end
    end

    # Tests model creation via the @reaction_network DSL.
    if GROUP == "All" || GROUP == "DSL"
        @testset "DSL" begin
            @time @safetestset "DSL Model Construction" begin include("dsl/dsl_basic_model_construction.jl") end
            @time @safetestset "DSL Basics" begin include("dsl/dsl_advanced_model_construction.jl") end
            @time @safetestset "DSL Options" begin include("dsl/dsl_options.jl") end
        end
    end

    # Tests compositional and hierarchical modelling.
    if GROUP == "All" || GROUP == "CompositionalModelling"
        @testset "CompositionalModelling" begin
            @time @safetestset "ReactionSystem Components Based Creation" begin include("compositional_modelling/component_based_model_creation.jl") end
        end
    end

    # Tests various miscellaneous features.
    if GROUP == "All" || GROUP == "Miscellaneous"
        @testset "Miscellaneous" begin
            @time @safetestset "API" begin include("miscellaneous_tests/api.jl") end
            @time @safetestset "Compound Species" begin include("miscellaneous_tests/compound_macro.jl") end
            @time @safetestset "Reaction Balancing" begin include("miscellaneous_tests/reaction_balancing.jl") end
            @time @safetestset "Units" begin include("miscellaneous_tests/units.jl") end
        end
    end

    # Tests reaction network analysis features.
    if GROUP == "All" || GROUP == "NetworkAnalysis"
        @testset "NetworkAnalysis" begin
            @time @safetestset "Conservation Laws" begin include("network_analysis/conservation_laws.jl") end
            @time @safetestset "Network Properties" begin include("network_analysis/network_properties.jl") end
        end
    end

    # Tests ODE, SDE, jump simulations, nonlinear solving, and steady state simulations.
    if GROUP == "All" || GROUP == "Simulation"
        @testset "Simulation" begin
            @time @safetestset "ODE System Simulations" begin include("simulation_and_solving/simulate_ODEs.jl") end
            @time @safetestset "Automatic Jacobian Construction" begin include("simulation_and_solving/jacobian_construction.jl") end
            @time @safetestset "SDE System Simulations" begin include("simulation_and_solving/simulate_SDEs.jl") end
            @time @safetestset "Jump System Simulations" begin include("simulation_and_solving/simulate_jumps.jl") end
            @time @safetestset "Nonlinear and SteadyState System Solving" begin include("simulation_and_solving/solve_nonlinear.jl") end
        end
    end

    # Tests upstream SciML and DiffEq stuff.
    if GROUP == "All" || GROUP == "Upstream"
        @testset "Upstream" begin
            @time @safetestset "MTK Structure Indexing" begin include("upstream/mtk_structure_indexing.jl") end
            @time @safetestset "MTK Problem Inputs" begin include("upstream/mtk_problem_inputs.jl") end
        end
    end

    # Tests spatial modelling and simulations.
    if GROUP == "All" || GROUP == "Spatial"
        @testset "Spatial" begin
            @time @safetestset "PDE Systems Simulations" begin include("spatial_modelling/simulate_PDEs.jl") end
            @time @safetestset "Lattice Reaction Systems" begin include("spatial_modelling/lattice_reaction_systems.jl") end
            @time @safetestset "ODE Lattice Systems Simulations" begin include("spatial_modelling/lattice_reaction_systems_ODEs.jl") end
        end
    end

    # Tests network visualization.
    if GROUP == "All" || GROUP == "Visualisation"
        @testset "Visualization" begin
            @time @safetestset "Latexify" begin include("visualisation/latexify.jl") end
            # Disable on Macs as can't install GraphViz via jll
            if !Sys.isapple()
                @time @safetestset "Graphs Visualisations" begin include("visualisation/graphs.jl") end
            end
        end
    end

    # Tests extensions.
    if GROUP == "All" || GROUP == "Extensions"
        @testset "Extensions" begin
            # @time @safetestset "BifurcationKit Extension" begin include("extensions/bifurcation_kit.jl") end
            # @time @safetestset "HomotopyContinuation Extension" begin include("extensions/homotopy_continuation.jl") end
            # @time @safetestset "Structural Identifiability Extension" begin include("extensions/structural_identifiability.jl") end
        end
    end

end # @time
