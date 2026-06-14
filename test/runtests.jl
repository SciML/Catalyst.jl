### Preparations ###

# Required for `@safetestset` and `@testset`, respectively.
using SafeTestsets, Test, SciMLTesting

### Run Tests ###
@time begin
    run_tests(;
        core = () -> nothing,
        groups = Dict(
            "Modeling" => () -> begin
                # Tests the `ReactionSystem` structure and its properties.
                @time @safetestset "Reaction Structure" begin include(joinpath(@__DIR__, "reactionsystem_core", "reaction.jl")) end
                @time @safetestset "ReactionSystem Structure" begin include(joinpath(@__DIR__, "reactionsystem_core", "reactionsystem.jl")) end
                @time @safetestset "Higher Order Reactions" begin include(joinpath(@__DIR__, "reactionsystem_core", "higher_order_reactions.jl")) end
                @time @safetestset "Symbolic Stoichiometry" begin include(joinpath(@__DIR__, "reactionsystem_core", "symbolic_stoichiometry.jl")) end
                @time @safetestset "Parameter Type Designation" begin include(joinpath(@__DIR__, "reactionsystem_core", "parameter_type_designation.jl")) end
                @time @safetestset "Custom CRN Functions" begin include(joinpath(@__DIR__, "reactionsystem_core", "custom_crn_functions.jl")) end
                @time @safetestset "Coupled CRN/Equation Systems" begin include(joinpath(@__DIR__, "reactionsystem_core", "coupled_equation_crn_systems.jl")) end
                @time @safetestset "Events" begin include(joinpath(@__DIR__, "reactionsystem_core", "events.jl")) end
                @time @safetestset "Functional Parameters" begin include(joinpath(@__DIR__, "reactionsystem_core", "functional_parameters.jl")) end

                # Tests model creation via the @reaction_network DSL.
                @time @safetestset "DSL Basic Model Construction" begin include(joinpath(@__DIR__, "dsl", "dsl_basic_model_construction.jl")) end
                @time @safetestset "DSL Advanced Model Construction" begin include(joinpath(@__DIR__, "dsl", "dsl_advanced_model_construction.jl")) end
                @time @safetestset "DSL Options" begin include(joinpath(@__DIR__, "dsl", "dsl_options.jl")) end

                # Tests compositional and hierarchical modelling.
                @time @safetestset "ReactionSystem Components Based Creation" begin include(joinpath(@__DIR__, "compositional_modelling", "component_based_model_creation.jl")) end # hierarchical modelling broken due to https://github.com/SciML/ModelingToolkit.jl/pull/4101
                @time @safetestset "Brownians, Jumps, and Poissonians Composition" begin include(joinpath(@__DIR__, "compositional_modelling", "brownians_and_jumps_composition.jl")) end

                # Tests various miscellaneous features.
                @time @safetestset "API" begin include(joinpath(@__DIR__, "miscellaneous_tests", "api.jl")) end
                @time @safetestset "Units" begin include(joinpath(@__DIR__, "miscellaneous_tests", "units.jl")) end # `_validate` currently no longer avaiable, awaiting advice.
                @time @safetestset "Compound Species" begin include(joinpath(@__DIR__, "miscellaneous_tests", "compound_macro.jl")) end
                @time @safetestset "Reaction Balancing" begin include(joinpath(@__DIR__, "miscellaneous_tests", "reaction_balancing.jl")) end

                # Tests reaction network analysis features.
                @time @safetestset "Conservation Laws" begin include(joinpath(@__DIR__, "network_analysis", "conservation_laws.jl")) end # Multiple issues. https://github.com/SciML/ModelingToolkit.jl/issues/4102 required to start debugging.
                @time @safetestset "Network Properties" begin include(joinpath(@__DIR__, "network_analysis", "network_properties.jl")) end
                @time @safetestset "CRN Theory" begin include(joinpath(@__DIR__, "network_analysis", "crn_theory.jl")) end
            end,
            "Simulation" => () -> begin
                # Tests ODE, SDE, jump simulations, nonlinear solving, and steady state simulations.
                @time @safetestset "ODE System Simulations" begin include(joinpath(@__DIR__, "simulation_and_solving", "simulate_ODEs.jl")) end
                @time @safetestset "Automatic Jacobian Construction" begin include(joinpath(@__DIR__, "simulation_and_solving", "jacobian_construction.jl")) end
                @time @safetestset "SDE System Simulations" begin include(joinpath(@__DIR__, "simulation_and_solving", "simulate_SDEs.jl")) end
                @time @safetestset "Jump System Simulations" begin include(joinpath(@__DIR__, "simulation_and_solving", "simulate_jumps.jl")) end
                @time @safetestset "Nonlinear and SteadyState System Solving" begin include(joinpath(@__DIR__, "simulation_and_solving", "solve_nonlinear.jl")) end

                # Tests upstream SciML and DiffEq stuff.
                @time @safetestset "MTK Structure Indexing" begin include(joinpath(@__DIR__, "upstream", "mtk_structure_indexing.jl")) end
                @time @safetestset "MTK Problem Inputs" begin include(joinpath(@__DIR__, "upstream", "mtk_problem_inputs.jl")) end # Required to fix lots of these: https://github.com/SciML/ModelingToolkit.jl/issues/4098
            end,
            "Hybrid" => () -> begin
                @time @safetestset "ReactionSystem Hybrid Solvers" begin include(joinpath(@__DIR__, "simulation_and_solving", "hybrid_models.jl")) end
            end,
            "Misc" => () -> begin
                @time @safetestset "ReactionSystem Serialisation" begin include(joinpath(@__DIR__, "miscellaneous_tests", "reactionsystem_serialisation.jl")) end
                # BROKEN
                #@time @safetestset "Latexify" begin include(joinpath(@__DIR__, "visualisation", "latexify.jl")) end # https://github.com/SciML/Catalyst.jl/issues/1352
            end,
            "Spatial" => () -> begin
                # Tests spatial modelling and simulations.
                @time @safetestset "PDE Systems Simulations" begin include(joinpath(@__DIR__, "spatial_modelling", "simulate_PDEs.jl")) end
                @time @safetestset "Spatial Reactions" begin include(joinpath(@__DIR__, "spatial_modelling", "spatial_reactions.jl")) end
                @time @safetestset "Discrete Space Reaction Systems" begin include(joinpath(@__DIR__, "spatial_modelling", "dspace_reaction_systems.jl")) end
                @time @safetestset "Spatial Discrete Space Variants" begin include(joinpath(@__DIR__, "spatial_modelling", "dspace_reaction_systems_space_types.jl")) end
                @time @safetestset "ODE Discrete Space Systems Simulations" begin include(joinpath(@__DIR__, "spatial_modelling", "dspace_reaction_systems_ODEs.jl")) end
                @time @safetestset "Jump Discrete Space Systems Simulations" begin include(joinpath(@__DIR__, "spatial_modelling", "dspace_reaction_systems_jumps.jl")) end
                @time @safetestset "Discrete Space Simulation Structure Interfacing" begin include(joinpath(@__DIR__, "spatial_modelling", "dspace_simulation_struct_interfacing.jl")) end
            end,
            "Extensions" => (; env = joinpath(@__DIR__, "extensions"),
                body = joinpath(@__DIR__, "runtests_extensions.jl")),
        ),
        # Quality assurance (Aqua, JET, ExplicitImports) via the SciMLTesting harness.
        # Runs only under `GROUP=QA` (never as part of `All`), in its own `test/qa`
        # sub-environment so the heavier QA tooling does not constrain the main test
        # environment's resolution.
        qa = (; env = joinpath(@__DIR__, "qa"),
            body = joinpath(@__DIR__, "qa", "qa.jl")),
        all = ["Modeling", "Simulation", "Hybrid", "Misc", "Spatial", "Extensions"],
    )
end # @time
