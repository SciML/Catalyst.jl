### Fetch the require packages ###
using SafeTestsets

### Run the tests ###
@time begin

    ### Tests the properties of ReactionSystems. ###
    @time @safetestset "ReactionSystem Structure" begin include("reactionsystem_structure/reactionsystem.jl") end
    @time @safetestset "Higher Order Reactions" begin include("reactionsystem_structure/higher_order_reactions.jl") end


    ### Tests model creation via the @reaction_network DSL. ###
    @time @safetestset "Basic DSL" begin include("dsl/dsl_basics.jl") end
    @time @safetestset "DSL Model Construction" begin include("dsl/dsl_model_construction.jl") end
    @time @safetestset "Custom CRN Functions" begin include("dsl/custom_functions.jl") end
    @time @safetestset "DSL Options" begin include("dsl/dsl_options.jl") end
    @time @safetestset "1.6 Arrows" begin include("dsl/newarrows.jl") end


    ### Non-DSL model cration and modication. ###
    @time @safetestset "ReactionSystem Components Based Creation" begin include("programmatic_model_creation/component_based_model_creation.jl") end
    @time @safetestset "Programmatic Model Expansion" begin include("programmatic_model_creation/model_modification.jl") end


    ### Reaction network analysis. ###
    @time @safetestset "Conservation Laws" begin include("network_analysis/conservation_laws.jl") end
    @time @safetestset "Network Properties" begin include("network_analysis/network_properties.jl") end


    ### Tests ODE, SDE, PDE, and Gillespie Simulations. ###
    @time @safetestset "ODE System Simulations" begin include("model_simulation/simulate_ODEs.jl") end
    @time @safetestset "Automatic Jacobian Construction" begin include("model_simulation/make_jacobian.jl") end
    @time @safetestset "U0 and Parameters Input Variants" begin include("model_simulation/u0_n_parameter_inputs.jl") end
    # @time @safetestset "DiffEq Steady State Solving" begin include("model_simulation/steady_state_problems.jl") end
    @time @safetestset "SDE System Simulations" begin include("model_simulation/simulate_SDEs.jl") end
    @time @safetestset "PDE Systems Simulations" begin include("model_simulation/simulate_pdes.jl") end
    @time @safetestset "Jump System Simulations" begin include("model_simulation/simulate_jumps.jl") end


    ### Tests network visualization. ###
    @time @safetestset "Latexify" begin include("visualization/latexify.jl") end
    # @time @safetestset "Basic Plotting" begin include("visualization/plotting.jl") end
    # Disable on Macs as can't install GraphViz via jll
    if !Sys.isapple()
        @time @safetestset "Graphs" begin include("visualization/graphs.jl") end
    end


    # Runs various miscellaneous tests.
    @time @safetestset "API" begin include("miscellaneous_tests/api.jl") end
    @time @safetestset "Units" begin include("miscellaneous_tests/units.jl") end
    @time @safetestset "Symbolic Stoichiometry" begin include("dsl/symbolic_stoichiometry.jl") end


end # @time
