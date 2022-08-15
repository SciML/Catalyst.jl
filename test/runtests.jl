### Fecth the require packages ###
using SafeTestsets

### Run the tests ###
@time begin

    # system type tests
    @time @safetestset "ReactionSystem Test" begin include("reactionsystem.jl") end
    @time @safetestset "ReactionSystem Components Test" begin include("reactionsystem_components.jl") end

    # Tests all features realted to constructing a model
    @time @safetestset "1.6 Arrows" begin include("newarrows.jl") end
    @time @safetestset "Basic DSL" begin include("dsl.jl") end
    @time @safetestset "Model Construction" begin include("make_model.jl") end
    @time @safetestset "Custom Functions" begin include("custom_functions.jl") end
    @time @safetestset "Model Modification" begin include("model_modification.jl") end
    @time @safetestset "Symbolic Stoichiometry" begin include("symbolic_stoich.jl") end
    @time @safetestset "Conservation Laws" begin include("conslaws.jl") end
    @time @safetestset "Network Properties" begin include("network_properties.jl") end

    # Test api
    @time @safetestset "API" begin include("api.jl") end

    # Tests various core properties of the package.
    @time @safetestset "Higher Order" begin include("higher_order_reactions.jl") end
    @time @safetestset "U0 and Parameters Input Types" begin include("u0_n_parameter_inputs.jl") end

    # Tests related to solving Ordinary Differential Equations.
    @time @safetestset "ODE System Solving" begin include("solve_ODEs.jl") end
    @time @safetestset "Make Jacobian" begin include("make_jacobian.jl") end
    #@time @safetestset "DiffEq Steady State Solving" begin include("steady_state_problems.jl") end

    # Tests related to solving Stochastic Differential Equations.
    @time @safetestset "SDE System Solving" begin include("solve_SDEs.jl") end

    # Tests related to solving Jump Systems.
    @time @safetestset "Jump System Solving" begin include("solve_jumps.jl") end

    # Miscellaneous tests
    #@time @safetestset "Basic Plotting" begin include("plotting.jl") end
    @time @safetestset "Latexify" begin include("latexify.jl") end

    # disable on Macs as can't install GraphViz via jll
    if !Sys.isapple()
        @time @safetestset "Graphs" begin include("graphs.jl") end
    end

    @time @safetestset "Units" begin include("units.jl") end
end # @time
