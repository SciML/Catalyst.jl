### Fecth the require packages ###
using SafeTestsets

### Run the tests ###
@time begin

# the following can't really be run until there is an artifact for Graphviz
@time @safetestset "Graphs" begin include("graphs.jl") end

# Tests all features realted to constructing a model
@time @safetestset "Model Construction" begin include("make_model.jl") end
@time @safetestset "Custom Functions" begin include("custom_functions.jl") end
@time @safetestset "Model Modification" begin include("model_modification.jl") end

# Test api
@time @safetestset "API" begin include("api.jl") end

# Tests various core properties of the package.
@time @safetestset "Higher Order" begin include("higher_order_reactions.jl") end
@time @safetestset "U0 and Parameters Input Types" begin include("u0_n_parameter_inputs.jl") end

# Tests related to solving Ordinary Differential Equations.
@time @safetestset "ODE System Solving" begin include("solve_ODEs.jl") end
@time @safetestset "Make Jacobian" begin include("make_jacobian.jl") end
@time @safetestset "DiffEq Steady State Solving" begin include("steady_state_problems.jl") end

# Tests related to solving Stochastic Differential Equations.
@time @safetestset "SDE System Solving" begin include("solve_SDEs.jl") end

# Tests related to solvingJump Systems.
@time @safetestset "Jump System Solving" begin include("solve_jumps.jl") end

# Miscellaneous tests
#@time @safetestset "Basic Plotting" begin include("plotting.jl") end
@time @safetestset "Latexify" begin include("latexify.jl") end

end # @time
