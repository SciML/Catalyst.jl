### Fecth the require packages ###
using SafeTestsets


### Decalre constants realted to groups.
const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = Sys.iswindows() && haskey(ENV,"APPVEYOR")
const is_TRAVIS = haskey(ENV,"TRAVIS")


### Ruyn the tests ###
@time begin

# Tests all features realted to constructing a model
if GROUP == "All" || GROUP == "Network Creation & Modification"
  @time @safetestset "Model Construction" begin include("make_model.jl") end
  @time @safetestset "Custom Functions" begin include("custom_functions.jl") end
  @time @safetestset "Model Modification" begin include("model_modification.jl") end
end

# Tests various core properties of the package.
if GROUP == "All" || GROUP == "Core Properties"
  @time @safetestset "Higher Order" begin include("higher_order_reactions.jl") end
  @time @safetestset "U0 and Parameters Input Types" begin include("u0_n_parameter_inputs.jl") end
  #@time @safetestset "Network Property Queries" begin include("property_query.jl") end
end

# Tests related to solving Ordinary Differential Equations.
if GROUP == "All" || GROUP == "Deterministic Tests"
  @time @safetestset "ODE System Solving" begin include("solve_ODEs.jl") end
  @time @safetestset "Make Jacobian" begin include("make_jacobian.jl") end
end

# Tests related to solving Stochastic Differential Equations.
if GROUP == "All" || GROUP == "SDE Problems"
  @time @safetestset "SDE System Solving" begin include("solve_SDEs.jl") end
end

# Tests related to solvingJump Systems.
if GROUP == "All" || GROUP == "Jump Problems"
  @time @safetestset "Jump System Solving" begin include("solve_jumps.jl") end
end

# Miscellaneous tests.
if GROUP == "All" || GROUP == "Miscellaneous"
  @time @safetestset "Basic Plotting" begin include("plotting.jl") end
  #@time @safetestset "Latexify" begin include("latexify.jl") end
end

end # @time
