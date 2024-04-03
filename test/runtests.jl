### Fetch the require packages ###
using SafeTestsets

### Run the tests ###
@time begin
    @time @safetestset "Structural Identifiability Extension" begin include("extensions/structural_identifiability.jl") end
    
end # @time
