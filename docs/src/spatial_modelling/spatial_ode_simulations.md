# [Spatial ODE simulations](@id spatial_dspace_ode_simulations)
Our [introduction to spatial space simulations](@ref spatial_dspace_modelling_intro) has already provided an extensive description of how to simulate [`DiscreteSpaceReactionSystem`](@ref)s using ODEs. Further tutorials have also shown how to [retrieve values from simulations](@ref dspace_simulation_structure_interaction_simulation_species) and or how to [plot them](@ref dspace_simulation_plotting). Here we will build on this, primarily discussing strategies for increasing ODE simulation performance. This is especially important for spatial simulations, as these typically are more computationally demanding as compared to non-spatial ones. While focusing on non-spatial simulations, this [ODE performance tutorial](@ref ode_simulation_performance) is also be useful to read.

## [Solver selection for spatial ODE simulations](@id spatial_dspace_ode_simulations_solvers)
Previously we have described [how to select ODE solvers, and how this can impact simulation performance](@ref ode_simulation_performance_solvers). This is especially relevant for spatial simulations. For stiff problems, `FBDF` is a good first solver to try. For non-stiff problems, `ROCK2` is instead a good first alternative. However, it is still worthwhile to explore a range of alternative solvers.

## [Jacobian options for spatial ODE simulations](@id spatial_dspace_ode_simulations_jacobians)
We have previously described how, when [implicit solvers are used to solve stiff ODEs](@ref ode_simulation_performance_stiffness), the [strategy for computing the system Jacobian](@ref ode_simulation_performance_jacobian) is important. This is especially the case for spatial simulations, where the Jacobian often is large and highly sparse. Catalyst implements special methods for spatial Jacobians. To utilise these, provide the `jac = true` argument to your `ODEProblem` when it is created (if `jac = false`, which is the default, [*automatic differentiation*](https://en.wikipedia.org/wiki/Automatic_differentiation) will be used for Jacobian computation). Here we simulate a [Brusselator](@ref basic_CRN_library_brusselator) while designating to use Catalyst's computed Jacobian:
```@example spatial_ode
using Catalyst, OrdinaryDiffEqBDF
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
diffusion_rx = @transport_reaction D X
space = CartesianGrid((20,20))
dsrs = DiscreteSpaceReactionSystem(brusselator, [diffusion_rx], space)

u0 = [:X => rand(20, 20), :Y => 10.0]
tspan = (0.0, 1.0)
ps = [:A => 1.0, :B => 4.0, :D => 0.2]
oprob = ODEProblem(dsrs, u0, tspan, ps; jac = true)
sol = solve(oprob, FBDF())
nothing # hide
``` 
For large systems, building a dense Jacobian can be problematic, in which case a [*sparse*](@ref ode_simulation_performance_sparse_jacobian) Jacobian can be designated using `sparse = true`: 
```@example spatial_ode
oprob = ODEProblem(dsrs, u0, tspan, ps; jac = true, sparse = true)
sol = solve(oprob, FBDF())
nothing # hide
``` 

It is possible to use `sparse = true` while `jac = false`, in which case a sparse Jacobian is computed using automatic differentiation. 
