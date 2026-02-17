# [Interfacing with Lattice Problems, Integrators, and Solutions](@id dspace_simulation_structure_interaction)
We have [previously described](@ref simulation_structure_interfacing) how to retrieve species and parameter values stored in non-spatial problems, integrators, and solutions. This section describes similar workflows for simulations based on [`DiscreteSpaceReactionSystem`](@ref)s.

Generally, while for non-spatial systems these operations can typically be done by indexing a structure directly, e.g. through
```julia
sol[:X]
```
there are no equally straightforward interfaces for spatial simulations. Typically, helper functions have to be used, e.g
```julia
spat_getu(sol, :X, dsrs)
```
Furthermore, there are some cases of interfacing which are currently not supported (e.g. updating parameter values in `JumpProblem`s). It is likely that these interfaces will be improved in the future (i.e. by introducing a similar syntax to the current non-spatial one). Finally, we note that many of the functions presented below have not been as extensively optimised for performance as other parts of the Catalyst package. Hence, you should take care when designing workflows which requires using them a large number of times.

!!! note
    Below we will describe various features using ODE simulations as examples. However, all interfaces (unless where else is stated) work identically for jump simulations.

## [Retrieving values from space simulations](@id dspace_simulation_structure_interaction_simulation_species)
Let us consider a simulation of a [`DiscreteSpaceReactionSystem`](@ref):
```@example dspace_struct_interaction_sims
using Catalyst, OrdinaryDiffEqDefault
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
diffusion_rx = @transport_reaction D X1
space = CartesianGrid((2,3))
dsrs = DiscreteSpaceReactionSystem(two_state_model, [diffusion_rx], space)

u0 = [:X1 => [0.0 0.0 0.0; 2.0 2.0 2.0], :X2 => 0.0]
tspan = (0.0, 1.0)
ps = [:k1 => 2.0, :k2 => 1.0, :D => 0.1]
oprob = ODEProblem(dsrs, u0, tspan, ps)
sol = solve(oprob)
nothing # hide
```
To retrieve the values of $X1$ across the simulation we use the `spat_getu` function. It takes three arguments:
- The solution objects from which we wish to retrieve values.
- The species which values we wish to retrieve.
- The [`DiscreteSpaceReactionSystem`](@ref) which was simulated.

```@example dspace_struct_interaction_sims
spat_getu(sol, :X1, dsrs)
```
Here, the output is a vector with $X1$'s value at each simulation time step. How the species's value is represented at each time step depends on the space which was originally used to create the [`DiscreteSpaceReactionSystem`](@ref):
- For Cartesian spaces, an array of the same size as the Cartesian space is used. Each array element corresponds to the species's value in the corresponding compartment.
- For masked spaces, a [sparse array](https://docs.julialang.org/en/v1/stdlib/SparseArrays/) of the same size as the masked space is used. Each filled array element corresponds to the species's value in the corresponding compartment. Unfilled array elements correspond to positions without compartments.
- For unstructured (graph) spaces, vectors are used. The i'th element in the vectors corresponds to the species's value in the i'th compartment.

Unlike for non-spatial simulations, `spat_getu` does not take vector (e.g. `spat_getu(sol, [:X1, :X2], dsrs)`) or symbolic expression (e.g. `spat_getu(sol, [X1 + X2], dsrs)`) inputs. However, it is possible to use symbolic variables as input (e.g. `spat_getu(sol, two_state_model.X1, dsrs)`).

### [Retrieving space simulations values at specific time points](@id dspace_simulation_structure_interaction_simulation_species_ts)
Just like for non-spatial solutions, it is possible to access the simulation's values at designated time points. This is possible even if the simulation did not stop at those specific time points (in which case an interpolated value is returned). To do this, the desired time points to sample are provided as a vector to `spat_getu` using the optional argument `t`. E.g. here we retrieve the simulation's (interpolated) values at time points `0.5` and `0.75`:

```@example dspace_struct_interaction_sims
spat_getu(sol, :X1, dsrs; t = [0.5, 0.75])
```

## [Retrieving and updating species values in problems and integrators](@id dspace_simulation_structure_interaction_prob_int_species)
Let us consider a spatial `ODEProblem`
```@example dspace_struct_interaction_prob_ints
using Catalyst, OrdinaryDiffEqDefault
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
diffusion_rx = @transport_reaction D X1
space = CartesianGrid((2,3))
dsrs = DiscreteSpaceReactionSystem(two_state_model, [diffusion_rx], space)

u0 = [:X1 => [0.0 0.0 0.0; 2.0 2.0 2.0], :X2 => 0.0]
tspan = (0.0, 1.0)
ps = [:k1 => 2.0, :k2 => 1.0, :D => 0.1]
oprob = ODEProblem(dsrs, u0, tspan, ps)
nothing # hide
```
We can retrieve the species values stored in `oprob` using the `spat_getu` function. It uses [identical syntax as for simulations](@ref dspace_simulation_structure_interaction_simulation_species) (except that you cannot specify a time point). However, it returns a single set of species values (while for simulations it returns a vector across different time steps):
```@example dspace_struct_interaction_prob_ints
spat_getu(oprob, :X1, dsrs)
```
Again, the format used corresponds to the space used to create the original [`DiscreteSpaceReactionSystem`](@ref). Here, even if a species has homogeneous values, the full format is used.
```@example dspace_struct_interaction_prob_ints
spat_getu(oprob, :X2, dsrs)
```

For both problems and integrators, species values can be updated using the `spat_setu!` function. It uses a similar syntax as `spat_getu`, but takes a fourth argument which is the new values to use for the designated species:
```@example dspace_struct_interaction_prob_ints
spat_setu!(oprob, :X1, dsrs, [1.0 2.0 3.0; 4.0 5.0 6.0])
```
Here, the same format (which depends on the used space) is used for the species's new values, as which is used when initially designating their initial conditions. I.e. to make $X1$'s initial condition values uniform we can call
```@example dspace_struct_interaction_prob_ints
spat_setu!(oprob, :X1, dsrs, 1.0)
```

!!! note
    It is currently not possible to change a species value at a single compartment only. To do so, you must first retrieve its values across all compartments using `spat_getu`, then modify this at the desired compartment, and then use the modified version as input to `spat_setu!`.

Species values in [integrators](@ref simulation_structure_interfacing_integrators) can be interfaced with using identical syntax as for problems.

## [Retrieving and updating parameter values in problems and integrators](@id dspace_simulation_structure_interaction_prob_int_parameters)
Retrieval and updating of parameter values for problems and integrators works similarly as for species, but with the following differences:
- The `spat_getp` and `spat_setp!` functions are used.
- It is currently not possible to interface with parameter values of `JumpProblem`s and their integrators.
- After parameter values are modified, the `rebuild_spat_internals!` function must be applied before the problem/integrator can be used for further analysis.
- Updating of [edge parameters](@ref spatial_dspace_modelling_intro_simulation_edge_parameters) is limited and uses a different interface. 

Let us consider the spatial `ODEProblem` we previously declared. We can check the value of $k1$ by using `spat_getp`
```@example dspace_struct_interaction_prob_ints
spat_getp(oprob, :k1, dsrs)
```
Next, we can update it using `spat_setp!` (here we also confirm that it now has the updated values):
```@example dspace_struct_interaction_prob_ints
spat_setp!(oprob, :k1, dsrs, [1.0 2.0 3.0; 4.0 5.0 6.0])
spat_getp(oprob, :k1, dsrs)
```

If we now were to simulate `oprob`, the simulation would not take the updated value of $k1$ into account. For our changes to take effect we might first need to call `rebuild_spat_internals!` with `oprob` as an input
```@example dspace_struct_interaction_prob_ints
rebuild_spat_internals!(oprob)
```
There are two different circumstances when `rebuild_spat_internals!` must be called:
- When modifying the value of an [edge parameter](@ref spatial_dspace_modelling_intro_simulation_edge_parameters).
- When changing a parameter from having spatially uniform values to spatially non-uniform values, or the other way around.

Parameter values of integrators can be interfaced with just like for problems (this is primarily relevant when using [*callbacks*](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)). Again, after doing so, the `rebuild_spat_internals!` function might need to be applied to the integrator.


### [Retrieving and updatingedge  parameter values in problems and integrators](@id dspace_simulation_structure_interaction_prob_int_parameters_edge_ps)
The `spat_getp` and `spat_setp!` functions cannot currently be applied to [edge parameters](@ref spatial_dspace_modelling_intro_simulation_edge_parameters). Instead, to access the value of an edge parameter, use 
```@example dspace_struct_interaction_prob_ints
oprob.ps[:D]
```
To update an edge parameter's value, use
```@example dspace_struct_interaction_prob_ints
oprob.ps[:D] = [0.2]
nothing # hide
```
This interface is somewhat limited, and the following aspects should be noted:
- Edge parameter values can only be interfaced with if the edge parameter's value is spatially uniform.
- When accessing an (spatially uniform) edge parameter's value, its single value will be encapsulated in a vector.
- When setting an (spatially uniform) edge parameter's value, you must encapsulate the new value in a vector.

