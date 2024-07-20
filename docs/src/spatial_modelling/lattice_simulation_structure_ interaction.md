# [Interfacing with Lattice Problems, Integrators, and Solutions](@id lattice_simulation_structure_interaction)
We have [previously described](@ref simulation_structure_interfacing) how to retrieve species and parameter values stored in non-spatial problems, integrators, and solutions. This section describes similar workflows for simulations based on [`LatticeReactionSystem`](@ref)s.

Generally, while for non-spatial systems these operations can typically be done by indexing a structure directly, e.g. through
```julia
sol[:X]
```
there are no equally straightforward interfaces for spatial simulations. Typically, helper functions have to be used, e.g
```julia
lat_getu(sol, :X, lrs)
```
Furthermore, there are some cases of interfacing which are currently not supported (e.g. updating parameter values in `JumpProblem`s). It is likely that these interfaces will be improved in the future (i.e. by introducing a similar syntax to the current non-spatial one). Finally, we note that many of the functions presented below have not been as extensively optimised for performance as other parts of the Catalyst package. Hence, you should take care when designing workflows which requires using them a large number of times.

!!! note
    Below we will describe various features using ODE simulations as examples. However, all interfaces (unless where else is stated) work identically for jump simulations.

## [Retrieving values from lattice simulations](@id lattice_simulation_structure_interaction_simulation_species)
Let us consider a simulation of a [`LatticeReactionSystem`](@ref):
```@example lattice_struct_interaction_sims
using Catalyst, OrdinaryDiffEq
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
diffusion_rx = @transport_reaction D X1
lattice = CartesianGrid((2,3))
lrs = LatticeReactionSystem(two_state_model, [diffusion_rx], lattice)

u0 = [:X1 => [0.0 0.0 0.0; 2.0 2.0 2.0], :X2 => 0.0]
tspan = (0.0, 1.0)
ps = [:k1 => 2.0, :k2 => 1.0, :D => 0.1]
oprob = ODEProblem(lrs, u0, tspan, ps)
sol = solve(oprob)
nothing # hide
```
To retrieve the values of $X1$ across the simulation we use the `lat_getu` function. It takes three arguments:
- The solution objects from which we wish to retrieve values.
- The species which values we wish to retrieve.
- The [`LatticeReactionSystem`](@ref) which was simulated.

```@example lattice_struct_interaction_sims
lat_getu(sol, :X1, lrs)
```
Here, the output is a vector with $X1$'s value at each simulation time step. How the species's value is represented at each time step depends on the lattice which was originally used to create the [`LatticeReactionSystem`](@ref):
- For Cartesian lattices, an array of the same size as the Cartesian lattice is used. Each array element corresponds to the species's value in the corresponding compartment.
- For masked lattices, a [sparse array](https://docs.julialang.org/en/v1/stdlib/SparseArrays/) of the same size as the masked lattice is used. Each filled array element corresponds to the species's value in the corresponding compartment. Unfilled array elements correspond to positions without compartments.
- For unstructured (graph) lattices, vectors are used. The i'th element in the vectors corresponds to the species's value in the i'th compartment.

Unlike for non-spatial simulations, `lat_getu` does not take vector (e.g. `lat_getu(sol, [:X1, :X2], lrs)`) or symbolic expression (e.g. `lat_getu(sol, [X1 + X2], lrs)`) inputs. However, it is possible to use symbolic variables as input (e.g. `lat_getu(sol, two_state_model.X1, lrs)`).

### [Retrieving lattice simulations values at specific time points](@id lattice_simulation_structure_interaction_simulation_species_ts)
Just like for non-spatial solutions, it is possible to access the simulation's values at designated time points. This is possible even if the simulation did not stop at those specific time points (in which case an interpolated value is returned). To do this, the desired time points to sample are provided as a vector to `lat_getu` using the optional argument `t`. E.g. here we retrieve the simulation's (interpolated) values at time points `0.5` and `0.75`:

```@example lattice_struct_interaction_sims
lat_getu(sol, :X1, lrs; t = [0.5, 0.75])
```

## [Retrieving and updating species values in problems and integrators](@id lattice_simulation_structure_interaction_prob_int_species)
Let us consider a spatial `ODEProblem`
```@example lattice_struct_interaction_prob_ints
using Catalyst, OrdinaryDiffEq
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
diffusion_rx = @transport_reaction D X1
lattice = CartesianGrid((2,3))
lrs = LatticeReactionSystem(two_state_model, [diffusion_rx], lattice)

u0 = [:X1 => [0.0 0.0 0.0; 2.0 2.0 2.0], :X2 => 0.0]
tspan = (0.0, 1.0)
ps = [:k1 => 2.0, :k2 => 1.0, :D => 0.1]
oprob = ODEProblem(lrs, u0, tspan, ps)
nothing # hide
```
We can retrieve the species values stored in `oprob` using the `lat_getu` function. It uses [identical syntax as for simulations](@ref lattice_simulation_structure_interaction_simulation_species) (except that you cannot specify a time point). However, it returns a single set of species values (while for simulations it returns a vector across different time steps):
```@example lattice_struct_interaction_prob_ints
lat_getu(oprob, :X1, lrs)
```
Again, the format used corresponds to the lattice used to create the original [`LatticeReactionSystem`](@ref). Here, even if a species has homogeneous values, the full format is used.
```@example lattice_struct_interaction_prob_ints
lat_getu(oprob, :X2, lrs)
```

For both problems and integrators, species values can be updated using the `lat_setu!` function. It uses a similar syntax as `lat_getu`, but takes a fourth argument which is the new values to use for the designated species:
```@example lattice_struct_interaction_prob_ints
lat_setu!(oprob, :X1, lrs, [1.0 2.0 3.0; 4.0 5.0 6.0])
```
Here, the same format (which depends on the used lattice) is used for the species's new values, as which is used when initially designating their initial conditions. I.e. to make $X1$'s initial condition values uniform we can call
```@example lattice_struct_interaction_prob_ints
lat_setu!(oprob, :X1, lrs, 1.0)
```

!!! note
    It is currently not possible to change a species value at a single compartment only. To do so, you must first retrieve its values across all compartments using `lat_getu`, then modify this at the desired compartment, and then use the modified version as input to `lat_setu!`.

Species values in [integrators](@ref simulation_structure_interfacing_integrators) can be interfaced with using identical syntax as for problems.

## [Retrieving and updating parameter values in problems and integrators](@id lattice_simulation_structure_interaction_prob_int_parameters)
Retrieval and updating of parameter values for problems and integrators works similarly as for species, but with the following differences:
- The `lat_getp` and `lat_setp!` functions are used.
- It is currently not possible to interface with parameter values of `JumpProblem`s and their integrators.
- After parameter values are modified, the `rebuild_lat_internals!` function must be applied before the problem/integrator can be used for further analysis.
- Updating of [edge parameters](@ref spatial_lattice_modelling_intro_simulation_edge_parameters) is limited and uses a different interface. 

Let us consider the spatial `ODEProblem` we previously declared. We can check the value of $k1$ by using `lat_getp`
```@example lattice_struct_interaction_prob_ints
lat_getp(oprob, :k1, lrs)
```
Next, we can update it using `lat_setp!` (here we also confirm that it now has the updated values):
```@example lattice_struct_interaction_prob_ints
lat_setp!(oprob, :k1, lrs, [1.0 2.0 3.0; 4.0 5.0 6.0])
lat_getp(oprob, :k1, lrs)
```

If we now were to simulate `oprob`, the simulation would not take the updated value of $k1$ into account. For our changes to take effect we might first need to call `rebuild_lat_internals!` with `oprob` as an input
```@example lattice_struct_interaction_prob_ints
rebuild_lat_internals!(oprob)
```
There are two different circumstances when `rebuild_lat_internals!` must be called:
- When modifying the value of an [edge parameter](@ref spatial_lattice_modelling_intro_simulation_edge_parameters).
- When changing a parameter from having spatially uniform values to spatially non-uniform values, or the other way around.

Parameter values of integrators can be interfaced with just like for problems (this is primarily relevant when using [*callbacks*](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)). Again, after doing so, the `rebuild_lat_internals!` function might need to be applied to the integrator.


### [Retrieving and updatingedge  parameter values in problems and integrators](@id lattice_simulation_structure_interaction_prob_int_parameters_edge_ps)
The `lat_getp` and `lat_setp!` functions cannot currently be applied to [edge parameters](@ref spatial_lattice_modelling_intro_simulation_edge_parameters). Instead, to access the value of an edge parameter, use 
```@example lattice_struct_interaction_prob_ints
oprob.ps[:D]
```
To update an edge parameter's value, use
```@example lattice_struct_interaction_prob_ints
oprob.ps[:D] = [0.2]
nothing # hide
```
This interface is somewhat limited, and the following aspects should be noted:
- Edge parameter values can only be interfaced with if the edge parameter's value is spatially uniform.
- When accessing an (spatially uniform) edge parameter's value, its single value will be encapsulated in a vector.
- When setting an (spatially uniform) edge parameter's value, you must encapsulate the new value in a vector.

