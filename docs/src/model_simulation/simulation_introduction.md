# [Model Simulation Introduction](@id simulation_intro)
Catalyst's core functionality is the creation of *chemical reaction network* (CRN) models that can be simulated using ODE, SDE, and jump simulations. How such simulations are carried out has already been described in [Catalyst's introduction](@ref introduction_to_catalyst). This page provides a deeper introduction, giving some additional background and introducing various simulation-related options. 

Here we will focus on the basics, with other sections of the simulation documentation describing various specialised features, or giving advice on performance. Anyone who plans on using Catalyst's simulation functionality extensively is recommended to also read the documentation on [solution plotting](@ref simulation_plotting), and on how to [interact with simulation problems, integrators, and solutions](@ref simulation_structure_interfacing). Anyone with an application for which performance is critical should consider reading the corresponding page on performance advice for [ODEs](@ref ode_simulation_performance), [SDEs](@ref ref), or [jump simulations](@ref ref).

### [Background to CRN simulations](@id simulation_intro_theory)
This section provides some brief theory on CRN simulations. For details on how to carry out these simulations in actual code, please skip to the following sections.

CRNs are defined by a set of *species* (with the amounts of these determining the system's state during simulations) and a set of *reaction events* (rules for how the state of the system changes). In real systems, the species amounts are *discrete copy-numbers*, describing the exact numbers of each species type present in the system (in systems biology this can e.g. be the number of a specific molecule present in a cell). Given rates for these reaction events, *stochastic chemical kinetics* provides a formula for simulating the system that recreates its real reaction process. During stochastic chemical kinetics simulations, the system's state is defined by discrete copy-numbers (denoting the number of each species present in the system). Next, at the occurrence of individual *reaction events*, the system's state is updated according to the occurred reaction. The result is a stochastic process. The most well-known approach for simulating stochastic chemical kinetics is [Gillespie's algorithm](https://en.wikipedia.org/wiki/Gillespie_algorithm).

In practice, these jump simulations are computationally expensive. In many cases, copy-numbers are so large that they can be approximated as *continuous concentrations*, and the time-development of the system as a *deterministic process*. This creates an ordinary differential equation (ODE), and is the chemical reaction network form most people are most familiar with. The rule for how ODEs are generated from CRNs is called the [*reaction rate equation*](https://en.wikipedia.org/wiki/Rate_equation) (RRE).

Here, the RRE enables fast, approximate, and deterministic simulations of CRNs, while stochastic chemical kinetics enables exact, stochastic, simulations of the true process. An intermediary approach is to use the [*chemical Langevin equation*](https://pubs.aip.org/aip/jcp/article/113/1/297/184125/The-chemical-Langevin-equation) (CLE) to formulate a stochastic differential equation (SDE). This approximates the system's state as continuous concentrations, but *does not* assume that its time development is deterministic. Generally, the CLE is used when copy-numbers are large enough that the continuous approximation holds, but not so large that the system's behaviour is deterministic. Generally, the advantage of SDE simulations (compared to jump ones) is that they are faster. Also, since the system state is continuous, interpretation of e.g. stability and steady state results from the deterministic (also continuous) domain is easier for SDEs (however one *should be careful* when making such interpretations).

These three different approaches are summed up in the following table:
```@raw html
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;}
.tg td{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
  overflow:hidden;padding:10px 5px;word-break:normal;}
.tg th{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
  font-weight:normal;overflow:hidden;padding:10px 5px;word-break:normal;}
.tg .tg-0pky{border-color:inherit;text-align:left;vertical-align:top}
</style>
<table class="tg">
<thead>
  <tr>
    <th class="tg-0pky">Interpretation</th>
    <th class="tg-0pky">Reaction rate equation</th>
    <th class="tg-0pky">Chemical Langevin equation</th>
    <th class="tg-0pky">Stochastic chemical kinetics</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky">Simulation form</td>
    <td class="tg-0pky">ODE simulations</td>
    <td class="tg-0pky">SDE simulations</td>
    <td class="tg-0pky">Jump simulations</td>
  </tr>
  <tr>
    <td class="tg-0pky">Example simulation methods</td>
    <td class="tg-0pky"><a href="https://en.wikipedia.org/wiki/Euler_method">Euler</a>, <a href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods">Runge-Kutta</a></td>
    <td class="tg-0pky"><a href="https://en.wikipedia.org/wiki/Euler%E2%80%93Maruyama_method">Euler-Maruyama</a>, <a href="https://en.wikipedia.org/wiki/Milstein_method">Milstein</a></td>
    <td class="tg-0pky"><a href="https://en.wikipedia.org/wiki/Gillespie_algorithm">Gillespie</a>, <a href="https://pubmed.ncbi.nlm.nih.gov/16321569/">Sorting direct</a></td>
  </tr>
  <tr>
    <td class="tg-0pky">Species units</td>
    <td class="tg-0pky">Concentration</td>
    <td class="tg-0pky">Concentration</td>
    <td class="tg-0pky">Copy-numbers</td>
  </tr>
  <tr>
    <td class="tg-0pky">Deterministic/Stochastic</td>
    <td class="tg-0pky">Deterministic</td>
    <td class="tg-0pky">Stochastic</td>
    <td class="tg-0pky">Stochastic</td>
  </tr>
  <tr>
    <td class="tg-0pky">Applicability</td>
    <td class="tg-0pky">Large species amounts</td>
    <td class="tg-0pky"><span style="font-weight:400;font-style:normal">Non-small species amounts</span></td>
    <td class="tg-0pky">Any species amounts</td>
  </tr>
  <tr>
    <td class="tg-0pky">Speed</td>
    <td class="tg-0pky">Typically fast</td>
    <td class="tg-0pky">Typically intermediate</td>
    <td class="tg-0pky">Typically slow</td>
  </tr>
  <tr>
    <td class="tg-0pky">Simulation package</td>
    <td class="tg-0pky"><a href="https://github.com/SciML/OrdinaryDiffEq.jl">OrdinaryDiffEq.jl</a></td>
    <td class="tg-0pky"><a href="https://github.com/SciML/StochasticDiffEq.jl">StochasticDiffEq.jl</a></td>
    <td class="tg-0pky"><a href="https://github.com/SciML/JumpProcesses.jl">JumpProcesses.jl</a></td>
  </tr>
</tbody>
</table>
```

## [Performing (ODE) simulations](@id simulation_intro_ODEs)
The following section gives a (more throughout than [previous]) introduction of how to simulate Catalyst models. This is exemplified using ODE simulations (some ODE-specific options will also be discussed). Later on, we will describe things specific to [SDE](@ref simulation_intro_SDEs) and [jump](@ref simulation_intro_jumps) simulations. All ODE simulations are performed using the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package, which full documentation can be found [here](https://docs.sciml.ai/OrdinaryDiffEq/stable/). A dedicated section giving advice on how to optimise ODE simulation performance can be found [here](@ref ode_simulation_performance)

To perform any simulation, we must first define our model, as well as the simulation's initial conditions, time span, and parameter values. Here we will use a simple [two-state model](@ref basic_CRN_library_two_states):
```@example simulation_intro_ode
using Catalyst
two_state_model = @reaction_network begin
  (k1,k2), X1 <--> X2
end
u0 = [:X1 => 100.0, :X2 => 200.0]
tspan = (0.0, 5.0)
ps = [:k1 => 2.0, :k2 => 5.0]
nothing # hide
```
To simulate the model we first bundle these up into an `ODEProblem`:
```@example simulation_intro_ode
oprob = ODEProblem(two_state_model, u0, tspan, ps)
nothing # hide
```
Next, we can simulate the model (requires loading the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package). Simulations are performed using the `solve` function.
```@example simulation_intro_ode
using OrdinaryDiffEq
sol = solve(oprob)
nothing # hide
```
Finally, the result can be plotted using the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package's `plot` function:
```@example simulation_intro_ode
using Plots
plot(sol)
```
More information on how to interact with solution structures is provided [here](@ref simulation_structure_interfacing) and on how to plot them [here](@ref simulation_plotting).

Some additional considerations:
- If a model without parameters has been declared, only the first three arguments must be provided to `ODEProblem`.
- While the first value of `tspan` will almost always be `0.0`, other starting times (both negative and positive) are possible.
- A discussion of various ways to represent species and parameters when designating their values in the `u0` and `ps` vectors can be found [here](@ref ref). 


### [Designating solvers and solver options](@id simulation_intro_solver_options)
While good defaults are generally selected, OrdinaryDiffEq enables the user to customise simulations through a long range of options that can be provided to the `solve` function. This includes specifying a [solver algorithm](https://en.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations), which can be provided as a second argument to `solve` (if none is provided, a suitable choice is automatically made). E.g. here we specify that the `Rodas5P` method should be used:
```@example simulation_intro_ode
sol = solve(oprob, Rodas5P())
nothing # hide
```
A full list of available solvers is provided [here](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/), and a discussion on optimal solver choices [here](@ref ode_simulation_performance_solvers).

Additional options can be provided as keyword arguments. E.g. the `adaptive` arguments determine whether adaptive time-stepping is used (for algorithms that permit this). This defaults to `true`, but can be disabled using
```@example simulation_intro_ode
sol = solve(oprob; adaptive = false)
nothing # hide
```

Here follows a list of solver options which might be of interest to the user.
- `adaptive`: Toggles adaptive time stepping for valid methods. Default to `true`.
- `dt`: For non-adaptive simulations, sets the step size (also sets the initial step size for adaptive methods).
- `saveat`: Determines the time points at which the simulation is saved. E.g. for `saveat = 2.0` the simulation is saved every second time unit. If not given, the solution is saved after each time step.
- `save_idxs`: Provides a vector of species whose values should be saved during the simulation. E.g. for `save_idxs = [:X1]`, only the value of species $X1$ is saved. 
- `maxiters`: The maximum number of time steps of the simulation. If this number is reached, the simulation is terminated.
- `seed`: Sets a seed for stochastic simulations. Stochastic simulations with the same seed generate identical results.

A full list of solver options can be found [here](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/).

### [Alternative problem input forms](@id simulation_intro_ODEs_input_forms)
Throughout Catalyst's documentation, we typically provide initial condition and parameter values as vectors. However, these can also be provided as tuples:
```@example simulation_intro_ode
u0 = (:X1 => 100.0, :X2 => 200.0)
tspan = (0.0, 5.0)
ps = (:k1 => 2.0, :k2 => 5.0)
oprob = ODEProblem(two_state_model, u0, tspan, ps)
nothing # hide
```
or dictionaries:
```@example simulation_intro_ode
u0 = Dict([:X1 => 100.0, :X2 => 200.0])
tspan = (0.0, 5.0)
ps = Dict([:k1 => 2.0, :k2 => 5.0])
oprob = ODEProblem(two_state_model, u0, tspan, ps)
nothing # hide
```
The forms used for `u0` and `ps` does not need to be the same (but can e.g. be a vector and a tuple). 

!!! note
    It is possible to [designate specific types for parameters](@ref dsl_advanced_options_parameter_types). When this is done, the tuple form for providing parameter values should be preferred.

Throughout Catalyst's documentation, we typically provide the time span as a tuple. However, if the first time point is `0.0` (which is typically the case), this can be omitted. Here, we supply only the simulation endpoint to our `ODEProblem`:
```@example simulation_intro_ode
tend = 5.0
oprob = ODEProblem(two_state_model, u0, tend, ps)
nothing # hide
```

## [Performing SDE simulations](@id simulation_intro_SDEs)
Catalyst uses the [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) package to perform SDE simulations. This section provides a brief introduction, with [StochasticDiffEq's documentation](https://docs.sciml.ai/StochasticDiffEq/stable/) providing a more extensive description. A dedicated section giving advice on how to optimise SDE simulation performance can be found [here](@ref ref). By default, Catalyst generates SDEs from CRN models using the chemical Langevin equation.

SDE simulations are performed in a similar manner to ODE simulations. The only exception is that an `SDEProblem` is created (rather than an `ODEProblem`). Furthermore, the [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) package (rather than the OrdinaryDiffEq package) is required for performing simulations. Here we simulate the two-state model for the same parameter set as previously used:
```@example simulation_intro_sde
using Catalyst, StochasticDiffEq
two_state_model = @reaction_network begin
  (k1,k2), X1 <--> X2
end
u0 = [:X1 => 100.0, :X2 => 200.0]
tspan = (0.0, 1.0)
ps = [:k1 => 2.0, :k2 => 5.0]

sprob = SDEProblem(two_state_model, u0, tspan, ps)
sol = solve(sprob, STrapezoid())
sol = solve(sprob, STrapezoid(); seed = 123) # hide
plot(sol)
```
we can see that while this simulation (unlike the ODE ones) exhibits some fluctuations. 

!!! note
    Unlike for ODE and jump simulations, there are no good heuristics for automatically selecting suitable SDE solvers. Hence, for SDE simulations a solver must be provided. `STrapezoid` will work for a large number of cases. When this is not the case, however, please check the list of [available SDE solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/) for a suitable alternative (making sure to select one compatible with non-diagonal noise and the [Ito interpretation]https://en.wikipedia.org/wiki/It%C3%B4_calculus).

### [Common SDE simulation pitfalls](@id simulation_intro_SDEs_pitfalls)
Next, let us reduce species amounts (using [`remake`](@ref simulation_structure_interfacing_problems_remake)), thereby also increasing the relative amount of noise, we encounter a problem when the model is simulated:
```@example simulation_intro_sde
sprob = remake(sprob; u0 = [:X1 => 0.33, :X2 => 0.66])
sol = solve(sprob, STrapezoid())
sol = solve(sprob, STrapezoid(); seed = 1234567) # hide
plot(sol)
```
Here, we receive a warning that the simulation was aborted. In the plot, we also see that it is incomplete. In this case we also note that species concentrations are very low (and sometimes, due to the relatively high amount of noise, even negative). This, combined with the early termination, suggests that we are simulating our model for too low species concentration for the assumptions of the CLE to hold. Instead, [jump simulations](@ref simulation_intro_jumps) should be used.

Next, let us consider a simulation for another parameter set:
```@example simulation_intro_sde
sprob = remake(sprob; u0 = [:X1 => 100.0, :X2 => 200.0], p = [:k1 => 200.0, :k2 => 500.0])
sol = solve(sprob, STrapezoid())
sol = solve(sprob, STrapezoid(); seed = 12345) # hide
plot(sol)
```
Again, the simulation is aborted. This time, however, species concentrations are relatively large, so the CLE might still hold. What has happened this time is that the accuracy of the simulations has not reached its desired threshold. This can be deal with [by reducing simulation tolerances](@ref ode_simulation_performance_error):
```@example simulation_intro_sde
sol = solve(sprob, STrapezoid(), abstol = 1e-1, reltol = 1e-1)
sol = solve(sprob, STrapezoid(); seed = 12345, abstol = 1e-1, reltol = 1e-1) # hide
plot(sol)
```

### [Scaling the noise in the chemical Langevin equation](@id simulation_intro_SDEs_noise_saling)
When using the CLE to generate SDEs from a CRN, it can sometimes be desirable to scale the magnitude of the noise. This can be done by introducing a *noise scaling term*, with each noise term generated by the CLE being multiplied with this term. A noise scaling term can be set using the `@default_noise_scaling` option:
```@example simulation_intro_sde
two_state_model = @reaction_network begin
  @default_noise_scaling 0.1
  (k1,k2), X1 <--> X2
end
```
Here, we set the noise scaling term to `0.1`, reducing the noise with a factor $10$ (noise scaling terms $>1.0$ increase the noise, while terms $<1.0$ reduce the noise). If we re-simulate the model using the low-concentration settings used previously, we see that the noise has been reduced (in fact by so much that the model can now be simulated without issues):
```@example simulation_intro_sde
u0 = [:X1 => 100.0, :X2 => 200.0]
tspan = (0.0, 1.0)
ps = [:k1 => 200.0, :k2 => 500.0]
sprob = SDEProblem(two_state_model, u0, tspan, ps)
sol = solve(sprob, STrapezoid())
sol = solve(sprob, STrapezoid(); seed = 123) # hide
plot(sol)
```

The `@default_noise_scaling` option can take any expression. This can be used to e.g. designate a *noise scaling parameter*:
```@example simulation_intro_sde
two_state_model = @reaction_network begin
  @parameters η
  @default_noise_scaling η
  (k1,k2), X1 <--> X2
end
```
Now we can tune the noise through $η$'s value. E.g. here we remove the noise entirely by setting $η = 0.0$ (thereby recreating an ODE simulation's behaviour):
```@example simulation_intro_sde
u0 = [:X1 => 0.33, :X2 => 0.66, :η => 0.0]
tspan = (0.0, 1.0)
ps = [:k1 => 2.0, :k2 => 5.0]
sprob = SDEProblem(two_state_model, u0, tspan, ps)
sol = solve(sprob, STrapezoid())
sol = solve(sprob, STrapezoid(); seed = 123) # hide
plot(sol)
```

!!! note
    Above, Catalyst is unable to infer that $η$ is a parameter from the `@default_noise_scaling η` option only. Hence, `@parameters η` is used to explicitly declare $η$ to be a parameter (as discussed in more detail [here](@ref dsl_advanced_options_declaring_species_and_parameters)).

It is possible to designate specific noise scaling terms for individual reactions through the `noise_scaling` [reaction metadata](@ref dsl_advanced_options_reaction_metadata). Here, CLE noise terms associated with a specific reaction are multiplied by that reaction's noise scaling term. Here we use this to turn off the noise in the $X1 \to X2$ reaction:
```@example simulation_intro_sde
two_state_model = @reaction_network begin
  k1, X1 <--> X2, [noise_scaling = 0.0]
  k2, X2 --> X1
end
nothing # hide
```
If the `@default_noise_scaling` option is used, that term is only applied to reactions *without* `noise_scaling` metadata.

While the `@default_noise_scaling` option is unavailable for [programmatically created models](@ref programmatic_CRN_construction), the [`remake_reactionsystem`](@ref) function can be used to achieve a similar effect.

## [Performing jump simulations using stochastic chemical kinetics](@id simulation_intro_jumps)

Catalyst uses the [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl) package to perform jump simulations. This section provides a brief introduction, with [JumpProcesses's documentation](https://docs.sciml.ai/JumpProcesses/stable/) providing a more extensive description. A dedicated section giving advice on how to optimise jump simulation performance can be found [here](@ref ref).

Jump simulations are performed using so-called `JumpProblem`s. Unlike ODEs and SDEs (for which the corresponding problem types can be created directly), jump simulations require first creating an intermediary `DiscreteProblem`. In this example, we first declare our two-state model and its initial conditions, time span, and parameter values.
```@example simulation_intro_jumps
using Catalyst
two_state_model = @reaction_network begin
  (k1,k2), X1 <--> X2
end
u0 = [:X1 => 5, :X2 => 10]
tspan = (0.0, 5.0)
ps = [:k1 => 2.0, :k2 => 5.0]
nothing # hide
```
!!! note
    Since jump simulations typically simulate the integer copy-numbers of each species present in the system, we designate our initial conditions for jump simulations as integers. Decimal-numbered initial conditions (and thus jump simulations) are, however, also possible. While ODE and SDE simulations accept integer initial conditions, these will be converted to decimal numbers.

Next, we bundle these into a `DiscreteProblem` (similarly to how `ODEProblem`s and `SDEProblem`s are created):
```@example simulation_intro_jumps
dprob = DiscreteProblem(two_state_model, u0, tspan, ps)
nothing # hide
```
This is then used as input to a `JumpProblem`. The `JumpProblem` also requires the CRN model as input.
```@example simulation_intro_jumps
jprob = JumpProblem(two_state_model, dprob, Direct())
nothing # hide
```
The `JumpProblem` can now be simulated using `solve` (just like any other problem type).
```@example simulation_intro_jumps
using JumpProcesses
sol = solve(jprob, SSAStepper())
nothing # hide
```
If we plot the solution we can see how the system's state does not change continuously, but instead in discrete jumps (due to the occurrence of the individual reactions of the system).
```@example simulation_intro_jumps
plot(sol)
```

### [Designating aggregators and simulation methods for jump simulations](@id simulation_intro_jumps_solver_designation)
Jump simulations (just like ODEs and SDEs) are performed using solver methods. Unlike ODEs and SDEs, jump simulations are carried out by two different types of methods acting in tandem. First, an *aggregator* method is used to (after each reaction) determine the time to, and type of, the next reaction. Next, a simulation method is used to actually carry out the simulation.

Several different aggregators are available (a full list is provided [here](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-for-Exact-Simulation)). To designate a specific one, provide it as the third argument to the `JumpProblem`. E.g. to designate that Gillespie's direct method (`Direct`) should be used, use:
```@example simulation_intro_jumps
jprob = JumpProblem(brusselator, dprob, Direct())
nothing # hide
```
Especially for large systems, the choice of aggregator is relevant to simulation performance. A guide for aggregator selection is provided [here](@ref ref).

Next, a simulation method can be provided (like for ODEs and SDEs) as the second argument to `solve`. Primarily two alternatives are available, `SSAStepper` and `FunctionMap` (other alternatives are only relevant when jump simulations are combined with ODEs/SDEs, which is described in more detail in JumpProcesses's documentation). Generally, `FunctionMap` is only used when a [continuous callback](@ref ref) is used (and `SSAStepper` otherwise). E.g. we can designate that the `FunctionMap` method should be used through:
```@example simulation_intro_jumps
sol = solve(jprob, FunctionMap())
nothing # hide
```

### [Jump simulations where some rate depends on time](@id simulation_intro_jumps_variableratejumps)
For some models, the rate of some reactions depend on time. E.g. consider the following [circadian model](https://en.wikipedia.org/wiki/Circadian_rhythm), where the production rate of some protein ($P$) depends on a sinusoid function:
```@example simulation_intro_jumps
circadian_model = @reaction_network begin
    A*(sin(2π*f*t - ϕ)+1)/2, 0 --> P
    d, P --> 0
end
```
This type of model will generate so called [*variable rate jumps*](@ref ref). Simulation of such model is non-trivial (and Catalyst currently lacks a good interface for this). A detailed description of how to carry out jump simulations for models with time-dependant rates can be found [here](https://docs.sciml.ai/JumpProcesses/stable/tutorials/simple_poisson_process/#VariableRateJumps-for-processes-that-are-not-constant-between-jumps).