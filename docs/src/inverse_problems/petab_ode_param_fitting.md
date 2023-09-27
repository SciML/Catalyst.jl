# [Parameter Fitting for ODEs using PEtab.jl](@id petab_parameter_fitting)
<<<<<<< refs/remotes/origin/master
<<<<<<< refs/remotes/origin/master
The [PEtab.jl package](https://github.com/sebapersson/PEtab.jl) implements the [PEtab format](https://petab.readthedocs.io/en/latest/) for fitting the parameters of deterministic CRN models to data ([please cite the corresponding papers if you use the PEtab approach in your research](@ref petab_citations)). PEtab.jl both implements methods for creating cost functions (determining parameter sets' fits to data), and for minimizing these cost functions. The PEtab approach covers most cases of fitting deterministic (ODE) models to data. When possible, it is recommended to use PEtab.jl for parameter fitting (and when not, an alternative, more general, approach will have to be used). This page describes how to combine PEtab.jl and Catalyst for parameter fitting, with a more extensive documentation available at https://sebapersson.github.io/PEtab.jl/stable/ (with this tutorial partially being an adaptation of this documentation).

## Introductory example
Let us consider a simple catalysis network, where an enzyme (*E*) turns a substrate (*S*) into a product (*P*):
```petab1
using Catalyst, PEtab
=======
The [PEtab.jl package](https://github.com/sebapersson/PEtab.jl) implements the [PEtab format](https://petab.readthedocs.io/en/latest/) for parameter fitting within systems biology. It implements both methods for creating a cost function (determening a parameter set's ability to fit a model to data), and for minimizing said cost function. This format covers most cases of fitting deterministic (ODE) models to data. For other cases, an alternative approach should be used. This page describes how to use PEtab for by Catalyst defined `ReactionSystem`, with more extesnive documentation being avaiable at https://sebapersson.github.io/PEtab.jl/stable/.
=======
The [PEtab.jl package](https://github.com/sebapersson/PEtab.jl) implements the [PEtab format](https://petab.readthedocs.io/en/latest/) for fitting the parameters of deterministic CRN models to data ([please cite the corresponding papers if you use the PEtab approach in you research](@ref petab_citations)). PEtab.jl both implements methods for creating cost functions (determining parameter sets' fits to data), and for minimizing these cost functions. The PEtab approach covers most cases of fitting deterministic (ODE) models to data. When possible, it is recommended to use PEtab.jl for parameter fitting (and when not, an alternative, more general, approach will have to be used). This page describes how to combine PEtab.jl and Catalyst for parameter fitting, with a more extensive documentation available at https://sebapersson.github.io/PEtab.jl/stable/ (with this tutorial partially being an adaptation of this documentation).
>>>>>>> init tutorial

## Introductory example
Let us consider a simple catalysis network, where an enzyme (*E*) turns a substrate (*S*) into a product (*P*):
```petab1
<<<<<<< refs/remotes/origin/master
using Catalyst
>>>>>>> initiate tutorial
=======
using Catalyst, PEtab
>>>>>>> init tutorial

rn = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
```
<<<<<<< refs/remotes/origin/master
<<<<<<< refs/remotes/origin/master
From some known initial condition, and a true parameter set (which we later want to recover from the data) we generate synthetic data (on which we will demonstrate the fitting process).
```petab1
# Define initial conditions and parameters.
u0 = [:S => 1.0, :E => 1.0, :SE => 0.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5]

# Generate synthetic data.
using OrdinaryDiffEq, Random
oprob = ODEProblem(rn, u0, (0.0, 10.0), p_true)
sol = solve(oprob, Tsit5(); saveat=0.1)
data = (0.8 .+ 0.4*rand(10)) .* sol[:P][10:10:end]

# Plots the true solutions and the (synthetic data) measurements.
using Plots
plot(sol; idxs=:P, label="True solution")
plot!(sol.t[10:10:end], data; label="Measurements", seriestype=:scatter)
```

Generally, PEtab takes five different inputs to define an optimisation problem (the minimiser of which generates a fitted parameter set):
1. **Model**: The model which we want to fit to the data (a `ReactionSystem`).
2. **Observables**: The possible observables that can be measured (a `Dict{String,PEtabObservable}`).
3. **Estimation parameters**: The parameters which we want to fit to the data (a `Vector{PEtabParameter}`). 
4. **Experimental (or simulation) conditions**: The simulations (each corresponding to a potential experiment) carried out during each step of the optimisation process (a `Dict{String,Dict{Symbol,Float64}}`).
5. **Measurements**: The measurements to which the model is fitted (a `DataFrame`). 

### Observables
The observables define the quantities that we may measure in our experiments. Typically, each corresponds to a single species, however, [more complicated observables are possible](@ref petab_observables_observables). For each observable, we also need a noise formula, defining the uncertainty in the measurement. By default, PEtab assumes normally distributed noise, with a mean equal to the true value and a standard deviation which we have to define. It is also possible to use [more advanced noise formulas](@ref petab_observables_noise_formula).

In our example, we only have a single possible observable, the `P` species. We will assume that the noise is normally distributed with a standard deviation `0.5` (in our case this is not true, however, typically the noise distribution is unknown). We combine this information in a `PEtabObservable` struct (to access the `P` species we must use [`@unpack`](@ref simulation_structure_interfacing_symbolic_representation)). Finally, we store all our observables in a dictionary, giving each an id tag (which is later used in the measurements input).
```petab1
@unpack P = rn
obs_P = PEtabObservable(P, 0.5)
observables = Dict("obs_P" => obs_P)
nothing # hide
```

### Parameters
Each parameter of the system can either be
1. Known (described [here](@ref petab_parameters_known)).
2. Depend on experimental/simulated conditions (described [here](@ref petab_simulation_conditions)). 
3. Be an unknown we wish to fit to data.

In our case, we wish to fit all three system parameters (*kB*, *kD*, and *kP*). For each, we create a single `PEtabParameter`, and then gather these into a single vector.
```petab1
par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
params = [par_kB, par_kD, par_kP]
nothing # hide
```
For each parameter, it is also possible to set [a lower or upper bound](@ref petab_parameters_bounds), set whether to use [logarithmic or linear scale](@ref petab_parameters_scales), or add a [prior distribution of its value](@ref petab_parameters_priors).

### Simulation conditions
Sometimes, several different experiments are performed on a system (each potentially generating several measurements). An experiment could be the time development of a system from a specific initial condition. Since each experimental condition (during the optimisation procedure, for a guess of the unknown parameters) generates a distinct simulation, these are also called simulation conditions. In our example, all data comes from a single experiment, and the simulation condition input is not required. How to define and use different experimental conditions is described [here](@ref petab_simulation_conditions).

### Measurements
Finally, we need to define the system measurements to which the parameters will be fitted. Each measurement combines:
1. The observable which is observed (here we use the id tag defined in the `observables` dictionary).
2. The time point of the measurement.
3. The measured value.
  
For cases where several simulation conditions are given, we also need to provide:

4. The simulation condition which generates the measurement ([here](@ref petab_simulation_conditions) is an example where this is used).

Measurements are provided in a `DataFrame` where each row corresponds to a measurement. The values are provided in the `obs_id`, `time`, `measurement`, and `simulation_id` columns. In our case we only need to fill in the first three:
```petab1
using DataFrames
measurements = DataFrame(obs_id="obs_P", time=sol.t[10:10:end], measurement=data)
```
Since, in our example, all measurements are of the same observable, we can set `obs_id="obs_P"`. If several different observables were measured, `obs_id` would be a vector of the same length as `time` and `measurement`.

### Creating a PEtabModel
Finally, we combine all inputs in a single `PEtabModel`. To it, we also pass the initial conditions of our simulations (using the `state_map` argument). It is also possible to have [initial conditions with uncertainty](@ref petab_simulation_initial_conditions_uncertainty), [that vary between different simulations](@ref petab_simulation_conditions), or [that we attempt to fit to the data](@ref petab_simulation_initial_conditions_fitted).
```petab1
petab_model = PEtabModel(rn, observables, measurements, params; state_map=u0)
nothing # hide
```

### Fitting parameters
We are now able to fit our model to the data. First, we create a `PEtabODEProblem`. Here, we use `petab_problem` as the only input, but it is also possible to set various [numeric solver and automatic differentiation option](@ref petab_simulation_conditions) (such as method or tolerance).
```petab1
petab_problem = PEtabODEProblem(petab_model)
nothing # hide
```
To fit a parameter set we use the `calibrate_model` function. In addition to our `PEtabODEProblem`, we must also provide an initial guess (which can be generated with the `generate_startguesses` function) and an optimisation algorithm (which needs to be imported specifically). PEtab.jl supports various [optimisation methods and options](@ref petab_optimisation_optimisers). 
```petab1
using Optim
p0 = generate_startguesses(petab_problem, 1)
res = calibrate_model(petab_problem, p0, IPNewton())
```

We can now simulate our model for the fitted parameter set, and compare the result to the measurements and true solution.
```petab1
oprob = ODEProblem(rn, u0, (0.0, 10.0), get_ps(res, petab_model))
sol = solve(oprob, Tsit5(); saveat=0.1)
plot!(sol; idxs=4, label="Fitted solution")
```
Here we use the `get_ps` function to retrieve a full parameter set using the optimal parameters. The calibration result can also be found in `res.xmin`, however, note that PEtab automatically ([unless a linear scale is selected](petab_parameters_scales)) converts parameters to logarithmic scale, so typically `10 .^res.xmin` are the values of interest. If you investigate the result from this example you might note that the found parameter set (while it fits the data well) does not correspond to the true parameter set. This phenomenon is related to *identifiability*, and is very important for parameter fitting.

### Final notes
PEtab.jl also supports [multistart optimisation](@ref petab_multistart_optimisation), [automatic equilibration before simulations](https://sebapersson.github.io/PEtab.jl/stable/Brannmark/), and [events](@ref petab_events). Various [plot recipes](@ref petab_plotting) exist for investigating the optimisation process. Please read the [PETab.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/) for a more complete description of the package's features.

## [Additional features: Observables](@id petab_observables)

### [Defining non-trivial observables](@id petab_observables_observables)
It is possible for observables to be any algebraic expression of species concentrations and parameters. E.g. in this example the total amount of `X` in the system is an observable:
```petab2
using Catalyst, PEtab # hide
rn = @reaction_network begin
    (k1,k2), X1 <--> X2
end
@unpack X1, X2 = rn
obs_X = PEtabObservable(X1 + X2, 0.5)
```

### [Advanced observables noise formulas](@id petab_observables_noise_formula)
In our basic example we assumed that the normally distributed noise had a standard deviation of `0.5`. However, this value may be a parameter (or indeed any algebraic expression). E.g, we could set
```petab1
@parameters σ
obs_P = PEtabObservable(P, σ)
```
and then let the parameter *σ* vary between different [simulation conditions](@ref petab_simulation_conditions).

PEtab.jl assumes that noise is normally distributed (with a standard deviation equal to the second argument provided to `PEtabObservable`). The only other (currently implemented) noise distribution is log-normally distributed noise, which is designated through the `transformation=:log` argument:
```petab1
obs_P = PEtabObservable(P, σ; transformation=:log)
```

## [Additional features: Parameters](@id petab_parameters)

### [Known parameters](@id petab_parameters_known)
In our previous example, all parameters were unknowns that we wished to fit to the data. If any parameters have known values, it is possible to provide these to `PEtabModel` through the `parameter_map` argument. E.g if we had known that *kB = 1.0*, then we would only define *kD* and *kP* as parameters we wish to fit:
```petab1
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
params = [par_kD, par_kP]
nothing # hide
```
We then provide `parameter_map=[:kB => 1.0]` when we assembly our model:
```petab1
petab_model = PEtabModel(rn, observables, measurements, params; state_map=u0, parameter_map=[:kB => 1.0])
nothing # hide
```

### [Parameter bounds](@id petab_parameters_bounds)
By default, when fitted, potential parameter values are assumed to be in the interval *[1e-3, 1e3]*. When declaring a `PEtabParameter` it is possible to change these values through the `lb` and `ub` arguments. E.g. we could use
```petab1
par_kB = PEtabParameter(:kB; lb=1e-2, ub=1e-2)
```
to achieve the more conservative bound *[1e-2, 1e2]* for the parameter *kB*.

### [Parameter scales](@id petab_parameters_scales)

Internally, parameters that are fitted are converted to a logarithmic scale (this has several advantages, e.g. automatically preventing negative parameter values). To prevent this conversion, the `scale=:lin` argument can be used. Here we change the scale of *kB* to linear:
```petab1
par_kB = PEtabParameter(:kB; scale=:lin)
```

### [Parameter priors](@id petab_parameters_priors)
If one has prior knowledge about the distribution of a parameter, it is possible to incorporate this into the model. The prior can be any continuous distribution from the [Distributions.jl package](https://github.com/JuliaStats/Distributions.jl). E.g we can use:
```petab1
using Distributions
par_kB = PEtabParameter(:kB; prior=Normal(1.0,0.2))
```
to set a normally distributed prior (with mean `1.0` and standard deviation `0.2`) on the value of *kB*. By default, the prior is assumed to be on the linear scale of the parameter (before any potential log transform). To specify that the prior is on the logarithmic scale, the `prior_on_linear_scale=false` argument can be provided:
```petab1
par_kB = PEtabParameter(:kB; prior=Normal(1.0,0.2), prior_on_linear_scale=false)
```

## [Simulation conditions](@id petab_simulation_conditions)
Sometimes, we have data from several different experimental conditions. Here, when a potential parameter set is evaluated during the fitting process, each experimental condition corresponds to one simulation condition (which produces one simulation). To account for this, PEtab permits us to define several different simulation conditions, with each condition being defined by specific values for some initial conditions and/or parameters. 

If, for our previous catalysis example, we had measured the system for two different initial values of *S* (*S(0)=1.0* and *S(0)=0.5*), these would correspond to two different simulation conditions. For each condition we define a `Dict` mapping the species to their initial condition (here, *S* is the only species in each `Dict`):
```petab1
c1 = Dict(:S => 1.0)
c2 = Dict(:S => 0.5)
nothing # hide
```
Similarly as for observables, we then gather the conditions in another `Dict`, giving each an id tag:
```petab1
simulation_conditions = Dict("c1" => c1, "c2" => c2)
nothing # hide
```
Again (like for observables), each measurement in the measurements `DataFrame` needs to be associated with a simulation condition id tag (describing which condition that measurements were taken from). Parameters, just like initial conditions, may vary between different conditions. If an initial condition (or parameter) occurs in one condition, it must occur in all of them.

Here follows a complete version of our basic example, but with measurements both for *S(0)=1.0* and *S(0)=0.5*.
```petab3
using Catalyst, PEtab

rn = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end

u0 = [:E => 1.0, :SE => 0.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5]

# Simulate data.
using OrdinaryDiffEq, Random
d1, t1 = let 
    oprob = ODEProblem(rn, [:S=>1.0; u0], (0.0, 10.0), p_true)
    sol = solve(oprob, Tsit5(); saveat=0.1)
    ((0.8 .+ 0.4*rand(10)) .* sol[:P][10:10:end], sol.t[10:10:end])
end
d2, t2 = let 
    oprob = ODEProblem(rn, [:S=>0.5; u0], (0.0, 10.0), p_true)
    sol = solve(oprob, Tsit5(); saveat=0.1)
    ((0.8 .+ 0.4*rand(10)) .* sol[:P][10:10:end], sol.t[10:10:end])
end

@unpack P = rn
obs_P = PEtabObservable(P, 0.5)
observables = Dict("obs_P" => obs_P)

par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
params = [par_kB, par_kD, par_kP]

c1 = Dict(:kB => 1.0)
c2 = Dict(:kB => 0.5)
simulation_conditions = Dict("c1" => c1, "c2" => c2)

using DataFrames
m1 = DataFrame(simulation_id="c1", obs_id="obs_P", time=t1, measurement=d1)
m2 = DataFrame(simulation_id="c2", obs_id="obs_P", time=t2, measurement=d2)
measurements = vcat(m1,m2)

petab_model = PEtabModel(rn, simulation_conditions, observables, measurements, params; state_map=u0)
nothing # hide
```
Note that the `u0` we pass into `PEtabModel` through the `state_map` argument no longer contains the value of *S* (as it is provided by the conditions). 

### [Varying parameters between different simulation conditions](@id petab_simulation_conditions)
Sometimes, the parameters that are used vary between the different conditions. Consider our catalysis example, if we had performed the experiment twice, using two different enzymes with different catalytic properties, this could have generated such conditions. Here, let us assume that the two enzymes yield different rates (*kP₁* and *kP₂*) for the `SE --> P + E` reaction (but are otherwise identical). 

[The author (Torkel) currently unsure how to carry this out, but will complete this tutorial when he knows.]

## [Additional features: Initial conditions](@id petab_simulation_initial_conditions)

### [Fitting initial conditions](@id petab_simulation_initial_conditions_fitted)
Sometimes, initial conditions are uncertain quantities which we wish to fit to the data. This is possible by [defining an initial condition as a parameter](dsl_description_parametric_initial_conditions): 
```petab4
using Catalyst, PEtab # hide
rn = @reaction_network begin
    @parameters E0
    @species E(t)=E0
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
nothing # hide
```
Here, the initial value of `E` is equal to the parameter `E0`. We modify our `u0` vector by removing `E` (which is no longer known):
```petab4
u0 = [:S => 1.0, :SE => 0.0, :P => 0.0]
nothing # hide
```
Next, we add `E0` to the parameters we wish to fit:
```petab4
par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
par_E0 = PEtabParameter(:E0)
params = [par_kB, par_kD, par_kP, par_E0]
nothing # hide
```
and we can use our updated `rn`, `u0`, and `params` as input to our `PEtabModel`.

### [Uncertain initial conditions](@id petab_simulation_initial_conditions_uncertainty)
Often, while an initial condition has been reported for an experiment, its exact value is uncertain. This can be modelled by making the initial condition a [parameter that is fitted to the data](@id petab_simulation_initial_conditions_fitted) and attach a prior to it corresponding to our certainty about its value. 

Let us consider our initial example, but where we want to add uncertainty to the initial conditions of `S` and `E`. We will add priors on these, assuming normal distributions with mean `1.0` and standard deviation `0.1`. For the synthetic measured data we will assume true values *S(0)=E(0)=1.0*.
```petab4
using Catalyst, PEtab

rn = @reaction_network begin
    @parameters S0 E0
    @species S(t)=S0 E(t)=E0
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end

u0 = [:SE => 0.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5]

using OrdinaryDiffEq, Random
oprob = ODEProblem(rn, [:S=>1.0; :E => 1.0; u0], (0.0, 10.0), p_true)
sol = solve(oprob, Tsit5(); saveat=0.1)
data = ((0.8 .+ 0.4*rand(10)) .* sol[:P][10:10:end], sol.t[10:10:end])

@unpack P = rn
obs_P = PEtabObservable(P, 0.5)
observables = Dict("obs_P" => obs_P)

par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
par_S0 = PEtabParameter(:S0)
par_E0 = PEtabParameter(:E0)
params = [par_kB, par_kD, par_kP, par_S0, par_E0]

using DataFrames
measurements = DataFrame(obs_id="obs_P", time=sol.t[10:10:end], measurement=data)

petab_model = PEtabModel(rn, observables, measurements, params; state_map=u0)
nothing # hide
```
Here, when we fit our data we will also gain values for `S0` and `E0`, however, unless we are explicitly interested in these, they can be ignored.

## [Additional features: Simulation options](@id petab_simulation_options)

While in our basic example, we do not provide any additional information to our `PEtabODEProblem`, this is an opportunity to specify how the model should be simulated, and what differentiation techniques to use for the optimisation procedure (if none are provided, appropriate defaults are selected). 

Here is an example, taken from the [more detailed PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/dev/Boehm/#Creating-a-PEtabODEProblem)
```petab1
petab_problem = PEtabODEProblem(petab_model, 
                                ode_solver=ODESolver(Rodas5P(), abstol=1e-8, reltol=1e-8), 
                                gradient_method=:ForwardDiff, 
                                hessian_method=:ForwardDiff)
nothing # hide
```
where we simulate our ODE model using the `Rodas5p` method (with absolute and relative tolerance both equal `1e-8`) and use [forward automatic differentiation](https://github.com/JuliaDiff/ForwardDiff.jl) for both gradient and hessian computation. More details on available ODE solver options can be found in [the PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/dev/API_choosen/#PEtab.ODESolver).

## [Additional features: Optimisation](@id petab_optimisation)

### [Optimisation methods and options](@id petab_optimisation_optimisers)
For our examples, we have used the `Optim.IPNewton` optimisation method. PEtab.jl supports [several additional optimisation methods](https://sebapersson.github.io/PEtab.jl/dev/Avaible_optimisers/#options_optimisers). Furthermore, the `calibrate_model` `options` argument permits the customisation of the options for any used optimiser. E.g. to increase the maximum number of iterations of the `Optim.IPNewton` method we would use:
```petab1
res = calibrate_model(petab_problem, p0, IPNewton(); options=Optim.Options(iterations = 10000))
nothing # hide
```

Please read the [PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/dev/Avaible_optimisers/#options_optimisers) to learn how to customise the various optimisers' properties. 

### Optimisation path recording
To record all the parameter sets evaluated (and their objective values) during the optimisation procedure, provide the `save_trace=true` argument to either `calibrate_model` or `calibrate_model_multistart`. This is required for the various [optimisation evaluation plots](@ref petab_plotting) provided by PEtab.jl.


## Objective function extraction
While PEtab.jl provides various tools for analysing the objective function generated by `PEtabODEProblem`, it is also possible to extract this function for customised analysis. Given a `PEtabODEProblem`
```petab1
petab_problem = PEtabODEProblem(petab_model)
nothing # hide
```
We can find the:
1. Objective function in the `petab_problem.compute_cost`. It takes a single argument (`p`) and returns the objective value.
2. Gradient in the `petab_problem.compute_gradient!` field. It takes two arguments (`g` and `p`) with `g` containing the updated gradient values.
3. Hessian in the `petab_problem.compute_hessian` field. It takes two arguments (`H` and `p`) with `H` containing the updated hessian values.

## [Multistart optimisation](@id petab_multistart_optimisation)
To avoid the optimisation process returning a local minimum, it is often advised to run it multiple, using different initial guesses. PEtab.jl supports this through the `calibrate_model_multistart` function. This is identical to the `calibrate_model` function, but takes two additional arguments:
1. `n_multistart`: The number of runs to perform.
2. `dir_save`: A location to which the output is saved. If `dir_save=nothing`, no output is saved.

And one additional optional argument:

3. `sampling_method`: Selects the sampling method with which to select the initial guesses (`QuasiMonteCarlo.LatinHypercubeSample()` used by default). 
 
Because `calibrate_model_multistart` handles initial guess sampling, unlike for `calibrate_model`, no initial guess has to be provided. 

Here, we fit parameters through 10 independent optimisation runs, using QuasiMonteCarlo's `SobolSample` method, and save the result to the "OptimizationRuns/runs_ms_1" folder:
```petab1
using Optim
using QuasiMonteCarlo
res_ms = calibrate_model_multistart(petab_problem, IPNewton(), 10, "OptimizationRuns/runs_ms_1"; sampling_method=QuasiMonteCarlo.SobolSample())
```
The best result across all runs can still be retrieved using `get_ps(res_ms, petab_model)`, with the results of the individual runs being stored in the `res_ms.runs` field.

## [Events](@id petab_events)
So far, we have assumed that all experiments, after initiation, run without interference. Experiments where conditions change, or where species are added/removed during the time course, can be represented through events (related to [callbacks](@ref advanced_simulations_callbacks)). In PEtab, an event is represented through the `PEtabEvent` structure. It takes three arguments:
1. The condition for triggering the event. This can either indicate a point in time, or a boolean condition.
2. A rule for updating the event's target
3. The event's target (either a species or parameter).

Here we create an event which adds `0.5` units of `S` to the system at time `5.0`:
```petab1
@unpack S = rn
event = PEtabEvent(5.0, S + 0.5, S)
```
Here, the first argument is evaluated to a scalar value, and which case it is interpreted as a time point at which the event happens. If we instead want the event to happen whenever the concentration of `S` falls below `0.5` we set the first argument to a boolean condition:
```petab1
event = PEtabEvent(S < 0.5, S + 0.5, S)
```

Whenever we have several events or not, we bundle them together in a single vector which is later passed to the `PEtabODEProblem` using the `events` argument:
```petab1
events = [event]
petab_model = PEtabModel(rn, observables, measurements, params; state_map=u0, events=events)
```

For more details on how to use events, including how t

## [Plot recipes](@id petab_plotting)
There exist various types of graphs that can be used to evaluate the parameter fitting process. These can be plotted using the `plot` command, where the input is either the result of a `calibrate_model` or a  `calibrate_model_multistart` run. To, for a single start calibration run, plot, for each iteration of the optimization process, the best objective value achieved so far, run:
```petab1
plot(res)
```
If the input instead is a multi-start calibration run, the output is a so-called waterfall plot:
```petab1
plot(res_ms)
```
Here, each dot shows the final objective value for a single run in the multi-start process. The runs are ordered by their objective values, and colours designate runs in the same local minimum.

To instead use the best objective value plot for a multi-start run (with one curve for each run), the `plot_type` argument is used:
```petab1
plot(res_ms; plot_type = :best_objective)
```
Generally, there exist several types of plots for both types of calibration results. More details of the types of available plots, and how to customise them, can be found [here](@ref https://sebapersson.github.io/PEtab.jl/stable/optimisation_output_plotting/).


## [Citations](@id petab_citations)
If you use this functionality in your research, [in addition to Catalyst](@ref catalyst_citation), please cite the following papers to support the authors of the PEtab.jl package (currently there is no article associated with this package) and the PEtab standard:
```
@misc{2023Petabljl,
  author       = {Ognissanti, Damiano and Arutjunjan, Rafael and Persson, Sebastian and Hasselgren, Viktor},
=======
From some known initial condition, and a parameter set (which we wish to fit) we generate synthetic data (on which we will demonstrate the fitting process).
=======
From some known initial condition, and a true parameter set (which we later want to recover from the data) we generate synthetic data (on which we will demonstrate the fitting process).
>>>>>>> init tutorial
```petab1
# Define initial conditions and parameters.
u0 = [:S => 1.0, :E => 1.0, :SE => 0.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5]

# Generate synthetic data.
using OrdinaryDiffEq, Random
oprob = ODEProblem(rn, u0, (0.0, 10.0), p_true)
sol = solve(oprob, Tsit5(); saveat=0.1)
data = (0.8 .+ 0.4*rand(10)) .* sol[:P][10:10:end]

# Plots the true solutions and the (synthetic data) measurements.
using Plots
plot(sol; idxs=:P, label="True solution")
plot!(sol.t[10:10:end], data; label="Measurements", seriestype=:scatter)
```

Generally, PEtab takes five different inputs to define an optimisation problem (the minimiser of which generates a fitted parameter set):
1. **Model**: The model which we want to fit to the data (a `ReactionSystem`).
2. **Observables**: The possible observables that can be measured (a `Dict{String,PEtabObservable}`).
3. **Estimation parameters**: The parameters which we want to fit to the data (a `Vector{PEtabParameter}`). 
4. **Experimental (or simulation) conditions**: The simulations (each corresponding to a potential experiment) carried out during each step of the optimisation process (a `Dict{String,Dict{Symbol,Float64}}`).
5. **Measurements**: The measurements to which the model is fitted (a `DataFrame`). 

### Observables
The observables define the quantities that we may measure in our experiments. Typically, each corresponds to a single species, however, [more complicated observables are possible](@ref petab_observables_observables). For each observable, we also need a noise formula, defining the uncertainty in the measurement. By default, PEtab assumes normally distributed noise, with a mean equal to the true value and a standard deviation which we have to define. It is also possible to use [more advanced noise formulas](@ref petab_observables_noise_formula).

In our example, we only have a single possible observable, the `P` species. We will assume that the noise is normal distributed with standard deviation `0.5` (in our case this is not true, however, typically the noise distribution is unknown). We combine this information in a `PEtabObservable` struct (to access the `P` species we must use [`@unpack`](@ref simulation_structure_interfacing_symbolic_representation)). Finally, we store all our observables in a dictionary, giving each an id tag (which is later used in the measurements input).
```petab1
@unpack P = rn
obs_P = PEtabObservable(P, 0.5)
observables = Dict("obs_P" => obs_P)
nothing # hide
```

### Parameters
Each parameter of the system can either be
1. Known (described [here](@ref petab_parameters_known)).
2. Depend on experimental/simulated conditions (described [here](@ref petab_simulation_conditions)). 
3. Be an unknown we wish to fit to data.

In our case we wish to fit all three system parameters (*kB*, *kD*, and *kP*). For each, we create a single `PEtabParameter`, and then gather these into a single vector.
```petab1
par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
parameters = [par_kB, par_kD, par_kP]
nothing # hide
```
For each parameter, it is also possible to set [a lower or upper bound](@ref petab_parameters_bounds), set whether to use [logarithmic or linear scale](@ref petab_parameters_scales), or add a [prior distribution of its value](@ref petab_parameters_priors).

### Simulation conditions
Sometimes, several different experiments are performed on a system (each potentially generating several measurements). An experiment could be the time development of a system from a specific initial condition. Since each experimental condition (during the optimisation procedure, for a guess of the unknown parameters) generate a distinct simulation, these are also called simulation conditions. In our example, all data comes from a single experiment, and the simulation condition input is not required. How to define and use different experimental conditions is described [here](@ref petab_simulation_conditions).

### Measurements
Finally, we need to define the system measurements to which the parameters will be fitted. Each measurement combine:
1. The observable which is observed (here we use the id tag defined in the `observables` dictionary).
2. The time point of the measurement.
3. The measured value.
  
For cases where several simulation conditions are given, we also need to provide:

4. The simulation condition which generates the measurement ([here](@ref petab_simulation_conditions) is an example where this is used).

Measurements are provided in a `DataFrame` where each row correspond to a measurement. The values are provided in the `obs_id`, `time`, `measurement`, and `simulation_id` columns. In our case we only need to fill in the first three:
```petab1
using DataFrames
measurements = DataFrame(obs_id="obs_P", time=sol.t[10:10:end], measurement=data)
```
Since, in our example, all measurements is of the same observable, we can set `obs_id="obs_P"`. If several different observables were measured, `obs_id` would be a vector of the same length as `time` and `measurement`.

### Creating a PEtabModel
Finally, we combine all inputs in a single `PEtabModel`. To it, we also pass the initial conditions of our simulations (using the `state_map` argument). It is also possible to have [initial conditions with uncertainty](@ref petab_simulation_initial_conditions_uncertainty), [that vary between different simulations](@ref petab_simulation_conditions), or [that we attempt to fit to the data](@ref petab_simulation_initial_conditions_fitted).
```petab1
petab_model = PEtabModel(rn, observables, measurements, parameters; state_map=u0)
nothing # hide
```

### Fitting parameters
We are now able to fit our model to the data. First we create a `PEtabODEProblem`. Here, we use `petab_problem` as the only input, but it is also possible to set various [numeric solver and automatic differentiation option](@ref petab_simulation_conditions) (such as method or tolerance).
```petab1
petab_problem = PEtabODEProblem(petab_model)
nothing # hide
```
To fit a parameters we use the `calibrate_model` function. In addition to our `PEtabODEProblem`, we must also provide an initial guess (which can be generated with the `generate_startguesses` function) and an optimisation algorithm (which needs to be imported specifically). PEtab.jl supports various [optimisation methods and options](@ref petab_optimisation_optimisers). 
```petab1
using Optim
p0 = generate_startguesses(petab_problem, 1)
res = calibrate_model(petab_problem, p0, IPNewton())
```

We can now simulate our model for the fitted parameter set, and compare the result to the measurements and true solution.
```petab1
oprob = ODEProblem(rn, u0, (0.0, 10.0), get_best_ps(res, petab_model))
sol = solve(oprob, Tsit5(); saveat=0.1)
plot!(sol; idxs=4, label="Fitted solution")
```
Here we use the `get_best_ps` function to retrieve a full parameter set using the optimal parameters. The calibration result can also be found in `res.xmin`, however, note that PEtab automatically ([unless linear scale is selected](petab_parameters_scales)) converts parameters to logarithmic scale, so typically `10 .^res.xmin` are the values of interest. If you investigate the result from this example you might note that the found parameter set (while it fits the data well) does not correspond to the true parameter set. This phenomena is related to *identifiability*, and is very important for parameter fitting.

### Final notes
PEtab.jl also supports [multistart optimisation](@ref petab_multistart_optimisation), [automatic equilibration before simulations](https://sebapersson.github.io/PEtab.jl/stable/Brannmark/), and [events](@ref petab_events). Various [plot recipes](@ref petab_plotting) exist for investigating the optimisation process. Please read the [PETab.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/) for a more complete description of the package's features.

## [Additional features: Observables](@id petab_observables)

### [Defining non-trivial observables](@id petab_observables_observables)
It is possible for observables to be any algebraic expression of species concentrations and parameters. E.g. in this example the total amount of `X` in the system is an observable:
```petab2
using Catalyst, PEtab # hide
rn = @reaction_network begin
    (k1,k2), X1 <--> X2
end
@unpack X1, X2 = rn
obs_X = PEtabObservable(X1 + X2, 0.5)
```

### [Advanced observables noise formulas](@id petab_observables_noise_formula)
In out basic example we assumed that the normally distributed noise had a standard deviation of `0.5`. However, this value may be a parameter (or indeed any algebraic expression). E.g, we could set
```petab1
@parameters σ
obs_P = PEtabObservable(P, σ)
```
and then let the parameter *σ* vary between different [simulation conditions](@ref petab_simulation_conditions).

PEtab.jl assumes that noise is normally distributed (with standard deviation equal to the second argument provided to `PEtabObservable`). The only other (currently implemented) noise distribution is log-normally distributed noise, which is designated through the `transformation=:log` argument:
```petab1
obs_P = PEtabObservable(P, σ; transformation=:log)
```

## [Additional features: Parameters](@id petab_parameters)

### [Known parameters](@id petab_parameters_known)
In our previous example, all parameters where unknowns that we wished to fit to the data. If any parameters have known values, it is possible to provide these to `PEtabModel` through the `parameter_map` argument. E.g if we had known that *kB = 1.0*, then we would only define *kD* and *kP* as parameters we wish to fit:
```petab1
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
parameters = [par_kD, par_kP]
nothing # hide
```
We then provide `parameter_map=[:kB => 1.0]` when we assembly our model:
```petab1
petab_model = PEtabModel(rn, observables, measurements, parameters; state_map=u0, parameter_map=[:kB => 1.0])
nothing # hide
```

### [Parameter bounds](@id petab_parameters_bounds)
By default, when fitted, potential parameter values are assumed to be in the interval *[1e-3, 1e3]*. When declaring a `PEtabParameter` is is possible to change these values through the `lb` and `ub` arguments. E.g. we could use
```petab1
par_kB = PEtabParameter(:kB; lb=1e-2, ub=1e-2)
```
to achieve the more conservative bound *[1e-2, 1e2]* for the parameter *kB*.

### [Parameter scales](@id petab_parameters_scales)

Internally, parameters that are fitted are converted to a logarithmic scale (this has several advantages, e.g. automatically preventing negative parameter values). To prevent this conversion, the `scale=:lin` argument can be used. Here we change the scale of *kB* to linear:
```petab1
par_kB = PEtabParameter(:kB; scale=:lin)
```

### [Parameter priors](@id petab_parameters_priors)
If one have prior knowledge about the distribution of a parameter, it is possible to incorporate this into the model. The prior can be any continuous distribution from the [Distributions.jl package](https://github.com/JuliaStats/Distributions.jl). E.g we can use:
```petab1
using Distributions
par_kB = PEtabParameter(:kB; prior=Normal(1.0,0.2))
```
to set a normally distributed prior (with mean `1.0` and standard deviation `0.2`) on the value of *kB*. By default, the prior is assumed to be on the linear scale of the parameter (before any potential log transform). To specify that the prior is on the logarithmic scale, the `prior_on_linear_scale=false` argument can be provided:
```petab1
par_kB = PEtabParameter(:kB; prior=Normal(1.0,0.2), prior_on_linear_scale=false)
```

## [Simulation conditions](@id petab_simulation_conditions)
Sometimes, we have data from several different experimental conditions. Here, when a potential parameter set is evaluated during the fitting process, each experimental condition corresponds to one simulation condition (which produce one simulation). To account for this, PEtab permits us to define several different simulation conditions, which each condition being defined by specific values for some initial conditions and/or parameters. 

If, for our previous catalysis example, we had measured the system for two different initial values of *S* (*S(0)=1.0* and *S(0)=0.5*), these would correspond to two different simulation conditions. For each condition we define a `Dict` mapping the species to their initial condition (here, *S* is the only species in each `Dict`):
```petab1
c1 = Dict(:S => 1.0)
c2 = Dict(:S => 0.5)
nothing # hide
```
Similarly as for observables, we then gather the conditions in another `Dict`, giving each an id tag:
```petab1
simulation_conditions = Dict("c1" => c1, "c2" => c2)
nothing # hide
```
Again (like for observables), each measurement in the measurements `DataFrame` needs to be associated with a simulation condition id tag (describing which condition that measurements was taken from). Parameters, just like initial conditions, may vary between different conditions. If an initial condition (or parameter) occurs in one condition, it must occur in all of them.

Here follows a complete version of our basic example, but with measurements both for *S(0)=1.0* and *S(0)=0.5*.
```petab3
using Catalyst, PEtab

rn = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end

u0 = [:E => 1.0, :SE => 0.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5]

# Simulate data.
using OrdinaryDiffEq, Random
d1, t1 = let 
    oprob = ODEProblem(rn, [:S=>1.0; u0], (0.0, 10.0), p_true)
    sol = solve(oprob, Tsit5(); saveat=0.1)
    ((0.8 .+ 0.4*rand(10)) .* sol[:P][10:10:end], sol.t[10:10:end])
end
d2, t2 = let 
    oprob = ODEProblem(rn, [:S=>0.5; u0], (0.0, 10.0), p_true)
    sol = solve(oprob, Tsit5(); saveat=0.1)
    ((0.8 .+ 0.4*rand(10)) .* sol[:P][10:10:end], sol.t[10:10:end])
end

@unpack P = rn
obs_P = PEtabObservable(P, 0.5)
observables = Dict("obs_P" => obs_P)

par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
parameters = [par_kB, par_kD, par_kP]

c1 = Dict(:kB => 1.0)
c2 = Dict(:kB => 0.5)
simulation_conditions = Dict("c1" => c1, "c2" => c2)

using DataFrames
m1 = DataFrame(simulation_id="c1", obs_id="obs_P", time=t1, measurement=d1)
m2 = DataFrame(simulation_id="c2", obs_id="obs_P", time=t2, measurement=d2)
measurements = vcat(m1,m2)

petab_model = PEtabModel(rn, simulation_conditions, observables, measurements, parameters; state_map=u0)
nothing # hide
```
Note that the `u0` we pass into `PEtabModel` through the `state_map` argument no longer contain the value of *S* (as it is provided by the conditions). 

### [Varying parameters between different simulation conditions](@id petab_simulation_conditions)
Sometimes, the parameters that are used vary between the different conditions. Consider our catalysis example, if we had performed the experiment twice, using two different enzymes with different catalytic properties, this could generate such conditions. Here, lets assume that the two enzymes yield different rates (*kP₁* and *kP₂*) for the `SE --> P + E` reaction (but are otherwise identical). 

[The author (Torkel) currently unsure how to carry this out, but will complete this tutorial when he knows.]

## [Additional features: Initial conditions](@id petab_simulation_initial_conditions)

### [Fitting initial conditions](@id petab_simulation_initial_conditions_fitted)
Sometimes, initial conditions are uncertain quantities which we wish to fit to the data. This is possible by [defining an initial condition as a parameter](dsl_description_parametric_initial_conditions): 
```petab4
using Catalyst, PEtab # hide
rn = @reaction_network begin
    @parameters E0
    @species E(t)=E0
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
nothing # hide
```
Here, the initial value of `E` is equal to the parameter `E0`. We modify our `u0` vector by removing `E` (which is no longer known):
```petab4
u0 = [:S => 1.0, :SE => 0.0, :P => 0.0]
nothing # hide
```
Next we add `E0` to the parameters we wish to fit:
```petab4
par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
par_E0 = PEtabParameter(:E0)
parameters = [par_kB, par_kD, par_kP, par_E0]
nothing # hide
```
and we can use our updated `rn`, `u0`, ad `parameters` as input to our `PEtabModel`.

### [Uncertain initial conditions](@id petab_simulation_initial_conditions_uncertainty)
Often, while an initial condition have been reported for an experiment, its exact value is uncertain. This can be modelled by making the initial condition a [parameter that is fitted to the data](@id petab_simulation_initial_conditions_fitted) and attach a prior to it corresponding to our certainty about its value. 

Let us consider our initial example, but where we want to add uncertainty to the initial conditions of `S` and `E`. We will add priors on these, assuming normal distributions with mean `1.0` and standard deviation `0.1`. For the synthetic measured data we will assume true values *S(0)=E(0)=1.0*.
```petab4
using Catalyst, PEtab

rn = @reaction_network begin
    @parameters S0 E0
    @species S(t)=S0 E(t)=E0
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end

u0 = [:SE => 0.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5]

# Simulate data.
using OrdinaryDiffEq, Random
oprob = ODEProblem(rn, [:S=>1.0; :E => 1.0; u0], (0.0, 10.0), p_true)
sol = solve(oprob, Tsit5(); saveat=0.1)
data = ((0.8 .+ 0.4*rand(10)) .* sol[:P][10:10:end], sol.t[10:10:end])

@unpack P = rn
obs_P = PEtabObservable(P, 0.5)
observables = Dict("obs_P" => obs_P)

par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
par_S0 = PEtabParameter(:S0)
par_E0 = PEtabParameter(:E0)
parameters = [par_kB, par_kD, par_kP, par_S0, par_E0]

using DataFrames
measurements = DataFrame(obs_id="obs_P", time=sol.t[10:10:end], measurement=data)

petab_model = PEtabModel(rn, observables, measurements, parameters; state_map=u0)
nothing # hide
```
Here, when we fit our data we will also gain values for `S0` and `E0`, however, unless we are explicitly interested in these, they can be ignored.

## [Additional features: Simulation options](@id petab_simulation_options)

While in our basic example, we do not provide any additional information to our `PEtabODEProblem`, this is an opportunity to specify how the model should be simulated, and what differentiation techniques to use for the optimisation procedure (if none are provided, appropriate defaults are selected). 

Here is an example, taken from the [more detailed PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/dev/Boehm/#Creating-a-PEtabODEProblem)
```petab1
petab_problem = PEtabODEProblem(petab_model, 
                                ode_solver=ODESolver(Rodas5P(), abstol=1e-8, reltol=1e-8), 
                                gradient_method=:ForwardDiff, 
                                hessian_method=:ForwardDiff)
nothing # hide
```
where we simulate our ODE model using the `Rodas5p` method (with absolute and relative tolerance both equal `1e-8`) and use [forward automatic differentiation](https://github.com/JuliaDiff/ForwardDiff.jl) for both gradient and hessian computation. More details on available ODE solver options can be found in [the PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/dev/API_choosen/#PEtab.ODESolver).

## [Additional features: Optimisation](@id petab_optimisation)

### [Optimisation methods and options](@id petab_optimisation_optimisers)
For our examples, we have used the `Optim.IPNewton` optimisation method. PEtab.jl supports [several additional optimisation methods](https://sebapersson.github.io/PEtab.jl/dev/Avaible_optimisers/#options_optimisers). Furthermore, the `calibrate_model` `options` argument permits the customisation of the options for any used optimiser. E.g. to increase the maximum number of iterations of the `Optim.IPNewton` method we would use:
```petab1
res = calibrate_model(petab_problem, p0, IPNewton(); options=Optim.Options(iterations = 10000))
nothing # hide
```

Please read the [PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/dev/Avaible_optimisers/#options_optimisers) to learn how to customise the various optimisers' properties. 

### Optimisation path recording
To record all the parameter sets evalauted (and their objective values) during the optimisation procedure, provide the `save_trace=true` argument to either `calibrate_model` or `calibrate_model_multistart`. This is required for the various [optimisation evaluation plots](@ref petab_plotting) provided by PEtab.jl.


## Objective function extraction
While PEtab.jl provides various tools for analysing the objective function generated by `PEtabODEProblem`, it is also possible to extract this function for customised analysis. Given a `PEtabODEProblem`
```petab1
petab_problem = PEtabODEProblem(petab_model)
nothing # hide
```
We can find the:
1. Objective function in the `petab_problem.compute_cost`. It takes a single argument (`p`) and return the objective value.
2. Gradient in the `petab_problem.compute_gradient!` field. It takes two arguments (`g` and `p`) with `g` containing the updated gradient values.
3. Hessian in the `petab_problem.compute_hessian` field. It takes two arguments (`H` and `p`) with `H` containing the updated hessian values.

## [Multistart optimisation](@id petab_multistart_optimisation)
To avoid the optimisation process returning a local minimum, it is often advice to run it multiple, using different initial guesses. PEtab.jl supports this through the `calibrate_model_multistart` function. This is identical to the `calibrate_model` function, but takes two additional arguments:
1. `n_multistart`: The number of runs to perform.
2. `dir_save`: A location to which the output is saved. If `dir_save=nothing`, no output is saved.

And one additional optional argument:

3. `sampling_method`: Selection the sampling method with which to select the initial guesses (`QuasiMonteCarlo.LatinHypercubeSample()` used by default). 
 
Because `calibrate_model_multistart` handles initial guess sampling, unlike for `calibrate_model`, no initial guess have to be provided. 

Here, we fits parameters by through 10 independent optimisation runs, using QuasiMonteCarlo's `SobolSample` method, and saves the result to the "OptimizationRuns/runs_ms_1" folder:
```petab1
using Optim
using QuasiMonteCarlo
res_ms = calibrate_model_multistart(petab_problem, IPNewton(), 10, "OptimizationRuns/runs_ms_1"; sampling_method=QuasiMonteCarlo.SobolSample())
```
The best result across all runs can still be retrieved using `get_best_ps(res_ms)`, with the results of the individual runs being stored in the `res_ms.runs` field.

## [Events](@id petab_events)
Support for events in PEtab and PEtab.jl is currently a work in progress, with an implementation in PEtab.jl expected soon.


## [Plot recipes](@id petab_plotting)
There exist various types of graphs that can be used to evaluate the parameter fitting process. Support for these in PEtab.jl is currently a work in progress, and they should be available soon.


## [Citations](@id petab_citations)
If you use this functionality in your research, please cite the following papers to support the authors of the PEtab.jl package (currently there is no article associated with this package) and the PEtab standard:
```
@misc{2023Petabljl,
<<<<<<< refs/remotes/origin/master
  author       = {Ognissanti, Damiano and Arutjunjan, Rafael and Persson, Sebastian and Hasselgren, Viktor.
    
    Isaacson, S. A. and Ilin, V. and Rackauckas, C. V.},
>>>>>>> initiate tutorial
=======
  author       = {Ognissanti, Damiano and Arutjunjan, Rafael and Persson, Sebastian and Hasselgren, Viktor},
>>>>>>> init tutorial
  title        = {{2023Petabljl.jl}},
  howpublished = {\url{https://github.com/sebapersson/PEtab.jl/}},
  year         = {2023}
}
<<<<<<< refs/remotes/origin/master
<<<<<<< refs/remotes/origin/master
```
```
=======
>>>>>>> initiate tutorial
=======
```
```
>>>>>>> update
@article{SchmiesterSch2021,
  author    = {Schmiester, Leonard AND Schälte, Yannik AND Bergmann, Frank T. AND Camba, Tacio AND Dudkin, Erika AND Egert, Janine AND Fröhlich, Fabian AND Fuhrmann, Lara AND Hauber, Adrian L. AND Kemmer, Svenja AND Lakrisenko, Polina AND Loos, Carolin AND Merkt, Simon AND Müller, Wolfgang AND Pathirana, Dilan AND Raimúndez, Elba AND Refisch, Lukas AND Rosenblatt, Marcus AND Stapor, Paul L. AND Städter, Philipp AND Wang, Dantong AND Wieland, Franz-Georg AND Banga, Julio R. AND Timmer, Jens AND Villaverde, Alejandro F. AND Sahle, Sven AND Kreutz, Clemens AND Hasenauer, Jan AND Weindl, Daniel},
  journal   = {PLOS Computational Biology},
  title     = {PEtab—Interoperable specification of parameter estimation problems in systems biology},
  year      = {2021},
  month     = {01},
  number    = {1},
  pages     = {1-10},
  volume    = {17},
  doi       = {10.1371/journal.pcbi.1008646},
  publisher = {Public Library of Science},
  url       = {https://doi.org/10.1371/journal.pcbi.1008646},
}
```
<<<<<<< refs/remotes/origin/master
<<<<<<< refs/remotes/origin/master

## References
[1]: [Schmiester, L et al. *PEtab—Interoperable specification of parameter estimation problems in systems biology*, PLOS Computational Biology (2021).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008646)
=======
=======

>>>>>>> init tutorial
## References
<<<<<<< refs/remotes/origin/master
[^1]: [Schmiester, L et al. *PEtab—Interoperable specification of parameter estimation problems in systems biology*, PLOS Computational Biology (2021).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008646)
>>>>>>> initiate tutorial
=======
[1]: [Schmiester, L et al. *PEtab—Interoperable specification of parameter estimation problems in systems biology*, PLOS Computational Biology (2021).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008646)
>>>>>>> update
