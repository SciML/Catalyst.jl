# [Parameter Fitting for ODEs using PEtab.jl](@id petab_parameter_fitting)
The [PEtab.jl package](https://github.com/sebapersson/PEtab.jl) implements the [PEtab format](https://petab.readthedocs.io/en/latest/) for fitting the parameters of deterministic CRN models to data [^1]. PEtab.jl both implements methods for creating cost functions (determining how well parameter sets fit to data), and for minimizing these cost functions. The PEtab approach covers most cases of fitting deterministic (ODE) models to data and is a good default choice when fitting reaction rate equation ODE models. This page describes how to combine PEtab.jl and Catalyst for parameter fitting, with the PEtab.jl package providing [a more extensive documentation](https://sebapersson.github.io/PEtab.jl/stable/) (this tutorial is partially an adaptation of this documentation). 

While PEtab's interface generally is very flexible, there might be specific use-cases where it cannot create an appropriate cost-function. Here, it is recommended to instead look at using [Optimization.jl](@ref optimization_parameter_fitting).

## Introductory example
Let us consider a simple catalysis network, where an enzyme ($E$) turns a substrate ($S$) into a product ($P$):
```@example petab1
using Catalyst, PEtab

rn = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
```
From some known initial condition, and a true parameter set (which we later want to recover from the data) we generate synthetic data (on which we will demonstrate the fitting process).
```@example petab1
# Define initial conditions and parameters.
u0 = [:S => 1.0, :E => 1.0, :SE => 0.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5]

# Generate synthetic data.
using OrdinaryDiffEqDefault
oprob_true = ODEProblem(rn, u0, (0.0, 10.0), p_true)
true_sol = solve(oprob_true)
data_sol = solve(oprob_true; saveat = 1.0)
data_ts = data_sol.t[2:end]
data_vals = (0.8 .+ 0.4*rand(10)) .* data_sol[:P][2:end]

# Plots the true solutions and the (synthetic data) measurements.
using Plots
default(bottom_margin = 4Plots.Measures.mm, left_margin = 4Plots.Measures.mm) # hide
plot(true_sol; idxs = :P, label = "True solution", lw = 4)
plot!(data_ts, data_vals; label = "Measurements", seriestype = :scatter, ms = 6, color = :blue)
```

Generally, PEtab takes five different inputs to define an optimisation problem (the minimiser of which generates a fitted parameter set):
1. **Model**: The model which we want to fit to the data (a `ReactionSystem`).
2. **Observables**: The possible observables that can be measured (a `Dict{String,PEtabObservable}`).
3. **Estimation parameters**: The parameters which we want to fit to the data (a `Vector{PEtabParameter}`).
4. **Experimental (or simulation) conditions**: The simulations (each corresponding to a potential experiment) carried out during each step of the optimisation process (a `Dict{String,Dict{Symbol,Float64}}`).
5. **Measurements**: The measurements to which the model is fitted (a `DataFrame`).

### Observables
The observables define the quantities that we may measure in our experiments. Typically, each corresponds to a single species, however, [more complicated observables are possible](@ref petab_observables_observables). For each observable, we also need a noise formula, defining the uncertainty in its measurements. By default, PEtab assumes normally distributed noise, with a mean equal to the true value and a standard deviation which we have to define. It is also possible to use [more advanced noise formulas](@ref petab_observables_noise_formula).

In our example, we only have a single possible observable, the `P` species. We will assume that the noise is normally distributed with a standard deviation `0.5` (in our case this is not true, however, typically the noise distribution is unknown and a guess must be made). We combine this information in a `PEtabObservable` struct (to access the `P` species we must use [`@unpack`](@ref simulation_structure_interfacing_symbolic_representation)). Finally, we store all our observables in a dictionary, giving each an id tag (which is later used in the measurements input).

```@example petab1
@unpack P = rn
obs_P = PEtabObservable(P, 0.5)
observables = Dict("obs_P" => obs_P)
nothing # hide
```

### Parameters
Each parameter of the system can either be
1. Known (described [here](@ref petab_parameters_known)).
2. Depend on experimental/simulation conditions (described [here](@ref petab_simulation_conditions)).
3. Be an unknown that we wish to fit to data.

In our case, we wish to fit all three system parameters ($kB$, $kD$, and $kP$). For each, we create a single `PEtabParameter`, and then gather these into a single vector.
```@example petab1
par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
params = [par_kB, par_kD, par_kP]
nothing # hide
```
For each parameter, it is also possible to set [a lower and/or upper bound](@ref petab_parameters_bounds) (by default, $(0.001,1000)$ is used), set whether to use [logarithmic or linear scale](@ref petab_parameters_scales), or add a [prior distribution of its value](@ref petab_parameters_priors).

### Simulation conditions
Sometimes, several different experiments are performed on a system (each potentially generating several measurements). An experiment could e.g. be the time development of a system from a specific initial condition. Since each experimental condition (during the optimisation procedure, for a guess of the unknown parameters) generates a distinct simulation, these are also called simulation conditions. In our example, all data comes from a single experiment, and the simulation condition input is not required. How to define and use different experimental conditions is described [here](@ref petab_simulation_conditions).

### Measurements
Finally, we need to define the system measurements to which the parameters will be fitted. Each measurement combines:
1. The observable which is observed (here we use the id tag defined in the `observables` dictionary).
2. The time point of the measurement.
3. The measured value.

For cases where several simulation conditions are given, we also need to provide:

4. The simulation condition which generates the measurement ([here](@ref petab_simulation_conditions) is an example where this is used).

*Note also, when [pre-equilibration](https://sebapersson.github.io/PEtab.jl/stable/Brannmark/) is used to initiate the system in a steady state, a fifth field is also required.*

Each individual measurement is provided as a row of a `DataFrame`. The values are provided in the `obs_id`, `time`, `measurement`, and `simulation_id` columns. In our case we only need to fill in the first three:
```@example petab1
using DataFrames
measurements = DataFrame(obs_id = "obs_P", time = data_ts, measurement = data_vals)
```
Since, in our example, all measurements are of the same observable, we can set `obs_id="obs_P"`. However, it is also possible to [include measurements from several different observables](@ref petab_simulation_measurements_several_observables).


### Creating a PEtabModel
Finally, we combine all inputs in a single `PEtabModel`. To it, we also pass the initial conditions of our simulations (using the `speciemap` argument). It is also possible to have [initial conditions with uncertainty](@ref petab_simulation_initial_conditions_uncertainty), [that vary between different simulations](@ref petab_simulation_conditions), or [that we attempt to fit to the data](@ref petab_simulation_initial_conditions_fitted).
```@example petab1
petab_model = PEtabModel(rn, observables, measurements, params; speciemap = u0)
nothing # hide
```

### Fitting parameters
We are now able to fit our model to the data. First, we create a `PEtabODEProblem`. Here, we use `petab_model` as the only input, but it is also possible to set various [numeric solver and automatic differentiation options](@ref petab_simulation_options) (such as method or tolerance).
```@example petab1
petab_problem = PEtabODEProblem(petab_model)
```
Since no additional input is given, default options are selected by PEtab.jl (and generally, its choices are good).

To fit a parameter set we use the `calibrate` function. In addition to our `PEtabODEProblem`, we must also provide an initial guess (which can be generated with the `generate_startguesses` function) and an optimisation algorithm (which needs to be imported specifically). PEtab.jl supports various [optimisation methods and options](@ref petab_optimisation_optimisers).
```@example petab1
using Optim
p0 = get_startguesses(petab_problem, 1)
p0 = [0.0, 0.0, 0.0] # hide
res = calibrate(petab_problem, p0, IPNewton()) # hide
res = calibrate(petab_problem, p0, IPNewton())
```

We can now simulate our model for the fitted parameter set, and compare the result to the measurements and true solution.
```@example petab1
oprob_fitted = remake(oprob_true; p = get_ps(res, petab_problem))
fitted_sol = solve(oprob_fitted)
plot!(fitted_sol; idxs = :P, label = "Fitted solution", linestyle = :dash, lw = 4, color = :lightblue)
```

Here we use the `get_ps` function to retrieve a full parameter set using the optimal parameters. Alternatively, the `ODEProblem` or fitted simulation can be retrieved directly using the `get_odeproblem` or `get_odesol` [functions](https://sebapersson.github.io/PEtab.jl/stable/API/), respectively (and the initial condition using the `get_u0` function). The calibration result can also be found in `res.xmin`, however, note that PEtab automatically ([unless a linear scale is selected](@ref petab_parameters_scales)) converts parameters to logarithmic scale, so typically `10 .^res.xmin` are the values of interest. If you investigate the result from this example you might note, that even if PEtab.jl has found the global optimum (which fits the data well), this does not actually correspond to the true parameter set. This phenomenon is related to the concept of *identifiability*, which is very important for parameter fitting.

### Final notes
PEtab.jl also supports [multistart optimisation](@ref petab_multistart_optimisation), [automatic pre-equilibration before simulations](https://sebapersson.github.io/PEtab.jl/stable/petab_preeq_simulations/), and [events](@ref petab_events). Various [plot recipes](@ref petab_plotting) exist for investigating the optimisation process. Please read the [PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/) for a more complete description of the package's features. Below follows additional details of various options and features (generally, PEtab is able to find good default values for most options that are not specified).

## [Additional features: Observables](@id petab_observables)

### [Defining non-trivial observables](@id petab_observables_observables)
It is possible for observables to be any algebraic expression of species concentrations and parameters. E.g. in this example the total amount of `X` in the system is an observable:
```@example petab2
using Catalyst, PEtab # hide
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
@unpack X1, X2 = two_state_model
obs_X = PEtabObservable(X1 + X2, 0.5)
```

A common application for this is to define an [*offset* and a *scale* for each observable](https://sebapersson.github.io/PEtab.jl/stable/petab_obs_noise/).

### [Advanced observables noise formulas](@id petab_observables_noise_formula)
In our basic example we assumed that the normally distributed noise had a standard deviation of `0.5`. However, this value may be a parameter (or indeed any algebraic expression). E.g, we could set
```@example petab1
@parameters σ
obs_P = PEtabObservable(P, σ)
```
and then let the parameter $σ$ vary between different [simulation conditions](@ref petab_simulation_conditions). Alternatively we could let the noise scale linearly with the species' concentration:
```@example petab1
obs_P = PEtabObservable(P, 0.05P)
```
It would also be possible to make $σ$ a *parameter that is fitted as a part of the parameter fitting process*.

PEtab.jl assumes that noise is normally distributed (with a standard deviation equal to the second argument provided to `PEtabObservable`). The only other (currently implemented) noise distribution is log-normally distributed noise, which is designated through the `transformation=:log` argument:
```@example petab1
obs_P = PEtabObservable(P, σ; transformation = :log)
```

## [Additional features: Parameters](@id petab_parameters)

### [Known parameters](@id petab_parameters_known)
In our previous example, all parameters were unknowns that we wished to fit to the data. If any parameters have known values, it is possible to provide these to `PEtabModel` through the `parameter_map` argument. E.g if we had known that $kB = 1.0$, then we would only define $kD$ and $kP$ as parameters we wish to fit:
```@example petab1
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
params = [par_kD, par_kP]
nothing # hide
```
We then provide `parameter_map=[:kB => 1.0]` when we assembly our model:
```@example petab1
petab_model_known_param = PEtabModel(rn, observables, measurements, params; speciemap = u0, parametermap = [:kB => 1.0])
nothing # hide
```

### [Parameter bounds](@id petab_parameters_bounds)
By default, when fitted, potential parameter values are assumed to be in the interval $(1e-3, 1e3)$. When declaring a `PEtabParameter` it is possible to change these values through the `lb` and `ub` arguments. E.g. we could use
```@example petab1
par_kB = PEtabParameter(:kB; lb = 1e-2, ub = 1e2)
```
to achieve the more conservative bound $(1e-2, 1e2)$ for the parameter $kB$.

### [Parameter scales](@id petab_parameters_scales)

Internally, parameters that are fitted are converted to a logarithmic scale (generally, this is considered advantageous[^2]). To prevent this conversion, the `scale=:lin` argument can be used. Here we change the scale of $kB$ to linear:

```@example petab1
par_kB = PEtabParameter(:kB; scale = :lin)
```

### [Parameter priors](@id petab_parameters_priors)
If we have prior knowledge about the distribution of a parameter, it is possible to incorporate this into the model. The prior can be any continuous, univariate, distribution from the [Distributions.jl package](https://github.com/JuliaStats/Distributions.jl). E.g we can use:

```@example petab1
using Distributions
par_kB = PEtabParameter(:kB; prior = Normal(1.0,0.2))
```
to set a normally distributed prior (with mean `1.0` and standard deviation `0.2`) on the value of $kB$. By default, the prior is assumed to be on the linear scale of the parameter (before any potential log transform). To specify that the prior is on the logarithmic scale, the `prior_on_linear_scale=false` argument can be provided:
```@example petab1
par_kB = PEtabParameter(:kB; prior = Normal(1.0,0.2), prior_on_linear_scale = false)
```
In this example, setting `prior_on_linear_scale=false` makes sense as a (linear) normal distribution is non-zero for negative values (an alternative is to use a log-normal distribution, e.g. `prior=LogNormal(3.0, 3.0)`).

## [Simulation conditions](@id petab_simulation_conditions)
Sometimes, we have data from different experimental conditions. Here, when a potential parameter set is evaluated during the fitting process, each experimental condition corresponds to one simulation condition (which produces one simulation). To account for this, PEtab permits the user to define different simulation conditions, with each condition being defined by specific values for some initial conditions and/or parameters.

If, for our previous catalysis example, we had measured the system for two different initial values of $S$ ($S(0)=1.0$ and $S(0)=\tfrac{1}{2}$), these would correspond to two different simulation conditions. For each condition we define a `Dict` mapping the species to their initial condition (here, $S$ is the only species in each `Dict`):
```@example petab1
c1 = Dict(:S => 1.0)
c2 = Dict(:S => 0.5)
nothing # hide
```
Similarly as for observables, we then gather the conditions in another `Dict`, giving each an id tag:
```@example petab1
simulation_conditions = Dict("c1" => c1, "c2" => c2)
nothing # hide
```
Again (like for observables), each measurement in the measurements `DataFrame` needs to be associated with a simulation condition id tag (describing which condition those measurements were taken from). Parameters, just like initial conditions, may vary between different conditions. If an initial condition (or parameter) occurs in one condition, it must occur in all of them.

Here follows a complete version of our basic example, but with measurements both for $S(0)=1.0$ and $S(0)=\tfrac{1}{2}$.
```@example petab3
using Catalyst, PEtab

rn = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end

u0 = [:E => 1.0, :SE => 0.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5]

# Simulate data.
using OrdinaryDiffEqDefault
t1, d1 = let
    oprob_true = ODEProblem(rn, [:S => 1.0; u0], (0.0, 10.0), p_true)
    data_sol = solve(oprob_true; saveat = 1.0)
    data_sol.t[2:end], (0.8 .+ 0.4*rand(10)) .* data_sol[:P][2:end]
end
t2, d2 = let
    oprob_true = ODEProblem(rn, [:S=>0.5; u0], (0.0, 10.0), p_true)
    data_sol = solve(oprob_true; saveat = 1.0)
    data_sol.t[2:end], (0.8 .+ 0.4*rand(10)) .* data_sol[:P][2:end]
end

@unpack P = rn
obs_P = PEtabObservable(P, 0.5)
observables = Dict("obs_P" => obs_P)

par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
params = [par_kB, par_kD, par_kP]

c1 = Dict(:S => 1.0)
c2 = Dict(:S => 0.5)
simulation_conditions = Dict("c1" => c1, "c2" => c2)

using DataFrames
m1 = DataFrame(simulation_id = "c1", obs_id = "obs_P", time = t1, measurement = d1)
m2 = DataFrame(simulation_id = "c2", obs_id = "obs_P", time = t2, measurement = d2)
measurements = vcat(m1,m2)

petab_model = PEtabModel(rn, observables, measurements, params; speciemap = u0, 
                         simulation_conditions = simulation_conditions)
nothing # hide
```
Note that the `u0` we pass into `PEtabModel` through the `speciemap` argument no longer contains the value of $S$ (as it is provided by the conditions).

## [Additional features: Measurements](@id petab_simulation_measurements)

### [Measurements of several observables](@id petab_simulation_measurements_several_observables)
In our previous example, all our measurements were from a single observable, `obs_P`. If we also had collected measurements of both $S$ and $P$:
```@example petab1
data_ts = data_sol.t[2:end]
data_vals_S = (0.8 .+ 0.4*rand(10)) .* data_sol[:S][2:end]
data_vals_P = (0.8 .+ 0.4*rand(10)) .* data_sol[:P][2:end]
nothing # hide
```
and then corresponding observables:
```@example petab1
@unpack S, P = rn
obs_S = PEtabObservable(S, 0.5)
obs_P = PEtabObservable(P, 0.5)
observables = Dict("obs_S" => obs_P, "obs_P" => obs_P)
nothing # hide
```
we are able to include all these measurements in the same `measurements` `DataFrame`:
```@example petab1
m1 = DataFrame(obs_id = "obs_P", time = data_ts, measurement = data_vals_S)
m2 = DataFrame(obs_id = "obs_S", time = data_ts, measurement = data_vals_P)
measurements = vcat(m1,m2)
```
which then can be used as input to `PEtabModel`.

### Varying parameters between different simulation conditions
Sometimes, the parameters that are used vary between the different conditions. Consider our catalysis example, if we had performed the experiment twice, using two different enzymes with different catalytic properties, this could have generated such conditions. The two enzymes could e.g. yield different rates ($kP_1$ and $kP_2$) for the `SE --> P + E` reaction, but otherwise be identical. Here, the parameters $kP_1$ and $kP_2$ are unique to their respective conditions. PEtab.jl provides support for cases such as this, and [its documentation](https://sebapersson.github.io/PEtab.jl/stable/petab_cond_specific/) provided instructions of how to handle them.

## [Additional features: Initial conditions](@id petab_simulation_initial_conditions)

### [Fitting initial conditions](@id petab_simulation_initial_conditions_fitted)
Sometimes, initial conditions are uncertain quantities which we wish to fit to the data. This is possible [by defining an initial condition as a parameter](@ref dsl_advanced_options_parametric_initial_conditions):
```@example petab4
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
```@example petab4
u0 = [:S => 1.0, :SE => 0.0, :P => 0.0]
nothing # hide
```
Next, we add `E0` to the parameters we wish to fit:
```@example petab4
par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
par_E0 = PEtabParameter(:E0)
params = [par_kB, par_kD, par_kP, par_E0]
nothing # hide
```
and we can use our updated `rn`, `u0`, and `params` as input to our `PEtabModel`.

### [Uncertain initial conditions](@id petab_simulation_initial_conditions_uncertainty)
Often, while an initial condition has been reported for an experiment, its exact value is uncertain. This can be modelled by making the initial condition a [parameter that is fitted to the data](@ref petab_simulation_initial_conditions_fitted) and attaching a prior to it corresponding to our certainty about its value.

Let us consider our initial example, but where we want to add uncertainty to the initial conditions of `S` and `E`. We will add priors on these, assuming normal distributions with mean `1.0` and standard deviation `0.1`. For the synthetic measured data we will use the true values $S(0) = E(0) = 1.0$.
```@example petab5
using Catalyst, Distributions, PEtab 

rn = @reaction_network begin
    @parameters S0 E0
    @species S(t)=S0 E(t)=E0
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end

u0 = [:SE => 0.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5, :S0=>1.0, :E0 => 1.0]

using OrdinaryDiffEqDefault
oprob_true = ODEProblem(rn, u0, (0.0, 10.0), p_true)
true_sol = solve(oprob_true)
data_sol = solve(oprob_true; saveat = 1.0)
data_ts = data_sol.t[2:end]
data_vals = (0.8 .+ 0.4*rand(10)) .* data_sol[:P][2:end]

@unpack P = rn
obs_P = PEtabObservable(P, 0.5)
observables = Dict("obs_P" => obs_P)

par_kB = PEtabParameter(:kB)
par_kD = PEtabParameter(:kD)
par_kP = PEtabParameter(:kP)
par_S0 = PEtabParameter(:S0; prior = Normal(1.0, 0.1))
par_E0 = PEtabParameter(:E0; prior = Normal(1.0, 0.1))
params = [par_kB, par_kD, par_kP, par_S0, par_E0]

using DataFrames
measurements = DataFrame(obs_id = "obs_P", time = data_ts, measurement = data_vals)

petab_model = PEtabModel(rn, observables, measurements, params; speciemap = u0)
nothing # hide
```
Here, when we fit our data we will also gain values for `S0` and `E0`, however, unless we are explicitly interested in these, they can be ignored.

## [Additional features: Simulation options](@id petab_simulation_options)

While in our basic example, we do not provide any additional information to our `PEtabODEProblem`, this is an opportunity to specify how the model should be simulated, and what automatic differentiation techniques to use for the optimisation procedure (if none are provided, appropriate defaults are selected).

Here is an example, adapted from the [more detailed PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/default_options/)
```@example petab1
using OrdinaryDiffEqRosenbrock
PEtabODEProblem(petab_model, 
                odesolver = ODESolver(Rodas5P(), abstol = 1e-8, reltol = 1e-8),
                gradient_method = :ForwardDiff, 
                hessian_method = :ForwardDiff)
nothing # hide
```
where we simulate our ODE model using the `Rodas5P` method (with absolute and relative tolerance both equal `1e-8`) and use [forward automatic differentiation](https://github.com/JuliaDiff/ForwardDiff.jl) for both gradient and hessian computation. More details on available ODE solver options can be found in [the PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/API/#PEtabODEProblem).

## [Additional features: Optimisation](@id petab_optimisation)

### [Optimisation methods and options](@id petab_optimisation_optimisers)
For our examples, we have used the `Optim.IPNewton` optimisation method. PEtab.jl supports [several additional optimisation methods](https://sebapersson.github.io/PEtab.jl/stable/pest_algs/). Furthermore, `calibrate`'s `options` argument permits the customisation of the options for any used optimiser. E.g. to designate the maximum number of iterations of the `Optim.IPNewton` method we would use:
```@example petab1
res = calibrate(petab_problem, p0, IPNewton(); options = Optim.Options(iterations = 10000))
nothing # hide
```

Please read the [PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/pest_algs/) to learn how to customize the various optimisers' properties.

### [Optimisation path recording](@id petab_optimisation_path_recording)

To record all the parameter sets evaluated (and their objective values) during the optimisation procedure, provide the `save_trace=true` argument to `calibrate` (or `calibrate_multistart`):
```@example petab1
res = calibrate(petab_problem, p0, IPNewton(); save_trace = true)
nothing # hide
```
This is required for the various [optimisation evaluation plots](@ref petab_plotting) provided by PEtab.jl. If desired, this information can be accessed in the calibration output's `.xtrace` and `.ftrace` fields.

## Objective function extraction
While PEtab.jl provides various tools for analysing the objective function generated by `PEtabODEProblem`, it is also possible to extract this function for customised analysis. Given a `PEtabODEProblem`
```@example petab1
petab_problem = PEtabODEProblem(petab_model)
nothing # hide
```
```julia
petab_problem = PEtabODEProblem(petab_model)
```
We can find the:
1. Objective function (negative log-likelihood) as the `petab_problem.nllh`. It takes a single argument (`p`) and returns the objective value.
2. Gradient as the `petab_problem.grad!` field. It takes two arguments (`g` and `p`) with the updated gradient values being written to `g`.
3. Hessian as the `petab_problem.hess!` field. It takes two arguments (`H` and `p`) with the updated hessian values being written to `H`.

## [Multi-start optimisation](@id petab_multistart_optimisation)
To avoid the optimisation process returning a local minimum, it is often advised to run it multiple times, using different initial guesses. PEtab.jl supports this through the `calibrate_multistart` function. This is identical to the `calibrate` function, but takes one additional arguments:

1. `nmultistarts`: The number of runs to perform.

And two additional optional argument:

2. `dirsave`: A location to which the output is automatically saved. If `dirsave=nothing`, no output is saved. It is recommended to save intermediate results for parameter estimation runs that take a long time, to not lose intermediate results if something goes wrong.
3. `sampling_method`: Selects the sampling method with which to select the initial guesses (`QuasiMonteCarlo.LatinHypercubeSample()` used by default).

Because `calibrate_multistart` handles initial guess sampling, unlike for `calibrate`, no initial guess has to be provided.

Here, we fit parameters through 10 independent optimisation runs, using QuasiMonteCarlo's `SobolSample` method, and save the result to the OptimisationRuns folder:
```@example petab1
using Optim
using QuasiMonteCarlo
mkdir("OptimisationRuns") # hide
res_ms = calibrate_multistart(petab_problem, IPNewton(), 10; dirsave = "OptimisationRuns",
    sampling_method = QuasiMonteCarlo.SobolSample())
res_ms = calibrate_multistart(petab_problem, IPNewton(), 10; dirsave = "OptimisationRuns", sampling_method = QuasiMonteCarlo.SobolSample()) # hide
nothing # hide
```
The best result across all runs can still be retrieved using `get_ps(res_ms, petab_problem)`, with the results of the individual runs being stored in the `res_ms.runs` field.

To load the result in a later session, we can call:
```@example petab1
res_ms = PEtabMultistartResult("OptimisationRuns")
nothing # hide
```
where `"OptimisationRuns"` is the name of the save directory (specified in `calibrate_multistart`). If the OptimisationRuns folder contains the output from several runs, we can designate which to load using the `which_run` argument. Here we load the second run to be saved in that folder:
```@example petab1
res_ms = PEtabMultistartResult("OptimisationRuns"; which_run = 2)
rm("OptimisationRuns", recursive = true) # hide
nothing # hide
```
By default, `which_run` loads the first run saved to that directory.

## [Events](@id petab_events)
So far, we have assumed that all experiments, after initiation, run without interference. Experiments where conditions change, or where species are added/removed during the time course, can be represented through events. In PEtab, an event is represented through the `PEtabEvent` structure. It takes three arguments:
1. The condition for triggering the event. This can either indicate a point in time, or a boolean condition.
2. A rule for updating the event's target
3. The event's target (either a species or parameter).

Here we create an event which adds `0.5` units of `S` to the system at time `5.0`:
```@example petab1
@unpack S = rn
event1 = PEtabEvent(5.0, S + 0.5, S)
```
Here, the first argument is evaluated to a scalar value, in which case it is interpreted as a time point at which the event happens. If we instead want the event to happen whenever the concentration of `S` falls below `0.5` we set the first argument to a boolean condition indicating this:
```@example petab1
event2 = PEtabEvent(S < 0.5, S + 0.5, S)
```
Here, the event only triggers whenever the condition changes from `false` to `true`, and not while it remains `true` (or when changing from `true` to `false`). E.g. this event only triggers when `S` concentration passes from more than $5.0$ to less that $5.0$.

Whenever we have several events or not, we bundle them together in a single vector which is later passed to the `PEtabODEProblem` using the `events` argument:
```@example petab1
params = [par_kB, par_kD, par_kP] # hide
events = [event1, event2]
petab_model = PEtabModel(rn, observables, measurements, params; speciemap = u0, events = events)
nothing # hide
```

More details on how to use events, including how to create events with multiple targets, can be found in [PEtab.jl's documentation](https://sebapersson.github.io/PEtab.jl/stable/petab_event/).

!!! note
    PEtab currently ignores events [created as a part of a Catalyst `ReactionSystem` model](@ref events), and does not support SciML-style events. Instead, events have to use the preceding interface.

## [Plot recipes](@id petab_plotting)
There exist various types of graphs that can be used to evaluate the parameter fitting process. These can be plotted using the `plot` command, where the input is either the result of a `calibrate` or a `calibrate_multistart` run. To be able to use this functionality, you have to ensure that PEtab.jl [records the optimisation process](@ref petab_optimisation_path_recording) by providing the `save_trace=true` argument to the calibration functions.

To, for a single start calibration run, plot, for each iteration of the optimization process, the best objective value achieved so far, run:
```@example petab1
res = calibrate(petab_problem, p0, IPNewton(); save_trace = true) # hide
plot(res)
```

For a multi-start calibration run, the default output is instead a so-called waterfall plot:
```@example petab1
res_ms = PEtabMultistartResult("../assets/boehm___for_petab_tutorial") # hide
plot(res_ms)
```

(for this, and the next plot, we use a multi-start optimisation result from a different model, which yields less trivial optimisation runs than our catalysis one)

In the waterfall plot, each dot shows the final objective value for a single run in the multi-start process. The runs are ordered by their objective values, and colours designate runs in the same local minimum. A common use of waterfall plots is to check whether a sufficient number of runs (typically $>5$) has converged to the same best local minimum (in which case it is assumed to be the *global* minimum).

To instead use the best objective value plot for a multi-start run (with one curve for each run), the `plot_type` argument is used:
```@example petab1
plot(res_ms; plot_type = :best_objective)
```

There exist several types of plots for both types of calibration results. More details of the types of available plots, and how to customise them, can be found [here](https://sebapersson.github.io/PEtab.jl/stable/optimisation_output_plotting/).

---
## [Citations](@id petab_citations)
If you use this functionality in your research, [in addition to Catalyst](@ref doc_index_citation), please cite the following papers to support the authors of the PEtab.jl package (currently there is no article associated with this package) and the PEtab standard:
```
@misc{2023Petabljl,
  author       = {Ognissanti, Damiano AND Arutjunjan, Rafael AND Persson, Sebastian AND Hasselgren, Viktor},
  title        = {{2023Petabljl.jl}},
  howpublished = {\url{https://github.com/sebapersson/PEtab.jl/}},
  year         = {2023}
}
```
```
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

---
## References
[^1]: [Schmiester, L et al. *PEtab—Interoperable specification of parameter estimation problems in systems biology*, PLOS Computational Biology (2021).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008646)
[^2]: [Hass, H et al. *Benchmark problems for dynamic modeling of intracellular processes*, Bioinformatics (2019).](https://academic.oup.com/bioinformatics/article/35/17/3073/5280731?login=false)
