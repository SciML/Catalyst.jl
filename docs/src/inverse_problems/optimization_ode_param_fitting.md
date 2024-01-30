# [Parameter Fitting for ODEs using SciML/Optimization.jl and DiffEqParamEstim.jl](@id optimization_parameter_fitting)
Fitting parameters to data involves solving an optimisation problem (that is, finding the parameter set that optimally fits your model to your data, typically by minimising a cost function). The SciML ecosystem's primary package for solving optimisation problems is [Optimization.jl](https://github.com/SciML/Optimization.jl). It provides access to a variety of solvers via a single common interface by wrapping a large number of optimisation libraries that have been implemented in Julia.

This tutorial demonstrates both how to create parameter fitting cost functions using the [DiffEqParamEstim.jl](https://github.com/SciML/DiffEqParamEstim.jl) package, and how to use Optimization.jl to minimise these. Optimization.jl can also be used in other contexts, such as finding parameter sets that maximise the magnitude of some system behavior. More details on how to use these packages can be found in their [respective](https://docs.sciml.ai/Optimization/stable/) [documentations](https://docs.sciml.ai/DiffEqParamEstim/stable/).

## Basic example

Let us consider a simple catalysis network, where an enzyme ($E$) turns a substrate ($S$) into a product ($P$):
```@example diffeq_param_estim_1 
using Catalyst
rn = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
```
From some known initial condition, and a true parameter set (which we later want to recover from the data) we generate synthetic data (on which we will demonstrate the fitting process).
```@example diffeq_param_estim_1 
# Define initial conditions and parameters.
u0 = [:S => 1.0, :E => 1.0, :SE => 0.0, :P => 0.0]
p_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5]

# Generate synthetic data.
using OrdinaryDiffEq
oprob_true = ODEProblem(rn, u0, (0.0, 10.0), p_true)
true_sol = solve(oprob_true, Tsit5())
data_sol = solve(oprob_true, Tsit5(); saveat=1.0)
data_ts = data_sol.t[2:end]
data_vals = (0.8 .+ 0.4*rand(10)) .* data_sol[:P][2:end]

# Plots the true solutions and the (synthetic data) measurements.
using Plots
plot(true_sol; idxs=:P, label="True solution", lw=8)
plot!(data_ts, data_vals; label="Measurements", seriestype=:scatter, ms=6, color=:blue)
```

Next, we will use DiffEqParamEstim.jl to build a loss function to measure how well our model's solutions fit the data.
```@example diffeq_param_estim_1 
using DiffEqParamEstim, Optimization
p_dummy = [:kB => 0.0, :kD => 0.0, :kP => 0.0]
oprob = ODEProblem(rn, u0, (0.0, 10.0), p_dummy)
loss_function = build_loss_objective(oprob, Tsit5(), L2Loss(data_ts, data_vals), Optimization.AutoForwardDiff(); maxiters=10000, verbose=false, save_idxs=4)
nothing # hide
```
To `build_loss_objective` we provide the following arguments:
- `oprob`: The `ODEProblem` with which we simulate our model (using some dummy parameter values, since we do not know these).
- `Tsit5()`: The numeric integrator we wish to simulate our model with.
- `L2Loss(data_ts, data_vals)`: Defines the loss function. While [other alternatives](https://docs.sciml.ai/DiffEqParamEstim/stable/getting_started/#Alternative-Cost-Functions-for-Increased-Robustness) are available, `L2Loss` is the simplest one (measuring the sum of squared distances between model simulations and data measurements). Its first argument is the time points at which the data is collected, and the second is the data's values.
- `Optimization.AutoForwardDiff()`: Our choice of [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) framework.

Furthermore, we can pass any number of additional optional arguments, these are then passed to the internal `solve()` function (which is used to solve our ODE). Here we provide the following additional arguments:
- `maxiters=10000`: If the ODE integrator takes a very large number of steps, that can be a sign of a very poor fit (or stiffness in the ODEs, but that is not a concern for our current example). Reducing the `maxiters` threshold reduces the time we waste on evaluating such models. 
- `verbose=false`: The simulation of models with highly unsuitable parameter sets typically generate various warnings (such as simulation terminations due to reaching the `maxiters` value). To avoid an overflow of such (here unnecessary) warnings as we evaluate a large number of parameter sets, we turn warnings off.
- `save_idxs=4`: The measured species (*P*) is the 4th species in our species vector (`species(rn)`). To ensure that the concentration of the right species is evaluated against the data, we set the numeric integrator to only save the value of this species.

Now we can create an `OptimizationProblem` using our `loss_function` and some initial guess of parameter values from which the optimiser will start:
```@example diffeq_param_estim_1 
optprob = OptimizationProblem(loss_function, [1.0, 1.0, 1.0])
nothing # hide
```

!!! note
    `OptimizationProblem` cannot currently accept parameter values in the form of a map (e.g. `[:kB => 1.0, :kD => 1.0, :kP => 1.0]`). These must be provided as individual values (using the same order as the parameters occur in in the `parameters(rs)` vector). Similarly, `build_loss_objective`'s `save_idxs` uses the species index, rather than the species directly. These inconsistencies should be remedied in future package releases.

Finally, we can optimise `optprob` to find the parameter set that best fits our data. Optimization.jl only provide a few optimisation methods natively. However, for each supported optimisation package, it provides a corresponding wrapper-package to import that optimisation package for use with Optimization.jl. E.g., if we wish to use [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)'s [Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) method, we must install and import the OptimizationOptimJL package. A summary of all, by Optimization.jl, supported optimisation packages can be found [here](https://docs.sciml.ai/Optimization/stable/#Overview-of-the-Optimizers). Here, we import the Optim.jl package and uses it to minimise our cost function (thus finding a parameter set that fits the data):
```@example diffeq_param_estim_1 
using OptimizationOptimJL
optsol = solve(optprob, Optim.NelderMead())
nothing # hide
```

We can now simulate our model for the corresponding parameter set, checking that it fits our data.
```@example diffeq_param_estim_1 
oprob_fitted = remake(oprob; p=optsol.u)
fitted_sol = solve(oprob_fitted, Tsit5())
plot!(fitted_sol; idxs=:P, label="Fitted solution", linestyle=:dash, lw=6, color=:lightblue)
```

!!! note
    Here, a good exercise is to check the resulting parameter set and note that, while it creates a good fit to the data, it does not actually correspond to the original parameter set. [Identifiability](https://www.sciencedirect.com/science/article/pii/S1364815218307278) is a concept that studies how to deal with this problem.

Say that we instead would like to use the [Broyden–Fletcher–Goldfarb–Shannon](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) algorithm, as implemented by the [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl) package. In this case we would run:
```@example diffeq_param_estim_1 
using OptimizationNLopt
sol = solve(optprob, NLopt.LD_LBFGS())
nothing # hide
```

## Optimisation problems with data for multiple species
Imagine that, in our previous example, we had measurements the concentration of both *S* and *P*:
```@example diffeq_param_estim_1 
data_vals_S = (0.8 .+ 0.4*rand(10)) .* data_sol[:S][2:end]
data_vals_P = (0.8 .+ 0.4*rand(10)) .* data_sol[:P][2:end]

plot(true_sol; idxs=[:S, :P], label=["True S" "True P"], lw=8)
plot!(data_ts, data_vals_S; label="Measured S", seriestype=:scatter, ms=6, color=:blue)
plot!(data_ts, data_vals_P; label="Measured P", seriestype=:scatter, ms=6, color=:red)
```

In this case we would have to use the `L2Loss(data_ts, hcat(data_vals_S, data_vals_P))` and `save_idxs=[1,4]` arguments in `loss_function`:
```@example diffeq_param_estim_1 
loss_function_S_P = build_loss_objective(oprob, Tsit5(), L2Loss(data_ts, Array(hcat(data_vals_S, data_vals_P)')), Optimization.AutoForwardDiff(); maxiters=10000, verbose=false, save_idxs=[1,4])
nothing # hide
```
Here, `Array(hcat(data_vals_S, data_vals_P)')` is required to put the data in the right form (in this case, a 2x10 matrix).

We can now fit our model to data and plot the results:
```@example diffeq_param_estim_1 
optprob_S_P = OptimizationProblem(loss_function_S_P, [1.0,1.0, 1.0])
optsol_S_P = solve(optprob_S_P, Optim.NelderMead())
oprob_fitted_S_P = remake(oprob; p=optsol_S_P.u)
fitted_sol_S_P = solve(oprob_fitted_S_P, Tsit5())
plot!(fitted_sol_S_P; idxs=[:S, :P], label="Fitted solution", linestyle=:dash, lw=6, color=[:lightblue :pink])
```

## Setting parameter constraints and boundaries
Sometimes, it is desired to set boundaries on parameter values. Indeed, this can speed up the optimisation process (by preventing searching through unfeasible parts of parameter space) or be a prerequisite of some optimisation methods. This can be done by passing the `lb` (lower bounds) and `up` (upper bounds) arguments to `OptimizationProblem`. These are vectors (of the same length as the number of parameters), with each argument corresponding to the boundary value of the parameter with the same index (as used in the `parameters(rn)` vector). If we wish to constrain each parameter to the interval *(0.1, 10.0)* this can be done through:
```@example diffeq_param_estim_1 
optprob = OptimizationProblem(loss_function, [1.0, 1.0, 1.0]; lb = [0.1, 0.1, 0.1], ub = [10.0, 10.0, 10.0])
nothing # hide
```

In addition to boundaries, Optimization.jl also supports setting [linear and non-linear constraints](https://docs.sciml.ai/Optimization/stable/tutorials/constraints/#constraints) on its output solution for some optimizers.

## Parameter fitting with known parameters
If from previous knowledge we know that $kD = 0.1$, and only want to fit the values of $kD$ and $kP$, this can be achieved through `build_loss_objective`'s `prob_generator` argument. First, we create a function (`fixed_p_prob_generator`) that modifies our `ODEProblem` to incorporate this knowledge:
```@example diffeq_param_estim_1 
fixed_p_prob_generator(prob, p) = remake(prob; p = vcat(p[1], 0.1, p[2]))
nothing # hide
```
Here, it takes the `ODEProblem` (`prob`) we simulate, and the parameter set used, during the optimisation process (`p`), and creates a modified `ODEProblem` (by setting a customised parameter vector [using `remake`](@ref simulation_structure_interfacing_remake)). Now we create our modified loss function:
```@example diffeq_param_estim_1 
loss_function_fixed_kD = build_loss_objective(oprob, Tsit5(), L2Loss(data_ts, data_vals), Optimization.AutoForwardDiff(); prob_generator = fixed_p_prob_generator, maxiters=10000, verbose=false, save_idxs=4)
nothing # hide
```

We can create an optimisation problem from this one like previously, but keep in mind that it (and its output results) only contains two parameter values (*kB* and *kP):
```@example diffeq_param_estim_1 
optprob_fixed_kD = OptimizationProblem(loss_function_fixed_kD, [1.0, 1.0])
optsol_fixed_kD = solve(optprob_fixed_kD, Optim.NelderMead())
nothing # hide
```

## Fitting parameters on the logarithmic scale
Often it can be advantageous to fit parameters on a [logarithmic, rather than linear, scale](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008646). The best way to proceed is to simply replace each parameters in the model definition by its logarithmic version:
```@example diffeq_param_estim_2
using Catalyst
rn = @reaction_network begin
    10^kB, S + E --> SE
    10^kD, SE --> S + E
    10^kP, SE --> P + E
end
```
And then proceeding, by keeping in mind that parameter values are logarithmic. Here, setting
```@example diffeq_param_estim_2
p_true = [:kB => 0.0, :kD => -1.0, :kP => 10^(0.5)]
nothing # hide
```
corresponds to the same true parameter values as used previously (`[:kB => 1.0, :kD => 0.1, :kP => 0.5]`).

## Parameter fitting to multiple experiments
Say that we had measured our model for several different initial conditions, and would like to fit our model to all these measurements simultaneously. This can be done by first creating a [corresponding `EnsembleProblem`](@ref advanced_simulations_ensemble_problems). How to then create loss functions for these are described in more detail [here](https://docs.sciml.ai/DiffEqParamEstim/stable/tutorials/ensemble/).

## Optimisation solver options
Optimization.jl supports various [optimisation solver options](https://docs.sciml.ai/Optimization/stable/API/solve/) that can be supplied to the `solve` command. For example, to set a maximum number of seconds (after which the optimisation process is terminated), you can use the `maxtime` argument:
```@example diffeq_param_estim_1 
optsol_fixed_kD = solve(optprob, Optim.NelderMead(); maxtime=100)
nothing # hide
```

---
## References
[^1]: [Alejandro F. Villaverde, Dilan Pathirana, Fabian Fröhlich, Jan Hasenauer, Julio R. Banga, *A protocol for dynamic model calibration*, Briefings in Bioinformatics (2023).](https://academic.oup.com/bib/article/23/1/bbab387/6383562?login=false)