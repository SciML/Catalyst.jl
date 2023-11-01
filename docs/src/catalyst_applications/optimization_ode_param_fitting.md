# [Parameter Fitting for ODEs using SciML/Optimization.jl and DiffEqParamEstim.jl](@id optimization_parameter_fitting)
Fitting parameters to data involves solving a optimisation problem (that is, finding the parameter set that optimally fits you model to your data, typically by minimising some cost function). The SciML ecosystem's primary package for solving optimisation problems is [Optimization.jl](https://github.com/SciML/Optimization.jl). While it does not implement any optimisation methods itself, it wraps the large number of optimisation methods that have been implemented in Julia in a single common interface. 

This tutorial both demonstrate how to use first create a parameter fitting cost function using the [DiffEqParamEstim.jl](https://github.com/SciML/DiffEqParamEstim.jl) package, and how to use Optimization.jl to optimise this (and potentially other) functions. More details on how to use these packages can be found in their [respective](https://docs.sciml.ai/Optimization/stable/) [documentations](https://docs.sciml.ai/DiffEqParamEstim/stable/).

## Basic example

Let us consider a simple catalysis network, where an enzyme (*E*) turns a substrate (*S*) into a product (*P*):
```@example diffeq_param_estim_1 
using Catalyst, PEtab
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
plot!(data_ts, data_vals; label="Measurements", seriestype=:scatter, ms=6)
```

Next, we will use DiffEqParamEstim to build a loss function, describing how well simulations of our model fit the data.
```@example diffeq_param_estim_1 
using DiffEqParamEstim, Optimization
p_dummy = [:kB => 0.0, :kD => 0.0, :kP => 0.0]
oprob = ODEProblem(rn, u0, (0.0, 10.0), p_dummy)
loss_function = build_loss_objective(oprob, Tsit5(), L2Loss(data_ts, data_vals), Optimization.AutoForwardDiff(); maxiters=10000, verbose=false, save_idxs=4)
nothing #
```
To `build_loss_objective` we provide the following arguments:
- `oprob`: The `ODEProblem` with which we simulate our model (using some dummy parameter values, since we do not know these).
- `Tsit5()`: The numeric integrator we wish to simulate our model with.
- `L2Loss(data_ts, data_vals)`: Defines the loss function. While [other alternatives](https://docs.sciml.ai/DiffEqParamEstim/stable/getting_started/#Alternative-Cost-Functions-for-Increased-Robustness) are available, `L2Loss` is the simplest one. Its first argument is the time points at which the data is collected, and the second is the data's values.
- `Optimization.AutoForwardDiff()`: Our choice of [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) framework.

Furthermore, we can pass any number of additional optional arguments, these are then passed to the internal `solve()` function (that are used to simulate our ODE). Here we provide the following additional arguments:
- `maxiters=10000`: If the ODE integrator takes a very large number of steps, that s typically a sign of a very poor fit. Reducing the `maxiters` threshold reduces the time we waste on evaluating such models. 
- `verbose=false`: The simulation of models with very unsuitable parameter sets typically generate warnings. To avoid an overflow of such (here unnecessary) warnings as we evaluate a large number of parameter sets, we turn warning off.
- `save_idxs=4`: The produce (*P*) is the 4th species in our species vector. To ensure that the concentration of the right species is evaluated against the data, we set the numeric integrator to only save teh value of this species.

Now we can create an `OptimizationProblem` using our `loss_function` and some initial guess of parameter values from which the optimizer will start:
```@example diffeq_param_estim_1 
optprob = OptimizationProblem(loss_function, [1.0, 1.0, 1.0])
nothing #
```

!!! note
    `OptimizationProblem` cannot currently accept parameter values in the form of a map. These must be provided as individual values (using the same order as the parameters occur in in the `parameters(rs)` vector).

Finally, we can use optimise `optprob` to find the parameter set that best fits our data. Optimization.jl does not provide any optimisation methods by default. Instead, for each supported optimisation package, it provides a corresponding wrapper-package to import that optimisation package for using with Optimization. E.g., if we wish to [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)'s [Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) method, we most install and import the OptimizationOptimJL package. A summary of all, by Optimization.jl, supported optimisation packages can be found [here](https://docs.sciml.ai/Optimization/stable/#Overview-of-the-Optimizers). Here, we import the Optim.jl package and uses it to minimise our cost function (finding a parameter set that fits the data):
```@example diffeq_param_estim_1 
using OptimizationOptimJL
sol = solve(optprob, NelderMead())
```

We can now simulate our model for the corresponding parameter set, checking that it fits our data.
```@example diffeq_param_estim_1 
oprob_fitted = remake(oprob; p=sol.u)
fitted_sol = solve(oprob, Tsit5())
plot!(fitted_sol; idxs=:P, label="Fitted solution", linestyle=:dash, lw=6, color=:purple)
```

!!! note
Here, a good exercise is to check the resulting parameter set and note that, while it creates a good fit to the data, it does not actually correspond to the original parameter set. [Identifiability](https://www.sciencedirect.com/science/article/pii/S1364815218307278) is a concept that studies how to deal with this problem.

Say that we instead would like to use the [CMA Evolution Strategy](https://en.wikipedia.org/wiki/CMA-ES) method, as implemented by the [CMAEvolutionStrategy.jl](https://github.com/jbrea/CMAEvolutionStrategy.jl) package, we would instead run
```@example diffeq_param_estim_1 
using OptimizationCMAEvolutionStrategy
sol = solve(optprob, CMAEvolutionStrategyOpt())
```

## Optimisation problems with data for multiple species

## Setting parameter constraints and boundaries