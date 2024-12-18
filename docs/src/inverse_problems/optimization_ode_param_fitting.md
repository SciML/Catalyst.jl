# [Parameter Fitting for ODEs using Optimization.jl](@id optimization_parameter_fitting)
Fitting parameters to data involves solving an optimisation problem (that is, finding the parameter set that optimally fits your model to your data, typically by minimising an objective function)[^1]. The SciML ecosystem's primary package for solving optimisation problems is [Optimization.jl](https://github.com/SciML/Optimization.jl). It provides access to a variety of solvers via a single common interface by wrapping a large number of optimisation libraries that have been implemented in Julia.

This tutorial demonstrates how to 
1. Create a custom objective function which minimiser corresponds to the parameter set optimally fitting the data.
2. Use Optimization.jl to minimize this objective function and find the parameter set providing the optimal fit.

For simple parameter fitting problems (such as the one outlined below), [PEtab.jl often provides a more straightforward parameter fitting interface](@ref petab_parameter_fitting). However, Optimization.jl provides additional flexibility in defining your objective function. Indeed, it can also be used in other contexts, such as [finding parameter sets that maximise the magnitude of some system behaviour](@ref behaviour_optimisation). More details on how to use Optimization.jl can be found in its [documentation](https://docs.sciml.ai/DiffEqOptimizationParamEstim/stable/).

## [Basic example](@id optimization_parameter_fitting_basics)

Let us consider a [Michaelis-Menten enzyme kinetics model](@ref basic_CRN_library_mm), where an enzyme ($E$) converts a substrate ($S$) into a product ($P$):
```@example optimization_paramfit_1 
using Catalyst
rn = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
```
From some known initial condition, and a true parameter set (which we later want to recover from the data) we generate synthetic data (on which we will demonstrate the fitting process).
```@example optimization_paramfit_1 
# Define initial conditions and parameters.
u0 = [:S => 1.0, :E => 1.0, :SE => 0.0, :P => 0.0]
ps_true = [:kB => 1.0, :kD => 0.1, :kP => 0.5]

# Generate synthetic data.
using OrdinaryDiffEqDefault
oprob_true = ODEProblem(rn, u0, (0.0, 10.0), ps_true)
true_sol = solve(oprob_true)
data_sol = solve(oprob_true; saveat = 1.0)
data_ts = data_sol.t[2:end]
data_vals = (0.8 .+ 0.4*rand(10)) .* data_sol[:P][2:end]

# Plots the true solutions and the (synthetic) data measurements.
using Plots
plot(true_sol; idxs = :P, label = "True solution", lw = 8)
plot!(data_ts, data_vals; label = "Measurements", seriestype=:scatter, ms = 6, color = :blue)
plt = plot(true_sol; idxs = :P, label = "True solution", lw = 8) # hide
plot!(plt, data_ts, data_vals; label = "Measurements", seriestype=:scatter, ms = 6, color = :blue) # hide
Catalyst.PNG(plot(plt; fmt = :png, dpi = 200)) # hide
```

Next, we will formulate an objective function which, for a single parameter set, simulates our model and computes the sum-of-square distance between the data and the simulation (non-sum-of-square approaches can be used, but this is the most straightforward one).
```@example optimization_paramfit_1 
ps_init = [:kB => 1.0, :kD => 1.0, :kP => 1.0]
oprob_base = ODEProblem(rn, u0, (0.0, 10.0), ps_init)
function objective_function(p, _)
    p = Pair.([:kB, :kD, :kP], p)
    oprob = remake(oprob_base; p)
    sol = solve(oprob; saveat = data_ts, save_idxs = :P, verbose = false, maxiters = 10000)
    SciMLBase.successful_retcode(sol) || return Inf
    return sum((sol .- data_vals) .^2)
end
```
When our optimisation algorithm searches parameter space it will likely consider many highly non-plausible parameter sets. To better handle this we:
1. Add `maxiters = 10000` to our `solve` command. As most well-behaved ODEs can be solved in relatively few timesteps, this speeds up the optimisation procedure by preventing us from spending too much time trying to simulate (for the model) unsuitable parameter sets.
2. Add `verbose = false` to our `solve` command. This prevents (potentially a very large number of) warnings from being printed to our output as unsuitable parameter sets are simulated.
3. Add the line `SciMLBase.successful_retcode(sol) || return Inf`, which returns an infinite value for parameter sets which does not lead to successful simulations.

To improve optimisation performance, rather than creating a new `ODEProblem` in each iteration, we pre-declare one which we [apply `remake` to](@ref simulation_structure_interfacing_problems_remake). We also use the `saveat = data_ts, save_idxs = :P` arguments to only save the values of the measured species at the measured time points.

We can now create an `OptimizationProblem` using our `objective_function` and some initial guess of parameter values from which the optimiser will start:
```@example optimization_paramfit_1 
p_guess = [1.0, 1.0, 1.0]
optprob = OptimizationProblem(objective_function, p_guess)
nothing # hide
```
!!! note
    `OptimizationProblem`s cannot currently accept parameter values in the form of a map (e.g. `[:kB => 1.0, :kD => 1.0, :kP => 1.0]`). These must be provided as individual values (using the same order as the parameters occur in in the `parameters(rs)` vector). This should hopefully be remedied in future Optimization releases.

!!! note
    Especially if you check Optimization.jl's documentation, you will note that objective functions have the `f(u,p)` form. This is because `OptimizationProblem`s (like e.g. `ODEProblem`s) can take both variables (which are varied during the optimisation procedure), but also parameters that are fixed. In our case, the *optimisation variables* correspond to our *model parameters*. Hence, our model parameter values (`p`) are the first argument (`u`). This is also why we find the optimisation solution (our optimised parameter set) in `opt_sol`'s `u` field. Our optimisation problem does not actually have any parameters, hence, the second argument of `objective_function` is unused (that is why we call it `_`, a name commonly indicating unused function arguments).

    There are several modifications to our problem where it would actually have parameters. E.g. we might want to run the optimisation where one parameter has a known fixed value. If we then would like to rerun this for alternative fixed values, this value could be encoded as an `OptimizationProblem` parameter.

Finally, we can solve `optprob` to find the parameter set that best fits our data. Optimization.jl only provides a few optimisation methods natively. However, for each supported optimisation package, it provides a corresponding wrapper package to import that optimisation package for use with Optimization.jl. E.g., if we wish to use [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl)'s [Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) method, we must install and import the OptimizationNLopt package. A summary of all, by Optimization.jl supported, optimisation packages can be found [here](https://docs.sciml.ai/Optimization/stable/#Overview-of-the-Optimizers). Here, we import the NLopt.jl package and uses it to minimise our objective function (thus finding a parameter set that fits the data):
```@example optimization_paramfit_1 
using OptimizationNLopt
optsol = solve(optprob, NLopt.LN_NELDERMEAD())
nothing # hide
```

We can now simulate our model for the found parameter set (stored in `optsol.u`), checking that it fits our data.
```@example optimization_paramfit_1 
oprob_fitted = remake(oprob_base; p = Pair.([:kB, :kD, :kP], optsol.u))
fitted_sol = solve(oprob_fitted)
plot!(fitted_sol; idxs = :P, label = "Fitted solution", linestyle = :dash, lw = 6, color = :lightblue)
plot!(plt, fitted_sol; idxs = :P, label = "Fitted solution", linestyle = :dash, lw = 6, color = :lightblue) # hide
Catalyst.PNG(plot(plt; fmt = :png, dpi = 200)) # hide
```

!!! note
    Here, a good exercise is to check the resulting parameter set and note that, while it creates a good fit to the data, it does not actually correspond to the original parameter set. Identifiability is a concept that studies how to deal with this problem.<!--NTS: re-add ref when identifiablity works again-->

Say that we instead would like to use a [genetic algorithm](https://en.wikipedia.org/wiki/Differential_evolution) approach, as implemented by the [Evolutionary.jl](https://github.com/wildart/Evolutionary.jl) package. In this case we can run:
```@example optimization_paramfit_1 
using OptimizationEvolutionary
sol = solve(optprob, Evolutionary.GA())
nothing # hide
```
to solve `optprob` for this combination of solve and implementation.

## [Utilising automatic differentiation](@id optimization_parameter_fitting_AD)
Optimisation methods can be divided into differentiation-free and differentiation-based optimisation methods. E.g. consider finding the minimum of the function $f(x) = x^2$, given some initial guess of $x$. Here, we can simply compute the differential and descend along it until we find $x=0$ (admittedly, for this simple problem the minimum can be computed directly). This principle forms the basis of optimisation methods such as [gradient descent](https://en.wikipedia.org/wiki/Gradient_descent), which utilises information of a function's differential to minimise it. When attempting to find a global minimum, to avoid getting stuck in local minimums, these methods are often augmented by additional routines. While the differentiation of most algebraic functions is trivial, it turns out that even complicated functions (such as the one we used above) can be differentiated computationally through the use of [*automatic differentiation* (AD)](https://en.wikipedia.org/wiki/Automatic_differentiation).

Through packages such as [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), [ReverseDiff.jl](https://github.com/JuliaDiff/ReverseDiff.jl), and [Zygote.jl](https://github.com/FluxML/Zygote.jl), Julia supports AD for most code. Specifically for code including simulation of differential equations, differentiation is supported by [SciMLSensitivity.jl](https://github.com/SciML/SciMLSensitivity.jl). Generally, AD can be used without specific knowledge from the user, however, it requires an additional step in the construction of our `OptimizationProblem`. Here, we create a [specialised `OptimizationFunction` from our objective function](https://docs.sciml.ai/Optimization/stable/API/optimization_function/#optfunction). To it, we will also provide our choice of AD method. There are [several alternatives](https://docs.sciml.ai/Optimization/stable/API/optimization_function/#Automatic-Differentiation-Construction-Choice-Recommendations), and in our case we will use `AutoForwardDiff()` (a good choice for small optimisation problems). We can then create a new `OptimizationProblem` using our updated objective function:
```@example behaviour_optimization
opt_func = OptimizationFunction(objective_function, AutoForwardDiff())
opt_prob = OptimizationProblem(opt_func, p_guess)
nothing # hide
```
Finally, we can find the optimum using some differentiation-based optimisation methods. Here we will use [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)'s [Broyden–Fletcher–Goldfarb–Shanno algorithm](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) implementation:
```@example behaviour_optimization
using OptimizationOptimJL
opt_sol = solve(opt_prob, OptimizationOptimJL.BFGS())
```

## [Optimisation problems with data for multiple species](@id optimization_parameter_fitting_multiple_species)
Imagine that, in our previous example, we had measurements of the concentration of both $S$ and $P$:
```@example optimization_paramfit_1 
data_vals_S = (0.8 .+ 0.4*rand(10)) .* data_sol[:S][2:end]
data_vals_P = (0.8 .+ 0.4*rand(10)) .* data_sol[:P][2:end]

plot(true_sol; idxs = [:S, :P], label = ["True S" "True P"], lw = 8)
plot!(data_ts, data_vals_S; label = "Measured S", seriestype=:scatter, ms = 6, color = :blue)
plot!(data_ts, data_vals_P; label = "Measured P", seriestype=:scatter, ms = 6, color = :red)
plt2 = plot(true_sol; idxs = [:S, :P], label = ["True S" "True P"], lw = 8)  # hide
plot!(plt2, data_ts, data_vals_S; label = "Measured S", seriestype = :scatter, ms = 6, color = :blue)  # hide
plot!(plt2, data_ts, data_vals_P; label = "Measured P", seriestype = :scatter, ms = 6, color = :red)  # hide
Catalyst.PNG(plot(plt2; fmt = :png, dpi = 200)) # hide
```

In this case we simply modify our objective function to take this into account:
```@example optimization_paramfit_1 
function objective_function_S_P(p, _)
    p = Pair.([:kB, :kD, :kP], p)
    oprob = remake(oprob_base; p)
    sol = solve(oprob; saveat = data_ts, save_idxs = [:S, :P], verbose = false, maxiters = 10000)
    SciMLBase.successful_retcode(sol) || return Inf
    return sum((sol[:S] .- data_vals_S) .^2 + (sol[:P] .- data_vals_P) .^2)
end
```
Here we do not normalise the contribution from each species to the objective function. However, if species are present at different concentration levels this might be necessary (or you might essentially only take the highest concentration species(s) into account).

We can now fit our model to data and plot the results:
```@example optimization_paramfit_1 
optprob_S_P = OptimizationProblem(objective_function_S_P, p_guess)
optsol_S_P = solve(optprob_S_P, Optim.NelderMead())
oprob_fitted_S_P = remake(oprob_base; p = optsol_S_P.u)
fitted_sol_S_P = solve(oprob_fitted_S_P)
plot!(fitted_sol_S_P; idxs=[:S, :P], label="Fitted solution", linestyle = :dash, lw = 6, color = [:lightblue :pink])
plot!(plt2, fitted_sol_S_P; idxs=[:S, :P], label="Fitted solution", linestyle = :dash, lw = 6, color = [:lightblue :pink]) # hide
Catalyst.PNG(plot(plt2; fmt = :png, dpi = 200)) # hide
```

## [Setting parameter constraints and boundaries](@id optimization_parameter_fitting_constraints)
Sometimes, it is desirable to set boundaries on parameter values. Indeed, this can speed up the optimisation process (by preventing searching through unfeasible parts of parameter space), and can also be a requirement for some optimisation methods. This can be done by passing the `lb` (lower bounds) and `up` (upper bounds) arguments to `OptimizationProblem`. These are vectors (of the same length as the number of parameters), with each argument corresponding to the boundary value of the parameter with the same index. If we wish to constrain each parameter to the interval $(0.1, 10.0)$ this can be done through:
```@example optimization_paramfit_1 
optprob = OptimizationProblem(objective_function, [1.0, 1.0, 1.0]; lb = [1e-1, 1e-1, 1e-1], ub = [1e1, 1e1, 1e1])
nothing # hide
```

In addition to boundaries, Optimization.jl also supports setting [linear and non-linear constraints](https://docs.sciml.ai/Optimization/stable/tutorials/constraints/#constraints) on its output solution (only available for some optimisers).

## [Parameter fitting with known parameters](@id optimization_parameter_fitting_known_parameters)
If we from previous knowledge know that $kD = 0.1$, and only want to fit the values of $kB$ and $kP$, this can be achieved by making corresponding changes to our objective function.
```@example optimization_paramfit_1 
function objective_function_known_kD(p, _)
    p = Pair.([:kB, :kD, :kP], [p[1], 0.1, p[2]])
    oprob = remake(oprob_base; p)
    sol = solve(oprob; saveat = data_ts, save_idxs = :P, verbose = false, maxiters = 10000)
    SciMLBase.successful_retcode(sol) || return Inf
    return sum((sol .- data_vals) .^2)
end
```
We can now create and solve the corresponding `OptimizationProblem`, but with only two parameters in the initial guess.
```@example optimization_paramfit_1 
optprob_known_kD = OptimizationProblem(objective_function_known_kD, [1.0, 1.0])
optsol_known_kD = solve(optprob_known_kD, Optim.NelderMead())
nothing # hide
```

## [Optimisation solver options](@id optimization_parameter_fitting_solver_options)
Optimization.jl supports various [optimisation solver options](https://docs.sciml.ai/Optimization/stable/API/solve/) that can be supplied to the `solve` command. For example, to set a maximum number of seconds (after which the optimisation process is terminated), you can use the `maxtime` argument:
```@example optimization_paramfit_1 
optsol_fixed_kD = solve(optprob, Optim.NelderMead(); maxtime = 100)
nothing # hide
```
It should be noted that not all solver options are available to all optimisation solvers.

## [Fitting parameters on the logarithmic scale](@id optimization_parameter_fitting_log_scale)
Often it can be advantageous to fit parameters on a [logarithmic scale (rather than on a linear scale)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008646). The most straightforward way to do this is to simply replace each parameter in the model definition by its logarithmic version:
```@example optimization_paramfit_2
using Catalyst
rn = @reaction_network begin
    10^kB, S + E --> SE
    10^kD, SE --> S + E
    10^kP, SE --> P + E
end
```
And then going forward, by keeping in mind that parameter values are logarithmic. Here, setting
```@example optimization_paramfit_2
p_true = [:kB => 0.0, :kD => -1.0, :kP => 10^(0.5)]
nothing # hide
```
corresponds to the same true parameter values as used previously (`[:kB => 1.0, :kD => 0.1, :kP => 0.5]`). Alternatively, we can provide the log-transform in the objective function:
```@example optimization_paramfit_1 
function objective_function_logtransformed(p, _)
    p = Pair.([:kB, :kD, :kP], 10.0 .^ p)
    oprob = remake(oprob_base; p)
    sol = solve(oprob; saveat = data_ts, save_idxs = :P, verbose = false, maxiters = 10000)
    SciMLBase.successful_retcode(sol) || return Inf
    return sum((sol .- data_vals) .^2)
end
```

---
## [Citation](@id optimization_parameter_fitting_citation)
If you use this functionality in your research, please cite the following paper to support the authors of the Optimization.jl package:
```
@software{vaibhav_kumar_dixit_2023_7738525,
	author = {Vaibhav Kumar Dixit and Christopher Rackauckas},
	month = mar,
	publisher = {Zenodo},
	title = {Optimization.jl: A Unified Optimization Package},
	version = {v3.12.1},
	doi = {10.5281/zenodo.7738525},
  	url = {https://doi.org/10.5281/zenodo.7738525},
	year = 2023
}
```

---
## References
[^1]: [Alejandro F. Villaverde, Dilan Pathirana, Fabian Fröhlich, Jan Hasenauer, Julio R. Banga, *A protocol for dynamic model calibration*, Briefings in Bioinformatics (2023).](https://academic.oup.com/bib/article/23/1/bbab387/6383562?login=false)
