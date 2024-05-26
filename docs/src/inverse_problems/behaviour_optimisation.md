# [Optimization for non-data fitting purposes](@id behaviour_optimisation)
In previous tutorials we have described how to use [PEtab.jl](@ref petab_parameter_fitting) and [Optimization.jl](@ref optimization_parameter_fitting) for parameter fitting. This involves solving an optimisation problem (to find the parameter set yielding the best model-to-data fit). There are, however, other situations that require solving optimisation problems. Typically, these involve the creation of a custom cost function, which optimum can then be found using Optimization.jl. In this tutorial we will describe this process, demonstrating how parameter space can be searched to find values that achieve a desired system behaviour. A more throughout description on how to solve these problems is provided by [Optimization.jl's documentation](https://docs.sciml.ai/Optimization/stable/) and the literature[^1]. 

## [Maximising the pulse amplitude of an incoherent feed forward loop](@id behaviour_optimisation_IFFL_example)
Incoherent feedforward loops (network motifs where a single component both activates and deactivates a downstream component) are able to generate pulses in response to step inputs[^2]. In this tutorial we will consider such an incoherent feedforward loop, attempting to generate a system with as prominent a response pulse as possible.

Our model consists of 3 species: $X$ (the input node), $Y$ (an intermediary), and $Z$ (the output node). In it, $X$ activates the production of both $Y$ and $Z$, with $Y$ also deactivating $Z$. When $X$ is activated, there will be a brief time window where $Y$ is still inactive, and $Z$ is activated. However, as $Y$ becomes active, it will turn $Z$ off. This creates a pulse of $Z$ activity. To trigger the system, we create [an event](@ref ref), which increases the production rate of $X$ ($pX$) by a factor of $10$ at time $t = 10$.
```@example behaviour_optimization
using Catalyst
incoherent_feed_forward = @reaction_network begin
    @discrete_events [10.0] => [pX ~ 10*pX]
    pX, 0 --> X
    pY*X, 0 --> Y
    pZ*X/Y, 0 --> Z
    1.0, (X,Y,Z) --> 0
end
```
To demonstrate this pulsing behaviour we will simulate the system for an example parameter set. We select an initial condition (`u0`) so the system begins in a steady state.
```@example behaviour_optimization
using OrdinaryDiffEq, Plots
example_p = [:pX => 0.1, :pY => 1.0, :pZ => 1.0]
tspan = (0.0, 50.0)
example_u0 = [:X => 0.1, :Y => 0.1, :Z => 1.0]

oprob = ODEProblem(incoherent_feed_forward, example_u0, tspan, example_p)
sol = solve(oprob, Tsit5())
plot(sol)
```
Here we note that, while $X$ and $Y$ reach new steady state levels in response to the increase in $pX$, $Z$ resumes to its initial concentration after the pulse.

We will now attempt to find the parameter set $(pX,pY,pZ)$ which maximises the response pulse amplitude (defined by the maximum activity of $Z$ subtracted by its steady state activity). To do this, we create a custom cost function:
```@example behaviour_optimization
function pulse_amplitude(p, _)
    ps = Dict([:pX => p[1], :pY => p[2], :pZ => p[2]])
    u0_new = [:X => ps[:pX], :Y => ps[:pX]*ps[:pY], :Z => ps[:pZ]/ps[:pY]^2]
    oprob_local = remake(oprob; u0=  u0_new, p = ps)
    sol = solve(oprob_local, verbose = false, maxiters = 10000)
    SciMLBase.successful_retcode(sol) || return Inf
    return -(maximum(sol[:Z]) - sol[:Z][1])
end
nothing # here
```
This cost function takes two arguments (a parameter value `p`, and an additional one which we will ignore here but discuss later). It first calculates the new initial steady state concentration for the given parameter set. Next, it creates an updated `ODEProblem` using the steady state as initial conditions and the, to the cost function provided, input parameter set. While we could create a new `ODEProblem` within the cost function, cost functions are often called a large number of times during the optimisation process (making performance important). Here, using [`remake` on a previously created `ODEProblem`](@ref simulation_structure_interfacing_problems_remake) is more performant than creating a new one. Just like [when using Optimization.jl to fit parameters to data](@ref optimization_parameter_fitting), we use the `verbose = false` option to prevent unnecessary simulation printouts, and a reduced `maxiters` value to reduce time spent simulating (for the model) unsuitable parameter sets. We also use `SciMLBase.successful_retcode(sol)` to check whether the simulation return code indicates a successful simulation (and if it did not, returns a large cost function value). Finally, Optimization.jl finds the function's *minimum value*, so to find the *maximum* relative pulse amplitude, we make our cost function return the negative pulse amplitude.

Just like for [parameter fitting](@ref optimization_parameter_fitting), we create a `OptimizationProblem` using our cost function, and some initial guess of the parameter value. We also set upper and lower bounds for each parameter using the `lb` and `ub` optional arguments (in this case limiting each parameter's value to the interval $(0.1,10.0)$).
```@example behaviour_optimization
using Optimization
initial_guess = [1.0, 1.0, 1.0]
opt_prob = OptimizationProblem(pulse_amplitude, initial_guess; lb = [1e-1, 1e-1, 1e-1], ub = [1e1, 1e1, 1e1])
```
!!! note
    As described in a [previous section on Optimization.jl](@ref optimization_parameter_fitting), `OptimizationProblem`s do not support setting parameter values using maps. We must instead set `initial_guess` values using a vector. Next, in the first line of our cost function, we reshape the parameter values to the common form used across Catalyst (e.g. `[:pX => p[1], :pY => p[2], :pZ => p[2]]`, however, here we use a dictionary to easier compute the steady state initial condition). We also note that the order used in this array corresponds to the order we give each parameter's bounds in `lb` and `ub`, and the order in which their values occur in the output solution.

As [previously described](@ref optimization_parameter_fitting), Optimization.jl supports a wide range of optimisation algorithms. Here we use one from [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl):
```@example behaviour_optimization
using OptimizationBBO
opt_sol = solve(opt_prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
```
Finally, we plot a simulation using the found parameter set (stored in `opt_sol.u`):
```@example behaviour_optimization
ps_res = Dict([:pX => opt_sol.u[1], :pY => opt_sol.u[2], :pZ => opt_sol.u[2]])
u0_res = [:X => ps_res[:pX], :Y => ps_res[:pX]*ps_res[:pY], :Z => ps_res[:pZ]/ps_res[:pY]^2]
oprob_res = remake(oprob; u0 = u0_res, p = ps_res)
sol_res = solve(oprob_res)
plot(sol_res; idxs=:Z)
```
For this model, it turns out that $Z$'s maximum pulse amplitude is equal to twice its steady state concentration. Hence, the maximisation of its pulse amplitude is equivalent to maximising its steady state concentration.

!!! note
    Especially if you check Optimization.jl's documentation, you will note that cost functions have the `f(u,p)` form. This is because `OptimizationProblem`s (like e.g. `ODEProblem`s) can take both variables (which can be varied in the optimisation problem), but also parameters that are fixed. In our case, the *optimisation variables* correspond to our *model parameters*. Hence, our model parameter values are the `u` input. This is also why we find the optimisation solution (our optimised parameter set) in `opt_sol`'s `u` field. Our optimisation problem does not actually have any parameters, hence, the second argument of `pulse_amplitude` is unused (that is why we call it `_`, a name commonly indicating unused function arguments). 
    
    There are several modifications to our problem where it would actually have parameters. E.g. our model might have had additional parameters (e.g. a degradation rate) which we would like to keep fixed throughout the optimisation process. If we then would like to run the optimisation process for several different values of these fixed parameters, we could have made them parameters to our `OptimizationProblem` (and their values provided as a third argument, after `initial_guess`).

## [Utilising automatic differentiation](@id behaviour_optimisation_AD)
Optimisation methods can be divided into differentiation-free and differentiation-based optimisation methods. E.g. consider finding the minimum of the function $f(x) = x^2$, given some initial guess of $x$. Here, we can simply compute the differential and descend along it until we find $x=0$ (admittedly, for this simple problem the minimum can be computed directly). This principle forms the basis of optimisation methods such as [gradient descent](https://en.wikipedia.org/wiki/Gradient_descent), which utilises information of a function's differential to minimise it. When attempting to find a global minimum, to avoid getting stuck in local minimums, these methods are often augmented by additional routines. While the differentiation of most algebraic functions is trivial, it turns out that even complicated functions (such as the one we used above) can be differentiated computationally through the use of [*automatic differentiation* (AD)](https://en.wikipedia.org/wiki/Automatic_differentiation).

Through packages such as [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), [ReverseDiff.jl](https://github.com/JuliaDiff/ReverseDiff.jl), and [Zygote.jl](https://github.com/FluxML/Zygote.jl), Julia supports AD for most code. Specifically for code including simulation of differential equations, differentiation is supported by [SciMLSensitivity.jl](https://github.com/SciML/SciMLSensitivity.jl). Generally, AD can be used without specific knowledge from the user, however, it requires an additional step in the construction of our `OptimizationProblem`. Here, we create a [specialised `OptimizationFunction` from our cost function](https://docs.sciml.ai/Optimization/stable/API/optimization_function/#optfunction). To it, we will also provide our choice of AD method. There are [several alternatives](https://docs.sciml.ai/Optimization/stable/API/optimization_function/#Automatic-Differentiation-Construction-Choice-Recommendations), and in our case we will use `AutoForwardDiff()` (a good choice for small optimisation problems). We can then create a new `OptimizationProblem` using our updated cost function:
```@example behaviour_optimization
opt_func = OptimizationFunction(pulse_amplitude, AutoForwardDiff())
opt_prob = OptimizationProblem(opt_func, initial_guess; lb = [1e-1, 1e-1, 1e-1], ub = [1e1, 1e1, 1e1])
``` 
Finally, we can find the optimum using some differentiation-based optimisation methods. Here we will use [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)'s `BFGS` method:
```@example behaviour_optimization
using OptimizationOptimJL
opt_sol = solve(opt_prob, OptimizationOptimJL.BFGS())
``` 

---
## [Citation](@id structural_identifiability_citation)
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
[^1]: [Mykel J. Kochenderfer, Tim A. Wheeler *Algorithms for Optimization*, The MIT Press (2019).](https://algorithmsbook.com/optimization/files/optimization.pdf)
[^2]: [Lea Goentoro, Oren Shoval, Marc W Kirschner, Uri Alon *The incoherent feedforward loop can provide fold-change detection in gene regulation*, Molecular Cell (2009).](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2896310/)