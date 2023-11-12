# Optimization for non-data fitting purposes
In previous tutorials we have described how to use [PEtab.jl](@ref petab_parameter_fitting) and [Optimization.jl](@ref optimization_parameter_fitting) for parameter fitting. This involves solving an optimisation problem (to find the parameter set best fitting the data). There are, however, other situations that requires solving optimisation problems. Typically, these involves the creation of a custom cost function, which optimum can then be found using Optimization.jl. In this tutorial we will describe this process, demonstrating how parameter space can be searched to find values that achieves a desired system behaviour.

## Maximising the pulse amplitude of a incoherent feed forward loops.
Incoherent feedforward loops (network motifs where a single components both activates and deactivates a downstream component) are able to generate pules in response to step inputs [^1]. In this tutorial we will consider such a incoherent feedforward loop, attempting to generate a system with as prominent response pulse as possible. Our model consists of 3 species: $X$ (the input node), $Y$ (an intermediary), and $Z$ (the output node). In it, $X$ activates the production of both $Y$ and $Z$, with $Y$ also deactivating $Z$. When $X$ is activated, their will be a brief time window where $Y$ is still inactive, and $Z$ is activated. However, as $Y$ becomes active, it will turn $Z$ off, creating a pulse.
```@example behaviour_optimization
using Catalyst
incoherent_feed_forward = @reaction_network begin
    pX, 0 --> X
    pY*X, 0 --> Y
    pZ*X/Y, 0 --> Z
    1.0, (X,Y,Z) --> 0
end
```
To demonstrate this pulsing behaviour we will simulate the system for an example parameter set. First, we will create the `ODEProblem` corresponding to our simulation. Here, we select an initial condition (`u0`) so the system begins in a steady state.
```@example behaviour_optimization
using DifferentialEquations
@unpack X, Y, Z, pX, pY, pZ = incoherent_feed_forward
example_p = [pX => 0.1, pY => 1.0, pZ => 1.0]
example_u0 = [X => 0.1, Y => 0.1, Z => 1.0]
ode_prob = ODEProblem(incoherent_feed_forward, example_u0, (0.0,50.0), example_p)
```
To trigger he activation of $X$ we will use a [callback](@ref advanced_simulations_callbacks), increasing its production rate ($pX$) by a factor of $10$ at the time $t=10.0$. We supply the callback to our simulation, and plot the result
```@example behaviour_optimization
using Plots
activation_cb = PresetTimeCallback([10.0], int -> int[pX] *= 10.0)
sol = solve(ode_prob, Tsit5(); callback=activation_cb)
plot(sol; lw=4)
```
Here we note that, while $X$ and $Y$ reach new steady state levels in response to the increase in $pX$, $Z$ resumes to its initial concentration after the pulse.

We will now attempt to find the parameter set $(pX,pY,pZ)$ which maximises the response pulse amplitude (defined by the maximum activity of $Z$ subtracted by its steady state activity). To do this, we create a custom cost function:
```@example behaviour_optimization
function pulse_amplitude(p, _)
    u0_new = [X => p[1], Y => p[1]*p[2], Z => p[3]/p[2]^2]
    ode_prob_local = remake(ode_prob; u0=u0_new, p=p)
    sol = solve(ode_prob_local, Tsit5(); callback=activation_cb, verbose=false, maxiters=10000)
    SciMLBase.successful_retcode(sol) || return Inf
    return -(maximum(sol[:Z])-sol[:Z][1])
end
nothing # here
```
The cost function takes two arguments (a parameter value `p`, and an additional one which we will ignore). It first calculates the new initial steady concentration (for the given parameter set), and then creates an updated `ODEProblem` using it as initial conditions (and the new parameter set as parameters). While we could create a new `ODEProblem` within the cost function, these functions are often called a large number of times during the optimization process (making performance important). Here, using [`remake` on a previously created `ODEProblem`](@id simulation_structure_interfacing_remake) is more performant than creating an new one. Next, we simulate our, remembering to use the activation callback. Just like when using Optimization.jl to fit parameters, we use the `verbose=false` options to prevent unnecessary print outs, and a reduced `maxiters` value to reduce time spent simulating (to the model) unsuitable parameter sets. We also use `SciMLBase.successful_retcode(sol)` to check whether the simulation return code indicate a successful simulation (and if it did not, returns a large cost function value). Finally, Optimization.jl finds the function's *minimum value*, so to find the *maximum* relative pulse amplitude, we make our cost function return the negation of that value.

Just like for [parameter fitting](@ref optimization_parameter_fitting), we create a `OptimizationProblem` using our cost function, and some initial guess of the parameter value. We also set upper and lower bounds for each parameter using the `lb` and `ub` optional arguments (in this case limiting each parameter's value to the interval $(0.1,10.0)$).
```@example behaviour_optimization
using Optimization
initial_guess = [1.0, 1.0, 1.0]
opt_prob = OptimizationProblem(cost_function, initial_guess; lb = [1e-1, 1e-1, 1e-1], ub = [1e1, 1e1, 1e1])
```
!!! note
    As described in a previous section on Optimization.jl, `OptimizationProblem`s do not support setting parameter values using maps. We must instead set `initial_guess`'s values using a vector. Here, the i'th value correspond to the value of the i'th parameter in the `parameters(incoherent_feed_forward)` vector.

Like described here, Optimization.jl supports a wide range of optimisation algorithms. Here we use one from BlackBoxOptim.jl
```@example behaviour_optimization
using OptimizationBBO
opt_sol = solve(opt_prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
```
Finally, we plot a simulation using the found parameter set (stored in `opt_sol.u`):
```@example behaviour_optimization
u0_result = [X => opt_sol.u[1], Y => opt_sol.u[1]*opt_sol.u[2], Z => opt_sol.u[3]/opt_sol.u[2]]
oprob_result = remake(ode_prob; u0=u0_result, p=opt_sol.u)
sol_result = solve(oprob_result, Tsit5(); callback=activation_cb)
plot(sol_result; lw=4, idxs=Z)
```
For this model, it actually turns out that $Z$'s maximum pulse amplitude is equal to twice its steady state concentration. Hence, the maximisation of its pulse amplitude is equivalent to maximising its steady state concentration.

!!! note
    Especially if you check Optimization.jl's documentation, you will note that cost functions have teh `f(u,p)` form. This is because `OptimizationProblem`s (like e.g. `ODEProblem`s) can take both variables (which can be varied in the optimisation problem), but also parameter that are fixed. In our case, the *optimisation variables* correspond to our *model parameters*. Hence, our model parameter values is the `u` input. This is also why we find the optimisation solution (our optimised parameter set) in `opt_sol`'s `u` field. Our optimisation problem does not actually have any parameters, hence, the second argument of `pulse_amplitude` is unused (that is why we call it `_`, a name commonly indicating unused function arguments). There are several modifications to our problem where it would actually have parameters. E.g. our model might have had additional parameters (e.g. a degradation rate) which we would like to keep fixed throughout the optimisation process. If we then would like to run the optimisation process for several different values of these fixed parameters, we could have made them parameters to our `OptimizationProblem` (and their values provided as a third argument, after `initial_guess`).


---
## References
[^1]: [Goentoro, L et al. *The incoherent feedforward loop can provide fold-change detection in gene regulation*, Molecular Cell (2009).](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2896310/)