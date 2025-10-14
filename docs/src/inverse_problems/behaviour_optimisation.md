# [Optimization for Non-data Fitting Purposes](@id behaviour_optimisation)

In previous tutorials we have described how to use [PEtab.jl](https://github.com/sebapersson/PEtab.jl) and [Optimization.jl](@ref optimization_parameter_fitting) for parameter fitting. This involves solving an optimisation problem (to find the parameter set yielding the best model-to-data fit). There are, however, other situations that require solving optimisation problems. Typically, these involve the creation of a custom objective function, which minimizer can then be found using Optimization.jl. In this tutorial we will describe this process, demonstrating how parameter space can be searched to find values that achieve a desired system behaviour. Many options used here are described in more detail in [the tutorial on using Optimization.jl for parameter fitting](@ref optimization_parameter_fitting). A more throughout description of how to solve these problems is provided by [Optimization.jl's documentation](https://docs.sciml.ai/Optimization/stable/) and the literature[^1].

## [Maximising the pulse amplitude of an incoherent feed forward loop](@id behaviour_optimisation_IFFL_example)

Incoherent feedforward loops (network motifs where a single component both activates and deactivates a downstream component) are able to generate pulses in response to step inputs[^2]. In this tutorial we will consider such an incoherent feedforward loop, attempting to generate a system with as prominent a response pulse as possible.

Our model consists of 3 species: $X$ (the input node), $Y$ (an intermediary), and $Z$ (the output node). In it, $X$ activates the production of both $Y$ and $Z$, with $Y$ also deactivating $Z$. When $X$ is activated, there will be a brief time window where $Y$ is still inactive, and $Z$ is activated. However, as $Y$ becomes active, it will turn $Z$ off. This creates a pulse of $Z$ activity. To trigger the system, we create [an event](@ref constraint_equations_events), which increases the production rate of $X$ ($pX$) by a factor of $10$ at time $t = 10$.
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
using OrdinaryDiffEqDefault, Plots
example_p = [:pX => 0.1, :pY => 1.0, :pZ => 1.0]
tspan = (0.0, 50.0)
example_u0 = [:X => 0.1, :Y => 0.1, :Z => 1.0]

oprob = ODEProblem(incoherent_feed_forward, example_u0, tspan, example_p)
sol = solve(oprob)
plot(sol)
```
Here we note that, while $X$ and $Y$ reach new steady state levels in response to the increase in $pX$, $Z$ resumes to its initial concentration after the pulse.

We will now attempt to find the parameter set $(pX,pY,pZ)$ which maximises the response pulse amplitude (defined by the maximum activity of $Z$ subtracted by its steady state activity). To do this, we create a custom objective function:
```@example behaviour_optimization
function pulse_amplitude(p, _)
    p = Dict([:pX => p[1], :pY => p[2], :pZ => p[2]])
    u0 = [:X => p[:pX], :Y => p[:pX]*p[:pY], :Z => p[:pZ]/p[:pY]^2]
    oprob_local = remake(oprob; u0, p)
    sol = solve(oprob_local; verbose = false, maxiters = 10000)
    SciMLBase.successful_retcode(sol) || return Inf
    return -(maximum(sol[:Z]) - sol[:Z][1])
end
nothing # hide
```
This objective function takes two arguments (a parameter value `p`, and an additional one which we will ignore but is discussed in a note [here](@ref optimization_parameter_fitting_basics)). It first calculates the new initial steady state concentration for the given parameter set. Next, it creates an updated `ODEProblem` using the steady state as initial conditions and the, to the objective function provided, input parameter set.  Finally, Optimization.jl finds the function's *minimum value*, so to find the *maximum* relative pulse amplitude, we make our objective function return the negative pulse amplitude.

As described [in our tutorial on parameter fitting using Optimization.jl](@ref optimization_parameter_fitting_basics) we use `remake`, `verbose = false`, `maxiters = 10000`, and a check on the simulations return code, all providing various advantages to the optimisation procedure (as explained in that tutorial).

Just like for [parameter fitting](@ref optimization_parameter_fitting_basics), we create an `OptimizationProblem` using our objective function, and some initial guess of the parameter values. We also [set upper and lower bounds](@ref optimization_parameter_fitting_constraints) for each parameter using the `lb` and `ub` optional arguments (in this case limiting each parameter's value to the interval $(0.1,10.0)$).
```@example behaviour_optimization
using Optimization
initial_guess = [1.0, 1.0, 1.0]
opt_prob = OptimizationProblem(pulse_amplitude, initial_guess; lb = [1e-1, 1e-1, 1e-1], ub = [1e1, 1e1, 1e1])
nothing # hide
```
!!! note
    As described in a [previous section on Optimization.jl](@ref optimization_parameter_fitting), `OptimizationProblem`s do not support setting parameter values using maps. We must instead set `initial_guess` values using a vector. Next, in the first line of our objective function, we reshape the parameter values to the common form used across Catalyst (e.g. `[:pX => p[1], :pY => p[2], :pZ => p[2]]`, however, here we use a dictionary to easier compute the steady state initial condition). We also note that the order used in this array corresponds to the order we give each parameter's bounds in `lb` and `ub`, and the order in which their values occur in the output solution.

As [previously described](@ref optimization_parameter_fitting), Optimization.jl supports a wide range of optimisation algorithms. Here we use one from [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl):
```@example behaviour_optimization
using OptimizationBBO
opt_sol = solve(opt_prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
nothing # hide
```
Finally, we plot a simulation using the found parameter set (stored in `opt_sol.u`):
```@example behaviour_optimization
ps_res = Dict([:pX => opt_sol.u[1], :pY => opt_sol.u[2], :pZ => opt_sol.u[2]])
u0_res = [:X => ps_res[:pX], :Y => ps_res[:pX]*ps_res[:pY], :Z => ps_res[:pZ]/ps_res[:pY]^2]
oprob_res = remake(oprob; u0 = u0_res, p = ps_res)
sol_res = solve(oprob_res)
plot(sol_res; idxs = :Z)
```
For this model, it turns out that $Z$'s maximum pulse amplitude is equal to twice its steady state concentration. Hence, the maximisation of its pulse amplitude is equivalent to maximising its steady state concentration.

## [Other optimisation options](@id behaviour_optimisation_options)

How to use Optimization.jl is discussed in more detail in [this tutorial](@ref optimization_parameter_fitting). This includes options such as using [automatic differentiation](@ref optimization_parameter_fitting_AD), [setting constraints](@ref optimization_parameter_fitting_constraints), and setting [optimisation solver options](@ref optimization_parameter_fitting_solver_options). Finally, it discusses the advantages of [carrying out the fitting in logarithmic space](@ref optimization_parameter_fitting_log_scale), something which can be advantageous for the problem described above as well.

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
