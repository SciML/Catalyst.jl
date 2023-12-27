# [Advice for performant ODE simulations](@id ode_simulation_performance)
We have previously described how to perform ODE simulations of *chemical reaction network* (CRN) models. These simulations are typically fast and require little additional consideration. However, when a model is simulated many times (e.g. as a part of solving an inverse problem), or is very large, simulation runtimes may become noticeable. Here we will give some advice on how to improve performance for these cases.

Generally, there are few good ways to, before a simulation, determine the best options. Hence, while we below provide several options, if you face an application for which reducing runtime is critical (e.g. if you need to simulate the same ODE many times), it might be required to manually trial these various options to see which yields the best performance ([BenchmarkTools.jl's](https://github.com/JuliaCI/BenchmarkTools.jl) `@btime` macro is useful for this purpose). It should be noted that most default options typically perform well, and it is primarily for large models where investigating these options is worthwhile. All ODE simulations of Catalyst models are performed using the OrdinaryDiffEq.jl package, [which documentation](https://docs.sciml.ai/DiffEqDocs/stable/) provides additional advice on performance.

Generally, this small checklist provides a quick guide for dealing with ODE performance:
1. If performance is not critical, use [the default solver choice](@ref ode_simulation_performance_solvers) and do not worry further about the issue.
2. Determine whether your problem is [non-stiff nor stiff](@ref ode_simulation_performance_stiffness), and use the `Tsit5` solver for the former case and `Rodas5P` for the latter.
3. If more performance would be useful, read about solver selection (both in [this tutorial](@ref ode_simulation_performance_solvers) and [OrdinaryDiffEq's documentation](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)) and then try a few different solvers to find one with good performance.
4. If you have a large ODE, try the [various options for efficient Jacobian computation](@ref ode_simulation_performance_jacobian) (noting that some are non-trivial to use, and should only be investigated if required).
5. If you plan to simulate your ODE many times, try [parallelise it on CPUs or GPUs](@ref investigating) (with preference for the former, which is easier to use).

## [Regarding stiff and non-stiff problems and solvers](@id ode_simulation_performance_stiffness)
Generally, ODE problems can be categorised into [*stiff ODEs* and *non-stiff ODEs*](https://en.wikipedia.org/wiki/Stiff_equation). This categorisation is important due to stiff ODEs requiring specialised solvers. A common cause of failure to simulate an ODE is the use of a non-stiff solver for a stiff problem. There is no exact way to determine whether a given ODE is stiff or not, however, systems with several different time scales (e.g. a CRN with both slow and fast reactions) typically generate stiff ODEs.

Here we simulate the (stiff) [brusselator](https://en.wikipedia.org/wiki/Brusselator) model using the `Tsit5` solver (which is designed for non-stiff ODEs):
```@example ode_simulation_performance_1
using Catalyst, OrdinaryDiffEq, Plots

brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end

u0 = [:X => 1.0, :Y => 0.0]
tspan = (0.0, 20.0)
p = [:A => 10.0, :B => 40.0]
oprob = ODEProblem(brusselator, u0, tspan, p)

sol1 = solve(oprob, Tsit5())
plot(sol1)
``` 
We note that we get a warning, indicating that an instability was detected (the typical indication of a non-stiff solver being used for a stiff ODE). Furthermore, the resulting plot ends at $t≈10$, meaning that the simulation was not completed. Indeed, we can confirm this by checking the *return code* of the solution object:
```@example ode_simulation_performance_1
sol1.retcode
```
Next, we instead try the `Rodas5P` solver (which is designed for stiff problems):
```@example ode_simulation_performance_1
sol2 = solve(oprob, Rodas5P())
plot(sol2)
``` 
This time the simulation was successfully completed, which can be confirmed by checking the return code:
```@example ode_simulation_performance_1
sol2.retcode
``` 

Generally, ODE solvers can be divided into [*explicit* and *implicit* solvers](https://en.wikipedia.org/wiki/Explicit_and_implicit_methods). Roughly, explicit solvers (which include `Tsit5`) are better for non-stiff problems, with implicit solvers (like `Rodas5P`) being required for stiff problems. While we could use implicit solvers for all problems (to guarantee successful simulations irrespective of stiffness), these are typically less performant on non-stiff problems (as compared to the explicit solvers). An important property of implicit solvers is that they require the *computation of a Jacobian* as part of their routine. This means that the various options for efficient Jacobian computation [described later in this tutorial](@ref ode_simulation_performance_jacobian) are only relevant to implicit solvers.

Finally, we should note that this is not an exact science, and sometimes explicit solvers can successfully solve a stiff problem. E.g. if we change the parameters values of our previous Brusselator model to `[:A => 1.0, :B => 4.0]`, the `Tsit5` will successfully be able to simulate it.


## [ODE solver selection](@id ode_simulation_performance_solvers)
OrdinaryDiffEq implements an unusually large number of ODE solvers, with the performance of the simulation heavily depending on which one is chosen. These are provided as the second argument to the `solve` command, e.g. here we use the `Tsit5` solver to simulate a simple production/degradation loop.
```@example ode_simulation_performance_2
using Catalyst, OrdinaryDiffEq

rn = @reaction_network begin
    (p,d), 0 <--> X
end
oprob = ODEProblem(rn, [:X => 0.1], (0.0,1.0), [:p => 1.0, :d => 0.2])
solve(oprob, Tsit5())
nothing # hide
``` 
Alternatively, if the full DifferentialEquations.jl package (rather than the OrdinaryDiffEq subpackage) is loaded, no solver is required (as one is automatically selected):
```@example ode_simulation_performance_2
using DifferentialEquations
solve(oprob)
nothing # hide
``` 
While the default choice is typically enough for most single simulations, if performance is important, it can be worthwhile exploring the available solvers to find one that is especially suited for the given problem. A complete list of possible ODE solvers, with advice on optimal selection, can be found [here](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/). This section will give some general advice.

The most important part of solver selection is to select one appropriate for [the problem's stiffness](@ref ode_simulation_performance_stiffness). Generally, the `Tsit5` solver is recommended for non-stiff problems, and `Rodas5P` for stiff problems. For large stiff problems (with many species), `QNDF` can be a good choice. We can illustrate the impact of these choices by simulating our production/degradation model using the `Tsit5`, `BS5` (an explicit solver yielding [low error in the solution](@ref ode_simulation_performance_error)), `Rodas5P`, and `QNDF` solvers (benchmarking their respective performance using [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl)):
```@example ode_simulation_performance_2
using BenchmarkTools
@btime solve(oprob, Tsit5())
@btime solve(oprob, BS5())
@btime solve(oprob, Rodas5P())
@btime solve(oprob, QNDF())
```
Here we note that the fastest solver is several times faster than the slowest one (`QNDF`, which is a poor choice for this ODE),

### [Simulation error, tolerance, and solver selection](@id ode_simulation_performance_error)
Numerical ODE simulations [approximate an ODE's continuous solution as a discrete vector](https://en.wikipedia.org/wiki/Discrete_time_and_continuous_time). This introduces errors in the computed solution. The magnitude of these errors can be controlled by setting solver *tolerances*. By reducing the tolerance, solution errors will be reduced, however, this will also increase simulation runtimes. The (absolute and relative) tolerance of a solver can be tuned through the `abstol` and `reltol` arguments. Here we see how runtime increases with larger tolerances:
```@example ode_simulation_performance_2
@btime solve(oprob, Tsit5(); abstol=1e-6, reltol=1e-6)
@btime solve(oprob, Tsit5(); abstol=1e-12, reltol=1e-12)
```
It should be noted, however, that the result of the second simulation will be a lot more accurate. Thus, ODE solver performance cannot be determined solely from runtime, but is a composite of runtime and error. Benchmarks comparing the performance (by plotting the runtime as a function of the error) of various ODE solvers, at various tolerances, for various CRN models, can be found in the [SciMLBenchmarks repository](https://docs.sciml.ai/SciMLBenchmarksOutput/stable/Bio/BCR/). 

Generally, whether solution error is a consideration depends on the application. If you want to compute the trajectory of an expensive space probe as it is sent from Earth, to slingshot Jupiter, and then reach Pluto a few years later, keeping error low will be essential. However, if you want to simulate a simple CRN to determine whether it oscillates for a given parameter set, a small error will not constitute a problem. An important aspect with regard to error is that it affects the selection of the optimal solver. E.g. if tolerance is low (generating larger errors) the `Rosenbrock23` method performs well for small, stiff, problems (again, more details can be found in [OrdinaryDiffEq's documentation](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)).


## [Jacobian computation options for implicit solvers](@id ode_simulation_performance_jacobian)
As [previously mentioned](@ref ode_simulation_performance_stiffness), implicit ODE solvers require the computation of the system's [Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant). The reason is that, roughly, in each time step, these solvers need to solve a non-linear equation to find the simulation's value at the next timestep (unlike explicit solvers, which can compute the value at the next time step directly). Typically this is done using [the Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method), which requires the Jacobian. Especially for large systems this can be computationally expensive (and a potential strain on available memory), in which case one might consider various Jacobian-computation options (as described below). A throughout tutorial on simulating a large, stiff, ODE can be found [here](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/advanced_ode_example/#stiff).

### [Building the Jacobian symbolically](@id ode_simulation_performance_symbolic_jacobian)
By default, OrdinaryDiffEq computes the Jacobian using [*automatic differentiation*](https://en.wikipedia.org/wiki/Automatic_differentiation) (however, using [*finite differences*](https://en.wikipedia.org/wiki/Finite_difference) is [also possible](https://docs.sciml.ai/DiffEqDocs/stable/features/performance_overloads/)). Since Catalyst builds its ODEs symbolically, it is able to *compute an analytic Jacobian symbolically*. Typically, this is advantageous when some rows of the Jacobian are highly dense (that is, some system species interact with most other species of the system).

To use this option, simply set `jac=true` when constructing an `ODEProblem`:
```@example ode_simulation_performance_3
using Catalyst, OrdinaryDiffEq

brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
u0 = [:X => 1.0, :Y => 0.0]
tspan = (0.0, 20.0)
p = [:A => 10.0, :B => 40.0]

oprob = ODEProblem(brusselator, u0, tspan, p; jac=true)
nothing # hide
```

### Using a sparse Jacobian
For a system with $n$ variables, the Jacobian will be an $n\times n$ matrix. This means that, as $n$ becomes large, the Jacobian can become *very* large, potentially causing a significant strain on memory. In these cases, most Jacobian entries are typically $0$. This means that a [*sparse*](https://en.wikipedia.org/wiki/Sparse_matrix) Jacobian (rather than a *dense* one, which is the default) can be advantageous. To designate sparse Jacobian usage, simply provide the `sparse=true` option when constructing an `ODEProblem`:
```@example ode_simulation_performance_3
oprob = ODEProblem(brusselator, u0, tspan, p; sparse=true)
nothing # hide
```

### [Linear solver selection](@id ode_simulation_performance_symbolic_jacobian_linear_solver)
When implicit solvers use e.g. the Newton-Raphson method to, in each simulation time step, solve a (typically non-linear) equation, they actually solve a linearised version of this equation. For this they use a linear solver, the choice of which can impact performance. To specify one, we use the `linsolve` option (given to the solver function, *not* the `solve` command). E.g. to use the `KLUFactorization` linear solver we run
```@example ode_simulation_performance_3
using DifferentialEquations
solve(oprob, Rodas5P(linsolve = KLUFactorization()))
nothing # hide
```
Please note that this requires the full DifferentialEquations package (rather than just the more lightweight OrdinaryDiffEq). A full list of potential linear solvers can be found [here](https://docs.sciml.ai/LinearSolve/dev/solvers/solvers/#Full-List-of-Methods), however, the default choice typically performs well. 

A unique approach to the linear solvers is to use a Jacobian-free Newton-Krylov method. These do not actually compute the Jacobian, but rather *the effect of multiplying it with a vector*. They are typically advantageous for large systems (with large Jacobians), and can be designated using the `KrylovJL_GMRES` linear solver:
```@example ode_simulation_performance_3
solve(oprob, Rodas5P(linsolve = KrylovJL_GMRES()))
nothing # hide
```
Since these methods do not depend on a Jacobian, certain Jacobian options (such as [computing it symbolically](@ref ode_simulation_performance_symbolic_jacobian)) are irrelevant to them. 

### [Designating preconditioners for Jacobian-free linear solvers](@ref ode_simulation_performance_preconditioners)
When an implicit method solves a linear equation through an iterative method, the rate of convergence depends on the numerical properties of the matrix defining the linear system. To speed up convergence, a [*preconditioner*](https://en.wikipedia.org/wiki/Preconditioner) can be applied to both sides of the linear equation, attempting to create an equation that converges faster. In practice, preconditioners are implemented as functions with a specific set of arguments. How to implement these is non-trivial, and we recommend reading OrdinaryDiffEq's documentation pages [here](https://docs.sciml.ai/DiffEqDocs/stable/features/linear_nonlinear/#Preconditioners:-precs-Specification) and [here](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/advanced_ode_example/#Adding-a-Preconditioner). In this example, we will define an [Incomplete LU](https://en.wikipedia.org/wiki/Incomplete_LU_factorization) preconditioner (which requires the [IncompleteLU.jl](https://github.com/haampie/IncompleteLU.jl) package):
```@example ode_simulation_performance_3
using IncompleteLU
function incompletelu(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = ilu(convert(AbstractMatrix, W), τ = 50.0)
    else
        Pl = Plprev
    end
    Pl, nothing
end
nothing # hide
```
Next, `incompletelu` can be supplied to our solver using the `precs` argument:
```@example ode_simulation_performance_3
solve(oprob, Rodas5P(precs = incompletelu))
nothing # hide
```
Finally, we note that since matrix-free linear solvers (like `KrylovJL_GMRES`) by default do not build a Jacobian. Hence, if we want to use them with a preconditioner we must tell them to build it. This can be done using the `concrete_jacobian=true` option:
```@example ode_simulation_performance_3
solve(oprob, Rodas5P(linsolve = KrylovJL_GMRES(), precs = incompletelu, concrete_jac=true))
nothing # hide
```

Generally, the use of preconditioners is only recommended for advanced users who are familiar with the concepts. However, for large systems, if performance is essential, they can still be worth looking into.

## [Parallelisation on CPUs and GPUs](@id ode_simulation_performance_parallelisation)
Whenever an ODE is simulated a large number of times (e.g. when investigating its behaviour for different parameter values), the best way to improve performance is to [parallelise the simulation over several processing units](https://en.wikipedia.org/wiki/Parallel_computing). Indeed, an advantage of the Julia programming language is that it was designed after the advent of parallel computing, making it well-suited for this task. Roughly, parallelisation can be divided into parallelisation on [CPUs](https://en.wikipedia.org/wiki/Central_processing_unit) and on [GPUs](https://en.wikipedia.org/wiki/General-purpose_computing_on_graphics_processing_units). CPU parallelisation is most straightforward, while GPU parallelisation requires specialised ODE solvers (which Catalyst have access to). 

Both CPU and GPU parallelisation require first building an `EnsembleProblem` (which defines the simulations you wish to perform) and then supplying this with the correct parallelisation options. These have [previously been introduced in Catalyst's documentation](@id advanced_simulations_montecarlo_simulations) (but in the context of convenient bundling of similar simulations, rather than to improve performance), with a more throughout description being found in [OrdinaryDiffEq's documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#ensemble). Finally, a general documentation of parallel computing in Julia is available [here](https://docs.julialang.org/en/v1/manual/parallel-computing/).

### [CPU parallelisation](@id ode_simulation_performance_parallelisation_CPU)
For this example (and the one for GPUs), we will consider a simple model of an enzyme ($E$) that converts a substrate ($S$) to a product ($P$):
```@example ode_simulation_performance_4
using Catalyst
rn = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
    d, S --> ∅
end
```
The model can be simulated, showing how $P$ is produced from $S$:
```@example ode_simulation_performance_3
using OrdinaryDiffEq, Plots
u0 = [:S => 1.0, :E => 1.0, :SE => 0.0, :P => 0.0]
p = [:kB => 1.0, :kD => 0.1, :kP => 0.5, :d => 0.1]
oprob = ODEProblem(rn, u0, (0.0, 50.0), p)
sol = solve(oprob, Tsit5())
plot(sol)
```
Due to the degradation of $S$, if the production rate is not high enough, the total amount of $P$ produced is reduced. For these tutorial, we will investigate this effect for a range of values of $kP$. This will require a large number of simulations (for various $kP$ values), which we will parallelise on CPUs (here) and GPUs ([later](@ref ode_simulation_performance_parallelisation_GPU)).

To parallelise our model simulations, we first need to create an `EnsembleProblem`. These describe which simulations we wish to perform. The input to this is:
- The `ODEProblem` corresponding to the model simulation (`SDEProblem` and `JumpProblem`s can also be supplied, enabling the parallelisation of these problem types).
- A function, `prob_func`, describing how to modify the problem for each simulation. If we wish to simulate the same, unmodified problem, in each simulation (primarily relevant for stochastic simulations), this argument is not required.

Here, `prob_func` takes 3 arguments:
- `prob`: The problem that it modifies at the start of each individual run (which will be the same as `EnsembleProblem`'s first argument).
- `i`: The index of the specific simulation.
- `repeat`: The repeat of the specific simulations. We will ignore this one for now, but more details (like how it is different from `i`) is provided [here](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#Building-a-Problem).

and output the `ODEProblem` simulated in the ith simulation.

Let us assume that we wish to simulate our model 100 times, for $kP = 0.01, 0.02, ..., 0.99, 1.0$. We can now define our `prob_func` as:
```@example ode_simulation_performance_4
function prob_func(prob, i, repeat) 
    prob[:kP] = 0.01*i
    return prob
end
nothing # hide
```
Next, we can create our `EnsembleProblem`:
```@example ode_simulation_performance_4
eprob = EnsembleProblem(oprob; prob_func=prob_func)
nothing # hide
```
We can now solve our `ODEProblem` using the same syntax we would use to solve the original `ODEProblem`, with the exception that an additional argument, `trajectories`, is required (which denotes how many simulations should be performed).
```@example ode_simulation_performance_4
esol = solve(eprob, Tsit5(); trajectories=100)
nothing # hide
```
to access the i'th solution we use `esol.u[i]`. To e.g. plot the 47'th solution we use:
```@example ode_simulation_performance_4
plot(esol.u[47])
```
To extract the amount of $P$ produced in each simulation, and plot this against the corresponding $kP$ value, we can use:
```@example ode_simulation_performance_4
plot(0.01:0.01:1.0, map(sol -> sol[:P][end], esol.u), xguide="kP", yguide="P produced", label="")
```

Above, we have simply used `EnsembleProblem` as a convenient interface to run a large number of similar simulations. However, these problems have the advantage that they allow the passing of an *ensemble algorithm* to the `solve` command, which describes a strategy for parallelising the simulations. By default, `EnsembleThreads` is used. This parallelises the simulations using [multithreading](https://en.wikipedia.org/wiki/Multithreading_(computer_architecture)) (parallelisation within a single process), which is typically advantageous for small problems. An alternative is `EnsembleDistributed` which instead parallelises the simulations using [multiprocessing](https://en.wikipedia.org/wiki/Multiprocessing) (parallelisation across multiple processes). To do this, we simply supply this additional solver to the solve command:
```@example ode_simulation_performance_4
esol = solve(eprob, Tsit5(), EnsembleDistributed(); trajectories=100)
nothing # hide
```
To utilise multiple processes, you must first give Julia access to these. First, you can check how many processes are available using the `nprocs` (which requires the [Distributed.jl](https://github.com/JuliaLang/Distributed.jl) package):
```@example ode_simulation_performance_4
using Distributed
nprocs()
```
Next, more processes can be added using `addprocs`, e.g. here we add an additional 4 processes:
```@example ode_simulation_performance_4
addprocs(4)
nothing # hide
```
Powerful personal computers and HPC clusters typically have a large number of available additional processes that can be added to improve performance.

While `EnsembleThreads` and `EnsembleDistributed` cover the main cases, additional ensemble algorithms exist. A more throughout description of these can be found [here](https://docs.sciml.ai/DiffEqDocs/dev/features/ensemble/#EnsembleAlgorithms).

Finally, it should be noted that OrdainryDiffEq, if additional processes are available, automatically parallelises the [linear solve part of implicit simulations](@ref ode_simulation_performance_symbolic_jacobian_linear_solver). It is thus possible to see performance improvements from adding additional processes on single simulations, even without running multiple simulations in parallel (this effect is primarily noticeable for large systems with many species).

### [GPU parallelisation](@id ode_simulation_performance_parallelisation_GPU)
GPUs are different from CPUs in that they are much more restricted in what computations they can carry out. However, unlike CPUs, they are typically available in far larger numbers. Their original purpose is for rendering graphics (which typically involves solving a large number of very simple computations, something CPUs with their few, but powerful, cores are unsuited for). Recently, they have also started to be applied to other problems, such as the simulation of ODEs. Generally, GPU parallelisation is only worthwhile when you have a very large number of parallel simulations to run (and access to good GPU resources, either locally or on a cluster).

Generally, we can parallelise `EnsembleProblem`s across several GPUs in a very similar manner to how we parallelised them across several CPUs, but by using a different ensemble algorithm (such as `EnsembleGPUArray`). However, there are some additional requirements:
- GPU parallelisation requires using the [DiffEqGPU.jl](https://github.com/SciML/DiffEqGPU.jl) package.
- Depending on which GPU hardware is used, a specific back-end package has to be installed and imported (e.g. CUDA for NVIDIA's GPUs or Metal for Apple's).
- For some cases, we must use a special ODE solver supporting simulations on GPUs.

Furthermore, to receive good performance, we should also make the following adaptations:
- By default, Julia's decimal numbers are implemented as `Float64`s, however, using `Float32`s is advantageous on GPUs. Ideally, all initial conditions and parameter values should be specified using these.
- We should designate all our vectors (i.e. initial conditions and parameter values) as [static vectors](https://github.com/JuliaArrays/StaticArrays.jl).

We will assume that we are using the CUDA GPU hardware, so we will first load the [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) backend package, as well as DiffEqGPU:
```@example ode_simulation_performance_5
using CUDA, DiffEqGPU
```
Which backend package you should use depends on your available hardware, with the alternative being listed [here](https://docs.sciml.ai/DiffEqGPU/stable/manual/backends/).

Next, we declare our model and `ODEProblem`. However, we make all values `Float64` (by appending `f0` to them) and all vectors static (by adding `@SVector` before their declaration, something which requires the [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) package).
```@example ode_simulation_performance_5
using Catalyst, OrdinaryDiffEq, StaticArrays

rn = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
    d, S --> ∅
end

using OrdinaryDiffEq, Plots
u0 = @SVector [:S => 1.0f0, :E => 1.0f0, :SE => 0.0f0, :P => 0.0f0]
p = @SVector [:kB => 1.0f0, :kD => 0.1f0, :kP => 0.5f0, :d => 0.1f0]
oprob = ODEProblem(rn, u0, (0.0f0, 50.0f0), p)
nothing # hide
```
When we declare our `prob_func` and `EnsembleProblem` we need to ensure that the updated `ODEProblem` uses `Float32`:
```@example ode_simulation_performance_5
function prob_func(prob, i, repeat) 
    prob[:kP] = 0.01f0*i
    return prob
end
eprob = EnsembleProblem(oprob; prob_func=prob_func)
nothing # hide
```

We can now simulate our model using a GPU-based ensemble algorithm. Currently, two such algorithms are available, `EnsembleGPUArray` and `EnsembleGPUKernel`. Their differences are that
* Only `EnsembleGPUKernel` requires arrays to be static arrays (although it is still advantageous for `EnsembleGPUArray`).
* While `EnsembleGPUArray` can use standard ODE solvers, `EnsembleGPUKernel` requires specialised versions (such as `GPUTsit5`). A list of available such solvers can be found [here](https://docs.sciml.ai/DiffEqGPU/dev/manual/ensemblegpukernel/#specialsolvers).

Generally, it is recommended to use `EnsembleGPUArray` for large models (that have at least $100$ variables), and `EnsembleGPUKernel` for smaller ones. Here we simulate our model using both approaches (noting that `EnsembleGPUKernel` requires `GPUTsit5`):
```@example ode_simulation_performance_5
esol1 = solve(eprob, Tsit5(), EnsembleGPUArray(CUDA.CUDABackend()); trajectories=100)
esol2 = solve(eprob, GPUTsit5(), EnsembleGPUKernel(CUDA.CUDABackend()); trajectories=100)
nothing # hide
```
Note that we have to provide the `CUDA.CUDABackend()` argument to our ensemble algorithms (to designate our GPU backend, in this case, CUDA).

Just like OrdinaryDiffEq is able to utilise parallel CPU processes to speed up the linear solve part of ODE simulations, GPUs can also be used. More details on this can be found [here](https://docs.sciml.ai/DiffEqGPU/stable/tutorials/within_method_gpu/). This is only recommended when ODEs are very large (at least 1,000 species), and typically not applicable to CRNs.

---
## References
