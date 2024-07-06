# [Advice for performant SDE simulations](@id sde_simulation_performance)
While there exist relatively straightforward approaches to manage performance for [ODE](@ref ode_simulation_performance) and [Jump](@ref ref) simulations, this is generally not the case for SDE simulations. Below, we briefly describe some options. However, as one starts to investigate these, one quickly reaches what is (or could be) active areas of research.

## [SDE solver selection](@id sde_simulation_performance_solvers)
We have previously described how [ODE solver selection](@ref ode_simulation_performance_solvers) can impact simulation performance. Again, it can be worthwhile to investigate solver selection's impact on performance for SDE simulations. Throughout this documentation, we generally use the `STrapezoid` solver as the default choice. However, if the `DifferentialEquations` package is loaded
```julia
using DifferentialEquations
```
automatic SDE solver selection enabled (just like is the case for ODEs by default). Generally, the automatic SDE solver choice enabled by `DifferentialEquations` is better than just using `STrapezoid`. Next, if performance is critical, it can be worthwhile to check the [list of available SDE solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/) to find one with advantageous performance for a given problem. When doing so, it is important to pick a solver compatible with *non-diagonal noise* and with [*Ito problems*](https://en.wikipedia.org/wiki/It%C3%B4_calculus). 

## [Options for Jacobian computation](@id sde_simulation_performance_jacobian)
In the section on ODE simulation performance, we describe various [options for computing the system Jacobian](@ref ode_simulation_performance_jacobian), and how these could be used to improve performance for [implicit solvers](@ref ode_simulation_performance_stiffness). These can be used in tandem with implicit SDE solvers (such as `STrapezoid`). However, due to additional considerations during SDE simulations, it is much less certain whether these will actually have any impact on performance. So while these options might be worth reading about and trialling, there is no guarantee that they will be beneficial.

## [Parallelisation on CPUs and GPUs](@id sde_simulation_performance_parallelisation)
We have previously described how simulation parallelisation can be used to [improve performance when multiple ODE simulations are carried out](@ref ode_simulation_performance_parallelisation). The same approaches can be used for SDE simulations. Indeed, it is often more relevant for SDEs, as these are often re-simulated using identical simulation conditions (to investigate their typical behaviour across many samples). CPU parallelisation of SDE simulations uses the [same approach as ODEs](@ref ode_simulation_performance_parallelisation_CPU). GPU parallelisation requires some additional considerations, which are described below.

### [GPU parallelisation of SDE simulations](@id sde_simulation_performance_parallelisation_GPU)
GPU parallelisation of SDE simulations uses a similar approach as that for [ODE simulations](@ref ode_simulation_performance_parallelisation_GPU). The main differences are that SDE parallelisation requires a GPU SDE solver (like `GPUEM`) and fixed time stepping.

We will assume that we are using the CUDA GPU hardware, so we will first load the [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) backend package, as well as DiffEqGPU:
```julia
using CUDA, DiffEqGPU
```
Which backend package you should use depends on your available hardware, with the alternatives being listed [here](https://docs.sciml.ai/DiffEqGPU/stable/manual/backends/).

Next, we create the `SDEProblem` which we wish to simulate. Like for ODEs, we ensure that all vectors are [static vectors](https://github.com/JuliaArrays/StaticArrays.jl) and that all values are `Float32`s. Here we prepare the parallel simulations of a simple [birth-death process](@ref basic_CRN_library_bd).
```@example sde_simulation_performance_gpu
using Catalyst
bd_model = @reaction_network begin
    (p,d), 0 <--> X
end
@unpack X, p, d = bd_model

u0 = @SVector [X => 20.0f0]
tspan = (0.0f0, 10.0f0)
ps = @SVector [p => 10.0f0, d => 1.0f0]
sprob = SDEProblem(bd_model, u0, tspan, ps)
nothing # hide
```
The `SDEProblem` is then used to [create an `EnsembleProblem`](@ref ensemble_simulations_monte_carlo).
```@example sde_simulation_performance_gpu
eprob = EnsembleProblem(sprob)
nothing # hide
```
Finally, we can solve our `EnsembleProblem` while:
- Using a valid GPU SDE solver (either [`GPUEM`](https://docs.sciml.ai/DiffEqGPU/stable/manual/ensemblegpukernel/#DiffEqGPU.GPUEM) or [`GPUSIEA`](https://docs.sciml.ai/DiffEqGPU/stable/manual/ensemblegpukernel/#DiffEqGPU.GPUSIEA)).
- Designating the GPU ensemble method, `EnsembleGPUKernel` (with the correct GPU backend as input).
- Designating the number of trajectories we wish to simulate.
- Designating a fixed time step size.

```julia
esol = solve(eprob, GPUEM(), EnsembleGPUKernel(CUDA.CUDABackend()); trajectories = 10000, dt = 0.01)
```

Above we parallelise GPU simulations with identical initial conditions and parameter values. However, [varying these](@ref ensemble_simulations_varying_conditions) is also possible.

### [Multilevel Monte Carlo](@id sde_simulation_performance_parallelisation_mlmc)
An approach for speeding up parallel stochastic simulations is so-called [*multilevel Monte Carlo approaches*](https://en.wikipedia.org/wiki/Multilevel_Monte_Carlo_method) (MLMC). These are used when a stochastic process is simulated repeatedly using identical simulation conditions. Here, instead of performing all simulations using identical [tolerance](@ref ode_simulation_performance_error), the ensemble is simulated using a range of tolerances (primarily lower ones, which yields faster simulations). Currently, [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) do not have a native implementation for performing MLMC simulations (this will hopefully be added in the future). However, if high performance of parallel SDE simulations is required, these approaches may be worth investigating.