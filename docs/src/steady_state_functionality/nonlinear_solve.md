# [Finding Steady States using NonlinearSolve.jl and SteadyStateDiffEq.jl](@id steady_state_solving)

Catalyst `ReactionSystem` models can be converted to ODEs (through [the reaction rate equation](@ref introduction_to_catalyst_ratelaws)). We have previously described how these ODEs' steady states can be found through [homotopy continuation](@ref homotopy_continuation). Generally, homotopy continuation (due to its ability to find *all* of a system's steady states) is the preferred approach. However, Catalyst supports two additional approaches for finding steady states:

- Through solving the nonlinear system produced by setting all ODE differentials to 0[^1].
- Through forward ODE simulation from an initial condition until a steady state has been reached.

While these approaches only find a single steady state, they offer two advantages as compared to homotopy continuation:

- They are typically much faster.
- They can find steady states for models that do not produce multivariate, rational, polynomial systems (which is a requirement for homotopy continuation to work). Examples include models with non-integer hill coefficients.

In practice, model steady states are found through [nonlinear system solving](@ref steady_state_solving_nonlinear) by creating a `NonlinearProblem`, and through forward ODE simulation by creating a `SteadyStateProblem`. These are then solved through solvers implemented in the [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl), package (with the latter approach also requiring the [SteadyStateDiffEq.jl](https://github.com/SciML/SteadyStateDiffEq.jl) package). This tutorial describes how to find steady states through these two approaches. More extensive descriptions of available solvers and options can be found in [NonlinearSolve's documentation](https://docs.sciml.ai/NonlinearSolve/stable/).

## [Steady state finding through nonlinear solving](@id steady_state_solving_nonlinear)

Let us consider a simple dimerisation system, where a protein ($P$) can exist in a monomer and a dimer form. The protein is produced at a constant rate from its mRNA, which is also produced at a constant rate.

```@example steady_state_solving_nonlinear
using Catalyst
dimer_production = @reaction_network begin
    pₘ, 0 --> mRNA
    pₚ, mRNA --> mRNA + P
    (k₁, k₂), 2P <--> P₂
    d, (mRNA, P, P₂) --> 0
end
```

This system corresponds to the following ODE:

```math
\begin{aligned}
\frac{dmRNA}{dt} &= pₘ - d \cdot mRNA \\
\frac{dP}{dt} &= pₚ \cdot mRNA - k₁ \cdot P + 2k₂ \cdot P₂ - d \cdot P \\
\frac{dP₂}{dt} &= k₁ \cdot P + 2k₂ \cdot P₂ \\
\end{aligned}
```

To find its steady states we need to solve:

```math
\begin{aligned}
0 &= pₘ - d \cdot mRNA \\
0 &= pₚ \cdot mRNA - k₁ \cdot P + 2k₂ \cdot P₂ - d \cdot P \\
0 &= k₁ \cdot P + 2k₂ \cdot P₂ \\
\end{aligned}
```

To solve this problem, we must first designate our parameter values, and also make an initial guess of the solution. Generally, for problems with a single solution (like this one), most arbitrary guesses will work fine (the exception typically being [systems with conservation laws](@ref steady_state_solving_nonlinear_conservation_laws)). Using these, we can create the `NonlinearProblem` that we wish to solve.

```@example steady_state_solving_nonlinear
p = [:pₘ => 0.5, :pₚ => 2.0, :k₁ => 5.0, :k₂ => 1.0, :d => 1.0]
u_guess = [:mRNA => 1.0, :P => 1.0, :P₂ => 1.0]
nlprob = NonlinearProblem(dimer_production, u_guess, p)
nothing # hide
```

Finally, we can solve it using the `solve` command, returning the steady state solution:

```@example steady_state_solving_nonlinear
using NonlinearSolve
sol = solve(nlprob)
```

Typically, a good default method is automatically selected for any problem. However, NonlinearSolve does provide [a wide range of potential solvers](https://docs.sciml.ai/NonlinearSolve/stable/solvers/NonlinearSystemSolvers/). If we wish to designate one, it can be supplied as a second argument to `solve`. Here, we use the [Newton Trust Region method](https://en.wikipedia.org/wiki/Trust_region), and then check that the solution is equal to the previous one.

```@example steady_state_solving_nonlinear
sol_ntr = solve(nlprob, TrustRegion())
sol ≈ sol_ntr
```

### [Systems with conservation laws](@id steady_state_solving_nonlinear_conservation_laws)

As described in the section on homotopy continuation, when finding the steady states of systems with conservation laws, [additional considerations have to be taken](@ref homotopy_continuation_conservation_laws). E.g. consider the following [two-state system](@ref basic_CRN_library_two_states):

```@example steady_state_solving_claws
using Catalyst, NonlinearSolve # hide
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
```

It has an infinite number of steady states. To make steady state finding possible, information of the system's conserved quantities (here $C = X1 + X2$) must be provided. Since these can be computed from system initial conditions (`u0`, i.e. those provided when performing ODE simulations), designating an `u0` is often the best way. There are two ways to do this. First, one can perform [forward ODE simulation-based steady state finding](@ref steady_state_solving_simulation), using the initial condition as the initial `u` guess. Alternatively, any conserved quantities can be eliminated when the `NonlinearProblem` is created. This feature is supported by [Catalyst's conservation law finding and elimination feature](@ref conservation_laws).

To eliminate conservation laws we simply provide the `remove_conserved = true` argument to `NonlinearProblem`:

```@example steady_state_solving_claws
p = [:k1 => 2.0, :k2 => 3.0]
u_guess = [:X1 => 3.0, :X2 => 1.0]
nl_prob = NonlinearProblem(two_state_model, u_guess, p; remove_conserved = true)
nothing # hide
```

here it is important that the quantities used in `u_guess` correspond to the conserved quantities we wish to use. E.g. here the conserved quantity $X1 + X2 = 3.0 + 1.0 = 4$ holds for the initial condition, and will hence also hold in the computed steady state as well. We can now find the steady states using `solve` like before:

```@example steady_state_solving_claws
sol = solve(nl_prob)
```

We note that the output only provides a single value. The reason is that the actual system solved only contains a single equation (the other being eliminated with the conserved quantity). To find the values of $X1$ and $X2$ we can [directly query the solution object for these species' values, using the species themselves as inputs](@ref simulation_structure_interfacing_solutions):

```@example steady_state_solving_claws
sol[[:X1, :X2]]
```

## [Finding steady states through ODE simulations](@id steady_state_solving_simulation)

The `NonlinearProblem`s generated by Catalyst corresponds to ODEs. A common method of solving these is to simulate the ODE from an initial condition until a steady state is reached. Here we do so for the dimerisation system considered in the previous section. First, we declare our model, initial condition, and parameter values.

```@example steady_state_solving_simulation
using Catalyst # hide
dimer_production = @reaction_network begin
    pₘ, 0 --> mRNA
    pₚ, mRNA --> mRNA + P
    (k₁, k₂), 2P <--> P₂
    d, (mRNA, P, P₂) --> 0
end
p = [:pₘ => 0.5, :pₚ => 2.0, :k₁ => 5.0, :k₂ => 1.0, :d => 1.0]
u0 = [:mRNA => 0.1, :P => 0.0, :P₂ => 0.0]
nothing # hide
```

Next, we provide these as an input to a `SteadyStateProblem`

```@example steady_state_solving_simulation
ssprob = SteadyStateProblem(dimer_production, u0, p)
nothing # hide
```

Finally, we can find the steady states using the `solver` command. Since `SteadyStateProblem`s are solved through forward ODE simulation, we must load the sub-library of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package that corresponds to the [selected ODE solver](@ref simulation_intro_solver_options). Any available ODE solver can be used, however, it has to be encapsulated by the `DynamicSS()` function. E.g. here we use the `Rodas5P` solver which is loaded from the `OrdinaryDiffEqRosenbrock` sub-library:

(which requires loading the SteadyStateDiffEq package).

```@example steady_state_solving_simulation
using SteadyStateDiffEq, OrdinaryDiffEqRosenbrock
solve(ssprob, DynamicSS(Rodas5P()))
```

Note that, unlike for nonlinear system solving, `u0` is not just an initial guess of the solution, but the initial conditions from which the steady state simulation is carried out. This means that, for a system with multiple steady states, we can determine the steady states associated with specific initial conditions (which is not possible when the nonlinear solving approach is used). This also permits us to easily [handle the presence of conservation laws](@ref steady_state_solving_nonlinear_conservation_laws). The forward ODE simulation approach (unlike homotopy continuation and nonlinear solving) cannot find unstable steady states.

Generally, `SteadyStateProblem`s can be solved using the [same options that are available for ODE simulations](@ref simulation_intro_solver_options). E.g. here we designate a specific `dt` step size:

```@example steady_state_solving_simulation
solve(ssprob, DynamicSS(Rodas5P()); dt = 0.01)
nothing # hide
```

It is possible to use solve `SteadyStateProblem`s using a nonlinear solver, and `NonlinearProblem`s using forward ODE simulation solvers:

```@example steady_state_solving_simulation
using NonlinearSolve
solve(ssprob, TrustRegion())
nothing # hide
```

```@example steady_state_solving_simulation
nlprob = NonlinearProblem(dimer_production, u0, p)
solve(nlprob, DynamicSS(Rodas5P()))
nothing # hide
```

However, especially when the forward ODE simulation approach is used, it is recommended to use the problem type which corresponds to the intended solver.

---

## [Citations](@id nonlinear_solve_citation)

If you use this functionality in your research, [in addition to Catalyst](@ref doc_index_citation), please cite the following paper to support the authors of the NonlinearSolve.jl package:

```bibtex
@article{pal2024nonlinearsolve,
  title={NonlinearSolve. jl: High-Performance and Robust Solvers for Systems of Nonlinear Equations in Julia},
  author={Pal, Avik and Holtorf, Flemming and Larsson, Axel and Loman, Torkel and Schaefer, Frank and Qu, Qingyu and Edelman, Alan and Rackauckas, Chris and others},
  journal={arXiv preprint arXiv:2403.16341},
  year={2024}
}
```

---

## References

[^1]: [J. Nocedal, S. J. Wright, *Numerical Optimization*, Springer (2006).](https://www.math.uci.edu/~qnie/Publications/NumericalOptimization.pdf)
