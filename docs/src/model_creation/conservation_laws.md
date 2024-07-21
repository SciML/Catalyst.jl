# [Working with conservation laws](@id conservation_laws)
Catalyst contain specialised functionality to work with *conserved quantities*, i.e. species quantities which are conserved throughout e.g. a simulation. This functionality is both automatically utilised by Catalyst (e.g. to [remove singular Jacobians during steady state computations](@ref homotopy_continuation_conservation_laws)), but also something that is available for the user to utilise directly (e.g. to [improve simulation performance](@ref ode_simulation_performance_conservation_laws)).

To illustrate conserved quantities, let us consider the following [two-state](@ref basic_CRN_library_two_states) model:
```@example conservation_laws
using Catalyst
rs = @reaction_network begin
 (k₁,k₂), X₁ <--> X₂
end
```
If we simulate it, we note that while the concentrations of $X₁$ and $X₂$ change throughout the simulation, the total concentration of $X$ ($= X₁ + X₂$) is constant:
```@example conservation_laws
using OrdinaryDiffEq, Plots
u0 = [:X₁ => 80.0, :X₂ => 20.0]
ps = [:k₁ => 10.0, :k₂ => 2.0]
oprob = ODEProblem(rs, u0, (0.0, 1.0), ps)
sol = solve(oprob)
plot(sol; idxs = [rs.X₁, rs.X₂, rs.X₁ + rs.X₂], label = ["X₁" "X₂" "X₁ + X₂ (a conserved quantity)"])
```
This makes sense, as while $X$ is converted between two different forms ($X₁$ and $X₂$), it is neither produced nor degraded. That is, it is a *conserved quantity*. Next, if we consider the ODE that our model generates:
```@example conservation_laws
using Latexify
latexify(rs; form = :ode)
```
we note that it essentially generates the same equation twice (i.e. $\frac{dX₁(t)}{dt} = -\frac{dX₂(t)}{dt}$). By designating our conserved quantity $X₁ + X₂ = Γ$, we can rewrite our differential equation model as a [differential-algebraic equation](https://en.wikipedia.org/wiki/Differential-algebraic_system_of_equations) (with a single differential equation and a single algebraic equation):
```math
\frac{dX₁(t)}{dt} = - k₁X₁(t) + k₂(-X₁(t) + Γ) \\
X₂(t) = -X₁(t) + Γ
```
Using Catalyst, it is possible to detect any such conserved quantities and eliminate them from the system. Here, when we convert our `ReactionSystem` to an `ODESystem`, we provide the `remove_conserved = true` argument to instruct Catalyst to perform this elimination:
```@example conservation_laws
osys = convert(ODESystem, rs; remove_conserved = true)
```
We note that the output system only contains a single (differential) equation and can hence be solved with an ODE solver. The second (algebraic) equation is stored as an [*observable*](@ref dsl_advanced_options_observables), and can be retrieved using the `observed` function:
```@example conservation_laws
observed(osys)
```
Using the `unknowns` function we can confirm that the ODE only has a single unknown variable:
```@example conservation_laws
unknowns(osys)
```
Next, using `parameters` we note that an additional parameter, `Γ[1]` has been added to the system:
```@example conservation_laws
parameters(osys)
```
Here, Catalyst encodes all conserved quantities in a single, [vector-valued](@ref dsl_advanced_options_vector_variables), parameter called `Γ`.

Practically, the `remove_conserved = true` argument can be provided when a `ReactionSystem` is converted to an `ODEProblem`:
```@example conservation_laws
using OrdinaryDiffEq, Plots
u0 = [:X₁ => 80.0, :X₂ => 20.0]
ps = [:k₁ => 10.0, :k₂ => 2.0]
oprob = ODEProblem(rs, u0, (0.0, 1.0), ps; remove_conserved = true)
nothing # hide
```
Here, while `Γ[1]` becomes a parameter of the new system, it has a [default value](@ref dsl_advanced_options_default_vals) equal to the corresponding conservation law. Hence, its value is computed from the initial condition `[:X₁ => 80.0, :X₂ => 20.0]`, and does not need to be provided in the parameter vector. Next, we can simulate and plot our model using normal syntax:
```@example conservation_laws
sol = solve(oprob)
plot(sol)
```
!!! note
    Any species eliminated using `remove_conserved = true` will not be plotted by default. However, it can be added to the plot using [the `idxs` plotting option](@ref simulation_plotting_options). E.g. here we would use `plot(sol; idxs = [:X₁, :X₂])` to plot both species.

!!! danger
    Currently, there is a bug in MTK where the values associated with conservation laws are not updated properly in response to [`remake`](@ref simulation_structure_interfacing_problems_remake) (or [other problem-updating functions, such as `getu`](@ref simulation_structure_interfacing_functions)). Hence, problems created using `remove_conserved = true` *should not* be modified, i.e. by changing initial conditions or the values of the `Γ`'s directly. Instead, to change these values new problems must be generated with the new initial condition map.

While `X₂` is an observable (and not unknown) of the ODE, we can [access it](@ref simulation_structure_interfacing_problems) just like if `remove_conserved = true` had not been used:
```@example conservation_laws
sol[:X₂]
```
!!! note
    Generally, `remove_conserved = true` should not change any model workflows. I.e. anything that works without this option should also work when an `ODEProblem` is created using `remove_conserved = true`.

!!! note
    The `remove_conserved = true` option is available when creating `SDEProblem`s, `NonlinearProblem`s, and `SteadyStateProblem`s (and their corresponding systems). However, it cannot be used when creating `JumpProblem`s.

## [Conservation law accessor functions](@id conservation_laws_accessors)

For any given `ReactionSystem` model, we can use `conservationlaw_constants` to compute all of a system's conserved quantities:
```@example conservation_laws
conservationlaw_constants(rs)
```
Next, the `conservedequations` function can be used to compute the algebraic equations produced when a system's conserved quantities are eliminated:
```@example conservation_laws
conservedequations(rs)
```
Finally, the `conservationlaws` function yields a $m \times n$ matrix, where $n$ is a system's number of species, $m$ its number of conservation laws, and element $(i,j)$ corresponds to the copy number of the $j$th species that occurs in the $i$th conserved quantity:
```@example conservation_laws
conservationlaws(rs)
```
I.e. in this case we have a single conserved quantity, which contains a single copy each of the system's two species.