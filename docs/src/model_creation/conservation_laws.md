# [Working with Conservation Laws](@id conservation_laws)
```@raw html
<details><summary><strong>Environment setup and package installation</strong></summary>
```
The following code sets up an environment for running the code on this page.
```julia
using Pkg
Pkg.activate(; temp = true) # Creates a temporary environment, which is deleted when the Julia session ends.
Pkg.add("Catalyst")
Pkg.add("Latexify")
Pkg.add("OrdinaryDiffEqDefault")
Pkg.add("Plots")
```
```@raw html
</details>
```
  \

Catalyst can detect, and eliminate for differential-equation based models, *conserved quantities*, i.e. linear combinations of species which are conserved via their chemistry. This functionality is both automatically utilised by Catalyst (e.g. to [remove singular Jacobians during steady state computations](@ref homotopy_continuation_conservation_laws)), but is also available for users to utilise directly (e.g. to potentially [improve simulation performance](@ref ode_simulation_performance_conservation_laws)).

To illustrate conserved quantities, let us consider the following [two-state](@ref basic_CRN_library_two_states) model:
```@example conservation_laws
using Catalyst
rs = @reaction_network begin
 (k₁,k₂), X₁ <--> X₂
end
```
If we simulate it, we note that while the concentrations of $X₁$ and $X₂$ change throughout the simulation, the total concentration of $X$ ($= X₁ + X₂$) is constant:
```@example conservation_laws
using OrdinaryDiffEqDefault, Plots
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
latexify(rs; form = :ode, math_delimiters = true) # hide
```
we note that it essentially generates the same equation twice (i.e. $\frac{dX₁(t)}{dt} = -\frac{dX₂(t)}{dt}$). By designating our conserved quantity $X₁ + X₂ = Γ$, we can rewrite our differential equation model as a [differential-algebraic equation](https://en.wikipedia.org/wiki/Differential-algebraic_system_of_equations) (with a single differential equation and a single algebraic equation):
```math
\frac{dX₁(t)}{dt} = - k₁X₁(t) + k₂(-X₁(t) + Γ) \\
X₂(t) = -X₁(t) + Γ
```
Using Catalyst, it is possible to detect any such conserved quantities and eliminate them from the system. Here, when we convert our `ReactionSystem` to an ODE `System`, we provide the `remove_conserved = true` argument to instruct Catalyst to perform this elimination:
```@example conservation_laws
osys = ode_model(rs; remove_conserved = true)
latexify(Catalyst.system_to_reactionsystem(osys; disable_forbidden_symbol_check = true); math_delimiters = true) # hide
```
We note that the output system only contains a single (differential) equation and can hence be solved with an ODE solver. The second (algebraic) equation is stored as an [*observable*](@ref dsl_advanced_options_observables), and can be retrieved using the `observed` function:
```@example conservation_laws
observed(osys)
```
Using the `unknowns` function we can confirm that the ODE only has a single unknown variable:
```@example conservation_laws
unknowns(osys)
```
Next, using `parameters` we note that an additional (vector) parameter, `Γ` has been added to the system:
```@example conservation_laws
parameters(osys)
```
Here, Catalyst encodes all conserved quantities in a single, [vector-valued](@ref dsl_advanced_options_vector_variables), parameter called `Γ`.

Practically, the `remove_conserved = true` argument can be provided when a `ReactionSystem` is converted to an `ODEProblem`:
```@example conservation_laws
using OrdinaryDiffEqDefault, Plots
u0 = [:X₁ => 80.0, :X₂ => 20.0]
ps = [:k₁ => 10.0, :k₂ => 2.0]
oprob = ODEProblem(rs, u0, (0.0, 1.0), ps; remove_conserved = true)
nothing # hide
```
Here, while `Γ` becomes a parameter of the new system, it has a default value equal to the corresponding conservation law. Hence, its value is computed from the initial condition `[:X₁ => 80.0, :X₂ => 20.0]`, and does not need to be provided in the parameter vector. Next, we can simulate and plot our model using normal syntax:
```@example conservation_laws
sol = solve(oprob)
plot(sol)
```
!!! note
    Any species eliminated using `remove_conserved = true` will not be plotted by default. However, it can be added to the plot using [the `idxs` plotting option](@ref simulation_plotting_options). E.g. here we would use `plot(sol; idxs = [:X₁, :X₂])` to plot both species.

While `X₂` is an observable (and not unknown) of the ODE, we can [access it](@ref simulation_structure_interfacing_problems) just like if `remove_conserved = true` had not been used:
```@example conservation_laws
sol[:X₂]
```
!!! note
    Generally, `remove_conserved = true` should not change any modelling workflows. I.e. anything that works without this option should also work when an `ODEProblem` is created using `remove_conserved = true`.

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

## [Updating conservation law values directly](@id conservation_laws_prob_updating)
Previously we noted that conservation law elimination adds the conservation law as a system parameter (named $Γ$). E.g. here we find it among the parameters of the ODE model generated by the previously used two-state system:
```@example conservation_laws
parameters(ode_model(rs; remove_conserved = true))
```
An advantage of this is that we can set conserved quantities's values directly in simulations. Here we simulate the model while specifying that $X₁ + X₂ = 10.0$ (and while also specifying an initial condition for $X₁$, but none for $X₂$).
```@example conservation_laws
u0 = [:X₁ => 6.0]
ps = [:k₁ => 1.0, :k₂ => 2.0, :Γ => [10.0]]
oprob = ODEProblem(rs, u0, 1.0, ps; remove_conserved = true)
sol = solve(oprob)
nothing # hide
```
Generally, for each conservation law, one can omit specifying either the conservation law constant, or one initial condition (whichever quantity is missing will be computed from the other ones). 
!!! note
    As previously mentioned, the conservation law parameter $Γ$ is a [*vector-valued* parameter](@ref dsl_advanced_options_vector_variables) with one value for each conservation law. That is why we above provide its value as a vector (`:Γ => [5.0]`). If there had been multiple conservation laws, we would provide the `:Γ` vector with one value for each of them (e.g. `:Γ => [10.0, 15.0]`).
!!! warn
    If you specify the value of a conservation law parameter, you *must not* specify the value of all species of that conservation law (this can result in an error). Instead, the value of exactly one species must be left unspecified.

Just like when we create a problem, if we [update the species (or conservation law parameter) values of `oprob`](@ref simulation_structure_interfacing_problems), the remaining ones will be recomputed to generate an accurate conservation law.
```@example conservation_laws
u0 = [:X₁ => 6.0, :X₂ => 4.0]
ps = [:k₁ => 1.0, :k₂ => 2.0]
oprob = ODEProblem(rs, u0, 10.0, ps; remove_conserved = true)
oprob = remake(oprob; u0 = [:X₁ => 16.0])
```
!!! warn
    If a problem with conservation laws have had values related to that conservation law updated using `remake`, the values stored within the problem will no longer appear correct. I.e. `oprob[:X₂]` will not necessarily yield the correct value. The correct values are, however, computed correctly during `solve`, and it is only when checking the content of the problem that erroneous values appears.

It is also possible to update the value of $Γ$. Here, as $X₂$ is the species eliminated by the conservation law (which can be checked using `conservedequations`), $X₂$'s value will be modified to ensure that $Γ$'s new value is correct. This, however, also requires designating `X₂ = nothing`
```@example conservation_laws
oprob = remake(oprob; u0 = [:X₂ => nothing], p = [:Γ => [30.0]] )
oprob[:X₂]
```

Generally, for a conservation law where we have:
- One (or more) non-eliminated species.
- One eliminated species.
- A conservation law parameter.
If the value of the conservation law parameter $Γ$'s value *has never been specified*, its value will be updated to accommodate changes in species values. If the conservation law parameter ($Γ$)'s value *has been specified* (either when the `ODEProblem` was created, or using `remake`), then the eliminated species's value will be updated to accommodate changes in the conservation law parameter or non-eliminated species's values. Furthermore, in this case, the value of the eliminated species *cannot be updated*. 

!!! warn
    When updating the values of problems with conservation laws it is additionally important to use `remake` (as opposed to direct indexing, e.g. setting `oprob[:X₁] = 16.0`).

### [Extracting the conservation law parameter's symbolic variable](@id conservation_laws_prob_updating_symvar)
Catalyst represents its models using the [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) computer algebraic system, something which allows the user to [form symbolic expressions of model components](@ref simulation_structure_interfacing_symbolic_representation). If you wish to extract and use the symbolic variable corresponding to a model's conserved quantities, you can use the following syntax:
```@example conservation_laws
Γ = Catalyst.get_networkproperties(rs).conservedconst
```