# [Coupling chemical reaction networks with differential or algebraic equations](@id coupled_equations)
Catalyst is primarily designed to model [chemical reaction networks (CRNs)](@ref ref), which are defined by a set of species and reaction events. These can then be used to generate e.g. ODEs which can be simulated. However, there are a large number phenomena, many relevant to the fields utilising CRNs, that cannot be described by this formalism. To model these, more general equations (either [differential](https://en.wikipedia.org/wiki/Differential_equation) or [algebraic](https://en.wikipedia.org/wiki/Algebraic_equation)) are needed. Here, Catalyst provides the functionality of coupling CRN models with algebraic or differential equations models. For creating pure non-CRN models, we recommend using [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl).

## [Basic example](@id coupled_equations_basic_example)
Let us consider a simple model of a growth factor ($G$), which controls a cell's volume ($V$). Here, we can model the production and degradation of the growth factor using a simple CRN in Catalyst:
```@example coupled_eqs_1
using Catalyst
crn = @reaction_network begin
    (p,d), 0 <--> G
end
```
CRNs, however, model quantities that exists in (potentially very large) discrete quantities (like number of protein molecules) and change in discrete reaction events (like production and degradation events). This does not hold for a cell's volume. Here we will instead model the cell volume using the following differential equation:
```math
\frac{d\mathbf{V}}{dt} = \mathbf{G(t)} - \mathbf{G(V)},
```
In this case, we can insert this differential equation directly into our model using the `@equations` option, and `D(V)` to describe the differential with respect to $V$:
```@example coupled_eqs_1
coupled_crn = @reaction_network begin
    @equations begin
        D(V) ~ G - V
    end
    (p,d), 0 <--> G
end
nothing # hide
```
We can check the ODE this model would generate [using Latexify](@ref ref):
```@example coupled_eqs_1
using Latexify
latexify(coupled_crn; form=:ode)
```
Finally, the model can be simualted using standard syntax (while providing initial conditions for $V$ as well as for $G$):
```@example coupled_eqs_1
using OrdinaryDiffEq, Plots
u0 = [:G => 0.1, :V => 1.0]
tspan = (0.0, 10.0)
ps = [:p => 1.0, :D => 0.2]
oprob = ODEProblem(coupled_crn, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

## Variables and species
The unknowns of CRN models are called *species*. These are unknown variables of a system which can partake in reaction events. Catalyst support another form of unknowns, *variables*. While variables (like species) may have values that change continuously throughout a simulation, they cannot be substrates or products of a reaction. They can, however, be the subject of differential equations. Species, *cannot* be part of differential equations.

## [Coupling CRNs with differential equation](@id coupled_equations_differential_equations)
Like [demonstrated in the previous example](@ref coupled_equations_basic_example), differential equations can be added to CRN models through the `@equations` option (in the DSL, non-DSL approaches for coupled CRN/equations [are described later](@ref ref)). These are then saved as part of the generated `ReactionSystem` model, and added to the system of ODEs generated from the CRN. The `@equations` option is followed by a `begin ... end` block, where each line correspond to a single equation. Here we create a simple growth model which describes the growth of a cell's volume ($V$), which depends on a growth factor ($G$), and also depletes a nutrient supply ($N$):
```@example coupled_eqs_2
growth_model = @reaction_network begin
    @equations begin
        D(V) ~ G - V
        D(N) ~ -V
    end
    (N*p,d), 0 <--> G
end
latexify(growth_model; form=:ode)
```
Here, while $G$ is inferred to be a species, $V$ and $N$ are inferred to be variables. We can check this using the `species` and `variables` functions:
```@example coupled_eqs_2
species(growth_model)
```
```@example coupled_eqs_2
variables(growth_model)
```
Previously, we have described how Catalyst automatically infers what is a species and what is a parameter depending on where in the model it occurs. Catalyst can also infer variables, however, the only situation where this is done is *terms that occur as the subject of the default differential on the left-hand side of a differential equation* E.g. above both $D$ and $V$ are inferred to be variables.Catalyst does not infer any other components to be species or variables. This means that if we want to include parameters in our differential equations, these must explicitly be declared as such:
```@example coupled_eqs_2
growth_model = @reaction_network begin
    @parameters g d
    @equations begin
        D(V) ~ g*G - V
        D(N) ~ -d*V
    end
    (N*p,d), 0 <--> G
end
parameters(growth_model)
```
Here, with the following exceptions:
- Variables can neither be substrates or products of reactions.
- Species cannot be subject of differentials in differential equations.
- Neither species nor variables can occur as stoichiometric coefficients within reactions

species and variables may occur throughout the model. E.g. in our example above, the production rate of $G$ depends on $N$, while the growth grate of $V$ depends on $G$.

Generally, differential equations are formed similarly as in the [ModelingToolkit.jl package](https://github.com/SciML/ModelingToolkit.jl) (on which [Catalyst is built](@ref ref)). In practise, they follow exactly the same rules as when forming expressions for [rates in Catalyst](@ref ref). The only noteworthy exceptions are:
- The `~` (not `=`!) is used to separate the left-hand and right-hand sides of the equations.
- By default, `D(...)` is used to represent the differential operator with respect to time (or the [default independent variable](@ref ref)).

Differential equations does not need to be written in a pure `D(...) ~ ...` form (with the differential isolated on the left-hand side). However, if this is not the case, `D` is not automatically inferred to be the differential operator. Here, a differential must be [manually declared](@ref ref). Higher-order differentials are also permitted, how to model these are described [here](@ref ref). In these cases, what is a variable is not inferred at all, and this must hence be declared explicitly.

## [Coupling CRNs with algebraic equation](@id coupled_equations_algebraic_equations)
Catalyst also permits the coupling of CRN models with *algebraic equations* (equations not containing differentials). In practise, these are handled similarly to how differential equations are handled, but with the following differences:
- The equations does not contain any differentials.
- The `structural_simplify = true` option must be provided to any `ODEProblem` (or other problem) generated from a problem containing an algebraic equation.
- Catalyst will not infer what is a variable

A special case of algebraic equations (where a new variable is trivially produced by an algebraic expression) are *observables*. These are described in a separate section [here](@ref ref).

## [Converting coupled crn/equations to non-ODE systems](@id coupled_equations_conversions)
As [previously described](@ref ref), Catalyst permits the conversion of `ReactionSystem`s to a range of model types. WIth the exception of `JumpProblem`s, all these can be generated from CRNs coupled to equation models. Below, we briefly describe each case.

### [Converting coupled models to SteadyState problems](@id coupled_equations_conversions_steady_state)

### [Converting coupled models to NonlinearProblems](@id coupled_equations_conversions_nonlinear)

### [Converting coupled models to SDEs](@id coupled_equations_conversions_SDEs)
Generally, Catalyst models containing variables (whenever these values' are given through differential/algebraic equations, or remain constant) can be converted to SDEs. Here, while the chemical Langevin equation is sued to generate noise for all species, equations for variables will remain deterministic. However, if the variables depend on species values, noise from these may still filter through.

!!! note
    Adding support for directly adding noise to variables during SDE conversion is a work in progress. If this is a feature you are interested in, please raise an issue.

### [Other applications of coupled CRN/equation models](@id coupled_equations_conversions_extensions)
In other sections, we describe how to perform [bifurcation analysis](@ref ref), [identifiability analysis](@ref ref), or [compute steady states using homotopy continuation](@ref ref) using chemical reaction network models. Generally, most types of analysis supported by Catalyst is supported for coupled models. below follows some exceptions and notes:
- Methods related to steady state analysis internally converts all systems to `NonlinearSystem`s. This means that all differential terms are set to `0`.
- Homotopy continuation can only be applied if all equations also corresponds to (rational) polynomials.


## [Special cases](@id coupled_equations_other)

### [Defining custom differentials](@id coupled_equations_other_custom_differentials)

### [Using higher-order differentials](@id coupled_equations_other_higher_order_differentials)


## [Coupled equations for programmatic modelling](@id coupled_equations_programmatic)