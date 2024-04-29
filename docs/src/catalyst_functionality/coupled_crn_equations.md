# [Coupling chemical reaction networks with differential or algebraic equations](@id coupled_equations)
Catalyst is primarily designed to model [chemical reaction networks (CRNs)](@ref ref), which are defined by a set of *species* and a set of *reaction events*. Once defined, CRNs can be converted to e.g. ODEs or SDEs. However, many phenomena relevant to fields utilising CRNs (e.g. systems biology) cannot be described using this formalism. To model these, more general equations (either [differential](https://en.wikipedia.org/wiki/Differential_equation) or [algebraic](https://en.wikipedia.org/wiki/Algebraic_equation)) are required. Here, Catalyst provides the functionality of coupling CRN models with algebraic or differential equations models. For creating pure non-CRN models, we recommend using [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl).

Example of phenomena which might be modelled with differential equations are:
- The volume of a cell, or its amount of available nutrients.
- The public's awayness of of a novel infectious disease.
- The current temperature in a chemical reactor.

Example of phenomena which might be modelled with algebraic equations are:
- The concentration of species which reactions are "fast enough" that they can be considered to be at steady state.
- The relations between pressure, volume, molecular numbers, and temperature given bye [*the ideal gas law*](https://en.wikipedia.org/wiki/Ideal_gas_law).

### [Basic example](@id coupled_equations_basic_example)
Let us consider a simple model of a growth factor ($G$), which controls a cell's volume ($V$). Here, we can model the production and degradation of the growth factor using a simple CRN in Catalyst:
```@example coupled_eqs_1
using Catalyst
crn = @reaction_network begin
    (p,d), 0 <--> G
end
```
CRNs, however, model quantities that exists in (potentially very large) discrete quantities (like number of protein molecules) and change in discrete reaction events (like production and degradation events). This does not hold for a cell's volume (which is a *continuous quantity*). Here we will instead model the cell volume using the following differential equation:
```math
\frac{d\mathbf{V}}{dt} = \mathbf{G(t)} - \mathbf{V(t)}
```
where the volume's growth increase with $G$ at saturates ones large enough. It is possible to insert this differential equation directly into our model using the DSL's `@equations` option (here, `D(V)` denotes the differential of $V$ with respect to time):
```@example coupled_eqs_basic_example
coupled_crn = @reaction_network begin
    @equations begin
        D(V) ~ G - 1.0
    end
    (p,d), 0 <--> G
end
nothing # hide
```
We can check the ODE this model would generate [using Latexify](@ref ref):
```@example coupled_eqs_basic_example
using Latexify
latexify(coupled_crn; form=:ode)
```
Finally, the model can be simulated using standard syntax (while providing initial conditions for $V$ as well as for $G$):
```@example coupled_eqs_1
using OrdinaryDiffEq, Plots
u0 = [:G => 0.1, :V => 1.0]
tspan = (0.0, 10.0)
ps = [:p => 1.0, :D => 0.2]
oprob = ODEProblem(coupled_crn, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

## [Variables and species](@id coupled_equations_variables)
Catalyst's normal CRN models depend on two types of quantities, *parameters* (which values are explicitly provided at the beginning of simulations) *species* (which values are implicitly inferred from a CRN model by e.g. simulating it). The coupling of CRNs and equations require the introduction of a third quantity, *variables*. Variables, just like species, have values which are implicitly given my a model. For many purposes, species and variables can be treated identically, and a model's set of species and variables is called its *unknowns*.

Variables and species, however, are different in that:
- Only species are permissible reaction reactants.
- Only variables are permitted subjects of differentials in coupled CRN/differential equations.

Previously, we have described how model species can [explicitly be declared using the `@species` DSL option](@ref ref). Model variables can be declared in exactly the same manner, but using the `@variables` option. While quantities declared using `@variables` are considered variables (not species), its syntax (for e.g. designating [default values](@ref ref) or [metadata](@ref ref)) is otherwise identical to `@species`. Below we declare a simple cell model (with a growth factor species $G$ and a volume variable $V$): 
```@example coupled_eqs_variable_intro
using Catalyst # hide
cell_model = @reaction_network begin
  @variables V(t)
  (p,d), 0 <--> G
end
```
While Catalyst, generally, can infer what quantities are species and parameters, this typically does not hold for variables. Catalyst only infers a quantity to be a variable if it is [the subject of a time differential](@ref coupled_equations_differential_equations). In all other situations, variables must explicitly be declared.

Note that variable (like species) are time-dependant, and (like species) must be declared as such (e.g. using `V(t)` rather than `V`) when the `@variables` option is used. We can now simulate our model (remembering to provide a initial conditions for $V$ as well as $G$):
```@example dsl_advanced_variables
using OrdinaryDiffEq, Plots
u0 = [:G => 0.1, :V => 1.0]
tspan = (0.0, 10.0)
p = [:p => 1.0, :d => 0.5]
oprob = ODEProblem(cell_model, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
Here, since we have not inserted $V$ into any equation describing how its value may change, it remains constant throughout the simulations. The following sections will describe how to declare such equation as part of Catalyst models.

Finally, we note that variables (like parameters and species) are valid components of the rate expressions. E.g. we can modify our cell module so that the reaction rates (inversely) depends on the volume variable:
```@example coupled_eqs_variable_intro
cell_model = @reaction_network begin
  @variables V(t)
  (p/V,d/V), 0 <--> G
end
```

It is possible to check a models species, variables, and unknowns through the `species`, `variables`, and `unknowns` functions:
```@example coupled_eqs_variable_intro
species(cell_model)
```
```@example coupled_eqs_variable_intro
variables(cell_model)
```
```@example coupled_eqs_variable_intro
unknowns(cell_model)
```


## [Coupling CRNs with differential equation](@id coupled_equations_differential_equations)
Differential equations can be added to CRN models through the DSL's `@equations` option ([programmatic](@ref ref) creation of couple CRN/equation model [is described later](@ref coupled_equations_programmatic)). These are then saved as part of the generated `ReactionSystem` model, and added to the system of ODEs generated from the CRN. 

Let us consider a simple cell growth model consisting of a single growth factor ($G$) that is produced/degraded according to the following reactions:
```math
\mathbf{N} \cdot \mathbf{p}, \mathbf{∅} \to \mathbf{G} \\
\mathbf{d}, \mathbf{G} \to \mathbf{∅}
```
Our model will also consider the variables $N$ (the amount of available nutrient) and $V$ (the cell's volume). These are described through the following differential equations:
```math
\frac{d\mathbf{N}}{dt} = - \mathbf{V(t)} \\
\frac{d\mathbf{V}}{dt} = \mathbf{G(t)} - \mathbf{V(t)}
```
The reactions can be modelled using Catalyst's normal reaction notation. The equations are added using the `@equations`, which is followed by a`begin ... end` block, where each line contain one equation. The time differential is written using `D(...)` and equality between the left and right-hand sides are denoted using `~` (*not* `=`!).
```@example coupled_eqs_diff_eq
growth_model = @reaction_network begin
    @equations begin
        D(N) ~ -V
        D(V) ~ G - V
    end
    (N*p,d), 0 <--> G
end
```
The equations left and right-hand sides are declared using the same rules as reaction rate (and can be any valid algebraic expression containing variables, species, parameters, numbers, and function calls). We can check the ODE our model generate [by using Latexify](@ref ref):
```@example coupled_eqs_diff_eq
latexify(growth_model; form=:ode)
```

Previously, we described how coupled CRN/equation models contain [both *species* and *variables*](@ref coupled_equations_variables), and how variables can be explicitly declared using the `@variables` option. Here, since both $N$ and $V$ are the subject of the default differential $D$, Catalyst automatically infers these to be variables. We which unknowns are species and variables through the corresponding functions:
```@example coupled_eqs_diff_eq
species(growth_model)
```
```@example coupled_eqs_diff_eq
variables(growth_model)
```

Currently, our differential equations depend on the species and variables $G$, $N$, and $V$. They may also depend on parameters, however, before we introduce additional quantities to these, we must be careful to ensure that Catalyst can correctly infer their type. Catalyst uses the following rules for automatically inferring what to consider a parameter, a species, and a variable:
- Quantities explicitly declared using the `@parameters`, `@species`, and `@variables` options have the corresponding types.
- Undeclared quantities occurring as reaction substrates or products are inferred to be species.
- Undeclared quantities occurring as the single subject of the default differential (`D`) on the left-hand side of an equation is inferred to be a variable.
- Remaining undeclared quantities that occur in either reaction rates or stoichiometries are inferred to be parameters.

Catalyst cannot infer the type of quantities that do not fulfil any of these conditions. This includes parameters that occur only within declared equations. Below we introduce such parameters by also explicitly declaring them as such using the `@parameters` option, and then confirm that they have been added properly using the `parameters` function:
```@example coupled_eqs_diff_eq
growth_model = @reaction_network begin
    @parameters dN dV
    @equations begin
        D(N) ~ - dN * V
        D(V) ~ G - dV * V
    end
    (N*p,d), 0 <--> G
end
parameters(growth_model)
```

Once declared, we can simulate our model using standard syntax (providing initial conditions for all unknowns, whether species or variables):
```@example coupled_eqs_diff_eq
using OrdinaryDiffEq, Plots # hide
u0 = [:G => 0.1, :V => 1.0, :N => 1.0]
tspan = (0.0, 10.0)
p = [:p => 1.0, :d => 0.5, :dN => 0.2, :dV => 0.1]
oprob = ODEProblem(growth_model, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```

Finally, a few additional notes regarding coupled CRN/differential equation models:
- The rules for declaring differential equations are identical to those used by the [ModelingToolkit.jl package](https://github.com/SciML/ModelingToolkit.jl) (on which [Catalyst is built](@ref ref)). Users familiar with this package can directly use this familiarity when declaring equations within Catalyst.
- While the default differential `D` is typically used, it is possible to [define custom differentials](@ref ref).
- Catalyst only infers a quantity `X` to be a variable for expressions on the form `D(X) ~ ...`, where the left-hand side contain no differentials. For equations with multiple, higher-order, or custom differentials, or where the left-hand side contain several terms, no variables are inferred (and these must be declared explicitly using `@variables`).
- Declared differential equations may contain higher-order differential equations, however, these have [some additional considerations](@ref ref).
- Generally `D` is inferred to be the differential with respect to the [system's independent variable](@ref ref) (typically time). However, this is only the case in the same situations when Catalyst infers variables (e.g. equations on the form `D(X) ~ ...`). For other situations, the differential must [explicitly be declared](@ref ref).

## [Coupling CRNs with algebraic equation](@id coupled_equations_algebraic_equations)
Catalyst also permits the coupling of CRN models with *algebraic equations* (equations not containing differentials). In practise, these are handled similarly to how differential equations are handled, but with the following differences:
- The `structural_simplify = true` option must be provided to any `ODEProblem` (or other problem) generated from a model containing algebraic equations.
- Catalyst cannot infer any variables, so these must explicitly be declared.

A special case of algebraic equations (where a new variable is trivially produced by an algebraic expression) are *observables*. These are described in a separate section [here](@ref ref).

## [Converting coupled CRN/equation models to non-ODE systems](@id coupled_equations_conversions)
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