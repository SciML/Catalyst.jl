# [Coupling chemical reaction networks with differential or algebraic equations](@id coupled_equations)
Catalyst is primarily designed to model [chemical reaction networks (CRNs)](@ref ref), which are defined by a set of *species* and a set of *reaction events*. Once defined, CRNs can be converted to e.g. ODEs or SDEs. However, many phenomena relevant to fields utilising CRNs (e.g. systems biology) cannot be described using this formalism. To model these, more general equations (either [differential](https://en.wikipedia.org/wiki/Differential_equation) or [algebraic](https://en.wikipedia.org/wiki/Algebraic_equation)) are required. Here, Catalyst provides the functionality of coupling CRN models with algebraic or differential equations models. For creating pure non-CRN models, we recommend using [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl).

Examples of phenomena which might be modelled with differential equations are:
- The volume of a cell, or its amount of available nutrients.
- The public's awareness of a novel infectious disease.
- The current temperature in a chemical reactor.

Examples of phenomena which might be modelled with algebraic equations are:
- The concentration of species whose reactions are "fast enough" that they can be considered to be at steady state.
- The relations between pressure, volume, molecular numbers, and temperature as given by [*the ideal gas law*](https://en.wikipedia.org/wiki/Ideal_gas_law).

### [Basic example](@id coupled_equations_basic_example)
Let us consider a simple model of a growth factor ($G$), which controls a cell's volume ($V$). Here, we can model the production and degradation of the growth factor using a simple CRN in Catalyst:
```@example coupled_eqs_1
using Catalyst
crn = @reaction_network begin
    (p,d), 0 <--> G
end
```
CRNs, however, model quantities that exist in (potentially very large) discrete quantities (like the number of protein molecules) and change in discrete reaction events (like production and degradation events). This does not hold for a cell's volume (which is a *continuous quantity*). Here we will instead model the cell volume using the following differential equation:
```math
\mathbf{\frac{dV}{dt} = G(t) - V(t)}
```
where the volume's growth increases with $G$ at saturates ones large enough. It is possible to insert this differential equation directly into our model using the DSL's `@equations` option (here, `D(V)` denotes the differential of $V$ with respect to time):
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
Finally, the model can be simulated using standard syntax (while providing initial conditions for both $V$ and $G$):
```@example coupled_eqs_1
using OrdinaryDiffEq, Plots
u0 = [:G => 0.1, :V => 1.0]
tspan = (0.0, 10.0)
ps = [:p => 1.0, :d => 0.2]
oprob = ODEProblem(coupled_crn, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

## [Variables and species](@id coupled_equations_variables)
Catalyst's normal CRN models depend on two types of quantities, *parameters* (which values are explicitly provided at the beginning of simulations) and *species* (which values are implicitly inferred from a CRN model and some initial condition through simulations). The coupling of CRNs and equations requires the introduction of a third quantity, *variables*. Variables, just like species, have values which are implicitly given by a model. For many purposes, species and variables can be treated identically. A model's set of species and variables is called its *unknowns*.

Variables and species, however, are different in that:
- Only species are permissible reaction reactants.
- Only variables are permitted subjects of differentials in coupled CRN/differential equations.

Previously, we have described how model species can [explicitly be declared using the `@species` DSL option](@ref ref). Model variables can be declared in exactly the same manner, but using the `@variables` option. While quantities declared using `@variables` are considered variables (not species), its syntax (for e.g. designating [default values](@ref ref) or [metadata](@ref ref)) is otherwise identical to `@species`. Below we declare a simple cell model (with a growth factor species $G$ and a volume variable $V$): 
```@example coupled_eqs_variable_intro
using Catalyst
cell_model = @reaction_network begin
  @variables V(t)
  (p,d), 0 <--> G
end
```
While Catalyst, generally, can infer what quantities are species and parameters, this typically does not hold for variables. Catalyst only infers a quantity to be a variable if it is [the subject of a time differential](@ref coupled_equations_differential_equations). In all other situations, variables must explicitly be declared.

Note that variables (like species) are time-dependent, and (like species) must be declared as such (e.g. using `V(t)` rather than `V`) when the `@variables` option is used. We can now simulate our model (remembering to provide initial conditions for $V$ as well as $G$):
```@example dsl_advanced_variables
using OrdinaryDiffEq, Plots
u0 = [:G => 0.1, :V => 1.0]
tspan = (0.0, 10.0)
p = [:p => 1.0, :d => 0.5]
oprob = ODEProblem(cell_model, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```
Here, since we have not inserted $V$ into any equation describing how its value may change, it remains constant throughout the simulations. The following sections will describe how to declare such equations as part of Catalyst models.

Finally, we note that variables (like parameters and species) are valid components of the rate expressions. E.g. we can modify our cell module so that the reaction rates (inversely) depend on the volume variable:
```@example coupled_eqs_variable_intro
cell_model = @reaction_network begin
  @variables V(t)
  (p/V,d/V), 0 <--> G
end
```

It is possible to check a model's species, variables, and unknowns through the `species`, `nonspecies`, and `unknowns` functions:
```@example coupled_eqs_variable_intro
species(cell_model)
```
```@example coupled_eqs_variable_intro
nonspecies(cell_model)
```
```@example coupled_eqs_variable_intro
unknowns(cell_model)
```


## [Coupling CRNs with differential equation](@id coupled_equations_differential_equations)
Differential equations can be added to CRN models through the DSL's `@equations` option ([programmatic](@ref ref) creation of couple CRN/equation model [is described later](@ref coupled_equations_programmatic)). These are then saved as part of the generated `ReactionSystem` model, and added to the system of ODEs generated from the CRN. 

Let us consider a simple cell growth model consisting of a single growth factor ($G$) that is produced/degraded according to the following reactions:
```math
\mathbf{N \cdot p, ∅ \to G} \\
\mathbf{d, G \to ∅}
```
Our model will also consider the variables $N$ (the amount of available nutrients) and $V$ (the cell's volume). These are described through the following differential equations:
```math
\mathbf{\frac{dN}{dt} = - V(t)} \\
\mathbf{\frac{dV}{dt} = G(t) - V(t)}
```
The reactions can be modelled using Catalyst's normal reaction notation. The equations are added using the `@equations` option, which is followed by a `begin ... end` block, where each line contains one equation. The time differential is written using `D(...)` and equality between the left and right-hand sides is denoted using `~` (*not* `=`!).
```@example coupled_eqs_diff_eq
using Catalyst # hide
growth_model = @reaction_network begin
    @equations begin
        D(N) ~ -V
        D(V) ~ G - V
    end
    (N*p,d), 0 <--> G
end
```
The equation's left and right-hand sides are declared using the same rules as used for reaction rates (and can be any valid algebraic expression containing variables, species, parameters, numbers, and function calls). We can check the ODE our model generates [using Latexify](@ref ref):
```@example coupled_eqs_diff_eq
latexify(growth_model; form=:ode)
```

Previously, we described how coupled CRN/equation models contain [both *species* and *variables*](@ref coupled_equations_variables), and how variables can be explicitly declared using the `@variables` option. Here, since both $N$ and $V$ are the subject of the default differential $D$, Catalyst automatically infers these to be variables. We can check which unknowns are species and variables through the corresponding functions:
```@example coupled_eqs_diff_eq
species(growth_model)
```
```@example coupled_eqs_diff_eq
nonspecies(growth_model)
```

Currently, our differential equations depend on the species and variables $G$, $N$, and $V$. They may also depend on parameters, however, before we introduce additional quantities to these, we must be careful to ensure that Catalyst can correctly infer their type. Catalyst uses the following rules for automatically inferring what to consider a parameter, a species, and a variable:
- Quantities explicitly declared using the `@parameters`, `@species`, and `@variables` options have the corresponding types.
- Undeclared quantities occurring as reaction substrates or products are inferred to be species.
- Undeclared quantities occurring as the single subject of the default differential (`D`) on the left-hand side of an equation are inferred to be variables.
- Remaining undeclared quantities that occur in either reaction rates or stoichiometries are inferred to be parameters.

Catalyst cannot infer the type of quantities that do not fulfil any of these conditions. This includes parameters that occur only within declared equations. Below we introduce such parameters by also explicitly declaring them as such [using the `@parameters` option](@ref ref), and then confirm that they have been added properly using the `parameters` function:
```@example coupled_eqs_diff_eq
growth_model = @reaction_network begin
    @parameters dN dV
    @equations begin
        D(N) ~ - dN * V * N
        D(V) ~ G - dV * V
    end
    (N+p,d), 0 <--> G
end
parameters(growth_model)
```

Once declared, we can simulate our model using standard syntax (providing initial conditions for all unknowns, whether species or variables):
```@example coupled_eqs_diff_eq
using OrdinaryDiffEq, Plots # hide
u0 = [:G => 0.0, :V => 1.0, :N => 1.0]
tspan = (0.0, 10.0)
p = [:p => 1.0, :d => 0.5, :dN => 0.2, :dV => 1.0]
oprob = ODEProblem(growth_model, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

Finally, a few additional notes regarding coupled CRN/differential equation models:
- If our model only contains a single equation, it can be declared on the same line as the `@equations` option (e.g. `@equations D(V) ~ G - dV * V`), without the use of a `begin ... end` block.
- The rules for declaring differential equations are identical to those used by the [ModelingToolkit.jl package](https://github.com/SciML/ModelingToolkit.jl) (on which [Catalyst is built](@ref ref)). Users familiar with this package can use their knowledge of it when declaring equations within Catalyst.
- While the default differential `D` is typically used, it is possible to [define custom differentials](@ref coupled_equations_other_custom_differentials).
- Catalyst only infers a quantity `X` to be a variable for expressions on the form `D(X) ~ ...`, where the left-hand side contains no differentials. For equations with multiple, higher-order, or custom differentials, or where the left-hand side contains several terms, no variables are inferred (and these must be declared explicitly using `@variables`).
- Declared differential equations may contain higher-order differential equations, however, these have [some additional considerations](@ref coupled_equations_other_higher_order_differentials).
- Generally, `D` is inferred to be the differential with respect to the [system's independent variable](@ref ref) (typically time). However, this is only the case in the same situations when Catalyst infers variables (e.g. equations on the form `D(X) ~ ...`). For other situations, the differential must [explicitly be declared](@ref coupled_equations_other_custom_differentials).

## [Coupling CRNs with algebraic equation](@id coupled_equations_algebraic_equations)
Catalyst also permits the coupling of CRN models with *algebraic equations* (equations not containing differentials). In practice, these are handled similarly to how differential equations are handled, but with the following differences:
- The `structural_simplify = true` option must be provided to any `ODEProblem` (or other problem) generated from a model containing algebraic equations.
- Catalyst cannot infer any variables, so these must explicitly be declared.

A special case of algebraic equations (where a new variable is trivially produced by an algebraic expression) is *observables*. These are described in a separate section [here](@ref ref).

Let us demonstrate the coupling of CRNs to algebraic equations by considering a system with a single component $X$, which can dimerise to form $X2$. $X$ is degraded at a constant rate, and its production is modelled as a [repressive Michaelis-Menten function](@ref ref) of the concentration of $X2$:
```@example coupled_eqs_alg_eq
using Catalyst # hide
dimerisation_model = @reaction_network begin
    (mmr(X2,v,K),d), 0 <--> X
    (kB,kD), 2X <--> X2
end
```
Here, if the dynamics of the binding/unbinding reactions is much faster than the production/degradation reactions (not an implausible assumption), the `2X <--> X2` reaction can be assumed to be at steady state at the time scale of the production/degradation reaction. We can then rephrase our model using the (algebraic) steady state equation to define the concentration of $X2$:
```@example coupled_eqs_alg_eq
algebraic_crn = @reaction_network begin
    @variables X2(t)
    @parameters kB kD
    @equations begin
        kB*X^2 ~ kD*X2
    end
    (mmr(X2,v,K),d), 0 <--> X
end
nothing # hide
```
We can now define initial conditions and parameter values, and use these to build an `ODEProblem`. Note that, since `X2`'s value can be computed by solving a nonlinear problem at the initial condition (using the supplied algebraic equation), its value *should not* need to be supplied to the initial condition vector.
```@example coupled_eqs_alg_eq
u0 = [:X => 1.0]
tspan = (0.0, 10.0)
ps = [:v => 2.0, :K => 1.0, :d => 0.2, :kB => 0.1, :kD => 0.4]
oprob = ODEProblem(algebraic_crn, u0, tspan, ps; structural_simplify = true)
nothing # hide
```
Here, since our model contains an algebraic equation, we *must* add the `structural_simplify = true` option to `ODEProblem` (this option is discussed in more detail [here](@ref ref)). This option eliminates any algebraic equations from the system, creating a pure ODE which can be solved using normal means. Next, our model can be simulated just like any other model:
```@example coupled_eqs_alg_eq
using OrdinaryDiffEq, Plots # hide
sol = solve(oprob)
plot(sol)
```
Here, we note that the variable $X2$ is not included in the plot (since it is internally eliminated from the ODE). However, it can still be plotted by [explicitly telling the `plot` function so](@ref ref):
```@example coupled_eqs_alg_eq
plot(sol; idxs = :X2)
```

## [Converting coupled CRN/equation models to non-ODE systems](@id coupled_equations_conversions)
As [previously described](@ref ref), Catalyst permits the conversion of `ReactionSystem`s to a range of model types. With the exception of `JumpProblem`s, all these can be generated from coupled CRNs/equation models. Below, we briefly describe each case.

### [Converting coupled models to SteadyState problems](@id coupled_equations_conversions_steady_state)
`SteadyStateProblems` are used to [find a system's steady state through forward ODE simulations](@ref ref). Internally, these consider the system as an ODE, and simulate it until it has reached a point where the solution's rate of change with time is negligible. Since these are (internally) based on `ODEProblem`s, they are handled for coupled CRN/equation models exactly like ODEs. Below we demonstrate how to find our growth model's steady state by creating and solving a `SteadyStateProblem`.
```@example coupled_eqs_steadystate_problem
using Catalyst, OrdinaryDiffEq # hide
using SteadyStateDiffEq

# Declare the model.
growth_model = @reaction_network begin
    @parameters dN dV
    @equations begin
        D(N) ~ - dN * V * N
        D(V) ~ G - dV * V
    end
    (N+p,d), 0 <--> G
end

# Create and solves a `SteadyStateProblem`.
u0 = [:G => 0.0, :V => 1.0, :N => 1.0]
p = [:p => 1.0, :d => 0.5, :dN => 0.2, :dV => 1.0]
ssprob = SteadyStateProblem(growth_model, u0, p)
sol = solve(ssprob, DynamicSS(Tsit5()))
```

### [Converting coupled models to NonlinearProblems](@id coupled_equations_conversions_nonlinear)
Like `SteadyStateProblem`s, `NonlinearProblem`s can be used [to find system steady states]@ref ref. However, rather than doing so through forward simulation, these set all differential terms equal to $0$, and then solve the resulting nonlinear system of equations. Again, the same approach can be used for coupled CRN/equation models. Below we demonstrate how to find our growth model's steady state by creating and solving a `NonlinearProblem`.
```@example coupled_eqs_nonlinear_problem
using Catalyst # hide
using NonlinearSolve

# Declare the model.
growth_model = @reaction_network begin
    @parameters dN dV
    @equations begin
        D(N) ~ - dN * V * N
        D(V) ~ G - dV * V
    end
    (N+p,d), 0 <--> G
end

# Create and solves a `NonlinearProblem`.
u0 = [:G => 1.0, :V => 1.0, :N => 1.0]
p = [:p => 1.0, :d => 0.5, :dN => 0.2, :dV => 1.0]
nprob = NonlinearProblem(growth_model, u0, p)
sol = solve(nprob)
```
Note that for `NonlinearProblem`s, `u0` is an initial guess of the steady state solution, rather than an initial condition for its simulation. Hence, we here use a trivial guess where all unknowns are set to $1.0$. 

!!! warn
    For coupled CRN/equation models with [higher-order differential equations](@ref coupled_equations_other_higher_order_differentials) or multiple differentials, Catalyst can still generate `NonlinearProblem`s (and does so by also setting higher-order differentials to $0$). This is, however, not officially supported. Users who use this should manually check whether this approach makes sense for their specific application. To enable this, supply the `all_differentials_permitted = true` option to your `NonlinearProblem`.

### [Converting coupled models to SDEs](@id coupled_equations_conversions_SDEs)
All Catalyst models containing variables (whenever these values are given through differential/algebraic equations, or remain constant) can be converted to SDEs. Here, while the [chemical Langevin equation is used to generate noise for all species](@ref ref), equations for variables will remain deterministic. However, if the variables depend on species values, noise from these may still filter through. Below we demonstrate this by creating and simulating an `SDEProblem`.
```@example coupled_eqs_sdes
using Catalyst, Plots # hide
using StochasticDiffEq

# Declare the model.
growth_model = @reaction_network begin
    @parameters dN dV
    @equations begin
        D(N) ~ - dN * V * N
        D(V) ~ G - dV * V
    end
    (N+p,d), 0 <--> G
end

# Create and simulate an `SDEProblem`.
u0 = [:G => 1.0, :V => 1.0, :N => 1.0]
tspan = (0.0, 10.0)
p = [:p => 0.2, :d => 0.1, :dN => 0.2, :dV => 1.0]
sprob = SDEProblem(growth_model, u0, tspan, p)
sol = solve(sprob, STrapezoid())
sol = solve(sprob, STrapezoid(); seed = 12345) # hide
plot(sol)
```
In the plot, we can confirm that while the equation for $V$ is deterministic, because it is coupled to $G$ (whose behaviour is stochastic), $V$ does exhibit some fluctuations. While the same also holds for $N$, in practise this is so far detached from the noise that its trajectory is mostly deterministic.

!!! note
    Adding support for directly adding noise to variables during SDE conversion is a work in progress. If this is a feature you are interested in, please [raise an issue]([@ref ref](https://github.com/SciML/Catalyst.jl/issues)).

### [Other applications of coupled CRN/equation models](@id coupled_equations_conversions_extensions)
In other sections, we describe how to perform [bifurcation analysis](@ref ref), [identifiability analysis](@ref ref), or [compute steady states using homotopy continuation](@ref ref) using chemical reaction network models. Generally, most types of analysis supported by Catalyst are supported for coupled models. Below follows some exceptions and notes:
- Methods related to steady state analysis internally convert all systems to `NonlinearSystem`s. This means that all differential terms are set to `0`.
- Homotopy continuation can only be applied if all equations also correspond to (rational) polynomials.


## [Special cases](@id coupled_equations_other)

### [Defining custom differentials](@id coupled_equations_other_custom_differentials)
Previously, differential equations have always been declared using the default differential, `D`. There are, however, two situations when non-default differentials are useful:
- When differentials with respect to independent variables other than the time independent variable are required.
- When a differential with a name other than D is required (e.g. because `D` represents another system quantity).

Non-default differentials can be declared using the `@differentials` option. E.g. here we declare our growth model, but call our differential `Δ` instead of `D`:
```@example coupled_eqs_nonlinear_problem
using Catalyst # hide
growth_model = @reaction_network begin
    @differential Δ = Differential(t)
    @variables N(t) V(t)
    @parameters dN dV
    @equations begin
        Δ(N) ~ - dN * V * N
        Δ(V) ~ G - dV * V
    end
    (N+p,d), 0 <--> G
end
```
Here, the left-hand side of the differential contains the differential's name only. The right-hand side is `Differential(t)`, where `t` is the default [time independent variable](@ref ref) used within Catalyst. If you wish to declare a differential with respect to another independent variable, replace `t` with its name. If you need to declare multiple differentials, this can be done by providing a `begin ... end` block (with one differential declaration on each line) to the `@differential` option. Note that since we are no longer using the default differential, $N$ and $V$ must be explicitly declared using the `@variables` option.

### [Using higher-order differentials](@id coupled_equations_other_higher_order_differentials)
Until now we have only declared first-order differential equations. However, higher-order differential equations can also be used. These are mostly handled identically to how first-order differential equations are handled, with three exceptions:
- Variables and differentials are not inferred for these, and must be explicitly declared.
- Like for [coupled CRN/algebraic equation models](@ref coupled_equations_algebraic_equations), the `structural_simplify = true` option must be provided to `ODEProblem`.
- Initial conditions must be provided for any lower-order differentials.

The second point will make things a bit more convoluted, as we will need to explicitly define the differential outside of the DSL to set its initial condition.

To demonstrate workflows using higher-order differentials we declare a simple model featuring these:
```@example coupled_eqs_higher_order_diffs
using Catalyst # hide
higher_order_model = @reaction_network begin
    @variables A(t)
    @parameters ω
    @differentials D = Differential(t)
    @equations begin
        D(D(A)) + ω*A / (X+1) ~ 0
    end
    (p,d), 0 <--> X
end
nothing # hide
```
Next, we must [@unpack](@ref ref) to retrieve all unknowns (as their initial conditions cannot be provided in `Symbol` form simultaneously as differential initial conditions are provided).
```@example coupled_eqs_higher_order_diffs
@unpack A, X = higher_order_model
nothing # hide
```
Next, we define the Differential `D` in our current scope (fetching the same independent variable used with `@reaction_network` from our model using `Catalyst.get_iv`):
```@example coupled_eqs_higher_order_diffs
D = Differential(Catalyst.get_iv(higher_order_model))
nothing # hide
```
Finally, we can define our `ODEProblem`, simulate it, and plot the solution:
```@example coupled_eqs_higher_order_diffs
u0 = [X => 0.1, A => 1.0, D(A) => 1.0]
tspan = (0.0, 10.0)
p = [:p => 1.0, :d => 0.5, :ω => 2.0]
oprob = ODEProblem(higher_order_model, u0, tspan, p; structural_simplify = true)
sol = solve(oprob)
plot(sol)
```
We note that while modelling of higher-order coupled CRN/differential equation models using the DSL is a bit convoluted (due to the requirement of re-declaring the differential outside of the DSL), this is not the case for [programmatic modelling](@ref coupled_equations_programmatic). Hence, for these models, this is instead the recommended approach.

## [Coupled CRN/equation models for programmatic modelling](@id coupled_equations_programmatic)
We have previously described how Catalyst's `ReactionSystem` models can be [created programmatically](@ref ref) (instead of through the DSL). It is possible to create coupled CRN/equation models programmatically in a very similar manner as through the DSL. Here, we modify our normal programmatic approach to:
- Declare any variables using the `@variables` macro.
- Fetching the default time differential using the `default_time_deriv` function. 
- Instead of declaring a reaction vector, we declare a combined reaction/equation vector, and use this one as input to our `ReactionSystem` constructor.

To demonstrate this we recreate our cell growth model, but using a programmatic approach. First, we declare the model's time independent variable, parameters, species, and variables. The variables are declared using the `@variables` macro (which works just like the `@species` macro, but declares its input as variables and not species).
```@example coupled_eqs_programmatic
using Catalyst
t = default_t()
@species G(t)
@variables N(t) V(t)
@parameters p d dN dV
nothing # hide
```
Next, we declare the differential with respect to the independent variable. Here, we do so through the `default_time_deriv` function (although custom ones can be defined using [similar syntax as for the DSL](@ref coupled_equations_other_custom_differentials), but without the `@differential` option).
```@example coupled_eqs_programmatic
D = default_time_deriv()
nothing # hide
```
Next, we declare our equation vector. This vector contains all the model's reactions *and* equations. Reactions are declared using normal programmatic syntax, while equations (whether differential or algebraic) are declared using syntax identical to that of the DSL:
```@example coupled_eqs_programmatic
eqs = [
    Reaction(p + N, [], [G]),
    Reaction(d, [G], []),
    D(N) ~ - dN * V * N,
    D(V) ~ G - dV * V
]
```
Finally, we use these as input to our `ReactionSystem`:
```@example coupled_eqs_programmatic
@named growth_model = ReactionSystem(eqs, t)
growth_model = complete(growth_model)
```
The resulting model can be simulated using standard syntax:
```@example coupled_eqs_programmatic
using OrdinaryDiffEq, Plots
u0 = [:G => 0.0, :V => 1.0, :N => 1.0]
tspan = (0.0, 10.0)
p = [:p => 1.0, :d => 0.5, :dN => 0.2, :dV => 1.0]
oprob = ODEProblem(growth_model, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

Finally, we note that the `reactions` function can be used to retrieve a `ReactionSystem`s reactions only, while `equations` retrieve both reactions and (algebraic and differential) equations
```@example coupled_eqs_programmatic
reactions(growth_model)
```
```@example coupled_eqs_programmatic
equations(growth_model)
```
Here, the vector returned by `equations(growth_model)` is reordered so that reactions occur first, and equations second. This is the case independently of the order these occur in in the input vector:
```@example coupled_eqs_programmatic
eqs = [
    D(N) ~ - dN * V * N,
    D(V) ~ G - dV * V,
    Reaction(p + N, [], [G]),
    Reaction(d, [G], [])
]
@named growth_model = ReactionSystem(eqs, t)
equations(growth_model)
```