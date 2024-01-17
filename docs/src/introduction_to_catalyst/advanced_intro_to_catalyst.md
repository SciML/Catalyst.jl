# [How Catalyst works and advanced introduction](@id advanced_intro_to_catalyst) 
We have [previously](@ref introduction_to_catalyst) introduced Catalyst, and how to use it to create and simulate simple chemical reaction network (CRN) models. Here, we provide a more in-depth introduction. We will explain how Catalyst is constructed (i.e. what packages it builds on and how), provide a comprehensive overview of the ways Catalyst models can be created, and describe how these are different.

This advanced introduction is not required for understanding the remaining documentation, nor for using any of Catalyst's features. It is not expected that new users will read it. Instead, it aims to provide additional context and understanding for more experienced users. Various section here are also used as references at later part of the documentation (for users who wants additional insight). For those who expect to use Catalyst extensively, we believe that, at some point, reading this tutorial will be useful.


## [How Catalyst is constructed](@id advanced_intro_to_catalyst_background)
Catalyst is built on top of the domain-agnostic modelling package [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl). ModelingToolkit is in turn built on the computer algebraic system [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl). Here we will introduce these packages, and how the feed into Catalyst.

### [The Symbolics.jl computer algebraic system](@id advanced_intro_to_catalyst_background_symbolics)
Symbolics is a Julia package for representing, and doing computations on, symbolic expressions. It allow us to declare symbolic variables that can be used to construct equations. Prinarily this is done through the `@variables` macro:
```@example advanced_intro_to_catalyst_symbolics
using Symbolics
@variables x y
nothing # hide
```
Here, both `x` and `y` becomes variables in the current scope, e.g. we can do
```@example advanced_intro_to_catalyst_symbolics
x
```
and it prints `x`, just like any other stored Julia value. Next, we can use these variables to create algebraic equations:
```@example advanced_intro_to_catalyst_symbolics
eq = x^2 - y^2
```
Symbolics primary feature is that it can do algebraic manipulation on algebraic expressions. E.g. we can use the `simplify` function to get a simpler (but equivalent) expression. E.g. it can correctly apply the [difference of two squares](https://en.wikipedia.org/wiki/Difference_of_two_squares) rule to our previous equation:
```@example advanced_intro_to_catalyst_symbolics
simplify(eq/(x+y))
```

!!! warning
    `@variables` (and similar macros introduced later on, like `@parameters` and `@species`) create new variables in the current scope. E.g. if you have previously defined `x = 5` and then run `@variables x`, this will overwrite the previous value of `x`.

Symbolics is used internally within Catalyst to represent all models. One advantage of this is that complicated algebraic expressions can be automatically simplified. E.g. if we create a model with a (trivially) unnecessarily complicated rate and check what it is, we see that it has been automatically simplified
```@example advanced_intro_to_catalyst_symbolics
rn = @reaction_network begin
    (k1*k2)/k2, X --> 0
end
reactions(rn)[1].rate
```
Especially when large ODEs are formed, this can remove redundant terms, making the ODE faster to evaluate, thus increasing simulation performance.

Another feature of Symbolics.jl is that it differentiate equations:
```@example advanced_intro_to_catalyst_symbolics
Symbolics.derivative(eq, x)
```
This is used internally to symbolically compute system Jacobians. As explained in a [later section]() this can sometimes be advantageous for simulations.


### [The ModelingToolkit.jl intermediary representation and modelling package](@id advanced_intro_to_catalyst_background_mtk)
ModelingToolkit.jl is a package with two separate functionalities:
- It provides an interface for building such mathematical models.
- It implements an intermediary representation for representing such mathematical models.

We will start with describing the first point, by showing how ModelingToolkit can be used to build a ODE model. First we need to designate the variables and parameters of the model. This is done in a similar manner as we defined variables using symbolics, however, with two differences:
- ODEs (typically) describe the evolution of its variables with time. The variables are thus *time dependent* and must be declared as functions of the time *independent variable*, $t$ (dependency on non-time variables is possible and [discussed later]()).
- Parameters must be declared using the `@parameters` macro (to distinguish them from teh variables).
```@example advanced_intro_to_catalyst_mtk
using ModelingToolkit
@variables X(t)
@parameters p d
```
Like previously, this creates `X`, `p`, and `d` as variables in our current scope:
```@example advanced_intro_to_catalyst_mtk
X
```
```@example advanced_intro_to_catalyst_mtk
p
```

Next, we wish to form the following differential equation:
```math
\frac{dX}{Xt} = p - d*X
```
Here, the derivative with anything with respect to $t$ is defined using the `D` function, i.e. we can form our differential equation as
```@example advanced_intro_to_catalyst_mtk
diff_eq = D(X) ~ p - d*X
```
Finally, to be able to use our differential equation for any purpose (e.g. to make simulations), we must first encapsulate it within an `ODESystem`. This takes vector listing all differential equations, variables, and parameters. It also requires designating the independent variable (here `t`) as the second argument:
```@example advanced_intro_to_catalyst_mtk
osys = ODESystem([diff_eq], t, [X], [p, d]; name = :simple_ODE)
```
The `ODESystem` also have a name (this is a required argument), which can be found using the `nameof` function:
```@example advanced_intro_to_catalyst_mtk
nameof(osys)
```
Sometimes the creation of an `ODESystem` is prepended with the `@named` macro. This automatically gives the system teh same name as the variable it is assigned to:
```@example advanced_intro_to_catalyst_mtk
@named osys = ODESystem([diff_eq], t, [X], [p, d])
nameof(osys)
```

Finally, `ODESystem`s can passed to e.g. `ODEProblems` to create simulations. This is very similar to how it is done for Catalyst, however, we can actually use the variables and parameter themselves to designate their values in the initial condition and parameter values vectors:

```@example advanced_intro_to_catalyst_mtk
using OrdinaryDiffEq, Plots

u0 = [X => 1.0]
tspan = (0.0, 10.0)
ps = [p => 1.0, d => 0.2]

oprob = ODEProblem(osys, u0, tspan, ps)
sol = solve(oprob, Tsit5())
plot(sol)
```

### [System types and their conversions](@id advanced_intro_to_catalyst_background_systems)
We have previously noted that Catalyst CRN models (generated by e.g. `@reaction_network`) are stored in so-called `ReactionSystem` structures. Above we also noted that ModelingToolkit's `ODESystem` can be used as input to `ODEProblem`s in the same manner as Catalyst's `ReactionSystem`. What happens is that `ReactionSystem`s can be converted to `ODESystem`:
```@example advanced_intro_to_catalyst_mtk
using Catalyst
rsys = @reaction_network begin
    (p,d), 0 <--> X
end
osys_catalyst = convert(ODESystem, rsys)
```
Here, when a `ReactionSystem` is supplied to a `ODEProblem`, it is internally converted to a `ODESystem`, which is then used to perform the simulation. Indeed, we can simulate our converted `REactionSystem` and get the same result as when we simulated an `ODESystem` directly.
```@example advanced_intro_to_catalyst_mtk
oprob_catalyst = ODEProblem(osys_catalyst, u0, tspan, ps)
sol_catalyst = solve(oprob_catalyst, Tsit5())
sol == sol_catalyst
```

ModelingToolkit implements a wide range of [system types](https://docs.sciml.ai/ModelingToolkit/stable/basics/AbstractSystem/), enabling to to model a wide range of systems. This include e.g.:
- [`SDESystem`s](https://docs.sciml.ai/ModelingToolkit/stable/systems/SDESystem/) to represent SDEs.
- [`JumpSystem`s](https://docs.sciml.ai/ModelingToolkit/stable/systems/JumpSystem/) to represent jump processes.
- [`NonlinearSystem`s](https://docs.sciml.ai/ModelingToolkit/stable/systems/NonlinearSystem/) to represent systems of nonlinear equations.

The advantage of representing a model as a CRN is that from this representation there is an ubiquitous conversion to ODEs (via the [reaction rate equation](@ref ref)), SDEs (via the [chemical Langevin equation](@ref ref)), and jump processes (via the [stochastic chemical kinetics](@ref ref)). In the same manner, the advantage of representing a model as a Catalyst `ReactionSystem` is that these can ubiquitous be converted to `ODESystem`s, `SDESystem`s, `JumpSystem`s, and `NonlinearSystem`s (the last one primarily used to find the steady states of the ODE form). I.e. all of these are possible:
```@example advanced_intro_to_catalyst_mtk
ssys = convert(SDESystem, rsys)
jsys = convert(JumpSystem, rsys)
nsys = convert(NonlinearSystem, rsys)
nothing # hide
```
with it next being possible to pass these into their respective problem types, e.g.
```@example advanced_intro_to_catalyst_mtk
using StochasticDiffEq
sprob = SDEProblem(ssys, u0, tspan, ps)
ssol = solve(sprob, ImplicitEM())
nothing # hide
```
Practically, in all these cases, the `ReactionSystem` can be passed directly into the desired problem (and the conversion being carried out automatically).

### [ModelingToolkit as an intermediary representation](@id advanced_intro_to_catalyst_background_ir)
Above we have demonstrated ModelingToolkit's second feature, that it provides an *intermediary representation* (IR) for representing mathematical models (as e.g. `ODESystem`s). Here, a large number of packages for building mathematical models in the Julia programming languages (e.g. Catalyst, [OrbitalTrajectories.jl](https://github.com/dpad/OrbitalTrajectories.jl), and [NBodySimulator.jl](https://github.com/SciML/NBodySimulator.jl)) builds their models using ModelingToolkit's IRs. Next, a large number of packages mode model simulation and analysis (e.g. [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl), [StructuralIdentifiability.jl](https://github.com/SciML/StructuralIdentifiability.jl), and [PEtab.jl](https://github.com/sebapersson/PEtab.jl)) ensures that their packages also are compatible with the IRs.

This is the primary way Catalyst is able to support such a large number of features. By utilising the ModelingToolkit IR, we automatically receive support from one of the largest ecosystems for mathematical modelling, while requiring only minimal implementation within Catalyst itself.

## [Types of modelling in Catalyst](@id advanced_intro_to_catalyst_modelling_approaches)
There are three different ways of building Catalyst `ReactionSystem` models. The first two are distinct from each other, while the third one can build on top of either (or both) of these. Here we will briefly describe all three approached (with more details of each described in separate sections of this documentation).

### [DSL-based modelling](@id advanced_intro_to_catalyst_modelling_approaches_DSL)
Catalyst implements a *domain-specific language* (DSL) through its `@reaction_network` macro. This enables the writing of CRN models in their native chemical reaction form. here we use it to create a simple [catalysis model]:
```@example advanced_intro_to_catalyst_modelling_approaches
using Catalyst
rs_dsl = @reaction_network rs begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
```
A throughout description of how to model CRNs using the DSL is provided [here](@ref ref), with some of its more advanced options being described [here](@ref ref).

### [Programmatic modelling](@id advanced_intro_to_catalyst_modelling_approaches_programmatic)
Previously we showed how to programmatically create `ODESystems` [using ModelingToolkit](@ref advanced_intro_to_catalyst_background_mtk). Catalyst have a similar interface for creating `ReactionSystem`s. Here, we first designate the species and parameters of our model (note that we use the `@species` macro instead of the `@variables` macro, which is similar, but specifically designate that the created variables are *species*).
```@example advanced_intro_to_catalyst_modelling_approaches
@species S(t) E(t) SE(t) P(t) 
@parameters kB kD kP
nothing # hide
```
Next, we create the reactions of the system. For each reaction we provide the rate, the substrates, and the products.
```@example advanced_intro_to_catalyst_modelling_approaches
rxs = [
    Reaction(kB, [S, E], [SE]), 
    Reaction(kD, [SE], [S, E]), 
    Reaction(kP, [SE], [P, E])
]
nothing # hide
```
Finally, we create a `ReactionSystem` by combing the reactions, and also designating the [independent variable](@ref advanced_intro_to_catalyst_indep_vars):
```@example advanced_intro_to_catalyst_modelling_approaches
rs_prog = ReactionSystem(rxs, t; name=:rs)
```
The resulting `ReactionSystem` is identical to the one we created directly using the DSL:
```@example advanced_intro_to_catalyst_modelling_approaches
rs_dsl == rs_prog
```
This equality requires the two systems to have the same `name`, here `:rs`. This is provided to the DSL using [the initial argument to @reaction_network](@ref ref) a to the programmatic construction using the `name` keyword.

Note that we here directly use symbolic variables. This can be used to form arbitrary expressions as rates (just like we did when we [described Symbolics.jl](@ref @advanced_intro_to_catalyst_background_symbolics)). I.e. to make a reaction directly dependent on `P` we can write
```@example advanced_intro_to_catalyst_modelling_approaches
Reaction(kP*P, [S], [P])
```
More details on programmatic model construction (e.g. how to use non-unitary stoichiometries) can be found [here](@ref programmatic_CRN_construction).

### [Hierarchical and compositional modelling](@id advanced_intro_to_catalyst_modelling_approaches_hierarchical)
The final type of model construction is hierarchical (or compositional) modelling. This involves first creating smaller models, and then composing them together to create a single, more complex model. The component models can be created using either DSL-based or programmatic modelling.

Here, we create a simple model of a single protein ($X$) which can exist in either the nucleus or the cytoplasm. It can be degraded both in the nucleus and the cytoplasm, but is only created in the nucleus. We will first create one model representing the nucleus and one the cytoplasm. 
```@example advanced_intro_to_catalyst_modelling_approaches
nucleus = @reaction_network begin
    p, 0 --> X
    d_n, X --> 0
end
```
```@example advanced_intro_to_catalyst_modelling_approaches
cytoplasm = @reaction_network begin
    d_c, X --> 0
end
```
Next, we create a top-level model with a transportation reaction (transporting $X$ from the nucleus to the cytoplasm). Here, we interpolate the $X$ variable into the top-level model (more details on interpolation is given [here](@ref advanced_intro_to_catalyst_interpolation)).
```@example advanced_intro_to_catalyst_modelling_approaches
model = @reaction_network begin
    e, $(nucleus.X) --> $(nucleus.X)
end
```
Now, we can use the `compose` function to combine all three models into a single one:
```@example advanced_intro_to_catalyst_modelling_approaches
@named model = compose(model, [nucleus, cytoplasm])
```

Hierarchical modelling can be used to model [system spread across several different compartments](@ref ref), to simplify the creation of [models with repetitive components](). More details on hierarchical and compositional modelling is provided [here](@ref compositional_modeling).

## [Designating species and parameters](@id advanced_intro_to_catalyst_designations)
Frequently when working with Catalyst, we need to represent a species or parameter. Primarily this happens when we set the initial conditions and parameter values of a simulation, but it can also be when we [designate which species to plot](@ref ref) or [updates a parameter during a simulation](@ref ref). There are a few different ways to represent species and parameters, each with their pros and cons. Hence, we will give a short overveiw of them here.

### [Symbol-based designation](@id advanced_intro_to_catalyst_designations_syms)

### [Symbolics-based designation](@id advanced_intro_to_catalyst_designations_symbolics)

### [Using @unpack to extract symbolics variables](@id advanced_intro_to_catalyst_designations_unpack)

### [Via systems symbolics-based designation](@id advanced_intro_to_catalyst_designations_sys_symb)


## Other low-level features
Finally, we describe some additional features of Catalyst. These are features that can be skipped, but for which some additional understanding of can be useful.

### [Interpolating variables into the DSL](@id advanced_intro_to_catalyst_interpolation)

### [System completeness](@id advanced_intro_to_catalyst_completeness)

### [Independent variables](@id advanced_intro_to_catalyst_indep_vars)

### [Variable and reaction metadata](@id advanced_intro_to_catalyst_metadata)