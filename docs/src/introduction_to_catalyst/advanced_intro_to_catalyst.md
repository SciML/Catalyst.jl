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
This is used internally to symbolically compute system Jacobians. As explained in a [later section](@ref ref) this can sometimes be advantageous for simulations.


### [The ModelingToolkit.jl intermediary representation and modelling package](@id advanced_intro_to_catalyst_background_mtk)
ModelingToolkit.jl is a package with two separate functionalities:
- It provides an interface for building such mathematical models.
- It implements an intermediary representation for representing such mathematical models.

We will start with describing the first point, by showing how ModelingToolkit can be used to build a ODE model. First we need to designate the variables and parameters of the model. This is done in a similar manner as we defined variables using symbolics, however, with two differences:
- ODEs (typically) describe the evolution of its variables with time. The variables are thus *time dependent* and must be declared as functions of the time *independent variable*, $t$ (dependency on non-time variables is possible and [discussed later]()).
- Parameters must be declared using the `@parameters` macro (to distinguish them from the variables).
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
Sometimes the creation of an `ODESystem` is prepended with the `@named` macro. This automatically gives the system the same name as the variable it is assigned to:
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
Frequently when working with Catalyst, we need to represent a species or parameter (or variables, when [mixed CRN-equation based models are used]()). Primarily this happens when we set the initial conditions and parameter values of a simulation, but it can also be when we [designate which species to plot](@ref ref) or [updates a parameter during a simulation](@ref ref). There are a three different ways to represent species and parameters, each with their pros and cons. Hence, we will give a short overview of them here.

### [Symbol-based designation](@id advanced_intro_to_catalyst_designations_syms)
Throughout this documentation, Symbol-based representation is used for parameters and species. Here, to represent $X$, the corresponding [`Symbol`](@ref ref), `:X` is used (a `Symbol` is created by pre-pending an expression with `:`). E.g. to simulate the following model:
```@example advanced_intro_to_catalyst_designations
using Catalyst
rs_dsl = @reaction_network begin
    (p,d), 0 <--> X
end
```
we can represent our initial conditions and parameters with their respective `Symbol` form:
```@example advanced_intro_to_catalyst_designations
using  OrdinaryDiffEq
u0_syms = [:X => 1.0]
tspan = (0.0, 10.0)
ps_syms = [:p => 1.0, :d => 0.2]
oprob_syms = ODEProblem(rs_dsl, u0_syms, tspan, ps_syms)
sol_syms = solve(oprob_syms, Tsit5())
nothing # hide
```
Next, to obtain a vector with the simulation's $X$ value we can again use $X$'s `Symbol` representation:
```@example advanced_intro_to_catalyst_designations
sol[:X]
```

There are two cases where `Symbol`-based designations cannot be used:
- When using [hierarchical/compositional modelling](advanced_intro_to_catalyst_modelling_approaches_hierarchical). Here, it cannot be inferred the value of $X$ in which submodel `:X` designates (i.e. its value in the nucleus or the cytosol). How to deal with this problem is described [three sections down](@ref advanced_intro_to_catalyst_designations_sys_symb).
- When needing to form an algebraic expression of species/parameters (as used e.g. [here](@ref petab_parameter_fitting_intro_observables)). Potential workaround is to use one of the following two approaches for designating species/parameters (often, the quantity of interest can also be [declared as an *observable*](@ref ref)).

### [Symbolics-based designation](@id advanced_intro_to_catalyst_designations_symbolics)
Previously, we showed how model creation using ModelingToolkit involves [explicitly declaring variables and parameters, which then can be used to designate their values](@ref advanced_intro_to_catalyst_background_mtk).We also showed how Catalyst `ReactionSystem`s can be [created programmatically, in a similar way to how ModelingToolkit models are created](@ref advanced_intro_to_catalyst_modelling_approaches_programmatic). Here, just like for ModelingToolkit, we can represent species and parameters using their respective Symbolics variables directly:
```@example advanced_intro_to_catalyst_designations
@variables t 
@species X(t) 
@parameters p d
rxs = [
    Reaction(p, nothing, [X]), 
    Reaction(d, [X], nothing)
]
rs_prog = ReactionSystem(rxs, t; name=:rs)

u0_symb = [X => 1.0]
ps_symb = [p => 1.0, d => 0.2]
oprob_symb = ODEProblem(rs_prog, u0_symb, tspan, ps_symb)
sol_symb = solve(oprob_symb, Tsit5())
nothing # hide
```
We can check that the resulting solution is the same as when we used `Symbol`-based designation (and DSL-based model creation):
```@example advanced_intro_to_catalyst_designations
sol_sym == sol_symb
```
We can also use the Symbolic form to e.g. index the output solution:
```@example advanced_intro_to_catalyst_designations
sol[X]
```

The advantage with the Symbolics for is that we can form algebraic expressions using this for of the variables/parameters. E.g. if we (for some reason) would like to evaluate the expression $2X + d$ across the solution vector, we can do this directly using the Symbolic form:
```@example advanced_intro_to_catalyst_designations
sol[2X + p]
```
Indeed, we can confirm that this returns the expected output:
```@example advanced_intro_to_catalyst_designations
sol[2X + p] == 2*sol[X] .+ 1.0
```

### [Using @unpack to extract symbolics variables](@id advanced_intro_to_catalyst_designations_unpack)
Because symbolic-based designation requires having access to the species and parameters in their symbolic form, it is primarily associated with programmatic modelling. I.e. you cannot run
```@example advanced_intro_to_catalyst_designations
sol[X]
nothing # hide
```
without first having run
```@example advanced_intro_to_catalyst_designations
@species X(t) 
```
Thus, one cannot directly combine symbolics-based designations with DSL-based modelling. This is, however, possible by either:
1. In addition to the DSL, declare species and parameters you wish to designate symbolically using `@species` and `@parameters`.
2. By using the `@unpack` macro.

Here, `@unpack` permits you to import variables from a `ReactionSystem` into your current scope. This is done by first listing what you wish to import, and then the system you wish to import it from, e.g. you can use:
```@example advanced_intro_to_catalyst_designations
@unpack X, p, d = rs_dsl
nothing # hide
```
to import `X`, `p`, and `d` to use for symbolic designation in your current scope. If you only wish to import a subset of a system's variables, you only list those you wish to import when using `@unpack`. Note that this will, like when using e.g. `@variables` overwrite any variable you currently have declared sharing a name with one of the imported ones.

### [Via systems symbolics-based designation](@id advanced_intro_to_catalyst_designations_sys_symb)
Previously we described [hierarchical/compositional modelling](@ref advanced_intro_to_catalyst_modelling_approaches_hierarchical). This allowed us to e.g. describe a single species, but which can exist in separate models (e.g. an protein $X$ which max exist both in the nucleus and the cytosol). Here, it is not possible to know whether `X` refers to $X$ in the nucleus or in the cytosol. To solve this, we can use *indirect symbolic representation*. Using this requires first designating a `ReactionSystem` as [*complete*](@ref advanced_intro_to_catalyst_completeness).
```@example advanced_intro_to_catalyst_designations
rs = @reaction_network begin
    (k1, k2), X1 <--> X2
end
rs = complete(rs)
nothing # hide
```
Next, we can designate a variable through `rs.X1` (first the system to which it belongs, next the variable we wish to designate). Using `unpack`, we can check that the retrieved variable is the same one as stored within the system:
```@example advanced_intro_to_catalyst_designations
X1 = @unpack X1 = rs
isequal(rs.X1, X1)
```
Here, indirect symbolic representation can be used exactly like symbolic representation to perform a simulation:
```@example advanced_intro_to_catalyst_designations
u0_inds = [rs.X => 1.0]
ps_inds = [rs.p => 1.0, rs.d => 0.2]
oprob_inds = ODEProblem(rs, u0_inds, tspan, ps_inds)
sol_inds = solve(oprob_inds, Tsit5())
sol_inds == sol_symb == sol_syms
```

Using indirect symbolic representation is required when using compositional modelling. E.g. let us consider our [nucleus/cytoplasm model from previously](@ref advanced_intro_to_catalyst_modelling_approaches_hierarchical):
```@example advanced_intro_to_catalyst_designations
nucleus = @reaction_network begin
    @incomplete
    p, 0 --> X
    d_n, X --> 0
end
cytoplasm = @reaction_network begin
    @incomplete
    d_c, X --> 0
end
nc_model = @reaction_network begin
    @incomplete
    e, $(nucleus.X) --> $(nucleus.X)
end
@named nc_model = compose(nc_model, [nucleus, cytoplasm])
nc_model = complete(nc_model)
```
We can see that this actually have two different $X$ species:
```@example advanced_intro_to_catalyst_designations
species(nc_model)
```
Here, if we try to access `X` directly in the top-level model, it fails (as it do not know which $X$ to access):
```@example advanced_intro_to_catalyst_designations
nc_model.X
```
Instead, we have to designate the exact submodel which $X$ we wish to designate:
```@example advanced_intro_to_catalyst_designations
nc_model.nucleus.X
```
If we wish to simulate our model, we have to do this for all components:
```@example advanced_intro_to_catalyst_designations
u0 = [nc_model.nucleus.X => 1.0, nc_model.cytoplasm.X => 0.0]
ps = [nc_model.nucleus.p => 1.0, nc_model.nucleus.X => 0.5, nc_model.cytoplasm.X => 0.2, nc_model.e => 0.5]
oprob_inds = ODEProblem(nc_model, u0, tspan, ps)
sol = solve(nc_model, Tsit5())
plot(sol)
```

Here, we showed how indirect symbolic indexing is required for compositional modelling (although alternatives exist to get around it in many cases, e.g. designating [default values](@ref ref) for all components). It has an alternative use, *to form algebraic expressions without using `@unpack`*. Previously, we showed how, when we declared out model using the DSL, [we could use `@unpack` to form an algebraic expression of components](@ref advanced_intro_to_catalyst_designations_unpack). Alternatively, we can access $X1$ and $X2$ using indirect symbolic representation:
```@example advanced_intro_to_catalyst_designations
rs = @reaction_network begin
    (k1, k2), X1 <--> X2
end
rs = complete(rs)

u0 = [:X => 1.0]
tspan = (0.0, 10.0)
ps = [:p => 1.0, :d => 0.2]
oprob = ODEProblem(rs, u0, tspan, ps)
sol = solve(oprob, Tsit5())
sol[rs.X1 + rs.X2]
```

## Other low-level features
Finally, we describe some additional features of Catalyst. These are features that can be skipped, but for which some additional understanding of can be useful. 

### [Interpolating variables into the DSL](@id advanced_intro_to_catalyst_interpolation)

The DSL allows Julia variables to be interpolated for the network name, within
rate constant expressions, or for species/stoichiometry within reactions. Using
the lower-level symbolic interface we can then define symbolic variables and
parameters outside of the macro, which can then be used within expressions in
the DSL (see the [Programmatic Construction of Symbolic Reaction Systems](@ref programmatic_CRN_construction)
tutorial for details on the lower-level symbolic interface). For example,
```@example advanced_intro_to_catalyst_interpolation
@parameters k α
@variables t
@species A(t)
spec = A
par = α
rate = k*A
name = :network
rn = @reaction_network $name begin
    $rate*B, 2*$spec + $par*B --> $spec + C
  end
```
As the parameters `k` and `α` were pre-defined and appeared via interpolation,
we did not need to declare them within the `@reaction_network` macro,
i.e. they are automatically detected as parameters:
```@example tut2
parameters(rn)
```
as are the species coming from interpolated variables
```@example tut2
species(rn)
```

!!! note
    When using interpolation, expressions like `2$spec` won't work; the
    multiplication symbol must be explicitly included like `2*$spec`.

### [System completeness](@id advanced_intro_to_catalyst_completeness)
Both Catalyst and ModelingToolkit creates various *systems*, which then are used as input to *problems*. Especially ModelingToolkit (but [sometimes also Catalyst](@ref ref)) constructs models by composing systems together. Here, *completeness* is a termed used to designate that a model is ready for analysis. The important part here is that:
- Complete models cannot be composed with other systems.
- Incomplete models should not be used for for analysis.

Catalyst models are by default *complete* when they are created. This means that they can be provided to e.g. a `ODESystem` directly. This, however, means that they cannot be used for composition. To do so, the component model must be designated as *incomplete* when it is created (described [here](@ref ref)). Once a model have been composed of incomplete blocks, it can be marked as complete by applying the `complete` function. It is worth noting that ModelingToolkit (which is primarily designed for compositional modelling) works the other way. Models are, by default, incomplete, and must always be designated as complete before they are used.


### [Independent variables](@id advanced_intro_to_catalyst_indep_vars)
A system consists of *dependent* and *independent* variables. Here, dependent variables are the typical variables of the system, while the independent are what they depend on. Most (CRN) systems only have a single independent variable, *time*. Next, its dependent variables (typically species) depend on (are function of) this independent variables (that is why we declare them as e.g. `X(t)` in programmatic modelling). Non-time independent variables are typically spatial variables. Here, a species's value could depend both on time and its spatial position in teh system (e.g. if you model the concentration of a chemical in a reactor with spatial patterns).

It is possible to designate non-time independent variables for both dsl-based and programmatic modelling:
```@example advanced_intro_to_catalyst_indep_vars
using Catalyst

# In the DSL.
@reaction_network begin
    @ivs s
    @species X(t,s)
    (p,d), 0 <--> X
end

# Programmatically.
@variables s
@species X(t,s)
@parameters p d
rxs = [
    Reaction(p, nothing, [X]),
    Reaction(d, [X], nothing)
]
ReactionSystem(rxs, t; sivs=s)
```
A `ReactionSystem`'s primary (time) independent variable can be retrieved using `get_iv(rsys)`, and its spatial ones using `get_sivs(rsys)`.

Currently, Catalyst primarily supports spatial modelling on discrete spatial structures, where spatial independent variables are not used. Hence, this feature is not required for any current Catalyst feature.

### [Variable and reaction metadata](@id advanced_intro_to_catalyst_metadata)

#### [Variable metadata](@id advanced_intro_to_catalyst_metadata_variables)
It is possible  to attach metadata to symbolic variables. The use of this depend on the metadata, and some may simply be a convenience for the user. E.g. the most general metadata, `description` allows you to attach an arbitrary string to a variable. Metadata is attached when a variable is created by adding it within square brackets after the variable in question.
```@example advanced_intro_to_catalyst_indep_metadata
using Catalyst
t = ModelingToolkit.time_iv
@species X(t) [description="My protein of interest"]
nothing # hide
```
When several variables are created, you only put square bracket with metadata after those variables you are interested in
```@example advanced_intro_to_catalyst_indep_metadata
@parameters p [description="The production rate of X"] d
nothing # hide
```
Metadata cannot be updated, however, it is possible to override variables with new ones with updated metadataL
```@example advanced_intro_to_catalyst_indep_metadata
@species X(t) [description="My protein of very high interest"]
nothing # hide
```
Each metadata have its own implementation, and it is not possible to use custom metadata (unless they implement it themselves). Each metadata also have a specific getter. I.e. to retieve the description we use
```@example advanced_intro_to_catalyst_indep_metadata
getdescription(X)
```

#### [Reaction metadata](@id advanced_intro_to_catalyst_metadata_reactions_)

Reactions also have metadata. When created using the DSL, a reaction's metadata is given after the reactions (separated from it by a comma):
```@example advanced_intro_to_catalyst_indep_metadata
@reaction_network begin
    p, ∅ --> X, [description="Production of X"]
    d, X --> ∅
end
```
Here, it is possible to attach metadata to [bundled reactions](@ref ref). Here, we add parentheses around the metadata blocks just as we do for the rates:
```@example advanced_intro_to_catalyst_indep_metadata
rs = @reaction_network begin
    (p,d), ∅ <--> X, ([description="Production of X"], [description="Degradation of X"])
    d, X --> ∅
end
nothing # hide
```

When using programmatic modelling, metadata is added as an optional argument
```@example advanced_intro_to_catalyst_indep_metadata
t = ModelinGToolkit.time_iv
@species X(t)
@parameters p d
rxs = [
    Reaction(p, nothing, [X]; metadata = "Production of X")
    Reaction(d, [X], nothing; metadata = "Degradation of X")
]
rs = ReactionSystem(rxs, t)
nothing # hide
```

Reaction metadata can be accessed using the `hasmetadata`, `getmetadata`, and `getallmetadata` functions, which are described [here](@ref ref). Unlike for variables, users can specify and use arbitrary metadata for their reactions.