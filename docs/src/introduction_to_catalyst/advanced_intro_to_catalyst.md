# [How Catalyst works and advanced introduction](@id advanced_intro_to_catalyst) 
We have [previously](@ref introduction_to_catalyst) introduced Catalyst, and how to use it to create and simulate simple chemical reaction network (CRN) models. Here, we provide a more in-depth introduction. We will explain how Catalyst is constructed (i.e. what packages it builds on and how), provide a comprehensive overview of the different ways Catalyst models can be created, and describe how these are different.

This advanced introduction is not required for understanding the remaining documentation, nor for using any of Catalyst's features. It is not expected that new users will read it. Instead, it aims to provide additional context and understanding for more experienced users. If you have just started exploring how to use Catalyst, reading this section is not needed. However, if you are using (or plan to use) Catalyst extensively, reading this page will be useful.

## [How Catalyst is constructed](@id advanced_intro_to_catalyst_background)
Catalyst is built on top of the domain-agnostic modelling package [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl). ModelingToolkit is in turn built on the computer algebraic system [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl). Here we will introduce these packages, and how the feed into Catalyst.

### [The Symbolics.jl computer algebraic system](@id advanced_intro_to_catalyst_background_symbolics)
Symbolics is a Julia package for representing, and doing computations on, symbolic expressions. It allow us to declare symbolic variables that can be used to construct equations. Primarily this is done through the `@variables` macro:
```@example advanced_intro_to_catalyst_symbolics
using Symbolics
@variables x y
nothing # hide
```
Here, both `x` and `y` becomes Julia variables in the current scope, e.g. we can do
```@example advanced_intro_to_catalyst_symbolics
x
```
to confirm that `x` is a Julia variable we can work with. Next, we can use these variables to create algebraic equations:
```@example advanced_intro_to_catalyst_symbolics
eq = x^2 - y^2
```
Symbolics's primary feature is that it enables algebraic manipulation on algebraic expressions. We can use its `simplify` function to get a simpler (but equivalent) expression. E.g. it can correctly apply the [difference of two squares](https://en.wikipedia.org/wiki/Difference_of_two_squares) rule to our previous equation:
```@example advanced_intro_to_catalyst_symbolics
simplify(eq / (x + y))
```

!!! warning
    `@variables` (and similar [macros introduced later on](@ref ref), like `@parameters` and `@species`) create new variables in the current scope. E.g. if you have previously defined `x = 5` and then run `@variables x`, this will overwrite the previous value of `x`.

Symbolics is used internally within Catalyst to represent all models. One advantage of this is that complicated algebraic expressions can be automatically simplified. E.g. if we create a model with a (trivially) unnecessarily complicated rate and check what it is, we note that it has been automatically simplified
```@example advanced_intro_to_catalyst_symbolics
using Catalyst
rn = @reaction_network begin
    (k1*k2)/k2, X --> 0
end
reactions(rn)[1].rate
```

Another feature of Symbolics is that it can differentiate equations:
```@example advanced_intro_to_catalyst_symbolics
Symbolics.derivative(eq, x)
```
This is used internally to symbolically compute system Jacobians. As explained in a [later section](@ref ode_simulation_performance_symbolic_jacobian) this can be used to improve simulation performance.


### [The ModelingToolkit.jl intermediary representation and modelling package](@id advanced_intro_to_catalyst_background_mtk)
ModelingToolkit.jl is a package with two separate functionalities:
- It provides an interface for building mathematical models from symbolic equations.
- It implements an intermediary representation for representing such mathematical models.

We will start with describing the first point, by showing how ModelingToolkit can be used to build a ODE model. First we need to pre-declare the model's variables and parameters. This is done in a similar manner as we defined variables using Symbolics, however, with two differences:
- ODEs (typically) describe the evolution of its variables with time. The variables are thus *time dependent* and must be declared as functions of the time [*independent variable*](@ref advanced_intro_to_catalyst_independent_vars), $t$.
- Parameters must be declared using the `@parameters` macro (to distinguish them from the variables).

Here, we import `t` (the time independent variable) and `D` (the differential with respect to `t`), declares the variable `X` and the parameters `p` and `d`:
```@example advanced_intro_to_catalyst_mtk
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
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
\frac{dX}{dt} = p - d*X
```
Here, the differential of anything with respect to $t$ is defined using the `D` function, i.e. we can form our differential equation as
```@example advanced_intro_to_catalyst_mtk
diff_eq = D(X) ~ p - d*X
```
!!! note
    When we form symbolic equations (e.g. here setting `D(X)` equal to `p - d*X`), `~` is used as the equivalence operator. Next, `=` is used in the normal Julia fashion to assign the equation `D(X) ~ p - d*X` to the Julia variable `diff_eq`.

Finally, to use our differential equation for any purpose (like simulations), we must first encapsulate it within an `ODESystem`. The `ODESystem` constructor's input is a vector listing all equations, and also the time independent variable:
```@example advanced_intro_to_catalyst_mtk
@mtkbuild osys = ODESystem([diff_eq], t)
```
Here, declarations of ModelingToolkit systems are prepended with an `@mtkbuild` statement. This have two effects:
- It sets a *system name* (which is equal to the variable it is stored in, in this case `osys`). System names are discussed in more detail [here](@ref ref).
- It applied the `structural_simplify` function to ths system (described in more detail [here](https://docs.sciml.ai/ModelingToolkit/stable/systems/NonlinearSystem/#Transformations)).

Finally, the `ODESystem` can used to construct an `ODEProblems`. This is done in a similar manner as for Catalyst models, however, we can use the actually variables and parameter themselves to designate their values in the initial condition and parameter values vectors:
```@example advanced_intro_to_catalyst_mtk
using OrdinaryDiffEq, Plots

u0 = [X => 1.0]
tspan = (0.0, 10.0)
ps = [p => 1.0, d => 0.2]

oprob = ODEProblem(osys, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```
[ModelingToolkit's documentation](https://docs.sciml.ai/ModelingToolkit/stable/) provide a more extensive description on how to create model's using it.

### [System types and their conversions](@id advanced_intro_to_catalyst_background_systems)
We have previously noted that Catalyst CRN models (generated by e.g. the `@reaction_network` macro) are stored in so-called `ReactionSystem` structures. These can (just like ModelingToolkit's `ODESystem`s) be used as input to `ODEProblem`s. Internally, what actually happens here is that Catalyst's `ReactionSystem`s are converted to ModelingToolkit `ODESystems`, which are used for the actual simulation. To demonstrate this we carry out the same process manually. First we use the `convert` function to convert our `ReactionSystem` to an `ODESystem`:
```@example advanced_intro_to_catalyst_mtk
using Catalyst
rsys = @reaction_network begin
    (p,d), 0 <--> X
end
osys_catalyst = convert(ODESystem, rsys)
osys_catalyst = complete(osys_catalyst)
```
Here we have implemented the same model (a [birth-death process](@ref basic_CRN_library_bd)) using both Catalyst and ModelingToolkit. Next, we simulate `osys_catalyst` in the exact same way as `osys` was simulated. 
```@example advanced_intro_to_catalyst_mtk
oprob_catalyst = ODEProblem(osys_catalyst, u0, tspan, ps)
sol_catalyst = solve(oprob_catalyst)
sol == sol_catalyst
```
Here we also confirm that the output simulations are identical.

ModelingToolkit implements a wide range of [system types](https://docs.sciml.ai/ModelingToolkit/stable/basics/AbstractSystem/), enabling us to model a wide range of systems. This include e.g.:
- [`ODESystem`s](https://docs.sciml.ai/ModelingToolkit/stable/systems/ODESystem/) to represent ODEs.
- [`SDESystem`s](https://docs.sciml.ai/ModelingToolkit/stable/systems/SDESystem/) to represent SDEs.
- [`JumpSystem`s](https://docs.sciml.ai/ModelingToolkit/stable/systems/JumpSystem/) to represent jump processes.
- [`NonlinearSystem`s](https://docs.sciml.ai/ModelingToolkit/stable/systems/NonlinearSystem/) to represent systems of nonlinear equations.

The advantage of representing a model as a CRN is that this representation can be unambiguously converted to ODEs (via the [reaction rate equation](@ref ref)), SDEs (via the [chemical Langevin equation](@ref ref)), and jump processes (via the [stochastic chemical kinetics](@ref ref)). In the same manner, the advantage of representing a model as a Catalyst `ReactionSystem` is that these can unambiguously be converted to `ODESystem`s, `SDESystem`s, `JumpSystem`s, and `NonlinearSystem`s (the last one is primarily used to find ODE steady states). I.e. all of these are possible:
```@example advanced_intro_to_catalyst_mtk
ssys = convert(SDESystem, rsys)
jsys = convert(JumpSystem, rsys)
nsys = convert(NonlinearSystem, rsys)
nothing # hide
```
Next, these can be passed their respective problem types, e.g.
```@example advanced_intro_to_catalyst_mtk
using StochasticDiffEq
ssys = complete(ssys)
sprob = SDEProblem(ssys, u0, tspan, ps)
ssol = solve(sprob, STrapezoid())
ssol = solve(sprob, STrapezoid(); seed = 12345) # hide
plot(ssol)
```

### [ModelingToolkit as an intermediary representation](@id advanced_intro_to_catalyst_background_ir)
Above we have indirectly demonstrated ModelingToolkit's second feature, that it provides an *intermediary representation* (IR) for representing mathematical models (as e.g. `ODESystem`s). Here, a large number of packages for building mathematical models in the Julia programming languages (e.g. Catalyst, [OrbitalTrajectories.jl](https://github.com/dpad/OrbitalTrajectories.jl), and [NBodySimulator.jl](https://github.com/SciML/NBodySimulator.jl)) builds their models using ModelingToolkit's IRs. Next, a large number of packages for model simulation and analysis (e.g. [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl), [StructuralIdentifiability.jl](https://github.com/SciML/StructuralIdentifiability.jl), and [PEtab.jl](https://github.com/sebapersson/PEtab.jl)) ensures that their packages also are compatible with these IRs.

This is the primary way Catalyst is able to support such a large number of features. By linking Catalyst to the ModelingToolkit IR, its models is automatically supported by a large ecosystem of packages implementing various modelling tools, while requiring only minimal implementation within Catalyst itself.


## [Types of Catalyst modelling approached](@id advanced_intro_to_catalyst_modelling_approaches)
There are three main ways of building Catalyst `ReactionSystem` models. The first two are distinct from each other, while the third one can build on top of either (or both) of these. Here we will briefly describe all three approached (with more details of each described in separate sections of this documentation).

### [DSL-based modelling](@id advanced_intro_to_catalyst_modelling_approaches_DSL)
Catalyst implements a *domain-specific language* (DSL) through its `@reaction_network` macro. This enables the writing of CRN models in their native chemical reaction form. here we use it to create a simple [Michaelis-Menten kinetics model](@ref basic_CRN_library_mm):
```@example advanced_intro_to_catalyst_modelling_approaches
using Catalyst
rs_dsl = @reaction_network rs begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
```
A throughout description of how to model CRNs using the DSL is provided [here](@ref dsl_description), with some of its more advanced options being described [here](@ref dsl_advanced_options). The DSL is the most frequent modelling approach used throughout this documentation. Its primary advantage is that it enables models to be created quickly, and using code which is highly human-interpretable.

### [Programmatic modelling](@id advanced_intro_to_catalyst_modelling_approaches_programmatic)
Previously we showed how to programmatically create `ODESystems` [using ModelingToolkit](@ref advanced_intro_to_catalyst_background_mtk). Catalyst have a similar interface for creating `ReactionSystem`s. Here, we first designate the species and parameters of our model (note that we use the `@species` macro instead of the `@variables` macro, which is similar, but specifically designate that the created [variables are *species*]).
```@example advanced_intro_to_catalyst_modelling_approaches
t = default_t()
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
Finally, we create a `ReactionSystem` by combing the reactions, and also designating the [independent variable](@ref advanced_intro_to_catalyst_independent_vars):
```@example advanced_intro_to_catalyst_modelling_approaches
rs_prog = ReactionSystem(rxs, t; name=:rs)
```
The resulting `ReactionSystem` is identical to the one we created directly using the DSL:
```@example advanced_intro_to_catalyst_modelling_approaches
rs_dsl == rs_prog
```

More details on programmatic model construction can be found [here](@ref programmatic_CRN_construction). While programmatic modelling typically is a bit more laborious than DSL-based, it also gives the user more control over the modelling process. Also, while DSL-based modelling requires the user to explicitly declare all reactions, programmatic modelling can be used to [automatically generate model reactions from simple rules](@ref programmatic_generative_linear_pathway).

!!! note
    In ModelingToolkit, the time independent variable is typically imported directly from the package through `using ModelingToolkit: t_nounits`. In Catalyst, the `default_t` importer function is instead used: `t = default_t()`.

### [Hierarchical and compositional modelling](@id advanced_intro_to_catalyst_modelling_approaches_hierarchical)
The final type of model construction is hierarchical (or compositional) modelling. This involves first creating smaller models, and then composing them together to create a single, more complex model. The component models can be created using either DSL-based or programmatic modelling.

Here, we create a simple model of a single protein ($X$) which is produced in the nucleus, from which it can be transported to the cytoplasm, where it is activated:
```@example advanced_intro_to_catalyst_modelling_approaches
# Declare sub models.
nucleus_sys = @network_component nucleus begin
    (p,d), 0 <--> Xᵢ
end
cytoplasm_sys = @network_component cytoplasm begin
    kₐ, Xᵢ --> Xₐ
    d, (Xᵢ, Xₐ) --> 0
end

# Assembly hierarchical model.
transport = @reaction kₜ, $(nucleus_sys.Xᵢ) --> $(cytoplasm_sys.Xᵢ)
@named rs = ReactionSystem([transport], default_t(); systems = [nucleus_sys, cytoplasm_sys])
rs = complete(rs)
```
This hierarchical model consists of a top-level model (which contains the transportation reaction only) which contains two sub-models (one describing the nucleus and one describing the cytoplasm).

Hierarchical modelling can be used to model system spread across several different compartments, or to create models from commonly occurring components. Hierarchical and compositional modelling is described in more detail [here](@ref compositional_modeling).


## [Designating species and parameters](@id advanced_intro_to_catalyst_designations)
Frequently when working with Catalyst, we need to represent a species or parameter (or variables, when [CRNs are coupled with equations](@ref ref)). Primarily this happens when we set the initial conditions and parameter values of a simulation, but it is also relevant when we e.g. [designate which species to plot](@ref simulation_plotting_options) or [updates a parameter during a simulation](@ref simulation_structure_interfacing_integrators). There are three different ways to represent species and parameters, each with their pros and cons. Here we will provide a short overview of each.

### [Symbol-based designation](@id advanced_intro_to_catalyst_designations_syms)
Throughout this documentation, we primarily represent parameters and species using symbols. Here, we represent $X$ through its corresponding [symbol](https://docs.julialang.org/en/v1/base/base/#Core.Symbol), `:X` is used (a symbol is created by pre-pending an expression with `:`). E.g. to simulate the following model:
```@example advanced_intro_to_catalyst_designations
using Catalyst
rs_dsl = @reaction_network begin
    (p,d), 0 <--> X
end
```
we can represent our initial conditions and parameters with their respective symbol forms:
```@example advanced_intro_to_catalyst_designations
using  OrdinaryDiffEq
u0_syms = [:X => 1.0]
tspan = (0.0, 10.0)
ps_syms = [:p => 1.0, :d => 0.2]
oprob_syms = ODEProblem(rs_dsl, u0_syms, tspan, ps_syms)
sol_syms = solve(oprob_syms)
nothing # hide
```
Next, to obtain a vector with the simulation's $X$ value we can again use $X$'s symbol representation:
```@example advanced_intro_to_catalyst_designations
sol_syms[:X]
```

There are two cases where symbol-based designations cannot be used:
- When using [hierarchical/compositional modelling](advanced_intro_to_catalyst_modelling_approaches_hierarchical). Here, if several sub-models contain species called `X`, it cannot be inferred which one `:X` refers to.
- When algebraic expression of species and/or parameters must be formed (as used e.g. [here](@ref ref)). Typically, these situations indicate that an [*observable*](@ref dsl_advanced_options_observables) should be created to designate the quantity of index.

Generally, symbol-based designates is used when [models are created using the DSL](@ref advanced_intro_to_catalyst_modelling_approaches_DSL).

### [Symbolics variables-based designation](@id advanced_intro_to_catalyst_designations_symbolics)
Previously, we showed how model creation using ModelingToolkit involves [explicitly declaring variables and parameters, which then can be used to designate their values](@ref advanced_intro_to_catalyst_background_mtk). We also showed how Catalyst `ReactionSystem`s can be [created programmatically, in a similar way to how ModelingToolkit models are created](@ref advanced_intro_to_catalyst_modelling_approaches_programmatic). Here, just like for ModelingToolkit, we can represent species and parameters using their respective symbolic variables directly:
```@example advanced_intro_to_catalyst_designations
t = default_t()
@species X(t) 
@parameters p d
rxs = [
    Reaction(p, nothing, [X]), 
    Reaction(d, [X], nothing)
]
rs_prog = ReactionSystem(rxs, t; name=:rs)
rs_prog = complete(rs_prog)

u0_symb = [X => 1.0]
ps_symb = [p => 1.0, d => 0.2]
oprob_symb = ODEProblem(rs_prog, u0_symb, tspan, ps_symb)
sol_symb = solve(oprob_symb)
nothing # hide
```
We can check that the resulting solution is the same as when we used symbol-based designation (and DSL-based model creation):
```@example advanced_intro_to_catalyst_designations
sol_syms == sol_symb
```
We can also use symbolic variables to e.g. index the output solution:
```@example advanced_intro_to_catalyst_designations
sol_symb[X]
```

The advantage with the Symbolics for is that we can form algebraic expressions from our symbolic variables. E.g. if we (for some reason) would like to evaluate the expression $2X + p$ across the solution vector, we can do this directly using the Symbolic form:
```@example advanced_intro_to_catalyst_designations
sol_symb[2X + p]
```
Indeed, we can confirm that this returns the expected output:
```@example advanced_intro_to_catalyst_designations
sol_symb[2X + p] == 2*sol[X] .+ 1.0
```

Symbolic variables-based representation is strictly superior to symbol-based representation. However, it is generally only used in combination with programmatic modelling (as during DSL-based modelling we never actually declare the symbolic variables, and thus cannot use them directly). To use symbolic variables-based representation in combination with the DSL, the [`@unpack` macro can be used](@ref ref).

### [System-based symbolic designation](@id advanced_intro_to_catalyst_designations_sys_symb)
However Catalyst models are created, their species and parameters are internally stored as symbolic variables. Here, it is possible to directly access these through the model itself. E.g. we can access species `X` through:
```@example advanced_intro_to_catalyst_designations
rs_prog.X
```
we can confirm that this is the same `X` which was used to create the model originally:
```@example advanced_intro_to_catalyst_designations
isequal(X, rs_prog.X)
```
Here, we can use the symbolic variables stored within the model to represent our initial conditions and parameter values:
```@example advanced_intro_to_catalyst_designations
u0_symb = [rs_prog.X => 1.0]
ps_symb = [rs_prog.p => 1.0, rs_prog.d => 0.2]
oprob_symb = ODEProblem(rs_prog, u0_symb, tspan, ps_symb)
sol_symb = solve(oprob_symb)
nothing # hide
```
or to index our solution:
```@example advanced_intro_to_catalyst_designations
sol_symb[rs_prog.X]
```

Since the models we created programmatically and through the DSL are identical, we can access the symbolic variables stored in the DSL-created model using this approach. This can be used to create algebraic expressions when DSL-based modelling is used (where we normally do not have access to the symbolic variables):
```@example advanced_intro_to_catalyst_designations
sol_syms[2rs_dsl.X + rs_dsl.p]
```

This approach for representing species and parameters is primarily used in combination with [hierarchical or compositional modelling](@ref compositional_modeling) (as it enables us to designate species that have identical names, but are stored in different subsystems). System-based symbolic accessing is discussed in more detail [here](@ref ref).

!!! warn
    When accessing species and parameters that are stored in a system using the above approach it is important to ensure that the model is *complete*. Model completeness is described [here](@ref ref).

## [Independent variables](@id advanced_intro_to_catalyst_independent_vars)
Above, we described how declared variables (so-called *dependent variables*, like `X(t)`) depend on a so-called *independent variable*. Here, dependent variables are the quantities which behaviour is described by our model, while the independent variables are variable these depend on. Most (CRN) systems only have a single independent variable, *time*. Next, its dependent variables (typically species) depend on (are function of) this independent variable (this is why we declare them as e.g. `X(t)` in programmatic modelling). Non-time independent variables exists, and these describe spatial dimensions. A species's value can depend on both time and one (or several) spatial dimension(s) (e.g. the concentration of a species `X` may depend on time and on its spatial location in a chemical reactor). They may also depend on a spatial independent variable only.

It is possible to designate non-time independent variables for both DSL-based and programmatic modelling:
```@example advanced_intro_to_catalyst_independent_vars
using Catalyst

# Using the DSL.
rs_dsl = @reaction_network begin
    @ivs t x
    @species X(t,x)
    (p,d), 0 <--> X
end

# Programmatically.
t = default_t()
@variables x
@species X(t,x)
@parameters p d
rxs = [
    Reaction(p, nothing, [X]),
    Reaction(d, [X], nothing)
]
@named rs_prog = ReactionSystem(rxs, t; spatial_ivs = [x])
nothing # hide
```
A `ReactionSystem`'s primary (time) independent variable can be retrieved using `Catalyst.get_iv(rsys)`, and its spatial ones using `Catalyst.get_sivs(rsys)`.

!!! note
    Above we uses Catalyst's `default_t` function to declare the time independent variable. No such function exists for spatial independent variables, as these have to be declare manually. Since these are just variables that do not depend on any other variables, they can simply be declared using e.g. `@variables x`.

Currently, Catalyst primarily supports spatial modelling on discrete spatial structures (where spatial independent variables are not used). Hence, the declaration of spatial independent variables are of limited use.
