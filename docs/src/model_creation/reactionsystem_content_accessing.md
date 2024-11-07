# [Accessing model properties](@id model_accessing)
Catalyst is based around the creation, analysis, and simulation of chemical reaction network models. Catalyst stores these models in [`ReactionSystem`](@ref) structures. This page describes some basic functions for accessing the content of these structures. This includes retrieving lists of species, parameters, or reactions that a model consists of. An extensive list of relevant functions for working with `ReactionSystem` models can be found in Catalyst's [API](@ref api).

!!! warning
    Generally, a field of a Julia structure can be accessed through `struct.fieldname`. E.g. a simulation's time vector can be retrieved using `simulation.t`. While Catalyst `ReactionSystem`s are structures, one should *never* access their fields using this approach, but rather using the accessor functions described below and in the [API](@ref api_accessor_functions) (direct accessing of fields can yield unexpected behaviours). E.g. to retrieve the species of a `ReactionsSystem` `rs`, use `Catalyst.get_species(rs)`, *not* `rs.species`. The reason is that, as shown [below](@ref model_accessing_symbolic_variables), Catalyst (and more generally any [ModelingToolkit]https://github.com/SciML/ModelingToolkit.jl system types) reserves this type of accessing for accessing symbolic variables stored in the system. I.e. `rs.X` refers to the `X` symbolic variable, not a field in `rs` named "X".

## [Direct accessing of symbolic model parameter and species](@id model_accessing_symbolic_variables)
Previously we have described how the parameters and species that Catalyst models contain are represented using so-called [*symbolic variables*](@ref introduction_to_catalyst) (and how these enable the forming of [*symbolic expressions*](@ref introduction_to_catalyst)). We have described how, during [programmatic modelling](@ref programmatic_CRN_construction), the user has [direct access to these](@ref programmatic_CRN_construction) and how this can be [taken advantage of](@ref programmatic_CRN_construction). We have also described how, during [DSL-based modelling](@ref dsl_description), the need for symbolic representation can be circumvented by [using `@unpack`](@ref dsl_advanced_options_symbolics_and_DSL_unpack) or by [creating an observable](@ref dsl_advanced_options_observables). However, sometimes, it is easier to *directly access a symbolic variable through the model itself*, something which we will describe here.

Let us consider the following [two-state model](@ref basic_CRN_library_two_states)
```@example model_accessing_symbolic_variables
using Catalyst
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
```
If we wish to access one of the symbolic variables stored in it (here `X1`, `X2`, `k1`, and `k2`), we simply write 
```@example model_accessing_symbolic_variables
rs.X1
```
to access e.g. `X1`. This symbolic variable can be used just like those [declared using `@parameters` and `@species`](@ref programmatic_CRN_construction):
```@example model_accessing_symbolic_variables
using OrdinaryDiffEq
u0 = [rs.X1 => 1.0, rs.X2 => 2.0]
ps = [rs.k1 => 2.0, rs.k2 => 4.0]
oprob = ODEProblem(rs, u0, (0.0, 10.0), ps)
sol = solve(oprob)
nothing # hide
```
We can also use them to form symbolic expressions:
```@example model_accessing_symbolic_variables
Xtot = rs.X1 + rs.X2
```
which can be used when we e.g. [plot our simulation](@ref simulation_plotting_options):
```@example model_accessing_symbolic_variables
using Plots
plot(sol; idxs = [rs.X1, rs.X2, Xtot])
```

Next we create our two-state model [programmatically](@ref programmatic_CRN_construction):
```@example model_accessing_symbolic_variables
t = default_t()
@species X1(t) X2(t)
@parameters k1 k2
rxs = [
    Reaction(k1, [X1], [X2]),
    Reaction(k2, [X2], [X1])
]
@named rs_prog = ReactionSystem(rxs, t)
rs_prog = complete(rs_prog)
nothing # hide
```
Here, we can confirm that the symbolic variables we access through the model are identical to those we used to create it:
```@example model_accessing_symbolic_variables
isequal(rs.k1, k1)
```

!!! warning
    When accessing model symbolic variables through the model (using e.g. `rs.X1`), it is important to first ensure that the [*model has been marked complete*](@ref programmatic_CRN_construction).

## [Accessing basic model properties](@id model_accessing_basics)

### [Accessing model parameter and species](@id model_accessing_basics_parameters_n_species)
Previously we showed how to access individual parameters or species of a `ReactionSystem` model. Next, the `parameters` and [`species`](@ref) functions allow us to retrieve *all* parameters and species as vectors:
```@example model_accessing_basics
using Catalyst # hide
sir = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end
parameters(sir)
```
```@example model_accessing_basics
species(sir)
```
These vectors contain the exact same symbolic variables that we would access through the system:
```@example model_accessing_basics
issetequal([sir.S, sir.I, sir.R], species(sir))
```

If we wish to count the number of parameters or species in a system, we can do this directly through the [`numparams`](@ref) and [`numspecies`](@ref) functions:
```@example model_accessing_basics
numparams(sir)
```
```@example model_accessing_basics
numspecies(sir)
```

### [Accessing model reactions](@id model_accessing_basics_reactions)
A vector containing all a model's [reactions](@ref programmatic_CRN_construction) can be retrieved using the [`reactions`](@ref) function:
```@example model_accessing_basics
reactions(sir)
```
We can count the number of reactions in a model using the [`numreactions`](@ref) function:
```@example model_accessing_basics
numreactions(sir)
```
Finally, a vector with all the reactions' rates can be retrieved using [`reactionrates`](@ref):
```@example model_accessing_basics
reactionrates(sir)
```

### [Accessing content of models coupled to equations](@id model_accessing_basics_reactions)
Previously, we have shown how to [couple equations to a chemical reaction network model](@ref constraint_equations_coupling_constraints), creating models containing [non-species unknowns (variables)](@ref constraint_equations_coupling_constraints). Here we create a birth-death model where some nutrient supply (modelled through the variable $N$) is depleted in the presence of $X$:
```@example model_accessing_basics
using Catalyst # hide
coupled_crn = @reaction_network begin
    @equations D(N) ~ -N * X
    (p/(1+N),d), 0 <--> X
end
```
Here, the `unknowns` function returns all unknowns (i.e. species and variables):
```@example model_accessing_basics
unknowns(coupled_crn)
```
Meanwhile, `species` returns the species only, while [`nonspecies`](@ref) returns the variables only:
```@example model_accessing_basics
species(coupled_crn)
```
```@example model_accessing_basics
nonspecies(coupled_crn)
```

Similarly, the `equations` function returns a vector with all reactions and equations of the model (ordered so that reactions occur first and equations thereafter):
```@example model_accessing_basics
equations(coupled_crn)
```
Meanwhile, [`reactions`](@ref) returns the reactions only, while [`nonreactions`](@ref) returns any algebraic or differential equations:
```@example model_accessing_basics
reactions(coupled_crn)
```
```@example model_accessing_basics
nonreactions(coupled_crn)
```

### [Accessing other model properties](@id model_accessing_basics_others)
There exist several other functions for accessing model properties. 

The `observed`, `continuous_events`, `discrete_events` functions can be used to access a model's [observables](@ref dsl_advanced_options_observables), [continuous events](@ref constraint_equations_events), and [discrete events](@ref constraint_equations_events), respectively.

The `ModelingToolkit.get_iv` function can be used to retrieve a [model's independent variable](@ref programmatic_CRN_construction):
```@example model_accessing_basics
ModelingToolkit.get_iv(sir)
```

## [Accessing properties of hierarchical models](@id model_accessing_hierarchical)
Previously, we have described how [compositional modelling can be used to create hierarchical models](@ref compositional_modeling). There are some special considerations when accessing the content of hierarchical models, which we will describe here.

First, we will create a simple hierarchical model. It describes a protein ($X$) which is created in its inactive form ($Xᵢ$) in the nucleus, from which it is transported to the cytoplasm, where it is activated.
```@example model_accessing_hierarchical
using Catalyst # hide
# Declare submodels.
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
This model consists of a top-level system, which contains the transportation reaction only, and two subsystems. We can retrieve all the subsystems of the top-level system through `Catalyst.get_systems`:
```@example model_accessing_hierarchical
Catalyst.get_systems(rs)
nothing # hide
```
!!! note
    If either of the subsystems had had further subsystems, these would *not* be retrieved by `Catalyst.get_systems` (which only returns the direct subsystems of the input system).

### [Accessing parameter and species of hierarchical models](@id model_accessing_hierarchical_symbolic_variables)
Our hierarchical model consists of a top-level system (`rs`) with two subsystems (`nucleus_sys` and `cytoplasm_sys`). Note that we have given our subsystems [names](@ref dsl_advanced_options_naming) (`nucleus` and `cytoplasm`, respectively). Above, we retrieved the subsystems by calling `Catalyst.get_systems` on our top-level system. We can also retrieve a subsystem directly by calling:
```@example model_accessing_hierarchical
rs.nucleus
```
```@example model_accessing_hierarchical
rs.cytoplasm
```
!!! note
    When accessing subsystems, we use the subsystems' [names](@ref dsl_advanced_options_naming), *not* the name of their variables (i.e. we call `rs.nucleus`, not `rs.nucleus_sys`).

Next, if we wish to access a species declared as a part of one of the subsystems, we do so through it. E.g. here we access `Xₐ` (which is part of the cytoplasm subsystem):
```@example model_accessing_hierarchical
rs.cytoplasm.Xₐ
```
Note that species contained in a subsystem have the subsystem's name prepended to their name when we access it.

Both subsystems contain a species `Xᵢ`. However, while they have the same name, *these are different species when accessed through their respective models*:
```@example model_accessing_hierarchical
isequal(rs.nucleus.Xᵢ, rs.cytoplasm.Xᵢ)
```

The same holds for the model parameters, i.e. while each subsystem contains a parameter `d`, these are considered different parameters:
```@example model_accessing_hierarchical
isequal(rs.nucleus.d, rs.cytoplasm.d)
```
The parameter `kₜ` is actually contained within the top-level model, and is accessed directly through it:
```@example model_accessing_hierarchical
rs.kₜ
```

Practically, when we simulate our hierarchical model, we use all of this to designate initial conditions and parameters. I.e. below we designate values for the two `Xᵢ` species and `d` parameters separately:
```@example model_accessing_hierarchical
using OrdinaryDiffEq, Plots
u0 = [rs.nucleus.Xᵢ => 0.0, rs.cytoplasm.Xᵢ => 0.0, rs.cytoplasm.Xₐ => 0.0]
ps = [rs.nucleus.p => 1.0, rs.nucleus.d => 0.2, rs.cytoplasm.kₐ => 5.0, rs.cytoplasm.d => 0.2, rs.kₜ => 0.1]
oprob = ODEProblem(rs, u0, (0.0, 10.0), ps)
sol = solve(oprob)
plot(sol)
```

!!! note
    When we access a symbolic variable through a subsystem (e.g. `rs.nucleus.Xᵢ`) that subsystem's name is prepended to the symbolic variable's name (we call this *namespacing*). This is also the case if we access it through the original model, i.e. `nucleus_sys.Xᵢ`. Namespacing is only performed when we access variables of [*incomplete systems*](@ref programmatic_CRN_construction). I.e. `isequal(nucleus_sys.d, cytoplasm_sys.d)` returns false (as the systems are incomplete and namespacing is performed). However, `isequal(complete(nucleus_sys).d, complete(cytoplasm_sys).d)` returns true  (as the systems are complete and namespacing is not performed). This is the reason that the system top-level system's name is never prepended when we do e.g. `rs.kₜ` (because here, `rs` is complete).

### [Accessing the content of hierarchical models](@id model_accessing_hierarchical_symbolic_variables)
In the last section, we noted that our hierarchical model contained several instances of the `Xᵢ` species. The [`species`](@ref) function, which retrieves all of a model's species shows that our model has three species (two types of `Xᵢ`, and one type of `Xₐ`)
```@example model_accessing_hierarchical
species(rs)
```
Similarly, `parameters` retrieves five different parameters. Here, we note that `kₜ` (which has no model name prepended) belongs to the top-level system (and not a subsystem):
```@example model_accessing_hierarchical
parameters(rs)
```

If we wish to retrieve the species (or parameters) that are specifically contained in the top-level system (and not only indirectly through its subsystems), we can use the `Catalyst.get_species` (or `Catalyst.get_ps`) functions:
```@example model_accessing_hierarchical
Catalyst.get_species(rs)
```
```@example model_accessing_hierarchical
Catalyst.get_ps(rs)
```
Here, our top-level model contains a single parameter (`kₜ`), and two the two versions of the `Xᵢ` species. These are all the symbolic variables that occur in the transportation reaction (`@kₜ, $(nucleus_sys.Xᵢ) --> $(cytoplasm_sys.Xᵢ)`), which is the only reaction of the top-level system. We can apply these functions to the systems as well. However, when we do so, the systems' names are not prepended:
```@example model_accessing_hierarchical
Catalyst.get_ps(rs.nucleus)
```

Generally, functions starting with `get_` retrieve only the components stored in the input system (and do not consider its subsystems), while the corresponding function without `get_` also retrieves the components stored in subsystems. I.e. compare the `Catalyst.get_rxs` and [`reactions`](@ref) functions:
```@example model_accessing_hierarchical
reactions(rs)
```
```@example model_accessing_hierarchical
Catalyst.get_rxs(rs)
```
Other examples of similar pairs of functions are `get_unknowns` and `unknowns`, and `get_observed` and `observed`.

!!! note
    Due to the large number of available accessor functions, most `get_` functions are not directly exported by Catalyst. Hence, these must be used as `Catalyst.get_rxs`, rather than `get_rxs` directly. Again, a full list of accessor functions and instructions on how to use them can be found in [Catalyst's API](@ref api).