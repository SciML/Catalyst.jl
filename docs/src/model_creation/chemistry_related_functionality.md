# [Chemistry-related functionality](@id chemistry_functionality)

While Catalyst has primarily been designed around the modelling of biological systems, reaction network models are also common across chemistry. This section describes two types of functionality, that while of general interest, should be especially useful in the modelling of chemical systems.
- The `@compound` option, which enables the user to designate that a specific species is composed of certain subspecies.
- The `balance_reaction` function, which enables the user to balance a reaction so the same number of components occur on both sides.

## Modelling with compound species

### Creating compound species programmatically
We will first show how to create compound species through [programmatic model construction](@ref programmatic_CRN_construction), and then demonstrate using the DSL. To create a compound species, use the `@compound` macro, first designating the compound, followed by its components (and their stoichiometries). In this example, we will create a CO₂ molecule, consisting of one C atom and two O atoms. First, we create species corresponding to the components:
```@example chem1
using Catalyst
t = default_t()
@species C(t) O(t) 
```
Next, we create the `CO2` compound species:
```@example chem1
@compound CO2 ~ C + 2O
```
Here, the compound is the first argument to the macro, followed by its component (with the left-hand and right-hand sides separated by a `~` sign). While non-compound species (such as `C` and `O`) have their independent variable (in this case `t`) designated, independent variables are generally not designated for compounds (these are instead directly inferred from their components). Components with non-unit stoichiometries have this value written before the component (generally, the rules for designating the components of a compound are identical to those of designating the substrates or products of a reaction). The created compound, `CO2`, is also a species, and can be used wherever e.g. `C` can be used:
```@example chem1
isspecies(CO2)
```
In its metadata, however, is stored information of its components, which can be retrieved using the `components` (returning a vector of its component species) and `coefficients` (returning a vector with each component's stoichiometry) functions:
```@example chem1
components(CO2)
```
```@example chem1
coefficients(CO2)
```
Alternatively, we can retrieve the components and their stoichiometric coefficients as a single vector using the `component_coefficients` function:
```@example chem1
component_coefficients(CO2)
```
Finally, it is possible to check whether a species is a compound or not using the `iscompound` function:
```@example chem1
iscompound(CO2)
```

Compound components that are also compounds are allowed, e.g. we can create a carbonic acid compound (H₂CO₃) that consists of CO₂ and H₂O:
```@example chem1
@species H(t)
@compound H2O ~ 2H + O
@compound H2CO3 ~ CO2 + H2O
```

When multiple compounds are created, they can be created simultaneously using the `@compounds` macro, e.g. the previous code-block can be re-written as:
```@example chem1
@species H(t)
@compounds begin
    H2O ~ 2H + O
    H2CO3 ~ CO2 + H2O
end
```

### Creating compound species within the DSL
It is also possible to declare species as compound species within the `@reaction_network` DSL, using the `@compounds` options:
```@example chem1
rn = @reaction_network begin
    @species C(t) H(t) O(t)
    @compounds begin
        C2O ~ C + 2O
        H2O ~ 2H + O
        H2CO3 ~ CO2 + H2O
    end
    (k1,k2), H2O + CO2 <--> H2CO3
end
```
When creating compound species using the DSL, it is important to note that *every component must be known to the system as a species, either by being declared using the `@species` or `@compound` options, or by appearing in a reaction*. E.g. the following is not valid
```julia 
rn = @reaction_network begin``
    @compounds begin
        C2O ~ C + 2O
        H2O ~ 2H + O
        H2CO3 ~ CO2 + H2O
    end
    (k1,k2), H2O+ CO2 <--> H2CO3
end
```
as the components `C`, `H`, and `O` are not declared as a species anywhere. Please also note that only `@compounds` can be used as an option in the DSL, not `@compound`.

### Designating metadata and default values for compounds
Just like for normal species, it is possible to designate metadata and default values for compounds. Metadata is provided after the compound name, but separated from it by a `,`:
```@example chem1
@compound (CO2, [unit="mol"]) ~ C + 2O
```
Default values are designated using `=`, and provided directly after the compound name.:
```@example chem1
@compound (CO2 = 2.0) ~ C + 2O
```
If both default values and meta data are provided, the metadata is provided after the default value:
```@example chem1
@compound (CO2 = 2.0, [unit="mol"]) ~ C + 2O
```
In all of these cases, the side to the left of the `~` must be enclosed within `()`.

### Compounds with multiple independent variables
While we generally do not need to specify independent variables for compound, if the components (together) have more than one independent variable, this have to be done:
```@example chem1
t = default_t()
@variables s
@species N(s) O(t) 
@compound NO2(t,s) ~ N + 2O
```
Here, `NO2` depend both on a spatial independent variable (`s`) and a time one (`t`). This is required since, while multiple independent variables can be inferred, their internal order cannot (and must hence be provided by the user).

## Balancing chemical reactions
One use of defining a species as a compound is that they can be used to balance reactions so that the number of components are the same on both sides. Catalyst provides the `balance_reaction` function, which takes a reaction, and returns a balanced version. E.g. let us consider a reaction when carbon dioxide is formed from carbon and oxide `C + O --> CO2`. Here, `balance_reaction` enables us to find coefficients creating a balanced reaction (in this case, where the number of carbon and oxygen atoms are the same on both sides). To demonstrate, we first created the unbalanced reactions:
```@example chem1
rx = @reaction k, C + O --> $CO2
```
Here, the reaction rate (`k`) is not involved in the reaction balancing. We use interpolation for `CO2`, ensuring that the `CO2` used in the reaction is the same one we previously defined as a compound of `C` and `O`. Next, we call the `balance_reaction` function
```@example chem1
balance_reaction(rx)
```
which correctly finds the (rather trivial) solution `C + 2O --> CO2`. Here we note that `balance_reaction` actually returns a vector. The reason is that the reaction balancing problem may have several solutions. Typically, there is only a single solution (in which case this is the vector's only element). No, or an infinite number of, solutions is also possible depending on the given reaction.

Let us consider a more elaborate example, the reaction between ammonia (NH₃) and oxygen (O₂) to form nitrogen monoxide (NO) and water (H₂O). Let us first create the components and the unbalanced reaction:
```@example chem2
using Catalyst # hide
t = default_t()
@species N(t) H(t) O(t) 
@compounds begin
    NH3 ~ N + 3H
    O2 ~ 2O
    NO ~ N + O
    H2O ~ 2H + O
end
unbalanced_reaction = @reaction k, $NH3 + $O2 --> $NO + $H2O
```
We can now create a balanced version (where the amount of H, N, and O is the same on both sides):
```@example chem2
balanced_reaction = balance_reaction(unbalanced_reaction)[1]
```

Reactions declared as a part of a `ReactionSystem` (e.g. using the DSL) can be retrieved for balancing using the [`reactions`](@ref) function. Please note that balancing these will not mutate the `ReactionSystem`, but a new reaction system will need to be created using the balanced reactions.

!!! note
    Reaction balancing is currently not supported for reactions involving compounds of compounds.

### Balancing full systems
It is possible to balance all the reactions of a reaction system simultaneously using the `balance_system` function. Here, the output is a new system, where all reactions are balanced. E.g. We can use it to balance this system of Methane formation/combustion:
```@example chem2
rs = @reaction_network begin
    @species C(t) O(t) H(t)
    @compounds begin
        H2(t) = 2H
        CH4(t) = C + 4H
        O2(t) = 2O
        CO2(t) = C + 2O
        H2O(t) = 2H + O
    end
    1.0, C + H2 --> CH4
    2.0, CH4 + O2 --> CO2 + H2O
end
rs_balanced = balance_system(rs)
```
Except for the modified reaction stoichiometries, the new system is identical to the previous one.