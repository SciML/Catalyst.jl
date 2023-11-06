# [Structural Identifiability Analysis](@id structural_identifiability)
During parameter fitting, parameter values are inferred from data. Identifiability is a concept describing to what extent the identification of parameter values for a certain model is actually feasible. Ideally, parameter fitting should always be accompanied with an identifiability analysis of the problem. Identifiability can be divided into *structural* and *practical* identifiability[^1]. Structural identifiability considers only  the system and what quantities we can measure to determine which quantities can be identified. Practical identifiability instead considers the available data, and determines what system quantities can be infeed from it. Generally, in the hypothetical case of an infinite amount of without noise, practical identifiability becomes identical to structural identifiability. Generally, structural identifiability is assessed before parameters are fitted, while practical identifiability is assessed afterwards.

Structural identifiability can be illustrated in the following example network:
${dx \over dt} = p1*p2*x(t)$
where, however much data is collected on *x*, it is impossible to determine the distinct values of *p1* and *p2* (these are non-identifiable).

Catalyst contains a special extension for carrying out structural identifiability analysis using the [StructuralIdentifiability.jl](https://github.com/SciML/StructuralIdentifiability.jl) package. This enables StructuralIdentifiability's `assess_identifiability`, `assess_local_identifiability`, and `find_identifiable_functions` functions to be called directly on Catalyst `ReactionSystems`. It also implements specialised routines to make these more efficient when applied to reaction network models. How to use these functions are described in the following tutorial, with [StructuralIdentifiability providing a more extensive documentation](https://docs.sciml.ai/StructuralIdentifiability/stable/). If you use this in your research, please [cite the StructuralIdentifiability.jl](@ref structural_identifiability_citation) and [Catalyst.jl](@ref catalyst_citation) publications.

Structural identifiability can be divided into *local* and *global* identifiability. If a model quantity (which can either be a parameter or an initial condition) is locally identifiable, it means that its true value can be determined down to a finite-number of possible options. This also means that there is some limited region around its true value, where the true value is the only possible value. Globally identifiable quantities' values, on the other hand, can be uniquely determined. Again, while parameter (and initial condition) identifiability can be confirmed structurally for a model, it does not necessarily mean that they are practically* identifiable for some given data.


## Global identifiability analysis

### Basic example

Local identifiability can be assessed using the `assess_identifiability` function. For each model quantity (parameters and initial conditions), it will asses whether they are:
- globally identifiable.
- locally identifiable.
- Unidentifiable.
  
To it, we provide our `ReactionSystem` model and a list of quantities that we are able to measure. Here, we consider a ... model. Let us say that we are able to measure teh values of ... and ..., we provide these at the `measured_quantities` argument. We can now assess identifiability in the following way:
```example si1
```
Here, ... are determined to be globally identifiable (and could theoretically be determined from data) and ... are locally identifiable (and for each, a finite number of candidate values can be determined from the data). Finally, ... are unidentifiable, and cannot be determined from data.

### Indicating known parameters
In the previous case we assumed that all parameters are unknown, however, this is not necessarily true. If there are parameters which value's are known, we can supply these using the `known_ps` argument. Indeed, this might turn other, previously unidentifiable, parameters identifiable. Let us consider this simple example:
```example si2
using Catalyst, StructuralIdentifiability # hide
rn = @reaction_network begin
    (p1+p2, d), 0 <--> X
end
```
Typically, the two production parameters ($p1$ and $p2$) are unidentifiable. However, we we already know the value of $p1$, then $p2$#s value becomes identifiable:
```example si2
assess_identifiability(rn; measured_quantities=[:X], known_ps=[:p1])
```

### Providing non-trivial measured quantities
Sometimes, we are not actually measuring species species, but rather some combinations of species (or possibly parameters). Here, any algebraic expression can be used in `measured_quantities`. If so, used species and parameters have to first be `@unpack`'ed from the model. Say that we have a model where an enzyme ($E$) is converted between an active and inactive form, which in turns activates the production of a product, $P$:
```example si3
using Catalyst, StructuralIdentifiability # hide
enzyme_activation = @reaction_network begin
    (kA,kD), Ei <--> Ea
    (Ea, d), 0 <-->P
end
```
and we can measure the total amount of $E$ ($=$Ei+Ea$), as well as the amount of $P$, we can use the following to assess identifiability:
```example si3
@unpack Ea, Ei = enzyme_activation
assess_identifiability(enzyme_activation; measured_quantities=[Ei+Ea, :P])
```

### Checking identifiability of specific quantities
By default, `asses_identifiability` assesses the identifiability of all parameters and species initial conditions. Sometimes, it is desirable to assess the identifiability of specific quantities. This can be done through the `funcs_to_check` argument. Let us consider our previous example, but say that we can only measure the amount of active enzyme ($Ea$), as well as the product ($P$). If we wish to determine whether the total amount of enzyme ($Ei+Ea$) is identifiable, we could use the following (again using `@unpack` to enable the formation of algebraic expression using the specific quantities):
```example si3
assess_identifiability(enzyme_activation; measured_quantities=[:Ei, :P], funcs_to_check=[Ei+Ea])
```

### Probability of correctness
The identifiabiltiy methods used can, in theory, produce erroneous results. However, it is possible to adjust the lower bound for the probability of correctness using the argument `p` (by default set to `0.99`, that is, at least a $99%$ chance of correctness). We can e.g. increase the bound through:
```example si2
assess_identifiability(rn; measured_quantities=[:X], p=0.999)
nothing # hide
```
giving a minimum bound of $99.9%$  chance of correctness. In practise, the bounds used by StructuralIdentifiability are very conservative, which means that while the minimum guaranteed probability of correctness in the default case is $99%$, in practise it is higher. While increasing teh value of `p` increases the certainty of correctness, it will also increase the time required to assess correctness.

## Local identifiability analysis
Local identifiability can be assessed using the `assess_local_identifiability` function. While this is already determined by the `assess_identifiability` function, local identifiability have the advantage that it is easier to compute. Hence, there might be models where global identifiability analysis fails (or takes prohibitively long time), where instead `assess_local_identifiability` can be used. This functions takes the same inputs as `assess_identifiability` and returns, for each quantity, `true` if iti is locally identifiable and `false` if it is not. Here we assesses local identifiability for the same model as used in the previous example:
```example si1
```
We note that all parameters that `assess_identifiability` determined as either globally or locally identifiable are determined to be locally identifiable, while teh remaining are considered unidentifiable.


## Finding identifiable functions
Finally, StructuralIdentifiability provides the `find_identifiable_functions` function. Rather than determining the identifiability of each parameter and initial condition of the model, it finds a minimal set of identifiable functions, such as any other identifiable expression of the model can be generated by these. Here, let us consider the following model ...
```example si5
```

The `find_identifiable_functions` functions tries to simplify its output functions to create nice expression. The degree to which it does this can be adjusted using the `simplify` keywords. Using the `:weak`, `:standard` (default), and `:strong` arguments, increased simplification can be forced (at the expense of longer runtimes).

## Creating StructuralIdentifiability compatible ODE models from Catalyst `ReactionSystem`s
While the functionality described above covers the vast majority of analysis as user might want to perform, the StructuralIdentifiability package supports several additional features . While these does not have inherent Catalyst support, we do provide the `make_si_ode` function to simplify their use. Similarly to the previous functions, it takes a `ReactionSystem` and lists of measured quantities and known parameter values. The output is a [ode of the standard form supported by StructuralIdentifiability](https://docs.sciml.ai/StructuralIdentifiability/stable/tutorials/creating_ode/#Defining-the-model-using-@ODEmodel-macro). It can be created using the following syntax:
```example si4
using Catalyst, StructuralIdentifiability # hide

```
and then used as input to various StructuralIdentifiability functions. In the following example we use the produced ode to
```example si4

```

---
## [Citation](@id structural_identifiability_citation)
If you use this functionality in your research, please cite the following paper to support the authors of the StructuralIdentifiability package:
```
@article{structidjl,
  author  = {Dong, R. and Goodbrake, C. and Harrington, H. and Pogudin G.},
  title   = {Differential Elimination for Dynamical Models via Projections with Applications to Structural Identifiability},
  journal = {SIAM Journal on Applied Algebra and Geometry},
  url     = {https://doi.org/10.1137/22M1469067},
  year    = {2023}
  volume  = {7},
  number  = {1},
  pages   = {194-235}
}
```

---
## References
[^1]: [Guillaume H.A. Joseph et al., *Introductory overview of identifiability analysis: A guide to evaluating whether you have the right type of data for your modeling purpose*, Environmental Modelling & Software (2019).](https://www.sciencedirect.com/science/article/pii/S1364815218307278)