# [Structural Identifiability Analysis](@id structural_identifiability)
During parameter fitting, parameter values are inferred from data. Identifiability is a concept describing to what extent the identification of parameter values for a certain model is actually feasible. Ideally, parameter fitting should always be accompanied with an identifiability analysis of the problem. Identifiability can be divided into *structural* and *practical* identifiability[^1]. Structural identifiability considers only  the system and what quantities we can measure to determine which quantities can be identified. Practical identifiability instead considers the available data, and determines what system quantities can be inferred from it. Generally, in the hypothetical case of an infinite amounts of noise-less data, practical identifiability becomes identical to structural identifiability. Generally, structural identifiability is assessed before parameters are fitted, while practical identifiability is assessed afterwards.

Structural identifiability (which is what this tutorial considers) can be illustrated by the following differential equation:
${dx \over dt} = p1*p2*x(t)$
where, however much data is collected on *x*, it is impossible to determine the distinct values of *p1* and *p2*. Hence, these parameters are these are non-identifiable.

Catalyst contains a special extension for carrying out structural identifiability analysis using the [StructuralIdentifiability.jl](https://github.com/SciML/StructuralIdentifiability.jl) package. This enables StructuralIdentifiability's `assess_identifiability`, `assess_local_identifiability`, and `find_identifiable_functions` functions to be called directly on Catalyst `ReactionSystem`s. It also implements specialised routines to make these more efficient when applied to reaction network models. How to use these functions are described in the following tutorial, with [StructuralIdentifiability providing a more extensive documentation](https://docs.sciml.ai/StructuralIdentifiability/stable/). If you use this in your research, please [cite the StructuralIdentifiability.jl](@ref structural_identifiability_citation) and [Catalyst.jl](@ref catalyst_citation) publications.

Structural identifiability can be divided into *local* and *global* identifiability. If a model quantity (which can either be a parameter or an initial condition) is locally identifiable, it means that its true value can be determined down to a finite-number of possible options. This also means that there is some limited region around its true value, where the true value is the only possible value. Globally identifiable quantities' values, on the other hand, can be uniquely determined. Again, while parameter (and initial condition) identifiability can be confirmed structurally for a model, it does not necessarily mean that they are practically identifiable for some given data.


## Global identifiability analysis

### Basic example

Local identifiability can be assessed using the `assess_identifiability` function. For each model quantity (parameters and initial conditions), it will asses whether they are:
- globally identifiable.
- locally identifiable.
- Unidentifiable.
  
To it, we provide our `ReactionSystem` model and a list of quantities that we are able to measure. Here, we consider a Goodwind oscillator (a simple 3-component model, where the three species $M$, $E$, and $P$ are produced and degraded, which may exhibit oscillations)[^2]. Let us say that we are able to measure the concentration of $M$, we then provide designate this using the `measured_quantities` argument. We can now assess identifiability in the following way:
```example si1
using StructuralIdentifiability, Catalyst
goodwind_oscillator = @reaction_network begin
    (pₘ/(1+P), dₘ), 0 <--> M
    (pₑ*M,dₑ), 0 <--> E
    (pₚ*E,dₚ), 0 <--> P
end
assess_identifiability(goodwind_oscillator; measured_quantities=[:M])
```
From the output, we find that `E`, `pₑ`, and `pₚ` (the initial concentration of $E$, and the production rates of $E$ and $P$, respectively) are non-identifiable. Next, `dₑ` and `dₚ` (the degradation rates of $E$ and $P$, respectively) are locally identifiable. Finally, `P`, `M`, `pₘ`, and `dₘ` (the initial concentrations of `P` and `M`, the production and degradation rate of `M`, respectively) are all globally identifiable. Next, we can assess identifiability in the case where we can measure all three species concentrations:
```example si1
assess_identifiability(goodwind_oscillator; measured_quantities=[:M, :P, :E])
```
in which case all initial conditions and parameters become identifiable.

### Indicating known parameters
In the previous case we assumed that all parameters are unknown, however, this is not necessarily true. If there are parameters which value's are known, we can supply these using the `known_p` argument. Indeed, this might turn other, previously unidentifiable, parameters identifiable. Let us consider the previous example, where we measure the concentration of $M$ only, but also happen to know the production rate of $E$ ($pₑ$):
```example si1
assess_identifiability(gwo; measured_quantities=[:M], known_p=[:pₑ])
```
Not only does this turn the previously non-identifiable `pₑ` (globally) identifiable (which is obvious, given that its value is now known), but this additional information improve identifiability for several other network components.

To, in a similar manner, indicate that certain initial conditions are known is a work in progress. Hopefully it should be an available feature in the near future.

### Providing non-trivial measured quantities
Sometimes, we are not actually measuring species species, but rather some combinations of species (or possibly parameters). Here, any algebraic expression can be used in `measured_quantities`. If so, used species and parameters have to first be `@unpack`'ed from the model. Say that we have a model where an enzyme ($E$) is converted between an active and inactive form, which in turns activates the production of a product, $P$:
```example si2
using Catalyst, StructuralIdentifiability # hide
enzyme_activation = @reaction_network begin
    (kA,kD), Eᵢ <--> Eₐ
    (Eₐ, d), 0 <-->P
end
```
and we can measure the total amount of $E$ ($=$Eᵢ+Eₐ$), as well as the amount of $P$, we can use the following to assess identifiability:
```example si2
@unpack Eᵢ, Eₐ = enzyme_activation
assess_identifiability(enzyme_activation; measured_quantities=[Eᵢ+Eₐ, :P])
nothing # hide
```

### Probability of correctness
The identifiability methods used can, in theory, produce erroneous results. However, it is possible to adjust the lower bound for the probability of correctness using the argument `p` (by default set to `0.99`, that is, at least a $99%$ chance of correctness). We can e.g. increase the bound through:
```example si2
assess_identifiability(goodwind_oscillator; measured_quantities=[:M], p=0.999)
nothing # hide
```
giving a minimum bound of $99.9%$  chance of correctness. In practise, the bounds used by StructuralIdentifiability are very conservative, which means that while the minimum guaranteed probability of correctness in the default case is $99%$, in practise it is higher. While increasing the value of `p` increases the certainty of correctness, it will also increase the time required to assess identifiability.

## Local identifiability analysis
Local identifiability can be assessed through the `assess_local_identifiability` function. While this is already determined by `assess_identifiability`, assessing local identifiability only have the advantage that it is easier to compute. Hence, there might be models where global identifiability analysis fails (or takes prohibitively long time), where instead `assess_local_identifiability` can be used. This functions takes the same inputs as `assess_identifiability` and returns, for each quantity, `true` if it is locally identifiable (and `false` if it is not). Here, for the Goodwind oscillator, we assesses it for local identifiability only:
```example si1
assess_local_identifiability(goodwind_oscillator; measured_quantities=[:M])
```
We note that the results are consistent with those produced by `assess_identifiability` (with globally or locally identifiable quantities here all being assessed as at least locally identifiable).

## Finding identifiable functions
Finally, StructuralIdentifiability provides the `find_identifiable_functions` function. Rather than determining the identifiability of each parameter and initial condition of the model, it finds a minimal set of identifiable functions, such as any other identifiable expression of the model can be generated by these. Let us again consider the Goodwind oscillator, using the `find_identifiable_functions` function we find that identifiability can be reduced two five globally identifiable expressions:
```example si1
find_identifiable_functions(goodwind_oscillator; measured_quantities=[:M])
```
Again, these results are consistent with those of produced by `assess_identifiability`. There, `pₑ` and `pₚ` where found to be globally identifiable. Here, they correspond directly to identifiable expressions. The remaining four parameters (`pₘ`, `dₘ`, `dₑ`, and `dₚ`) occurs as part of more complicated composite expressions.

The `find_identifiable_functions` functions tries to simplify its output functions to create nice expression. The degree to which it does this can be adjusted using the `simplify` keywords. Using the `:weak`, `:standard` (default), and `:strong` arguments, increased simplification can be forced (at the expense of longer runtimes).

## Creating StructuralIdentifiability compatible ODE models from Catalyst `ReactionSystem`s
While the functionality described above covers the vast majority of analysis as user might want to perform, the StructuralIdentifiability package supports several additional features. While these does not have inherent Catalyst support, we do provide the `make_si_ode` function to simplify their use. Similarly to the previous functions, it takes a `ReactionSystem` and lists of measured quantities and known parameter values. The output is a [ode of the standard form supported by StructuralIdentifiability](https://docs.sciml.ai/StructuralIdentifiability/stable/tutorials/creating_ode/#Defining-the-model-using-@ODEmodel-macro). It can be created using the following syntax:
```example si1
si_ode = make_si_ode(goodwind_oscillator; measured_quantities=[:M])
nothing # hide
```
and then used as input to various StructuralIdentifiability functions. In the following example we uses StructuralIdentifiability's `print_for_DAISY` function, printing the model as an expression that can be used by the [DAISY](https://daisy.dei.unipd.it/) software for identifiability analysis[^3].
```example si1
print_for_DAISY(si_ode)
nothing # hide
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
[^2]: [Goodwin B.C., *Oscillatory Behavior in Enzymatic Control Processes*, Advances in Enzyme Regulation (1965).](https://www.sciencedirect.com/science/article/pii/0065257165900671?via%3Dihub)
[^3]: [Bellu G., et al., *DAISY: A new software tool to test global identifiability of biological and physiological systems*, Computer Methods and Programs in Biomedicine (2007).](https://www.sciencedirect.com/science/article/abs/pii/S0169260707001605)