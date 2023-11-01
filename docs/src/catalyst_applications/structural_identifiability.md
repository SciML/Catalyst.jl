# [Structural Identifiability Analysis](@id structural_identifiability)
During parameter fitting, parameter values are inferred from data. Identifiability is a concept describing to what extent the identification of parameter values for a certain model is actually feasible. Ideally, parameter fitting should always be accompanied with an identifiability analysis of the problem. Identifiability can be divided into *structural* and *practical* identifiability[^1]. Structural identifiability considers only  the system and what quantities we can measure to determine which quantities can be identified. Practical identifiability instead considers the available data, and determines what system quantities can be infeed from it. Generally, in the hypothetical case of an infinite amount of without noise, practical identifiability becomes identical to structural identifiability. Generally, structural identifiability is assessed before parameters are fitted, while practical identifiability is assessed afterwards.

Structural identifiability can be illustrated in the following example network:
${dx \over dt} = p1*p2*x(t)$
where, however much data I have on *x*, it is impossible to determine the values of *p1* and *p2* (these are non-identifiable).

Catalyst contains a special extension for carrying out structural identifiability analysis using the [StructuralIdentifiability.jl](https://github.com/SciML/StructuralIdentifiability.jl) package. This enables StructuralIdentifiability's `assess_identifiability` and `assess_local_identifiability` functions to be called directly on Catalyst `ReactionSystems`. How to use these functions are described in the following tutorial, with [StructuralIdentifiability providing a more extensive documentation](https://docs.sciml.ai/StructuralIdentifiability/stable/). If you use this in your research, please [cite the StructuralIdentifiability.jl](@ref structural_identifiability_citation) and [Catalyst.jl](@ref catalyst_citation) publications.

## Global sensitivity analysis

## Local sensitivity analysis

## Creating StructuralIdentifiability compatible ODE models from Catalyst `ReactionSystem`s


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