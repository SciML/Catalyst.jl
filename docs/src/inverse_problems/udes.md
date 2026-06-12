# [Fitting universal differential equations](@id udes)
```@raw html
<details><summary><strong>Environment setup and package installation</strong></summary>
```
The following code sets up an environment for running the code on this page.
```julia
using Pkg
Pkg.activate(; temp = true) # Creates a temporary environment, which is deleted when the Julia session ends.
Pkg.add("Catalyst")
Pkg.add("DataFrames")
Pkg.add("Distributions")
Pkg.add("Lux")
Pkg.add("ModelingToolkitNeuralNets")
Pkg.add("OrdinaryDiffEq")
Pkg.add("Optimisers")
Pkg.add("PEtab")
Pkg.add("Plots")
```
```@raw html
</details>
```
```@raw html
<details><summary><strong>Quick-start example</strong></summary>
```
The following code provides a brief example of how to perform parameter fitting using the [Optimization.jl](https://github.com/SciML/Optimization.jl) package.
```julia
# Declare the model (a mutual activation loop).
using Catalyst

# Generate synthetic data to which we will apply the fitting procedure.
using Distributions, OrdinaryDiffEq, Plots

# Create a UDE using Lu and ModelingToolkitNeuralNets.
using ModelingToolkitNeuralNets, Lux

# Fit the model using PEtab.
using PEtab, Optimisers

# Plot the fitted simulation and recovered activation functions.
plot()
plot()
```
```@raw html
</details>
```
  \
Default stuff open meues with env and copy example.

## Fitting a rate-based UDE

Short introduction to rate-based UDEs goes here.

Generate synthetic data.

Declare UDE

Fit using PEtab

Plotting the results



---

## [Citations](@id petab_citations)

If you use this functionality in your research, [in addition to Catalyst](@ref doc_index_citation), please cite the following papers to support the authors of the Lux.jl and PEtab.jl packages:

```bibtex
@software{pal2025lux,
  author    = {Pal, Avik},
  title     = {{Lux: Explicit Parameterization of Deep Neural Networks in Julia}},
  publisher = {Zenodo},
  year      = {2025},
  doi       = {10.5281/zenodo.7808903},
  url       = {https://doi.org/10.5281/zenodo.7808903}
}
```
```bibtex
@article{PEtabBioinformatics2025,
  title={PEtab.jl: advancing the efficiency and utility of dynamic modelling},
  author={Persson, Sebastian and Fr{\"o}hlich, Fabian and Grein, Stephan and Loman, Torkel and Ognissanti, Damiano and Hasselgren, Viktor and Hasenauer, Jan and Cvijovic, Marija},
  journal={Bioinformatics},
  volume={41},
  number={9},
  pages={btaf497},
  year={2025},
  publisher={Oxford University Press}
}
```

---

## References

[^1]:
    [Rackauckas, R et al. _Universal differential equations for scientific machine learning_, arXiv (2020).](https://arxiv.org/abs/2001.04385)

[^2]:
    [Loman, T and R, Baker. _Functional and parametric identifiability for universal differential equations applied to chemical reaction networks_, arXiv (2025).](https://arxiv.org/abs/2510.14140)
