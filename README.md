# Catalyst.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://docs.sciml.ai/Catalyst/stable/)
[![API Stable](https://img.shields.io/badge/API-stable-blue.svg)](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/)
<!---  [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://docs.sciml.ai/Catalyst/dev/)
[![API Dev](https://img.shields.io/badge/API-dev-blue.svg)](https://docs.sciml.ai/Catalyst/dev/api/catalyst_api/) -->

[![Build Status](https://github.com/SciML/Catalyst.jl/workflows/CI/badge.svg)](https://github.com/SciML/Catalyst.jl/actions?query=workflow%3ACI)
[![codecov.io](https://codecov.io/gh/SciML/Catalyst.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/Catalyst.jl)
[![Coverage Status](https://coveralls.io/repos/github/SciML/Catalyst.jl/badge.svg?branch=master)](https://coveralls.io/github/SciML/Catalyst.jl?branch=master)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

Catalyst.jl is a symbolic modeling package for analysis and high-performance
simulation of chemical reaction networks. Catalyst defines symbolic
[`ReactionSystem`](https://docs.sciml.ai/Catalyst/stable/catalyst_functionality/programmatic_CRN_construction/)s,
which can be created programmatically or easily
specified using Catalyst's domain-specific language (DSL). Leveraging
[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), Catalyst enables
large-scale simulations through auto-vectorization and parallelism. Symbolic
`ReactionSystem`s can be used to generate ModelingToolkit-based models, allowing
the easy simulation and parameter estimation of mass action ODE models, Chemical
Langevin SDE models, stochastic chemical kinetics jump process models, and more.
Generated models can be used with solvers throughout the broader
[SciML](https://sciml.ai) ecosystem, including higher-level SciML packages (e.g.
for sensitivity analysis, parameter estimation, machine learning applications,
etc).

## Breaking changes and new features

**NOTE:** Version 14 is a breaking release, prompted by the release of ModelingToolkit.jl version 9. This caused several breaking changes in how Catalyst models are represented and interfaced with.

Breaking changes and new functionality are summarized in the [HISTORY.md](HISTORY.md) file. This also includes a special migration guide for version 14.

## Tutorials and documentation

The latest tutorials and information on using the package are available in the [stable
documentation](https://docs.sciml.ai/Catalyst/stable/). The [in-development
documentation](https://docs.sciml.ai/Catalyst/dev/) describes unreleased features in
the current master branch.

Several YouTube video tutorials and overviews are also available (however, these use older versions of Catalyst, and some notation may be out-of-date):
- From JuliaCon 2023: A short 15-minute overview of Catalyst (version 13) is
available in the talk [Catalyst.jl, Modeling Chemical Reaction Networks](https://www.youtube.com/watch?v=yreW94n98eM&ab_channel=TheJuliaProgrammingLanguage).
- From JuliaCon 2022: A 3-hour tutorial workshop overviewing how to use
 Catalyst (version 12.1) and its more advanced features. [Workshop
 video](https://youtu.be/tVfxT09AtWQ), [Workshop Pluto.jl
 Notebooks](https://github.com/SciML/JuliaCon2022_Catalyst_Workshop). 
- From SIAM CSE 2021: A short 15-minute overview of Catalyst (version 6) is
available in the talk [Modeling Biochemical Systems with
Catalyst.jl](https://www.youtube.com/watch?v=5p1PJE5A5Jw).
- From JuliaCon 2018: A short 13-minute overview of Catalyst (when it was known
 as DiffEqBiological) is available in the talk [Efficient
 Modelling of Biochemical Reaction
 Networks](https://www.youtube.com/watch?v=s1e72k5XD6s)

Finally, an overview of the package and its features (as of version 13) can also be found in its corresponding research paper, [Catalyst: Fast and flexible modeling of reaction networks](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011530).

## Features

#### Features of Catalyst

- [The Catalyst DSL](@ref ref) provides a simple and readable format for manually specifying reaction 
 network models using chemical reaction notation.
- Catalyst `ReactionSystem`s provides a symbolic representation of reaction networks,
 built on [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) and
 [Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/).
- The [Catalyst.jl API](http://docs.sciml.ai/Catalyst/stable/api/catalyst_api) provides functionality 
 for extending networks, building networks programmatically, and for composing 
 multiple networks together.
- Generated systems can be simulated using any
 [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/)
 [ODE/SDE/jump solver](@ref ref), and can be used within `EnsembleProblem`s for carrying
 out [parallelized parameter sweeps and statistical sampling](@ref ref). Plot recipes
 are available for [visualizing of all solutions](@ref ref).
- Non-integer (e.g. `Float64`) stoichiometric coefficients [are supported](@ref ref) for generating
 ODE models, and symbolic expressions for stoichiometric coefficients [are supported](@ref ref) for
 all system types.
- A [network analysis suite](@ref ref) permits the computation of linkage classes, deficiencies, and
 reversibilities.
- [Conservation laws can be detected and utilised](@ref ref) to reduce system sizes, and to generate
 non-singular Jacobians (e.g. during conversion to ODEs, SDEs, and steady state equations).
- Catalyst reaction network models can be [coupled with differential and algebraic equations](@ref ref)
 (which are then incorporated during conversion to ODEs, SDEs, and steady state equations).
- Models can be [coupled with events](@ref ref) that affect the system and its state during simulations.
- By leveraging ModelingToolkit, users have a variety of options for generating
 optimized system representations to use in solvers. These include construction
 of [dense or sparse Jacobians](@ref ref), [multithreading or parallelization of generated
 derivative functions](@ref ref), [automatic classification of reactions into optimized
 jump types for Gillespie type simulations](@ref ref), [automatic construction of
 dependency graphs for jump systems](@ref ref), and more.
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) symbolic
 expressions and Julia `Expr`s can be obtained for all rate laws and functions determining the 
 deterministic and stochastic terms within resulting ODE, SDE or jump models.
- [Steady states](@ref ref) (and their [stabilities](@ref ref)) can be computed for model ODE representations.


#### Features of Catalyst composing with other packages
- [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) Can be used to [perform model ODE 
 simulations](@ref ref).
- [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) Can be used to [perform model 
 SDE simulations](@ref ref).
- [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl) Can be used to [model jump 
 simulations](@ref ref).
- Support for [parallelisation of all simulations]((@ref ref)), including parallelisation of 
 [ODE simulations on GPUs](@ref ref) using
 [DiffEqGPU.jl](https://github.com/SciML/DiffEqGPU.jl).
- [Latexify](https://korsbo.github.io/Latexify.jl/stable/) can be used to [generate LaTeX 
 expressions](@ref ref) corresponding to generated mathematical models or the
 underlying set of reactions.
- [Graphviz](https://graphviz.org/) can be used to generate and [visualize reaction network graphs](@ref ref) 
 (reusing the Graphviz interface created in [Catlab.jl](https://algebraicjulia.github.io/Catlab.jl/stable/).)
- Models steady states can be computed through homotopy continuation using [HomotopyContinuation.jl](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl)
 (which can find *all* steady states of systems with multiple ones), by forward ODE simulations using
 [SteadyStateDiffEq.jl)](https://github.com/SciML/SteadyStateDiffEq.jl), or by nonlinear systems 
 solving using [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl).
- [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) can be used to [compute 
 bifurcation diagrams](@ref ref) of models' steady states (including finding periodic orbits).
- [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl) can be used to compute 
 model [basins of attraction](@ref ref) and [Lyapunov spectrums](@ref ref).
- [StructuralIdentifiability.jl](https://github.com/SciML/StructuralIdentifiability.jl) can be used 
 to [perform structural identifiability analysis](@ref ref).
- [Optimization.jl](https://github.com/SciML/Optimization.jl), [DiffEqParamEstim.jl](https://github.com/SciML/DiffEqParamEstim.jl), 
 and [PEtab.jl](https://github.com/sebapersson/PEtab.jl) can all be used to [fit model parameters to data](@ref ref).
- [GlobalSensitivity.jl](https://github.com/SciML/GlobalSensitivity.jl) can be used to perform 
 [global sensitivity analysis](@ref ref) of model behaviours.
 
#### Features of packages built upon Catalyst
- Catalyst [`ReactionSystem`](@ref)s can be [imported from SBML files](@ref ref) via
 [SBMLImporter.jl](https://github.com/SciML/SBMLImporter.jl) and [SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl), 
 and [from BioNetGen .net files](@ref ref) and various stoichiometric matrix network representations
 using [ReactionNetworkImporters.jl](https://github.com/SciML/ReactionNetworkImporters.jl).
- [MomentClosure.jl](https://github.com/augustinas1/MomentClosure.jl) allows generation of symbolic 
 ModelingToolkit `ODESystem`s that represent moment closure approximations to moments of the 
 Chemical Master Equation, from reaction networks defined in Catalyst.
- [FiniteStateProjection.jl](https://github.com/kaandocal/FiniteStateProjection.jl)
 allows the construction and numerical solution of Chemical Master Equation
 models from reaction networks defined in Catalyst.
- [DelaySSAToolkit.jl](https://github.com/palmtree2013/DelaySSAToolkit.jl) can
 augment Catalyst reaction network models with delays, and can simulate the
 resulting stochastic chemical kinetics with delays models.
- [BondGraphs.jl](https://github.com/jedforrest/BondGraphs.jl) a package for
 constructing and analyzing bond graphs models, which can take Catalyst models as input.
- [PEtab.jl](https://github.com/sebapersson/PEtab.jl) a package that implements the PEtab format for 
 fitting reaction network ODEs to data. Input can be provided either as SBML files or as Catalyst 
 `ReactionSystem`s.


## Illustrative example

#### Deterministic ODE simulation of Michaelis-Menten enzyme kinetics
Here we show a simple example where a model is created using the Catalyst DSL, and then simulated as 
an ordinary differential equation.

```julia
# Fetch required packages.
using Catalyst, DifferentialEquations, Plots

# Create model.
model = @reaction_network begin
 kB, S + E --> SE
 kD, SE --> S + E
 kP, SE --> P + E
end

# Create an ODE that can be simulated.
u0 = [:S => 50, :E => 10, :SE => 0, :P => 0]
tspan = (0., 200.)
ps = (:kB => 0.01, :kD => 0.1, :kP => 0.1)
ode = ODEProblem(model, u0, tspan, ps)

# Simulate ODE and plot results.
sol = solve(ode)
plot(sol; lw = 5)
```

![](https://user-images.githubusercontent.com/1814174/87864114-3bf9dd00-c932-11ea-83a0-58f38aee8bfb.png)

#### Stochastic jump simulations
The same model can be used as input to other types of simulations. E.g. here we instead perform a 
jump simulation
```julia
# Create and simulate a jump process (here using Gillespie's direct algorithm).
dprob = DiscreteProblem(model, u0, tspan, ps)
jprob = JumpProblem(model, dprob, Direct())
jump_sol = solve(jprob, SSAStepper())
plot(jump_sol; lw = 2)
```


## Elaborate example

In the above example, we used basic Catalyst-based workflows to simulate a simple model. Here we instead show how various Catalyst features can compose to create a much more advanced model. Our model describes how the volume of a cell ($V$) is affected by a growth factor ($G$). Typically the growth factor is inactive ($Gi$), but it is activated ($Ga$) by the presence of sunlight (modeled as the cyclic sinusoid $kA*(sin(t)+1)$). When the cell reaches a critical volume ($V$) it goes through cell division. First, we declare our model:
```julia
cell_model = @reaction_network begin
 @parameters V_thres g
 @equations begin
 D(V) ~ g*Ga
 end
 @continuous_events begin
 [V ~ V_thres] => [V ~ V/2]
 end
 (kA*(sin(t)+1), kI), Gi <--> Ga
end
```
Next, we can use [Latexify.jl](https://korsbo.github.io/Latexify.jl/stable/) to show the ordinary differential equations associated with this model:
```julia
using Latexify
latexify(cell_model; form = :ode)
```
In this case we would like to perform stochastic simulations, so we transform our model to an SDE:
```julia
u0 = [:V => 0.5, :Gi => 1.0, :Ga => 0.0]
tspan = (0.0, 10.0)
ps = [:V_thres => 1.0, :g => 0.5, :kA => 5.0, :kI => 2.0]
sprob = SDEProblem(cell_model, u0, tspan, ps)
```
Finally, we simulate it and plot the result.
```julia
sol = solve(oprob, Tsit5())
plot(sol)
```

## Getting help or getting involved
Catalyst developers are active on the [Julia Discourse](https://discourse.julialang.org/), 
the [Julia Slack](https://julialang.slack.com) channels \#sciml-bridged and \#sciml-sysbio, and the 
[Julia Zulip sciml-bridged channel](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged). 
For bugs or feature requests, [open an issue](https://github.com/SciML/Catalyst.jl/issues).

If you are interested in participating in the development of Catalyst, or integrating your package(s) 
with it, developer documentation can be found [here](@ref ref). Are you a student (or similar) who 
wishes to do a Google Summer of Code (or similar) project tied to Catalyst? Information on how to get 
involved, including good first issues to get familiar with working on the package, can be found [here](@ref ref).

## Supporting and citing Catalyst.jl
The software in this ecosystem was developed as part of academic research. If you would like to help 
support it, please star the repository as such metrics may help us secure funding in the future. If 
you use Catalyst as part of your research, teaching, or other activities, we would be grateful if you 
could cite our work:
```
@article{CatalystPLOSCompBio2023,
 doi = {10.1371/journal.pcbi.1011530},
 author = {Loman, Torkel E. AND Ma, Yingbo AND Ilin, Vasily AND Gowda, Shashi AND Korsbo, Niklas AND Yewale, Nikhil AND Rackauckas, Chris AND Isaacson, Samuel A.},
 journal = {PLOS Computational Biology},
 publisher = {Public Library of Science},
 title = {Catalyst: Fast and flexible modeling of reaction networks},
 year = {2023},
 month = {10},
 volume = {19},
 url = {https://doi.org/10.1371/journal.pcbi.1011530},
 pages = {1-19},
 number = {10},
}
```

We also maintain a user survey, asking basic questions about how users utilise the package. The survey 
is available [here](ref), and only takes about 5 minutes to fill out. We are grateful to those who 
fill out the survey, as this helps us further develop the package.