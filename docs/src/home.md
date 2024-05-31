# [Catalyst.jl for Reaction Network Modeling](@id doc_home)

Catalyst.jl is a symbolic modeling package for analysis and high-performance
simulation of chemical reaction networks. Catalyst defines symbolic
[`ReactionSystem`](https://docs.sciml.ai/Catalyst/stable/catalyst_functionality/programmatic_CRN_construction/)s,
which can be created programmatically or easily
specified using Catalyst's domain-specific language (DSL). Leveraging
[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), Catalyst enables
large-scale simulations through auto-vectorization and parallelism. Symbolic
`ReactionSystem`s can be used to generate ModelingToolkit-based models, allowing
the easy simulation and parameter estimation of mass action ODE models, chemical
Langevin SDE models, stochastic chemical kinetics jump process models, and more.
Generated models can be used with solvers throughout the broader
[SciML](https://sciml.ai) ecosystem, including higher-level SciML packages (e.g.
for sensitivity analysis, parameter estimation, machine learning applications,
etc).

## [Features](@id doc_home_features)

#### [Features of Catalyst](@id doc_home_features_catalyst)
- [The Catalyst DSL](@ref ref) provides a simple and readable format for manually specifying reaction 
 network models using chemical reaction notation.
- Catalyst `ReactionSystem`s provides a symbolic representation of reaction networks,
 built on [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) and
 [Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/).
- The [Catalyst.jl API](http://docs.sciml.ai/Catalyst/stable/api/catalyst_api) provides functionality 
 for extending networks, building networks programmatically, and for composing 
 multiple networks together.
- Generated models can be simulated using any
 [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/)
 [ODE/SDE/jump solver](@ref ref), and can be used within `EnsembleProblem`s for carrying
 out [parallelized parameter sweeps and statistical sampling](@ref ref). Plot recipes
 are available for [visualization of all solutions](@ref ref).
- Non-integer (e.g. `Float64`) stoichiometric coefficients [are supported](@ref ref) for generating
 ODE models, and symbolic expressions for stoichiometric coefficients [are supported](@ref ref) for
 all system types.
- A [network analysis suite](@ref ref) permits the computation of linkage classes, deficiencies, and
 reversibilities.
- [Conservation laws can be detected and utilized](@ref ref) to reduce system sizes, and to generate
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

#### [Features of Catalyst composing with other packages](@id doc_home_features_composed)
- [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) Can be used to [perform model ODE 
 simulations](@ref ref).
- [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) Can be used to [perform model 
 SDE simulations](@ref ref).
- [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl) Can be used to [model jump 
 simulations](@ref ref).
- Support for [parallelization of all simulations](@ref ref), including parallelization of 
 [ODE simulations on GPUs](@ref ref) using
 [DiffEqGPU.jl](https://github.com/SciML/DiffEqGPU.jl).
- [Latexify](https://korsbo.github.io/Latexify.jl/stable/) can be used to [generate LaTeX 
 expressions](@ref ref) corresponding to generated mathematical models or the
 underlying set of reactions.
- [Graphviz](https://graphviz.org/) can be used to generate and [visualize reaction network graphs](@ref ref) 
 (reusing the Graphviz interface created in [Catlab.jl](https://algebraicjulia.github.io/Catlab.jl/stable/).)
- Model steady states can be computed through homotopy continuation using [HomotopyContinuation.jl](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl)
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
 [global sensitivity analysis](@ref ref) of model behaviors.
 
#### [Features of packages built upon Catalyst](@id doc_home_features_other_packages)
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
- [BondGraphs.jl](https://github.com/jedforrest/BondGraphs.jl), a package for
 constructing and analyzing bond graphs models, which can take Catalyst models as input.
- [PEtab.jl](https://github.com/sebapersson/PEtab.jl), a package that implements the PEtab format for 
 fitting reaction network ODEs to data. Input can be provided either as SBML files or as Catalyst 
 `ReactionSystem`s.

## [How to read this documentation](@id doc_home_documentation)
The Catalyst documentation is separated into sections describing Catalyst's various features. Where appropriate, some sections will also give advice on best practices for various modeling workflows, and provide links with further reading. Each section also contains a set of relevant example workflows. Finally, the [API](@ref api) section contains a list of all functions exported by Catalyst (as well as descriptions of them and their inputs and outputs).

New users are recommended to start with either the [Introduction to Catalyst and Julia for New Julia users](@ref catalyst_for_new_julia_users) or [Introduction to Catalyst](@ref introduction_to_catalyst) sections (depending on whether they are familiar with Julia programming or not). This should be enough to carry out many basic Catalyst workflows. Next, [The Catalyst DSL](@ref ref) section gives a more throughout introduction to model creation, while the [Introduction to model simulation](@ref ref) section more through describes how simulations are carried out. Once you have gotten started using Catalyst, you can read whichever sections are relevant to your work (they should all be clearly labelled).

This documentation contains code which is dynamically run whenever it is built. If you copy the code and run it in your Julia environment it should work. The exact Julia environment that is used in this documentation can be found [here](@ref doc_home_reproducibility).

For most code blocks in this documentation, the output of the last line of code is printed at the of the block, e.g.
```@example home1
1 + 2
```
and
```@example home1
@reaction_network begin
    (p,d), 0 <--> X
end
```
However, in some situations (e.g. when output is extensive, or irrelevant to what is currently being described) we have disabled this, e.g. like here:
```@example home1
1 + 2
nothing # hide
```
and here:
```@example home1
@reaction_network begin
    (p,d), 0 <--> X
end
nothing # hide
```

## [Installation](@id doc_home_installation)
Catalyst is an officially registered Julia package, which can be installed through the Julia package manager:
```julia
using Pkg
Pkg.add("Catalyst")
```

Many Catalyst features require the installation of additional packages. E.g. for ODE-solving and simulation plotting
```julia
Pkg.add("OrdinaryDiffEq")
Pkg.add("Plots")
```
is also needed.

A more throughout guide for setting up Catalyst and installing Julia packages can be found [here](@ref catalyst_for_new_julia_users_packages).

## [Illustrative example](@id doc_home_example)

#### [Deterministic ODE simulation of Michaelis-Menten enzyme kinetics](@id doc_home_example_ode)
Here we show a simple example where a model is created using the Catalyst DSL, and then simulated as 
an ordinary differential equation.

```@example home_simple_example
# Fetch required packages.
using Catalyst, OrdinaryDiffEq, Plots

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

#### [Stochastic jump simulations](@id doc_home_example_jump)
The same model can be used as input to other types of simulations. E.g. here we instead perform a 
jump simulation
```@example home_simple_example
# Create and simulate a jump process (here using Gillespie's direct algorithm).
using JumpProcesses
dprob = DiscreteProblem(model, u0, tspan, ps)
jprob = JumpProblem(model, dprob, Direct())
jump_sol = solve(jprob, SSAStepper())
jump_sol = solve(jprob, SSAStepper(); seed = 1234) # hide
plot(jump_sol; lw = 2)
```

## [Elaborate example](@id doc_home_elaborate_example)
In the above example, we used basic Catalyst-based workflows to simulate a simple model. Here we 
instead show how various Catalyst features can compose to create a much more advanced model. Our 
model describes how the volume of a cell ($V$) is affected by a growth factor ($G$). The growth 
factor only promotes growth while in its phosphorylated form ($Gᴾ$). The phosphorylation of $G$ 
($G \to Gᴾ$) is promoted by sunlight (modeled as the cyclic sinusoid $kₐ*(sin(t)+1)$), which
phosphorylates the growth factor (producing $Gᴾ$). When the cell reaches a critical volume ($V$)
it undergoes through cell division. First, we declare our model:
```@example home_elaborate_example
using Catalyst
cell_model = @reaction_network begin
    @parameters Vₘₐₓ g Ω
    @default_noise_scaling Ω
    @equations begin
        D(V) ~ g*Gᴾ
    end
    @continuous_events begin
        [V ~ Vₘₐₓ] => [V ~ V/2]
    end
    kₚ*(sin(t)+1)/V, G --> Gᴾ
    kᵢ/V, Gᴾ --> G
end
```
Next, we can use [Latexify.jl](https://korsbo.github.io/Latexify.jl/stable/) to show the ordinary differential equations associated with this model:
```@example home_elaborate_example
using Latexify
latexify(cell_model; form = :ode)
```
In this case, we would instead like to perform stochastic simulations, so we transform our model to an SDE:
```@example home_elaborate_example
u0 = [:V => 0.5, :G => 1.0, :Gᴾ => 0.0]
tspan = (0.0, 20.0)
ps = [:Vₘₐₓ => 1.0, :g => 0.2, :kₚ => 5.0, :kᵢ => 2.0, :Ω => 0.1]
sprob = SDEProblem(cell_model, u0, tspan, ps)
```
Finally, we simulate it and plot the result.
```@example home_elaborate_example
using StochasticDiffEq
sol = solve(sprob, STrapezoid())
sol = solve(sprob, STrapezoid(); seed = 1234) # hide
plot(sol; xguide = "Time (au)", lw = 2)
```

## [Getting Help](@id doc_home_help)
Catalyst developers are active on the [Julia Discourse](https://discourse.julialang.org/), 
the [Julia Slack](https://julialang.slack.com) channels \#sciml-bridged and \#sciml-sysbio, and the 
[Julia Zulip sciml-bridged channel](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged). 
For bugs or feature requests, [open an issue](https://github.com/SciML/Catalyst.jl/issues).

## [Supporting and Citing Catalyst.jl](@id doc_home_citation)
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

## [Reproducibility](@id doc_home_reproducibility)
```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```
```@example
using Pkg # hide
Pkg.status() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>and using this machine and Julia version.</summary>
```
```@example
using InteractiveUtils # hide
versioninfo() # hide
```
```@raw html
</details>
```