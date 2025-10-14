# [Catalyst.jl for Reaction Network Modeling](@id doc_index)

Catalyst.jl is a symbolic modeling package for analysis and high-performance
simulation of chemical reaction networks. Catalyst defines symbolic
[`ReactionSystem`](@ref)s, which can be created programmatically or easily
specified using Catalyst's domain-specific language (DSL). Leveraging
[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), Catalyst enables
large-scale simulations through auto-vectorization and parallelism. Symbolic
`ReactionSystem`s can be used to generate ModelingToolkit-based models, allowing
the easy simulation and parameter estimation of mass action ODE models, Chemical
Langevin SDE models, stochastic chemical kinetics jump process models, and more.
Generated models can be used with solvers throughout the broader Julia and
[SciML](https://sciml.ai) ecosystems, including higher-level SciML packages (e.g.
for sensitivity analysis, parameter estimation, machine learning applications,
etc).

## [Features](@id doc_index_features)

#### [Features of Catalyst](@id doc_index_features_catalyst)

- [The Catalyst DSL](@ref dsl_description) provides a simple and readable format for manually specifying reaction network models using chemical reaction notation.
- Catalyst `ReactionSystem`s provides a symbolic representation of reaction networks, built on [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) and [Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/).
- The [Catalyst.jl API](@ref api) provides functionality for building networks programmatically and for composing multiple networks together.
- Leveraging ModelingToolkit, generated models can be converted to symbolic reaction rate equation ODE models, symbolic Chemical Langevin Equation models, and symbolic stochastic chemical kinetics (jump process) models. These can be simulated using any [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) [ODE/SDE/jump solver](@ref simulation_intro), and can be used within `EnsembleProblem`s for carrying out [parallelized parameter sweeps and statistical sampling](@ref ensemble_simulations). Plot recipes are available for [visualization of all solutions](@ref simulation_plotting).
- Non-integer (e.g. `Float64`) stoichiometric coefficients [are supported](@ref dsl_description_stoichiometries_decimal) for generating ODE models, and symbolic expressions for stoichiometric coefficients [are supported](@ref parametric_stoichiometry) for all system types.
- A [network analysis suite](@ref network_analysis_structural_aspects) permits the computation of linkage classes, deficiencies, reversibility, and other network properties.
- [Conservation laws can be detected and utilized](@ref conservation_laws) to reduce system sizes, and to generate non-singular Jacobians (e.g. during conversion to ODEs, SDEs, and steady state equations).
- Catalyst reaction network models can be [coupled with differential and algebraic equations](@ref constraint_equations_coupling_constraints) (which are then incorporated during conversion to ODEs, SDEs, and steady state equations).
- Models can be [coupled with events](@ref constraint_equations_events) that affect the system and its state during simulations.
- By leveraging ModelingToolkit, users have a variety of options for generating optimized system representations to use in solvers. These include construction of [dense or sparse Jacobians](@ref ode_simulation_performance_sparse_jacobian), [multithreading or parallelization of generated derivative functions](@ref ode_simulation_performance_parallelisation), [automatic classification of reactions into optimized jump types for Gillespie type simulations](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#jump_types), [automatic construction of dependency graphs for jump systems](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-Requiring-Dependency-Graphs), and more.
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) symbolic expressions and Julia `Expr`s can be obtained for all rate laws and functions determining the deterministic and stochastic terms within resulting ODE, SDE, or jump models.
- [Steady states](@ref homotopy_continuation) (and their [stabilities](@ref steady_state_stability)) can be computed for model ODE representations.

#### [Features of Catalyst composing with other packages](@id doc_index_features_composed)

- [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) Can be used to numerically solve generated reaction rate equation ODE models.
- [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) can be used to numerically solve generated Chemical Langevin Equation SDE models.
- [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl) can be used to numerically sample generated Stochastic Chemical Kinetics Jump Process models.
- Support for [parallelization of all simulations](@ref ode_simulation_performance_parallelisation), including parallelization of [ODE](@ref ode_simulation_performance_parallelisation_GPU) and [SDE](@ref sde_simulation_performance_parallelisation_GPU) simulations on GPUs using [DiffEqGPU.jl](https://github.com/SciML/DiffEqGPU.jl).
- [Latexify](https://korsbo.github.io/Latexify.jl/stable/) can be used to [generate LaTeX expressions](@ref visualisation_latex) corresponding to generated mathematical models or the underlying set of reactions.
- [GraphMakie](https://docs.makie.org/stable/) recipes are provided that can be used to generate and [visualize reaction network graphs](@ref visualisation_graphs)
- Model steady states can be [computed through homotopy continuation](@ref homotopy_continuation) using [HomotopyContinuation.jl](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl) (which can find *all* steady states of systems with multiple ones), by [forward ODE simulations](@ref steady_state_solving_simulation) using [SteadyStateDiffEq.jl](https://github.com/SciML/SteadyStateDiffEq.jl), or by [numerically solving steady-state nonlinear equations](@ref steady_state_solving_nonlinear) using [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl).
- [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) can be used to compute bifurcation diagrams of model steady states (including finding periodic orbits).
- [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl) can be used to compute model [basins of attraction](@ref dynamical_systems_basins_of_attraction), [Lyapunov spectrums](@ref dynamical_systems_lyapunov_exponents), and other dynamical system properties.
- [Optimization.jl](https://github.com/SciML/Optimization.jl) and [PEtab.jl](https://github.com/sebapersson/PEtab.jl) can all be used to [fit model parameters to data](https://sebapersson.github.io/PEtab.jl/stable/Define_in_julia/).
- [GlobalSensitivity.jl](https://github.com/SciML/GlobalSensitivity.jl) can be used to perform [global sensitivity analysis](@ref global_sensitivity_analysis) of model behaviors.
- [SciMLSensitivity.jl](https://github.com/SciML/SciMLSensitivity.jl) can be used to compute local sensitivities of functions containing forward model simulations.
- [StructuralIdentifiability.jl](https://github.com/SciML/StructuralIdentifiability.jl) can be used to perform structural identifiability analysis.

#### [Features of packages built upon Catalyst](@id doc_index_features_other_packages)

- Catalyst [`ReactionSystem`](@ref)s can be [imported from SBML files](@ref model_file_import_export_sbml) via [SBMLImporter.jl](https://github.com/sebapersson/SBMLImporter.jl) and [SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl), and [from BioNetGen .net files](@ref model_file_import_export_sbml_rni_net) and various stoichiometric matrix network representations using [ReactionNetworkImporters.jl](https://github.com/SciML/ReactionNetworkImporters.jl).
- [MomentClosure.jl](https://github.com/augustinas1/MomentClosure.jl) allows generation of symbolic ModelingToolkit `ODESystem`s that represent moment closure approximations to moments of the Chemical Master Equation, from reaction networks defined in Catalyst.
- [FiniteStateProjection.jl](https://github.com/kaandocal/FiniteStateProjection.jl) allows the construction and numerical solution of Chemical Master Equation models from reaction networks defined in Catalyst.
- [DelaySSAToolkit.jl](https://github.com/palmtree2013/DelaySSAToolkit.jl) can augment Catalyst reaction network models with delays, and can simulate the resulting stochastic chemical kinetics with delays models.

## [How to read this documentation](@id doc_index_documentation)

The Catalyst documentation is separated into sections describing Catalyst's various features. Where appropriate, some sections will also give advice on best practices for various modeling workflows, and provide links with further reading. Each section also contains a set of relevant example workflows. Finally, the [API](@ref api) section contains a list of all functions exported by Catalyst (as well as descriptions of them and their inputs and outputs).

New users are recommended to start with either the [Introduction to Catalyst and Julia for New Julia users](@ref catalyst_for_new_julia_users) or [Introduction to Catalyst](@ref introduction_to_catalyst) sections (depending on whether they are familiar with Julia programming or not). This should be enough to carry out many basic Catalyst workflows.

This documentation contains code which is dynamically run whenever it is built. If you copy the code and run it in your Julia environment it should work. The exact Julia environment that is used in this documentation can be found [here](@ref doc_index_reproducibility).

For most code blocks in this documentation, the output of the last line of code is printed at the of the block, e.g.

```@example home_display
1 + 2
```

and

```@example home_display
using Catalyst # hide
@reaction_network begin
    (p,d), 0 <--> X
end
```

However, in some situations (e.g. when output is extensive, or irrelevant to what is currently being described) we have disabled this, e.g. like here:

```@example home_display
1 + 2
nothing # hide
```

and here:

```@example home_display
@reaction_network begin
    (p,d), 0 <--> X
end
nothing # hide
```

## [Installation](@id doc_index_installation)

Catalyst is an officially registered Julia package, which can be installed through the Julia package manager:

```julia
using Pkg
Pkg.add("Catalyst")
```

Many Catalyst features require the installation of additional packages. E.g. for ODE-solving and simulation plotting

```julia
Pkg.add("OrdinaryDiffEqDefault")
Pkg.add("Plots")
```

is also needed.

It is **strongly** recommended to install and use Catalyst in its own environment with the
minimal set of needed packages. For example, to install Catalyst and Plots in a
new environment named `catalyst_project` (saved in the current directory) one
can say

```julia
Pkg.activate("catalyst_project")
Pkg.add("Catalyst")
Pkg.add("Plots")
```

If one would rather just create a temporary environment that is not saved when
exiting Julia you can say

```julia
Pkg.activate(; temp = true)
Pkg.add("Catalyst")
Pkg.add("Plots")
```

After installation, we suggest running

```julia
Pkg.status("Catalyst")
```

to confirm that the latest version of Catalyst was installed (and not an older version).
If you have installed into a new environment this should always be the case. However, if you
installed into an existing environment, such as the default Julia global environment, the presence
of incompatible versions of other pre-installed packages could lead to older versions of Catalyst
being installed. In this case we again recommend creating a new environment and installing Catalyst
there to obtain the latest version.

A more thorough guide for setting up Catalyst and installing Julia packages can be found [here](@ref catalyst_for_new_julia_users_packages).

## [Illustrative example](@id doc_index_example)

#### [Deterministic ODE simulation of Michaelis-Menten enzyme kinetics](@id doc_index_example_ode)

Here we show a simple example where a model is created using the Catalyst DSL, and then simulated as
an ordinary differential equation.

```@example home_simple_example
# Fetch required packages.
using Catalyst, OrdinaryDiffEqDefault, Plots

# Create model.
model = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end

# Create an ODE that can be simulated.
u0 = [:S => 50.0, :E => 10.0, :SE => 0.0, :P => 0.0]
tspan = (0., 200.)
ps = [:kB => 0.01, :kD => 0.1, :kP => 0.1]
ode = ODEProblem(model, u0, tspan, ps)

# Simulate ODE and plot results.
sol = solve(ode)
plot(sol; lw = 5)
```

#### [Stochastic jump simulations](@id doc_index_example_jump)

The same model can be used as input to other types of simulations. E.g. here we instead generate and simulate a stochastic chemical kinetics jump process model.

```@example home_simple_example
# Create and simulate a jump process (here using Gillespie's direct algorithm).
# The initial conditions are now integers as we track exact populations for each species.
using JumpProcesses
u0_integers = [:S => 50, :E => 10, :SE => 0, :P => 0]
jinput = JumpInputs(model, u0_integers, tspan, ps)
jprob = JumpProblem(jinput)
jump_sol = solve(jprob)
jump_sol = solve(jprob; seed = 1234) # hide
plot(jump_sol; lw = 2)
```

## [More elaborate example](@id doc_index_elaborate_example)

In the above example, we used basic Catalyst workflows to simulate a simple
model. Here we instead show how various Catalyst features can compose to create
a much more advanced model. Our model describes how the volume of a cell ($V$)
is affected by a growth factor ($G$). The growth factor only promotes growth
while in its phosphorylated form ($G^P$). The phosphorylation of $G$ ($G \to G^P$)
is promoted by sunlight, which is modeled as the cyclic sinusoid $k_a (\sin(t) + 1)$.
When the cell reaches a critical volume ($V_m$) it undergoes cell division. First, we
declare our model:

```@example home_elaborate_example
using Catalyst
cell_model = @reaction_network begin
    @parameters Vₘ g
    @equations begin
        D(V) ~ g*Gᴾ
    end
    @continuous_events begin
        [V ~ Vₘ] => [V ~ V/2]
    end
    kₚ*(sin(t)+1)/V, G --> Gᴾ
    kᵢ/V, Gᴾ --> G
end
```

We now study the system as a Chemical Langevin Dynamics SDE model, which can be generated as follows

```@example home_elaborate_example
u0 = [:V => 25.0, :G => 50.0, :Gᴾ => 0.0]
tspan = (0.0, 20.0)
ps = [:Vₘ => 50.0, :g => 0.3, :kₚ => 100.0, :kᵢ => 60.0]
sprob = SDEProblem(cell_model, u0, tspan, ps)
```

This problem encodes the following stochastic differential equation model:

```math
\begin{align*}
dG(t) &=  - \left( \frac{k_p(\sin(t)+1)}{V(t)} G(t) + \frac{k_i}{V(t)} G^P(t) \right) dt - \sqrt{\frac{k_p (\sin(t)+1)}{V(t)} G(t)} \, dW_1(t) + \sqrt{\frac{k_i}{V(t)} G^P(t)} \, dW_2(t) \\
dG^P(t) &= \left( \frac{k_p(\sin(t)+1)}{V(t)} G(t) - \frac{k_i}{V(t)} G^P(t) \right) dt + \sqrt{\frac{k_p (\sin(t)+1)}{V(t)} G(t)} \, dW_1(t) - \sqrt{\frac{k_i}{V(t)} G^P(t)} \, dW_2(t) \\
dV(t) &= \left(g \, G^P(t)\right) dt
\end{align*}
```

where the $dW_1(t)$ and $dW_2(t)$ terms represent independent Brownian Motions, encoding the noise added by the Chemical Langevin Equation. Finally, we can simulate and plot the results.

```@example home_elaborate_example
using StochasticDiffEq, Plots
sol = solve(sprob, EM(); dt = 0.05)
sol = solve(sprob, EM(); dt = 0.05, seed = 1234) # hide
plot(sol; xguide = "Time (au)", lw = 2)
```

## [Getting Help](@id doc_index_help)

Catalyst developers are active on the [Julia Discourse](https://discourse.julialang.org/) and
the [Julia Slack](https://julialang.slack.com) channels \#sciml-bridged and \#sciml-sysbio.
For bugs or feature requests, [open an issue](https://github.com/SciML/Catalyst.jl/issues).

## [Supporting and Citing Catalyst.jl](@id doc_index_citation)

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

## [Reproducibility](@id doc_index_reproducibility)

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
