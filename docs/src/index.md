# [Catalyst.jl for Reaction Network Modeling](@id doc_index)

Catalyst.jl is a symbolic modeling package for analysis and high-performance simulation of chemical reaction networks and related dynamical systems. Models can be specified using an intuitive [domain-specific language (DSL)](@ref dsl_description) or constructed [programmatically](@ref programmatic_CRN_construction). Catalyst supports ODE, steady-state ODE, SDE, stochastic chemical kinetics (jump), and [hybrid](@ref simulation_intro) simulations, including models that couple reactions with differential equations, events, and external noise (via Brownian Motions and/or Poisson Processes).

Built on [ModelingToolkitBase.jl](https://github.com/SciML/ModelingToolkit.jl/tree/master/lib/ModelingToolkitBase) and [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), Catalyst leverages symbolic computation for sparsity exploitation, Jacobian construction, dependency graph analysis, and parallelism. Generated models integrate with the broader Julia and [SciML](https://sciml.ai) ecosystems for sensitivity analysis, parameter estimation, bifurcation analysis, and more.

## [Features](@id doc_index_features)

- [**DSL for reaction networks**](@ref dsl_description) — a readable, concise format for specifying models using chemical reaction notation.
- [**Multiple simulation types**](@ref simulation_intro) — generate and simulate ODE, steady-state ODE, SDE, jump, and hybrid models from a single [`ReactionSystem`](@ref).
- [**Coupled models**](@ref constraint_equations_coupling_constraints) — combine reactions with differential equations, [events](@ref events), Brownian noise (`@brownians`), and Poisson jumps (`@poissonians`).
- [**Network analysis**](@ref network_analysis_structural_aspects) — compute linkage classes, deficiencies, reversibility, and other network properties.
- [**Compositional modeling**](@ref compositional_modeling) — build models hierarchically using [`@network_component`](@ref), [`compose`](@ref), and [`extend`](@ref).
- [**Spatial modeling**](@ref spatial_modelling) — simulate reaction networks on discrete spatial domains.
- [**Steady state analysis**](@ref homotopy_continuation) — find and analyze steady states, stability, and bifurcation diagrams.
- [**Inverse problems**](@ref optimization_parameter_fitting) — parameter estimation, sensitivity analysis, and structural identifiability.
- [**Model I/O**](@ref model_file_import_export) — import from SBML and BioNetGen `.net` files, export to LaTeX and other formats.
- [**Visualization**](@ref visualisation_graphs) — reaction network graphs and LaTeX rendering.

## [Ecosystem](@id doc_index_ecosystem)

Catalyst integrates with a wide range of Julia packages:

| Category | Packages |
|----------|----------|
| ODE/SDE/Jump solving | [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl), [StochasticDiffEq](https://github.com/SciML/StochasticDiffEq.jl), [JumpProcesses](https://github.com/SciML/JumpProcesses.jl) |
| GPU parallelism | [DiffEqGPU](https://github.com/SciML/DiffEqGPU.jl) |
| Steady states & bifurcations | [HomotopyContinuation](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl), [SteadyStateDiffEq](https://github.com/SciML/SteadyStateDiffEq.jl), [NonlinearSolve](https://github.com/SciML/NonlinearSolve.jl), [BifurcationKit](https://github.com/bifurcationkit/BifurcationKit.jl) |
| Parameter estimation | [Optimization](https://github.com/SciML/Optimization.jl), [PEtab](https://github.com/sebapersson/PEtab.jl), [Turing](https://github.com/TuringLang/Turing.jl) |
| Sensitivity & identifiability | [GlobalSensitivity](https://github.com/SciML/GlobalSensitivity.jl), [SciMLSensitivity](https://github.com/SciML/SciMLSensitivity.jl), [StructuralIdentifiability](https://github.com/SciML/StructuralIdentifiability.jl) |
| Dynamical systems | [DynamicalSystems](https://github.com/JuliaDynamics/DynamicalSystems.jl) |
| Visualization | [Plots](https://github.com/JuliaPlots/Plots.jl), [Makie](https://github.com/MakieOrg/Makie.jl), [GraphMakie](https://github.com/MakieOrg/GraphMakie.jl), [Latexify](https://korsbo.github.io/Latexify.jl/stable/) |
| Model import | [SBMLImporter](https://github.com/sebapersson/SBMLImporter.jl), [SBMLToolkit](https://github.com/SciML/SBMLToolkit.jl), [ReactionNetworkImporters](https://github.com/SciML/ReactionNetworkImporters.jl) |
| Stochastic extensions | [MomentClosure](https://github.com/augustinas1/MomentClosure.jl), [FiniteStateProjection](https://github.com/kaandocal/FiniteStateProjection.jl), [DelaySSAToolkit](https://github.com/palmtree2013/DelaySSAToolkit.jl) |

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
jprob = JumpProblem(model, u0_integers, tspan, ps)
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
is promoted by sunlight, which is modeled as the cyclic sinusoid $k_p (\sin(t) + 1)$.
When the cell reaches a critical volume ($V_m$) it undergoes cell division.
Environmental stochasticity ($\sigma\,dW$) is added to the volume equation via
the `@brownians` DSL option. First, we declare our model:
```@example home_elaborate_example
using Catalyst
cell_model = @reaction_network begin
    @parameters Vₘ g σ
    @brownians W
    @equations begin
        D(V) ~ g*Gᴾ + σ*W
    end
    @continuous_events begin
        [V ~ Vₘ] => [V => V/2]
    end
    kₚ*(sin(t)+1)/V, G --> Gᴾ
    kᵢ/V, Gᴾ --> G
end
```
We now study the system as a Chemical Langevin Dynamics SDE model:
```@example home_elaborate_example
u0 = [:V => 25.0, :G => 50.0, :Gᴾ => 0.0]
tspan = (0.0, 20.0)
ps = [:Vₘ => 50.0, :g => 0.3, :kₚ => 100.0, :kᵢ => 60.0, :σ => 0.5]
sprob = SDEProblem(cell_model, u0, tspan, ps; mtkcompile = true)
```
This problem encodes the following stochastic differential equation model:
```math
\begin{align*}
dG(t) &=  - \left( \frac{k_p(\sin(t)+1)}{V(t)} G(t) + \frac{k_i}{V(t)} G^P(t) \right) dt - \sqrt{\frac{k_p (\sin(t)+1)}{V(t)} G(t)} \, dW_1(t) + \sqrt{\frac{k_i}{V(t)} G^P(t)} \, dW_2(t) \\
dG^P(t) &= \left( \frac{k_p(\sin(t)+1)}{V(t)} G(t) - \frac{k_i}{V(t)} G^P(t) \right) dt + \sqrt{\frac{k_p (\sin(t)+1)}{V(t)} G(t)} \, dW_1(t) - \sqrt{\frac{k_i}{V(t)} G^P(t)} \, dW_2(t) \\
dV(t) &= \left(g \, G^P(t)\right) dt + \sigma \, dW(t)
\end{align*}
```
where $dW_1(t)$ and $dW_2(t)$ are the Chemical Langevin Equation noise terms from the reactions, and $dW(t)$ is an independent Brownian motion representing environmental stochasticity. Finally, we can simulate and plot the results:
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
