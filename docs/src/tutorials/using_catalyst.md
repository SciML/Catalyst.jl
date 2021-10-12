# Using Catalyst
In this tutorial we'll provide an introduction to using Catalyst to specify
chemical reaction networks, and then to solve ODE, jump, and SDE models generated
from them. Let's start by using the Catalyst [`@reaction_network`](@ref) macro
to specify a simple chemical reaction network: the well-known repressilator.

We first import the basic packages we'll need:

```example ex
# If not already installed, first hit "]" within a Julia REPL. Then type:
# add Catalyst DifferentialEquations Plots Latexify

using Catalyst, DifferentialEquations, Plots, Latexify
```

We now construct the reaction network. The basic types of arrows and predefined
rate laws one can use are discussed in detail within the tutorial, [The
Reaction DSL](@ref). Here, we use a mix of first order, zero order, and repressive
Hill function rate laws. Note, $\varnothing$ corresponds to the empty state, and
is used for zeroth order production and first order degradation reactions:

```julia
repressilator = @reaction_network Repressilator begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₂
    hillr(P₂,α,K,n), ∅ --> m₃
    (δ,γ), m₁ <--> ∅
    (δ,γ), m₂ <--> ∅
    (δ,γ), m₃ <--> ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    β, m₃ --> m₃ + P₃
    μ, P₁ --> ∅
    μ, P₂ --> ∅
    μ, P₃ --> ∅
end α K n δ γ β μ
```
which gives
```
Model Repressilator with 15 equations
States (6):
  m₁(t)
  m₂(t)
  m₃(t)
  P₁(t)
  P₂(t)
  P₃(t)
Parameters (7):
  α
  K
  n
  δ
  γ
  β
  μ
```
showing that we've created a new network model named `Repressilator` with the
listed chemical species and states. [`@reaction_network`](@ref) returns a
[`ReactionSystem`](@ref), which we saved in the `repressilator` variabled. It
can be converted to a variety of other mathematical models represented as
`ModelingToolkit.AbstractSystem`s, or analyzed in various ways using Catalyst's
network API. 

To see the reactions in the system we use

```julia
reactions(repressilator)
```

which gives

```
15-element Vector{Reaction}:
 Catalyst.hillr(P₃(t), α, K, n), ∅ --> m₁
 Catalyst.hillr(P₁(t), α, K, n), ∅ --> m₂
 Catalyst.hillr(P₂(t), α, K, n), ∅ --> m₃
 δ, m₁ --> ∅
 γ, ∅ --> m₁
 δ, m₂ --> ∅
 γ, ∅ --> m₂
 δ, m₃ --> ∅
 γ, ∅ --> m₃
 β, m₁ --> m₁ + P₁
 β, m₂ --> m₂ + P₂
 β, m₃ --> m₃ + P₃
 μ, P₁ --> ∅
 μ, P₂ --> ∅
 μ, P₃ --> ∅
```

We can also use Latexify to see the corresponding reactions, which shows what
the `hillr` terms correspond to mathematically

```julia
latexify(repressilator, starred=true)
```
```math
\begin{align*}
\require{mhchem}
\ce{ \varnothing &->[\frac{\alpha K^{n}}{K^{n} + P{_3}^{n}}] m{_1}}\\
\ce{ \varnothing &->[\frac{\alpha K^{n}}{K^{n} + P{_1}^{n}}] m{_2}}\\
\ce{ \varnothing &->[\frac{\alpha K^{n}}{K^{n} + P{_2}^{n}}] m{_3}}\\
\ce{ m{_1} &<=>[\delta][\gamma] \varnothing}\\
\ce{ m{_2} &<=>[\delta][\gamma] \varnothing}\\
\ce{ m{_3} &<=>[\delta][\gamma] \varnothing}\\
\ce{ m{_1} &->[\beta] m{_1} + P{_1}}\\
\ce{ m{_2} &->[\beta] m{_2} + P{_2}}\\
\ce{ m{_3} &->[\beta] m{_3} + P{_3}}\\
\ce{ P{_1} &->[\mu] \varnothing}\\
\ce{ P{_2} &->[\mu] \varnothing}\\
\ce{ P{_3} &->[\mu] \varnothing}
\end{align*}
```

Assuming [Graphviz](https://graphviz.org/) is installed and commandline
accessible, within a Jupyter notebook we can also graph the reaction network by
```julia
g = Graph(repressilator)
```
giving

![Repressilator solution](../assets/repressilator.svg)

The network graph shows a variety of information, representing each species as a
blue node, and each reaction as an orange dot. Black arrows from species to
reactions indicate reactants, and are labelled with their input stoichiometry.
Similarly, black arrows from reactions to species indicate products, and are
labelled with their output stoichiometry. In contrast, red arrows from a species
to reactions indicate the species is used within the reactions' rate
expressions. For the repressilator, the reactions
```julia
hillr(P₃,α,K,n), ∅ --> m₁
hillr(P₁,α,K,n), ∅ --> m₂
hillr(P₂,α,K,n), ∅ --> m₃
```
have rates that depend on the proteins, and hence lead to red arrows from each
`Pᵢ`.

Note, from the REPL or scripts one can always use [`savegraph`](@ref) to save
the graph (assuming `Graphviz` is installed).

## Mass Action ODE Models
Let's now use our `ReactionSystem` to generate and solve a corresponding mass
action ODE model. We first convert the system to a `ModelingToolkit.ODESystem`
by
```julia
odesys = convert(ODESystem, repressilator)
```
We can once again use Latexify to look at the corresponding ODE model
```julia; results="hidden";
latexify(odesys)
```
```math
\begin{aligned}
\frac{dm_1(t)}{dt} =& \frac{\alpha K^{n}}{K^{n} + \left( \mathrm{P_3}\left( t \right) \right)^{n}} - \delta \mathrm{m_1}\left( t \right) + \gamma \\
\frac{dm_2(t)}{dt} =& \frac{\alpha K^{n}}{K^{n} + \left( \mathrm{P_1}\left( t \right) \right)^{n}} - \delta \mathrm{m_2}\left( t \right) + \gamma \\
\frac{dm_3(t)}{dt} =& \frac{\alpha K^{n}}{K^{n} + \left( \mathrm{P_2}\left( t \right) \right)^{n}} - \delta \mathrm{m_3}\left( t \right) + \gamma \\
\frac{dP_1(t)}{dt} =& \beta \mathrm{m_1}\left( t \right) - \mu \mathrm{P_1}\left( t \right) \\
\frac{dP_2(t)}{dt} =& \beta \mathrm{m_2}\left( t \right) - \mu \mathrm{P_2}\left( t \right) \\
\frac{dP_3(t)}{dt} =& \beta \mathrm{m_3}\left( t \right) - \mu \mathrm{P_3}\left( t \right)
\end{aligned}
```
(Note, there is currently a Latexify bug that causes different fonts to be used
for the species symbols on each side of the equations.)

Before we can solve the ODEs, we need to specify the values of the parameters in
the model, the initial condition, and the time interval to solve the model on.
To do this we need to build mappings from the parameters and the species to the
corresponding numerical values for parameters and initial conditions. We can do
this by using the `species` and `parameters` commands to get the corresponding
symbolic variables:
```julia
species(repressilator)
```
```julia
6-element Array{Term{Real},1}:
 m₁(t)
 m₂(t)
 m₃(t)
 P₁(t)
 P₂(t)
 P₃(t)
```
```julia
parameters(repressilator)
```
```julia
7-element Array{Sym{ModelingToolkit.Parameter{Real}},1}:
 α
 K
 n
 δ
 γ
 β
 μ
```
Knowing these orderings, we can create parameter and initial condition vectors,
and then set up the `ODEProblem` we want to solve:

```julia
# parameters [α,K,n,δ,γ,β,μ]
p = (.5, 40, 2, log(2)/120, 5e-3, 20*log(2)/120, log(2)/60)

# initial condition [m₁,m₂,m₃,P₁,P₂,P₃]
u₀ = [0.,0.,0.,20.,0.,0.]

# mappings from symbolic variables to their values
u₀map = species(repressilator) .=> u₀
pmap  = parameters(repressilator) .=> p

# time interval to solve on
tspan = (0., 10000.)

# create the ODEProblem we want to solve
oprob = ODEProblem(repressilator, u₀map, tspan, pmap)
```
Here 
```julia
pmap
```
gives 
```
7-element Vector{Pair{Sym{Real, Base.ImmutableDict{DataType, Any}}, B} where B}:
 α => 0.5
 K => 40
 n => 2
 δ => 0.0057762265046662105
 γ => 0.005
 β => 0.11552453009332422
 μ => 0.011552453009332421
```
a vector of `Pair`s, mapping each symbolic parameter to its numerical value.

Note, by passing `repressilator` directly to the `ODEProblem`, Catalyst has to
(internally) call `convert(ODESystem, repressilator)` again to generate the
symbolic ODEs. We could instead pass `odesys` directly like
```julia
oprob2 = ODEProblem(odesys, u₀map, tspan, pmap)
```
`oprob` and `oprob2` are functionally equivalent, each representing the same
underlying problem.

At this point we are all set to solve the ODEs. We can now use any ODE solver
from within the
[DifferentialEquations.jl](https://diffeq.sciml.ai/latest/solvers/ode_solve/)
package. We'll use the recommended default explicit solver, `Tsit5()`, and then
plot the solutions:

```julia
sol = solve(oprob, Tsit5(), saveat=10.)
plot(sol)
```
![Repressilator ODE Solutions](../assets/repressilator_odes.svg)

We see the well-known oscillatory behavior of the repressilator! For more on the
choices of ODE solvers, see the [DifferentialEquations.jl
documentation](https://diffeq.sciml.ai/dev/solvers/ode_solve/).

---

## Stochastic Simulation Algorithms (SSAs) for Stochastic Chemical Kinetics
Let's now look at a stochastic chemical kinetics model of the repressilator,
modeling it with jump processes. Here, we will construct a
[DiffEqJump](https://github.com/SciML/DiffEqJump.jl) `JumpProblem` that uses
Gillespie's `Direct` method, and then solve it to generate one realization of
the jump process:

```julia
# redefine the initial condition to be integer valued
u₀ = [0,0,0,20,0,0]

# recreate the symbolic variable to value mapping since u₀ is now integer-valued
u₀map = species(repressilator) .=> u₀

# next we create a discrete problem to encode that our species are integer valued:
dprob = DiscreteProblem(repressilator, u₀map, tspan, pmap)

# now, we create a JumpProblem, and specify Gillespie's Direct Method as the solver:
jprob = JumpProblem(repressilator, dprob, Direct(), save_positions=(false,false))

# now, let's solve and plot the jump process:
sol = solve(jprob, SSAStepper(), saveat=10.)
plot(sol)
```
![Repressilator SSA Solutions](../assets/repressilator_jumps.svg)

We see that oscillations remain, but become much noisier. Note, in constructing
the `JumpProblem` we could have used any of the SSAs that are part of DiffEqJump
instead of the `Direct` method, see the list of SSAs (i.e., constant rate jump
aggregators) in the
[documentation](https://diffeq.sciml.ai/dev/types/jump_types/#Constant-Rate-Jump-Aggregators-1).

Common questions that arise in using the DiffEqJump SSAs (i.e. Gillespie methods) 
are collated in the [DiffEqJump FAQ](https://diffeq.sciml.ai/latest/tutorials/discrete_stochastic_example/#FAQ).

---
## Chemical Langevin Equation (CLE) Stochastic Differential Equation (SDE) Models
At an intermediate physical scale between macroscopic ODE models and microscopic
stochastic chemical kinetics models lies the CLE, given by a system of SDEs that
add to each ODE above a noise term. As the repressilator has species that get
very close to zero in size, it is not a good candidate to model with the CLE
(where solutions can then go negative and become unphysical). Let's create a
simpler reaction network for a birth-death process that will stay non-negative:

```julia
bdp = @reaction_network begin
  c₁, X --> 2X
  c₂, X --> 0
  c₃, 0 --> X
end c₁ c₂ c₃
p = (1.0,2.0,50.)
u₀ = [5.]
tspan = (0.,4.)
```

The corresponding Chemical Langevin Equation SDE is then

```math
dX(t) = \left( c_1 X\left( t \right) - c_2 X\left( t \right) + c_3 \right) dt + \sqrt{c_1 X(t)} dW_1(t) - \sqrt{c_2 X(t)} dW_2(t) + \sqrt{c_3} dW_3(t)
```

where each $W_i(t)$ denotes an independent Brownian Motion. We can solve the CLE
model by creating an `SDEProblem` and solving it similarly to what we did for ODEs
above:

```julia
# SDEProblem for CLE
sprob = SDEProblem(bdp, species(bdp) .=> u₀, tspan, parameters(bdp) .=> p)

# solve and plot, tstops is used to specify enough points
# that the plot looks well-resolved
sol = solve(sprob, LambaEM(), tstops=range(0., step=4e-3, length=1001))
plot(sol)
```
![CLE Solution](../assets/birthdeath_cle.svg)

We again have complete freedom to select any of the
StochasticDiffEq.jl SDE solvers, see the
[documentation](https://diffeq.sciml.ai/dev/solvers/sde_solve/).

---
## Notes
1. For each of the preceding models we converted the `ReactionSystem` to, i.e.,
   ODEs, jumps, or SDEs, we had two paths for conversion:

    a. Convert to the corresponding ModelingToolkit system and then use it in
       creating the corresponding problem.

    b. Directly create the desired problem type from the `ReactionSystem`.

   The latter is more convenient, however, the former will be more efficient if
   one needs to repeatedly create the associated `Problem`.
2. ModelingToolkit offers many options for optimizing the generated ODEs and
   SDEs, including options to build functions for evaluating Jacobians and/or
   multithreaded versions of derivative evaluation functions. See the options
   for
   [`ODEProblem`s](https://mtk.sciml.ai/dev/systems/ODESystem/#DiffEqBase.ODEProblem)
   and
   [`SDEProblem`s](https://mtk.sciml.ai/dev/systems/SDESystem/#DiffEqBase.SDEProblem).
