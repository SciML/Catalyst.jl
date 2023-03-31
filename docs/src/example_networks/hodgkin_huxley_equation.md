# [Hodgkin-Huxley Equation](@id hodgkin_huxley_equation)

This tutorial shows how to programmatically construct a [`ReactionSystem`](@ref) corresponding to the chemistry underlying the [Hodgkin–Huxley model](https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model) using [ModelingToolkit](http://docs.sciml.ai/ModelingToolkit/stable/)/[Catalyst](http://docs.sciml.ai/Catalyst/stable/).

The Hodgkin–Huxley model, or conductance-based model, is a mathematical model that describes how action potentials in neurons are initiated and propagated. It is a set of nonlinear differential equations that approximates the electrical characteristics of excitable cells such as neurons and muscle cells. It is a continuous-time dynamical system.

We begin by importing some necessary packages.
```julia
using ModelingToolkit, Catalyst, NonlinearSolve
using DifferentialEquations, IfElse
using Plots, GraphRecipes
```

Let's build a simple Hodgkin-Huxley model for a single neuron, with the voltage, V(t), included as a constraint ODESystem.

We first specify the transition rates for the three gating variables, $m(t)$, $n(t)$ and $h(t)$.

$$s' \xleftrightarrow[\beta_s(V(t))]{\alpha_s(V(t))} s, \quad s \in \{m,n,h\}$$

where $m$, $n$ and $h$, are gating variables that determine the fraction of active(open) or inactive ($m' = 1 - m, n' = 1 -n, h' = 1 - h$) receptors.

The transition rate functions, which depend on the voltage, $V(t)$, are then

```julia
begin 
	function αₘ(V) 
		theta = (V + 45) / 10
		IfElse.ifelse(theta == 0.0, 1.0, theta/(1 - exp(-theta)))
	end
	βₘ(V) = 4*exp(-(V + 70)/18)
	
	αₕ(V) = .07 * exp(-(V + 70)/20)
	βₕ(V) = 1/(1 + exp(-(V + 40)/10))
	
	function αₙ(V)
		theta = (V + 60) / 10
		IfElse.ifelse(theta == 0.0, .1, .1*theta / (1 - exp(-theta)))
	end
	βₙ(V) = .125 * exp(-(V + 70)/80)
end
```
* We now declare the symbolic variable, `V(t)`, that will represent voltage.

* We tell Catalyst not to generate an equation for it from the reactions we list, using the `isbcspecies` metadata.

* This label tells Catalyst an ODE or nonlinear equation for `V(t)` will be provided in a constraint system.

Aside: `bcspecies` means a boundary condition species, a terminology from SBML.

$\begin{equation}V\left( t \right)\end{equation}$

```julia
@variables V(t) [isbcspecies=true]
```


