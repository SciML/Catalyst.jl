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

We first specify the transition rates for the three gating variables, $$m(t)$$, $$n(t)$$, $$h(t)$$

$$s' \xleftrightarrow[\beta_s(V(t))]{\alpha_s(V(t))} s, \quad s \in \{m,n,h\}$$