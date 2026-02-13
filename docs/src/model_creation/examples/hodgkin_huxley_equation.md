# [Hodgkin-Huxley Equation](@id hodgkin_huxley_equation)

This tutorial shows how to construct a
[Catalyst](http://docs.sciml.ai/Catalyst/stable/) [`ReactionSystem`](@ref) that
includes a coupled ODE, corresponding to the [Hodgkin–Huxley
model](https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model) for an
excitable cell. The Hodgkin–Huxley model is a mathematical model that describes
how action potentials in neurons are initiated and propagated. It is a
continuous-time dynamical system given by a coupled system of nonlinear
differential equations that model the electrical characteristics of excitable
cells such as neurons and muscle cells.

We'll present two different approaches for constructing the model. The first
will show how it can be built entirely within a single DSL model, while the
second illustrates another work flow, showing how separate models containing the
chemistry and the dynamics of the transmembrane potential can be combined into a
complete model.

We begin by importing some necessary packages:
```@example hh1
using ModelingToolkitBase, Catalyst, NonlinearSolveFirstOrder, Plots, OrdinaryDiffEqRosenbrock
```

## Building the model via the Catalyst DSL
Let's build a simple Hodgkin-Huxley model for a single neuron, with the voltage,
$V(t)$, included as a coupled ODE. We first specify the transition rates for
three gating variables, $m(t)$, $n(t)$ and $h(t)$.

$$s' \xleftrightarrow[\beta_s(V(t))]{\alpha_s(V(t))} s, \quad s \in \{m,n,h\}$$

Here each gating variable is used in determining the fraction of active (i.e.
open) or inactive ($m' = 1 - m$, $n' = 1 -n$, $h' = 1 - h$) sodium ($m$ and $h$)
and potassium ($n$) channels.

The transition rate functions, which depend on the voltage, $V(t)$, are given by

```@example hh1
function αₘ(V)
    theta = (V + 45) / 10
    ifelse(theta == 0.0, 1.0, theta/(1 - exp(-theta)))
end
βₘ(V) = 4*exp(-(V + 70)/18)

αₕ(V) = .07 * exp(-(V + 70)/20)
βₕ(V) = 1/(1 + exp(-(V + 40)/10))

function αₙ(V)
    theta = (V + 60) / 10
    ifelse(theta == 0.0, .1, .1*theta / (1 - exp(-theta)))
end
βₙ(V) = .125 * exp(-(V + 70)/80)
nothing # hide
```

We also declare a function to represent an applied current in our model, which we
will use to perturb the system and create action potentials. 
```@example hh1
Iapp(t,I₀) = I₀ * sin(2*pi*t/30)^2
```

We now declare a `ReactionSystem` that encompasses the Hodgkin-Huxley model.
Note, we will also give the (default) values for our parameters as part of
constructing the model to avoid having to specify them later on via parameter
maps.

```@example hh1
hhmodel = @reaction_network hhmodel begin
    @parameters begin
        C = 1.0 
        ḡNa = 120.0 
        ḡK = 36.0 
        ḡL = .3 
        ENa = 45.0 
        EK = -82.0 
        EL = -59.0 
        I₀ = 0.0
    end

    @variables V(t)

    (αₙ(V), βₙ(V)), n′ <--> n
    (αₘ(V), βₘ(V)), m′ <--> m
    (αₕ(V), βₕ(V)), h′ <--> h
    
    @equations begin
        D(V) ~ -1/C * (ḡK*n^4*(V-EK) + ḡNa*m^3*h*(V-ENa) + ḡL*(V-EL)) + Iapp(t,I₀)
    end
end
```
For now we turn off the applied current by setting its amplitude, `I₀`, to zero.

`hhmodel` is now a `ReactionSystem` that is coupled to an internal constraint
ODE for $dV/dt$. Let's now solve to steady-state, as we can then use these
resting values as an initial condition before applying a current to create an
action potential.

```@example hh1
tspan = (0.0, 50.0)
u₀ = [:V => -70, :m => 0.0, :h => 0.0, :n => 0.0,
	  :m′ => 1.0, :n′ => 1.0, :h′ => 1.0]
oprob = ODEProblem(hhmodel, u₀, tspan)
hhsssol = solve(oprob, Rosenbrock23())
nothing # hide
```

From the artificial initial condition we specified, the solution approaches the
physiological steady-state via firing one action potential:

```@example hh1
plot(hhsssol, idxs = hhmodel.V)
```

We now save this steady-state to use as the initial condition for simulating how
a resting neuron responds to an applied current. We save the steady-state values
as a mapping from the symbolic variables to their steady-states that we can
later use as an initial condition:

```@example hh1
u_ss = unknowns(hhmodel) .=> hhsssol(tspan[2], idxs = unknowns(hhmodel))
nothing # hide
```

Finally, starting from this resting state let's solve the system when the
amplitude of the stimulus is non-zero and see if we get action potentials

```@example hh1
tspan = (0.0, 50.0)
@unpack V,I₀ = hhmodel
oprob = ODEProblem(hhmodel, u_ss, tspan, [I₀ => 10.0])
sol = solve(oprob)
plot(sol, idxs = V, legend = :outerright)
```

We observe three action potentials due to the steady applied current.

## Building the model via composition of separate systems for the ion channel and transmembrane voltage dynamics 

As an illustration of how one can construct models from individual components,
we now separately construct and compose the model components.

We start by defining systems to model each ionic current. Note we now use
`@network_component` instead of `@reaction_network` as we want the models to be
composable and not marked as finalized.
```@example hh1
IKmodel = @network_component IKmodel begin
    @parameters ḡK = 36.0 EK = -82.0 
    @variables V(t) Iₖ(t)
    (αₙ(V), βₙ(V)), n′ <--> n
    @equations Iₖ ~ ḡK*n^4*(V-EK)
end

INamodel = @network_component INamodel begin
    @parameters ḡNa = 120.0 ENa = 45.0 
    @variables V(t) Iₙₐ(t)
    (αₘ(V), βₘ(V)), m′ <--> m
    (αₕ(V), βₕ(V)), h′ <--> h
    @equations Iₙₐ ~ ḡNa*m^3*h*(V-ENa) 
end

ILmodel = @network_component ILmodel begin
    @parameters ḡL = .3 EL = -59.0 
    @variables V(t) Iₗ(t)
    @equations Iₗ ~ ḡL*(V-EL)
end
nothing # hide
```

We next define the voltage dynamics with unspecified values for the currents
```@example hh1
hhmodel2 = @network_component hhmodel2 begin
    @parameters C = 1.0 I₀ = 0.0
    @variables V(t) Iₖ(t) Iₙₐ(t) Iₗ(t)
    @equations D(V) ~ -1/C * (Iₖ + Iₙₐ + Iₗ) + Iapp(t,I₀)
end
nothing # hide
```
Finally, we extend the `hhmodel` with the systems defining the ion channel currents
```@example hh1
@named hhmodel2 = extend(IKmodel, hhmodel2)
@named hhmodel2 = extend(INamodel, hhmodel2)
@named hhmodel2 = extend(ILmodel, hhmodel2)
hhmodel2 = complete(hhmodel2)
```
Let's again solve the system starting from the previously calculated resting
state, using the same applied current as above (to verify we get the same
figure). Note, we now run `structural_simplify` from ModelingToolkit to
eliminate the algebraic equations for the ionic currents when constructing the
`ODEProblem`:

```@example hh1
@unpack I₀,V = hhmodel2
oprob = ODEProblem(hhmodel2, u_ss, tspan, [I₀ => 10.0]; structural_simplify = true)
sol = solve(oprob)
plot(sol, idxs = V, legend = :outerright)
```

We observe the same solutions as from our original model.
