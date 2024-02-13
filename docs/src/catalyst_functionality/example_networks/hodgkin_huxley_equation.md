# [Hodgkin-Huxley Equation](@id hodgkin_huxley_equation)

This tutorial shows how to programmatically construct a
[Catalyst](http://docs.sciml.ai/Catalyst/stable/) [`ReactionSystem`](@ref) that
is coupled to a constraint ODE, corresponding to the [Hodgkin–Huxley
model](https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model) for an
excitable cell. The Hodgkin–Huxley model is a mathematical model that describes
how action potentials in neurons are initiated and propagated. It is a
continuous-time dynamical system given by a coupled system of nonlinear
differential equations that model the electrical characteristics of excitable
cells such as neurons and muscle cells.

We begin by importing some necessary packages.
```@example hh1
using ModelingToolkit, Catalyst, NonlinearSolve
using DifferentialEquations, Symbolics
using Plots
import Catalyst: t_nounits as t, D_nounits as D
```

We'll build a simple Hodgkin-Huxley model for a single neuron, with the voltage,
V(t), included as a constraint ODE. We first specify the transition rates for
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

We now declare the symbolic variable, `V(t)`, that will represent the
transmembrane potential

```@example hh1
@variables V(t)
nothing # hide
```

and a `ReactionSystem` that models the opening and closing of receptors

```@example hh1
hhrn = @reaction_network hhmodel begin
    (αₙ($V), βₙ($V)), n′ <--> n
    (αₘ($V), βₘ($V)), m′ <--> m
    (αₕ($V), βₕ($V)), h′ <--> h
end
nothing # hide
```

Next we create a `ModelingToolkit.ODESystem` to store the equation for `dV/dt`

```@example hh1
@parameters C=1.0 ḡNa=120.0 ḡK=36.0 ḡL=.3 ENa=45.0 EK=-82.0 EL=-59.0 I₀=0.0
I = I₀ * sin(2*pi*t/30)^2

# get the gating variables to use in the equation for dV/dt
@unpack m,n,h = hhrn

eqs = [Dₜ(V) ~ -1/C * (ḡK*n^4*(V-EK) + ḡNa*m^3*h*(V-ENa) + ḡL*(V-EL)) + I/C]
@named voltageode = ODESystem(eqs, t)
nothing # hide
```

Notice, we included an applied current, `I`, that we will use to perturb the system and create action potentials. For now we turn this off by setting its amplitude, `I₀`, to zero.

Finally, we add this ODE into the reaction model as

```@example hh1
@named hhmodel = extend(voltageode, hhrn)
nothing # hide
```

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
plot(hhsssol, idxs = V)
```

We now save this steady-state to use as the initial condition for simulating how
a resting neuron responds to an applied current.

```@example hh1
u_ss = hhsssol.u[end]
nothing # hide
```

Finally, starting from this resting state let's solve the system when the
amplitude of the stimulus is non-zero and see if we get action potentials

```@example hh1
tspan = (0.0, 50.0)
@unpack I₀ = hhmodel
oprob = ODEProblem(hhmodel, u_ss, tspan, [I₀ => 10.0])
sol = solve(oprob)
plot(sol, vars = V, legend = :outerright)
```

We observe three action potentials due to the steady applied current.