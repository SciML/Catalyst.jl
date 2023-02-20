# [Finding Steady States through Homotopy Continuation](@id homotopy_continuation)

The steady states of a dynamical system ${dx \over dt} = f(x)$ can be found by
solving $0 = f(x)$. This is typically a hard problem, and generally, there is no
method guaranteed to find all steady states for a system that has multiple ones.
However, many chemical reaction networks generate polynomial systems (for
example those which are purely mass action or have only have Hill functions with
integer Hill exponents). The roots of these can reliably be found through a
*homotopy continuation* algorithm. This is implemented in Julia through the
[HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/) package.
In this tutorial, we will demonstrate how homotopy continuation can be used to
find the steady states of mass action chemical reaction networks implemented in
Catalyst.

## Basic example
For this tutorial, we will use a model from Wilhem (2009)[^1] (which
demonstrates bistability in a small chemical reaction network). We declare the
model and the parameter set for which we want to find the steady states:
```@example hc1
using Catalyst, ModelingToolkit
import HomotopyContinuation
const MT = ModelingToolkit
const HC = HomotopyContinuation

wilhelm_2009_model = @reaction_network begin
    k1, Y --> 2X
    k2, 2X --> X + Y
    k3, X + Y --> Y
    k4, X --> 0
end

# add default parameters values to model
setdefaults!(wilhelm_2009_model, [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5])
nothing   # hide
```
Next, we will need to extract the actual equations from our model. In addition,
we will substitute in our parameter values to these equations.
```@example hc1
ns = convert(NonlinearSystem, wilhelm_2009_model)

# this gets the parameter values ordered consistent with parameters(ns)
pvals = MT.varmap_to_vars([], MT.parameters(ns); defaults = MT.defaults(ns))

subs = Dict(MT.parameters(ns) .=> pvals)
neweqs = map(eq -> substitute(eq.rhs, subs), equations(ns))
```
Finally, we use Catalyst's `to_multivariate_poly` function to reinterpret our
symbolic equations in a polynomial representation that is compatible with
HomotopyContinuation. We can then apply HomotopyContinuation's `solve` command
to find the roots, using `real_solutions` to filter our non-physical complex
steady-states:
```@example hc1
polyeqs = Catalyst.to_multivariate_poly(neweqs)
sols = HC.real_solutions(HC.solve(polyeqs))
```
Note that HomotopyContinuation orders variables lexicographically, so this will
be the ordering present in each steady-state solution vector (i.e. `[X1, X2]` is
the ordering here).

While it is not the case for this CRN, we note that solutions with negative
species concentrations can be valid (unphysical) steady-states for certain
systems. These will need to be filtered out as well.

## Systems with conservation laws
Finally, some systems are under-determined, and have an infinite number of
possible steady states. These are typically systems containing a conservation
law, e.g.
```@example hc3
using Catalyst
import HomotopyContinuation
const MT = ModelingToolkit
const HC = HomotopyContinuation

two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
```
Catalyst allows the conservation laws to be computed using the
`conservationlaws` function. By using these to reduce the dimensionality of the
system, as well specifying the initial amount of each species,
HomotopyContinuation can again be used to find steady states. First, we set the
default values of the system's initial conditions and parameter values. This
will allow the system to automatically find the conserved amounts.
```@example hc3
setdefaults!(two_state_model, [:X1 => 1.0, :X2 => 1.0, :k1 => 2.0, :k2 => 1.0])
```
Next, we create a `NonlinearSystem`, while also removing one species via the
conservation equation.
```@example hc3
ns = convert(NonlinearSystem, two_state_model; remove_conserved = true)
```
Again, we next create a dictionary for parameter values that we substitute in to
give our final equation.
```@example hc3
pvals = MT.varmap_to_vars([], MT.parameters(ns); defaults = MT.defaults(ns))
subs = Dict(MT.parameters(ns) .=> pvals)
neweqs = map(eq -> substitute(eq.rhs, subs), equations(ns))
```
Notice, our equations are just for `X1` as `X2` was eliminated via the conservation law.

Finally, we convert to polynomial form and solve for the steady-states
```@example hc3
polyeqs = Catalyst.to_multivariate_poly(neweqs)
sols = HC.real_solutions(HC.solve(polyeqs))
```

If we also want the corresponding value for `X2`, we can substitute
into the equation for it from the conservation laws:
```@example hc3
# get the X2 symbolic variable
@unpack X2 = two_state_model

# get its algebraic formula in terms of X1 and parameters
ceqs = conservedequations(two_state_model)
X2eqidx = findfirst(eq -> isequal(eq.lhs, X2), ceqs)
X2eq = ceqs[X2eqidx].rhs

# for each SS, set X1's value in the subs map and calculate X2
@unpack X1 = two_state_model
X2 = map(sols) do x
    X1val = x[1]
    subs[MT.value(X1)] = X1val
    substitute(X2eq, subs)
end
```
giving that the steady-state for `X2` is about `1.33333`.

As an alternative, we could have coupled `neweqs` with the conservation law
relations to have HomotopyContinuation find the steady-state simultaneously:
```@example hc3
# move all the terms in the conserved equations to one side
# and substitute in the parameter values
subs = Dict(MT.parameters(ns) .=> pvals)
conservedrelations = map(eq -> substitute(eq.rhs - eq.lhs, subs), ceqs)
neweqs = vcat(neweqs, conservedrelations)

# calculate the steady-states
polyeqs = Catalyst.to_multivariate_poly(neweqs)
sols = HC.real_solutions(HC.solve(polyeqs))
```
---
## References
[^1]: [Wilhelm, T. *The smallest chemical reaction system with bistability*, BMC Systems Biology (2009).](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-050Wilhelm-3-90)
