# [Finding Steady States through Homotopy Continuation](@id homotopy_continuation)

The steady states of a dynamical system ${dx \over dt} = f(x)$ can be found by
solving $0 = f(x)$. This is typically a hard problem, and generally, there is no
method guaranteed to find all steady states for a system that has multiple ones.
However, most chemical reaction networks generate polynomial systems (the main
exception is when Hill functions with non-integer exponents are used). The roots
of these can reliably be found through a *homotopy continuation* algorithm. This
is implemented in Julia through the
[HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/) package.
In this tutorial, we will demonstrate how homotopy continuation can be used to
find the steady states of a chemical reaction network implemented in  Catalyst.

## Basic example
For this tutorial, we will use a model from Wilhem (2009)[^1] (which
demonstrates bistability in a small chemical reaction network). We declare the
model and the parameter set for which we want to find the steady states:
```@example hc1
using Catalyst, ModelingToolkit
const MT = ModelingToolkit

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
new_eqs = map(eq -> substitute(eq.rhs, subs), equations(ns))

```
Finally, we use the `as_polynomial` function to read our symbolic expression as
a polynomial, within it, we can apply homotopy continuation's `solve` command to
find the roots. In addition, we use the `real_solutions` to filter away
imaginary roots (as CRNs' states typically are non-imaginary):
```@example hc1
using HomotopyContinuation
sols = real_solutions(as_polynomial((f, x...) -> HomotopyContinuation.solve(collect(x)), new_eqs...))
```
While it is not the case for this CRN, we note that some solutions with negative
species concentrations may still appear. Typically, these will need to be
filtered away as well.

## Rational polynomial systems
It is not uncommon for CRNs to generate systems corresponding to rational
multivariate polynomials (e.g. through Hill functions). The roots of these can
also be found using homotopy continuation. An expanded tutorial for this will be
published once some awaited improvements to the `as_polynomial` function are
completed.

## Systems with conservation laws
Finally, some systems are underdetermined, and have an infinite number of
possible steady states. These are typically systems containing a conservation
law, e.g.
```@example hc3
using Catalyst
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
```
However, the conservation laws can be computed using the `conservationlaws`
function. By supplying these, as well as fixed concentrations (in this case the
total amount of *X*, that is *X1+X2*), steady states can be found. First, we set
the default values of the system's initial conditions and parameter values. This
will allow the system to automatically find the conserved amounts.
```@example hc3
setdefaults!(two_state_model, [:X1 => 1.0, :X2 => 1.0, :k1 => 2.0, :k2 => 1.0])
```
Next, we create a `NonlinearSystem`, while also removing the redundant equation.
```@example hc3
ns = convert(NonlinearSystem, two_state_model; remove_conserved=true)
```
Again, we will create the dictionary for parameter values that we will sub in.
However, we will do it slightly differently so that the conserved quantities are
accounted for.
```@example hc3
const MT = ModelingToolkit
subs = Dict(MT.parameters(ns) .=> MT.varmap_to_vars([], MT.parameters(ns); defaults=MT.defaults(ns)))
```
We now extract the equation produced by the conservation law, and then sub in
the parameter values creating a final set of equations (like previously). Unlike
previously, we have to do `eq.rhs-eq.lhs`, as `cons_eq` may contain important
information on both the lhs and rhs.
```@example hc3
cons_eq = conservedequations(two_state_model)
new_eqs = map(eq -> substitute(eq.rhs - eq.lhs, subs), [equations(ns)...; cons_eq...])
```
Finally, we compute the solution:
```@example hc3
using HomotopyContinuation
sols = real_solutions(as_polynomial((f, x...) -> HomotopyContinuation.solve(collect(x)), new_eqs...))
```

---
## References
[^1]: [Wilhelm, T. *The smallest chemical reaction system with bistability*, BMC Systems Biology (2009).](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-050Wilhelm-3-90)
