# Parametric Stoichiometry
Catalyst supports stoichiometric coefficients that involve parameters, species, or even general expressions. In this tutorial we show several examples of how to use parametric stoichiometry, and discuss several caveats to be aware of.

Let's first consider a simple reversible reaction where the number of reactants is a parameter, and the number of products is the product of two parameters. Note, currently Catalyst's `@reaction_network` macro does not support symbolic stoichiometry, so the model needs to be specified through the symbolic API:
```@example s1
using Catalyst, Latexify, DifferentialEquations, ModelingToolkit, Plots
@parameters k₊,k₋,m,n
@variables t, A(t), B(t)
rxs = [Reaction(k₊,[A],[B],[m],[m*n]),
       Reaction(k₋,[B],[A])]
@named revsys = ReactionSystem(rxs,t)
reactions(revsys)
```
Let's now convert the system to ODEs and look at the resulting equations:
```@example s1
osys = convert(ODESystem, revsys)
equations(osys)
show(stdout, MIME"text/plain"(), equations(osys)) # hide
```
Notice, as described in the [Reaction rate laws used in simulations](@ref) section, the default rate laws involve factorials in the stoichiometric coefficients. As ModelingToolkit currently converts numeric parameters to a common type, this can lead to difficulties since the `factorial` function only accepts integer input, i.e. the integer parameters in `p`
```@example s1
p  = (k₊ => 1.0, k₋ => 1.0, m => 2, n => 2)
u₀ = [A => 1.0, B => 1.0]
oprob = ODEProblem(osys, u₀, (0.0,1.0), p)
```
are converted to floating point:
```@example s1
oprob.p
```
Calling
```julia
sol = solve(oprob, Tsit5())
```
will now give an error that
```julia
    MethodError: no method matching factorial(::Float64)
```

There are two ways around this problem. First, we can rebuild `oprob` to use a parameter tuple of the correct type. This is complicated slightly as we need to know the parameter ordering used internally by ModelingToolkit. A robust way to do this when the parameter ordering is not known is the following:
```@example s1
pmap = Dict(p)
pcorrect = Tuple(pmap[psym] for psym in parameters(osys))
oprob = remake(oprob, p=pcorrect)
oprob.p
```
now `oprob.p` has the correct type to use in solving the system
```@example s1
sol = solve(oprob, Tsit5())
plot(sol)
```
*Note, to allow for mixed parameter types (i.e. integers and floats in this example), it is necessary to use a tuple to store parameters.*

An alternative approach to avoid the issues of using mixed floating point and integer variables is to disable the rescaling of rate laws described in [Reaction rate laws used in simulations](@ref) section. This requires passing the `combinatoric_ratelaws=false` keyword to `convert` or to `ODEProblem` (if directly building the problem from a `ReactionSystem` instead of first converting to an `ODESystem`). For the previous example this gives the following (different) system of ODEs
```@example s1
osys = convert(ODESystem, revsys; combinatoric_ratelaws=false)
equations(osys)
show(stdout, MIME"text/plain"(), equations(osys)) # hide
```
Since we no longer have factorial functions appearing, our example will now run when ModelingToolkit converts `m` and `n` to be floating point:
```@example s1
oprob = ODEProblem(osys, u₀, (0.0,1.0), p)
sol = solve(oprob, Tsit5())
plot(sol)
```