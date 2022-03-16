# Parametric Stoichiometry
Catalyst supports stoichiometric coefficients that involve parameters, species, or even general expressions. In this tutorial we show several examples of how to use parametric stoichiometry, and discuss several caveats to be aware of.

Let's first consider a simple reversible reaction where the number of reactants is a parameter and the number of products is the product of two parameters. Note, currently Catalyst's `@reaction_network` macro does not support stoichiometric parameters, so they need to be specified through the symbolic API interface:
```@example s1
using Catalyst, Latexify, DifferentialEquations, ModelingToolkit, Plots
@parameters k₊,k₋,m,n
@variables t, A(t), B(t)
rxs = [Reaction(k₊,[A],[B],[m],[m*n]),
       Reaction(k₋,[B],[A])]
@named revsys = ReactionSystem(rxs,t)
reactions(revsys)
```

Let's now convert the system to ODEs and look at the resulting equations. 
```@example s1
osys = convert(ODESystem, revsys; combinatori_ratelaws=false)
equations(osys)
show(stdout, MIME"text/plain"(), equations(osys)) # hide
```
Notice, as described in the [Reaction rate laws used in simulations](@ref) section, the default rate law involves factorials in the stoichiometric coefficients. As ModelingToolkit currently converts numeric parameters to a common type, this can lead to difficulties since the `factorial` function only accepts integer input, i.e. this doesn't work:
```@example s1
p  = (k₊ => 1.0, k₋ => 1.0, m => 2, n => 2)
u₀ = [A => 1.0, B => 1.0]
oprob = ODEProblem(osys, u₀, (0.0,1.0), p)
```
notice that
```@example s1
oprob.p
```
has converted all parameters to floating point. Calling
```@example s1
sol = solve(oprob, Tsit5())
```
will now give an error that
```julia
    MethodError: no method matching factorial(::Float64)
```

There are two aways around this problem. First, we can remake `oprob` to use the correct parameter values. This is complicated slightly as we need to know the parameter ordering used internally by ModelingToolkit. We can do this like: