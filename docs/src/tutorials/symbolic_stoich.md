# Parametric Stoichiometry
Catalyst supports stoichiometric coefficients that involve parameters, species, or even general expressions. In this tutorial we show several examples of how to use parametric stoichiometry, and discuss several caveats to be aware of.

## Using Symbolic Stoichiometry
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

An alternative approach to avoid the issues of using mixed floating point and integer variables is to disable the rescaling of rate laws as described in [Reaction rate laws used in simulations](@ref) section. This requires passing the `combinatoric_ratelaws=false` keyword to `convert` or to `ODEProblem` (if directly building the problem from a `ReactionSystem` instead of first converting to an `ODESystem`). For the previous example this gives the following (different) system of ODEs
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

## Gene expression with randomly produced amounts of protein
As a second example, let's build the negative feedback model from [MomentClosure.jl](https://augustinas1.github.io/MomentClosure.jl/dev/tutorials/geometric_reactions+conditional_closures/) that involves a bursty reaction that produces a random amount of protein. First, let's define our chemical species and parameters
```@example s1
@parameters k₊, k₋, kₚ, γₚ, b
@variables t, G₋(t), G₊(t), P(t)
nothing # hide
```
Here `G₋` denotes the repressed state, and `G₊` the active state where the gene can transcribe. `P` denotes the protein product of the gene. We will assume that proteins are produced in bursts that produce `m` proteins, where `m` is a (shifted) geometric random variable with mean `b`. To define `m` we must register the `Distributions.Geometric` distribution from Distributions.jl with Symbolics.jl, after which we can use it in symbolic expressions:
```@example s1
using Distributions: Geometric
@register Geometric(b)
m = rand(Geometric(1/b)) + 1
nothing # hide
```
Note, as we require the shifted geometric distribution, we add one to Distributions.jl's `Geometric` random variable (which includes zero).

We can now define our model
```@example s1
rxs = [Reaction(k₊, [G₋], [G₊]),
       Reaction(k₋*P^2, [G₊], [G₋]),
       Reaction(kₚ, [G₊], [G₊,P], [1], [1,m]),
       Reaction(γₚ, [P], nothing)]
@named burstyrn = ReactionSystem(rxs, t)
reactions(burstyrn)
show(stdout, MIME"text/plain"(), reactions(burstyrn)) # hide
```
and convert it to a jump process representation
```@example s1
jsys = convert(JumpSystem, burstyrn; combinatoric_ratelaws=false)
equations(jsys)
show(stdout, MIME"text/plain"(), equations(jsys)) # hide
```
Notice, the `equations` of `jsys` have three `MassActionJump`s for the first three reactions, and one `ConstantRateJump` for the last reaction. If we examine the `ConstantRateJump` more closely we can see the generated `rate` and `affect!` functions for the bursty reaction that makes protein
```@example s1
equations(jsys)[4].rate
show(stdout, MIME"text/plain"(), equations(jsys)[4].rate) # hide
```
```@example s1
equations(jsys)[4].affect!
show(stdout, MIME"text/plain"(), equations(jsys)[4].affect!) # hide
```
Finally, we can now simulate our jumpsystem
```@example s1
pmean = 200
bval = 70
γₚval = 1
k₋val = 0.001
k₊val = 0.05
kₚval = pmean * γₚval * (k₋val * pmean^2 + k₊val) / (k₊val * bval)
p = (k₊ => k₊val, k₋ => k₋val, kₚ => kₚval, γₚ => γₚval, b => bval)
u₀ = [G₊ => 1, G₋ => 0, P => 1]
tspan = (0., 6.0)   # time interval to solve over
dprob = DiscreteProblem(jsys, u₀, tspan, p)
jprob = JumpProblem(jsys, dprob, Direct())
sol = solve(jprob, SSAStepper())
plot(sol.t, sol[P], legend=false, xlabel="time", ylabel="P(t)")
```
To double check our results are consistent with MomentClosure.jl, let's calculate and plot the average amount of protein
```@example s1
function getmean(jprob, Nsims, tv)
    Pmean = zeros(length(tv))
    @variables t, P(t)
    for n in 1:Nsims
        sol = solve(jprob, SSAStepper())        
        Pmean .+= sol(tv, idxs=P)
    end
    Pmean ./= Nsims
end
tv = range(tspan[1],tspan[2],step=.1)
psim_mean = getmean(jprob, 20000, tv)
plot(tv, psim_mean, ylabel="average of P(t)", xlabel="time", xlim=(0.0,6.0), legend=false)
```