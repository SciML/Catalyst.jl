# [Symbolic Stochiometries](@id parametric_stoichiometry)
Catalyst supports stoichiometric coefficients that involve parameters, species,
or even general expressions. In this tutorial we show several examples of how to
use symbolic stoichiometries, and discuss several caveats to be aware of.

*Note, this tutorial requires ModelingToolkit v8.5.4 or greater to work properly.*

## Using symbolic stoichiometry
Let's first consider a simple reversible reaction where the number of reactants
is a parameter, and the number of products is the product of two parameters.
```@example s1
using Catalyst, Latexify, DifferentialEquations, ModelingToolkit, Plots
revsys = @reaction_network revsys begin
    k₊, m*A --> (m*n)*B
    k₋, B --> A
end
reactions(revsys)
```
Note, as always the `@reaction_network` macro defaults to setting all symbols
neither used as a reaction substrate nor a product to be parameters. Hence, in
this example we have two species (`A` and `B`) and four parameters (`k₊`, `k₋`,
`m`, and `n`). In addition, the stoichiometry is applied to the rightmost symbol
in a given term, i.e. in the first equation the substrate `A` has stoichiometry
`m` and the product `B` has stoichiometry `m*n`. For example, in
```@example s1
rn = @reaction_network begin
    k, A*C --> 2B
    end
reactions(rn)
```
we see two species, `(B,C)`, with `A` treated as a parameter representing the
stoichiometric coefficient of `C`, i.e.
```@example s1
rx = reactions(rn)[1]
rx.substrates[1],rx.substoich[1]
```
We could have equivalently specified our systems directly via the Catalyst
API. For example, for `revsys` we would could use
```@example s1
@parameters k₊, k₋, m, n
@variables t
@species A(t), B(t)
rxs = [Reaction(k₊, [A], [B], [m], [m*n]),
       Reaction(k₋, [B], [A])]
revsys2 = ReactionSystem(rxs,t; name=:revsys)
revsys2 == revsys
```
which can be simplified using the `@reaction` macro to
```@example s1
rxs2 = [(@reaction k₊, m*A --> (m*n)*B),
        (@reaction k₋, B --> A)]
revsys3 = ReactionSystem(rxs2,t; name=:revsys)
revsys3 == revsys
```
Note, the `@reaction` macro again assumes all symbols are parameters except the
substrates or reactants (i.e. `A` and `B`). For example, in
`@reaction k, F*A + 2(H*G+B) --> D`, the substrates are `(A,G,B)` with
stoichiometries `(F,2*H,2)`.

Let's now convert `revsys` to ODEs and look at the resulting equations:
```@example s1
osys = convert(ODESystem, revsys)
equations(osys)
show(stdout, MIME"text/plain"(), equations(osys)) # hide
```
Notice, as described in the [Reaction rate laws used in simulations](@ref)
section, the default rate laws involve factorials in the stoichiometric
coefficients. For this reason we must specify `m` and `n` as integers, and hence
*use a tuple for the parameter mapping*
```@example s1
p  = (k₊ => 1.0, k₋ => 1.0, m => 2, n => 2)
u₀ = [A => 1.0, B => 1.0]
oprob = ODEProblem(osys, u₀, (0.0, 1.0), p)
nothing # hide
```
We can now solve and plot the system
```@example s1
sol = solve(oprob, Tsit5())
plot(sol)
```
*If we had used a vector to store parameters, `m` and `n` would be converted to
floating point giving an error when solving the system.*

An alternative approach to avoid the issues of using mixed floating point and
integer variables is to disable the rescaling of rate laws as described in
[Reaction rate laws used in simulations](@ref) section. This requires passing
the `combinatoric_ratelaws=false` keyword to `convert` or to `ODEProblem` (if
directly building the problem from a `ReactionSystem` instead of first
converting to an `ODESystem`). For the previous example this gives the following
(different) system of ODEs
```@example s1
osys = convert(ODESystem, revsys; combinatoric_ratelaws = false)
equations(osys)
show(stdout, MIME"text/plain"(), equations(osys)) # hide
```
Since we no longer have factorial functions appearing, our example will now run
even with floating point values for `m` and `n`:
```@example s1
p  = (k₊ => 1.0, k₋ => 1.0, m => 2.0, n => 2.0)
oprob = ODEProblem(osys, u₀, (0.0, 1.0), p)
sol = solve(oprob, Tsit5())
plot(sol)
```

## Gene expression with randomly produced amounts of protein
As a second example, let's build the negative feedback model from
[MomentClosure.jl](https://augustinas1.github.io/MomentClosure.jl/dev/tutorials/geometric_reactions+conditional_closures/)
that involves a bursty reaction that produces a random amount of protein.

In our model `G₋` will denote the repressed state, and `G₊` the active state
where the gene can transcribe. `P` will denote the protein product of the gene.
We will assume that proteins are produced in bursts that produce `m` proteins,
where `m` is a (shifted) geometric random variable with mean `b`. To define `m`
we must register the `Distributions.Geometric` distribution from
Distributions.jl with Symbolics.jl, after which we can use it in symbolic
expressions:
```@example s1
using Distributions: Geometric
@register_symbolic Geometric(b)
@parameters b
m = rand(Geometric(1/b)) + 1
nothing # hide
```
Note, as we require the shifted geometric distribution, we add one to
Distributions.jl's `Geometric` random variable (which includes zero).

We can now define our model
```@example s1
burstyrn = @reaction_network burstyrn begin
    k₊, G₋ --> G₊
    k₋*P^2, G₊ --> G₋
    kₚ, G₊ --> G₊ + $m*P
    γₚ, P --> ∅
end
reactions(burstyrn)
show(stdout, MIME"text/plain"(), reactions(burstyrn)) # hide
```
The parameter `b` does not need to be explicitly declared in the
`@reaction_network` macro as it is detected when the expression
`rand(Geometric(1/b)) + 1` is substituted for `m`.

We next convert our network to a jump process representation
```@example s1
jsys = convert(JumpSystem, burstyrn; combinatoric_ratelaws = false)
equations(jsys)
show(stdout, MIME"text/plain"(), equations(jsys)) # hide
```
Notice, the `equations` of `jsys` have three `MassActionJump`s for the first
three reactions, and one `ConstantRateJump` for the last reaction. If we examine
the `ConstantRateJump` more closely we can see the generated `rate` and
`affect!` functions for the bursty reaction that makes protein
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
p = symmap_to_varmap(jsys, (:k₊ => k₊val, :k₋ => k₋val, :kₚ => kₚval,
                            :γₚ => γₚval, :b => bval))
u₀ = symmap_to_varmap(jsys, [:G₊ => 1, :G₋ => 0, :P => 1])
tspan = (0., 6.0)   # time interval to solve over
dprob = DiscreteProblem(jsys, u₀, tspan, p)
jprob = JumpProblem(jsys, dprob, Direct())
sol = solve(jprob, SSAStepper())
plot(sol.t, sol[jsys.P], legend = false, xlabel = "time", ylabel = "P(t)")
```
To double check our results are consistent with MomentClosure.jl, let's
calculate and plot the average amount of protein (which is also plotted in the
MomentClosure.jl
[tutorial](https://augustinas1.github.io/MomentClosure.jl/dev/tutorials/geometric_reactions+conditional_closures/)).
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
plot(tv, psim_mean; ylabel = "average of P(t)", xlabel = "time",
                    xlim = (0.0,6.0), legend = false)
```
Comparing, we see similar averages for `P(t)`.