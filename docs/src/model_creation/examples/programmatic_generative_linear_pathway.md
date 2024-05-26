# [Programmatic, generative, modelling of a linear pathway](@id programmatic_generative_linear_pathway)
This example will show how to use programmatic, generative, modelling to model a system implicitly. I.e. rather than listing all system reactions explicitly, the reactions are implicitly generated from a simple set of rules. This example is specifically designed to show how [programmatic modelling](@ref ref) enables *generative workflows* (demonstrating one of its advantages as compared to [DSL-based modelling](@ref dsl_description)). In our example, we will model linear pathways, so we will first introduce these. Next, we will model them first using the DSL, and then using a generative programmatic workflow.

## [Linear pathways](@id programmatic_generative_linear_pathway_intro)
Linear pathways consists of a series of species ($X_0$, $X_1$, $X_2$, ..., $X_n$) where each activates the subsequent one. These are often modelled through the following reaction system:
```math
X_{i-1}/\tau,\hspace{0.33cm} ∅ \to X_{i}\\
1/\tau,\hspace{0.33cm} X_{i} \to ∅
```
for $i = 1, ..., n$, where the activation of $X_1$ depends on some input species $X_0$. 

A common use of these linear pathways is the implementation of *time delays*. Consider a species $X(t)$ which is activated by species $X_0(t)$. This can be modelled by making the production rate of $X(t)$ a function of the *time-delayed* value of $X_0(t)$:
```math
f(X_0(t-\tau)),\hspace{0.33cm} ∅ \to X
```
This is a so-called *discrete-delay* (which will generate a *delay differential equation*). However, in reality, $X(t)$ probably does not depend on only $f(X_0(t-\tau))$, but rather *a distribution of previous $X_0(t)$ values*. This can be modelled through a *distributed delay*s
```math
f(\int_{0}^{\inf} X_0(t-\tau)g(\tau) d\tau),\hspace{0.33cm} ∅ \to X
```
for some kernel $g(\tau)$. Here, a common kernel is a [gamma distribution](https://en.wikipedia.org/wiki/Gamma_distribution), which generates a gamma-distributed delay:
```math
g(\tau; \alpha, \beta) = \frac{\beta^{\alpha}\tau^{\alpha-1}}{\Gamma(\alpha)}e^{-\beta\tau}
```
When this is converted to an ODE, this generates an integro-differential equation. These (as well as the simpler delay differential equations) can be difficult to solve and analyse (especially when SDE or jump simulations are desired). Here, *the linear chain trick* can be used to instead model the delay as a linear pathway of the form described above[^1]. A result by Fargue shows that this is equivalent to a gamma-distributed delay, where $\alpha$ is equivalent to $n$ (the number of species in our linear pathway) and $\beta$ to %\tau$ (the delay length term)[^2]. While modelling time delays using the linear chain trick introduces additional system species, it is often advantageous as it enables simulations using standard ODE, SDE, and Jump methods.

## [Modelling linear pathways using the DSL](@id programmatic_generative_linear_pathway_dsl)
It is known that two linear pathways have similar delays if the following equality holds:
```math
\frac{1}{\tau_1 n_1} = \frac{1}{\tau_2 n_2}
```
However, the shape of the delay depends on the number of intermediaries ($n$). Here we wish to investigate this shape for two choices of $n$ ($n = 3$ and $n = 10$). We do so by implementing two models using the DSL, one for each $n$. 
```@example programmatic_generative_linear_pathway_dsl
using Catalyst

lp_n3 = @reaction_network begin
    1/τ, X0 --> 0
    (X0/τ, 1/τ), 0 <--> X1
    (X1/τ, 1/τ), 0 <--> X2
    (X2/τ, 1/τ), 0 <--> X3
end

lp_n10 = @reaction_network begin
    1/τ, X0 --> 0
    (X0/τ, 1/τ), 0 <--> X1
    (X1/τ, 1/τ), 0 <--> X2
    (X2/τ, 1/τ), 0 <--> X3
    (X3/τ, 1/τ), 0 <--> X4
    (X4/τ, 1/τ), 0 <--> X5
    (X5/τ, 1/τ), 0 <--> X6
    (X6/τ, 1/τ), 0 <--> X7
    (X7/τ, 1/τ), 0 <--> X8
    (X8/τ, 1/τ), 0 <--> X9
    (X9/τ, 1/τ), 0 <--> X10
end
nothing # hide
```
Next, we prepare an ODE for each model (scaling the initial concentration of $X_0$ and the value of $\tau$ appropriately for each model).
```@example programmatic_generative_linear_pathway_dsl
using OrdinaryDiffEq, Plots
u0_n3 = [:X0 => 3*1.0, :X1 => 0.0, :X2 => 0.0, :X3 => 0.0]
ps_n3 = [:τ => 1.0/3]
oprob_n3 = ODEProblem(lp_n3, u0_n3, (0.0, 5.0), ps_n3)

u0_n10 = [:X0 => 10*1.0, :X1 => 0.0, :X2 => 0.0, :X3 => 0.0, :X4 => 0.0, :X5 => 0.0, :X6 => 0.0, :X7 => 0.0, :X8 => 0.0, :X9 => 0.0, :X10 => 0.0]
ps_n10 = [:τ => 1.0/10.0]
oprob_n10 = ODEProblem(lp_n10, u0_n10, (0.0, 5.0), ps_n10)
nothing # hide
```
Finally, we plot the concentration of the final species in each linear pathway, noting that while the two pulses both peak at $t = 1.0$, their shapes depend on $n$.
```@example programmatic_generative_linear_pathway_dsl
sol_n3 = solve(oprob_n3)
sol_n10 = solve(oprob_n10)
plot(sol_n3; idxs = :X3, label = "n = 3")
plot!(sol_n10; idxs = :X10, label = "n = 10")
```

## [Modelling linear pathways using programmatic, generative, modelling](@id programmatic_generative_linear_pathway_generative)
Above, we investigated the impact of linear pathways' lengths on their behaviours. Since the models were implemented using the DSL, we had to implement a new model for each pathway (in each case writing out all reactions). Here, we will instead show how [programmatic modelling](@ref ref) can be used to generate pathways of arbitrary lengths.

First, we create a function, `generate_lp`, which creates a linear pathway model of length `n`. It utilises [*vector variables*](@ref ref) to create an arbitrary number of species, and also creates an [observable](@ref ref) for the final species of the chain.
```@example programmatic_generative_linear_pathway_generative
using Catalyst # hide
t = default_t()
@parameters τ
function generate_lp(n)
    # Creates a vector `X` with n+1 species.
    @species X(t)[1:n+1]
    @species Xend(t)

    # Generate
    #     (1) A degradation reaction for the input species.
    #     (2) Production reactions for all intermediary species.
    #     (2) Degradation reactions for all intermediary species.
    rxs = [
        Reaction(1/τ, [X[1]], []);
        [Reaction(X[i]/τ, [], [X[i+1]]) for i in 1:n];
        [Reaction(1/τ, [X[i+1]], []) for i in 1:n]     
    ]

    # Assembly and return a complete `ReactionSystem` (including an observable for the final species).
    @named lp_n = ReactionSystem(rxs, t; observed = [Xend ~ X[end]])
    return complete(lp_n)
end
nothing # hide
```
Next, we create a function that generates an `ODEProblem` (with appropriate initial conditions and parameter values) for arbitrarily lengthed linear pathway models.
```@example programmatic_generative_linear_pathway_generative
function generate_oprob(n)
    lp = generate_lp(n)
    X_init = fill(0.0, n + 1)
    X_init[1] = n * 1.0
    u0 = [lp.X => X_init]
    ps = [τ => 1.0 / n]
    return ODEProblem(lp, u0, (0.0, 5.0), ps)
end
nothing # hide
```
We can now simulate linear pathways of arbitrary lengths using a simple syntax. We use this to recreate our previous result from the DSL:
```@example programmatic_generative_linear_pathway_generative
sol_n3 = solve(generate_oprob(3))
sol_n10 = solve(generate_oprob(10))
plot(sol_n3; idxs = :Xend, label = "n = 3")
plot!(sol_n10; idxs = :Xend, label = "n = 10")
```
If we wish to investigate the behaviour of a pathway with a different length, we can easily add this to the plot
```@example programmatic_generative_linear_pathway_generative
sol_n20 = solve(generate_oprob(20))
plot!(sol_n20; idxs = :Xend, label = "n = 20")
```


---
## References
[^1]: [J. Metz, O. Diekmann *The Abstract Foundations of Linear Chain Trickery* (1991).](https://ir.cwi.nl/pub/1559/1559D.pdf)
[^2]: D. Fargue *Reductibilite des systemes hereditaires a des systemes dynamiques (regis par des equations differentielles aux derivees partielles)*, Comptes rendus de l'Académie des Sciences (1973).
[^3]: [N. Korsbo, H. Jönsson *It’s about time: Analysing simplifying assumptions for modelling multi-step pathways in systems biology*, PLoS Computational Biology (2020).](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007982)