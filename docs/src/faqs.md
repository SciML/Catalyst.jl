# FAQs

## How to index solution objects using symbolic variables and observables?
One can directly use symbolic variables to index into SciML solution objects.
Moreover, observables can also be evaluated in this way. For example,
consider the system
```@example faq1
using Catalyst, OrdinaryDiffEqTsit5, Plots
rn = @reaction_network ABtoC begin
  (k₊,k₋), A + B <--> C
end
nothing    # hide
```
Let's convert it to a system of ODEs, using the conservation laws of the system
to eliminate two of the species:
```@example faq1
osys = ode_model(rn; remove_conserved = true)
osys = complete(osys)
```
Notice the resulting ODE system has just one ODE, while algebraic observables
have been added for the two removed species (in terms of the conservation law
constants, `Γ[1]` and `Γ[2]`)
```@example faq1
observed(osys)
```
Let's solve the system and see how to index the solution using our symbolic
variables
```@example faq1
u0 = [osys.A => 1.0, osys.B => 2.0, osys.C => 0.0]
ps = [osys.k₊ => 1.0, osys.k₋ => 1.0]
oprob = ODEProblem(osys, u0, (0.0, 10.0), ps)
sol = solve(oprob, Tsit5())
```
Suppose we want to plot just species `C`, without having to know its integer
index in the unknown vector. We can do this using the symbolic variable `C`, which
we can get at in several ways
```@example faq1
sol[osys.C]
```
or
```@example faq1
@unpack C = osys
sol[C]
```
To evaluate `C` at specific times and plot it we can just do
```@example faq1
t = range(0.0, 10.0, length = 101)
plot(sol(t, idxs = C), label = "C(t)", xlabel = "t")
```
If we want to get multiple variables we can just do
```@example faq1
@unpack A, B = osys
sol(t, idxs = [A, B])
```
Plotting multiple variables using the SciML plot recipe can be achieved
like
```@example faq1
plot(sol; idxs = [A, B])
```

## [How to disable rescaling of reaction rates in rate laws?](@id faq_combinatoric_ratelaws)
As explained in the [Reaction rate laws used in simulations](@ref introduction_to_catalyst_ratelaws) section, for
a reaction such as `k, 2X --> 0`, the generated rate law will rescale the rate
constant, giving `k*X^2/2` instead of `k*X^2` for ODEs and `k*X*(X-1)/2` instead
of `k*X*(X-1)` for jumps. This can be disabled when directly `convert`ing a
[`ReactionSystem`](@ref). If `rn` is a generated [`ReactionSystem`](@ref), we can
do
```@example faq1
osys = ode_model(rn; combinatoric_ratelaws=false)
```
Disabling these rescalings should work for all conversions of `ReactionSystem`s
to other `ModelingToolkit.AbstractSystem`s.

When creating a [`ReactionSystem`](@ref) using the DSL, combinatoric rate laws can be disabled (for 
the created system, and all systems derived from it) using the `@combinatoric_ratelaws` option (providing `false` as its only input):
```@example faq1
rn = @reaction_network begin
    @combinatoric_ratelaws false
    k, 2X --> 0
end
nothing # hide
```

## How to use non-integer stoichiometric coefficients?
```@example faq2
using Catalyst
rn = @reaction_network begin
  k, 2.5*A --> 3*B
end
```
or directly via
```@example faq2
t = default_t()
@parameters k b
@species A(t) B(t) C(t) D(t)
rx1 = Reaction(k,[B,C],[B,D], [2.5,1],[3.5, 2.5])
rx2 = Reaction(2*k, [B], [D], [1], [2.5])
rx3 = Reaction(2*k, [B], [D], [2.5], [2])
@named mixedsys = ReactionSystem([rx1, rx2, rx3], t, [A, B, C, D], [k, b])
mixedsys = complete(mixedsys)
osys = ode_model(mixedsys; combinatoric_ratelaws = false)
osys = complete(osys)
```
Note, when using `ode_model(mixedsys; combinatoric_ratelaws=false)` the
`combinatoric_ratelaws=false` parameter must be passed. This is also true when
calling `ODEProblem(mixedsys,...; combinatoric_ratelaws=false)`. As described
above, this disables Catalyst's standard rescaling of reaction rates when
generating reaction rate laws, see also the [Reaction rate laws used in
simulations](@ref introduction_to_catalyst_ratelaws) section. Leaving this keyword out for systems with floating
point stoichiometry will give an error message.

For a more extensive documentation of using non-integer stoichiometric
coefficients, please see the [Symbolic Stochiometries](@ref
parametric_stoichiometry) section.

## How to set default values for initial conditions and parameters?
How to set defaults when using the `@reaction_network` macro is described in
more detail [here](@ref dsl_advanced_options_default_vals). There are several ways to do
this. Using the DSL, one can use the `@species` and `@parameters` options:
```@example faq3
using Catalyst
sir = @reaction_network sir begin
    @species S(t)=999.0 I(t)=1.0 R(t)=0.0
    @parameters β=1e-4 ν=0.01
    β, S + I --> 2I
    ν, I --> R
end
show(stdout, MIME"text/plain"(), sir) # hide
```

When directly constructing a `ReactionSystem`, we can set the symbolic values to
have the desired default values, and this will automatically be propagated
through to the equation solvers:
```@example faq3
using Catalyst, Plots, OrdinaryDiffEqTsit5
t = default_t()
@parameters β=1e-4 ν=.01
@species S(t)=999.0 I(t)=1.0 R(t)=0.0
rx1 = Reaction(β, [S, I], [I], [1,1], [2])
rx2 = Reaction(ν, [I], [R])
@named sir = ReactionSystem([rx1, rx2], t)
sir = complete(sir)
oprob = ODEProblem(sir, [], (0.0, 250.0))
sol = solve(oprob, Tsit5())
plot(sol)
```

One can also build a mapping from symbolic parameter/species to value/initial
condition and pass these to the `ReactionSystem` via the `initial_conditions` keyword
argument:
```@example faq3
@parameters β ν
@species S(t) I(t) R(t)
rx1 = Reaction(β, [S,I], [I], [1,1], [2])
rx2 = Reaction(ν, [I], [R])
ics = [β => 1e-4, ν => .01, S => 999.0, I => 1.0, R => 0.0]
@named sir = ReactionSystem([rx1, rx2], t; initial_conditions = ics)
nothing # hide
```

## How to specify initial conditions and parameters values for `ODEProblem` and other problem types?
To explicitly pass initial conditions and parameters we can use mappings from
Julia `Symbol`s corresponding to each variable/parameter to their values, or
from ModelingToolkit symbolic variables/parameters to their values. Using
`Symbol`s we have
```@example faq4
using Catalyst, OrdinaryDiffEqTsit5
rn = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end
u0 = [:S => 999.0, :I => 1.0, :R => 0.0]
p  = (:α => 1e-4, :β => .01)
op1  = ODEProblem(rn, u0, (0.0, 250.0), p)
nothing  # hide
```
while using ModelingToolkit symbolic variables we have
```@example faq4
t = default_t()
u0 = [rn.S => 999.0, rn.I => 1.0, rn.R => 0.0]
p  = (rn.α => 1e-4, rn.β => .01)
op2  = ODEProblem(rn, u0, (0.0, 250.0), p)
nothing  # hide
```

## How to include non-reaction terms in equations for a chemical species?
One method to add non-reaction terms into an ODE or algebraic equation for a
chemical species is to add a new (non-species) unknown variable that represents
those terms, let it be the rate of zero order reaction, and add a constraint
equation. I.e., to add a force of `(1 + sin(t))` to ``dA/dt`` in a system with
the reaction `k, A --> 0`, we can do
```@example faq5
using Catalyst
t = default_t()
@variables f(t)
rx1 = @reaction k, A --> 0
rx2 = @reaction $f, 0 --> A
eq = f ~ (1 + sin(t))
@named rs = ReactionSystem([rx1, rx2, eq], t)
rs = complete(rs)
osys = ode_model(rs)
```
In the final ODE model, `f` can be eliminated by using
`ModelingToolkitBase.mtkcompile`
```@example faq5
osyss = mtkcompile(osys)
full_equations(osyss)
```

## How to modify generated ODEs?
Conversion to other `ModelingToolkit.AbstractSystem`s allows the possibility to
modify the system with further terms that are difficult to encode as a chemical
reaction or a constraint equation. For example, an alternative method to the
previous question for adding a forcing function, $1 + \sin(t)$, to the ODE for
`dA/dt` is
```@example faq6
using Catalyst
rn = @reaction_network begin
    k, A --> 0
end
osys = ode_model(rn)
dAdteq = equations(osys)[1]
t      = ModelingToolkitBase.get_iv(osys)
dAdteq = Equation(dAdteq.lhs, dAdteq.rhs + 1 + sin(t))

# create a new ODESystem with the modified equation
@named osys2  = ODESystem([dAdteq], t)
```

## How to override mass action kinetics rate laws?
While generally one wants the reaction rate law to use the law of mass action,
so the reaction
```@example faq7
using Catalyst
using Latexify # hide
rn = @reaction_network begin
    k, X --> ∅
end
ode_model(rn)
latexify(rn; form = :ode, math_delimiters = true) # hide
```
occurs at the (ODE) rate ``d[X]/dt = -k[X]``, it is possible to override this by
using any of the following non-filled arrows when declaring the reaction: `<=`,
`⇐`, `⟽`, `=>`, `⇒`, `⟾`, `⇔`, `⟺`. This means that the reaction
```@example faq7
rn = @reaction_network begin
    k, X => ∅
end
ode_model(rn)
latexify(rn; form = :ode, math_delimiters = true) # hide
```
will occur at rate ``d[X]/dt = -k`` (which might become a problem since ``[X]``
will be degraded at a constant rate even when very small or equal to 0).

Note, stoichiometric coefficients are still included, i.e. the reaction
```@example faq7
rn = @reaction_network begin
    k, 2*X ⇒ ∅
end
ode_model(rn)
latexify(rn; form = :ode, math_delimiters = true) # hide
```
has rate ``d[X]/dt = -2 k``.

## [How to specify user-defined functions as reaction rates?](@id user_functions)
The reaction network DSL can "see" user-defined functions that work with
ModelingToolkit. e.g., this is should work
```@example faq8
using Catalyst
myHill(x) = 2*x^3/(x^3+1.5^3)
rn = @reaction_network begin
    myHill(X), ∅ --> X
end
```
In some cases, it may be necessary or desirable to register functions with
Symbolics.jl before their use in Catalyst, see the discussion
[here](https://symbolics.juliasymbolics.org/stable/manual/functions/).

## [How does the Catalyst DSL (`@reaction_network`) infer what different symbols represent?](@id faq_dsl_sym_inference)
When declaring a model using the Catalyst DSL, e.g.
```@example faq_dsl_inference
using Catalyst
rn = @reaction_network begin
 (p,d), 0 <--> X
end
```
Catalyst can automatically infer that `X` is a species and `p` and `d` are parameters. In total, Catalyst can infer the following quantities:
- Species (from reaction reactants).
- Parameters (from reaction rates and stoichiometries).
- (non-species) Variables (from the `@equations` option).
- Differential functions (from the `@equations` option).
- Observables (from the [`@observables` option](@ref dsl_advanced_options_observables)).
- Compound species (from the [`@compounds` option](@ref chemistry_functionality_compounds_DSL)).

Inference of species, variables, and parameters follows the following steps:
1. Every symbol [explicitly declared](@ref dsl_advanced_options_declaring_species_and_parameters) using the `@species`, `@variables`, and `@parameters` options are assigned to the corresponding category.
2. Every symbol not declared in (1) that occurs as a reaction reactant is inferred as a species.
3. Every symbol not declared in (1) or (2) that occurs in an expression provided after `@equations` is inferred as a variable.
4. Every symbol not declared in (1), (2), or (3) that occurs either as a reaction rate or stoichiometric coefficient is inferred to be a parameter.

Here, in 
```@example faq_dsl_inference
using Catalyst
rn = @reaction_network begin
    @parameters p1
    @equations V ~ X + p1
 X + V + p1 + p2, 0 --> X
end
```
`p` is first set as a parameter (as it is explicitly declared as such). Next, `X` is inferred as a species. Next, `V` is inferred as a variable. Finally, `p2` is inferred as a parameter.

Next, if any expression `D(...)` (where `...` can be anything) is encountered within the `@equations` option, `D` is inferred to be the differential with respect to the default independent variable (typically `t`). Note that using  `D` in this way, while also using it in another form (e.g. in a reaction rate) will produce an error.

Any symbol used as the left-hand side within the `@observables` option is inferred to be an observable. These are by default assumed to be *variables*. It is possible to simultaneously explicitly declare an observable using the `@species` or `@variables` options (in the former case, the observable will be treated as a species instead). Using observables within most other expressions (e.g. as a reactant) will produce an error.

Any symbol declared as a compound using the `@compound` option is automatically inferred to be a system species.

Symbols occurring within other expressions will not be inferred as anything. These must either occur in one of the forms described above (which enables Catalyst to infer what they are) or be explicitly declared. E.g. having a parameter which only occurs in an event:
```julia
using Catalyst
rn_error = @reaction_network begin
    @discrete_events 1.0 => [X ~ X + Xadd] 
 d, X --> 0
end
```
is not permitted. E.g. here `Xadd` must be explicitly declared as a parameter using `@parameters`:
```@example faq_dsl_inference
using Catalyst
rn = @reaction_network begin
    @parameters Xadd
    @discrete_events 1.0 => [X ~ X + Xadd] 
 d, X --> 0
end
```

It is possible to turn off all inference (requiring all symbols to be declared using `@parameters`, `@species`, and `@variables`) through the [`@require_declaration` option](@ref faq_require_declaration).

## [How can I turn off automatic inferring of species and parameters when using the DSL?](@id faq_require_declaration)
This option can be set using the `@require_declaration` option inside `@reaction_network`. In this case all the species, parameters, and variables in the system must be pre-declared using one of the `@species`, `@parameters`, or `@variables` macros. For more information about what is inferred automatically and not, please see the section on [`@require_declaration`](@ref dsl_advanced_options_require_dec).

```@example faq9
using Catalyst
rn = @reaction_network begin
    @require_declaration
    @species A(t) B(t)
    @parameters k1 k2
    (k1, k2), A <--> B
end
```

## [What to be aware of when using `remake` with conservation law elimination and NonlinearProblems?](@id faq_remake_nonlinprob)

When constructing `NonlinearSystem`s or `NonlinearProblem`s with `remove_conserved = true`, i.e.
```julia
# for rn a ReactionSystem
nsys = ss_ode_model(rn; remove_conserved = true)

# or 
nprob = NonlinearProblem(rn, u0, p; remove_conserved = true)
```
`remake` is currently unable to correctly update all `u0` values when the
conserved constant(s), `Γ`, are updated. As an example consider the following
```@example faq_remake
using Catalyst, NonlinearSolveFirstOrder
rn = @reaction_network begin
    (k₁,k₂), X₁ <--> X₂
    (k₃,k₄), X₁ + X₂ <--> 2X₃
end
u0 = [:X₁ => 1.0, :X₂ => 2.0, :X₃ => 3.0]
ps = [:k₁ => 0.1, :k₂ => 0.2, :k₃ => 0.3, :k₄ => 0.4]
nlsys = ss_ode_model(rn; remove_conserved = true, conseqs_remake_warn = false)
nlsys = complete(nlsys)
equations(nlsys)
```
If we generate a `NonlinearProblem` from this system the conservation constant,
`Γ[1]`, is automatically set to `X₁ + X₂ + X₃ = 6` and the initial values are
those in `u0`. i.e if
```@example faq_remake
nlprob1 = NonlinearProblem(nlsys, u0, ps)
```
then
```@example faq_remake
nlprob1[(:X₁, :X₂, :X₃)] == (1.0, 2.0, 3.0)
```
and
```@example faq_remake
nlprob1.ps[:Γ][1] == 6.0
```
If we now try to change a value of `X₁`, `X₂`, or `X₃` using `remake`, the
conserved constant will be recalculated. i.e. if
```@example faq_remake
nlprob2 = remake(nlprob1; u0 = [:X₂ => 3.0])
```
compare  
```@example faq_remake
println("Correct u0 is: ", (1.0, 3.0, 3.0), "\n", "remade value is: ", nlprob2[(:X₁, :X₂, :X₃)]) 
```
and
```@example faq_remake
println("Correct Γ is: ", 7.0, "\n", "remade value is: ", nlprob2.ps[:Γ][1])
```
However, if we try to directly change the value of `Γ` it is not always the case
that a `u0` value will correctly update so that the conservation law is
conserved. Consider
```@example faq_remake
nlprob3 = remake(nlprob1; u0 = [:X₂ => nothing], p = [:Γ => [4.0]])
```
Setting `[:X₂ => nothing]` for other problem types communicates that the
`u0` value for `X₂` should be solved for. However, if we examine the values we
find
```@example faq_remake
println("Correct u0 is: ", (1.0, 0.0, 3.0), "\n", "remade value is: ", nlprob3[(:X₁, :X₂, :X₃)]) 
```
and
```@example faq_remake
println("Correct Γ is: ", 4.0, "\n", "remade value is: ", nlprob3.ps[:Γ][1])
```
As such, the `u0` value for `X₂` has not updated, and the conservation law is
now violated by the `u0` values, i.e,
```@example faq_remake
(nlprob3[:X₁] + nlprob3[:X₂] + nlprob3[:X₃]) == nlprob3.ps[:Γ][1] 
```
Currently, the only way to avoid this issue is to manually specify updated
values for the `u0` components, which will ensure that `Γ` updates appropriately
as in the first example. i.e. we manually set `X₂` to the value it should be and
`Γ` will be updated accordingly:
```@example faq_remake
nlprob4 = remake(nlprob1; u0 = [:X₂ => 0.0])
```
so that
```@example faq_remake
println("Correct u0 is: ", (1.0, 0.0, 3.0), "\n", "remade value is: ", nlprob4[(:X₁, :X₂, :X₃)]) 
```
and
```@example faq_remake
println("Correct Γ is: ", 4.0, "\n", "remade value is: ", nlprob4.ps[:Γ][1])
```

*Summary:* it is not recommended to directly update `Γ` via `remake`, but to
instead update values of the initial guesses in `u0` to obtain a desired `Γ`. At
this time the behavior when updating `Γ` can result in `u0` values that do not
satisfy the conservation law defined by `Γ` as illustrated above. 