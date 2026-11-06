# [Common Mistakes](@id common_mistakes)
New Catalyst users often run into a handful of recurring pitfalls. This page
collects the most common ones and shows how to avoid them. If you are just
getting started, skim this page after the introductory tutorials; if something
mysteriously "does not connect", check here first.

## [Declaring `@species` / `@parameters` outside `@reaction_network`](@id common_mistakes_outside_declarations)
A frequent misunderstanding is that declaring symbolic species or parameters
*outside* of `@reaction_network` somehow registers them into the next network
you create. It does not. The DSL creates its own symbolic variables from the
names that appear in the reaction expressions. Metadata such as default values
attached to the *outside* symbols is **not** copied into the network unless you
either:

1. declare them *inside* the DSL with `@species` / `@parameters`, or
2. interpolate the outside symbols with `$` (see
   [Interpolating variables into the DSL](@ref dsl_advanced_options_symbolics_and_DSL_interpolation)).

### What goes wrong
Here we attach defaults to outside symbols, then build a network that *mentions*
the same names but does not interpolate them:
```@example common_mistakes_outside
using Catalyst
t = default_t()
@species X(t)=5.0
@parameters p=1.0 d=0.1

rn_plain = @reaction_network begin
    (p, d), 0 <--> X
end
nothing # hide
```
Those defaults live on the outside Julia variables. The network itself has no
stored defaults for `X`, `p`, or `d`, so building an `ODEProblem` without
providing values fails:
```@example common_mistakes_outside
using OrdinaryDiffEq # hide
try
    ODEProblem(rn_plain, [], (0.0, 1.0), [])
catch e
    println(sprint(showerror, e)[1:min(end, 180)], "…")
end
```

### What to do instead
Prefer declaring defaults **inside** the DSL (usually the simplest option):
```@example common_mistakes_outside
rn_inside = @reaction_network begin
    @species X(t)=5.0
    @parameters p=1.0 d=0.1
    (p, d), 0 <--> X
end
oprob_inside = ODEProblem(rn_inside, [], (0.0, 1.0), [])
oprob_inside.u0
```
Or interpolate the outside symbols so their metadata is carried through:
```@example common_mistakes_outside
rn_interp = @reaction_network begin
    ($p, $d), 0 <--> $X
end
oprob_interp = ODEProblem(rn_interp, [], (0.0, 1.0), [])
oprob_interp.u0
```
!!! tip
    If you do not need programmatic construction or shared symbolic
    expressions, you usually do **not** need `@species` / `@parameters` outside
    the DSL at all. Let `@reaction_network` infer names, or declare them inside
    the `begin ... end` block.

## [Defining `t` outside `@reaction_network`](@id common_mistakes_independent_variable)
Related to the previous point: creating `t = default_t()` (or `@variables t`)
outside a network does not by itself wire that independent variable into the
DSL. `@reaction_network` always introduces its own independent variable (also
named `t` by default).

In practice `isequal` often holds between an outside `t` and the network's
independent variable, so indexing and many workflows still work. The painful
part is again **metadata and identity**: `@species X(t)=…` declared with an
outside `t` will not magically give the network default values unless you
interpolate `$X` (or declare `X` inside the DSL).

```@example common_mistakes_t
using Catalyst
t = default_t()
@species Y(t)=1.0
@parameters k=1.0

rn_no_dollar = @reaction_network begin
    k, 0 --> Y
end

rn_with_dollar = @reaction_network begin
    $k, 0 --> $Y
end

using OrdinaryDiffEq # hide
oprob = ODEProblem(rn_with_dollar, [], (0.0, 1.0), [])
oprob.u0
```
!!! note
    For ordinary DSL models you can ignore outside `t` entirely. Use
    `default_t()` primarily when you are building models
    [programmatically](@ref programmatic_CRN_construction) or interpolating
    pre-built symbolic expressions.

## [Interpolation: when you need `$` (and when you do not)](@id common_mistakes_interpolation)
Catalyst's `$` interpolation splices an existing Julia/symbolic value into the
DSL. It is powerful, but easy to misuse.

### Interpolation is usually unnecessary
For a self-contained network, plain names are enough:
```@example common_mistakes_interp
using Catalyst
rn = @reaction_network begin
    (p, d), 0 <--> X
end
species(rn), parameters(rn)
```

### When interpolation *is* needed
Use `$` when you want to reuse an already-created symbolic variable, parameter,
or algebraic expression (for example a shared Michaelis–Menten term). See
[Interpolating variables into the DSL](@ref dsl_advanced_options_symbolics_and_DSL_interpolation).

### Incomplete interpolation
A subtle bug is interpolating an *expression* that mentions outside symbols,
but forgetting to interpolate other species that should be the *same* objects
(or assuming outside declarations alone are enough — see the first section).
The expression's symbols are inserted as-is; names written without `$` are
still created by the DSL.

```@example common_mistakes_interp
t = default_t()
@species X1(t) E(t)
@parameters d
d_rate = d / (1 + E)

# Interpolate the shared rate expression; also interpolate X1 if you need the
# outside X1 object (e.g. same defaults / identity) to be the network species.
rn = @reaction_network begin
    $d_rate, $X1 --> 0
end
species(rn), parameters(rn)
```
!!! warning
    Writing `2$spec` is invalid Julia parsing in this context; use `2*$spec`.
    See the note in the
    [interpolation docs](@ref dsl_advanced_options_symbolics_and_DSL_interpolation).

## [Callbacks and `VariableRateJump` saving fill up memory](@id common_mistakes_save_positions)
Long stochastic simulations can quietly save enormous trajectories. Two common
causes are documented in
[this Discourse thread](https://discourse.julialang.org/t/memory-issues-for-long-ssa-simulations-with-mass-action-jumps/133530):

1. **Callbacks default to `save_positions = (true, true)`**, saving the state
   both before and after the affect. An always-true `DiscreteCallback` used to
   aggregate statistics therefore doubles (or more) the stored trajectory.
2. **`VariableRateJump`s** (reactions whose rates depend explicitly on time, or
   other non-mass-action forms that become variable-rate jumps) also save around
   jump events; pass `save_positions = (false, false)` to `JumpProblem` when you
   do not need the full path.

### Turn off callback saving
```@example common_mistakes_saving
using Catalyst, JumpProcesses, SciMLBase

two_state_model = @reaction_network begin
    (k1, k2), X1 <--> X2
end

p = (:k1 => 1.0, :k2 => 2.0)
u0 = [:X1 => 20, :X2 => 0]
tspan = (0.0, 50.0)

mutable struct TimeAggregator
    last_t::Float64
end
(a::TimeAggregator)(integrator) = (a.last_t = integrator.t)

# BAD for memory: default save_positions=(true, true) on an always-true callback
agg_default = TimeAggregator(0.0)
cb_default = DiscreteCallback((u, t, integrator) -> true, agg_default)
jp_default = JumpProblem(two_state_model, u0, tspan, p; save_positions = (false, false))
sol_default = solve(jp_default; callback = cb_default)
length(sol_default.t)
```
```@example common_mistakes_saving
# GOOD: disable callback saving when you only need side effects / aggregates
agg_off = TimeAggregator(0.0)
cb_off = DiscreteCallback((u, t, integrator) -> true, agg_off;
                          save_positions = (false, false))
jp_off = JumpProblem(two_state_model, u0, tspan, p; save_positions = (false, false))
sol_off = solve(jp_off; callback = cb_off)
length(sol_off.t)
```
With the same number of jumps, the second solution only keeps the start and end
points (length `2`), while the first stores thousands of duplicates.

### `VariableRateJump` models need an ODE timestepper
If any rate depends explicitly on `t` (or otherwise becomes a
`VariableRateJump`), `solve` requires an ODE algorithm, and
`save_positions` on `JumpProblem` still matters for how many jump times are
recorded:
```@example common_mistakes_saving
using OrdinaryDiffEq

circadian_model = @reaction_network begin
    A * (sin(2π * f * t - ϕ) + 1) / 2, 0 --> P
    d, P --> 0
end

u0map = [:P => 0]
pmap = [:f => 1.0, :A => 2.0, :ϕ => 0.0, :d => 1.0]
tspan = (0.0, 2.0)

jprob = JumpProblem(circadian_model, u0map, tspan, pmap;
                    save_positions = (false, false))
sol = solve(jprob, Tsit5())
length(sol.t)
```
See also [Jump simulations where some rate depends on time](@ref simulation_intro_jumps_variableratejumps)
and the DifferentialEquations.jl docs on
[Event Handling and Callback Functions](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/).
