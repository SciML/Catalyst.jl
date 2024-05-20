# Steady state stability computation
After system steady states have been found using [HomotopyContinuation.jl](@ref homotopy_continuation), [NonlinearSolve.jl](@ref nonlinear_solve), or other means, their stability can be computed using Catalyst's `steady_state_stability` function. Systems with conservation laws will automatically have these removed, permitting stability computation on systems with singular Jacobian.

!!! warn 
    Catalyst currently computes steady state stabilities using the naive approach of checking whether a system's largest eigenvalue is negative. While more advanced stability computation methods exist (and would be a welcome addition to Catalyst), there is no direct plans to implement these.

## Basic examples
Let us consider the following basic example:
```@example stability_1
using Catalyst
rn = @reaction_network begin 
    (p,d), 0 <--> X
end
```
It has a single (stable) steady state at $X = p/d$. We can confirm stability using the `steady_state_stability` function, to which we provide the steady state, the reaction system, and the parameter values:
```@example stability_1
ps = [:p => 2.0, :d => 0.5]
steady_state = [:X => 4.0]
steady_state_stability(steady_state, rn, ps)
```

Next, let us consider the following self-activation loop:
```@example stability_1
sa_loop = @reaction_network begin 
    (hill(X,v,K,n),d), 0 <--> X
end
```
For certain parameter choices, this system exhibits multi-stability. Here, we can find the steady states [using homotopy continuation](@ref homotopy_continuation):
```@example stability_1
import HomotopyContinuation
ps = [:v => 2.0, :K => 0.5, :n => 3, :d => 1.0]
steady_states = hc_steady_states(sa_loop, ps)
```
Next, we can apply `steady_state_stability` directly to this steady state vector, receiving the stability for each:
```@example stability_1
steady_state_stability(steady_states, sa_loop, ps)
```

!!! note
    For systems with [conservation laws](@ref homotopy_continuation_conservation_laws), `steady_state_jac` must be supplied a `u0` vector (indicating species concentrations for conservation law computation). This is required to eliminate the conserved quantities, preventing a singular Jacobian. These are supplied using the `u0` optional argument.

## Pre-computing the Jacobian to increase performance when computing stability for many steady states
Catalyst uses the system Jacobian to compute steady state stability, and the Jacobian is computed once for each call to `steady_state_stability`. If you repeatedly compute stability for steady states of the same system, pre-computing the Jacobian and supplying it to the `steady_state_stability` function can improve performance. 

In this example we use the self-activation loop from previously, pre-computes the Jacobian, and uses it to multiple `steady_state_stability` calls:
```@example stability_1
ss_jac = steady_state_jac(sa_loop)

ps_1 = [:v => 2.0, :K => 0.5, :n => 3, :d => 1.0]
steady_states_1 = hc_steady_states(sa_loop, ps)
stability_1 = steady_state_stability(steady_states_1, sa_loop, ps_1; ss_jac=ss_jac)

ps_2 = [:v => 4.0, :K => 1.5, :n => 2, :d => 1.0]
steady_states_2 = hc_steady_states(sa_loop, ps)
stability_2 = steady_state_stability(steady_states_2, sa_loop, ps_2; ss_jac=ss_jac)
nothing # hide
```

It is possible to designate that a [sparse Jacobian](@ref ref) should be used using the `sparse = true` option (either to `steady_state_jac` or directly to `steady_state_stability`):
```@example stability_1
ss_jac = steady_state_jac(sa_loop; sparse = true)
nothing # hide
```