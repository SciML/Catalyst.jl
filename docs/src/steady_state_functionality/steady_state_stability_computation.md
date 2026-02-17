# [Steady state stability computation](@id steady_state_stability)
```@raw html
<details><summary><strong>Environment setup and package installation</strong></summary>
```
The following code sets up an environment for running the code on this page.
```julia
using Pkg
Pkg.activate(; temp = true) # Creates a temporary environment, which is deleted when the Julia session ends.
Pkg.add("Catalyst")
Pkg.add("HomotopyContinuation")
```
```@raw html
</details>
```
```@raw html
<details><summary><strong>Quick-start example</strong></summary>
```
The following code provides a brief example of how a steady state's stability can be determine using the `steady_state_stability` function.
```julia
using Catalyst
rn = @reaction_network begin 
    (p,d), 0 <--> X
end
ps = [:p => 2.0, :d => 0.5]
steady_state = [:X => 4.0] # This must be a true steady state of the system.
steady_state_stability(steady_state, rn, ps)
```
```@raw html
</details>
```
  \
  

After system steady states have been found using [HomotopyContinuation.jl](@ref homotopy_continuation), [NonlinearSolve.jl](@ref steady_state_solving), or other means, their stability can be computed using Catalyst's `steady_state_stability` function. Systems with conservation laws will automatically have these removed, permitting stability computation on systems with singular Jacobian.

!!! warning 
    Catalyst currently computes steady state stabilities using the naive approach of checking whether a system's largest eigenvalue real part is negative. Furthermore, Catalyst uses a tolerance `tol = 10*sqrt(eps())` to determine whether a computed eigenvalue is far away enough from 0 to be reliably considered non-zero. This threshold can be changed through the `tol` keyword argument.

## [Basic examples](@id steady_state_stability_basics)
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

Next, let us consider the following [self-activation loop](@ref basic_CRN_library_self_activation):
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
Next, we can apply `steady_state_stability` to each steady state yielding a vector of their stabilities:
```@example stability_1
[steady_state_stability(sstate, sa_loop, ps) for sstate in steady_states]
```

Finally, as described above, Catalyst uses an optional argument, `tol`, to determine how strict to make the stability check.  I.e. below we set the tolerance to `1e-6` (a larger value, that is stricter, than the default of `10*sqrt(eps())`)
```@example stability_1
[steady_state_stability(sstate, sa_loop, ps; tol = 1e-6) for sstate in steady_states]
nothing# hide
```

## [Pre-computing the Jacobian to increase performance when computing stability for many steady states](@id steady_state_stability_jacobian)
Catalyst uses the system Jacobian to compute steady state stability, and the Jacobian is computed once for each call to `steady_state_stability`. If you repeatedly compute stability for steady states of the same system, pre-computing the Jacobian and supplying it to the `steady_state_stability` function can improve performance. 

In this example we use the self-activation loop from previously, pre-computes its Jacobian, and uses it to multiple `steady_state_stability` calls:
```@example stability_1
ss_jac = steady_state_jac(sa_loop)

ps_1 = [:v => 2.0, :K => 0.5, :n => 3, :d => 1.0]
steady_states_1 = hc_steady_states(sa_loop, ps)
stabs_1 = [steady_state_stability(st, sa_loop, ps_1; ss_jac) for st in steady_states_1]

ps_2 = [:v => 4.0, :K => 1.5, :n => 2, :d => 1.0]
steady_states_2 = hc_steady_states(sa_loop, ps)
stabs_2 = [steady_state_stability(st, sa_loop, ps_2; ss_jac) for st in steady_states_2]
nothing # hide
```

!!! warning
    For systems with [conservation laws](@ref homotopy_continuation_conservation_laws), `steady_state_jac` must be supplied a `u0` vector (indicating species concentrations for conservation law computation). This is required to eliminate the conserved quantities, preventing a singular Jacobian. These are supplied using the `u0` optional argument.
