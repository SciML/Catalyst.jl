# [Hybrid simulations](@id hybrid_simulations)
```@raw html
<details><summary><strong>Environment setup and package installation</strong></summary>
```
The following code sets up an environment for running the code on this page.
```julia
using Pkg
Pkg.activate(; temp = true) # Creates a temporary environment, which is deleted when the Julia session ends.
Pkg.add("Catalyst")
Pkg.add("JumpProcesses")
Pkg.add("OrdinaryDiffEqDefault")
Pkg.add("Plots")
Pkg.add("StochasticDiffEq")
```
```@raw html
</details>
```
```@raw html
<details><summary><strong>Quick-start example</strong></summary>
```
The following code provides a brief example of how to run a hybrid simulation (here we run a hybrid ODE/Jump simulation).
```julia
using Catalyst, JumpProcesses, OrdinarDiffEqDefault, Plots

# First we create the normal problem we wish to run an ensemble simulation of.
rn = @reaction_network begin
    (kOn, kOff), Gi <--> Ga
    p*Ga, 0 --> X
    d, X --> 0
end
u0 = [:X => 100.0]
ps = [:p => 50.0, :d => 0.25]
hprob = HybridProblem(rn, u0, 10.0, ps)
sol = solve(hprob)
plot(esol)
```
```@raw html
</details>
```
  \

