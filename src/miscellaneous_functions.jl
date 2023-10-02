# Contains various functions that does not clearly belong to any other files. Hopefully these can eventually be moved to a better place.

### Steady State Stability ###

"""
    steady_state_stability(s_state, ps, sys)

Computes the stability of steady states. For `Vector{Float64}` input, returns a `Bool` corresponding the state's stability. For `Vector{Vector{Float64}}` input, returns a `Vector{Bool}` corresponding the states' stability. 

Arguments:
- `s_state`: The steate(s) for which we want to compute the stability. These are assumed to be system steady states.
- `ps` the parameter values for which we wish to compute stability.
- `sys`: The system for which we wish to comptue stability. Can be a `ReactionSystem` or an `ODEFunction`.

Example
```@repl
rs = @reaction_network begin
    (p,d), 0 <--> X
end
ps = [:p => 1.0, :d => 0.2]
steady_state_stability([5.0], ps, rs)
```
gives
```
true
```

Notes:
- If a `ReactionSystem` is given, it is coverted to a `ODEFunction` internally. If you plan to compute the stability for a large number of steady states. First covert your ReactionSystem to a ODEFunction using `ODEFunction(convert(ODESystem, rs); jac=true)` and then use this in your input. 


"""
steady_state_stability(s_states::Vector{Float64}, ps, rs::ReactionSystem) = steady_state_stability(s_states, ODEFunction(convert(ODESystem, rs); jac=true))
# If a ODESystem is given directyl, no conversion is needed.
function steady_state_stability(s_states::Vector{Float64}, ps, ofun::ODEFunction)
    maximum(eigen(ofun.jac(ss, last.(ps), Inf)).values) < 0.0
end
steady_state_stability(s_states::Vector{Vector{Float64}}, args...) =  steady_state_stability.(s_states, args...)
