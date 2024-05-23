### Stability Analysis ###

"""
    stability(u::Vector{T}, rs::ReactionSystem, p; tol = 10*sqrt(eps())
                ss_jac = steady_state_jac(u, rs, p))

Compute the stability of a steady state (Returned as a `Bool`, with `true` indicating stability).

Arguments:
- `u`: The steady state for which we want to compute stability.
- `rs`: The `ReactionSystem` model for which we want to compute stability.
- `p`: The parameter set for which we want to compute stability.
- `tol = 10*sqrt(eps())`: The tolerance used for determining whether computed eigenvalues real 
parts can reliably be considered negative/positive. In practise, a steady state is considered 
stable if the corresponding Jacobian's maximum eigenvalue real part is < 0. However, if this maximum 
eigenvalue is in the range `-tol< eig < tol`, and error is throw, as we do not deem that stability
can be ensured with enough certainty. The choice `tol = 10*sqrt(eps())` has *not* been subject
to much analysis.
- `ss_jac = steady_state_jac(u, rs)`: It is possible to pre-compute the 
Jacobian used for stability computation using `steady_state_jac`. If stability is computed 
for many states, precomputing the Jacobian may speed up evaluation.

Example:
```julia
# Creates model.
rn = @reaction_network begin
    (p,d), 0 <--> X
end
p = [:p => 1.0, :d => 0.5]

# Finds (the trivial) steady state, and computes its stability.
steady_state = [2.0]
steady_state_stability(steady_state, rn, p)

Notes:
- The input `u` can (currently) be both a vector of values (like `[1.0, 2.0]`) and a vector map 
(like `[X => 1.0, Y => 2.0]`). The reason is that most methods for computing stability only
produces vectors of values. However, if possible, it is recommended to work with values on the
form of maps.
- Catalyst currently computes steady state stabilities using the naive approach of checking whether
a system's largest eigenvalue real part is negative. While more advanced stability computation 
methods exist (and would be a welcome addition to Catalyst), there is no direct plans to implement 
these. Furthermore, Catalyst uses a arbitrary tolerance tol ~ 1.5*10^-7 to determine whether a 
computed eigenvalue is far away enough from 0 to be reliably used. This selected threshold, 
however, have not been subject to further analysis (and can be changed through the `tol` argument).
```
"""
function steady_state_stability(u::Vector, rs::ReactionSystem, ps; tol = 10*sqrt(eps(ss_val_type(u))),
                                ss_jac = steady_state_jac(rs; u0 = u))
    # Warning checks.
    if !is_autonomous(rs) 
        error("Attempting to compute stability for a non-autonomous system (e.g. where some rate depend on $(rs.iv)). This is not possible.")
    end

    # If `u` is a vector of values, we convert it to a map. Also, if there are conservation laws,
    # corresponding values must be eliminated from `u`.
    u = steady_state_u_conversion(u, rs)
    if length(u) > length(unknowns(ss_jac.f.sys))
        u = filter(u_v -> any(isequal(u_v[1]), unknowns(ss_jac.f.sys)), u)
    end

    # Generates the Jacobian at the steady state (technically, `ss_jac` is an `ODEProblem` with dummy values for `u0` and `p`).
    J = zeros(length(u), length(u))
    ss_jac = remake(ss_jac; u0 = u, p = ps)
    ss_jac.f.jac(J, ss_jac.u0, ss_jac.p, Inf)
    
    # Computes stability (by checking that the real part of all eigenvalues is negative).
    max_eig = maximum(real(ev) for ev in eigvals(J))
    if abs(max_eig) < tol
        error("The system Jacobian's maximum eigenvalue at the steady state is within the tolerance range (abs($max_eig) < $tol). Hence, stability could not be reliably determined. If you still wish to compute the stability, reduce the `tol` argument, but note that floating point error in the eigenvalue calculation could lead to incorrect results.")
    end
    return max_eig < 0
end

# Used to determine the type of the steady states values, which is then used to set the tolerance's
# type.
ss_val_type(u::Vector{T}) where {T} = T
ss_val_type(u::Vector{Pair{S,T}}) where {S,T} = T
ss_val_type(u::Dict{S,T}) where {S,T} = T

"""
    steady_state_jac(rs::ReactionSystem; u0 = [])

Creates the Jacobian function which can be used as input to `steady_state_stability`. Useful when 
a large number of stability computation has to be carried out in a performant manner.

Arguments:
    - `rs`: The reaction system model for which we want to compute stability.
    - `u0 = []`: For systems with conservation laws, a `u` is required to compute the conserved quantities.
 
Example:
```julia
# Creates model.
rn = @reaction_network begin
    (p,d), 0 <--> X
end
p = [:p => 1.0, :d => 0.5]

# Creates the steady state jacobian.
ss_jac = steady_state_jacobian(rn)

# Finds (the trivial) steady state, and computes stability.
steady_state = [2.0]
steady_state_stability(steady_state, rn, p; ss_jac)

Notes:
    - In practise, this function returns an `ODEProblem` (with a computed Jacobian) set up in 
    such a way that it can be used by the `steady_state_stability` function.
```
"""
function steady_state_jac(rs::ReactionSystem; u0 = [sp => 0.0 for sp in unknowns(rs)], 
                          combinatoric_ratelaws = get_combinatoric_ratelaws(rs))
    # If u0 is a vector of values, must be converted to something MTK understands.

    # Converts u0 to values MTK understands, and checks that potential conservation laws are accounted for.
    u0 = steady_state_u_conversion(u0, rs)
    conservationlaw_errorcheck(rs, u0)

    # Creates an `ODEProblem` with a Jacobian. Dummy values for `u0` and `ps` must be provided.
    ps = [p => 0.0 for p in parameters(rs)]
    return ODEProblem(rs, u0, 0, ps; jac = true, remove_conserved = true, 
                      combinatoric_ratelaws = combinatoric_ratelaws)
end

# Converts a `u` vector from a vector of values to a map.
function steady_state_u_conversion(u, rs::ReactionSystem)
    if (u isa Vector{<:Number})
        if length(u) == length(unknowns(rs))
            u = [sp => v for (sp,v) in zip(unknowns(rs), u)]
        else
            error("You are trying to generate a stability Jacobian, providing u0 to compute conservation laws. Your provided u0 vector has length < the number of system states. If you provide a u0 vector, these have to be identical.")
        end
    end
    return symmap_to_varmap(rs, u)
end
