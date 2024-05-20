### Stability Analysis ###

"""
    stability(u::Vector{T}, rs::ReactionSystem, p; sparse = false, 
              ss_jac = steady_state_jac(u, rs, p; sparse=sparse))

Compute the stability of a steady state (Returned as a `Bool`, with `true` indicating stability).

Arguments:
    - `u`: The steady state for which we want to compute stability.
    - `rs`: The `ReactionSystem` model for which we want to compute stability.
    - `p`: The parameter set for which we want to compute stability.
    - `sparse = false`: If we wish to create a sparse Jacobian for the stability computation.
    - `ss_jac = steady_state_jac(u, rs; sparse = sparse)`: It is possible to pre-compute the 
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
```
"""
function steady_state_stability(u::Vector, rs::ReactionSystem, ps; sparse = false, 
                                ss_jac = steady_state_jac(rs; u0 = u, sparse = sparse))
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

    # Computes stability (by checking that the real part of all eigenvalues are negative).
    # Here, `ss_jac` is a `ODEProblem` with dummy values for `u0` and `p`.
    J = zeros(length(u), length(u))
    ss_jac = remake(ss_jac; u0 = u, p = ps)
    ss_jac.f.jac(J, ss_jac.u0, ss_jac.p, Inf)
    return maximum(real(ev) for ev in (eigvals(J))) < 0
end

"""
    steady_state_jac(rs::ReactionSystem; u0 = [], sparse = false)

Creates the Jacobian function which can be used as input to `steady_state_stability`. Useful when 
a large number of stability computation has to be carried out in a performant manner.

Arguments:
    - `rs`: The reaction system model for which we want to compute stability.
    - `u0 = []`: For systems with conservation laws, a `u` is required to compute the conserved quantities.
    - `sparse = false`: If we wish to create a sparse Jacobian for the stability computation.
 
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
                          sparse = false, combinatoric_ratelaws = get_combinatoric_ratelaws(rs))
    # If u0 is a vector of values, must be converted to something MTK understands.

    # Converts u0 to values MTK understands, and checks that potential conservation laws are accounted for.
    u0 = steady_state_u_conversion(u0, rs)
    conservationlaw_errorcheck(rs, u0)

    # Creates an `ODEProblem` with a Jacobian. Dummy values for `u0` and `ps` must be provided.
    ps = [p => 0.0 for p in parameters(rs)]
    return ODEProblem(rs, u0, 0, ps; jac = true, remove_conserved = true, 
                      combinatoric_ratelaws = get_combinatoric_ratelaws(rs))
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