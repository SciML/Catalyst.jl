### Stability Analysis ###

"""
    stability(u::Vector{T}, rs::ReactionSystem, p; sparse=false, ss_jac = steady_state_jac(u, rs, p; sparse=sparse), t = Inf, non_autonomous_warn = true)

Compute the stability of a steady state (Returned as a `Bool`, with `true` indicating stability).

Arguments:
    - `u`: The steady state for which we want to compute stability.
    - `rs`: The reaction system model for which we want to compute stability.
    - `p`: The parameter set for which we want to compute stability.
    - `sparse=false`: If we wish to create a sparse Jacobian for the stability computation.
    - `ss_jac = steady_state_jac(u, rs; sparse=sparse)`: It is possible to pre-compute the Jacobian used for stability computation using `steady_state_jac`. If stability is computed for many states, precomputing the Jacobian may speed up evaluation.
    - `t = Inf`: The time point at which stability is computed. 
    - `non_autonomous_warn = true`: If the system is non-autonomous (e.g. a rate depends on t), a warning will be given. Set this to false to remove that. Alternatively, specify a nonInf value for `t`.

Example:
```julia
# Creates model.
rn = @reaction_network begin
    (p,d), 0 <--> X
end
p = [:p => 1.0, :d => 0.5]

# Finds (the trivial) steady state, and computes stability.
steady_state = [2.0]
steady_state_stability(steady_state, rn, p)
```

Notes:
    y states (each being a `Vector`) is provided as `u`, stability is computed for each state separately.
"""
function steady_state_stability(u::Vector{T}, rs::ReactionSystem, p; 
                                sparse=false, ss_jac = steady_state_jac(rs; u0=u, sparse=sparse), 
                                t = Inf, non_autonomous_warn =true) where T
    # Warning checks.
    !is_autonomous(rs) && non_autonomous_warn && @warn "Attempting to compute stability for a non-autonomous system. Set `non_autonomous_warn = false` to disable this warning."

    # Because Jacobian currently requires u and p to be a normal vector.
    # Can be removed once this get fixed in MTK.
    if (u isa Vector{<:Pair}) || (u isa Dict) 
        u_dict = Dict(symmap_to_varmap(rs, u))
        u = [u_dict[var] for var in states(rs)]        
    end
    if (p isa Vector{<:Pair}) || (p isa Dict)
        p_dict = Dict(symmap_to_varmap(rs, p))
        p = [p_dict[var] for var in parameters(rs)]
    end

    # Computes stability (by checking that the real part of all eigenvalues are negative).
    jac = ss_jac(u, p, Inf)
    return maximum(real.(eigvals(jac))) < 0
end
# Computes the stability for a vector of steady states.
function steady_state_stability(us::Vector{Vector{T}}, rs::ReactionSystem, p; 
                                sparse=false, ss_jac = steady_state_jac(rs; u0=us[1], sparse=sparse)) where T
    return [steady_state_stability(u, rs, p) for u in us]
end

"""
    steady_state_jac(rs::ReactionSystem; u0=[], sparse=false)

Creates the Jacobian function which can be used as input to `steady_state_stability`. Useful when a large number of stability computation has to be carried out in a performant manner.
 
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
```

Arguments:
    - `rs`: The reaction system model for which we want to compute stability.
    - `u0=[]`: For systems with conservation laws, a `u` is required to compute the conserved quantities.
    - `sparse=false`: If we wish to create a sparse Jacobian for the stability computation.
"""
function steady_state_jac(rs::ReactionSystem; u0=[], sparse=false) where T
    # If u0 is a vector of values, must be converted to something MTK understands.
    if (u0 isa Vector{<:Number})
        if length(u0) == length(states(rs))
            u0 = Pair.(states(rs), u0)
        elseif !isempty(u0)
            error("You are trying to compute a stability matrix, providing u0 to compute conservation laws. Your provided u0 vector has length < the number of system states. If you provide a u0 vector, these have to be identical.")
        end
    end

    # Converts u0 to values MTK understands, and checks that potential conservation laws are accounted for.
    u0 = symmap_to_varmap(rs, u0)
    conservationlaw_errorcheck(rs, u0)

    # Creates a ODESystem and jacobian.
    osys = convert(ODESystem, rs; remove_conserved = true, defaults=Dict(u0))
    return ODEFunction(osys; jac=true, sparse=sparse).jac
end