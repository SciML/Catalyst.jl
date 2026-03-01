### Lattice Simulation Structure Species Getters/Setters ###

"""
    spat_setp!(sim_struct, p, dsrs::DiscreteSpaceReactionSystem, p_val)

For a problem or integrators, update its `p` vector with the input `p_val`. For non-discrete space models,
this is can be done through direct interfacing (e.g. `prob[p] = 1.0`). However, for
`DiscreteSpaceReactionSystem`-based problems and integrators, this function must be used instead.

Arguments:
- `sim_struct`: The simulation structure which `u` value we wish to update. Can be either a `ODEProblem`, `JumpProblem`, or an integrator derived from either of these.
- `p`: The species which value we wish to update. Can be provided either in its symbolic form (e.g. `k`) or as a symbol (e.g. `:k`).
- `dsrs`: The `DiscreteSpaceReactionSystem` which was used to generate the structure we wish to modify.
- `p_val`: The parameter's new values. Must be given in a form which is also a valid initial input to the `ODEProblem`/`JumpProblem`.

Example:
```julia
# Prepare `DiscreteSpaceReactionSystem`s.
using Catalyst
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
dsrs = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,3)))

# Prepares a corresponding ODEProblem.
u0 = [:X1 => 1.0, :X2 => 2.0]
tspan = (0.0, 50.0)
ps = [:k1 => [1.0 2.0 3.0; 4.0 5.0 6.0], :k2 => 1.0, :D => 0.01]
oprob = ODEProblem(dsrs, u0, tspan, ps)

# Updates the `ODEProblem`.
spat_setp!(oprob, :k1, dsrs, 0.0) # Sets `k1` to uniformly 0 across the discrete space.
spat_setp!(oprob, :k2, dsrs, [1.0 0.0 0.0; 0.0 0.0 0.0]) # Sets `k2` to `1.0` in one vertex, and 0 elsewhere.
```
"""
function spat_setp!(sim_struct, p, dsrs::DiscreteSpaceReactionSystem, p_vals)
    # Error checks.
    (p_vals isa Number) || check_dspace_format(extract_dspace(dsrs), p_vals)
    edge_param_check(p, dsrs)

    # Converts symbol parameter to symbolic and find the correct species index and numbers.
    (p isa Symbol) && (p = _symbol_to_var(dsrs, p))
    (p isa Num) && (p = unwrap(p))
    p_idx, p_tot = get_p_idxs(p, dsrs)

    # Reshapes the values to a vector of the correct form, and calls spat_setp! on the input structure.
    p_vals_reshaped = vertex_value_form(p_vals, dsrs, p)
    spat_setp!(sim_struct, p_idx, p_vals_reshaped, num_verts(dsrs))
end

# Note: currently, `spat_setp!(oprob::ODEProblem, ...`) and `spat_setp!(SciMLBase.AbstractODEIntegrator, ...`)
# are identical and could be merged to a single function.
function spat_setp!(oprob::ODEProblem, p_idx::Int64, p_vals, num_verts)
    if length(p_vals) == 1
        foreach(idx -> (oprob.p[p_idx][idx] = p_vals[1]), 1:num_verts)
    elseif length(p_vals) == length(oprob.p[p_idx])
        foreach(idx -> (oprob.p[p_idx][idx] = p_vals[idx]), 1:num_verts)
    elseif length(oprob.p[p_idx]) == 1
        oprob.p[p_idx][1] = p_vals[1]
        foreach(idx -> (push!(oprob.p[p_idx], p_vals[idx])), 2:num_verts)
    end
end
function spat_setp!(jprob::JumpProblem, p_idx::Int64, p_vals, num_verts)
    error("The `spat_setp!` function is currently not supported for `JumpProblem`s.")
end
function spat_setp!(oint::SciMLBase.AbstractODEIntegrator, p_idx::Int64, p_vals, num_verts)
    if length(p_vals) == 1
        foreach(idx -> (oint.p[p_idx][idx] = p_vals[1]), 1:num_verts)
    elseif length(p_vals) == length(oint.p[p_idx])
        foreach(idx -> (oint.p[p_idx][idx] = p_vals[idx]), 1:num_verts)
    elseif length(oint.p[p_idx]) == 1
        oint.p[p_idx][1] = p_vals[1]
        foreach(idx -> (push!(oint.p[p_idx], p_vals[idx])), 2:num_verts)
    end
end
function spat_setp!(jint::JumpProcesses.SSAIntegrator, p_idx::Int64, p_vals, num_verts)
    error("The `spat_setp!` function is currently not supported for jump simulation integrators.")
end

"""
    spat_getp(sim_struct, p, dsrs::DiscreteSpaceReactionSystem)

For a problem or integrators, retrieves its `p` values. For non-discrete space models,
this is can be done through direct interfacing (e.g. `prob[p]`). However, for
`DiscreteSpaceReactionSystem`-based problems and integrators, this function must be used instead. The
output format depends on the discrete space (a dense array for cartesian grid discrete spaces, a sparse array for
masked grid discrete spaces, and a vector for graph discrete spaces). This format is similar to what is used to
designate parameter initial values.

Arguments:
- `sim_struct`: The simulation structure which `p` value we wish to retrieve. Can be either a `ODEProblem`,
`JumpProblem`, or an integrator derived from either of these.
- `p`: The species which value we wish to update. Can be provided either in its symbolic form (e.g. `k`) or as a symbol (e.g. `:k`).
- `dsrs`: The `DiscreteSpaceReactionSystem` which was used to generate the structure we wish to modify.

Notes:
- Even if the parameter is spatially uniform, a full array with its values across all vertices will be retrieved.

Example:
```julia
# Prepare `DiscreteSpaceReactionSystem`s.
using Catalyst
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
dsrs = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,3)))

# Prepares a corresponding ODEProblem.
u0 = [:X1 => 1.0, :X2 => 2.0]
tspan = (0.0, 50.0)
ps = [:k1 => [1.0 2.0 3.0; 4.0 5.0 6.0], :k2 => 1.0, :D => 0.01]
oprob = ODEProblem(dsrs, u0, tspan, ps)

# Updates the `ODEProblem`.
spat_getp(oprob, :k1, dsrs) # Retrieves the value of `k1`.
```
"""
function spat_getp(sim_struct, p, dsrs::DiscreteSpaceReactionSystem)
    edge_param_check(p, dsrs)
    p_idx, _ = get_p_idxs(p, dsrs)
    spat_getp(sim_struct, p_idx, extract_dspace(dsrs), num_verts(dsrs))
end

# Retrieves the discrete space values for problem or integrator structures.
function spat_getp(oprob::ODEProblem, p_idx::Int64, dspace, num_verts)
    vals = oprob.p[p_idx]
    (length(vals) == 1) && (vals = fill(vals[1], num_verts))
    return reshape_vals(vals, dspace)
end
function spat_getp(jprob::JumpProblem, p_idx::Int64, dspace, num_verts)
    error("The `spat_getp` function is currently not supported for `JumpProblem`s.")
end
function spat_getp(oint::SciMLBase.AbstractODEIntegrator, p_idx::Int64, dspace, num_verts)
    vals = oint.p[p_idx]
    (length(vals) == 1) && (vals = fill(vals[1], num_verts))
    return reshape_vals(vals, dspace)
end
function spat_getp(jint::JumpProcesses.SSAIntegrator, p_idx::Int64, dspace, num_verts)
    error("The `spat_getp` function is currently not supported for jump simulation integrators.")
end

"""
    spat_setu!(sim_struct, sp, dsrs::DiscreteSpaceReactionSystem, u)

For a problem or integrators, update its `u` vector with the input `u`. For non-discrete space models,
this is can be done through direct interfacing (e.g. `prob[X] = 1.0`). However, for
`DiscreteSpaceReactionSystem`-based problems and integrators, this function must be used instead.

Arguments:
- `sim_struct`: The simulation structure which `u` value we wish to update. Can be either a `ODEProblem`, `JumpProblem`, or an integrator derived from either of these.
- `sp`: The species which value we wish to update. Can be provided either in its symbolic form (e.g. `X`) or as a symbol (e.g. `:X`).
- `dsrs`: The `DiscreteSpaceReactionSystem` which was used to generate the structure we wish to modify.
- `u`: The species's new values. Must be given in a form which is also a valid initial input to the `ODEProblem`/`JumpProblem`.

Example:
```julia
# Prepare `DiscreteSpaceReactionSystem`s.
using Catalyst
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
dsrs = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,3)))

# Prepares a corresponding ODEProblem.
u0 = [:X1 => [1.0 2.0 3.0; 4.0 5.0 6.0], :X2 => 2.0]
tspan = (0.0, 50.0)
ps = [:k1 => 2.0, :k2 => 1.0, :D => 0.01]
oprob = ODEProblem(dsrs, u0, tspan, ps)

# Updates the `ODEProblem`.
spat_setu!(oprob, :X1, dsrs, 0.0) # Sets `X1` to uniformly 0 across the discrete space.
spat_setu!(oprob, :X2, dsrs, [1.0 0.0 0.0; 0.0 0.0 0.0]) # Sets `X2` to `1.0` in one vertex, and 0 elsewhere.
```
"""
function spat_setu!(sim_struct, sp, dsrs::DiscreteSpaceReactionSystem, u)
    # Checks that if u is non-uniform, it has the correct format for the system's discrete space.
    (u isa Number) || check_dspace_format(extract_dspace(dsrs), u)

    # Converts symbol species to symbolic and finds correct species index and numbers.
    (sp isa Symbol) && (sp = _symbol_to_var(dsrs, sp))
    (sp isa Num) && (sp = unwrap(sp))
    sp_idx, sp_tot = get_sp_idxs(sp, dsrs)

    # Reshapes the values to a vector of the correct form, and calls spat_setu! on the input structure.
    u_reshaped = vertex_value_form(u, dsrs, sp)
    spat_setu!(sim_struct, sp_idx, sp_tot, u_reshaped, num_verts(dsrs))
end

function spat_setu!(oprob::ODEProblem, sp_idx::Int64, sp_tot::Int64, u, num_verts)
    if length(u) == 1
        foreach(idx -> (oprob.u0[sp_idx + (idx - 1) * sp_tot] = u[1]), 1:num_verts)
    else
        foreach(idx -> (oprob.u0[sp_idx + (idx - 1) * sp_tot] = u[idx]), 1:num_verts)
    end
end
function spat_setu!(jprob::JumpProblem, sp_idx::Int64, sp_tot::Int64, u, num_verts)
    if length(u) == 1
        foreach(idx -> (jprob.prob.u0[sp_idx, idx] = u[1]), 1:num_verts)
    else
        foreach(idx -> (jprob.prob.u0[sp_idx, idx] = u[idx]), 1:num_verts)
    end
end
function spat_setu!(oint::SciMLBase.AbstractODEIntegrator, sp_idx::Int64, sp_tot::Int64,
        u, num_verts)
    if length(u) == 1
        foreach(idx -> (oint.u[sp_idx + (idx - 1) * sp_tot] = u[1]), 1:num_verts)
    else
        foreach(idx -> (oint.u[sp_idx + (idx - 1) * sp_tot] = u[idx]), 1:num_verts)
    end
end
function spat_setu!(jint::JumpProcesses.SSAIntegrator, sp_idx::Int64, sp_tot::Int64,
        u, num_verts)
    if length(u) == 1
        foreach(idx -> (jint.u[sp_idx, idx] = u[1]), 1:num_verts)
    else
        foreach(idx -> (jint.u[sp_idx, idx] = u[idx]), 1:num_verts)
    end
end

"""
    spat_getu(sim_struct, sp, dsrs::DiscreteSpaceReactionSystem)

For a problem or integrators, retrieves its `u` values. For non-discrete space models,
this is can be done through direct interfacing (e.g. `prob[X]`). However, for
`DiscreteSpaceReactionSystem`-based problems and integrators, this function must be used instead. The
output format depends on the discrete space (a dense array for cartesian grid discrete spaces, a sparse array for
masked grid discrete spaces, and a vector for graph discrete spaces). This format is similar to which is used to
designate species initial conditions.

Arguments:
- `sim_struct`: The simulation structure which `u` value we wish to retrieve. Can be either a `ODEProblem`, `JumpProblem`, or an integrator derived from either of these.
- `sp`: The species which value we wish to update. Can be provided either in its symbolic form (e.g. `X`) or as a symbol (e.g. `:X`).
- `dsrs`: The `DiscreteSpaceReactionSystem` which was used to generate the structure we wish to modify.

Notes:
- Even if the species is spatially uniform, a full array with its values across all vertices will be retrieved.

Example:
```julia
# Prepare `DiscreteSpaceReactionSystem`s.
using Catalyst
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
dsrs = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,3)))

# Prepares a corresponding ODEProblem.
u0 = [:X1 => [1.0 2.0 3.0; 4.0 5.0 6.0], :X2 => 2.0]
tspan = (0.0, 50.0)
ps = [:k1 => 2.0, :k2 => 1.0, :D => 0.01]
oprob = ODEProblem(dsrs, u0, tspan, ps)

# Updates the `ODEProblem`.
spat_getu(oprob, :X1, dsrs) # Retrieves the value of `X1`.
```
"""
function spat_getu(sim_struct, sp, dsrs::DiscreteSpaceReactionSystem)
    sp_idx, sp_tot = get_sp_idxs(sp, dsrs)
    spat_getu(sim_struct, sp_idx, sp_tot, extract_dspace(dsrs))
end

# Retrieves the discrete space values for problem or integrator structures.
function spat_getu(oprob::ODEProblem, sp_idx, sp_tot, dspace)
    return reshape_vals(oprob.u0[sp_idx:sp_tot:end], dspace)
end
function spat_getu(jprob::JumpProblem, sp_idx, sp_tot, dspace)
    return reshape_vals(jprob.prob.u0[sp_idx, :], dspace)
end
function spat_getu(oint::SciMLBase.AbstractODEIntegrator, sp_idx, sp_tot, dspace)
    return reshape_vals(oint.u[sp_idx:sp_tot:end], dspace)
end
function spat_getu(jint::JumpProcesses.SSAIntegrator, sp_idx, sp_tot, dspace)
    return reshape_vals(jint.u[sp_idx, :], dspace)
end

# A single function, `spat_getu`, which contains all interfacing functionality. However,
# long-term it should be replaced with a sleeker interface. Ideally as MTK-wide support for
# discrete space problems and solutions is introduced. Note that SciML considers jump simulation solutions
# as `ODESolution`, hence that type is specified.
"""
    spat_getu(sol, sp, dsrs::DiscreteSpaceReactionSystem; t = nothing)

A function for retrieving the solution of a `DiscreteSpaceReactionSystem`-based simulation on various
desired forms. Generally, for `DiscreteSpaceReactionSystem`s, the values in `sol` is ordered in a
way which is not directly interpretable by the user. Furthermore, the normal Catalyst interface
for solutions (e.g. `sol[:X]`) does not work for these solutions. Hence this function is used instead.

The output is a vector, which in each position contains sp's value (either at a time step of time,
depending on the input `t`). Its shape depends on the discrete space (using a similar form as heterogeneous
initial conditions). I.e. for a NxM cartesian grid, the values are NxM matrices. For a masked grid,
the values are sparse matrices. For a graph discrete space, the values are vectors (where the value in
the n'th position corresponds to sp's value in the n'th vertex).

Arguments:
- `sol`: The solution from which we wish to retrieve some values.
- `sp`: The species which value we wish to update. Can be provided either in its symbolic form (e.g. `X`) or as a symbol (e.g. `:X`).
- `dsrs`: The `DiscreteSpaceReactionSystem` which was simulated to generate the solution.
- `t = nothing`: If `nothing`, we simply return the solution across all saved time steps (default). If `t` instead is a vector (or range of values), returns the solution interpolated at these time points.

Notes:
- The `spat_getu` is not optimised for performance. However, it should still be quite performant, but there might be some limitations if called a very large number of times.
- Long-term it is likely that this function gets replaced with a sleeker interface.

Example:
```julia
using Catalyst, OrdinaryDiffEqDefault

# Prepare `DiscreteSpaceReactionSystem`s.
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
dsrs = DiscreteSpaceReactionSystem(rs, [tr], CartesianGrid((2,2)))

# Create problems.
u0 = [:X1 => 1, :X2 => 2]
tspan = (0.0, 10.0)
ps = [:k1 => 1, :k2 => 2.0, :D => 0.1]

oprob = ODEProblem(dsrs1, u0, tspan, ps)
osol = solve(oprob)
spat_getu(osol, :X1, dsrs) # Returns the value of X1 at each time step.
spat_getu(osol, :X1, dsrs; t = 0.0:10.0) # Returns the value of X1 at times 0.0, 1.0, ..., 10.0
```
"""
function spat_getu(sol::ODESolution, sp, dsrs::DiscreteSpaceReactionSystem; t = nothing)
    (t isa Number) && error("The input `t` to `spat_getu` must be a `AbstractVector`.")
    sp_idx, sp_tot = get_sp_idxs(sp, dsrs)
    spat_getu(sol, extract_dspace(dsrs), t, sp_idx, sp_tot)
end

# Function which handles the input in the case where `t` is `nothing` (i.e. return `sp`s value
# across all sample points).
function spat_getu(sol::ODESolution, dspace, t::Nothing, sp_idx::Int64, sp_tot::Int64)
    # ODE simulations contain, in each data point, all values in a single vector. Jump simulations
    # instead in a matrix (NxM, where N is the number of species and M is the number of vertices). We
    # must consider each case separately.
    if sol.prob isa ODEProblem
        return [reshape_vals(vals[sp_idx:sp_tot:end], dspace) for vals in sol.u]
    elseif sol.prob isa DiscreteProblem
        return [reshape_vals(vals[sp_idx, :], dspace) for vals in sol.u]
    else
        error("Unknown type of solution provided to `spat_getu`. Only ODE or Jump solutions are supported.")
    end
end

# Function which handles the input in the case where `t` is a range of values (i.e. return `sp`s
# value at all designated time points.
function spat_getu(sol::ODESolution, dspace, t::AbstractVector{T}, sp_idx::Int64,
        sp_tot::Int64) where {T <: Number}
    # Checks that an appropriate `t` is provided (however, DiffEq does permit out-of-range `t`s).
    if (minimum(t) < sol.t[1]) || (maximum(t) > sol.t[end])
        error("The range of the t values provided for sampling, ($(minimum(t)),$(maximum(t))) is not fully within the range of the simulation time span ($(sol.t[1]),$(sol.t[end])).")
    end

    # ODE simulations contain, in each data point, all values in a single vector. Jump simulations
    # instead in a matrix (NxM, where N is the number of species and M is the number of vertices). We
    # must consider each case separately.
    if sol.prob isa ODEProblem
        return [reshape_vals(sol(ti)[sp_idx:sp_tot:end], dspace) for ti in t]
    elseif sol.prob isa DiscreteProblem
        return [reshape_vals(sol(ti)[sp_idx, :], dspace) for ti in t]
    else
        error("Unknown type of solution provided to `spat_getu`. Only ODE or Jump solutions are supported.")
    end
end

### Lattice Internal Rebuilding ###

"""
    rebuild_spat_internals!(sciml_struct)

Rebuilds the internal functions for simulating a DiscreteSpaceReactionSystem. Whenever a problem or
integrator has had its parameter values updated, this function should be called for the update to
be taken into account. For ODE simulations, `rebuild_spat_internals!` needs only to be called when
- An edge parameter has been updated.
- When a parameter with spatially homogeneous values has been given spatially heterogeneous values (or vice versa).

Arguments:
- `sciml_struct`: The problem (e.g. an `ODEProblem`) or an integrator which we wish to rebuild.

Notes:
- Currently does not work for `DiscreteProblem`s, `JumpProblem`s, or their integrators.
- The function is not built with performance in mind, so avoid calling it multiple times in performance-critical applications.

Example:
```julia
# Creates an initial `ODEProblem`
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
grid = CartesianGrid((2,2))
dsrs = DiscreteSpaceReactionSystem(rs, [tr], grid)

u0 = [:X1 => 2, :X2 => [5 6; 7 8]]
tspan = (0.0, 10.0)
ps = [:k1 => 1.5, :k2 => [1.0 1.5; 2.0 3.5], :D => 0.1]

oprob = ODEProblem(dsrs, u0, tspan, ps)

# Updates parameter values.
oprob.ps[:ks] = [2.0 2.5; 3.0 4.5]
oprob.ps[:D] = 0.05

# Rebuilds `ODEProblem` to make changes have an effect.
rebuild_spat_internals!(oprob)
```
"""
function rebuild_spat_internals!(oprob::ODEProblem)
    rebuild_spat_internals!(oprob.f.f, oprob.p, oprob.f.f.dsrs)
end

# Function for rebuilding a `DiscreteSpaceReactionSystem` integrator after it has been updated.
function rebuild_spat_internals!(integrator::SciMLBase.AbstractODEIntegrator)
    rebuild_spat_internals!(integrator.f.f, integrator.p, integrator.f.f.dsrs)
end

# Function which rebuilds a `LatticeTransportODEFunction` functor for a new parameter set.
function rebuild_spat_internals!(lt_ofun::LatticeTransportODEFunction, ps_new,
        dsrs::DiscreteSpaceReactionSystem)
    # Computes Jacobian properties.
    jac = !isnothing(lt_ofun.jac_transport)
    sparse = lt_ofun.sparse

    # Recreates the new parameters on the requisite form.
    ps_new = [(length(p) == 1) ? p[1] : p for p in deepcopy(ps_new)]
    ps_new = [p => p_val for (p, p_val) in zip(parameters(dsrs), deepcopy(ps_new))]
    vert_ps,
    edge_ps = dspace_process_p(ps_new, vertex_parameters(dsrs),
        edge_parameters(dsrs), dsrs)
    ps_new = [vert_ps; edge_ps]

    # Creates the new transport rates and transport Jacobian part.
    transport_rates = make_sidxs_to_transrate_map(vert_ps, edge_ps, dsrs)
    if !isnothing(lt_ofun.jac_transport)
        lt_ofun.jac_transport .= 0.0
        set_jac_transport_values!(lt_ofun.jac_transport, transport_rates, dsrs)
    end

    # Computes new field values.
    heterogeneous_vert_p_idxs = make_heterogeneous_vert_p_idxs(ps_new, dsrs)
    mtk_ps, p_setters = make_mtk_ps_structs(ps_new, dsrs, heterogeneous_vert_p_idxs)
    t_rate_idx_types, leaving_rates = make_t_types_and_leaving_rates(transport_rates, dsrs)

    # Updates functor fields.
    replace_vec!(lt_ofun.heterogeneous_vert_p_idxs, heterogeneous_vert_p_idxs)
    replace_vec!(lt_ofun.p_setters, p_setters)
    replace_vec!(lt_ofun.transport_rates, transport_rates)
    replace_vec!(lt_ofun.t_rate_idx_types, t_rate_idx_types)
    lt_ofun.leaving_rates .= leaving_rates

    # Updating the `MTKParameters` structure is a bit more complicated.
    p_dict = Dict(ps_new)
    osys = complete(ode_model(reactionsystem(dsrs)))
    for p in parameters(osys)
        MT.setp(osys, p)(lt_ofun.mtk_ps, (p_dict[p] isa Number) ? p_dict[p] : p_dict[p][1])
    end

    return nothing
end

# Specialised function which replaced one vector in another in a mutating way.
# Required to update the vectors in the `LatticeTransportODEFunction` functor.
function replace_vec!(vec1, vec2)
    l1 = length(vec1)
    l2 = length(vec2)

    # Updates the fields, then deletes superfluous fields, or additional ones.
    for (i, v) in enumerate(vec2[1:min(l1, l2)])
        vec1[i] = v
    end
    foreach(idx -> deleteat!(vec1, idx), l1:-1:(l2 + 1))
    foreach(val -> push!(vec1, val), vec2[(l1 + 1):l2])
end

# Currently not implemented.
function rebuild_spat_internals!(dprob::DiscreteProblem)
    error("Modification and/or rebuilding of `DiscreteProblem`s is currently not supported. Please create a new problem instead.")
end
function rebuild_spat_internals!(jprob::JumpProblem)
    error("Modification and/or rebuilding of `JumpProblem`s is currently not supported. Please create a new problem instead.")
end
function rebuild_spat_internals!(jprob::JumpProcesses.SSAIntegrator)
    error("Modification and/or rebuilding of `JumpProblem` integrators is currently not supported. Please create a new problem instead.")
end

### Utility ###

# Function which extracts a discrete space reaction system's space, taking extra consideration for
# masked grid (converts to sparse array template and checks that it is not 3d).
function extract_dspace(dsrs::DiscreteSpaceReactionSystem)
    dspace = Catalyst.dspace(dsrs)
    if has_masked_dspace(dsrs)
        if grid_dims(dsrs) == 3
            error("The `spat_getu` function is not defined for systems based on 3d sparse arrays. Please raise an issue at the Catalyst GitHub site if this is something which would be useful to you.")
        end
        dspace = sparse(dspace)
    end
    return dspace
end

# Functions which in each sample point reshape the vector of values to the correct form (depending
# on the type of space used).
function reshape_vals(vals, dspace::CartesianGridRej{N, T}) where {N, T}
    return reshape(vals, dspace.dims...)
end
function reshape_vals(vals, dspace::AbstractSparseArray{Bool, Int64, 1})
    return SparseVector(dspace.n, dspace.nzind, vals)
end
function reshape_vals(vals, dspace::AbstractSparseArray{Bool, Int64, 2})
    return SparseMatrixCSC(dspace.m, dspace.n, dspace.colptr, dspace.rowval, vals)
end
function reshape_vals(vals, dspace::DiGraph)
    return vals
end

# Get a species index and the total number of species. Also handles different symbolic forms.
function get_sp_idxs(sp, dsrs::DiscreteSpaceReactionSystem)
    (sp isa Symbol) && (sp = _symbol_to_var(dsrs, sp))
    sp_idx = findfirst(isequal(sp), species(dsrs))
    sp_tot = length(species(dsrs))
    return sp_idx, sp_tot
end

# Get a parameter index and the total number of parameters. Also handles different symbolic forms.
function get_p_idxs(p, dsrs::DiscreteSpaceReactionSystem)
    (p isa Symbol) && (p = _symbol_to_var(dsrs, p))
    p_idx = findfirst(isequal(p), parameters(dsrs))
    p_tot = length(parameters(dsrs))
    return p_idx, p_tot
end

# Checks that a provided input correspond to the correct dspace format.
function check_dspace_format(dspace::CartesianGridRej, u)
    (u isa AbstractArray) ||
        error("The input u should be an AbstractArray. It is a $(typeof(u)).")
    (size(u) == dspace.dims) ||
        error("The input u should have size $(dspace.dims), but has size $(size(u)).")
end
function check_dspace_format(dspace::AbstractSparseArray, u)
    (u isa AbstractArray) ||
        error("The input u should be an AbstractArray. It is a $(typeof(u)).")
    (size(u) == size(dspace)) ||
        error("The input u should have size $(size(dspace)), but has size $(size(u)).")
end
function check_dspace_format(dspace::DiGraph, u)
    (u isa AbstractArray) ||
        error("The input u should be an AbstractVector. It is a $(typeof(u)).")
    (length(u) == nv(dspace)) ||
        error("The input u should have length $(nv(dspace)), but has length $(length(u)).")
end

# Throws an error when interfacing with an edge parameter.
function edge_param_check(p, dsrs)
    (p isa Symbol) && (p = _symbol_to_var(dsrs, p))
    if isedgeparameter(p)
        throw(ArgumentError("The `spat_getp` and `spat_setp!` functions currently does not support edge parameter updating. If you require this functionality, please raise an issue on the Catalyst GitHub page and we can add this feature."))
    end
end
