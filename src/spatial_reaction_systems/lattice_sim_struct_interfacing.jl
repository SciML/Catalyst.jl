### Lattice Simulation Structure Species Getters/Setters ###

"""
    lat_setp!(sim_struct, p, lrs::LatticeReactionSystem, p_val)

For a problem or integrators, update its `p` vector with the input `p_val`. For non-lattice models,
this is can be done through direct interfacing (e.g. `prob[p] = 1.0`). However, for
`LatticeReactionSystem`-based problems and integrators, this function must be used instead.

Arguments:
- `sim_struct`: The simulation structure which `u` value we wish to update. Can be either a `ODEProblem`, `JumpProblem`, or an integrator derived from either of these.
- `p`: The species which value we wish to update. Can be provided either in its symbolic form (e.g. `k`) or as a symbol (e.g. `:k`).
- `lrs`: The `LatticeReactionSystem` which was used to generate the structure we wish to modify.
- `p_val`: The parameter's new values. Must be given in a form which is also a valid initial input to the `ODEProblem`/`JumpProblem`.

Example:
```julia
# Prepare `LatticeReactionSystem`s.
using Catalyst
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
lrs = LatticeReactionSystem(rs, [tr], CartesianGrid((2,3)))

# Prepares a corresponding ODEProblem.
u0 = [:X1 => 1.0, :X2 => 2.0]
tspan = (0.0, 50.0)
ps = [:k1 => [1.0 2.0 3.0; 4.0 5.0 6.0], :k2 => 1.0, :D => 0.01]
oprob = ODEProblem(lrs, u0, tspan, ps)

# Updates the `ODEProblem`.
lat_setp!(oprob, :k1, lrs, 0.0) # Sets `k1` to uniformly 0 across the lattice.
lat_setp!(oprob, :k2, lrs, [1.0 0.0 0.0; 0.0 0.0 0.0]) # Sets `k2` to `1.0` in one vertex, and 0 elsewhere.
```
"""
function lat_setp!(sim_struct, p, lrs::LatticeReactionSystem, p_vals)
    # Error checks.
    (p_vals isa Number) || check_lattice_format(extract_lattice(lrs), p_vals)
    edge_param_check(p, lrs)

    # Converts symbol parameter to symbolic and find the correct species index and numbers.
    (p isa Symbol) && (p = _symbol_to_var(lrs, p))
    (p isa Num) && (p = unwrap(p))
    p_idx, p_tot = get_p_idxs(p, lrs)

    # Reshapes the values to a vector of the correct form, and calls lat_setp! on the input structure.
    p_vals_reshaped = vertex_value_form(p_vals, lrs, p)
    return lat_setp!(sim_struct, p_idx, p_vals_reshaped, num_verts(lrs))
end

# Note: currently, `lat_setp!(oprob::ODEProblem, ...`) and `lat_setp!(SciMLBase.AbstractODEIntegrator, ...`)
# are identical and could be merged to a single function.
function lat_setp!(oprob::ODEProblem, p_idx::Int64, p_vals, num_verts)
    return if length(p_vals) == 1
        foreach(idx -> (oprob.p[p_idx][idx] = p_vals[1]), 1:num_verts)
    elseif length(p_vals) == length(oprob.p[p_idx])
        foreach(idx -> (oprob.p[p_idx][idx] = p_vals[idx]), 1:num_verts)
    elseif length(oprob.p[p_idx]) == 1
        oprob.p[p_idx][1] = p_vals[1]
        foreach(idx -> (push!(oprob.p[p_idx], p_vals[idx])), 2:num_verts)
    end
end
function lat_setp!(jprob::JumpProblem, p_idx::Int64, p_vals, num_verts)
    error("The `lat_setp!` function is currently not supported for `JumpProblem`s.")
end
function lat_setp!(oint::SciMLBase.AbstractODEIntegrator, p_idx::Int64, p_vals, num_verts)
    return if length(p_vals) == 1
        foreach(idx -> (oint.p[p_idx][idx] = p_vals[1]), 1:num_verts)
    elseif length(p_vals) == length(oint.p[p_idx])
        foreach(idx -> (oint.p[p_idx][idx] = p_vals[idx]), 1:num_verts)
    elseif length(oint.p[p_idx]) == 1
        oint.p[p_idx][1] = p_vals[1]
        foreach(idx -> (push!(oint.p[p_idx], p_vals[idx])), 2:num_verts)
    end
end
function lat_setp!(jint::JumpProcesses.SSAIntegrator, p_idx::Int64, p_vals, num_verts)
    error("The `lat_setp!` function is currently not supported for jump simulation integrators.")
end

"""
    lat_getp(sim_struct, p, lrs::LatticeReactionSystem)

For a problem or integrators, retrieves its `p` values. For non-lattice models,
this is can be done through direct interfacing (e.g. `prob[p]`). However, for
`LatticeReactionSystem`-based problems and integrators, this function must be used instead. The
output format depends on the lattice (a dense array for cartesian grid lattices, a sparse array for
masked grid lattices, and a vector for graph lattices). This format is similar to what is used to
designate parameter initial values.

Arguments:
- `sim_struct`: The simulation structure which `p` value we wish to retrieve. Can be either a `ODEProblem`,
`JumpProblem`, or an integrator derived from either of these.
- `p`: The species which value we wish to update. Can be provided either in its symbolic form (e.g. `k`) or as a symbol (e.g. `:k`).
- `lrs`: The `LatticeReactionSystem` which was used to generate the structure we wish to modify.

Notes:
- Even if the parameter is spatially uniform, a full array with its values across all vertices will be retrieved.

Example:
```julia
# Prepare `LatticeReactionSystem`s.
using Catalyst
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
lrs = LatticeReactionSystem(rs, [tr], CartesianGrid((2,3)))

# Prepares a corresponding ODEProblem.
u0 = [:X1 => 1.0, :X2 => 2.0]
tspan = (0.0, 50.0)
ps = [:k1 => [1.0 2.0 3.0; 4.0 5.0 6.0], :k2 => 1.0, :D => 0.01]
oprob = ODEProblem(lrs, u0, tspan, ps)

# Updates the `ODEProblem`.
lat_getp(oprob, :k1, lrs) # Retrieves the value of `k1`.
```
"""
function lat_getp(sim_struct, p, lrs::LatticeReactionSystem)
    edge_param_check(p, lrs)
    p_idx, _ = get_p_idxs(p, lrs)
    return lat_getp(sim_struct, p_idx, extract_lattice(lrs), num_verts(lrs))
end

# Retrieves the lattice values for problem or integrator structures.
function lat_getp(oprob::ODEProblem, p_idx::Int64, lattice, num_verts)
    vals = oprob.p[p_idx]
    (length(vals) == 1) && (vals = fill(vals[1], num_verts))
    return reshape_vals(vals, lattice)
end
function lat_getp(jprob::JumpProblem, p_idx::Int64, lattice, num_verts)
    error("The `lat_getp` function is currently not supported for `JumpProblem`s.")
end
function lat_getp(oint::SciMLBase.AbstractODEIntegrator, p_idx::Int64, lattice, num_verts)
    vals = oint.p[p_idx]
    (length(vals) == 1) && (vals = fill(vals[1], num_verts))
    return reshape_vals(vals, lattice)
end
function lat_getp(jint::JumpProcesses.SSAIntegrator, p_idx::Int64, lattice, num_verts)
    error("The `lat_getp` function is currently not supported for jump simulation integrators.")
end

"""
    lat_setu!(sim_struct, sp, lrs::LatticeReactionSystem, u)

For a problem or integrators, update its `u` vector with the input `u`. For non-lattice models,
this is can be done through direct interfacing (e.g. `prob[X] = 1.0`). However, for
`LatticeReactionSystem`-based problems and integrators, this function must be used instead.

Arguments:
- `sim_struct`: The simulation structure which `u` value we wish to update. Can be either a `ODEProblem`, `JumpProblem`, or an integrator derived from either of these.
- `sp`: The species which value we wish to update. Can be provided either in its symbolic form (e.g. `X`) or as a symbol (e.g. `:X`).
- `lrs`: The `LatticeReactionSystem` which was used to generate the structure we wish to modify.
- `u`: The species's new values. Must be given in a form which is also a valid initial input to the `ODEProblem`/`JumpProblem`.

Example:
```julia
# Prepare `LatticeReactionSystem`s.
using Catalyst
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
lrs = LatticeReactionSystem(rs, [tr], CartesianGrid((2,3)))

# Prepares a corresponding ODEProblem.
u0 = [:X1 => [1.0 2.0 3.0; 4.0 5.0 6.0], :X2 => 2.0]
tspan = (0.0, 50.0)
ps = [:k1 => 2.0, :k2 => 1.0, :D => 0.01]
oprob = ODEProblem(lrs, u0, tspan, ps)

# Updates the `ODEProblem`.
lat_setu!(oprob, :X1, lrs, 0.0) # Sets `X1` to uniformly 0 across the lattice.
lat_setu!(oprob, :X2, lrs, [1.0 0.0 0.0; 0.0 0.0 0.0]) # Sets `X2` to `1.0` in one vertex, and 0 elsewhere.
```
"""
function lat_setu!(sim_struct, sp, lrs::LatticeReactionSystem, u)
    # Checks that if u is non-uniform, it has the correct format for the system's lattice.
    (u isa Number) || check_lattice_format(extract_lattice(lrs), u)

    # Converts symbol species to symbolic and finds correct species index and numbers.
    (sp isa Symbol) && (sp = _symbol_to_var(lrs, sp))
    (sp isa Num) && (sp = unwrap(sp))
    sp_idx, sp_tot = get_sp_idxs(sp, lrs)

    # Reshapes the values to a vector of the correct form, and calls lat_setu! on the input structure.
    u_reshaped = vertex_value_form(u, lrs, sp)
    return lat_setu!(sim_struct, sp_idx, sp_tot, u_reshaped, num_verts(lrs))
end

function lat_setu!(oprob::ODEProblem, sp_idx::Int64, sp_tot::Int64, u, num_verts)
    return if length(u) == 1
        foreach(idx -> (oprob.u0[sp_idx + (idx - 1) * sp_tot] = u[1]), 1:num_verts)
    else
        foreach(idx -> (oprob.u0[sp_idx + (idx - 1) * sp_tot] = u[idx]), 1:num_verts)
    end
end
function lat_setu!(jprob::JumpProblem, sp_idx::Int64, sp_tot::Int64, u, num_verts)
    return if length(u) == 1
        foreach(idx -> (jprob.prob.u0[sp_idx, idx] = u[1]), 1:num_verts)
    else
        foreach(idx -> (jprob.prob.u0[sp_idx, idx] = u[idx]), 1:num_verts)
    end
end
function lat_setu!(
        oint::SciMLBase.AbstractODEIntegrator, sp_idx::Int64, sp_tot::Int64,
        u, num_verts
    )
    return if length(u) == 1
        foreach(idx -> (oint.u[sp_idx + (idx - 1) * sp_tot] = u[1]), 1:num_verts)
    else
        foreach(idx -> (oint.u[sp_idx + (idx - 1) * sp_tot] = u[idx]), 1:num_verts)
    end
end
function lat_setu!(
        jint::JumpProcesses.SSAIntegrator, sp_idx::Int64, sp_tot::Int64,
        u, num_verts
    )
    return if length(u) == 1
        foreach(idx -> (jint.u[sp_idx, idx] = u[1]), 1:num_verts)
    else
        foreach(idx -> (jint.u[sp_idx, idx] = u[idx]), 1:num_verts)
    end
end

"""
    lat_getu(sim_struct, sp, lrs::LatticeReactionSystem)

For a problem or integrators, retrieves its `u` values. For non-lattice models,
this is can be done through direct interfacing (e.g. `prob[X]`). However, for
`LatticeReactionSystem`-based problems and integrators, this function must be used instead. The
output format depends on the lattice (a dense array for cartesian grid lattices, a sparse array for
masked grid lattices, and a vector for graph lattices). This format is similar to which is used to
designate species initial conditions.

Arguments:
- `sim_struct`: The simulation structure which `u` value we wish to retrieve. Can be either a `ODEProblem`, `JumpProblem`, or an integrator derived from either of these.
- `sp`: The species which value we wish to update. Can be provided either in its symbolic form (e.g. `X`) or as a symbol (e.g. `:X`).
- `lrs`: The `LatticeReactionSystem` which was used to generate the structure we wish to modify.

Notes:
- Even if the species is spatially uniform, a full array with its values across all vertices will be retrieved.

Example:
```julia
# Prepare `LatticeReactionSystem`s.
using Catalyst
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
lrs = LatticeReactionSystem(rs, [tr], CartesianGrid((2,3)))

# Prepares a corresponding ODEProblem.
u0 = [:X1 => [1.0 2.0 3.0; 4.0 5.0 6.0], :X2 => 2.0]
tspan = (0.0, 50.0)
ps = [:k1 => 2.0, :k2 => 1.0, :D => 0.01]
oprob = ODEProblem(lrs, u0, tspan, ps)

# Updates the `ODEProblem`.
lat_getu(oprob, :X1, lrs) # Retrieves the value of `X1`.
```
"""
function lat_getu(sim_struct, sp, lrs::LatticeReactionSystem)
    sp_idx, sp_tot = get_sp_idxs(sp, lrs)
    return lat_getu(sim_struct, sp_idx, sp_tot, extract_lattice(lrs))
end

# Retrieves the lattice values for problem or integrator structures.
function lat_getu(oprob::ODEProblem, sp_idx, sp_tot, lattice)
    return reshape_vals(oprob.u0[sp_idx:sp_tot:end], lattice)
end
function lat_getu(jprob::JumpProblem, sp_idx, sp_tot, lattice)
    return reshape_vals(jprob.prob.u0[sp_idx, :], lattice)
end
function lat_getu(oint::SciMLBase.AbstractODEIntegrator, sp_idx, sp_tot, lattice)
    return reshape_vals(oint.u[sp_idx:sp_tot:end], lattice)
end
function lat_getu(jint::JumpProcesses.SSAIntegrator, sp_idx, sp_tot, lattice)
    return reshape_vals(jint.u[sp_idx, :], lattice)
end

# A single function, `lat_getu`, which contains all interfacing functionality. However,
# long-term it should be replaced with a sleeker interface. Ideally as MTK-wide support for
# lattice problems and solutions is introduced. Note that SciML considers jump simulation solutions
# as `ODESolution`, hence that type is specified.
"""
    lat_getu(sol, sp, lrs::LatticeReactionSystem; t = nothing)

A function for retrieving the solution of a `LatticeReactionSystem`-based simulation on various
desired forms. Generally, for `LatticeReactionSystem`s, the values in `sol` is ordered in a
way which is not directly interpretable by the user. Furthermore, the normal Catalyst interface
for solutions (e.g. `sol[:X]`) does not work for these solutions. Hence this function is used instead.

The output is a vector, which in each position contains sp's value (either at a time step of time,
depending on the input `t`). Its shape depends on the lattice (using a similar form as heterogeneous
initial conditions). I.e. for a NxM cartesian grid, the values are NxM matrices. For a masked grid,
the values are sparse matrices. For a graph lattice, the values are vectors (where the value in
the n'th position corresponds to sp's value in the n'th vertex).

Arguments:
- `sol`: The solution from which we wish to retrieve some values.
- `sp`: The species which value we wish to update. Can be provided either in its symbolic form (e.g. `X`) or as a symbol (e.g. `:X`).
- `lrs`: The `LatticeReactionSystem` which was simulated to generate the solution.
- `t = nothing`: If `nothing`, we simply return the solution across all saved time steps (default). If `t` instead is a vector (or range of values), returns the solution interpolated at these time points.

Notes:
- The `lat_getu` is not optimised for performance. However, it should still be quite performant, but there might be some limitations if called a very large number of times.
- Long-term it is likely that this function gets replaced with a sleeker interface.

Example:
```julia
using Catalyst, OrdinaryDiffEqDefault

# Prepare `LatticeReactionSystem`s.
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1
lrs = LatticeReactionSystem(rs, [tr], CartesianGrid((2,2)))

# Create problems.
u0 = [:X1 => 1, :X2 => 2]
tspan = (0.0, 10.0)
ps = [:k1 => 1, :k2 => 2.0, :D => 0.1]

oprob = ODEProblem(lrs1, u0, tspan, ps)
osol = solve(oprob)
lat_getu(osol, :X1, lrs) # Returns the value of X1 at each time step.
lat_getu(osol, :X1, lrs; t = 0.0:10.0) # Returns the value of X1 at times 0.0, 1.0, ..., 10.0
```
"""
function lat_getu(sol::ODESolution, sp, lrs::LatticeReactionSystem; t = nothing)
    (t isa Number) && error("The input `t` to `lat_getu` must be a `AbstractVector`.")
    sp_idx, sp_tot = get_sp_idxs(sp, lrs)
    return lat_getu(sol, extract_lattice(lrs), t, sp_idx, sp_tot)
end

# Function which handles the input in the case where `t` is `nothing` (i.e. return `sp`s value
# across all sample points).
function lat_getu(sol::ODESolution, lattice, t::Nothing, sp_idx::Int64, sp_tot::Int64)
    # ODE simulations contain, in each data point, all values in a single vector. Jump simulations
    # instead in a matrix (NxM, where N is the number of species and M is the number of vertices). We
    # must consider each case separately.
    if sol.prob isa ODEProblem
        return [reshape_vals(vals[sp_idx:sp_tot:end], lattice) for vals in sol.u]
    elseif sol.prob isa DiscreteProblem
        return [reshape_vals(vals[sp_idx, :], lattice) for vals in sol.u]
    else
        error("Unknown type of solution provided to `lat_getu`. Only ODE or Jump solutions are supported.")
    end
end

# Function which handles the input in the case where `t` is a range of values (i.e. return `sp`s
# value at all designated time points.
function lat_getu(
        sol::ODESolution, lattice, t::AbstractVector{T}, sp_idx::Int64,
        sp_tot::Int64
    ) where {T <: Number}
    # Checks that an appropriate `t` is provided (however, DiffEq does permit out-of-range `t`s).
    if (minimum(t) < sol.t[1]) || (maximum(t) > sol.t[end])
        error("The range of the t values provided for sampling, ($(minimum(t)),$(maximum(t))) is not fully within the range of the simulation time span ($(sol.t[1]),$(sol.t[end])).")
    end

    # ODE simulations contain, in each data point, all values in a single vector. Jump simulations
    # instead in a matrix (NxM, where N is the number of species and M is the number of vertices). We
    # must consider each case separately.
    if sol.prob isa ODEProblem
        return [reshape_vals(sol(ti)[sp_idx:sp_tot:end], lattice) for ti in t]
    elseif sol.prob isa DiscreteProblem
        return [reshape_vals(sol(ti)[sp_idx, :], lattice) for ti in t]
    else
        error("Unknown type of solution provided to `lat_getu`. Only ODE or Jump solutions are supported.")
    end
end

### Lattice Internal Rebuilding ###

"""
    rebuild_lat_internals!(sciml_struct)

Rebuilds the internal functions for simulating a LatticeReactionSystem. Whenever a problem or
integrator has had its parameter values updated, this function should be called for the update to
be taken into account. For ODE simulations, `rebuild_lat_internals!` needs only to be called when
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
lrs = LatticeReactionSystem(rs, [tr], grid)

u0 = [:X1 => 2, :X2 => [5 6; 7 8]]
tspan = (0.0, 10.0)
ps = [:k1 => 1.5, :k2 => [1.0 1.5; 2.0 3.5], :D => 0.1]

oprob = ODEProblem(lrs, u0, tspan, ps)

# Updates parameter values.
oprob.ps[:ks] = [2.0 2.5; 3.0 4.5]
oprob.ps[:D] = 0.05

# Rebuilds `ODEProblem` to make changes have an effect.
rebuild_lat_internals!(oprob)
```
"""
function rebuild_lat_internals!(oprob::ODEProblem)
    return rebuild_lat_internals!(oprob.f.f, oprob.p, oprob.f.f.lrs)
end

# Function for rebuilding a `LatticeReactionSystem` integrator after it has been updated.
function rebuild_lat_internals!(integrator::SciMLBase.AbstractODEIntegrator)
    return rebuild_lat_internals!(integrator.f.f, integrator.p, integrator.f.f.lrs)
end

# Function which rebuilds a `LatticeTransportODEFunction` functor for a new parameter set.
function rebuild_lat_internals!(
        lt_ofun::LatticeTransportODEFunction, ps_new,
        lrs::LatticeReactionSystem
    )
    # Computes Jacobian properties.
    jac = !isnothing(lt_ofun.jac_transport)
    sparse = lt_ofun.sparse

    # Recreates the new parameters on the requisite form.
    ps_new = [(length(p) == 1) ? p[1] : p for p in deepcopy(ps_new)]
    ps_new = [p => p_val for (p, p_val) in zip(parameters(lrs), deepcopy(ps_new))]
    vert_ps,
        edge_ps = lattice_process_p(
        ps_new, vertex_parameters(lrs),
        edge_parameters(lrs), lrs
    )
    ps_new = [vert_ps; edge_ps]

    # Creates the new transport rates and transport Jacobian part.
    transport_rates = make_sidxs_to_transrate_map(vert_ps, edge_ps, lrs)
    if !isnothing(lt_ofun.jac_transport)
        lt_ofun.jac_transport .= 0.0
        set_jac_transport_values!(lt_ofun.jac_transport, transport_rates, lrs)
    end

    # Computes new field values.
    heterogeneous_vert_p_idxs = make_heterogeneous_vert_p_idxs(ps_new, lrs)
    mtk_ps, p_setters = make_mtk_ps_structs(ps_new, lrs, heterogeneous_vert_p_idxs)
    t_rate_idx_types, leaving_rates = make_t_types_and_leaving_rates(transport_rates, lrs)

    # Updates functor fields.
    replace_vec!(lt_ofun.heterogeneous_vert_p_idxs, heterogeneous_vert_p_idxs)
    replace_vec!(lt_ofun.p_setters, p_setters)
    replace_vec!(lt_ofun.transport_rates, transport_rates)
    replace_vec!(lt_ofun.t_rate_idx_types, t_rate_idx_types)
    lt_ofun.leaving_rates .= leaving_rates

    # Updating the `MTKParameters` structure is a bit more complicated.
    p_dict = Dict(ps_new)
    osys = complete(ode_model(reactionsystem(lrs)))
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
    return foreach(val -> push!(vec1, val), vec2[(l1 + 1):l2])
end

# Currently not implemented.
function rebuild_lat_internals!(dprob::DiscreteProblem)
    error("Modification and/or rebuilding of `DiscreteProblem`s is currently not supported. Please create a new problem instead.")
end
function rebuild_lat_internals!(jprob::JumpProblem)
    error("Modification and/or rebuilding of `JumpProblem`s is currently not supported. Please create a new problem instead.")
end
function rebuild_lat_internals!(jprob::JumpProcesses.SSAIntegrator)
    error("Modification and/or rebuilding of `JumpProblem` integrators is currently not supported. Please create a new problem instead.")
end

### Utility ###

# Function which extracts a lattice reaction system's lattice, taking extra consideration for
# masked grid (converts to sparse array template and checks that it is not 3d).
function extract_lattice(lrs::LatticeReactionSystem)
    lattice = Catalyst.lattice(lrs)
    if has_masked_lattice(lrs)
        if grid_dims(lrs) == 3
            error("The `lat_getu` function is not defined for systems based on 3d sparse arrays. Please raise an issue at the Catalyst GitHub site if this is something which would be useful to you.")
        end
        lattice = sparse(lattice)
    end
    return lattice
end

# Functions which in each sample point reshape the vector of values to the correct form (depending
# on the type of lattice used).
function reshape_vals(vals, lattice::CartesianGridRej{N, T}) where {N, T}
    return reshape(vals, lattice.dims...)
end
function reshape_vals(vals, lattice::AbstractSparseArray{Bool, Int64, 1})
    return SparseVector(lattice.n, lattice.nzind, vals)
end
function reshape_vals(vals, lattice::AbstractSparseArray{Bool, Int64, 2})
    return SparseMatrixCSC(lattice.m, lattice.n, lattice.colptr, lattice.rowval, vals)
end
function reshape_vals(vals, lattice::DiGraph)
    return vals
end

# Get a species index and the total number of species. Also handles different symbolic forms.
function get_sp_idxs(sp, lrs::LatticeReactionSystem)
    (sp isa Symbol) && (sp = _symbol_to_var(lrs, sp))
    sp_idx = findfirst(isequal(sp), species(lrs))
    sp_tot = length(species(lrs))
    return sp_idx, sp_tot
end

# Get a parameter index and the total number of parameters. Also handles different symbolic forms.
function get_p_idxs(p, lrs::LatticeReactionSystem)
    (p isa Symbol) && (p = _symbol_to_var(lrs, p))
    p_idx = findfirst(isequal(p), parameters(lrs))
    p_tot = length(parameters(lrs))
    return p_idx, p_tot
end

# Checks that a provided input correspond to the correct lattice format.
function check_lattice_format(lattice::CartesianGridRej, u)
    (u isa AbstractArray) ||
        error("The input u should be an AbstractArray. It is a $(typeof(u)).")
    return (size(u) == lattice.dims) ||
        error("The input u should have size $(lattice.dims), but has size $(size(u)).")
end
function check_lattice_format(lattice::AbstractSparseArray, u)
    (u isa AbstractArray) ||
        error("The input u should be an AbstractArray. It is a $(typeof(u)).")
    return (size(u) == size(lattice)) ||
        error("The input u should have size $(size(lattice)), but has size $(size(u)).")
end
function check_lattice_format(lattice::DiGraph, u)
    (u isa AbstractArray) ||
        error("The input u should be an AbstractVector. It is a $(typeof(u)).")
    return (length(u) == nv(lattice)) ||
        error("The input u should have length $(nv(lattice)), but has length $(length(u)).")
end

# Throws an error when interfacing with an edge parameter.
function edge_param_check(p, lrs)
    (p isa Symbol) && (p = _symbol_to_var(lrs, p))
    return if isedgeparameter(p)
        throw(ArgumentError("The `lat_getp` and `lat_setp!` functions currently does not support edge parameter updating. If you require this functionality, please raise an issue on the Catalyst GitHub page and we can add this feature."))
    end
end
