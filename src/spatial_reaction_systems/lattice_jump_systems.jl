### JumpProblem ###

# Builds a spatial DiscreteProblem from a Lattice Reaction System.
function DiffEqBase.DiscreteProblem(lrs::LatticeReactionSystem, u0_in, tspan,
        p_in = DiffEqBase.NullParameters(), args...; kwargs...)
    if !is_transport_system(lrs)
        error("Currently lattice Jump simulations only supported when all spatial reactions are transport reactions.")
    end

    # Converts potential symmaps to varmaps.
    u0_in = symmap_to_varmap(lrs, u0_in)
    p_in = symmap_to_varmap(lrs, p_in)

    # Converts u0 and p to their internal forms.
    # u0 is simply a vector with all the species' initial condition values across all vertices.
    # u0 is [spec 1 at vert 1, spec 2 at vert 1, ..., spec 1 at vert 2, ...].
    u0 = lattice_process_u0(u0_in, species(lrs), lrs)
    # vert_ps and `edge_ps` are vector maps, taking each parameter's Symbolics representation to its value(s).
    # vert_ps values are vectors. Here, index (i) is a parameter's value in vertex i.
    # edge_ps values are sparse matrices. Here, index (i,j) is a parameter's value in the edge from vertex i to vertex j.
    # Uniform vertex/edge parameters store only a single value (a length 1 vector, or size 1x1 sparse matrix).
    vert_ps,
    edge_ps = lattice_process_p(p_in, vertex_parameters(lrs),
        edge_parameters(lrs), lrs)

    # Returns a DiscreteProblem (which basically just stores the processed input).
    return DiscreteProblem(u0, tspan, [vert_ps; edge_ps], args...; kwargs...)
end

# Builds a spatial JumpProblem from a DiscreteProblem containing a `LatticeReactionSystem`.
function JumpProcesses.JumpProblem(lrs::LatticeReactionSystem, dprob, aggregator, args...;
        combinatoric_ratelaws = get_combinatoric_ratelaws(reactionsystem(lrs)),
        name = nameof(reactionsystem(lrs)), kwargs...)
    # Error checks.
    if !isnothing(dprob.f.sys)
        throw(ArgumentError("Unexpected `DiscreteProblem` passed into `JumpProblem`. Was a `LatticeReactionSystem` used as input to the initial `DiscreteProblem`?"))
    end

    # Computes hopping constants and mass action jumps (requires some internal juggling).
    # Currently, the resulting JumpProblem does not depend on parameters (no way to incorporate these).
    # Hence the parameters of this one do not actually matter. If at some point JumpProcess can
    # handle parameters this can be updated and improved.
    # The non-spatial DiscreteProblem have a u0 matrix with entries for all combinations of species and vertices.
    hopping_constants = make_hopping_constants(dprob, lrs)
    sma_jumps = make_spatial_majumps(dprob, lrs)
    non_spat_dprob = DiscreteProblem(reshape(dprob.u0, num_species(lrs), num_verts(lrs)),
        dprob.tspan, first.(dprob.p[1]))

    # Creates and returns a spatial JumpProblem (masked lattices are not supported by these).
    spatial_system = has_masked_lattice(lrs) ? get_lattice_graph(lrs) : lattice(lrs)
    return JumpProblem(non_spat_dprob, aggregator, sma_jumps;
        hopping_constants, spatial_system, name, kwargs...)
end

# Creates the hopping constants from a discrete problem and a lattice reaction system.
function make_hopping_constants(dprob::DiscreteProblem, lrs::LatticeReactionSystem)
    # Creates the all_diff_rates vector, containing for each species, its transport rate across all edges.
    # If the transport rate is uniform for one species, the vector has a single element, else one for each edge.
    spatial_rates_dict = Dict(compute_all_transport_rates(Dict(dprob.p), lrs))
    all_diff_rates = [haskey(spatial_rates_dict, s) ? spatial_rates_dict[s] : [0.0]
                      for s in species(lrs)]

    # Creates an array (of the same size as the hopping constant array) containing all edges.
    # First the array is a NxM matrix (number of species x number of vertices). Each element is a
    # vector containing all edges leading out from that vertex (sorted by destination index).
    edge_array = [Pair{Int64, Int64}[] for _1 in 1:num_species(lrs), _2 in 1:num_verts(lrs)]
    for e in edge_iterator(lrs), s_idx in 1:num_species(lrs)

        push!(edge_array[s_idx, e[1]], e)
    end
    foreach(e_vec -> sort!(e_vec; by = e -> e[2]), edge_array)

    # Creates the hopping constants array. It has the same shape as the edge array, but each
    # element is that species transportation rate along that edge
    hopping_constants = [[Catalyst.get_edge_value(all_diff_rates[s_idx], e)
                          for e in edge_array[s_idx, src_idx]]
                         for s_idx in 1:num_species(lrs), src_idx in 1:num_verts(lrs)]
    return hopping_constants
end

# Creates a SpatialMassActionJump struct from a (spatial) DiscreteProblem and a LatticeReactionSystem.
# Could implement a version which, if all reactions' rates are uniform, returns a MassActionJump.
# Not sure if there is any form of performance improvement from that though. Likely not the case.
function make_spatial_majumps(dprob, lrs::LatticeReactionSystem)
    # Creates a vector, storing which reactions have spatial components.
    is_spatials = [has_spatial_vertex_component(rx.rate, dprob.p)
                   for rx in reactions(reactionsystem(lrs))]

    # Creates templates for the rates (uniform and spatial) and the stoichiometries.
    # We cannot fetch reactant_stoich and net_stoich from a (non-spatial) MassActionJump.
    # The reason is that we need to re-order the reactions so that uniform appears first, and spatial next.
    num_rxs = length(reactions(reactionsystem(lrs)))
    u_rates = Vector{Float64}(undef, num_rxs - count(is_spatials))
    s_rates = Matrix{Float64}(undef, count(is_spatials), num_verts(lrs))
    reactant_stoich = Vector{Vector{Pair{Int64, Int64}}}(undef, num_rxs)
    net_stoich = Vector{Vector{Pair{Int64, Int64}}}(undef, num_rxs)

    # Loops through reactions with non-spatial rates, computes their rates and stoichiometries.
    cur_rx = 1
    for (is_spat, rx) in zip(is_spatials, reactions(reactionsystem(lrs)))
        is_spat && continue
        u_rates[cur_rx] = compute_vertex_value(rx.rate, lrs; ps = dprob.p)[1]
        substoich_map = Pair.(rx.substrates, rx.substoich)
        reactant_stoich[cur_rx] = int_map(substoich_map, reactionsystem(lrs))
        net_stoich[cur_rx] = int_map(rx.netstoich, reactionsystem(lrs))
        cur_rx += 1
    end
    # Loops through reactions with spatial rates, computes their rates and stoichiometries.
    for (is_spat, rx) in zip(is_spatials, reactions(reactionsystem(lrs)))
        is_spat || continue
        s_rates[cur_rx - length(u_rates), :] .= compute_vertex_value(rx.rate, lrs;
            ps = dprob.p)
        substoich_map = Pair.(rx.substrates, rx.substoich)
        reactant_stoich[cur_rx] = int_map(substoich_map, reactionsystem(lrs))
        net_stoich[cur_rx] = int_map(rx.netstoich, reactionsystem(lrs))
        cur_rx += 1
    end
    # SpatialMassActionJump expects empty rate containers to be nothing.
    isempty(u_rates) && (u_rates = nothing)
    (count(is_spatials) == 0) && (s_rates = nothing)

    return SpatialMassActionJump(u_rates, s_rates, reactant_stoich, net_stoich, nothing)
end

### Extra ###

# Temporary. Awaiting implementation in SII, or proper implementation within Catalyst (with
# more general functionality).
function int_map(map_in, sys)
    return [MT.variable_index(sys, pair[1]) => pair[2] for pair in map_in]
end

# Currently unused. If we want to create certain types of MassActionJumps (instead of SpatialMassActionJumps) we can take this one back.
# Creates the (non-spatial) mass action jumps from a (non-spatial) DiscreteProblem (and its Reaction System of origin).
# function make_majumps(non_spat_dprob, rs::ReactionSystem)
#     # Computes various required inputs for assembling the mass action jumps.
#     js = make_sck_jump(rs)
#     statetoid = Dict(MT.value(state) => i for (i, state) in enumerate(unknowns(rs)))
#     eqs = equations(js)
#     invttype = non_spat_dprob.tspan[1] === nothing ? Float64 : typeof(1 / non_spat_dprob.tspan[2])
#
#     # Assembles the non-spatial mass action jumps.
#     p = (non_spat_dprob.p isa DiffEqBase.NullParameters || non_spat_dprob.p === nothing) ? Num[] : non_spat_dprob.p
#     majpmapper = MT.JumpSysMajParamMapper(js, p; jseqs = eqs, rateconsttype = invttype)
#     return MT.assemble_maj(eqs.x[1], statetoid, majpmapper)
# end
