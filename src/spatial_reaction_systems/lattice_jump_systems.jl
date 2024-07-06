### JumpProblem ###

# Builds a spatial DiscreteProblem from a Lattice Reaction System.
function DiffEqBase.DiscreteProblem(lrs::LatticeReactionSystem, u0_in, tspan,
        p_in = DiffEqBase.NullParameters(), args...; kwargs...)
    if !is_transport_system(lrs)
        error("Currently lattice Jump simulations only supported when all spatial reactions are transport reactions.")
    end

    # Converts potential symmaps to varmaps
    # Vertex and edge parameters may be given in a tuple, or in a common vector, making parameter case complicated.
    u0_in = symmap_to_varmap(lrs, u0_in)
    p_in = (p_in isa Tuple{<:Any, <:Any}) ?
           (symmap_to_varmap(lrs, p_in[1]), symmap_to_varmap(lrs, p_in[2])) :
           symmap_to_varmap(lrs, p_in)

    # Converts u0 and p to their internal forms.
    # u0 is [spec 1 at vert 1, spec 2 at vert 1, ..., spec 1 at vert 2, ...].
    u0 = lattice_process_u0(u0_in, species(lrs), lrs.num_verts)
    # Both vert_ps and edge_ps becomes vectors of vectors. Each have 1 element for each parameter.
    # These elements are length 1 vectors (if the parameter is uniform),
    # or length num_verts/nE, with unique values for each vertex/edge (for vert_ps/edge_ps, respectively).
    vert_ps, edge_ps = lattice_process_p(p_in, vertex_parameters(lrs),
        edge_parameters(lrs), lrs)

    # Returns a DiscreteProblem.
    # Previously, a Tuple was used for (vert_ps, edge_ps), but this was converted to a Vector internally.
    return DiscreteProblem(u0, tspan, [vert_ps, edge_ps], args...; kwargs...)
end

# Builds a spatial JumpProblem from a DiscreteProblem containing a Lattice Reaction System.
function JumpProcesses.JumpProblem(lrs::LatticeReactionSystem, dprob, aggregator,
        args...; name = nameof(lrs.rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(lrs.rs), kwargs...)
    # Error checks.
    if !isnothing(dprob.f.sys)
        error("Unexpected `DiscreteProblem` passed into `JumpProblem`. Was a `LatticeReactionSystem` used as input to the initial `DiscreteProblem`?")
    end

    # Computes hopping constants and mass action jumps (requires some internal juggling).
    # Currently, JumpProcesses requires uniform vertex parameters (hence `p=first.(dprob.p[1])`).
    # Currently, the resulting JumpProblem does not depend on parameters (no way to incorporate these).
    # Hence the parameters of this one does nto actually matter. If at some point JumpProcess can
    # handle parameters this can be updated and improved.
    # The non-spatial DiscreteProblem have a u0 matrix with entries for all combinations of species and vertexes.
    hopping_constants = make_hopping_constants(dprob, lrs)
    sma_jumps = make_spatial_majumps(dprob, lrs)
    non_spat_dprob = DiscreteProblem(
        reshape(dprob.u0, lrs.num_species, lrs.num_verts), dprob.tspan, first.(dprob.p[1]))

    return JumpProblem(non_spat_dprob, aggregator, sma_jumps;
        hopping_constants, spatial_system = lrs.lattice, name, kwargs...)
end

# Creates the hopping constants from a discrete problem and a lattice reaction system.
function make_hopping_constants(dprob::DiscreteProblem, lrs::LatticeReactionSystem)
    # Creates the all_diff_rates vector, containing for each species, its transport rate across all edges.
    # If transport rate is uniform for one species, the vector have a single element, else one for each edge.
    spatial_rates_dict = Dict(compute_all_transport_rates(dprob.p[1], dprob.p[2], lrs))
    all_diff_rates = [haskey(spatial_rates_dict, s) ? spatial_rates_dict[s] : [0.0]
                      for s in species(lrs)]

    # Creates the hopping constant Matrix. It contains one element for each combination of species and vertex.
    # Each element is a Vector, containing the outgoing hopping rates for that species, from that vertex, on that edge.
    hopping_constants = [Vector{Float64}(undef, length(lrs.lattice.fadjlist[j]))
                         for i in 1:(lrs.num_species), j in 1:(lrs.num_verts)]

    # For each edge, finds each position in `hopping_constants`.
    for (e_idx, e) in enumerate(edges(lrs.lattice))
        dst_idx = findfirst(isequal(e.dst), lrs.lattice.fadjlist[e.src])
        # For each species, sets that hopping rate.
        for s_idx in 1:(lrs.num_species)
            hopping_constants[s_idx, e.src][dst_idx] = get_component_value(
                all_diff_rates[s_idx], e_idx)
        end
    end

    return hopping_constants
end

# Creates a SpatialMassActionJump struct from a (spatial) DiscreteProblem and a LatticeReactionSystem.
# Could implementation a version which, if all reaction's rates are uniform, returns a MassActionJump.
# Not sure if there is any form of performance improvement from that though. Possibly is not the case.
function make_spatial_majumps(dprob, lrs::LatticeReactionSystem)
    # Creates a vector, storing which reactions have spatial components.
    is_spatials = [Catalyst.has_spatial_vertex_component(rx.rate, lrs;
                       vert_ps = dprob.p[1]) for rx in reactions(lrs.rs)]

    # Creates templates for the rates (uniform and spatial) and the stoichiometries.
    # We cannot fetch reactant_stoich and net_stoich from a (non-spatial) MassActionJump.
    # The reason is that we need to re-order the reactions so that uniform appears first, and spatial next.
    u_rates = Vector{Float64}(undef, length(reactions(lrs.rs)) - count(is_spatials))
    s_rates = Matrix{Float64}(undef, count(is_spatials), lrs.num_verts)
    reactant_stoich = Vector{Vector{Pair{Int64, Int64}}}(undef, length(reactions(lrs.rs)))
    net_stoich = Vector{Vector{Pair{Int64, Int64}}}(undef, length(reactions(lrs.rs)))

    # Loops through reactions with non-spatial rates, computes their rates and stoichiometries.
    cur_rx = 1
    for (is_spat, rx) in zip(is_spatials, reactions(lrs.rs))
        is_spat && continue
        u_rates[cur_rx] = compute_vertex_value(rx.rate, lrs; vert_ps = dprob.p[1])[1]
        substoich_map = Pair.(rx.substrates, rx.substoich)
        reactant_stoich[cur_rx] = int_map(substoich_map, lrs.rs)
        net_stoich[cur_rx] = int_map(rx.netstoich, lrs.rs)
        cur_rx += 1
    end
    # Loops through reactions with spatial rates, computes their rates and stoichiometries.
    for (is_spat, rx) in zip(is_spatials, reactions(lrs.rs))
        is_spat || continue
        s_rates[cur_rx - length(u_rates), :] = compute_vertex_value(rx.rate, lrs;
            vert_ps = dprob.p[1])
        substoich_map = Pair.(rx.substrates, rx.substoich)
        reactant_stoich[cur_rx] = int_map(substoich_map, lrs.rs)
        net_stoich[cur_rx] = int_map(rx.netstoich, lrs.rs)
        cur_rx += 1
    end
    # SpatialMassActionJump expects empty rate containers to be nothing.
    isempty(u_rates) && (u_rates = nothing)
    (count(is_spatials) == 0) && (s_rates = nothing)

    return SpatialMassActionJump(u_rates, s_rates, reactant_stoich, net_stoich)
end

### Extra ###

# Temporary. Awaiting implementation in SII, or proper implementation withinCatalyst (with more general functionality).
function int_map(map_in, sys)
    return [ModelingToolkit.variable_index(sys, pair[1]) => pair[2] for pair in map_in]
end

# Currently unused. If we want to create certain types of MassActionJumps (instead of SpatialMassActionJumps) we can take this one back.
# Creates the (non-spatial) mass action jumps from a (non-spatial) DiscreteProblem (and its Reaction System of origin).
# function make_majumps(non_spat_dprob, rs::ReactionSystem)
#     # Computes various required inputs for assembling the mass action jumps.
#     js = convert(JumpSystem, rs)
#     statetoid = Dict(ModelingToolkit.value(state) => i for (i, state) in enumerate(states(rs)))
#     eqs = equations(js)
#     invttype = non_spat_dprob.tspan[1] === nothing ? Float64 : typeof(1 / non_spat_dprob.tspan[2])
#
#     # Assembles the non-spatial mass action jumps.
#     p = (non_spat_dprob.p isa DiffEqBase.NullParameters || non_spat_dprob.p === nothing) ? Num[] : non_spat_dprob.p
#     majpmapper = ModelingToolkit.JumpSysMajParamMapper(js, p; jseqs = eqs, rateconsttype = invttype)
#     return ModelingToolkit.assemble_maj(eqs.x[1], statetoid, majpmapper)
# end
