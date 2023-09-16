### JumpProblem ###

# Builds a spatial DiscreteProblem from a Lattice Reaction System.

# Creates a DiscreteProblem from a LatticeReactionSystem.
function DiffEqBase.DiscreteProblem(lrs::LatticeReactionSystem, u0_in, tspan, p_in = DiffEqBase.NullParameters(), args...; kwargs...)
    is_transport_system(lrs) || error("Currently lattice Jump simulations only supported when all spatial reactions are transport reactions.")

    # Converts potential symmaps to varmaps.
    u0_in = symmap_to_varmap(lrs, u0_in)
    p_in = (p_in isa Tuple{<:Any,<:Any}) ? (symmap_to_varmap(lrs, p_in[1]),symmap_to_varmap(lrs, p_in[2])) : symmap_to_varmap(lrs, p_in)

    # Converts u0 and p to Vector{Vector{Float64}} form.
    u0 = lattice_process_u0(u0_in, species(lrs), lrs.nV)
    pC, pD = lattice_process_p(p_in, vertex_parameters(lrs), edge_parameters(lrs), lrs)
    
    # Creates DiscreteProblem.
    return DiscreteProblem(lrs.rs, u0, tspan, (pC, pD), args...; kwargs...)
end

# Builds a spatial JumpProblem from a DiscreteProblem containg a Lattice Reaction System.
function JumpProcesses.JumpProblem(lrs::LatticeReactionSystem, dprob, aggregator, args...; name = nameof(lrs.rs), combinatoric_ratelaws = get_combinatoric_ratelaws(lrs.rs),checks = false, kwargs...)
    # Error checks.
    dprob.p isa Tuple{Vector{Vector{Float64}}, Vector{Vector{Float64}}} || error("Parameters in input DiscreteProblem is of an unexpected type: $(typeof(dprob.p)). Was a LatticeReactionProblem passed into the DiscreteProblem when it was created?")
    any(length.(dprob.p[1]) .> 1) && error("Spatial reaction rates are currently not supported in lattice jump simulations.")
    
    # Creates JumpProblem.
    hopping_constants = make_hopping_constants(dprob, lrs)
    ___dprob = DiscreteProblem(reshape(dprob.u0, lrs.nS, lrs.nV), dprob.tspan, first.(dprob.p[1]))
    majumps_ = make_majumps(___dprob, lrs.rs)
    return JumpProblem(___dprob, aggregator, majumps_, hopping_constants = hopping_constants, spatial_system = lrs.lattice, save_positions = (true, false))
end

# Creates the hopping constants from a discrete problem and a lattice reaction system.
function make_hopping_constants(dprob::DiscreteProblem, lrs::LatticeReactionSystem)
    spatial_rates_dict = Dict(compute_all_spatial_rates(dprob.p[1], dprob.p[2], lrs))
    all_diff_rates = [haskey(spatial_rates_dict, s) ? spatial_rates_dict[s] : [0.0] for s in species(lrs)]
    hopping_constants = [Vector{Float64}(undef, length(lrs.lattice.fadjlist[j])) for i in 1:(lrs.nS), j in 1:(lrs.nV)]
    for (e_idx, e) in enumerate(edges(lrs.lattice)), s_idx in 1:(lrs.nS)
        dst_idx = findfirst(isequal(e.dst), lrs.lattice.fadjlist[e.src])
        hopping_constants[s_idx, e.src][dst_idx] = get_component_value(all_diff_rates[s_idx], e_idx)
    end
    return hopping_constants
end

# Creates the mass action jumps from a discrete problem and a reaction system.
function make_majumps(non_spat_dprob, rs::ReactionSystem)
    js = convert(JumpSystem, rs)
    statetoid = Dict(ModelingToolkit.value(state) => i for (i, state) in enumerate(states(rs)))
    eqs = equations(js)
    invttype = non_spat_dprob.tspan[1] === nothing ? Float64 : typeof(1 / non_spat_dprob.tspan[2])
    p = (non_spat_dprob.p isa DiffEqBase.NullParameters || non_spat_dprob.p === nothing) ? Num[] : non_spat_dprob.p    
    majpmapper = ModelingToolkit.JumpSysMajParamMapper(js, p; jseqs = eqs, rateconsttype = invttype)
    return ModelingToolkit.assemble_maj(eqs.x[1], statetoid, majpmapper)
end