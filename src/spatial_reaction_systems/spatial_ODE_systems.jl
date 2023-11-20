### Spatial ODE Functor Structures ###

# Functor structure containing the information for the forcing function of a spatial ODE with spatial movement on a lattice.
struct LatticeTransportODEf{R,S,T}
    """The ODEFunction of the (non-spatial) reaction system which generated this function."""
    ofunc::R
    """The number of vertices."""
    num_verts::Int64
    """The number of species."""
    num_species::Int64
    """The values of the parameters which values are tied to vertexes."""
    vert_ps::Vector{Vector{S}}
    """Temporary vector. For parameters which values are identical across the lattice, at some point these have to be converted of a length num_verts vector. To avoid re-allocation they are written to this vector."""
    work_vert_ps::Vector{S}    
    """For each parameter in vert_ps, its value is a vector with length either num_verts or 1. To know whenever a parameter's value need expanding to the work_vert_ps array, its length needs checking. This check is done once, and the value stored to this array. This field (specifically) is an enumerate over that array."""
    v_ps_idx_types::Vector{Bool}
    """A vector of pairs, with a value for each species with transportation. The first value is the species index (in the species(::ReactionSystem) vector), and the second is a vector with its transport rate values. If the transport rate is uniform (across all edges), that value is the only value in the vector. Else, there is one value for each edge in the lattice."""
    transport_rates::Vector{Pair{Int64, Vector{S}}}
    """A matrix, NxM, where N is the number of species with transportation and M the number of vertexes. Each value is the total rate at which that species leaves that vertex (e.g. for a species with constant diffusion rate D, in a vertex with n neighbours, this value is n*D)."""
    leaving_rates::Matrix{S}
    """An (enumerate'ed) iterator over all the edges of the lattice."""
    edges::Graphs.SimpleGraphs.SimpleEdgeIter{SimpleDiGraph{Int64}}
    """The edge parameters used to create the spatial ODEProblem. Currently unused, but will be needed to support changing these (e.g. due to events). Contain one vector for each edge parameter (length one if uniform, else one value for each edge)."""
    edge_ps::Vector{Vector{T}}
    
    function LatticeTransportODEf(ofunc::R, vert_ps::Vector{Vector{S}}, transport_rates::Vector{Pair{Int64, Vector{S}}}, edge_ps::Vector{Vector{T}}, lrs::LatticeReactionSystem) where {R,S,T}
        leaving_rates = zeros(length(transport_rates), lrs.num_verts)
        for (s_idx, trpair) in enumerate(transport_rates)
            rates = last(trpair)
            for (e_idx, e) in enumerate(edges(lrs.lattice))    
                # Updates the exit rate for species s_idx from vertex e.src
                leaving_rates[s_idx, e.src] += get_component_value(rates, e_idx)    
            end
        end
        work_vert_ps = zeros(lrs.num_verts)
        # 1 if ps are constant across the graph, 0 else.
        v_ps_idx_types = map(vp -> length(vp) == 1, vert_ps)
        eds = edges(lrs.lattice)                  
        new{R,S,T}(ofunc, lrs.num_verts, lrs.num_species, vert_ps, work_vert_ps, v_ps_idx_types, transport_rates, leaving_rates, eds, edge_ps)
    end
end

# Functor structure containing the information for the forcing function of a spatial ODE with spatial movement on a lattice.
struct LatticeTransportODEjac{Q,R,S,T}
    """The ODEFunction of the (non-spatial) reaction system which generated this function."""
    ofunc::Q
    """The number of vertices."""
    num_verts::Int64
    """The number of species."""
    num_species::Int64
    """The values of the parameters which values are tied to vertexes."""
    vert_ps::Vector{Vector{R}}
    """Temporary vector. For parameters which values are identical across the lattice, at some point these have to be converted of a length(num_verts) vector. To avoid re-allocation they are written to this vector."""
    work_vert_ps::Vector{R} 
    """For each parameter in vert_ps, it either have length num_verts or 1. To know whenever a parameter's value need expanding to the work_vert_ps array, its length needs checking. This check is done once, and the value stored to this array. This field (specifically) is an enumerate over that array."""
    v_ps_idx_types::Vector{Bool}
    """Whether the Jacobian is sparse or not."""
    sparse::Bool
    """The transport rates. Can be a dense matrix (for non-sparse) or as the "nzval" field if sparse."""
    jac_transport::S
    """The edge parameters used to create the spatial ODEProblem. Currently unused, but will be needed to support changing these (e.g. due to events). Contain one vector for each edge parameter (length one if uniform, else one value for each edge)."""
    edge_ps::Vector{Vector{T}}

    function LatticeTransportODEjac(ofunc::R, vert_ps::Vector{Vector{S}}, lrs::LatticeReactionSystem, jac_transport::Union{Nothing, SparseMatrixCSC{Float64, Int64}}, edge_ps::Vector{Vector{T}}, sparse::Bool) where {R,S,T}
        work_vert_ps = zeros(lrs.num_verts)
        v_ps_idx_types = map(vp -> length(vp) == 1, vert_ps)
        new{R,S,typeof(jac_transport),T}(ofunc, lrs.num_verts, lrs.num_species, vert_ps, work_vert_ps, v_ps_idx_types, sparse, jac_transport, edge_ps)
    end
end

### ODEProblem ###

# Creates an ODEProblem from a LatticeReactionSystem.
function DiffEqBase.ODEProblem(lrs::LatticeReactionSystem, u0_in, tspan,
                               p_in = DiffEqBase.NullParameters(), args...;
                               jac = false, sparse = false, 
                               name = nameof(lrs), include_zero_odes = true,
                               combinatoric_ratelaws = get_combinatoric_ratelaws(lrs.rs),
                               remove_conserved = false, checks = false, kwargs...)
    is_transport_system(lrs) || error("Currently lattice ODE simulations are only supported when all spatial reactions are TransportReactions.")
    
    # Converts potential symmaps to varmaps (parameter conversion is more involved since the vertex and edge parameters may be given in a tuple, or in a common vector).
    u0_in = symmap_to_varmap(lrs, u0_in)
    p_in = (p_in isa Tuple{<:Any,<:Any}) ? (symmap_to_varmap(lrs, p_in[1]),symmap_to_varmap(lrs, p_in[2])) : symmap_to_varmap(lrs, p_in)    

    # Converts u0 and p to their internal forms.
    # u0 is [spec 1 at vert 1, spec 2 at vert 1, ..., spec 1 at vert 2, ...].
    u0 = lattice_process_u0(u0_in, species(lrs), lrs.num_verts)                                   
    # Both vert_ps and edge_ps becomes vectors of vectors. Each have 1 element for each parameter. 
    # These elements are length 1 vectors (if the parameter is uniform), or length num_verts/nE, with unique values for each vertex/edge (for vert_ps/edge_ps, respectively).
    vert_ps, edge_ps = lattice_process_p(p_in, vertex_parameters(lrs), edge_parameters(lrs), lrs)   

    # Creates ODEProblem.
    ofun = build_odefunction(lrs, vert_ps, edge_ps, jac, sparse, name, include_zero_odes, combinatoric_ratelaws, remove_conserved, checks)
    return ODEProblem(ofun, u0, tspan, vert_ps, args...; kwargs...) 
end

# Builds an ODEFunction for a spatial ODEProblem.
function build_odefunction(lrs::LatticeReactionSystem, vert_ps::Vector{Vector{T}},
                           edge_ps::Vector{Vector{T}}, jac::Bool, sparse::Bool,
                           name, include_zero_odes, combinatoric_ratelaws, remove_conserved, checks) where {T}
    println()
    remove_conserved && error("Removal of conserved quantities is currently not supported for `LatticeReactionSystem`s")

    # Creates a map, taking (the index in species(lrs) each species (with transportation) to its transportation rate (uniform or one value for each edge).
    transport_rates = make_sidxs_to_transrate_map(vert_ps, edge_ps, lrs)    

    # Prepares the Jacobian and forcing functions (depending on jacobian and sparsity selection).
    if jac
        ofunc_dense = ODEFunction(convert(ODESystem, lrs.rs; name, combinatoric_ratelaws, include_zero_odes, checks); jac = true, sparse = false) # Always used for build_jac_prototype.
        ofunc_sparse = ODEFunction(convert(ODESystem, lrs.rs; name, combinatoric_ratelaws, include_zero_odes, checks); jac = true, sparse = true) # Always used for LatticeTransportODEjac.
        jac_vals = build_jac_prototype(ofunc_sparse.jac_prototype, transport_rates, lrs; set_nonzero = true)
        if sparse
            f = LatticeTransportODEf(ofunc_sparse, vert_ps, transport_rates, edge_ps, lrs)
            jac_vals = build_jac_prototype(ofunc_sparse.jac_prototype, transport_rates, lrs; set_nonzero = true)
            J = LatticeTransportODEjac(ofunc_dense, vert_ps, lrs, jac_vals, edge_ps, true)
            jac_prototype = jac_vals
        else
            f = LatticeTransportODEf(ofunc_dense, vert_ps, transport_rates, edge_ps, lrs)
            J = LatticeTransportODEjac(ofunc_dense, vert_ps, lrs, jac_vals, edge_ps, false)
            jac_prototype = nothing
        end
    else
        if sparse
            ofunc_sparse = ODEFunction(convert(ODESystem, lrs.rs; name, combinatoric_ratelaws, include_zero_odes, checks); jac = false, sparse = true)
            f = LatticeTransportODEf(ofunc_sparse, vert_ps, transport_rates, edge_ps, lrs)
            jac_prototype = build_jac_prototype(ofunc_sparse.jac_prototype, transport_rates, lrs; set_nonzero = false)
        else
            ofunc_dense = ODEFunction(convert(ODESystem, lrs.rs; name, combinatoric_ratelaws, include_zero_odes, checks); jac = false, sparse = false)
            f = LatticeTransportODEf(ofunc_dense, vert_ps, transport_rates, edge_ps, lrs)
            jac_prototype = nothing
        end
        J = nothing
    end

    return ODEFunction(f; jac = J, jac_prototype = jac_prototype)           
end

# Builds a jacobian prototype. If requested, populate it with the Jacobian's (constant) values as well.
function build_jac_prototype(ns_jac_prototype::SparseMatrixCSC{Float64, Int64}, trans_rates, lrs::LatticeReactionSystem;
                             set_nonzero = false)
    # Finds the indexes of the transport species, and the species with transport only (and no non-spatial dynamics).
    trans_species = first.(trans_rates)
    trans_only_species = filter(s_idx -> !Base.isstored(ns_jac_prototype, s_idx, s_idx), trans_species)

    # Finds the indexes of all terms in the non-spatial jacobian.
    ns_jac_prototype_idxs = findnz(ns_jac_prototype)
    ns_i_idxs = ns_jac_prototype_idxs[1]
    ns_j_idxs = ns_jac_prototype_idxs[2]

    # List the indexes of all non-zero Jacobian terms.
    non_spat_terms = [[get_index(vert, s_i, lrs.num_species), get_index(vert, s_j, lrs.num_species)] for vert in 1:(lrs.num_verts) for (s_i, s_j) in zip(ns_i_idxs,ns_j_idxs)]        # Indexes of elements due to non-spatial dynamics.
    trans_only_leaving_terms = [[get_index(e.src, s_idx, lrs.num_species), get_index(e.src, s_idx, lrs.num_species)] for e in edges(lrs.lattice) for s_idx in trans_only_species]     # Indexes due to terms for a species leaves its current vertex (but does not have non-spatial dynamics). If the non-spatial Jacobian is fully dense, these would already be accounted for. 
    trans_arriving_terms = [[get_index(e.src, s_idx, lrs.num_species), get_index(e.dst, s_idx, lrs.num_species)] for e in edges(lrs.lattice) for s_idx in trans_species]              # Indexes due to terms for species arriving into a new vertex.
    all_terms = [non_spat_terms; trans_only_leaving_terms; trans_arriving_terms]

    # Creates a jacobian prototype with 0 values in all positions).
    jac_prototype = sparse(first.(all_terms), last.(all_terms), fill(0.0, length(all_terms)))

    # Set element values.
    if set_nonzero
        for (s, rates) in trans_rates, (e_idx, e) in enumerate(edges(lrs.lattice))
            jac_prototype[get_index(e.src, s, lrs.num_species), get_index(e.src, s, lrs.num_species)] -= get_component_value(rates, e_idx)    # Term due to species leaving source vertex.
            jac_prototype[get_index(e.src, s, lrs.num_species), get_index(e.dst, s, lrs.num_species)] += get_component_value(rates, e_idx)    # Term due to species arriving to destination vertex.
        end
    end

    return jac_prototype
end

# Defines the forcing functor's effect on the (spatial) ODE system.
function (f_func::LatticeTransportODEf)(du, u, p, t)
    # Updates for non-spatial reactions.
    for vert_i in 1:(f_func.num_verts)
        # gets the indices of species at vertex i
        idxs = get_indexes(vert_i, f_func.num_species)
        
        # vector of vertex ps at vert_i
        vert_i_ps = view_vert_ps_vector!(f_func.work_vert_ps, p, vert_i, enumerate(f_func.v_ps_idx_types))
        
        # evaluate reaction contributions to du at vert_i
        f_func.ofunc((@view du[idxs]),  (@view u[idxs]), vert_i_ps, t)
    end

    # s_idx is species index among transport species, s is index among all species
    # rates are the species' transport rates
    for (s_idx, (s, rates)) in enumerate(f_func.transport_rates)  
        # Rate for leaving vert_i
        for vert_i in 1:(f_func.num_verts)                                 
            idx = get_index(vert_i, s, f_func.num_species)
            du[idx] -= f_func.leaving_rates[s_idx, vert_i] * u[idx]
        end
        # Add rates for entering a given vertex via an incoming edge
        for (e_idx, e) in enumerate(f_func.edges)
            idx_dst = get_index(e.dst, s, f_func.num_species)
            idx_src = get_index(e.src, s, f_func.num_species)
            du[idx_dst] += get_component_value(rates, e_idx) * u[idx_src] 
        end
    end
end

# Defines the jacobian functor's effect on the (spatial) ODE system.
function (jac_func::LatticeTransportODEjac)(J, u, p, t)
    J .= 0.0 

    # Update the Jacobian from reaction terms
    for vert_i in 1:(jac_func.num_verts)
        idxs = get_indexes(vert_i, jac_func.num_species)
        vert_ps = view_vert_ps_vector!(jac_func.work_vert_ps, p, vert_i, enumerate(jac_func.v_ps_idx_types))
        jac_func.ofunc.jac((@view J[idxs, idxs]), (@view u[idxs]), vert_ps, t)
    end

    # Updates for the spatial reactions (adds the Jacobian values from the diffusion reactions).
    J .+= jac_func.jac_transport
end