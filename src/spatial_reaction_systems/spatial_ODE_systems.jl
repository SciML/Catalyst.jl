### Spatial ODE Functor Structures ###

# Functor with information for the forcing function of a spatial ODE with spatial movement on a lattice.
struct LatticeTransportODEf{S,T}
    """The ODEFunction of the (non-spatial) reaction system which generated this function."""
    ofunc::S
    """The number of vertices."""
    num_verts::Int64
    """The number of species."""
    num_species::Int64
    """The values of the parameters which values are tied to vertexes."""
    vert_ps::Vector{Vector{T}}
    """
    Temporary vector. For parameters which values are identical across the lattice, 
    at some point these have to be converted of a length num_verts vector. 
    To avoid re-allocation they are written to this vector.
    """
    work_vert_ps::Vector{T}    
    """
    For each parameter in vert_ps, its value is a vector with length either num_verts or 1. 
    To know whenever a parameter's value need expanding to the work_vert_ps array, its length needs checking. 
    This check is done once, and the value stored to this array. True means a uniform value.
    """
    v_ps_idx_types::Vector{Bool}
    """
    A vector of sparse, with a value for each species with transportation. 
    The first value is the species index (in the species(::ReactionSystem) vector), 
    and the second is a vector with its transport rate values. 
    If the transport rate is uniform (across all edges), that value is the only value in the vector. 
    Else, there is one value for each edge in the lattice.
    """
    transport_rates::Vector{Pair{Int64,SparseMatrixCSC{T, Int64}}}
    """
    For each transport rate in transport_rates, its value is a (sparse) matrix with size either 
    (num_verts,num_verts) or (1,1). In the second case, that transportation rate is uniform across all edges. 
    To know how to access transport rate's value (without checking sizes), we can use this vector directly.
    True means a uniform value.
    """
    t_rate_idx_types::Vector{Bool}
    """
    A matrix, NxM, where N is the number of species with transportation and M the number of vertexes. 
    Each value is the total rate at which that species leaves that vertex 
    (e.g. for a species with constant diffusion rate D, in a vertex with n neighbours, this value is n*D).
    """
    leaving_rates::Matrix{T}
    """An iterator over all the edges of the lattice."""
    edge_iterator::Vector{Pair{Int64, Int64}}
    
    function LatticeTransportODEf(ofunc::S, vert_ps::Vector{Pair{BasicSymbolic{Real},Vector{T}}}, transport_rates::Vector{Pair{Int64, SparseMatrixCSC{T, Int64}}}, 
                                  lrs::LatticeReactionSystem) where {S,T}        
        # Records. which parameters and rates are uniform and not.
        v_ps_idx_types = map(vp -> length(vp[2]) == 1, vert_ps)
        t_rate_idx_types = map(tr -> size(tr[2]) == (1,1), transport_rates)
        
        # Input `vert_ps` is a vector map taking each parameter symbolic to its value (potentially a vector).
        # This vector is sorted according to the parameters order. Here, we simply extract its values only.
        vert_ps = [vp[2] for vp in vert_ps]

        # Computes the leaving rates.
        leaving_rates = zeros(length(transport_rates), lrs.num_verts)
        for (s_idx, trpair) in enumerate(transport_rates)
            t_rate = trpair[2]
            for e in lrs.edge_iterator 
                # Updates the exit rate for species s_idx from vertex e.src
                leaving_rates[s_idx, e[1]] += get_transport_rate(t_rate, e, t_rate_idx_types[s_idx]) 
            end
        end

        # Declares `work_vert_ps` (used as storage during computation) and an iterator over the edges.
        work_vert_ps = zeros(length(vert_ps))
        edge_iterator = lrs.edge_iterator               
        new{S,T}(ofunc, lrs.num_verts, lrs.num_species, vert_ps, work_vert_ps, 
                               v_ps_idx_types, transport_rates, t_rate_idx_types, leaving_rates, edge_iterator)
    end
end

# Functor with information for the Jacobian function of a spatial ODE with spatial movement on a lattice.
struct LatticeTransportODEjac{R,S,T}
    """The ODEFunction of the (non-spatial) reaction system which generated this function."""
    ofunc::R
    """The number of vertices."""
    num_verts::Int64
    """The number of species."""
    num_species::Int64
    """The values of the parameters which values are tied to vertexes."""
    vert_ps::Vector{Vector{S}}
    """
    Temporary vector. For parameters which values are identical across the lattice, 
    at some point these have to be converted of a length(num_verts) vector. 
    To avoid re-allocation they are written to this vector.
    """
    work_vert_ps::Vector{S} 
    """
    For each parameter in vert_ps, it either have length num_verts or 1. 
    To know whenever a parameter's value need expanding to the work_vert_ps array, 
    its length needs checking. This check is done once, and the value stored to this array. 
    This field (specifically) is an enumerate over that array.
    """
    v_ps_idx_types::Vector{Bool}
    """Whether the Jacobian is sparse or not."""
    sparse::Bool
    """The transport rates. This is a dense or sparse matrix (depending on what type of Jacobian is used)."""
    jac_transport::T

    function LatticeTransportODEjac(ofunc::R, vert_ps::Vector{Pair{BasicSymbolic{Real},Vector{S}}}, lrs::LatticeReactionSystem, 
                                    jac_transport::Union{Nothing, SparseMatrixCSC{Float64, Int64}}, 
                                    sparse::Bool) where {R,S}
        # Input `vert_ps` is a vector map taking each parameter symbolic to its value (potentially a vector).
        # This vector is sorted according to the parameters order. Here, we simply extract its values only.
        vert_ps = [vp[2] for vp in vert_ps]

        work_vert_ps = zeros(lrs.num_verts)
        v_ps_idx_types = map(vp -> length(vp) == 1, vert_ps)
        new{R,S,typeof(jac_transport)}(ofunc, lrs.num_verts, lrs.num_species, vert_ps, 
                                        work_vert_ps, v_ps_idx_types, sparse, jac_transport)
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
    if !is_transport_system(lrs)
        error("Currently lattice ODE simulations are only supported when all spatial reactions are TransportReactions.")
    end
    
    # Converts potential symmaps to varmaps
    # Vertex and edge parameters may be given in a tuple, or in a common vector, making parameter case complicated.
    u0_in = symmap_to_varmap(lrs, u0_in)
    p_in = symmap_to_varmap(lrs, p_in)    

    # Converts u0 and p to their internal forms.
    # u0 is simply a vector with all the species initial condition values across all vertexes.
    # u0 is [spec 1 at vert 1, spec 2 at vert 1, ..., spec 1 at vert 2, ...].
    u0 = lattice_process_u0(u0_in, species(lrs), lrs)                                       
    # vert_ps and `edge_ps` are vector maps, taking each parameter's Symbolics to its value(s).
    # vert_ps values are vectors. Here, index (i) is a parameters value in vertex i.
    # edge_ps becomes sparse matrix. Here, index (i,j) is a parameters value in the edge from vertex i to vertex j.
    # Uniform vertex/edge parameters stores only a single value (in a length 1 vector, or size 1x1 sparse matrix).
    # This is the parameters single value.
    # In the `ODEProblem` vert_ps and edge_ps are merged (but for building the ODEFunction, they are separate).
    vert_ps, edge_ps = lattice_process_p(p_in, vertex_parameters(lrs), edge_parameters(lrs), lrs)  

    # Creates ODEProblem.
    ofun = build_odefunction(lrs, vert_ps, edge_ps, jac, sparse, name, include_zero_odes, 
                                combinatoric_ratelaws, remove_conserved, checks)

    # Combines `vert_ps` and `edge_ps` to a single vector with values only (not a map). Creates ODEProblem.
    ps = [p[2] for p in [vert_ps; edge_ps]]
    return ODEProblem(ofun, u0, tspan, ps, args...; kwargs...) 
end

# Builds an ODEFunction for a spatial ODEProblem.
function build_odefunction(lrs::LatticeReactionSystem, vert_ps::Vector{Pair{A,Vector{T}}},
                           edge_ps::Vector{Pair{B,SparseMatrixCSC{T, Int64}}}, jac::Bool, sparse::Bool,
                           name, include_zero_odes, combinatoric_ratelaws, remove_conserved, checks) where {A,B,T}
    if remove_conserved 
        error("Removal of conserved quantities is currently not supported for `LatticeReactionSystem`s")
    end

    # Creates a map, taking (the index in species(lrs) each species (with transportation)
    # to its transportation rate (uniform or one value for each edge).
    transport_rates = make_sidxs_to_transrate_map(vert_ps, edge_ps, lrs)

    # Prepares the Jacobian and forcing functions (depending on jacobian and sparsity selection).
    osys = complete(convert(ODESystem, lrs.rs; name, combinatoric_ratelaws, include_zero_odes, checks))
    if jac
        # `build_jac_prototype` currently assumes a sparse (non-spatial) Jacobian. Hence compute this.
        # `LatticeTransportODEjac` currently assumes a dense (non-spatial) Jacobian. Hence compute this.
        # Long term we could write separate version of these functions for generic input.
        ofunc_dense = ODEFunction(osys; jac = true, sparse = false)
        ofunc_sparse = ODEFunction(osys; jac = true, sparse = true)
        jac_vals = build_jac_prototype(ofunc_sparse.jac_prototype, transport_rates, lrs; set_nonzero = true)
        if sparse
            f = LatticeTransportODEf(ofunc_sparse, vert_ps, transport_rates, lrs)
            jac_vals = build_jac_prototype(ofunc_sparse.jac_prototype, transport_rates, lrs; set_nonzero = true)
            J = LatticeTransportODEjac(ofunc_dense, vert_ps, lrs, jac_vals, true)
            jac_prototype = jac_vals
        else
            f = LatticeTransportODEf(ofunc_dense, vert_ps, transport_rates, lrs)
            J = LatticeTransportODEjac(ofunc_dense, vert_ps, lrs, jac_vals, false)
            jac_prototype = nothing
        end
    else
        if sparse
            ofunc_sparse = ODEFunction(osys; jac = false, sparse = true)
            f = LatticeTransportODEf(ofunc_sparse, vert_ps, transport_rates, lrs)
            jac_prototype = build_jac_prototype(ofunc_sparse.jac_prototype, transport_rates, lrs; set_nonzero = false)
        else
            ofunc_dense = ODEFunction(osys; jac = false, sparse = false)
            f = LatticeTransportODEf(ofunc_dense, vert_ps, transport_rates, lrs)
            jac_prototype = nothing
        end
        J = nothing
    end

    return ODEFunction(f; jac = J, jac_prototype = jac_prototype)           
end

# Builds a jacobian prototype. If requested, populate it with the Jacobian's (constant) values as well.
function build_jac_prototype(ns_jac_prototype::SparseMatrixCSC{Float64, Int64}, transport_rates::Vector{Pair{Int64,SparseMatrixCSC{T, Int64}}}, 
                                lrs::LatticeReactionSystem; set_nonzero = false) where T
    # Finds the indexes of the transport species, and the species with transport only (and no non-spatial dynamics).
    trans_species = [tr[1] for tr in transport_rates]
    trans_only_species = filter(s_idx -> !Base.isstored(ns_jac_prototype, s_idx, s_idx), trans_species)

    # Finds the indexes of all terms in the non-spatial jacobian.
    ns_jac_prototype_idxs = findnz(ns_jac_prototype)
    ns_i_idxs = ns_jac_prototype_idxs[1]
    ns_j_idxs = ns_jac_prototype_idxs[2]

    # Prepares vectors to store i and j indexes of Jacobian entries.
    idx = 1
    num_entries = lrs.num_verts * length(ns_i_idxs) + 
                  lrs.num_edges * (length(trans_only_species) + length(trans_species))
    i_idxs = Vector{Int}(undef, num_entries)
    j_idxs = Vector{Int}(undef, num_entries)

    # Indexes of elements due to non-spatial dynamics.
    for vert in 1:lrs.num_verts
        for n in 1:length(ns_i_idxs)
            i_idxs[idx] = get_index(vert, ns_i_idxs[n], lrs.num_species)
            j_idxs[idx] = get_index(vert, ns_j_idxs[n], lrs.num_species)
            idx += 1
        end
    end

    # Indexes of elements due to spatial dynamics.
    for e in lrs.edge_iterator
        # Indexes due to terms for a species leaves its current vertex (but does not have
        # non-spatial dynamics). If the non-spatial Jacobian is fully dense, these would already
        # be accounted for.
        for s_idx in trans_only_species
            i_idxs[idx] = get_index(e[1], s_idx, lrs.num_species)
            j_idxs[idx] = i_idxs[idx]
            idx += 1
        end
        # Indexes due to terms for species arriving into a new vertex.
        for s_idx in trans_species
            i_idxs[idx] = get_index(e[1], s_idx, lrs.num_species)
            j_idxs[idx] = get_index(e[2], s_idx, lrs.num_species)
            idx += 1
        end
    end

    # Create sparse jacobian prototype with 0-valued entries.
    jac_prototype = sparse(i_idxs, j_idxs, zeros(num_entries))

    # Set element values.
    if set_nonzero
        for (s, rates) in transport_rates, e in lrs.edge_iterator 
            idx_src = get_index(e[1], s, lrs.num_species)
            idx_dst = get_index(e[2], s, lrs.num_species)
            val = get_transport_rate(rates, e, size(rates)==(1,1))

            # Term due to species leaving source vertex.
            jac_prototype[idx_src, idx_src] -= val
            
            # Term due to species arriving to destination vertex.   
            jac_prototype[idx_src, idx_dst] += val
        end
    end

    return jac_prototype
end

# Defines the forcing functor's effect on the (spatial) ODE system.
function (f_func::LatticeTransportODEf)(du, u, p, t)
    # Updates for non-spatial reactions.
    for vert_i in 1:(f_func.num_verts)
        # Gets the indices of species at vertex i.
        idxs = get_indexes(vert_i, f_func.num_species)
        
        # Updates the vector which contains the vertex parameter values for vertex vert_i.
        update_work_vert_ps!(f_func, p, vert_i)
        
        # Evaluate reaction contributions to du at vert_i.
        f_func.ofunc((@view du[idxs]),  (@view u[idxs]), f_func.work_vert_ps, t)
    end

    # s_idx is species index among transport species, s is index among all species.
    # rates are the species' transport rates.
    for (s_idx, (s, rates)) in enumerate(f_func.transport_rates)  
        # Rate for leaving vert_i
        for vert_i in 1:(f_func.num_verts)                                 
            idx_src = get_index(vert_i, s, f_func.num_species)
            du[idx_src] -= f_func.leaving_rates[s_idx, vert_i] * u[idx_src]
        end
        # Add rates for entering a given vertex via an incoming edge.
        for e in f_func.edge_iterator            
            idx_src = get_index(e[1], s, f_func.num_species)
            idx_dst = get_index(e[2], s, f_func.num_species)
            du[idx_dst] += get_transport_rate(s_idx, f_func, e) * u[idx_src]
        end
    end
end

# Defines the jacobian functor's effect on the (spatial) ODE system.
function (jac_func::LatticeTransportODEjac)(J, u, p, t)
    J .= 0.0 

    # Update the Jacobian from reaction terms.
    for vert_i in 1:(jac_func.num_verts)
        idxs = get_indexes(vert_i, jac_func.num_species)
        update_work_vert_ps!(jac_func, p, vert_i)
        jac_func.ofunc.jac((@view J[idxs, idxs]), (@view u[idxs]), jac_func.work_vert_ps, t)
    end

    # Updates for the spatial reactions (adds the Jacobian values from the diffusion reactions).
    J .+= jac_func.jac_transport
end
