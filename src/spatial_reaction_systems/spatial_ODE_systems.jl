### Spatial ODE Functor Structures ###

# Functor with information for the forcing function of a spatial ODE with spatial movement on a lattice.
struct LatticeTransportODEf{S,T}
    """The ODEFunction of the (non-spatial) ReactionSystem that generated this LatticeTransportODEf instance."""
    ofunc::S
    """The number of vertices."""
    num_verts::Int64
    """The number of species."""
    num_species::Int64
    """The indexes of the vertex parameters in the parameter vector (`parameters(lrs)`)."""
    vert_p_idxs::Vector{Int64}
    """The indexes of the edge parameters in the parameter vector (`parameters(lrs)`)."""
    edge_p_idxs::Vector{Int64}
    """
    The non-spatial `ReactionSystem` which was used to create the `LatticeReactionSystem` contain
    a set of parameters (either identical to, or a sub set of, `parameters(lrs)`). This vector
    contain the indexes of the non-spatial system's parameters in `parameters(lrs)`. These are
    required to manage the non-spatial ODEFunction in the spatial call.
    """
    nonspatial_rs_p_idxs::Vector{Int64}
    """The values of the parameters that are tied to vertices."""
    vert_ps::Vector{Vector{T}}
    """
    Vector for storing temporary values. Repeatedly during simulations, we need to retrieve the 
    parameter values in a certain vertex. However, since most parameters (likely) are uniform
    (and hence only have 1 value stored), we need to create a new vector each time we need to retrieve
    the parameter values in a new vertex. To avoid relocating these values repeatedly, we write them
    to this vector.
    """
    work_ps::Vector{T}    
    """
    For each parameter in vert_ps, its value is a vector with a length of either num_verts or 1. 
    To know whenever a parameter's value needs expanding to the work_ps array, its length needs checking. 
    This check is done once, and the value is stored in this array. True means a uniform value.
    """
    v_ps_idx_types::Vector{Bool}
    """
    A vector that stores, for each species with transportation, its transportation rate(s). 
    Each entry is a pair from (the index of) the transported species (in the `species(lrs)` vector)
    to its transportation rate (each species only has a single transportation rate, the sum of all
    its transportation reactions' rates). If the transportation rate is uniform across all edges, 
    stores a single value (in a size (1,1) sparse matrix). Otherwise, stores these in a sparse matrix 
    where value (i,j) is the species transportation rate from vertex i to vertex j.
    """
    transport_rates::Vector{Pair{Int64,SparseMatrixCSC{T, Int64}}}
    """
    For each transport rate in transport_rates, its value is a (sparse) matrix with a size of either 
    (num_verts,num_verts) or (1,1). In the second case, the transportation rate is uniform across 
    all edges. To avoid having to check which case holds for each transportation rate, we store the
    corresponding case in this value. `true` means that a species has a uniform transportation rate.
    """
    t_rate_idx_types::Vector{Bool}
    """
    A matrix, NxM, where N is the number of species with transportation and M is the number of vertices. 
    Each value is the total rate at which that species leaves that vertex 
    (e.g. for a species with constant diffusion rate D, in a vertex with n neighbours, this value is n*D).
    """
    leaving_rates::Matrix{T}
    """An iterator over all the edges of the lattice."""
    edge_iterator::Vector{Pair{Int64, Int64}}
    
    function LatticeTransportODEf(ofunc::S, vert_ps::Vector{Pair{BasicSymbolic{Real},Vector{T}}}, 
                                  transport_rates::Vector{Pair{Int64, SparseMatrixCSC{T, Int64}}}, 
                                  lrs::LatticeReactionSystem) where {S,T}        
        # Records which parameters and rates are uniform and which are not.
        v_ps_idx_types = map(vp -> length(vp[2]) == 1, vert_ps)
        t_rate_idx_types = map(tr -> size(tr[2]) == (1,1), transport_rates)

        # Computes the indexes of various parameters in in the `parameters(lrs)` vector.
        vert_p_idxs = subset_indexes_of(vertex_parameters(lrs), parameters(lrs))
        edge_p_idxs = subset_indexes_of(edge_parameters(lrs), parameters(lrs))
        nonspatial_rs_p_idxs = subset_indexes_of(parameters(reactionsystem(lrs)), parameters(lrs))

        # Computes the indexes of the vertex parameters in the vector of parameters.
        # Input `vert_ps` is a vector map taking each parameter symbolic to its value (potentially a 
        # vector). This vector is already sorted according to the order of the parameters. Here, we extract 
        # its values only and put them into `vert_ps`.
        vert_ps = [vp[2] for vp in vert_ps]

        # Computes the leaving rate matrix.
        leaving_rates = zeros(length(transport_rates), num_verts(lrs))
        for (s_idx, tr_pair) in enumerate(transport_rates)
            for e in Catalyst.edge_iterator(lrs)
                # Updates the exit rate for species s_idx from vertex e.src.
                leaving_rates[s_idx, e[1]] += get_transport_rate(tr_pair[2], e, t_rate_idx_types[s_idx]) 
            end
        end

        # Declares `work_ps` (used as storage during computation) and the edge iterator.
        work_ps = zeros(length(parameters(lrs)))
        edge_iterator = Catalyst.edge_iterator(lrs) 
        new{S,T}(ofunc, num_verts(lrs), num_species(lrs), vert_p_idxs, edge_p_idxs, 
                 nonspatial_rs_p_idxs, vert_ps, work_ps, v_ps_idx_types, transport_rates, 
                 t_rate_idx_types, leaving_rates, edge_iterator)
    end
end

# Functor with information for the Jacobian function of a spatial ODE with spatial movement on a lattice.
struct LatticeTransportODEjac{R,S,T}
    """The ODEFunction of the (non-spatial) ReactionSystem that generated this LatticeTransportODEf instance."""
    ofunc::R
    """The number of vertices."""
    num_verts::Int64
    """The number of species."""
    num_species::Int64
    """The indexes of the vertex parameters in the parameter vector (`parameters(lrs)`)."""
    vert_p_idxs::Vector{Int64}
    """The indexes of the edge parameters in the parameter vector (`parameters(lrs)`)."""
    edge_p_idxs::Vector{Int64}
    """
    The non-spatial `ReactionSystem` which was used to create the `LatticeReactionSystem` contain
    a set of parameters (either identical to, or a sub set of, `parameters(lrs)`). This vector
    contain the indexes of the non-spatial system's parameters in `parameters(lrs)`. These are
    required to manage the non-spatial ODEFunction in the spatial call.
    """
    nonspatial_rs_p_idxs::Vector{Int64}
    """The values of the parameters that are tied to vertices."""
    vert_ps::Vector{Vector{S}}
    """
    Vector for storing temporary values. Repeatedly during simulations, we need to retrieve the 
    parameter values in a certain vertex. However, since most parameters (likely) are uniform
    (and hence only have 1 value stored), we need to create a new vector each time we need to retrieve
    the parameter values in a new vertex. To avoid relocating these values repeatedly, we write them
    to this vector.
    """
    work_ps::Vector{S}    
    """
    For each parameter in vert_ps, its value is a vector with a length of either num_verts or 1. 
    To know whenever a parameter's value needs expanding to the work_ps array, its length needs checking. 
    This check is done once, and the value is stored in this array. True means a uniform value.
    """
    v_ps_idx_types::Vector{Bool}
    """Whether the Jacobian is sparse or not."""
    sparse::Bool
    """The transport rates. This is a dense or sparse matrix (depending on what type of Jacobian is used)."""
    jac_transport::T

    function LatticeTransportODEjac(ofunc::R, vert_ps::Vector{Pair{BasicSymbolic{Real},Vector{S}}}, 
                                    jac_transport::Union{Nothing, SparseMatrixCSC{Float64, Int64}}, 
                                    lrs::LatticeReactionSystem, sparse::Bool) where {R,S}

        # Computes the indexes of various parameters in in the `parameters(lrs)` vector.
        vert_p_idxs = subset_indexes_of(vertex_parameters(lrs), parameters(lrs))
        edge_p_idxs = subset_indexes_of(edge_parameters(lrs), parameters(lrs))
        nonspatial_rs_p_idxs = subset_indexes_of(parameters(reactionsystem(lrs)), parameters(lrs))

        # Input `vert_ps` is a vector map taking each parameter symbolic to its value (potentially a 
        # vector). This vector is already sorted according to the order of the parameters. Here, we extract 
        # its values only and put them into `vert_ps`.
        vert_ps = [vp[2] for vp in vert_ps]

        work_ps = zeros(length(parameters(lrs)))
        v_ps_idx_types = map(vp -> length(vp) == 1, vert_ps)
        new{R,S,typeof(jac_transport)}(ofunc, num_verts(lrs), num_species(lrs) , vert_p_idxs, 
                                       edge_p_idxs, nonspatial_rs_p_idxs, vert_ps, 
                                       work_ps, v_ps_idx_types, sparse, jac_transport)
    end
end

# For each symbolic in syms1, returns a vector with their indexes in syms2.
function subset_indexes_of(syms1, syms2)
    [findfirst(isequal(sym1, sym2) for sym2 in syms2) for sym1 in syms1]
end

### ODEProblem ###

# Creates an ODEProblem from a LatticeReactionSystem.
function DiffEqBase.ODEProblem(lrs::LatticeReactionSystem, u0_in, tspan,
                               p_in = DiffEqBase.NullParameters(), args...;
                               jac = false, sparse = false, 
                               name = nameof(lrs), include_zero_odes = true,
                               combinatoric_ratelaws = get_combinatoric_ratelaws(reactionsystem(lrs)),
                               remove_conserved = false, checks = false, kwargs...)
    if !is_transport_system(lrs)
        error("Currently lattice ODE simulations are only supported when all spatial reactions are TransportReactions.")
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
    # edge_ps becomes a sparse matrix. Here, index (i,j) is a parameter's value in the edge from vertex i to vertex j.
    # Uniform vertex/edge parameters store only a single value (a length 1 vector, or size 1x1 sparse matrix).
    # In the `ODEProblem` vert_ps and edge_ps are merged (but for building the ODEFunction, they are separate).
    vert_ps, edge_ps = lattice_process_p(p_in, vertex_parameters(lrs), edge_parameters(lrs), lrs)  

    # Creates the ODEFunction.
    ofun = build_odefunction(lrs, vert_ps, edge_ps, jac, sparse, name, include_zero_odes, 
                             combinatoric_ratelaws, remove_conserved, checks)

    # Combines `vert_ps` and `edge_ps` to a single vector with values only (not a map). Creates ODEProblem.
    pval_dict = Dict([vert_ps; edge_ps])
    ps = [pval_dict[p] for p in parameters(lrs)]
    return ODEProblem(ofun, u0, tspan, ps, args...; kwargs...) 
end

# Builds an ODEFunction for a spatial ODEProblem.
function build_odefunction(lrs::LatticeReactionSystem, vert_ps::Vector{Pair{BasicSymbolic{Real},Vector{T}}},
                           edge_ps::Vector{Pair{BasicSymbolic{Real},SparseMatrixCSC{T, Int64}}}, 
                           jac::Bool, sparse::Bool, name, include_zero_odes, combinatoric_ratelaws, 
                           remove_conserved, checks) where {T}
    if remove_conserved 
        error("Removal of conserved quantities is currently not supported for `LatticeReactionSystem`s")
    end

    # Creates a map, taking (the index in species(lrs) of) each species (with transportation)
    # to its transportation rate (uniform or one value for each edge). The rates are sparse matrices.
    transport_rates = make_sidxs_to_transrate_map(vert_ps, edge_ps, lrs)

    # Prepares the Jacobian and forcing functions (depending on jacobian and sparsity selection).
    osys = complete(convert(ODESystem, reactionsystem(lrs); name, combinatoric_ratelaws, include_zero_odes, checks))
    if jac
        # `build_jac_prototype` currently assumes a sparse (non-spatial) Jacobian. Hence compute this.
        # `LatticeTransportODEjac` currently assumes a dense (non-spatial) Jacobian. Hence compute this.
        # Long term we could write separate versions of these functions for generic input.
        ofunc_dense = ODEFunction(osys; jac = true, sparse = false)
        ofunc_sparse = ODEFunction(osys; jac = true, sparse = true)
        jac_transport = build_jac_prototype(ofunc_sparse.jac_prototype, transport_rates, lrs; set_nonzero = true)
        if sparse
            f = LatticeTransportODEf(ofunc_sparse, vert_ps, transport_rates, lrs)
            jac_transport = build_jac_prototype(ofunc_sparse.jac_prototype, transport_rates, lrs; set_nonzero = true)
            J = LatticeTransportODEjac(ofunc_dense, vert_ps, jac_transport, lrs, true)
            jac_prototype = jac_transport
        else
            f = LatticeTransportODEf(ofunc_dense, vert_ps, transport_rates, lrs)
            J = LatticeTransportODEjac(ofunc_dense, vert_ps, jac_transport, lrs, false)
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

# Builds a jacobian prototype.
# If requested, populate it with the constant values of the Jacobian's transportation part.
function build_jac_prototype(ns_jac_prototype::SparseMatrixCSC{Float64, Int64}, 
                             transport_rates::Vector{Pair{Int64,SparseMatrixCSC{T, Int64}}}, 
                             lrs::LatticeReactionSystem; set_nonzero = false) where {T}
    # Finds the indices of  both the transport species,
    # and the species with transport only (that is, with no non-spatial dynamics but with spatial dynamics).
    trans_species = [tr[1] for tr in transport_rates]
    trans_only_species = filter(s_idx -> !Base.isstored(ns_jac_prototype, s_idx, s_idx), trans_species)

    # Finds the indices of all terms in the non-spatial jacobian.
    ns_jac_prototype_idxs = findnz(ns_jac_prototype)
    ns_i_idxs = ns_jac_prototype_idxs[1]
    ns_j_idxs = ns_jac_prototype_idxs[2]

    # Prepares vectors to store i and j indices of Jacobian entries.
    idx = 1
    num_entries = num_verts(lrs) * length(ns_i_idxs) + 
                  num_edges(lrs) * (length(trans_only_species) + length(trans_species))
    i_idxs = Vector{Int}(undef, num_entries)
    j_idxs = Vector{Int}(undef, num_entries)

    # Indices of elements caused by non-spatial dynamics.
    for vert in 1:num_verts(lrs)
        for n in 1:length(ns_i_idxs)
            i_idxs[idx] = get_index(vert, ns_i_idxs[n], num_species(lrs))
            j_idxs[idx] = get_index(vert, ns_j_idxs[n], num_species(lrs))
            idx += 1
        end
    end

    # Indices of elements caused by spatial dynamics.
    for e in edge_iterator(lrs)
        # Indexes due to terms for a species leaving its source vertex (but does not have
        # non-spatial dynamics). If the non-spatial Jacobian is fully dense, these would already
        # be accounted for.
        for s_idx in trans_only_species
            i_idxs[idx] = get_index(e[1], s_idx, num_species(lrs))
            j_idxs[idx] = i_idxs[idx]
            idx += 1
        end
        # Indexes due to terms for species arriving into a destination vertex.
        for s_idx in trans_species
            i_idxs[idx] = get_index(e[1], s_idx, num_species(lrs))
            j_idxs[idx] = get_index(e[2], s_idx, num_species(lrs))
            idx += 1
        end
    end

    # Create a sparse Jacobian prototype with 0-valued entries.
    jac_prototype = sparse(i_idxs, j_idxs, zeros(num_entries))

    # Set element values.
    if set_nonzero
        for (s, rates) in transport_rates, e in edge_iterator(lrs)
            idx_src = get_index(e[1], s, num_species(lrs))
            idx_dst = get_index(e[2], s, num_species(lrs))
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
        # Gets the indices of all the species at vertex i.
        idxs = get_indexes(vert_i, f_func.num_species)
        
        # Updates the work vector to contain the vertex parameter values for vertex vert_i.
        update_work_vert_ps!(f_func, p, vert_i)
        
        # Evaluate reaction contributions to du at vert_i.
        f_func.ofunc((@view du[idxs]),  (@view u[idxs]), nonspatial_ps(f_func), t)
    end

    # s_idx is the species index among transport species, s is the index among all species.
    # rates are the species' transport rates.
    for (s_idx, (s, rates)) in enumerate(f_func.transport_rates)  
        # Rate for leaving source vertex vert_i. 
        for vert_i in 1:(f_func.num_verts)  
            idx_src = get_index(vert_i, s, f_func.num_species)                       
            du[idx_src] -= f_func.leaving_rates[s_idx, vert_i] * u[idx_src]
        end
        # Add rates for entering a destination vertex via an incoming edge.
        for e in f_func.edge_iterator   
            idx_src = get_index(e[1], s, f_func.num_species)
            idx_dst = get_index(e[2], s, f_func.num_species)     
            du[idx_dst] += get_transport_rate(s_idx, f_func, e) * u[idx_src]
        end
    end
end

# Defines the Jacobian functor's effect on the (spatial) ODE system.
function (jac_func::LatticeTransportODEjac)(J, u, p, t)
    J .= 0.0 

    # Update the Jacobian from non-spatial reaction terms.
    for vert_i in 1:(jac_func.num_verts)
        idxs = get_indexes(vert_i, jac_func.num_species)
        update_work_vert_ps!(jac_func, p, vert_i)
        jac_func.ofunc.jac((@view J[idxs, idxs]), (@view u[idxs]), nonspatial_ps(jac_func), t)
    end

    # Updates for the spatial reactions (adds the Jacobian values from the transportation reactions).
    J .+= jac_func.jac_transport
end
