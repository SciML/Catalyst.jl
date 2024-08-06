### Spatial ODE Functor Structure ###

# Functor with information about a spatial Lattice Reaction ODEs forcing and Jacobian functions.
# Also used as ODE Function input to corresponding `ODEProblem`.
struct LatticeTransportODEFunction{P, Q, R, S, T}
    """
    The ODEFunction of the (non-spatial) ReactionSystem that generated this
    LatticeTransportODEFunction instance.
    """
    ofunc::P
    """The lattice's number of vertices."""
    num_verts::Int64
    """The system's number of species."""
    num_species::Int64
    """
    Stores an index for each heterogeneous vertex parameter (i.e. vertex parameter which value is 
    not identical across the lattice). Each index corresponds to its position in the full parameter
    vector (`parameters(lrs)`).
    """
    heterogeneous_vert_p_idxs::Vector{Int64}
    """
    The MTKParameters structure which corresponds to the non-spatial `ReactionSystem`. During 
    simulations, as we loop through each vertex, this is updated to correspond to the vertex 
    parameters of that specific vertex.
    """
    mtk_ps::Q
    """
    Stores a SymbolicIndexingInterface `setp` function for each heterogeneous vertex parameter (i.e. 
    vertex parameter whose value is not identical across the lattice). The `setp` function at index
    i of `p_setters` corresponds to the parameter in index i of `heterogeneous_vert_p_idxs`.
    """
    p_setters::R
    """
    A vector that stores, for each species with transportation, its transportation rate(s). 
    Each entry is a pair from (the index of) the transported species (in the `species(lrs)` vector)
    to its transportation rate (each species only has a single transportation rate, the sum of all
    its transportation reactions' rates). If the transportation rate is uniform across all edges, 
    stores a single value (in a size (1,1) sparse matrix). Otherwise, stores these in a sparse
    matrix  where value (i,j) is the species transportation rate from vertex i to vertex j.
    """
    transport_rates::Vector{Pair{Int64, SparseMatrixCSC{S, Int64}}}
    """
    For each transport rate in transport_rates, its value is a (sparse) matrix with a size of either 
    (num_verts,num_verts) or (1,1). In the second case, the transportation rate is uniform across 
    all edges. To avoid having to check which case holds for each transportation rate, we store the
    corresponding case in this value. `true` means that a species has a uniform transportation rate.
    """
    t_rate_idx_types::Vector{Bool}
    """
    A matrix, NxM, where N is the number of species with transportation and M is the number of
    vertices. Each value is the total rate at which that species leaves that vertex (e.g. for a 
    species with constant diffusion rate D, in a vertex with n neighbours, this value is n*D).
    """
    leaving_rates::Matrix{S}
    """An iterator over all the edges of the lattice."""
    edge_iterator::Vector{Pair{Int64, Int64}}
    """
    The transport rates. This is a dense or sparse matrix (depending on what type of Jacobian is
    used).
    """
    jac_transport::T
    """ Whether sparse jacobian representation is used. """
    sparse::Bool
    """Remove when we add this as problem metadata"""
    lrs::LatticeReactionSystem

    function LatticeTransportODEFunction(ofunc::P, ps::Vector{<:Pair},
            lrs::LatticeReactionSystem, sparse::Bool,
            jac_transport::Union{Nothing, Matrix{S}, SparseMatrixCSC{S, Int64}},
            transport_rates::Vector{Pair{Int64, SparseMatrixCSC{S, Int64}}}) where {P, S}
        # Computes `LatticeTransportODEFunction` functor fields.
        heterogeneous_vert_p_idxs = make_heterogeneous_vert_p_idxs(ps, lrs)
        mtk_ps, p_setters = make_mtk_ps_structs(ps, lrs, heterogeneous_vert_p_idxs)
        t_rate_idx_types, leaving_rates = make_t_types_and_leaving_rates(transport_rates,
            lrs)

        # Creates and returns the `LatticeTransportODEFunction` functor. 
        new{P, typeof(mtk_ps), typeof(p_setters), S, typeof(jac_transport)}(ofunc,
            num_verts(lrs), num_species(lrs), heterogeneous_vert_p_idxs, mtk_ps, p_setters,
            transport_rates, t_rate_idx_types, leaving_rates, Catalyst.edge_iterator(lrs),
            jac_transport, sparse, lrs)
    end
end

# `LatticeTransportODEFunction` helper functions (re-used by rebuild function later on).

# Creates a vector with the heterogeneous vertex parameters' indexes in the full parameter vector.
function make_heterogeneous_vert_p_idxs(ps, lrs)
    p_dict = Dict(ps)
    return findall((p_dict[p] isa Vector) && (length(p_dict[p]) > 1)
    for p in parameters(lrs))
end

# Creates the MTKParameters structure and `p_setters` vector (which are used to manage
# the vertex parameter values during the simulations).
function make_mtk_ps_structs(ps, lrs, heterogeneous_vert_p_idxs)
    p_dict = Dict(ps)
    nonspatial_osys = complete(convert(ODESystem, reactionsystem(lrs)))
    p_init = [p => p_dict[p][1] for p in parameters(nonspatial_osys)]
    mtk_ps = MT.MTKParameters(nonspatial_osys, p_init)
    p_setters = [MT.setp(nonspatial_osys, p)
                 for p in parameters(lrs)[heterogeneous_vert_p_idxs]]
    return mtk_ps, p_setters
end

# Computes the transport rate type vector and leaving rate matrix.
function make_t_types_and_leaving_rates(transport_rates, lrs)
    t_rate_idx_types = [size(tr[2]) == (1, 1) for tr in transport_rates]
    leaving_rates = zeros(length(transport_rates), num_verts(lrs))
    for (s_idx, tr_pair) in enumerate(transport_rates)
        for e in Catalyst.edge_iterator(lrs)
            # Updates the exit rate for species s_idx from vertex e.src.
            leaving_rates[s_idx, e[1]] += get_transport_rate(tr_pair[2], e,
                t_rate_idx_types[s_idx])
        end
    end
    return t_rate_idx_types, leaving_rates
end

### Spatial ODE Functor Functions ###

# Defines the functor's effect when applied as a forcing function.
function (lt_ofun::LatticeTransportODEFunction)(du::AbstractVector, u, p, t)
    # Updates for non-spatial reactions.
    for vert_i in 1:(lt_ofun.num_verts)
        # Gets the indices of all the species at vertex i.
        idxs = get_indexes(vert_i, lt_ofun.num_species)

        # Updates the functors vertex parameter tracker (`mtk_ps`) to contain the vertex parameter 
        # values for vertex vert_i. Then evaluates the reaction contributions to du at vert_i.
        update_mtk_ps!(lt_ofun, p, vert_i)
        lt_ofun.ofunc((@view du[idxs]), (@view u[idxs]), lt_ofun.mtk_ps, t)
    end

    # s_idx is the species index among transport species, s is the index among all species.
    # rates are the species' transport rates.
    for (s_idx, (s, rates)) in enumerate(lt_ofun.transport_rates)
        # Rate for leaving source vertex vert_i. 
        for vert_i in 1:(lt_ofun.num_verts)
            idx_src = get_index(vert_i, s, lt_ofun.num_species)
            du[idx_src] -= lt_ofun.leaving_rates[s_idx, vert_i] * u[idx_src]
        end
        # Add rates for entering a destination vertex via an incoming edge.
        for e in lt_ofun.edge_iterator
            idx_src = get_index(e[1], s, lt_ofun.num_species)
            idx_dst = get_index(e[2], s, lt_ofun.num_species)
            du[idx_dst] += get_transport_rate(rates, e, lt_ofun.t_rate_idx_types[s_idx]) *
                           u[idx_src]
        end
    end
end

# Defines the functor's effect when applied as a Jacobian.
function (lt_ofun::LatticeTransportODEFunction)(J::AbstractMatrix, u, p, t)
    # Resets the Jacobian J's values.
    J .= 0.0

    # Update the Jacobian from non-spatial reaction terms.
    for vert_i in 1:(lt_ofun.num_verts)
        # Gets the indices of all the species at vertex i.
        idxs = get_indexes(vert_i, lt_ofun.num_species)

        # Updates the functors vertex parameter tracker (`mtk_ps`) to contain the vertex parameter 
        # values for vertex vert_i. Then evaluates the reaction contributions to J at vert_i.
        update_mtk_ps!(lt_ofun, p, vert_i)
        lt_ofun.ofunc.jac((@view J[idxs, idxs]), (@view u[idxs]), lt_ofun.mtk_ps, t)
    end

    # Updates for the spatial reactions (adds the Jacobian values from the transportation reactions).
    J .+= lt_ofun.jac_transport
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
        error("Currently lattice ODE simulations are only supported when all spatial reactions are `TransportReaction`s.")
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
    # In the `ODEProblem` vert_ps and edge_ps are merged (but for building the ODEFunction, they are separate).
    vert_ps, edge_ps = lattice_process_p(p_in, vertex_parameters(lrs),
        edge_parameters(lrs), lrs)

    # Creates the ODEFunction.
    ofun = build_odefunction(lrs, vert_ps, edge_ps, jac, sparse, name, include_zero_odes,
        combinatoric_ratelaws, remove_conserved, checks)

    # Combines `vert_ps` and `edge_ps` to a single vector with values only (not a map). Creates ODEProblem.
    pval_dict = Dict([vert_ps; edge_ps])
    ps = [pval_dict[p] for p in parameters(lrs)]
    return ODEProblem(ofun, u0, tspan, ps, args...; kwargs...)
end

# Builds an ODEFunction for a spatial ODEProblem.
function build_odefunction(lrs::LatticeReactionSystem, vert_ps::Vector{Pair{R, Vector{T}}},
        edge_ps::Vector{Pair{S, SparseMatrixCSC{T, Int64}}},
        jac::Bool, sparse::Bool, name, include_zero_odes, combinatoric_ratelaws,
        remove_conserved, checks) where {R, S, T}
    # Error check.
    if remove_conserved
        throw(ArgumentError("Removal of conserved quantities is currently not supported for `LatticeReactionSystem`s"))
    end

    # Prepares the inputs to the `LatticeTransportODEFunction` functor. 
    osys = complete(convert(ODESystem, reactionsystem(lrs);
        name, combinatoric_ratelaws, include_zero_odes, checks))
    ofunc_dense = ODEFunction(osys; jac = true, sparse = false)
    ofunc_sparse = ODEFunction(osys; jac = true, sparse = true)
    transport_rates = make_sidxs_to_transrate_map(vert_ps, edge_ps, lrs)

    # Depending on Jacobian and sparsity options, compute the Jacobian transport matrix and prototype.
    if !sparse && !jac
        jac_transport = nothing
        jac_prototype = nothing
    else
        jac_sparse = build_jac_prototype(ofunc_sparse.jac_prototype, transport_rates, lrs;
            set_nonzero = jac)
        jac_dense = Matrix(jac_sparse)
        jac_transport = (jac ? (sparse ? jac_sparse : jac_dense) : nothing)
        jac_prototype = (sparse ? jac_sparse : nothing)
    end

    # Creates the `LatticeTransportODEFunction` functor (if `jac`, sets it as the Jacobian as well).
    f = LatticeTransportODEFunction(ofunc_dense, [vert_ps; edge_ps], lrs, sparse,
        jac_transport, transport_rates)
    J = (jac ? f : nothing)

    # Extracts the `Symbol` form for parameters (but not species). Creates and returns the `ODEFunction`.
    paramsyms = [MT.getname(p) for p in parameters(lrs)]
    sys = SciMLBase.SymbolCache([], paramsyms, [])
    return ODEFunction(f; jac = J, jac_prototype, sys)
end

# Builds a Jacobian prototype.
# If requested, populate it with the constant values of the Jacobian's transportation part.
function build_jac_prototype(ns_jac_prototype::SparseMatrixCSC{Float64, Int64},
        transport_rates::Vector{Pair{Int64, SparseMatrixCSC{T, Int64}}},
        lrs::LatticeReactionSystem; set_nonzero = false) where {T}
    # Finds the indices of  both the transport species,
    # and the species with transport only (that is, with no non-spatial dynamics but with spatial dynamics).
    trans_species = [tr[1] for tr in transport_rates]
    trans_only_species = filter(s_idx -> !Base.isstored(ns_jac_prototype, s_idx, s_idx),
        trans_species)

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

    # Create a sparse Jacobian prototype with 0-valued entries. If requested,
    # updates values with non-zero entries.
    jac_prototype = sparse(i_idxs, j_idxs, zeros(T, num_entries))
    set_nonzero && set_jac_transport_values!(jac_prototype, transport_rates, lrs)

    return jac_prototype
end

# For a Jacobian prototype with zero-valued entries. Set entry values according to a set of
# transport reaction values.
function set_jac_transport_values!(jac_prototype, transport_rates, lrs)
    for (s, rates) in transport_rates, e in edge_iterator(lrs)
        idx_src = get_index(e[1], s, num_species(lrs))
        idx_dst = get_index(e[2], s, num_species(lrs))
        val = get_transport_rate(rates, e, size(rates) == (1, 1))

        # Term due to species leaving source vertex.
        jac_prototype[idx_src, idx_src] -= val

        # Term due to species arriving to destination vertex.   
        jac_prototype[idx_src, idx_dst] += val
    end
end
