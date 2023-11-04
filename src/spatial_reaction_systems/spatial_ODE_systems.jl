### Spatial ODE Functor Structures ###

# Functor structure containing the information for the forcing function of a spatial ODE with spatial movement on a lattice.
struct LatticeDiffusionODEf{S,T}
    """The ODEFunction of the (non-spatial) reaction system which generated this function."""
    ofunc::S
    """The number of vertices."""
    num_verts::Int64
    """The number of species."""
    num_species::Int64
    """The values of the parameters which values are tied to vertexes."""
    vert_ps::Vector{Vector{Float64}}
    """Temporary vector. For parameters which values are identical across the lattice, at some point these have to be converted of a length num_verts vector. To avoid re-allocation they are written to this vector."""
    work_vert_ps::Vector{Float64}    
    """For each parameter in vert_ps, its value is a vector with length either num_verts or 1. To know whenever a parameter's value need expanding to the work_vert_ps array, its length needs checking. This check is done once, and the value stored to this array. This field (specifically) is an enumerate over that array."""
    enum_v_ps_idx_types::Base.Iterators.Enumerate{BitVector}
    """A vector of pairs, with a value for each species with transportation. The first value is the species index (in the species(::ReactionSystem) vector), and the second is a vector with its transport rate values. If the transport rate is uniform (across all edges), that value is the only value in the vector. Else, there is one value for each edge in the lattice."""
    transport_rates::Vector{Pair{Int64, Vector{Float64}}}
    """A matrix, NxM, where N is the number of species with transportation and M the number of vertexes. Each value is the total rate at which that species leaves that vertex (e.g. for a species with constant diffusion rate D, in a vertex with n neighbours, this value is n*D)."""
    leaving_rates::Matrix{Float64}
    """An (enumerate'ed) iterator over all the edges of the lattice."""
    enum_edges::T
    
    function LatticeDiffusionODEf(ofunc::S, vert_ps, transport_rates::Vector{Pair{Int64, Vector{Float64}}}, lrs::LatticeReactionSystem) where {S}
        leaving_rates = zeros(length(transport_rates), lrs.num_verts)           # Initialises the leaving rates matrix with zeros.
        for (s_idx, rates) in enumerate(last.(transport_rates)),
            (e_idx, e) in enumerate(edges(lrs.lattice))                         # Iterates through all edges, and all transport rates (map from each diffusing species to its rates across edges).
    
            leaving_rates[s_idx, e.src] += get_component_value(rates, e_idx)    # Updates the leaving rate for that combination of vertex and species. RHS finds the value of edge "e_idx" in the vector of diffusion rates ("rates").
        end
        work_vert_ps = zeros(lrs.num_verts)                                     # Initialises the work vert_ps vector to be empty.
        enum_v_ps_idx_types = enumerate(length.(vert_ps) .== 1)                 # Creates a Boolean vector whether each vertex parameter need expanding or (and enumerates it, since it always appear in this form).
        enum_edges = deepcopy(enumerate(edges(lrs.lattice)))                    # Creates an iterator over all the edges. Again, this is always used in the enumerated form. 
        new{S,typeof(enum_edges)}(ofunc, lrs.num_verts, lrs.num_species, vert_ps, work_vert_ps, enum_v_ps_idx_types, transport_rates, leaving_rates, enum_edges)
    end
end

# Functor structure containing the information for the forcing function of a spatial ODE with spatial movement on a lattice.
struct LatticeDiffusionODEjac{S,T}
    """The ODEFunction of the (non-spatial) reaction system which generated this function."""
    ofunc::S
    """The number of vertices."""
    num_verts::Int64
    """The number of species."""
    num_species::Int64
    """The values of the parameters which values are tied to vertexes."""
    vert_ps::Vector{Vector{Float64}}
    """Temporary vector. For parameters which values are identical across the lattice, at some point these have to be converted of a length(num_verts) vector. To avoid re-allocation they are written to this vector."""
    work_vert_ps::Vector{Float64} 
    """For each parameter in vert_ps, it either have length num_verts or 1. To know whenever a parameter's value need expanding to the work_vert_ps array, its length needs checking. This check is done once, and the value stored to this array. This field (specifically) is an enumerate over that array."""
    enum_v_ps_idx_types::Base.Iterators.Enumerate{BitVector}
    """Whether the Jacobian is sparse or not."""
    sparse::Bool
    """The values of the Jacobian. All the diffusion rates. Eitehr in matrix form (for non-sparse, in this case with potential zeros) or as the "nzval" field of the sparse jacobian matrix."""
    jac_values::T

    function LatticeDiffusionODEjac(ofunc::S, vert_ps, lrs::LatticeReactionSystem, jac_prototype::Union{Nothing, SparseMatrixCSC{Float64, Int64}}, sparse::Bool) where {S}
        work_vert_ps = zeros(lrs.num_verts)                                 # Initialises the work vert_ps vector to be empty.
        enum_v_ps_idx_types = enumerate(length.(vert_ps) .== 1)             # Creates a Boolean vector whether each vertex parameter need expanding or (an enumerates it, since it always appear in this form).
        jac_values = sparse ? jac_prototype.nzval : Matrix(jac_prototype)   # Retrieves the diffusion values (form depending on Jacobian sparsity).
        new{S,typeof(jac_values)}(ofunc, lrs.num_verts, lrs.num_species, vert_ps, work_vert_ps, enum_v_ps_idx_types, sparse, jac_values)
    end
end

### ODEProblem ###

# Creates an ODEProblem from a LatticeReactionSystem.
function DiffEqBase.ODEProblem(lrs::LatticeReactionSystem, u0_in, tspan,
                               p_in = DiffEqBase.NullParameters(), args...;
                               jac = true, sparse = jac, kwargs...)
    is_transport_system(lrs) || error("Currently lattice ODE simulations only supported when all spatial reactions are transport reactions.")
    
    # Converts potential symmaps to varmaps.
    u0_in = symmap_to_varmap(lrs, u0_in)
    p_in = (p_in isa Tuple{<:Any,<:Any}) ? (symmap_to_varmap(lrs, p_in[1]),symmap_to_varmap(lrs, p_in[2])) : symmap_to_varmap(lrs, p_in)    # Parameters can be given in Tuple form (where the first element is the vertex parameters and the second the edge parameters). In this case, we have to covert each element separately.

    # Converts u0 and p to their internal forms.
    u0 = lattice_process_u0(u0_in, species(lrs), lrs.num_verts)                                     # u0 becomes a vector ([species 1 at vertex 1, species 2 at vertex 1, ..., species 1 at vertex 2, ...]).
    vert_ps, edge_ps = lattice_process_p(p_in, vertex_parameters(lrs), edge_parameters(lrs), lrs)   # Both vert_ps and edge_ps becomes vectors of vectors. Each have 1 element for each parameter. These elements are length 1 vectors (if the parameter is uniform), or length num_verts/nE, with unique values for each vertex/edge (for vert_ps/edge_ps, respectively).

    # Creates ODEProblem.
    ofun = build_odefunction(lrs, vert_ps, edge_ps, jac, sparse)        # Builds the ODEFunction.
    return ODEProblem(ofun, u0, tspan, vert_ps, args...; kwargs...)     # Creates a normal ODEProblem.
end

# Builds an ODEFunction for a spatial ODEProblem.
function build_odefunction(lrs::LatticeReactionSystem, vert_ps::Vector{Vector{Float64}},
                           edge_ps::Vector{Vector{Float64}}, use_jac::Bool, sparse::Bool)
    # Prepares (non-spatial) ODE functions and list of spatially moving species and their rates.
    ofunc = ODEFunction(convert(ODESystem, lrs.rs); jac = use_jac, sparse = false)                  # Creates the (non-spatial) ODEFunction corresponding to the (non-spatial) reaction network.
    ofunc_sparse = ODEFunction(convert(ODESystem, lrs.rs); jac = use_jac, sparse = true)            # Creates the same function, but sparse. Could insert so this is only computed for sparse cases.
    transport_rates_speciesmap = compute_all_transport_rates(vert_ps, edge_ps, lrs)                 # Creates a map (Vector{Pair}), mapping each species that is transported to a vector with its transportation rate. If the rate is uniform across all edges, the vector will be length 1 (with this value), else there will be a separate value for each edge.
    transport_rates = [findfirst(isequal(spat_rates[1]), species(lrs)) => spat_rates[2]
                       for spat_rates in transport_rates_speciesmap]                                # Remakes "transport_rates_speciesmap". Rates are identical, but the species are represented as their index (in the species(::ReactionSystem) vector). In "transport_rates_speciesmap" they instead were Symbolics.

    f = LatticeDiffusionODEf(ofunc, vert_ps, transport_rates, lrs)                                  # Creates a functor for the ODE f function (incorporating spatial and non-spatial reactions).
    jac_prototype = (use_jac || sparse) ?
                    build_jac_prototype(ofunc_sparse.jac_prototype, transport_rates,
                                        lrs; set_nonzero = use_jac) : nothing                       # Computes the Jacobian prototype (nothing if `jac=false`). 
    jac = use_jac ? LatticeDiffusionODEjac(ofunc, vert_ps, lrs, jac_prototype, sparse) : nothing    # (Potentially) Creates a functor for the ODE Jacobian function (incorporating spatial and non-spatial reactions).
    return ODEFunction(f; jac = jac, jac_prototype = (sparse ? jac_prototype : nothing))            # Creates the ODEFunction used in the ODEProblem.
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
    non_spat_terms = [[get_index(vert, s_i, lrs.num_species), get_index(vert, s_j, lrs.num_species)] for vert in 1:(lrs.num_verts) for (s_i, s_j) in zip(ns_i_idxs,ns_j_idxs)]               # Indexes of elements due to non-spatial dynamics.
    trans_only_leaving_terms = [[get_index(e.src, s_idx, lrs.num_species), get_index(e.src, s_idx, lrs.num_species)] for e in edges(lrs.lattice) for s_idx in trans_only_species]     # Indexes due to terms for a species leaves its current vertex (but does not have non-spatial dynamics). If the non-spatial Jacobian is fully dense, these would already be accounted for. 
    trans_arriving_terms = [[get_index(e.src, s_idx, lrs.num_species), get_index(e.dst, s_idx, lrs.num_species)] for e in edges(lrs.lattice) for s_idx in trans_species]              # Indexes due to terms for species arriving into a new vertex.
    all_terms = [non_spat_terms; trans_only_leaving_terms; trans_arriving_terms]

    # Creates a jacobian prototype with 0 values in all positions).
    jac_prototype = sparse(first.(all_terms), last.(all_terms), fill(0.0, length(all_terms)))

    # Set element values.
    if set_nonzero
        for (s, rates) in trans_rates, (e_idx, e) in enumerate(edges(lrs.lattice))    # Loops through all species with transportation and all edges along which the can be transported.
            jac_prototype[get_index(e.src, s, lrs.num_species), get_index(e.src, s, lrs.num_species)] -= get_component_value(rates, e_idx)    # Term due to species leaving source vertex.
            jac_prototype[get_index(e.src, s, lrs.num_species), get_index(e.dst, s, lrs.num_species)] += get_component_value(rates, e_idx)    # Term due to species arriving to destination vertex.
        end
    end

    return jac_prototype
end

# Defines the forcing functor's effect on the (spatial) ODE system.
function (f_func::LatticeDiffusionODEf)(du, u, p, t)
    # Updates for non-spatial reactions.
    for comp_i::Int64 in 1:(f_func.num_verts)                                # Loops through each vertex of the lattice. Applies the (non-spatial) ODEFunction to the species in that vertex.
        f_func.ofunc((@view du[get_indexes(comp_i, f_func.num_species)]),    # Get the indexes of the i'th vertex (current one in the loop) in the u vector. Uses this to create a view of the vector to which the new species concentrations are written.  
                     (@view u[get_indexes(comp_i, f_func.num_species)]),     # Same as above, but reads the current species concentrations.
                      view_vert_ps_vector!(f_func.work_vert_ps, p, comp_i, f_func.enum_v_ps_idx_types),   # Gets a vector with the values of the (vertex) parameters in the current vertex.
                      t)                                                     # Time.
    end

    # Updates for spatial reactions.
    for (s_idx, (s, rates)) in enumerate(f_func.transport_rates)            # Loops through all species with transportation. Here: s_idx is its index among the species with transportations. s is its index among all species (in the species(::ReactionSystem) vector). rates is its rates values (vector length 1 if uniform, else same length as the number of edges).
        for comp_i::Int64 in 1:(f_func.num_verts)                           # Loops through all vertexes.
            du[get_index(comp_i, s, f_func.num_species)] -= f_func.leaving_rates[s_idx, comp_i] *
                                                u[get_index(comp_i, s,
                                                f_func.num_species)]        # Finds the leaving rate of this species in this vertex. Updates the du vector at that vertex/species combination with the corresponding rate (leaving rate times concentration).
        end
        for (e_idx::Int64, edge::Graphs.SimpleGraphs.SimpleEdge{Int64}) in f_func.enum_edges  # Loops through all edges.
            du[get_index(edge.dst, s, f_func.num_species)] += get_component_value(rates, e_idx) *
                                                  u[get_index(edge.src, s,
                                                  f_func.num_species)]      # For the destination of this edge, we want to add the influx term to du. This is ["rates" value for this edge]*[the species concentration in the source vertex]. 
        end
    end
end

# Defines the jacobian functor's effect on the (spatial) ODE system.
function (jac_func::LatticeDiffusionODEjac)(J, u, p, t)
    # Because of weird stuff where the Jacobian is not reset that I don't understand properly.
    reset_J_vals!(J)    # Sets all Jacobian values to 0 (because they are not by default, this is weird but and I could not get it to work otherwise, tried to get Chris to explain but he wouldn't. Hopefully this can be improved once I get him to explain).

    # Updates for non-spatial reactions.
    for comp_i::Int64 in 1:(jac_func.num_verts)                                # Loops through all vertexes and applies the (non-spatial) Jacobian to the species in that vertex.
        jac_func.ofunc.jac((@view J[get_indexes(comp_i, jac_func.num_species),
                           get_indexes(comp_i, jac_func.num_species)]),
                           (@view u[get_indexes(comp_i, jac_func.num_species)]),
                           view_vert_ps_vector!(jac_func.work_vert_ps, p, comp_i, jac_func.enum_v_ps_idx_types), t) # These inputs are the same as when f_func.ofunc was applied in the previous block.
    end

    # Updates for the spatial reactions.
    add_spat_J_vals!(J, jac_func)       # Adds the Jacobian values from the diffusion reactions.
end
# Resets the jacobian matrix within a jac call. Separate for spatial and non-spatial cases.
reset_J_vals!(J::Matrix) = (J .= 0.0)                   
reset_J_vals!(J::SparseMatrixCSC) = (J.nzval .= 0.0)
# Updates the jacobian matrix with the difusion values. Separate for spatial and non-spatial cases.
add_spat_J_vals!(J::SparseMatrixCSC, jac_func::LatticeDiffusionODEjac) = (J.nzval .+= jac_func.jac_values)
add_spat_J_vals!(J::Matrix, jac_func::LatticeDiffusionODEjac) = (J .+= jac_func.jac_values)