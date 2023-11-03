### Spatial ODE Functor Structures ###

# Functor structure containg the information for the forcing function of a spatial ODE with spatial movement on a lattice.
struct LatticeDiffusionODEf{R,S,T}
    """The ODEFunction of the (non-spatial) reaction system which generated this function."""
    ofunc::R
    """The number of vertices."""
    nV::Int64
    """The number of species."""
    nS::Int64
    """The values of the parameters which values are tied to vertexes."""
    pV::Vector{Vector{Float64}}
    """Temporary vector. For parameters which values are identical across the lattice, at some point these have to be converted of a length(nV) vector. To avoid re-allocation they are written to this vector."""
    work_pV::Vector{Float64}    
    """For each parameter in pV, it either have length nV or 1. To know whenever a parameter's value need expanding to the work_pV array, its length needs checking. This check is done once, and the value stored to this array. This field (specifically) is an enumerate over that array."""
    enumerated_pV_idx_types::Base.Iterators.Enumerate{BitVector}
    """A vector of pairs, with a value for each species with transportation. The first value is the species index (in the species(::ReactionSystem) vector), and the second is a vector with its diffusion rate values. If the diffusion rate is uniform, that value is the only value in the vector. Else, there is one value for each edge in the lattice."""
    spatial_rates::Vector{Pair{Int64, Vector{Float64}}}
    """A matrix, NxM, where N is the number of species with transportation and M the number of vertexes. Each value is the total rate at which that species leaves that vertex (e.g. for a species with constant diffusion rate D, in a vertex with n neighbours, this value is n*D)."""
    leaving_rates::Matrix{Float64}
    """An (enumerate'ed) itarator over all the edges of the lattice."""
    enumerated_edges::T
    
    function LatticeDiffusionODEf(ofunc::R, pV, spatial_rates::Vector{Pair{Int64, Vector{Float64}}}, lrs::LatticeReactionSystem) where {R, S}
        leaving_rates = zeros(length(spatial_rates), lrs.nV)            # Initialises the leaving rates matrix with zeros.
        for (s_idx, rates) in enumerate(last.(spatial_rates)),
            (e_idx, e) in enumerate(edges(lrs.lattice))                 # Iterates through all edges, and all spatial rates (map from each diffusing species to its rates across edges).
    
            leaving_rates[s_idx, e.src] += get_component_value(rates, e_idx)    # Updates the leaving rate for that combination of vertex and species. RHS finds the value of edge "e_idx" in the vector of diffusion rates ("rates").
        end
        work_pV = zeros(lrs.nV)                                         # Initialises the work pV vector to be empty.
        enumerated_pV_idx_types = enumerate(length.(pV) .== 1)          # Creates a Boolean vector whether each vertex parameter need expanding or (an enumerates it, since it always appear in this form).
        enumerated_edges = deepcopy(enumerate(edges(lrs.lattice)))      # Creates an iterator over all the edges. Again, this is always used in the enumerated form. 
        new{R,S,typeof(enumerated_edges)}(ofunc, lrs.nV, lrs.nS, pV, work_pV, enumerated_pV_idx_types, spatial_rates, leaving_rates, enumerated_edges)
    end
end

# Functor structure containing the information for the forcing function of a spatial ODE with spatial movement on a lattice.
struct LatticeDiffusionODEjac{S,T}
    """The ODEFunction of the (non-spatial) reaction system which generated this function."""
    ofunc::S
    """The number of vertices."""
    nV::Int64
    """The number of species."""
    nS::Int64
    """The values of the parameters which values are tied to vertexes."""
    pV::Vector{Vector{Float64}}
    """Temporary vector. For parameters which values are identical across the lattice, at some point these have to be converted of a length(nV) vector. To avoid re-allocation they are written to this vector."""
    work_pV::Vector{Float64} 
    """For each parameter in pV, it either have length nV or 1. To know whenever a parameter's value need expanding to the work_pV array, its length needs checking. This check is done once, and the value stored to this array. This field (specifically) is an enumerate over that array."""
    enumerated_pV_idx_types::Base.Iterators.Enumerate{BitVector}
    """Whether the Jacobian is sparse or not."""
    sparse::Bool
    """The values of the Jacobian. All the diffusion rates. Eitehr in matrix form (for non-sparse, in this case with potential zeros) or as the "nzval" field of the sparse jacobian matrix."""
    jac_values::T

    function LatticeDiffusionODEjac(ofunc::S, pV, lrs::LatticeReactionSystem, jac_prototype::Union{Nothing, SparseMatrixCSC{Float64, Int64}}, sparse::Bool) where {S, T}
        work_pV = zeros(lrs.nV)                                             # Initialises the work pV vector to be empty.
        enumerated_pV_idx_types = enumerate(length.(pV) .== 1)              # Creates a Boolean vector whether each vertex parameter need expanding or (an enumerates it, since it always appear in this form).
        jac_values = sparse ? jac_prototype.nzval : Matrix(jac_prototype)   # Retrieves the diffusion values (form depending on Jacobian sparsity).
        new{S,typeof(jac_values)}(ofunc, lrs.nV, lrs.nS, pV, work_pV, enumerated_pV_idx_types, sparse, jac_values)
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
    u0 = lattice_process_u0(u0_in, species(lrs), lrs.nV)                                    # u0 becomes a vector ([species 1 at vertex 1, species 2 at vertex 1, ..., species 1 at vertex 2, ...])
    pV, pE = lattice_process_p(p_in, vertex_parameters(lrs), edge_parameters(lrs), lrs)     # Both pV and pE becomes vectors of vectors. Each have 1 element for each parameter. These elements are length 1 vectors (if the parameter is uniform), or length nV/nE, with unique values for each vertex/edge (for pV/pE, respectively).

    # Creates ODEProblem.
    ofun = build_odefunction(lrs, pV, pE, jac, sparse)          # Builds the ODEFunction.
    return ODEProblem(ofun, u0, tspan, pV, args...; kwargs...)  # Creates a normal ODEProblem.
end

# Builds an ODEFunction for a spatial ODEProblem.
function build_odefunction(lrs::LatticeReactionSystem, pV::Vector{Vector{Float64}},
                           pE::Vector{Vector{Float64}}, use_jac::Bool, sparse::Bool)
    # Prepares (non-spatial) ODE functions and list of spatially moving species and their rates.
    ofunc = ODEFunction(convert(ODESystem, lrs.rs); jac = use_jac, sparse = false)              # Creates the (non-spatial) ODEFunction corresponding to the (non-spatial) reaction network.
    ofunc_sparse = ODEFunction(convert(ODESystem, lrs.rs); jac = use_jac, sparse = true)        # Creates the same function, but sparse. Could insert so this is only computed for sparse cases.
    spatial_rates_speciesmap = compute_all_spatial_rates(pV, pE, lrs)                           # Creates a map (Vector{Pair}), mapping each species that is transported to a vector with its transportation rate. If the rate is uniform across all edges, the vector will be length 1 (with this value), else there will be a separate value for each edge.
    spatial_rates = [findfirst(isequal(spat_rates[1]), species(lrs)) => spat_rates[2]
                       for spat_rates in spatial_rates_speciesmap]                              # Remakes "spatial_rates_speciesmap". Rates are identical, but the species are represented as their index (in the species(::ReactionSystem) vector). In "spatial_rates_speciesmap" they instead were Symbolics.

    f = LatticeDiffusionODEf(ofunc, pV, spatial_rates, lrs)                                     # Creates a functor for the ODE f function (incorporating spatial and non-spatial reactions).
    jac_prototype = (use_jac || sparse) ?
                    build_jac_prototype(ofunc_sparse.jac_prototype, spatial_rates,
                                        lrs; set_nonzero = use_jac) : nothing                   # Computes the Jacobian prototype (nothing if `jac=false`). 
    jac = use_jac ? LatticeDiffusionODEjac(ofunc, pV, lrs, jac_prototype, sparse) : nothing     # (Potentially) Creates a functor for the ODE Jacobian function (incorporating spatial and non-spatial reactions).
    return ODEFunction(f; jac = jac, jac_prototype = (sparse ? jac_prototype : nothing))        # Creates the ODEFunction used in the ODEProblem.
end

# Builds a jacobian prototype. If requested, populate it with the Jacobian's (constant) values as well.
function build_jac_prototype(ns_jac_prototype::SparseMatrixCSC{Float64, Int64}, spatial_rates, lrs::LatticeReactionSystem;
                             set_nonzero = false)
    spat_species = first.(spatial_rates)        # Gets a list of species with transportation (As their index in the full species vector. If you have species [X(t), Y(t)], with transportation for Y only, this becomes [2]).

    # Gets list of indexes for species that move spatially, but are involved in no other reaction.
    only_spat = [(s in spat_species) && !Base.isstored(ns_jac_prototype, s, s)
                 for s in 1:(lrs.nS)]           # Probably rare, but these creates weird special cases where there block diagonal part of the Jacobian have empty spaces.

    # Declares sparse array content.
    J_colptr = fill(1, lrs.nV * lrs.nS + 1)     # Initiates the column values of the (sparse) jacobian. 
    J_nzval = fill(0.0,                         # Initiates the row values of the (sparse) jacobian. Has one value for each value in teh sparse jacobian.            
                   lrs.nV * (nnz(ns_jac_prototype) +                        # The number of values due to the non-spatial jacobian (its number of values, and multiplied by the number of vertexes).
                   count(only_spat)) +                                      # 
                   length(edges(lrs.lattice)) * length(spatial_rates))      # 
    J_rowval = fill(0, length(J_nzval))                                      

    # Finds filled elements.
    for comp in 1:(lrs.nV), s in 1:(lrs.nS)                                 # Loops through all vertexes and species.                
        col_idx = get_index(comp, s, lrs.nS)                                # For the current vertex+species, finds its column (same as its index in the `u` vector).

        # Column values.
        local_elements = in(s, spat_species) *
                         (length(lrs.lattice.fadjlist[comp]) + only_spat[s])
        spatial_elements = -(ns_jac_prototype.colptr[(s + 1):-1:s]...)
        J_colptr[col_idx + 1] = J_colptr[col_idx] + local_elements + spatial_elements

        # Row values.
        rows = ns_jac_prototype.rowval[ns_jac_prototype.colptr[s]:(ns_jac_prototype.colptr[s + 1] - 1)] .+
               (comp - 1) * lrs.nS
        if in(s, spat_species)
            # Finds the location of the spatial_elements, and inserts the elements from the non-spatial part into this.
            spatial_rows = (lrs.lattice.fadjlist[comp] .- 1) .* lrs.nS .+ s
            split_idx = isempty(rows) ? 1 : findfirst(spatial_rows .> rows[1])
            isnothing(split_idx) && (split_idx = length(spatial_rows) + 1)
            rows = vcat(spatial_rows[1:(split_idx - 1)], rows,
            spatial_rows[split_idx:end])
            if only_spat[s]
                split_idx = findfirst(rows .> get_index(comp, s, lrs.nS))
                isnothing(split_idx) && (split_idx = length(rows) + 1)
                insert!(rows, split_idx, get_index(comp, s, lrs.nS))
            end
        end
        J_rowval[J_colptr[col_idx]:(J_colptr[col_idx + 1] - 1)] = rows
    end
    
    # Set element values.
    if !set_nonzero
        J_nzval .= 1.0
    else
        for (s_idx, (s, rates)) in enumerate(spatial_rates),
            (e_idx, edge) in enumerate(edges(lrs.lattice))

            col_start = J_colptr[get_index(edge.src, s, lrs.nS)]
            col_end = J_colptr[get_index(edge.src, s, lrs.nS) + 1] - 1
            column_view = @view J_rowval[col_start:col_end]

            # Updates the source value.
            val_idx_src = col_start +
                          findfirst(column_view .== get_index(edge.src, s, lrs.nS)) - 1
            J_nzval[val_idx_src] -= get_component_value(rates, e_idx)

            # Updates the destination value.
            val_idx_dst = col_start +
                          findfirst(column_view .== get_index(edge.dst, s, lrs.nS)) - 1
            J_nzval[val_idx_dst] += get_component_value(rates, e_idx)
        end
    end

    return SparseMatrixCSC(lrs.nS * lrs.nV, lrs.nS * lrs.nV, J_colptr, J_rowval, J_nzval)
end

# Defines the forcing functor's effect on the (spatial) ODE system.
function (f_func::LatticeDiffusionODEf)(du, u, p, t)
    # Updates for non-spatial reactions.
    for comp_i::Int64 in 1:(f_func.nV)                              # Loops through each vertex of the lattice. Applies the (non-spatial) ODEFunction to the species in that vertex.
        f_func.ofunc((@view du[get_indexes(comp_i, f_func.nS)]),    # Get the indexes of the i'th vertex (current one in the loop) in the u vector. Uses this to create a view of the vector to which the new species concentrations are written.  
                     (@view u[get_indexes(comp_i, f_func.nS)]),     # Same as above, but reads the current species concentrations.
                      view_pV_vector!(f_func.work_pV, p, comp_i, f_func.enumerated_pV_idx_types),   # Gets a vector with the values of the (vertex) parameters in the current vertex.
                      t)                                            # Time.
    end

    # Updates for spatial reactions.
    for (s_idx, (s, rates)) in enumerate(f_func.spatial_rates)      # Loops through all species with transportation. Here: s_idx is its index among the species with transportations. s is its index among all species (in the species(::ReactionSystem) vector). rates is its rates values (vector length 1 if uniform, else same length as the number of edges).
        for comp_i::Int64 in 1:(f_func.nV)                          # Loops through all vertexes.
            du[get_index(comp_i, s, f_func.nS)] -= f_func.leaving_rates[s_idx, comp_i] *
                                                u[get_index(comp_i, s,
                                                f_func.nS)]         # Finds the leaving rate of this species in this vertex. Updates the du vector at that vertex/species combination with the corresponding rate (leaving rate times concentration).
        end
        for (e_idx::Int64, edge::Graphs.SimpleGraphs.SimpleEdge{Int64}) in f_func.enumerated_edges  # Loops through all edges.
            du[get_index(edge.dst, s, f_func.nS)] += get_component_value(rates, e_idx) *
                                                  u[get_index(edge.src, s,
                                                  f_func.nS)]                                       # For the destination of this edge, we want to add the influx term to du. This is ["rates" value for this edge]*[the species concentration in the source vertex]. 
        end
    end
end

# Defines the jacobian functor's effect on the (spatial) ODE system.
function (jac_func::LatticeDiffusionODEjac)(J, u, p, t)
    # Because of weird stuff where the Jacobian is not reset that I don't understand properly.
    reset_J_vals!(J)    # Sets all Jacobian values to 0 (because they are not by default, this is weird but and I could not get it to work otherwise, tried to get Chris to explain but he wouldn't. Hopefully this can be improved once I get him to explain).

    # Updates for non-spatial reactions.
    for comp_i::Int64 in 1:(jac_func.nV)                                # Loops through all vertexes and applies the (non-spatial) Jacobian to the species in that vertex.
        jac_func.ofunc.jac((@view J[get_indexes(comp_i, jac_func.nS),
                           get_indexes(comp_i, jac_func.nS)]),
                  (@view u[get_indexes(comp_i, jac_func.nS)]),
                  view_pV_vector!(jac_func.work_pV, p, comp_i, jac_func.enumerated_pV_idx_types), t)# These inputs are the same as when f_func.ofunc was applied in the previous block.
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