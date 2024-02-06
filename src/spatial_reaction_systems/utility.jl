### Processes Input u0 & p ###

# Defines _symbol_to_var, but where the system is a LRS. Required to make symmapt_to_varmap to work.
function _symbol_to_var(lrs::LatticeReactionSystem, sym)
    # Checks if sym is a parameter.
    p_idx = findfirst(sym==p_sym for p_sym in ModelingToolkit.getname.(parameters(lrs)))    
    isnothing(p_idx) || return parameters(lrs)[p_idx]

    # Checks if sym is a species.
    s_idx = findfirst(sym==s_sym for s_sym in ModelingToolkit.getname.(species(lrs)))       
    isnothing(s_idx) || return species(lrs)[s_idx]

    error("Could not find property parameter/species $sym in lattice reaction system.")
end

# From u0 input, extracts their values and store them in the internal format.
# Internal format: a vector  on the form [spec 1 at vert 1, spec 2 at vert 1, ..., spec 1 at vert 2, ...]).
function lattice_process_u0(u0_in, u0_syms, lrs::LatticeReactionSystem)
    # u0 values can be given in various forms. This converts it to a Vector{Pair{Symbolics,...}} form.
    # Top-level vector: Maps each species to its value(s).
    u0 = lattice_process_input(u0_in, u0_syms)   

    # Species' initial condition values can be given in different forms (also depending on the lattice).
    # This converts each species's values to a Vector. For species with uniform initial conditions,
    # The value vector holds that value only. For spatially heterogeneous initial conditions,
    # the vector have teh same length as the number of vertexes (with one value for each).
    u0 = vertex_value_map(u0, lrs)

    # Converts the initial condition to a single Vector (with one values for each species and vertex).
    return expand_component_values([entry[2] for entry in u0], lrs.num_verts)                  
end

# From p input, splits it into diffusion parameters and compartment parameters.
# Store these in the desired internal format.
function lattice_process_p(ps_in, ps_vertex_syms, ps_edge_syms, lrs::LatticeReactionSystem)
    # p values can be given in various forms. This converts it to a Vector{Pair{Symbolics,...}} form.
    # Top-level vector: Maps each parameter to its value(s).
    # Second-level: Contains either a vector (vertex parameters) or a sparse matrix (edge parameters).
    # For uniform parameters these have size 1/1x1. Else, they have size num_verts/num_vertsxnum_verts.
    ps = lattice_process_input(ps_in, [ps_vertex_syms; ps_edge_syms])   

    # Split the parameter vector into one for vertex parameters and one for edge parameters. 
    # Next, converts the values to the correct form (vectors for vert_ps and sparse matrices for edge_ps).
    vert_ps, edge_ps = split_parameters(ps, ps_vertex_syms, ps_edge_syms)     
    vert_ps = vertex_value_map(vert_ps, lrs)
    edge_ps = edge_value_map(edge_ps, lrs)

    return vert_ps, edge_ps
end

# The input (parameters or initial conditions) may either be a dictionary (symbolics to value(s).)
# or a map (in vector or tuple form) from symbolics to value(s). This converts the input to a
# (Vector) map from symbolics to value(s), where the entries have the same order as `syms`.
function lattice_process_input(input::Dict{BasicSymbolic{Real}, <:Any}, syms::Vector{BasicSymbolic{Real}})
    # Error checks
    if !isempty(setdiff(keys(input), syms))
        error("You have provided values for the following unrecognised parameters/initial conditions: $(setdiff(keys(input), syms)).")
    end
    if !isempty(setdiff(syms, keys(input)))
        error("You have not provided values for the following parameters/initial conditions: $(setdiff(syms, keys(input))).")
    end

    return [sym => input[sym] for sym in syms]
end
function lattice_process_input(input, syms::Vector{BasicSymbolic{Real}}) 
    lattice_process_input(Dict(input), syms)
end

# Splits parameters into vertex and edge parameters.
#function split_parameters(ps::Vector{<: Pair}, p_vertex_syms::Vector, p_edge_syms::Vector)  
function split_parameters(ps, p_vertex_syms, p_edge_syms)  
    vert_ps = [p for p in ps if any(isequal(p[1]), p_vertex_syms)]     
    edge_ps = [p for p in ps if any(isequal(p[1]), p_edge_syms)]   
    return vert_ps, edge_ps
end

# Converts the values for initial condition/vertex parameters to the correct form:
# Map from symbolics to to vectors of length either 1 (for uniform values) or num_verts.
function vertex_value_map(values, lrs)
    isempty(values) && (return Pair{BasicSymbolic{Real}, Vector{Float64}}[])
    return [entry[1] => vertex_value_form(entry[2], lrs, entry[1]) for entry in values]
end
# Converts the values for a specific component (species/parameter) to the correct vector form.
function vertex_value_form(values, lrs::LatticeReactionSystem, sym)
    (values isa AbstractArray) || (return [values])
    if values isa Vector
        if has_grid_lattice(lrs) && (size(values) == grid_size(lrs))
            vertex_value_form(values, lrs.num_verts, lrs.lattice, sym)
        end
        if (length(values) != lrs.num_verts) 
            error("You have provided ($(length(values))) values for $sym. This is not equal to the number of vertexes ($(lrs.num_verts)).")
        end
        return values
    end
    return vertex_value_form(values, lrs.num_verts, lrs.lattice, sym)
end
# Converts values to correct vector form for a Cartesian grid lattice.
function vertex_value_form(values::AbstractArray, num_verts::Int64, lattice::CartesianGridRej{S,T}, sym) where {S,T}
    if size(values) != lattice.dims
        error("The values for $sym did not have the same format as the lattice. Expected a $(lattice.dims) array, got one of size $(size(values))")
    end
    if (length(values) != num_verts) 
        error("You have provided ($(length(values))) values for $sym. This is not equal to the number of vertexes ($(num_verts)).")
    end
    return [values[flat_idx] for flat_idx in 1:num_verts]
end
# Converts values to correct vector form for a masked grid lattice.
function vertex_value_form(values::AbstractArray, num_verts::Int64, lattice::Array{Bool, T}, sym) where {T}
    if size(values) != size(lattice)
        error("The values for $sym did not have the same format as the lattice. Expected a $(size(lattice)) array, got one of size $(size(values))")
    end
    return_values = Vector{typeof(values[1])}(undef, num_verts)
    cur_idx = 0
    for i = 1:length(lattice)
        lattice[i] || continue
        return_values[cur_idx += 1] = values[i]
    end
    if (length(return_values) != num_verts) 
        error("You have provided ($(length(return_values))) values for $sym. This is not equal to the number of vertexes ($(num_verts)).")
    end
    return return_values
end

# Converts the values for initial condition/vertex parameters to the correct form:
# Map from symbolics to to vectors of length either 1 (for uniform values) or num_verts.
function edge_value_map(values, lrs)
    isempty(values) && (return Pair{BasicSymbolic{Real}, SparseMatrixCSC{Float64, Int64}}[])
    return [entry[1] => edge_value_form(entry[2], lrs, entry[1]) for entry in values]
end
# Converts the values for a specific component (species/parameter) to the correct vector form.
function edge_value_form(values, lrs::LatticeReactionSystem, sym)
    # If a scalar have been given, converts it to a size (1,1) sparse matrix.
    (values isa SparseMatrixCSC) || (return sparse([1], [1], [values]))
    
    # Error checks.
    if nnz(values) != lrs.num_edges
        error("You have provided ($(nnz(values))) values for $sym. This is not equal to the number of edges ($(lrs.num_edges)).")
    end
    if !all(Base.isstored(values, e[1], e[2]) for e in lrs.edge_iterator)
        error("Values was not provided for some edges for edge parameter $sym.")
    end

    return values
end

# Creates a map, taking each species (with transportation) to its transportation rate.
# The species is represented by its index (in species(lrs). 
# If the rate is uniform across all edges, the transportation rate will be a size (1,1) sparse matrix.
# Else, the rate will be a size (num_verts,num_verts) sparse matrix.
# In the first step, computes a map from species symbolics form to value(s).
# Second step converts to map from species index to value(s).
function make_sidxs_to_transrate_map(vert_ps::Vector{Pair{BasicSymbolic{Real},Vector{T}}}, 
                                     edge_ps::Vector{Pair{BasicSymbolic{Real},SparseMatrixCSC{T, Int64}}},
                                     lrs::LatticeReactionSystem) where T
    p_val_dict = Dict(vcat(vert_ps, edge_ps))
    transport_rates_speciesmap = compute_all_transport_rates(p_val_dict, lrs)
    return Pair{Int64,SparseMatrixCSC{T, Int64}}[
        speciesmap(lrs.rs)[spat_rates[1]] => spat_rates[2] for spat_rates in transport_rates_speciesmap
    ]
end

# Computes the transport rates for all species with transportation rates. Output is a map
# taking each species; symbolics form to its transportation rates across all edges.
function compute_all_transport_rates(p_val_dict, lrs::LatticeReactionSystem)
    # For all species with transportation, compute their transportation rate (across all edges). 
    # This is a vector, pairing each species to these rates.
    unsorted_rates = [s => compute_transport_rates(get_transport_rate_law(s, lrs), p_val_dict, lrs) 
                        for s in spatial_species(lrs)] 
    
    # Sorts all the species => rate pairs according to their species index in species(::ReactionSystem).
    return sort(unsorted_rates; by = rate -> findfirst(isequal(rate[1]), species(lrs)))   
end
# For a species, retrieves the symbolic expression for its transportation rate
# (likely only a single parameter, such as `D`, but could be e.g. L*D, where L and D are parameters).
# If there are several transportation reactions for the species, their sum is used.
#function get_transport_rate_law(s::BasicSymbolic{Real}, lrs::LatticeReactionSystem)
function get_transport_rate_law(s, lrs)
    rates = filter(sr -> isequal(s, sr.species), lrs.spatial_reactions)
    return sum(getfield.(rates, :rate))
end
# For the numeric expression describing the rate of transport (likely only a single parameter, e.g. `D`), 
# and the values of all our parameters, computes the transport rate(s).
# If all parameters the rate depend on are uniform all edges, this becomes a length 1 vector.
# Else a vector with each value corresponding to the rate at one specific edge.
#function compute_transport_rates(rate_law::Num, p_val_dict, lrs::LatticeReactionSystem)
function compute_transport_rates(rate_law, p_val_dict, lrs)
    # Finds parameters involved in rate and create a function evaluating the rate law.
    relevant_ps = Symbolics.get_variables(rate_law)
    rate_law_func = drop_expr(@RuntimeGeneratedFunction(build_function(rate_law, relevant_ps...)))

    # If all these parameters are spatially uniform, the rates becomes a size (1,1) sparse matrix.
    # Else, the rates becomes a size (num_verts,num_verts) sparse matrix.
    if all(size(p_val_dict[p]) == (1,1) for p in relevant_ps)  
        relevant_p_vals = [get_edge_value(p_val_dict[p], [1 => 1]) for p in relevant_ps]
        return sparse([1],[1],rate_law_func(relevant_p_vals...))
    else
        transport_rates = spzeros(lrs.num_verts, lrs.num_verts)
        for e in lrs.edge_iterator
            relevant_p_vals = [get_edge_value(p_val_dict[p], e) for p in relevant_ps]
            transport_rates[e...] = rate_law_func(relevant_p_vals...)[1]
        end
        return transport_rates
    end
end

# Produces a dictionary with all parameters' values. Vertex parameters have their values converted to
# a sparse matrix (with one value for each edge, always using the source vertex's value)
function param_dict(vert_ps, edge_ps, lrs)
    return merge(Dict(zip(vertex_parameters(lrs), vert_ps)), Dict(zip(edge_parameters(lrs), edge_ps)))
end

### Accessing Unknown & Parameter Array Values ###


# Converts a vector of vectors to a long vector.
# These are used when the initial condition is converted to a single vector (from vector of vector form).
function expand_component_values(values, num_verts)
    vcat([get_vertex_value.(values, vert) for vert in 1:num_verts]...)
end

# Gets the index in the u array of species s in vertex vert (when their are num_species species).
get_index(vert::Int64, s::Int64, num_species::Int64) = (vert - 1) * num_species + s
# Gets the indexes in the u array of all species in vertex vert (when their are num_species species).
get_indexes(vert::Int64, num_species::Int64) = ((vert - 1) * num_species + 1):(vert * num_species)

# Returns the value of a parameter in an edge. For vertex parameters, uses their values in the source.
function get_edge_value(values::Vector{T}, edge) where {T}
    return (length(values) == 1) ? values[1] : values[edge[1]]
end
function get_edge_value(values::SparseMatrixCSC{T, Int64}, edge) where {T}
    return (size(values) == (1,1)) ? values[1,1] : values[edge[1],edge[2]]
end

# Returns the value of an initial condition of parameter in a vertex.
function get_vertex_value(values::Vector{T}, vert_idx) where {T}
    return (length(values) == 1) ? values[1] : values[vert_idx]
end









# Finds the transport rate of a parameter going from a source vertex to a destination vertex.
function get_transport_rate(transport_rate, edge::Pair{Int64,Int64}, t_rate_idx_types::Bool)
    return t_rate_idx_types ? transport_rate[1,1] : transport_rate[edge[1],edge[2]]
end
# Finds the transportation rate for a specific species and a `LatticeTransportODEf` struct.
function get_transport_rate(trans_s_idx::Int64, f_func::LatticeTransportODEf, edge::Pair{Int64,Int64})
    get_transport_rate(f_func.transport_rates[trans_s_idx][2], edge, f_func.t_rate_idx_types[trans_s_idx])
end



# Updates the internal work_vert_ps vector for a given location.
# To this vector, we write the systems parameter values at a specific vertex.
function update_work_vert_ps!(work_vert_ps, vert_ps, comp, vert_ps_idx_types)
    # Loops through all parameters.
    for (idx,loc_type) in enumerate(vert_ps_idx_types)
        # If the parameter is uniform across the spatial structure, it will have a length-1 value vector
        # (which value we write to the work vector).
        # Else, we extract it value at the specific location.
        work_vert_ps[idx] = (loc_type ? vert_ps[idx][1] : vert_ps[idx][comp])   
    end
end
# Input is always either a LatticeTransportODEf or LatticeTransportODEjac function (which fields we then pass on).
function update_work_vert_ps!(lt_ode_func, vert_ps, comp)
    return update_work_vert_ps!(lt_ode_func.work_vert_ps, vert_ps, comp, lt_ode_func.v_ps_idx_types)
end

# Expands a u0/p information stored in Vector{Vector{}} for to Matrix form
# (currently only used in Spatial Jump systems).
function matrix_expand_component_values(values::Vector{<:Vector}, n)
    reshape(expand_component_values(values, n), length(values), n)
end

# For an expression, computes its values using the provided state and parameter vectors.
# The expression is assumed to be valid in edges (and can have edges parameter components).
# If some component is non-uniform, output is a vector of length equal to the number of vertexes.
# If all components are uniform, the output is a length one vector.
function compute_edge_value(exp, lrs::LatticeReactionSystem, edge_ps)
    # Finds the symbols in the expression. Checks that all correspond to edge parameters.
    relevant_syms = Symbolics.get_variables(exp)
    if !all(any(isequal(sym, p) for p in edge_parameters(lrs)) for sym in relevant_syms) 
        error("An non-edge parameter was encountered in expressions: $exp. Here, only edge parameters are expected.")
    end

    # Creates a Function tha computes the expressions value for a parameter set.
    exp_func = drop_expr(@RuntimeGeneratedFunction(build_function(exp, relevant_syms...)))
    # Creates a dictionary with the value(s) for all edge parameters.
    sym_val_dict = vals_to_dict(edge_parameters(lrs), edge_ps)

    # If all values are uniform, compute value once. Else, do it at all edges.
    if !has_spatial_edge_component(exp, lrs, edge_ps)
        return [exp_func([sym_val_dict[sym][1] for sym in relevant_syms]...)]
    end
    return [exp_func([get_component_value(sym_val_dict[sym], idxE) for sym in relevant_syms]...) 
                                                                            for idxE in 1:lrs.num_edges]
end

# For an expression, computes its values using the provided state and parameter vectors.
# The expression is assumed to be valid in vertexes (and can have vertex parameter and state components).
# If at least one component is non-uniform, output is a vector of length equal to the number of vertexes.
# If all components are uniform, the output is a length one vector.
function compute_vertex_value(exp, lrs::LatticeReactionSystem; u=nothing, vert_ps=nothing)
    # Finds the symbols in the expression. Checks that all correspond to states or vertex parameters.
    relevant_syms = Symbolics.get_variables(exp)
    if any(any(isequal(sym) in edge_parameters(lrs)) for sym in relevant_syms) 
        error("An edge parameter was encountered in expressions: $exp. Here, on vertex-based components are expected.")
    end
    # Creates a Function tha computes the expressions value for a parameter set.
    exp_func = drop_expr(@RuntimeGeneratedFunction(build_function(exp, relevant_syms...)))
    # Creates a dictionary with the value(s) for all edge parameters.
    if !isnothing(u) && !isnothing(vert_ps)
        all_syms = [species(lrs); vertex_parameters(lrs)]
        all_vals = [u; vert_ps]
    elseif !isnothing(u) && isnothing(vert_ps)
        all_syms = species(lrs)
        all_vals = u

    elseif isnothing(u) && !isnothing(vert_ps)
        all_syms = vertex_parameters(lrs)
        all_vals = vert_ps
    else
        error("Either u or vertex_ps have to be provided to has_spatial_vertex_component.")
    end
    sym_val_dict = vals_to_dict(all_syms, all_vals)
    
    # If all values are uniform, compute value once. Else, do it at all edges.
    if !has_spatial_vertex_component(exp, lrs; u, vert_ps)
        return [exp_func([sym_val_dict[sym][1] for sym in relevant_syms]...)]
    end
    return [exp_func([get_component_value(sym_val_dict[sym], idxV) for sym in relevant_syms]...) 
                                                                            for idxV in 1:lrs.num_verts]
end

### System Property Checks ###

# For a Symbolic expression, a LatticeReactionSystem, and a parameter list of the internal format:
# Checks if any edge parameter in the expression have a spatial component (that is, is not uniform).
function has_spatial_edge_component(exp, lrs::LatticeReactionSystem, edge_ps)
    # Finds the edge parameters in the expression. Computes their indexes.
    exp_syms = Symbolics.get_variables(exp)
    exp_edge_ps = filter(sym -> any(isequal(sym), edge_parameters(lrs)), exp_syms)
    p_idxs = [findfirst(isequal(sym, edge_p) for edge_p in edge_parameters(lrs)) for sym in exp_syms]
    # Checks if any of the corresponding value vectors have length != 1 (that is, is not uniform).
    return any(length(edge_ps[p_idx]) != 1 for p_idx in p_idxs)
end

# For a Symbolic expression, a LatticeReactionSystem, and a parameter list of the internal format (vector of vectors):
# Checks if any vertex parameter in the expression have a spatial component (that is, is not uniform).
function has_spatial_vertex_component(exp, lrs::LatticeReactionSystem; u=nothing, vert_ps=nothing)
    # Finds all the symbols in the expression.
    exp_syms = Symbolics.get_variables(exp)

    # If vertex parameter values where given, checks if any of these have non-uniform values.
    if !isnothing(vert_ps)
        exp_vert_ps = filter(sym -> any(isequal(sym), vertex_parameters(lrs)), exp_syms)
        p_idxs = [ModelingToolkit.parameter_index(lrs.rs, sym) for sym in exp_vert_ps]
        any(length(vert_ps[p_idx]) != 1 for p_idx in p_idxs) && return true
    end

    # If states values where given, checks if any of these have non-uniform values.
    if !isnothing(u)
        exp_u = filter(sym -> any(isequal(sym), species(lrs)), exp_syms)
        u_idxs = [ModelingToolkit.variable_index(lrs.rs, sym) for sym in exp_u]
        any(length(u[u_idx]) != 1 for u_idx in u_idxs) && return true
    end

    return false
end