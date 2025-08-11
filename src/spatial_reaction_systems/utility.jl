### Processes Input u0 & p ###

# Defines _symbol_to_var, but where the system is a LRS. Required to make symmapt_to_varmap to work.
function _symbol_to_var(lrs::LatticeReactionSystem, sym)
    # Checks if sym is a parameter.
    p_idx = findfirst(sym == p_sym for p_sym in ModelingToolkit.getname.(parameters(lrs)))
    isnothing(p_idx) || return parameters(lrs)[p_idx]

    # Checks if sym is a species.
    s_idx = findfirst(sym == s_sym for s_sym in ModelingToolkit.getname.(species(lrs)))
    isnothing(s_idx) || return species(lrs)[s_idx]

    error("Could not find property parameter/species $sym in lattice reaction system.")
end

# From u0 input, extract their values and store them in the internal format.
# Internal format: a vector on the form [spec 1 at vert 1, spec 2 at vert 1, ..., spec 1 at vert 2, ...]).
function lattice_process_u0(u0_in, u0_syms::Vector, lrs::LatticeReactionSystem)
    # u0 values can be given in various forms. This converts it to a Vector{Pair{Symbolics,...}} form.
    # Top-level vector: Maps each species to its value(s).
    u0 = lattice_process_input(u0_in, u0_syms)

    # Species' initial condition values can be given in different forms (also depending on the lattice).
    # This converts each species's values to a Vector. In it, for species with uniform initial conditions,
    # it holds that value only. For spatially heterogeneous initial conditions,
    # the vector has the same length as the number of vertices (storing one value for each).
    u0 = vertex_value_map(u0, lrs)

    # Converts the initial condition to a single Vector (with one value for each species and vertex).
    return expand_component_values([entry[2] for entry in u0], num_verts(lrs))
end

# From a parameter input, split it into vertex parameters and edge parameters.
# Store these in the desired internal format.
function lattice_process_p(ps_in, ps_vertex_syms::Vector,
        ps_edge_syms::Vector, lrs::LatticeReactionSystem)
    # p values can be given in various forms. This converts it to a Vector{Pair{Symbolics,...}} form.
    # Top-level vector: Maps each parameter to its value(s).
    # Second-level: Contains either a vector (vertex parameters) or a sparse matrix (edge parameters).
    # For uniform parameters these have size 1/(1,1). Else, they have size num_verts/(num_verts,num_verts).
    ps = lattice_process_input(ps_in, [ps_vertex_syms; ps_edge_syms])

    # Split the parameter vector into one for vertex parameters and one for edge parameters.
    # Next, convert their values to the correct form (vectors for vert_ps and sparse matrices for edge_ps).
    vert_ps, edge_ps = split_parameters(ps, ps_vertex_syms, ps_edge_syms)
    vert_ps = vertex_value_map(vert_ps, lrs)
    edge_ps = edge_value_map(edge_ps, lrs)

    return vert_ps, edge_ps
end

# The input (parameters or initial conditions) may either be a dictionary (symbolics to value(s).)
# or a map (in vector or tuple form) from symbolics to value(s). This converts the input to a
# (Vector) map from symbolics to value(s), where the entries have the same order as `syms`.
function lattice_process_input(input::Dict{<:Any, T}, syms::Vector) where {T}
    # Error checks
    if !isempty(setdiff(keys(input), syms))
        throw(ArgumentError("You have provided values for the following unrecognised parameters/initial conditions: $(setdiff(keys(input), syms))."))
    end
    if !isempty(setdiff(syms, keys(input)))
        throw(ArgumentError("You have not provided values for the following parameters/initial conditions: $(setdiff(syms, keys(input)))."))
    end

    return [sym => input[sym] for sym in syms]
end
function lattice_process_input(input, syms::Vector)
    if ((input isa Vector) || (input isa Tuple)) && all(entry isa Pair for entry in input)
        return lattice_process_input(Dict(input), syms)
    end
    throw(ArgumentError("Input parameters/initial conditions have the wrong format ($(typeof(input))). These should either be a Dictionary, or a Tuple or a Vector (where each entry is a Pair taking a parameter/species to its value)."))
end

# Splits parameters into vertex and edge parameters.
function split_parameters(ps, p_vertex_syms::Vector, p_edge_syms::Vector)
    vert_ps = [p for p in ps if any(isequal(p[1]), p_vertex_syms)]
    edge_ps = [p for p in ps if any(isequal(p[1]), p_edge_syms)]
    return vert_ps, edge_ps
end

# Converts the values for the initial conditions/vertex parameters to the correct form:
# A map vector from symbolics to vectors of either length 1 (for uniform values) or num_verts.
function vertex_value_map(values, lrs::LatticeReactionSystem)
    isempty(values) && (return Pair{BasicSymbolic{Real}, Vector{Float64}}[])
    return [entry[1] => vertex_value_form(entry[2], lrs, entry[1]) for entry in values]
end

# Converts the values for an individual species/vertex parameter to its correct vector form.
function vertex_value_form(values, lrs::LatticeReactionSystem, sym::BasicSymbolic)
    # If the value is a scalar (i.e. uniform across the lattice), return it in vector form.
    (values isa AbstractArray) || (return [values])

    # If the value is a vector (something all three lattice types accept).
    if values isa Vector
        # For the case where we have a 1d (Cartesian or masked) grid, and the vector's values
        # correspond to individual grid points.
        if has_grid_lattice(lrs) && (size(values) == grid_size(lrs))
            return vertex_value_form(values, num_verts(lrs), lattice(lrs), sym)
        end

        # For the case where the i'th value of the vector corresponds to the value in the i'th vertex.
        # This is the only (non-uniform) case possible for graph grids.
        if (length(values) != num_verts(lrs))
            throw(ArgumentError("You have provided ($(length(values))) values for $sym. This is not equal to the number of vertices ($(num_verts(lrs)))."))
        end
        return values
    end

    # (2d and 3d) Cartesian and masked grids can take non-vector, non-scalar, values input.
    return vertex_value_form(values, num_verts(lrs), lattice(lrs), sym)
end

# Converts values to the correct vector form for a Cartesian grid lattice.
function vertex_value_form(values::AbstractArray, num_verts::Int64,
        lattice::CartesianGridRej{N, T}, sym::BasicSymbolic) where {N, T}
    if size(values) != lattice.dims
        throw(ArgumentError("The values for $sym did not have the same format as the lattice. Expected a $(lattice.dims) array, got one of size $(size(values))"))
    end
    if (length(values) != num_verts)
        throw(ArgumentError("You have provided ($(length(values))) values for $sym. This is not equal to the number of vertices ($(num_verts))."))
    end
    return [values[flat_idx] for flat_idx in 1:num_verts]
end

# Converts values to the correct vector form for a masked grid lattice.
function vertex_value_form(values::AbstractArray, num_verts::Int64,
        lattice::Array{Bool, T}, sym::BasicSymbolic) where {T}
    if size(values) != size(lattice)
        throw(ArgumentError("The values for $sym did not have the same format as the lattice. Expected a $(size(lattice)) array, got one of size $(size(values))"))
    end

    # Pre-declares a vector with the values in each vertex (return_values).
    # Loops through the lattice and the values, adding these to the return_values.
    return_values = Vector{typeof(values[1])}(undef, num_verts)
    cur_idx = 0
    for (idx, val) in enumerate(values)
        lattice[idx] || continue
        return_values[cur_idx += 1] = val
    end

    # Checks that the correct number of values was provided, and returns the values.
    if (length(return_values) != num_verts)
        throw(ArgumentError("You have provided ($(length(return_values))) values for $sym. This is not equal to the number of vertices ($(num_verts))."))
    end
    return return_values
end

# Converts the values for the edge parameters to the correct form:
# A map vector from symbolics to sparse matrices of size either (1,1) or (num_verts,num_verts).
function edge_value_map(values, lrs::LatticeReactionSystem)
    isempty(values) && (return Pair{BasicSymbolic{Real}, SparseMatrixCSC{Float64, Int64}}[])
    return [entry[1] => edge_value_form(entry[2], lrs, entry[1]) for entry in values]
end

# Converts the values for an individual edge parameter to its correct sparse matrix form.
function edge_value_form(values, lrs::LatticeReactionSystem, sym)
    # If the value is a scalar (i.e. uniform across the lattice), return it in sparse matrix form.
    (values isa SparseMatrixCSC) || (return sparse([1], [1], [values]))

    # Error checks.
    if nnz(values) != num_edges(lrs)
        throw(ArgumentError("You have provided ($(nnz(values))) values for $sym. This is not equal to the number of edges ($(num_edges(lrs)))."))
    end
    if !all(Base.isstored(values, e[1], e[2]) for e in edge_iterator(lrs))
        throw(ArgumentError("Values was not provided for some edges for edge parameter $sym."))
    end

    # Unlike initial conditions/vertex parameters, (unless uniform) edge parameters' values  are
    # always provided in the same (sparse matrix) form.
    return values
end

# Creates a map, taking each species (with transportation) to its transportation rate.
# The species is represented by its index (in species(lrs).
# If the rate is uniform across all edges, the transportation rate will be a size (1,1) sparse matrix.
# Else, the rate will be a size (num_verts,num_verts) sparse matrix.
function make_sidxs_to_transrate_map(vert_ps::Vector{Pair{R, Vector{T}}},
        edge_ps::Vector{Pair{S, SparseMatrixCSC{T, Int64}}},
        lrs::LatticeReactionSystem) where {R, S, T}
    # Creates a dictionary with each parameter's value(s).
    p_val_dict = Dict(vcat(vert_ps, edge_ps))

    # First, compute a map from species in their symbolics form to their values.
    # Next, convert to map from species index to values.
    transport_rates_speciesmap = compute_all_transport_rates(p_val_dict, lrs)
    return Pair{Int64, SparseMatrixCSC{T, Int64}}[speciesmap(reactionsystem(lrs))[spat_rates[1]] => spat_rates[2]
                                                  for spat_rates in
                                                      transport_rates_speciesmap]
end

# Computes the transport rates for all species with transportation rates. Output is a map
# taking each species' symbolics form to its transportation rates across all edges.
function compute_all_transport_rates(p_val_dict, lrs::LatticeReactionSystem)
    # For all species with transportation, compute their transportation rate (across all edges).
    # This is a vector, pairing each species to these rates.
    unsorted_rates = [s => compute_transport_rates(s, p_val_dict, lrs)
                      for s in spatial_species(lrs)]

    # Sorts all the species => rate pairs according to their species index in species(lrs).
    return sort(unsorted_rates; by = rate -> findfirst(isequal(rate[1]), species(lrs)))
end

# For the expression describing the rate of transport (likely only a single parameter, e.g. `D`),
# and the values of all our parameters, compute the transport rate(s).
# If all parameters that the rate depends on are uniform across all edges, this becomes a length-1 vector.
# Else it becomes a vector where each value corresponds to the rate at one specific edge.
function compute_transport_rates(s::BasicSymbolic, p_val_dict, lrs::LatticeReactionSystem)
    # Find parameters involved in the rate and create a function evaluating the rate law.
    rate_law = get_transport_rate_law(s, lrs)
    relevant_ps = Symbolics.get_variables(rate_law)
    rate_law_func = drop_expr(@RuntimeGeneratedFunction(build_function(
        rate_law, relevant_ps...)))

    # If all these parameters are spatially uniform, the rates become a size (1,1) sparse matrix.
    # Else, the rates become a size (num_verts,num_verts) sparse matrix.
    if all(size(p_val_dict[p]) == (1, 1) for p in relevant_ps)
        relevant_p_vals = [get_edge_value(p_val_dict[p], 1 => 1) for p in relevant_ps]
        return sparse([1], [1], rate_law_func(relevant_p_vals...))
    else
        transport_rates = spzeros(num_verts(lrs), num_verts(lrs))
        for e in edge_iterator(lrs)
            relevant_p_vals = [get_edge_value(p_val_dict[p], e) for p in relevant_ps]
            transport_rates[e...] = rate_law_func(relevant_p_vals...)[1]
        end
        return transport_rates
    end
end

# For a species, retrieve the symbolic expression for its transportation rate
# (likely only a single parameter, such as `D`, but could be e.g. L*D, where L and D are parameters).
# If there are several transportation reactions for the species, their sum is used.
function get_transport_rate_law(s::BasicSymbolic, lrs::LatticeReactionSystem)
    rates = filter(sr -> isequal(s, sr.species), spatial_reactions(lrs))
    return sum(getfield.(rates, :rate))
end

### Accessing Unknown & Parameter Array Values ###

# Converts a vector of vectors to a single, long, vector.
# These are used when the initial condition is converted to a single vector (from vector of vector form).
function expand_component_values(values::Vector{<:Vector}, num_verts::Int64)
    vcat([get_vertex_value.(values, vert) for vert in 1:num_verts]...)
end

# Gets the index in the u array of species s in vertex vert (when there are num_species species).
get_index(vert::Int64, s::Int64, num_species::Int64) = (vert - 1) * num_species + s
# Gets the indices in the u array of all species in vertex vert (when there are num_species species).
function get_indexes(vert::Int64, num_species::Int64)
    return ((vert - 1) * num_species + 1):(vert * num_species)
end

# Returns the value of a parameter in an edge. For vertex parameters, use their values in the source.
function get_edge_value(values::Vector{T}, edge::Pair{Int64, Int64}) where {T}
    return (length(values) == 1) ? values[1] : values[edge[1]]
end
function get_edge_value(values::SparseMatrixCSC{T, Int64},
        edge::Pair{Int64, Int64}) where {T}
    return (size(values) == (1, 1)) ? values[1, 1] : values[edge[1], edge[2]]
end

# Returns the value of an initial condition of vertex parameter in a vertex.
function get_vertex_value(values::Vector{T}, vert_idx::Int64) where {T}
    return (length(values) == 1) ? values[1] : values[vert_idx]
end

# Finds the transport rate of a parameter along a specific edge.
function get_transport_rate(transport_rate::SparseMatrixCSC{T, Int64},
        edge::Pair{Int64, Int64}, t_rate_idx_types::Bool) where {T}
    return t_rate_idx_types ? transport_rate[1, 1] : transport_rate[edge[1], edge[2]]
end

# For a `LatticeTransportODEFunction`, update its stored parameters (in `mtk_ps`) so that they
# the heterogeneous parameters' values correspond to the values in the specified vertex.
function update_mtk_ps!(lt_ofun::LatticeTransportODEFunction, all_ps::Vector{T},
        vert::Int64) where {T}
    for (setp, idx) in zip(lt_ofun.p_setters, lt_ofun.heterogeneous_vert_p_idxs)
        setp(lt_ofun.mtk_ps, all_ps[idx][vert])
    end
end

# For an expression, compute its values using the provided state and parameter vectors.
# The expression is assumed to be valid in vertices (and can have vertex parameter and state components).
# If at least one component is non-uniform, output is a vector of length equal to the number of vertices.
# If all components are uniform, the output is a length one vector.
function compute_vertex_value(exp, lrs::LatticeReactionSystem; u = [], ps = [])
    # Finds the symbols in the expression. Checks that all correspond to unknowns or vertex parameters.
    relevant_syms = Symbolics.get_variables(exp)
    if any(any(isequal(sym) in edge_parameters(lrs)) for sym in relevant_syms)
        throw(ArgumentError("An edge parameter was encountered in expressions: $exp. Here, only vertex-based components are expected."))
    end

    # Creates a Function that computes the expression value for a parameter set.
    exp_func = drop_expr(@RuntimeGeneratedFunction(build_function(exp, relevant_syms...)))

    # Creates a dictionary with the value(s) for all edge parameters.
    value_dict = Dict(vcat(u, ps))

    # If all values are uniform, compute value once. Else, do it at all edges.
    if all(length(value_dict[sym]) == 1 for sym in relevant_syms)
        return [exp_func([value_dict[sym][1] for sym in relevant_syms]...)]
    end
    return [exp_func([get_vertex_value(value_dict[sym], vert_idx) for sym in relevant_syms]...)
            for vert_idx in 1:num_verts(lrs)]
end

### System Property Checks ###

# For a Symbolic expression, and a parameter set, check if any relevant parameters have a
# spatial component. Filters out any parameters that are edge parameters.
function has_spatial_vertex_component(exp, ps)
    relevant_syms = Symbolics.get_variables(exp)
    value_dict = Dict(filter(p -> p[2] isa Vector, ps))
    return any(length(value_dict[sym]) > 1 for sym in relevant_syms)
end

### Utilities previously in ModelingToolkit.jl ###

function todict(d)
    eltype(d) <: Pair || throw(ArgumentError("The variable-value mapping must be a Dict."))
    d isa Dict ? d : Dict(d)
end

_merge(d1, d2) = merge(todict(d1), todict(d2))

MT.refreshed_metadata(::Nothing) = MT.MetadataT() # FIXME: Type piracy
