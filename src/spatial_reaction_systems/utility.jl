### Processes Input u0 & p ###

# Defines _symbol_to_var, but where the system is a LRS. Required to make symmapt_to_varmap to work.
function _symbol_to_var(lrs::LatticeReactionSystem, sym)
    p_idx = findfirst(sym==p_sym for p_sym in ModelingToolkit.getname.(parameters(lrs)))    # Checks if sym is a parameter.
    isnothing(p_idx) || return parameters(lrs)[p_idx]
    s_idx = findfirst(sym==s_sym for s_sym in ModelingToolkit.getname.(species(lrs)))       # Checks if sym is a species.
    isnothing(s_idx) || return species(lrs)[s_idx]
    error("Could not find property parameter/species $sym in lattice reaction system.")
end

# From u0 input, extracts their values and store them in the internal format (a vector  on the form [species 1 at vertex 1, species 2 at vertex 1, ..., species 1 at vertex 2, ...]).
function lattice_process_u0(u0_in, u0_syms, num_verts)
    u0 = lattice_process_input(u0_in, u0_syms, num_verts)   # u0 values can be given in various forms. This converts it to a Vector{Vector{}} form (Vector containing one vector for each species. Second vector has one value if species uniform across lattice, else one value for each vertex).
    check_vector_lengths(u0, length(u0_syms), num_verts)    # Perform various error checks on the (by the user provided) initial conditions.         
    expand_component_values(u0, num_verts)                  # Converts the Vector{Vector{}} format to a single Vector (with one values for each species and vertex).
end

# From p input, splits it into diffusion parameters and compartment parameters, and store these in the desired internal format.
function lattice_process_p(p_in, p_vertex_syms, p_edge_syms, lrs::LatticeReactionSystem)
    vert_ps_in, edge_ps_in = split_parameters(p_in, p_vertex_syms, p_edge_syms)     # If the user provided parameters as a single map (mixing vertex and edge parameters) , these are split into two separate vectors.
    vert_ps = lattice_process_input(vert_ps_in, p_vertex_syms, lrs.num_verts)       # Parameter values can be given in various forms. This converts it to the Vector{Vector{}} form.
    edge_ps = lattice_process_input(edge_ps_in, p_edge_syms, lrs.num_edges)         # Parameter values can be given in various forms. This converts it to the Vector{Vector{}} form.
    lrs.init_digraph || duplicate_trans_params!(edge_ps, lrs)                        # If the lattice was defined as an undirected graph (with N edges), but then provides N/2 values for some edge parameter, we presume they want to expand that parameters value so it has the same value in both directions. 
    check_vector_lengths(vert_ps, length(p_vertex_syms), lrs.num_verts)             # Perform various error checks on the (by the user provided) vertex parameters.   
    check_vector_lengths(edge_ps, length(p_edge_syms), lrs.num_edges)               # Perform various error checks on the (by the user provided) edge parameters. 
    return vert_ps, edge_ps
end

# Splits parameters into those for the vertexes and those for the edges.
split_parameters(ps::Tuple{<:Any, <:Any}, args...) = ps     # If they are already split, return that.
function split_parameters(ps::Vector{<:Number}, args...)    # Providing parameters to a spatial reaction system as a single vector of values (e.g. [1.0, 4.0, 0.1]) is not allowed, Either use tuple (e.g. ([1.0, 4.0], [0.1])) or map format (e.g. [A => 1.0, B => 4.0, D => 0.1]). 
    error("When providing parameters for a spatial system as a single vector, the paired form (e.g :D =>1.0) must be used.")
end
function split_parameters(ps::Vector{<: Pair}, p_vertex_syms::Vector, p_edge_syms::Vector)  # Splitting is only done for Vectors of Pairs (where the first value is a Symbols, and the second a value).
    vert_ps_in = [p for p in ps if any(isequal(p[1]), p_vertex_syms)]   # The input vertex parameters.   
    edge_ps_in = [p for p in ps if any(isequal(p[1]), p_edge_syms)]     # The input edge parameters.
    (sum(length.([vert_ps_in, edge_ps_in])) != length(ps)) && error("These input parameters are not recognised: $(setdiff(first.(ps), vcat(first.([vert_ps_in, edge_ps_in]))))") # Error check, in case some input parameters where neither recognised as vertex or edge parameters.
    return vert_ps_in, edge_ps_in
end

# Input is allowed on the following forms (after potential Symbol maps have been converted to Symbolic maps:
    # - A vector of values, where the i'th value corresponds to the value of the i'th: initial condition value (for u0_in), vertex parameter value (for vert_ps_in), or edge parameter value (for edge_ps_in).
    # - A vector of vectors of values. The same as previously, but here the species/parameter can have different values across the spatial structure.
    # - A map of symbols to their values. These can either be a single value (if uniform across the spatial structure) or a vector (with different values for each vertex/edge). These can be mixed (e.g. [X => 1.0, Y => [1.0, 2.0, 3.0, 4.0]] is allowed).
    # - A matrix. E.g. for initial conditions you can have a num_species * num_vertex matrix, indicating the value of each species at each vertex.
# The lattice_process_input function takes input initial conditions/vertex parameters/edge parameters of whichever form the user have used, and converts them to the Vector{Vector{}} form used internally. E.g. for parameters the top vector contain one vector for each parameter (using the same order as in parameters(::ReactionSystem)). If a parameter is uniformly-values across the spatial structure, its vector has a single value. Else, it has a number of values corresponding to the number of vertexes/edges (for edge/vertex parameters). Initial conditions works similarly.

# If the input is given in a map form, the vector needs sorting and the first value removed. The creates a Vector{Vector{Value}} or Vector{value} form, which is then again sent to lattice_process_input for reprocessing.
function lattice_process_input(input::Vector{<:Pair}, syms::Vector{BasicSymbolic{Real}}, args...)
    isempty(setdiff(first.(input), syms)) || error("Some input symbols are not recognised: $(setdiff(first.(input), syms)).")
    sorted_input = sort(input; by = p -> findfirst(isequal(p[1]), syms))
    return lattice_process_input(last.(sorted_input), syms, args...)
end
# If the input is a matrix: Processes the input and gives it in a form where it is a vector of vectors (some of which may have a single value). Sends it back to lattice_process_input for reprocessing. 
function lattice_process_input(input::Matrix{<:Number}, args...)
    lattice_process_input([vec(input[i, :]) for i in 1:size(input, 1)], args...)
end
function lattice_process_input(input::Array{<:Number, 3}, args...)
    error("3 dimensional array parameter input currently not supported.")    # Supposedly we want to support this type of input at some point.
end
# If the input is a Vector containing both vectors and single values, converts it to the Vector{<:Vector} form. Technically this last lattice_process_input is probably not needed.
function lattice_process_input(input::Vector{<:Any}, args...)
    isempty(input) ? Vector{Vector{Float64}}() :
    lattice_process_input([(val isa Vector{<:Number}) ? val : [val] for val in input],
                          args...)
end
# If the input is of the correct form already, return it.
lattice_process_input(input::Vector{<:Vector}, syms::Vector{BasicSymbolic{Real}}, n::Int64) = input

# Checks that a value vector have the right length, as well as that of all its sub vectors. Error check if e.g. the user does not provide values for all species/parameters, or for one, provides a vector of values, but that has the wrong length (e.g providing 7 values for one species, but there are 8 vertexes). 
function check_vector_lengths(input::Vector{<:Vector}, n_syms, n_locations)
    (length(input)==n_syms) || error("Missing values for some initial conditions/parameters. Expected $n_syms values, got $(length(input)).")
    isempty(setdiff(unique(length.(input)), [1, n_locations])) || error("Some inputs where given values of inappropriate length.")
end

# For transport parameters, if the lattice was given as an undirected graph of size n, this is converted to a directed graph of size 2n.
# If transport parameters are given with n values, we want to use the same value for both directions.
# Since the order of edges in the new graph is non-trivial, this function distributes the n input values to a 2n length vector, putting the correct value in each position.
function duplicate_trans_params!(edge_ps::Vector{Vector{Float64}}, lrs::LatticeReactionSystem)
    cum_adjacency_counts = [0;cumsum(length.(lrs.lattice.fadjlist[1:end-1]))]
    for idx in 1:length(edge_ps)    # Loops through all edge parameter.
        (2length(edge_ps[idx]) == lrs.num_edges) || continue # If the edge parameter already have one value for each directed edge, it is fine and we can continue.
        
        # This entire thing depends on the fact that, in the edges(lattice) iterator, the edges are sorted by (1) Their source node, (2) Their destination node.
        new_vals = Vector{Float64}(undef,lrs.num_edges)     # A vector where we will put the edge parameters new values. Has the correct length (the number of directed edges in the lattice).
        original_edge_count = 0                             # As we loop through the edges of the di-graph, this keeps track of each edge's index in the original graph. 
        for edge in edges(lrs.lattice)                      # For each edge.
            (edge.src < edge.dst) ? (original_edge_count += 1) : continue    # The digraph conversion only adds edges so that src > dst.
            idx_fwd = cum_adjacency_counts[edge.src] + findfirst(isequal(edge.dst),lrs.lattice.fadjlist[edge.src]) # For original edge i -> j, finds the index of i -> j in DiGraph.
            idx_bwd = cum_adjacency_counts[edge.dst] + findfirst(isequal(edge.src),lrs.lattice.fadjlist[edge.dst]) # For original edge i -> j, finds the index of j -> i in DiGraph.
            new_vals[idx_fwd] = edge_ps[idx][original_edge_count]
            new_vals[idx_bwd] = edge_ps[idx][original_edge_count]
        end
        edge_ps[idx] = new_vals # Replaces the edge parameters values with the updated value vector.
    end
end

# For a set of input values on the given forms, and their symbolics, convert into a dictionary.
vals_to_dict(syms::Vector, vals::Vector{<:Vector}) = Dict(zip(syms, vals))
# Produces a dictionary with all parameter values.
function param_dict(vert_ps, edge_ps, lrs)
    merge(vals_to_dict(vertex_parameters(lrs), vert_ps),
          vals_to_dict(edge_parameters(lrs), edge_ps))
end

# Computes the transport rates and stores them in a desired format (a Dictionary from species index to rates across all edges).
function compute_all_transport_rates(vert_ps::Vector{Vector{Float64}}, edge_ps::Vector{Vector{Float64}}, lrs::LatticeReactionSystem)
    param_value_dict = param_dict(vert_ps, edge_ps, lrs)    # Creates a dict, allowing us to access the values of wll parameters.
    unsorted_rates = [s => Symbolics.value.(compute_transport_rates(get_transport_rate_law(s, lrs), param_value_dict, lrs.num_edges)) for s in spatial_species(lrs)] # For all species with transportation, compute their transportation rate (across all edges). This is a vector, pairing each species to these rates.
    return sort(unsorted_rates; by=rate -> findfirst(isequal(rate[1]), species(lrs)))   # Sorts all the species => rate pairs according to their species index in species(::ReactionSystem),
end
# For a species, retrieves the symbolic expression for its transportation rate (likely only a single parameter, such as `D`, but could be e.g. L*D, where L and D are parameters).
function get_transport_rate_law(s::BasicSymbolic{Real}, lrs::LatticeReactionSystem)
    rates = filter(sr -> isequal(s, sr.species), lrs.spatial_reactions)
    (length(rates) > 1) && error("Species $s have more than one diffusion reaction.")    # We could allows several and simply sum them though, easy change.
    return rates[1].rate
end
# For the numeric expression describing the rate of transport (likely only a single parameter, such as `D`), and the values of all our parameters, computes the transport rate(s). If all parameters the rate depend on are uniform all edges, this becomes a length 1 vector. Else a vector with each value corresponding to the rate at one specific edge.
function compute_transport_rates(rate_law::Num,
                                 param_value_dict::Dict{SymbolicUtils.BasicSymbolic{Real}, Vector{Float64}}, num_edges::Int64)
    relevant_parameters = Symbolics.get_variables(rate_law)                 # Extracts the parameters that this rate depend on.
    if all(length(param_value_dict[P]) == 1 for P in relevant_parameters)   # If all these parameters are spatially uniform.
        return [
            substitute(rate_law,
                       Dict(p => param_value_dict[p][1] for p in relevant_parameters)), # Computes a length 1 vector with this value.
        ]
    end
    return [substitute(rate_law,                                            # If at least on parameter the rate depends on have a value varying across all edges, we have to compute one rate value for each edge.
                       Dict(p => get_component_value(param_value_dict[p], idxE)
                            for p in relevant_parameters)) for idxE in 1:num_edges]
end

### Accessing State & Parameter Array Values ###

# Gets the index in the u array of species s in vertex vert (when their are num_species species).
get_index(vert::Int64, s::Int64, num_species::Int64) = (vert - 1) * num_species + s
# Gets the indexes in the u array of all species in vertex vert (when their are num_species species).
get_indexes(vert::Int64, num_species::Int64) = ((vert - 1) * num_species + 1):(vert * num_species)

# We have many vectors of length 1 or n, for which we want to get value idx (or the one value, if length is 1), this function gets that. Here:
    # - values is the vector with the values of teh component across all locations (where the internal vectors may or may not be of size 1).
    # - component_idx is the initial condition species/vertex parameter/edge parameters's index. This is predominantly used for parameters, for initial conditions, it is only used once (at initialisation) to re-process the input vector.
    # - location_idx is the index of the vertex or edge for which we wish to access a initial condition or parameter values
# The first two function takes the full value vector, and call the function of at the components specific index.
function get_component_value(values::Vector{<:Vector}, component_idx::Int64,
                             location_idx::Int64)
    get_component_value(values[component_idx], location_idx)
end
function get_component_value(values::Vector{<:Vector}, component_idx::Int64,
                             location_idx::Int64, location_types::Vector{Bool})     # Sometimes we have pre-computed, for each component, whether it's vector is length 1 or not. This is stored in location_types.
    get_component_value(values[component_idx], location_idx, location_types[component_idx])
end
# For a components value (which is a vector of either length 1 or some other length), retrieves its value.
function get_component_value(values::Vector{<:Number}, location_idx::Int64)
    get_component_value(values, location_idx, length(values) == 1)
end
function get_component_value(values::Vector{<:Number}, location_idx::Int64,
                             location_type::Bool)
    location_type ? values[1] : values[location_idx]                                # Again, the location type (length of teh value vector) may be pre-computed.
end

# Converts a vector of vectors to a long vector. These are used when the initial condition is converted to a single vector (from vector of vector form).
function expand_component_values(values::Vector{<:Vector}, n)
    vcat([get_component_value.(values, comp) for comp in 1:n]...)
end
function expand_component_values(values::Vector{<:Vector}, n, location_types::Vector{Bool})
    vcat([get_component_value.(values, comp, location_types) for comp in 1:n]...)
end

# Creates a view of the vert_ps vector at a given location. Provides a work vector to which the converted vector is written. 
function view_vert_ps_vector!(work_vert_ps, vert_ps, comp, enumerated_vert_ps_idx_types)
    for (idx,loc_type) in enumerated_vert_ps_idx_types      # Loops through all parameters.
        work_vert_ps[idx] = (loc_type ? vert_ps[idx][1] : vert_ps[idx][comp])   # If the parameter is uniform across the spatial structure, it will have a length-1 value vector (which value we write to the work vector). Else, we extract it value at the specific location.
    end
    return work_vert_ps
end

# Expands a u0/p information stored in Vector{Vector{}} for to Matrix form (currently only used in Spatial Jump systems).
function matrix_expand_component_values(values::Vector{<:Vector}, n)
    reshape(expand_component_values(values, n), length(values), n)
end
