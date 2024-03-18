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
function lattice_process_u0(u0_in, u0_syms, num_verts)
    # u0 values can be given in various forms. This converts it to a Vector{Vector{}} form.
    # Top-level vector: Contains one vector for each species. 
    # Second-level vector: contain one value if species uniform across lattice, else one value for each vertex).
    u0 = lattice_process_input(u0_in, u0_syms, num_verts)   

    # Perform various error checks on the (by the user provided) initial conditions.
    check_vector_lengths(u0, length(u0_syms), num_verts)             

    # Converts the Vector{Vector{}} format to a single Vector (with one values for each species and vertex).
    expand_component_values(u0, num_verts)                  
end

# From p input, splits it into diffusion parameters and compartment parameters.
# Store these in the desired internal format.
function lattice_process_p(p_in, p_vertex_syms, p_edge_syms, lrs::LatticeReactionSystem)
    # If the user provided parameters as a single map (mixing vertex and edge parameters):
    # Split into two separate vectors.
    vert_ps_in, edge_ps_in = split_parameters(p_in, p_vertex_syms, p_edge_syms)     

    # Parameter values can be given in various forms. This converts it to the Vector{Vector{}} form.
    vert_ps = lattice_process_input(vert_ps_in, p_vertex_syms, lrs.num_verts)       

    # Parameter values can be given in various forms. This converts it to the Vector{Vector{}} form.
    edge_ps = lattice_process_input(edge_ps_in, p_edge_syms, lrs.num_edges)         

    # If the lattice defined as (N edge) undirected graph, and we provides N/2 values for some edge parameter:
    # Presume they want to expand that parameters value so it has the same value in both directions. 
    lrs.init_digraph || duplicate_trans_params!(edge_ps, lrs)   
    
    # Perform various error checks on the (by the user provided) vertex and edge parameters.
    check_vector_lengths(vert_ps, length(p_vertex_syms), lrs.num_verts)                
    check_vector_lengths(edge_ps, length(p_edge_syms), lrs.num_edges) 

    return vert_ps, edge_ps
end

# Splits parameters into those for the vertexes and those for the edges.

# If they are already split, return that.
split_parameters(ps::Tuple{<:Any, <:Any}, args...) = ps    
# Providing parameters to a spatial reaction system as a single vector of values (e.g. [1.0, 4.0, 0.1]) is not allowed.
# Either use tuple (e.g. ([1.0, 4.0], [0.1])) or map format (e.g. [A => 1.0, B => 4.0, D => 0.1]).  
function split_parameters(ps::Vector{<:Number}, args...)    
    error("When providing parameters for a spatial system as a single vector, the paired form (e.g :D =>1.0) must be used.")
end
# Splitting is only done for Vectors of Pairs (where the first value is a Symbols, and the second a value).
function split_parameters(ps::Vector{<: Pair}, p_vertex_syms::Vector, p_edge_syms::Vector)  
    vert_ps_in = [p for p in ps if any(isequal(p[1]), p_vertex_syms)]     
    edge_ps_in = [p for p in ps if any(isequal(p[1]), p_edge_syms)]     

    # Error check, in case some input parameters where neither recognised as vertex or edge parameters.
    if (sum(length.([vert_ps_in, edge_ps_in])) != length(ps))
        error("These input parameters are not recognised: $(setdiff(first.(ps), vcat(first.([vert_ps_in, edge_ps_in]))))") 
    end
    
    return vert_ps_in, edge_ps_in
end

# Input may have the following forms (after potential Symbol maps to Symbolic maps conversions):
    # - A vector of values, where the i'th value corresponds to the value of the i'th
    #   initial condition value (for u0_in), vertex parameter value (for vert_ps_in), or edge parameter value (for edge_ps_in).
    # - A vector of vectors of values. The same as previously, 
    #   but here the species/parameter can have different values across the spatial structure.
    # - A map of Symbols to values. These can either be a single value (if uniform across the spatial structure)
    #   or a vector (with different values for each vertex/edge). 
    #   These can be mixed (e.g. [X => 1.0, Y => [1.0, 2.0, 3.0, 4.0]] is allowed).
    # - A matrix. E.g. for initial conditions you can have a num_species * num_vertex matrix,
    #   indicating the value of each species at each vertex.

# The lattice_process_input function takes input initial conditions/vertex parameters/edge parameters
# of whichever form the user have used, and converts them to the Vector{Vector{}} form used internally.
# E.g. for parameters the top-level vector contain one vector for each parameter (same order as in parameters(::ReactionSystem)).
# If a parameter is uniformly-values across the spatial structure, its vector has a single value. 
# Else, it has a number of values corresponding to the number of vertexes/edges (for edge/vertex parameters). 
# Initial conditions works similarly.

# If the input is given in a map form, the vector needs sorting and the first value removed. 
# The creates a Vector{Vector{Value}} or Vector{value} form, which is then again sent to lattice_process_input for reprocessing.
function lattice_process_input(input::Vector{<:Pair}, syms::Vector{BasicSymbolic{Real}}, args...)
    if !isempty(setdiff(first.(input), syms)) 
        error("Some input symbols are not recognised: $(setdiff(first.(input), syms)).")
    end
    sorted_input = sort(input; by = p -> findfirst(isequal(p[1]), syms))
    return lattice_process_input(last.(sorted_input), syms, args...)
end
# If the input is a matrix: Processes the input and gives it in a form where it is a vector of vectors
# (some of which may have a single value). Sends it back to lattice_process_input for reprocessing. 
function lattice_process_input(input::Matrix{<:Number}, args...)
    lattice_process_input([vec(input[i, :]) for i in 1:size(input, 1)], args...)
end
# Possibly we want to support this type of input at some point.
function lattice_process_input(input::Array{<:Number, 3}, args...)
    error("3 dimensional array parameter input currently not supported.")    
end
# If the input is a Vector containing both vectors and single values, converts it to the Vector{<:Vector} form.
# Technically this last lattice_process_input is probably not needed.
function lattice_process_input(input::Vector{<:Any}, args...)
    isempty(input) ? Vector{Vector{Float64}}() :
    lattice_process_input([(val isa Vector{<:Number}) ? val : [val] for val in input],
                          args...)
end
# If the input is of the correct form already, return it.
lattice_process_input(input::Vector{<:Vector}, syms::Vector{BasicSymbolic{Real}}, n::Int64) = input

# Checks that a value vector have the right length, as well as that of all its sub vectors. 
# Error check if e.g. the user does not provide values for all species/parameters,
# or for one: provides a vector of values, but that has the wrong length 
# (e.g providing 7 values for one species, but there are 8 vertexes). 
function check_vector_lengths(input::Vector{<:Vector}, n_syms, n_locations)
    if (length(input)!=n_syms) 
        error("Missing values for some initial conditions/parameters. Expected $n_syms values, got $(length(input)).")
    end
    if !isempty(setdiff(unique(length.(input)), [1, n_locations])) 
        error("Some inputs where given values of inappropriate length.")
    end
end

# For transport parameters, if the lattice was given as an undirected graph of size n:
# this is converted to a directed graph of size 2n.
# If transport parameters are given with n values, we want to use the same value for both directions.
# Since the order of edges in the new graph is non-trivial, this function 
# distributes the n input values to a 2n length vector, putting the correct value in each position.
function duplicate_trans_params!(edge_ps::Vector{Vector{Float64}}, lrs::LatticeReactionSystem)
    cum_adjacency_counts = [0;cumsum(length.(lrs.lattice.fadjlist[1:end-1]))]
    for idx in 1:length(edge_ps)
        #  If the edge parameter already values for each directed edge, we can continue.
        (2length(edge_ps[idx]) == lrs.num_edges) || continue #
        
        # This entire thing depends on the fact that, in the edges(lattice) iterator, the edges are sorted by:
        # (1) Their source node
        # (2) Their destination node.

        # A vector where we will put the edge parameters new values. 
        # Has the correct length (the number of directed edges in the lattice).
        new_vals = Vector{Float64}(undef, lrs.num_edges)     
        # As we loop through the edges of the di-graph, this keeps track of each edge's index in the original graph. 
        original_edge_count = 0                             
        for edge in edges(lrs.lattice)                      # For each edge.
            # The digraph conversion only adds edges so that src > dst.
            (edge.src < edge.dst) ? (original_edge_count += 1) : continue    
            # For original edge i -> j, finds the index of i -> j in DiGraph.
            idx_fwd = cum_adjacency_counts[edge.src] + findfirst(isequal(edge.dst),lrs.lattice.fadjlist[edge.src]) 
            # For original edge i -> j, finds the index of j -> i in DiGraph.
            idx_bwd = cum_adjacency_counts[edge.dst] + findfirst(isequal(edge.src),lrs.lattice.fadjlist[edge.dst]) 
            new_vals[idx_fwd] = edge_ps[idx][original_edge_count]
            new_vals[idx_bwd] = edge_ps[idx][original_edge_count]
        end
        # Replaces the edge parameters values with the updated value vector.
        edge_ps[idx] = new_vals 
    end
end

# For a set of input values on the given forms, and their symbolics, convert into a dictionary.
vals_to_dict(syms::Vector, vals::Vector{<:Vector}) = Dict(zip(syms, vals))
# Produces a dictionary with all parameter values.
function param_dict(vert_ps, edge_ps, lrs)
    merge(vals_to_dict(vertex_parameters(lrs), vert_ps),
          vals_to_dict(edge_parameters(lrs), edge_ps))
end

# Computes the transport rates and stores them in a desired format 
# (a Dictionary from species index to rates across all edges).
function compute_all_transport_rates(vert_ps::Vector{Vector{Float64}}, edge_ps::Vector{Vector{Float64}}, lrs::LatticeReactionSystem)
    # Creates a dict, allowing us to access the values of wll parameters.
    p_val_dict = param_dict(vert_ps, edge_ps, lrs)    

    # For all species with transportation, compute their transportation rate (across all edges). 
    # This is a vector, pairing each species to these rates.
    unsorted_rates = [s => compute_transport_rates(get_transport_rate_law(s, lrs), p_val_dict, lrs.num_edges) 
                        for s in spatial_species(lrs)] 
    
    # Sorts all the species => rate pairs according to their species index in species(::ReactionSystem).
    return sort(unsorted_rates; by=rate -> findfirst(isequal(rate[1]), species(lrs)))   
end
# For a species, retrieves the symbolic expression for its transportation rate
# (likely only a single parameter, such as `D`, but could be e.g. L*D, where L and D are parameters).
# We could allows several transportation reactions for one species and simply sum them though, easy change.
function get_transport_rate_law(s::BasicSymbolic{Real}, lrs::LatticeReactionSystem)
    rates = filter(sr -> isequal(s, sr.species), lrs.spatial_reactions)
    (length(rates) > 1) && error("Species $s have more than one diffusion reaction.")    
    return rates[1].rate
end
# For the numeric expression describing the rate of transport (likely only a single parameter, e.g. `D`), 
# and the values of all our parameters, computes the transport rate(s).
# If all parameters the rate depend on are uniform all edges, this becomes a length 1 vector.
# Else a vector with each value corresponding to the rate at one specific edge.
function compute_transport_rates(rate_law::Num,
                                p_val_dict::Dict{SymbolicUtils.BasicSymbolic{Real}, Vector{Float64}}, num_edges::Int64)
    # Finds parameters involved in rate and create a function evaluating teh rate law.
    relevant_ps = Symbolics.get_variables(rate_law)
    rate_law_func = drop_expr(@RuntimeGeneratedFunction(build_function(rate_law, relevant_ps...)))

    # If all these parameters are spatially uniform. `rates` becomes a vector with 1 value.
    if all(length(p_val_dict[P]) == 1 for P in relevant_ps)  
        return [rate_law_func([p_val_dict[p][1] for p in relevant_ps]...)]
    # If at least on parameter the rate depends on have a value varying across all edges,
    # we have to compute one rate value for each edge.
    else
        return [rate_law_func([get_component_value(p_val_dict[p], idxE) for p in relevant_ps]...) 
                    for idxE in 1:num_edges]
    end
end

# Creates a map, taking each species (with transportation) to its transportation rate.
# The species is represented by its index (in species(lrs). 
# If the rate is uniform across all edges, the vector will be length 1 (with this value),
# else there will be a separate value for each edge.
# Pair{Int64, Vector{T}}[] is required in case vector is empty (otherwise it becomes Any[], causing type error later).
function make_sidxs_to_transrate_map(vert_ps::Vector{Vector{Float64}}, edge_ps::Vector{Vector{T}}, 
                                                                lrs::LatticeReactionSystem) where T
    transport_rates_speciesmap = compute_all_transport_rates(vert_ps, edge_ps, lrs)
    return Pair{Int64, Vector{T}}[
        speciesmap(lrs.rs)[spat_rates[1]] => spat_rates[2] for spat_rates in transport_rates_speciesmap
    ]
end

### Accessing Unknown & Parameter Array Values ###

# Gets the index in the u array of species s in vertex vert (when their are num_species species).
get_index(vert::Int64, s::Int64, num_species::Int64) = (vert - 1) * num_species + s
# Gets the indexes in the u array of all species in vertex vert (when their are num_species species).
get_indexes(vert::Int64, num_species::Int64) = ((vert - 1) * num_species + 1):(vert * num_species)

# For vectors of length 1 or n, we want to get value idx (or the one value, if length is 1).
# This function gets that. Here:
# - values is the vector with the values of the component across all locations
#   (where the internal vectors may or may not be of size 1).
# - component_idx is the initial condition species/vertex parameter/edge parameters's index.
#   This is predominantly used for parameters, for initial conditions,
#   it is only used once (at initialisation) to re-process the input vector.
# - location_idx is the index of the vertex or edge for which we wish to access a initial condition or parameter values.
# The first two function takes the full value vector, and call the function of at the components specific index.
function get_component_value(values::Vector{<:Vector}, component_idx::Int64,
                             location_idx::Int64)
    get_component_value(values[component_idx], location_idx)
end
# Sometimes we have pre-computed, for each component, whether it's vector is length 1 or not.
# This is stored in location_types.
function get_component_value(values::Vector{<:Vector}, component_idx::Int64,
                             location_idx::Int64, location_types::Vector{Bool})     
    get_component_value(values[component_idx], location_idx, location_types[component_idx])
end
# For a components value (which is a vector of either length 1 or some other length), retrieves its value.
function get_component_value(values::Vector{<:Number}, location_idx::Int64)
    get_component_value(values, location_idx, length(values) == 1)
end
# Again, the location type (length of the value vector) may be pre-computed.
function get_component_value(values::Vector{<:Number}, location_idx::Int64,
                             location_type::Bool)
    location_type ? values[1] : values[location_idx]                                
end

# Converts a vector of vectors to a long vector.
# These are used when the initial condition is converted to a single vector (from vector of vector form).
function expand_component_values(values::Vector{<:Vector}, n)
    vcat([get_component_value.(values, comp) for comp in 1:n]...)
end
function expand_component_values(values::Vector{<:Vector}, n, location_types::Vector{Bool})
    vcat([get_component_value.(values, comp, location_types) for comp in 1:n]...)
end

# Creates a view of the vert_ps vector at a given location.
# Provides a work vector to which the converted vector is written. 
function view_vert_ps_vector!(work_vert_ps, vert_ps, comp, enumerated_vert_ps_idx_types)
    # Loops through all parameters.
    for (idx,loc_type) in enumerated_vert_ps_idx_types      
        # If the parameter is uniform across the spatial structure, it will have a length-1 value vector
        # (which value we write to the work vector).
        # Else, we extract it value at the specific location.
        work_vert_ps[idx] = (loc_type ? vert_ps[idx][1] : vert_ps[idx][comp])   
    end
    return work_vert_ps
end

# Expands a u0/p information stored in Vector{Vector{}} for to Matrix form
# (currently only used in Spatial Jump systems).
function matrix_expand_component_values(values::Vector{<:Vector}, n)
    reshape(expand_component_values(values, n), length(values), n)
end
