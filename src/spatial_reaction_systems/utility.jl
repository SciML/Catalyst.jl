### Processes Input u0 & p ###

# Required to make symmapt_to_varmap to work.
function _symbol_to_var(lrs::LatticeReactionSystem, sym)
    p_idx = findfirst(sym==p_sym for p_sym in ModelingToolkit.getname.(parameters(lrs)))
    isnothing(p_idx) || return parameters(lrs)[p_idx]
    s_idx = findfirst(sym==s_sym for s_sym in ModelingToolkit.getname.(species(lrs)))
    isnothing(s_idx) || return species(lrs)[s_idx]
    error("Could not find property parameter/species $sym in lattice reaction system.")
end

# From u0 input, extracts their values and store them in the internal format.
function lattice_process_u0(u0_in, u0_syms, nV)
    u0 = lattice_process_input(u0_in, u0_syms, nV)
    check_vector_lengths(u0, length(u0_syms), nV)
    expand_component_values(u0, nV)
end

# From p input, splits it into diffusion parameters and compartment parameters, and store these in the desired internal format.
function lattice_process_p(p_in, p_vertex_syms, p_edge_syms, lrs::LatticeReactionSystem)
    pV_in, pE_in = split_parameters(p_in, p_vertex_syms, p_edge_syms)
    pV = lattice_process_input(pV_in, p_vertex_syms, lrs.nV)
    pE = lattice_process_input(pE_in, p_edge_syms, lrs.nE)
    lrs.init_digraph || duplicate_spat_params!(pE, lrs)
    check_vector_lengths(pV, length(p_vertex_syms), lrs.nV)
    check_vector_lengths(pE, length(p_edge_syms), lrs.nE)
    return pV, pE
end

# Splits parameters into those for the compartments and those for the connections.
split_parameters(ps::Tuple{<:Any, <:Any}, args...) = ps
function split_parameters(ps::Vector{<:Number}, args...)
    error("When providing parameters for a spatial system as a single vector, the paired form (e.g :D =>1.0) must be used.")
end
function split_parameters(ps::Vector{<:Pair}, p_vertex_syms::Vector, p_edge_syms::Vector)
    pV_in = [p for p in ps if any(isequal(p[1]), p_vertex_syms)]
    pE_in = [p for p in ps if any(isequal(p[1]), p_edge_syms)]
    (sum(length.([pV_in, pE_in])) != length(ps)) && error("These input parameters are not recongised: $(setdiff(first.(ps), vcat(first.([pV_in, pE_in]))))")
    return pV_in, pE_in
end

# If the input is given in a map form, the vector needs sorting and the first value removed.
function lattice_process_input(input::Vector{<:Pair}, syms::Vector{BasicSymbolic{Real}}, args...)
    isempty(setdiff(first.(input), syms)) || error("Some input symbols are not recognised: $(setdiff(first.(input), syms)).")
    sorted_input = sort(input; by = p -> findfirst(isequal(p[1]), syms))
    return lattice_process_input(last.(sorted_input), syms, args...)
end
# Processes the input and gvies it in a form where it is a vector of vectors (some of which may have a single value).
function lattice_process_input(input::Matrix{<:Number}, args...)
    lattice_process_input([vec(input[i, :]) for i in 1:size(input, 1)], args...)
end
function lattice_process_input(input::Array{<:Number, 3}, args...)
    error("3 dimensional array parameter inpur currently not supported.")
end
function lattice_process_input(input::Vector{<:Any}, args...)
    isempty(input) ? Vector{Vector{Float64}}() :
    lattice_process_input([(val isa Vector{<:Number}) ? val : [val] for val in input],
                          args...)
end
lattice_process_input(input::Vector{<:Vector}, syms::Vector{BasicSymbolic{Real}}, n::Int64) = input

# Checks that a value vector have the right length, as well as that of all its sub vectors.
function check_vector_lengths(input::Vector{<:Vector}, n_syms, n_locations)
    (length(input)==n_syms) || error("Missing values for some initial conditions/parameters. Expected $n_syms values, got $(length(input)).")
    isempty(setdiff(unique(length.(input)), [1, n_locations])) || error("Some inputs where given values of inappropriate length.")
end

# For spatial parameters, if the lattice was given as an undirected graph of size n this is converted to a directed graph of size 2n.
# If spatial parameters are given with n values, we want to sue the same value for both directions.
# Since the order of edges in the new graph is non-trivial, this function distributes the n input values to a 2n length vector, putting teh correct value in each position.
function duplicate_spat_params!(pE::Vector{Vector{Float64}}, lrs::LatticeReactionSystem)
    cum_adjacency_counts = [0;cumsum(length.(lrs.lattice.fadjlist[1:end-1]))]
    for idx in 1:length(pE)
        (2length(pE[idx]) == lrs.nE) || continue # Only apply the function if the parameter have n values.
        
        new_vals = Vector{Float64}(undef,lrs.nE) # The new values.
        original_edge_count = 0                  # As we loop through the edges of the di-graph, this keeps track of each edge's index in the original graph. 
        for edge in edges(lrs.lattice)
            (edge.src < edge.dst) ? (original_edge_count += 1) : continue    # The digraph convertion only adds edges so that src > dst.
            idx_fwd = cum_adjacency_counts[edge.src] + findfirst(isequal(edge.dst),lrs.lattice.fadjlist[edge.src]) # For original edge i -> j, finds the index of i -> j in DiGraph.
            idx_bwd = cum_adjacency_counts[edge.dst] + findfirst(isequal(edge.src),lrs.lattice.fadjlist[edge.dst]) # For original edge i -> j, finds the index of j -> i in DiGraph.
            new_vals[idx_fwd] = pE[idx][original_edge_count]
            new_vals[idx_bwd] = pE[idx][original_edge_count]
        end
        pE[idx] = new_vals
    end
end

# For a set of input values on the given forms, and their symbolics, convert into a dictionary.
vals_to_dict(syms::Vector, vals::Vector{<:Vector}) = Dict(zip(syms, vals))
# Produces a dictionary with all parameter values.
function param_dict(pV, pE, lrs)
    merge(vals_to_dict(vertex_parameters(lrs), pV),
          vals_to_dict(edge_parameters(lrs), pE))
end

# Computes the spatial rates and stores them in a format (Dictionary of species index to rates across all edges).
function compute_all_spatial_rates(pV::Vector{Vector{Float64}}, pE::Vector{Vector{Float64}}, lrs::LatticeReactionSystem)
    param_value_dict = param_dict(pV, pE, lrs)
    unsorted_rates = [s => Symbolics.value.(compute_spatial_rates(get_spatial_rate_law(s, lrs), param_value_dict, lrs.nE)) for s in spatial_species(lrs)]
    return sort(unsorted_rates; by=rate -> findfirst(isequal(rate[1]), species(lrs)))
end
function get_spatial_rate_law(s::BasicSymbolic{Real}, lrs::LatticeReactionSystem)
    rates = filter(sr -> isequal(s, sr.species), lrs.spatial_reactions)
    (length(rates) > 1) && error("Species $s have more than one diffusion reaction.")    # We could allows several and simply sum them though, easy change.
    return rates[1].rate
end
function compute_spatial_rates(rate_law::Num,
                                 param_value_dict::Dict{SymbolicUtils.BasicSymbolic{Real}, Vector{Float64}}, nE::Int64)
    relevant_parameters = Symbolics.get_variables(rate_law)
    if all(length(param_value_dict[P]) == 1 for P in relevant_parameters)
        return [
            substitute(rate_law,
                       Dict(p => param_value_dict[p][1] for p in relevant_parameters)),
        ]
    end
    return [substitute(rate_law,
                       Dict(p => get_component_value(param_value_dict[p], idxE)
                            for p in relevant_parameters)) for idxE in 1:nE]
end

### Accessing State & Parameter Array Values ###

# Gets the index in the u array of species s in compartment comp (when their are nS species).
get_index(comp::Int64, s::Int64, nS::Int64) = (comp - 1) * nS + s
# Gets the indexes in the u array of all species in comaprtment comp (when their are nS species).
get_indexes(comp::Int64, nS::Int64) = ((comp - 1) * nS + 1):(comp * nS)

# We have many vectors of length 1 or n, for which we want to get value idx (or the one value, if length is 1), this function gets that.
function get_component_value(values::Vector{<:Vector}, component_idx::Int64,
                             location_idx::Int64)
    get_component_value(values[component_idx], location_idx)
end
function get_component_value(values::Vector{<:Vector}, component_idx::Int64,
                             location_idx::Int64, location_types::Vector{Bool})
    get_component_value(values[component_idx], location_idx, location_types[component_idx])
end
function get_component_value(values::Vector{<:Number}, location_idx::Int64)
    get_component_value(values, location_idx, length(values) == 1)
end
function get_component_value(values::Vector{<:Number}, location_idx::Int64,
                             location_type::Bool)
    location_type ? values[1] : values[location_idx]
end
# Converts a vector of vectors to a long vector.
function expand_component_values(values::Vector{<:Vector}, n)
    vcat([get_component_value.(values, comp) for comp in 1:n]...)
end
function expand_component_values(values::Vector{<:Vector}, n, location_types::Vector{Bool})
    vcat([get_component_value.(values, comp, location_types) for comp in 1:n]...)
end
# Creates a view of the pV vector at a given comaprtment.
function view_pV_vector!(work_pV, pV, comp, enumerated_pV_idx_types)
    for (idx,loc_type) in enumerated_pV_idx_types
        work_pV[idx] = (loc_type ? pV[idx][1] : pV[idx][comp])
    end
    return work_pV
end
# Expands a u0/p information stored in Vector{Vector{}} for to Matrix form (currently used in Spatial Jump systems).
function matrix_expand_component_values(values::Vector{<:Vector}, n)
    reshape(expand_component_values(values, n), length(values), n)
end
