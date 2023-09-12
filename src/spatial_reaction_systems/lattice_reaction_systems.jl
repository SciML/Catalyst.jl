### Lattice Reaction Network Structure ###
# Desribes a spatial reaction network over a graph.
struct LatticeReactionSystem{S,T} # <: MT.AbstractTimeDependentSystem # Adding this part messes up show, disabling me from creating LRSs
    """The reaction system within each comaprtment."""
    rs::ReactionSystem{S}
    """The spatial reactions defined between individual nodes."""
    spatial_reactions::Vector{T}
    """The graph on which the lattice is defined."""
    lattice::SimpleDiGraph{Int64}

    # Derrived values.
    """The number of compartments."""
    nV::Int64
    """The number of edges."""
    nE::Int64
    """The number of species."""
    nS::Int64
    """Whenever the initial input was a digraph."""
    init_digraph::Bool
    """Species that may move spatially."""
    spat_species::Vector{BasicSymbolic{Real}}
    """Parameters which values are tied to edges (adjacencies).."""
    edge_parameters::Vector{BasicSymbolic{Real}}

    function LatticeReactionSystem(rs::ReactionSystem{S},
                                   spatial_reactions::Vector{T},
                                   lattice::DiGraph; init_digraph = true) where {S, T}
        (T <: AbstractSpatialReaction) || error("The secodn argument must be a vector of AbstractSpatialReaction subtypes.") # There probably some better way to acertain that T has that type. Not sure how.
        spat_species = unique(vcat(spatial_species.(spatial_reactions)...))
        rs_edge_parameters = filter(isedgeparameter, parameters(rs))
        srs_edge_parameters = setdiff(vcat(parameters.(spatial_reactions)...), parameters(rs))
        edge_parameters = unique([rs_edge_parameters; srs_edge_parameters])


        foreach(sr -> check_spatial_reaction_validity(rs, sr), spatial_reactions)   
        return new{S,T}(rs, spatial_reactions, lattice, nv(lattice), ne(lattice), length(species(rs)), init_digraph, spat_species, edge_parameters)
    end
end
function LatticeReactionSystem(rs, srs, lat::SimpleGraph)
    return LatticeReactionSystem(rs, srs, graph_to_digraph(lat); init_digraph = false)
end
function LatticeReactionSystem(rs, sr::AbstractSpatialReaction, lat)
    return LatticeReactionSystem(rs, [sr], lat)
end
function LatticeReactionSystem(rs, sr::AbstractSpatialReaction, lat::SimpleGraph)
    return LatticeReactionSystem(rs, [sr], graph_to_digraph(lat); init_digraph = false)
end

# Converts a graph to a digraph (in a way where we know where the new edges are in teh edge vector).
function graph_to_digraph(g1)
    g2 = Graphs.SimpleDiGraphFromIterator(reshape(permutedims(hcat(collect(edges(g1)),
                                                       reverse.(edges(g1)))), :, 1)[:])
    add_vertices!(g2, nv(g1) - nv(g2))
    return g2
end

### Lattice ReactionSystem Getters ###

# Get all species.
species(lrs::LatticeReactionSystem) = unique([species(lrs.rs); lrs.spat_species])
# Get all species that may be transported.
spatial_species(lrs::LatticeReactionSystem) = lrs.spat_species

# Get all parameters.
ModelingToolkit.parameters(lrs::LatticeReactionSystem) = unique([species(lrs.rs); lrs.edge_parameters])
# Get all parameters which values are tied to vertexes (compartments).
vertex_parameters(lrs::LatticeReactionSystem) = setdiff(parameters(lrs), edge_parameters(lrs))
# Get all parameters which values are tied to edges (adjacencies).
edge_parameters(lrs::LatticeReactionSystem) = lrs.edge_parameters

# Checks if a lattice reaction system is a pure (linear) transport reaction system.
is_transport_system(lrs::LatticeReactionSystem) = all(typeof.(lrs.spatial_reactions) .== TransportReaction)

### Processes Input u0 & p ###

# From u0 input, extracts their values and store them in the internal format.
function lattice_process_u0(u0_in, u0_symbols, nV)
    u0 = lattice_process_input(u0_in, u0_symbols, nV)
    check_vector_lengths(u0, nV)
    expand_component_values(u0, nV)
end

# From p input, splits it into diffusion parameters and compartment parameters, and store these in the desired internal format.
function lattice_process_p(p_in, p_comp_symbols, p_diff_symbols, lrs::LatticeReactionSystem)
    pV_in, pE_in = split_parameters(p_in, p_comp_symbols, p_diff_symbols)
    pV = lattice_process_input(pV_in, p_comp_symbols, lrs.nV)
    pE = lattice_process_input(pE_in, p_diff_symbols, lrs.nE)
    lrs.init_digraph || foreach(idx -> duplicate_spat_params!(pE, idx, lrs), 1:length(pE))
    check_vector_lengths(pV, lrs.nV)
    check_vector_lengths(pE, lrs.nE)
    return pV, pE
end

# Splits parameters into those for the compartments and those for the connections.
split_parameters(ps::Tuple{<:Any, <:Any}, args...) = ps
function split_parameters(ps::Vector{<:Number}, args...)
    error("When providing parameters for a spatial system as a single vector, the paired form (e.g :D =>1.0) must be used.")
end
function split_parameters(ps::Vector{<:Pair}, p_comp_symbols::Vector,
                          p_diff_symbols::Vector)
    pV_in = [p for p in ps if Symbol(p[1]) in p_comp_symbols]
    pE_in = [p for p in ps if Symbol(p[1]) in p_diff_symbols]
    (sum(length.([pV_in, pE_in])) != length(ps)) &&
        error("These input parameters are not recongised: $(setdiff(first.(ps), vcat(first.([pV_in, pE_in]))))")
    return pV_in, pE_in
end

# If the input is given in a map form, the vector needs sorting and the first value removed.
function lattice_process_input(input::Vector{<:Pair}, symbols::Vector{Symbol}, args...)
    (length(setdiff(Symbol.(first.(input)), symbols)) != 0) &&
        error("Some input symbols are not recognised: $(setdiff(Symbol.(first.(input)), symbols)).")
    sorted_input = sort(input;
                        by = p -> findfirst(ModelingToolkit.getname(p[1]) .== symbols))
    return lattice_process_input(last.(sorted_input), symbols, args...)
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
lattice_process_input(input::Vector{<:Vector}, symbols::Vector{Symbol}, n::Int64) = input
function check_vector_lengths(input::Vector{<:Vector}, n)
    isempty(setdiff(unique(length.(input)), [1, n])) ||
        error("Some inputs where given values of inappropriate length.")
end

# For spatial parameters, if the graph was given as an undirected graph of length n, and the paraemter have n values, expand so that the same value are given for both values on the edge.
function duplicate_spat_params!(pE::Vector{Vector{Float64}}, idx::Int64,
                                lrs::LatticeReactionSystem)
    (2length(pE[idx]) == lrs.nE) && (pE[idx] = [p_val for p_val in pE[idx] for _ in 1:2])
end

# For a set of input values on the given forms, and their symbolics, convert into a dictionary.
vals_to_dict(syms::Vector, vals::Vector{<:Vector}) = Dict(zip(syms, vals))
# Produces a dictionary with all parameter values.
function param_dict(pV, pE, lrs)
    merge(vals_to_dict(compartment_parameters(lrs), pV),
          vals_to_dict(diffusion_parameters(lrs), pE))
end

# Computes the spatial rates and stores them in a format (Dictionary of species index to rates across all edges).
function compute_all_spatial_rates(pV::Vector{Vector{Float64}},
                                     pE::Vector{Vector{Float64}},
                                     lrs::LatticeReactionSystem)
    param_value_dict = param_dict(pV, pE, lrs)
    return [s => Symbolics.value.(compute_spatial_rates(get_spatial_rate_law(s, lrs),
                                                          param_value_dict, lrs.nE))
            for s in spatial_species(lrs)]
end
function get_spatial_rate_law(s::Symbolics.BasicSymbolic, lrs::LatticeReactionSystem)
    rates = filter(sr -> isequal(ModelingToolkit.getname(s), sr.species_sym),
                   lrs.spatial_reactions)
    (length(rates) > 1) && error("Species $s have more than one diffusion reaction.")    # We could allows several and simply sum them though, easy change.
    return rates[1].rate
end
function compute_spatial_rates(rate_law::Num,
                                 param_value_dict::Dict{Any, Vector{Float64}}, nE::Int64)
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
