### Diffusion Reaction Structure. ###

# Implements the diffusionparameter metadata field.
struct DiffusionParameter end
Symbolics.option_to_metadata_type(::Val{:diffusionparameter}) = DiffusionParameter

isdiffusionparameter(x::Num, args...) = isdiffusionparameter(Symbolics.unwrap(x), args...)
function isdiffusionparameter(x, default = false)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, DiffusionParameter, default)
end

# Abstract spatial reaction structures.
abstract type AbstractSpatialReaction end

# A diffusion reaction. These are simple to hanlde, and should cover most types of spatial reactions.
# Currently only permit constant rates.
struct DiffusionReaction <: AbstractSpatialReaction
    """The rate function (excluding mass action terms). Currently only constants supported"""
    rate::Num
    """The species that is subject to difusion."""
    species::Num
    """A symbol representation of the species that is subject to difusion."""
    species_sym::Symbol    # Required for identification in certain cases.

    # Creates a diffusion reaction.
    function DiffusionReaction(rate::Num, species::Num)
        new(rate, species, ModelingToolkit.getname(species))
    end
    function DiffusionReaction(rate::Number, species::Num)
        new(Num(rate), species, ModelingToolkit.getname(species))
    end
    function DiffusionReaction(rate::Symbol, species::Num)
        new(Symbolics.variable(rate), species, ModelingToolkit.getname(species))
    end
    function DiffusionReaction(rate::Num, species::Symbol)
        new(rate, Symbolics.variable(species), species)
    end
    function DiffusionReaction(rate::Number, species::Symbol)
        new(Num(rate), Symbolics.variable(species), species)
    end
    function DiffusionReaction(rate::Symbol, species::Symbol)
        new(Symbolics.variable(rate), Symbolics.variable(species), species)
    end
end
# Creates a vector of DiffusionReactions.
function diffusion_reactions(diffusion_reactions)
    [DiffusionReaction(dr[1], dr[2]) for dr in diffusion_reactions]
end
# Gets the parameters in a diffusion reaction.
ModelingToolkit.parameters(dr::DiffusionReaction) = Symbolics.get_variables(dr.rate)

### Lattice Reaction Network Structure ###
# Desribes a spatial reaction network over a graph.
struct LatticeReactionSystem # <: MT.AbstractTimeDependentSystem # Adding this part messes up show, disabling me from creating LRSs
    """The reaction system within each comaprtment."""
    rs::ReactionSystem
    """The spatial reactions defined between individual nodes."""
    spatial_reactions::Vector{<:AbstractSpatialReaction}
    """The graph on which the lattice is defined."""
    lattice::DiGraph

    # Derrived values.
    """The number of compartments."""
    nC::Int64
    """The number of edges."""
    nE::Int64
    """The number of species."""
    nS::Int64
    """Whenever the initial input was a di graph."""
    init_digraph::Bool

    function LatticeReactionSystem(rs::ReactionSystem,
                                   spatial_reactions::Vector{<:AbstractSpatialReaction},
                                   lattice::DiGraph; init_digraph = true)
        return new(rs, spatial_reactions, lattice, nv(lattice), ne(lattice),
                   length(species(rs)), init_digraph)
    end
    function LatticeReactionSystem(rs::ReactionSystem,
                                   spatial_reactions::Vector{<:AbstractSpatialReaction},
                                   lattice::SimpleGraph)
        return LatticeReactionSystem(rs, spatial_reactions, graph_to_digraph(lattice);
                                     init_digraph = false)
    end
    function LatticeReactionSystem(rs::ReactionSystem,
                                   spatial_reaction::AbstractSpatialReaction,
                                   lattice::Graphs.AbstractGraph)
        return LatticeReactionSystem(rs, [spatial_reaction], lattice)
    end
end
# Covnerts a graph to a digraph (in a way where we know where the new edges are in teh edge vector).
function graph_to_digraph(g1)
    g2 = Graphs.SimpleDiGraphFromIterator(reshape(permutedims(hcat(collect(edges(g1)),
                                                       reverse.(edges(g1)))), :, 1)[:])
    add_vertices!(g2, nv(g1) - nv(g2))
    return g2
end
# Gets the species of a lattice reaction system.
species(lrs::LatticeReactionSystem) = species(lrs.rs)
function diffusion_species(lrs::LatticeReactionSystem)
    filter(s -> ModelingToolkit.getname(s) in getfield.(lrs.spatial_reactions, :species_sym),
           species(lrs.rs))
end

# Gets the parameters in a lattice reaction system.
function ModelingToolkit.parameters(lrs::LatticeReactionSystem)
    unique(vcat(parameters(lrs.rs),
                Symbolics.get_variables.(getfield.(lrs.spatial_reactions, :rate))...))
end
function compartment_parameters(lrs::LatticeReactionSystem)
    filter(p -> !is_spatial_param(p, lrs), parameters(lrs))
end
function diffusion_parameters(lrs::LatticeReactionSystem)
    filter(p -> is_spatial_param(p, lrs), parameters(lrs))
end

# Checks whenever a parameter is a spatial parameter or not. 
function is_spatial_param(p, lrs)
    hasmetadata(p, DiffusionParameter) && getmetadata(p, DiffusionParameter) &&
        (return true)    # Wanted to just depend on metadata, but seems like we cannot implement that trivially.
    return (any(isequal(p), parameters(lrs.rs)) ? false : true)
end

### Processes Input u0 & p ###

# From u0 input, extracts their values and store them in the internal format.
function lattice_process_u0(u0_in, u0_symbols, nC)
    u0 = lattice_process_input(u0_in, u0_symbols, nC)
    check_vector_lengths(u0, nC)
    expand_component_values(u0, nC)
end

# From p input, splits it into diffusion parameters and compartment parameters, and store these in the desired internal format.
function lattice_process_p(p_in, p_comp_symbols, p_diff_symbols, lrs::LatticeReactionSystem)
    pC_in, pD_in = split_parameters(p_in, p_comp_symbols, p_diff_symbols)
    pC = lattice_process_input(pC_in, p_comp_symbols, lrs.nC)
    pD = lattice_process_input(pD_in, p_diff_symbols, lrs.nE)
    lrs.init_digraph || foreach(idx -> duplicate_diff_params!(pD, idx, lrs), 1:length(pD))
    check_vector_lengths(pC, lrs.nC)
    check_vector_lengths(pD, lrs.nE)
    return pC, pD
end

# Splits parameters into those for the compartments and those for the connections.
split_parameters(ps::Tuple{<:Any, <:Any}, args...) = ps
function split_parameters(ps::Vector{<:Number}, args...)
    error("When providing parameters for a spatial system as a single vector, the paired form (e.g :D =>1.0) must be used.")
end
function split_parameters(ps::Vector{<:Pair}, p_comp_symbols::Vector,
                          p_diff_symbols::Vector)
    pC_in = [p for p in ps if Symbol(p[1]) in p_comp_symbols]
    pD_in = [p for p in ps if Symbol(p[1]) in p_diff_symbols]
    (sum(length.([pC_in, pD_in])) != length(ps)) &&
        error("These input parameters are not recongised: $(setdiff(first.(ps), vcat(first.([pC_in, pE_in]))))")
    return pC_in, pD_in
end

# If the input is given in a map form, teh vector needs sorting and the first value removed.
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

# For diffusion parameters, if the graph was given as an undirected graph of length n, and the paraemter have n values, exapnd so that the same value are given for both values on the edge.
function duplicate_diff_params!(pD::Vector{Vector{Float64}}, idx::Int64,
                                lrs::LatticeReactionSystem)
    (2length(pD[idx]) == lrs.nE) && (pD[idx] = [p_val for p_val in pD[idx] for _ in 1:2])
end

# For a set of input values on the given forms, and their symbolics, convert into a dictionary.
vals_to_dict(syms::Vector, vals::Vector{<:Vector}) = Dict(zip(syms, vals))
# Produces a dictionary with all parameter values.
function param_dict(pC, pD, lrs)
    merge(vals_to_dict(compartment_parameters(lrs), pC),
          vals_to_dict(diffusion_parameters(lrs), pD))
end

# Computes the diffusion rates and stores them in a format (Dictionary of species index to rates across all edges).
function compute_all_diffusion_rates(pC::Vector{Vector{Float64}},
                                     pD::Vector{Vector{Float64}},
                                     lrs::LatticeReactionSystem)
    param_value_dict = param_dict(pC, pD, lrs)
    return [s => Symbolics.value.(compute_diffusion_rates(get_diffusion_rate_law(s, lrs),
                                                          param_value_dict, lrs.nE))
            for s in diffusion_species(lrs)]
end
function get_diffusion_rate_law(s::Symbolics.BasicSymbolic, lrs::LatticeReactionSystem)
    rates = filter(sr -> isequal(ModelingToolkit.getname(s), sr.species_sym),
                   lrs.spatial_reactions)
    (length(rates) > 1) && error("Species $s have more than one diffusion reaction.")    # We could allows several and simply sum them though, easy change.
    return rates[1].rate
end
function compute_diffusion_rates(rate_law::Num,
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
# Creates a view of the pC vector at a given comaprtment.
function view_pC_vector!(work_pC, pC, comp, enumerated_pC_idx_types)
    for (idx,loc_type) in enumerated_pC_idx_types
        work_pC[idx] = (loc_type ? pC[idx][1] : pC[idx][comp])
    end
    return work_pC
end
# Expands a u0/p information stored in Vector{Vector{}} for to Matrix form (currently used in Spatial Jump systems).
function matrix_expand_component_values(values::Vector{<:Vector}, n)
    reshape(expand_component_values(values, n), length(values), n)
end
