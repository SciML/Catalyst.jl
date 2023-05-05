### Spatial Reaction Structure. ###
# Describing a spatial reaction that involves species from two neighbouring compartments.
# Currently only permit constant rate.
struct SpatialReaction
    """The rate function (excluding mass action terms). Currentl only cosntants supported"""
    rate::Symbol
    """Reaction substrates (source and destination)."""
    substrates::Tuple{Vector{Symbol}, Vector{Symbol}}
    """Reaction products (source and destination)."""
    products::Tuple{Vector{Symbol}, Vector{Symbol}}
    """The stoichiometric coefficients of the reactants (source and destination)."""
    substoich::Tuple{Vector{Int64}, Vector{Int64}}
    """The stoichiometric coefficients of the products (source and destination)."""
    prodstoich::Tuple{Vector{Int64}, Vector{Int64}}
    """The net stoichiometric coefficients of all species changed by the reaction (source and destination)."""
    netstoich::Tuple{Vector{Pair{Symbol,Int64}}, Vector{Pair{Symbol,Int64}}}
    """
    `false` (default) if `rate` should be multiplied by mass action terms to give the rate law.
    `true` if `rate` represents the full reaction rate law.
    Currently only `false`, is supported.
    """
    only_use_rate::Bool

    """These are similar to substrates, products, and netstoich, but ses species index (instead ) """
    function SpatialReaction(rate, substrates::Tuple{Vector, Vector}, products::Tuple{Vector, Vector}, substoich::Tuple{Vector{Int64}, Vector{Int64}}, prodstoich::Tuple{Vector{Int64}, Vector{Int64}}; only_use_rate = false)
        new(rate, substrates, products, substoich, prodstoich, get_netstoich.(substrates, products, substoich, prodstoich), only_use_rate)
    end
end

# As a spatial reaction, but replaces the species (and parameter) symbols with their index.
# For internal use only (to avoid having to constantly look up species indexes).
struct SpatialReactionIndexed
    rate::Int64
    substrates::Tuple{Vector{Int64}, Vector{Int64}}
    products::Tuple{Vector{Int64}, Vector{Int64}}
    substoich::Tuple{Vector{Int64}, Vector{Int64}}
    prodstoich::Tuple{Vector{Int64}, Vector{Int64}}
    netstoich::Tuple{Vector{Pair{Int64,Int64}}, Vector{Pair{Int64,Int64}}}
    only_use_rate::Bool

    function SpatialReactionIndexed(sr::SpatialReaction, species_list::Vector{Symbol}, param_list::Vector{Symbol})
        get_s_idx(species::Symbol) = findfirst(species .== (species_list))
        rate = findfirst(sr.rate .== (param_list))
        substrates = Tuple([get_s_idx.(sr.substrates[i]) for i in 1:2])
        products = Tuple([get_s_idx.(sr.products[i]) for i in 1:2])
        netstoich = Tuple([Pair.(get_s_idx.(first.(sr.netstoich[i])), last.(sr.netstoich[i])) for i in 1:2])
        new(rate, substrates, products, sr.substoich, sr.prodstoich, netstoich, sr.only_use_rate)
    end
end

"""
    DiffusionReaction(rate,species)

Simple function to create a diffusion spatial reaction. 
    Equivalent to SpatialReaction(rate,([species],[]),([],[species]),([1],[]),([],[1]))
"""
DiffusionReaction(rate,species) = SpatialReaction(rate,([species],Symbol[]),(Symbol[],[species]),([1],Int64[]),(Int64[],[1]))

"""
    OnewaySpatialReaction(rate, substrates, products, substoich, prodstoich)

Simple function to create a spatial reactions where all substrates are in teh soruce compartment, and all products in the destination.
Equivalent to SpatialReaction(rate,(substrates,[]),([],products),(substoich,[]),([],prodstoich))
"""
OnewaySpatialReaction(rate, substrates::Vector, products::Vector, substoich::Vector{Int64}, prodstoich::Vector{Int64}) = SpatialReaction(rate,(substrates,Symbol[]),(Symbol[],products),(substoich,Int64[]),(Int64[],prodstoich))



### Lattice Reaction Network Structure ###
# Couples:
# A reaction network (that is simulated within each compartment).
# A set of spatial reactions (denoting interaction between comaprtments).
# A network of compartments (a meta graph that can contain some additional infro for each compartment).
# The lattice is a DiGraph, normals graphs are converted to DiGraphs (with one edge in each direction).
struct LatticeReactionSystem # <: MT.AbstractTimeDependentSystem # Adding this part messes up show, disabling me from creating LRSs
    """The spatial reactions defined between individual nodes."""
    rs::ReactionSystem
    """The spatial reactions defined between individual nodes."""
    spatial_reactions::Vector{SpatialReaction}
    """A list of parameters that occur in the spatial reactions."""
    spatial_params::Vector{Symbol}
    """The graph on which the lattice is defined."""
    lattice::DiGraph
    """Dependent (state) variables representing amount of each species. Must not contain the
    independent variable."""

    function LatticeReactionSystem(rs, spatial_reactions, lattice::DiGraph)
        return new(rs, spatial_reactions, unique(getfield.(spatial_reactions, :rate)), lattice)
    end
    function LatticeReactionSystem(rs, spatial_reactions, lattice::SimpleGraph)
        return new(rs, spatial_reactions, unique(getfield.(spatial_reactions, :rate)), DiGraph(lattice))
    end
end

### ODEProblem ###
# Creates an ODEProblem from a LatticeReactionSystem.
function DiffEqBase.ODEProblem(lrs::LatticeReactionSystem, u0, tspan,
                               p = DiffEqBase.NullParameters(), args...;
                               jac=true, sparse=true, kwargs...)
    @unpack rs, spatial_reactions, lattice = lrs

    pV_in, pE_in = split_parameters(p, lrs.spatial_params)
    u_idxs = Dict(reverse.(enumerate(Symbolics.getname.(states(rs)))))
    pV_idxes = Dict(reverse.(enumerate(Symbol.(parameters(rs)))))
    pE_idxes = Dict(reverse.(enumerate(lrs.spatial_params)))

    nS,nV,nE = length.([states(lrs.rs), vertices(lrs.lattice), edges(lrs.lattice)])
    u0 = Vector(reshape(matrix_form(u0, nV, u_idxs),1:nS*nV))

    pV = matrix_form(pV_in, nV, pV_idxes)
    pE = matrix_form(pE_in, nE, pE_idxes)

    ofun = build_odefunction(lrs, jac, sparse)
    return ODEProblem(ofun, u0, tspan, (pV, pE), args...; kwargs...)
end

# Splits parameters into those for the compartments and those for the connections.
split_parameters(parameters::Tuple, spatial_params::Vector{Symbol}) = parameters
function split_parameters(parameters::Vector, spatial_params::Vector{Symbol})
    filter(p -> !in(p[1], spatial_params), parameters),
    filter(p -> in(p[1], spatial_params), parameters)
end

# Converts species and parameters to matrices form.
matrix_form(input::Matrix{Float64}, args...) = input
function matrix_form(input::Vector{Pair{Symbol, Vector{Float64}}}, n, index_dict)
    isempty(input) && return zeros(0,n)
    mapreduce(permutedims, vcat, last.(sort(input, by = i -> index_dict[i[1]])))
end
function matrix_form(input::Vector, n, index_dict)
    isempty(input) && return zeros(0,n)
    matrix_form(map(i -> (i[2] isa Vector) ? i[1] => i[2] : i[1] => fill(i[2], n), input),
                n, index_dict)
end

# Builds an ODEFunction.
function build_odefunction(lrs::LatticeReactionSystem, use_jac::Bool, sparse::Bool)
    ofunc = ODEFunction(convert(ODESystem, lrs.rs); jac=true)
    nS,nV = length.([states(lrs.rs), vertices(lrs.lattice)])
    spatial_reactions = [SpatialReactionIndexed(sr, map(s -> Symbol(s.f), species(lrs.rs)), lrs.spatial_params) for sr in lrs.spatial_reactions]

    f = build_f(ofunc, nS, nV, spatial_reactions, lrs.lattice.fadjlist)
    jac = (use_jac ? build_jac(ofunc,nS,nV,spatial_reactions, lrs.lattice.fadjlist) : nothing)
    jac_prototype = (use_jac ? build_jac_prototype(nS,nV,spatial_reactions, lrs, sparse) : nothing)

    return ODEFunction(f; jac=jac, jac_prototype=jac_prototype)
end

# Creates a function for simulating the spatial ODE with spatial reactions.
function build_f(ofunc::SciMLBase.AbstractODEFunction{true}, nS::Int64, nV::Int64, spatial_reactions::Vector{SpatialReactionIndexed}, adjlist::Vector{Vector{Int64}})
    return function(du, u, p, t)
        # Updates for non-spatial reactions.
        for comp_i::Int64 in 1:nV
            ofunc((@view du[get_indexes(comp_i,nS)]), (@view u[get_indexes(comp_i,nS)]), (@view p[1][:,comp_i]), t)
        end
    
        # Updates for spatial reactions.
        for comp_i::Int64 in 1:nV
            for comp_j::Int64 in adjlist[comp_i], sr::SpatialReactionIndexed in spatial_reactions
                rate::Float64 = get_rate(sr, p[2], (@view u[get_indexes(comp_i,nS)]), (@view u[get_indexes(comp_j,nS)]))
                for stoich::Pair{Int64,Int64} in sr.netstoich[1]
                    du[get_index(comp_i,stoich[1],nS)] += rate * stoich[2]
                end
                for stoich::Pair{Int64,Int64} in sr.netstoich[2]
                    du[get_index(comp_j,stoich[1],nS)] += rate * stoich[2]
                end
            end
        end
    end
end

# Get the rate of a specific reaction.
function get_rate(sr::SpatialReactionIndexed, pE::Matrix{Float64}, u_src, u_dst)
    product::Float64 = pE[sr.rate]
    !isempty(sr.substrates[1]) && for (sub::Int64,stoich::Int64) in zip(sr.substrates[1], sr.substoich[1])
        product *= u_src[sub]^stoich / factorial(stoich)
    end
    !isempty(sr.substrates[2]) && for (sub::Int64,stoich::Int64) in zip(sr.substrates[2], sr.substoich[2])
        product *= u_dst[sub]^stoich / factorial(stoich)
    end
    return product
end

function build_jac(ofunc::SciMLBase.AbstractODEFunction{true}, nS::Int64, nV::Int64, spatial_reactions::Vector{SpatialReactionIndexed}, adjlist::Vector{Vector{Int64}})
    return function(J, u, p, t)
        J .= 0

        # Updates for non-spatial reactions.
        for comp_i::Int64 in 1:nV
            ofunc.jac((@view J[get_indexes(comp_i,nS),get_indexes(comp_i,nS)]), (@view u[get_indexes(comp_i,nS)]), (@view p[1][:,comp_i]), t)
        end
    
        # Updates for spatial reactions.
        for comp_i::Int64 in 1:nV
            for comp_j::Int64 in adjlist[comp_i], sr::SpatialReactionIndexed in spatial_reactions
                for sub::Int64 in sr.substrates[1]
                    rate::Float64 = get_rate_differential(sr, p[2], sub, (@view u[get_indexes(comp_i,nS)]), (@view u[get_indexes(comp_j,nS)]))
                    for stoich::Pair{Int64,Int64} in sr.netstoich[1]
                        J[get_index(comp_i,stoich[1],nS),get_index(comp_i,sub,nS)] += rate * stoich[2]
                    end
                    for stoich::Pair{Int64,Int64} in sr.netstoich[2]
                        J[get_index(comp_j,stoich[1],nS),get_index(comp_i,sub,nS)] += rate * stoich[2]
                    end
                end
                for sub::Int64 in sr.substrates[2]
                    rate::Float64 = get_rate_differential(sr, p[2], sub, (@view u[get_indexes(comp_j,nS)]), (@view u[get_indexes(comp_i,nS)]))
                    for stoich::Pair{Int64,Int64} in sr.netstoich[1]
                        J[get_index(comp_i,stoich[1],nS),get_index(comp_j,sub,nS)] += rate * stoich[2]
                    end
                    for stoich::Pair{Int64,Int64} in sr.netstoich[2]
                        J[get_index(comp_j,stoich[1],nS),get_index(comp_j,sub,nS)] += rate * stoich[2]
                    end
                end
            end
        end
    end
end

function get_rate_differential(sr::SpatialReactionIndexed, pE::Matrix{Float64}, diff_species::Int64, u_src, u_dst)
    product::Float64 = pE[sr.rate]
    !isempty(sr.substrates[1]) && for (sub::Int64,stoich::Int64) in zip(sr.substrates[1], sr.substoich[1])
        if diff_species==sub
            product *= stoich*u_src[sub]^(stoich-1) / factorial(stoich)
        else
            product *= u_src[sub]^stoich / factorial(stoich)
        end
    end
    !isempty(sr.substrates[2]) && for (sub::Int64,stoich::Int64) in zip(sr.substrates[2], sr.substoich[2])
        product *= u_dst[sub]^stoich / factorial(stoich)
    end
    return product
end

function build_jac_prototype(nS::Int64, nV::Int64, spatial_reactions::Vector{SpatialReactionIndexed}, lrs::LatticeReactionSystem, sparse::Bool)
    jac_prototype = (sparse ? spzeros(nS*nV,nS*nV) : zeros(nS*nV,nS*nV)::Matrix{Float64})

    # Sets non-spatial reactions.
    # Tries to utilise sparsity within each comaprtment.
    for comp_i in 1:nV, reaction in reactions(lrs.rs)        
        for substrate in reaction.substrates, ns in reaction.netstoich
            sub_idx = findfirst(isequal(substrate), states(lrs.rs))
            spec_idx = findfirst(isequal(ns[1]), states(lrs.rs))
            jac_prototype[get_index(comp_i,spec_idx,nS), get_index(comp_i,sub_idx,nS)] = 1
        end
    end

    # Updates for spatial reactions.
    for comp_i::Int64 in 1:nV
        for comp_j::Int64 in lrs.lattice.fadjlist[comp_i]::Vector{Int64}, sr::SpatialReactionIndexed in spatial_reactions
            for sub::Int64 in sr.substrates[1]
                for stoich::Pair{Int64,Int64} in sr.netstoich[1]
                    jac_prototype[get_index(comp_i,stoich[1],nS),get_index(comp_i,sub,nS)] = 1.0
                end
                for stoich::Pair{Int64,Int64} in sr.netstoich[2]
                    jac_prototype[get_index(comp_j,stoich[1],nS),get_index(comp_i,sub,nS)] = 1.0
                end
            end
            for sub::Int64 in sr.substrates[2]
                for stoich::Pair{Int64,Int64} in sr.netstoich[1]
                    jac_prototype[get_index(comp_i,stoich[1],nS),get_index(comp_j,sub,nS)] = 1.0
                end
                for stoich::Pair{Int64,Int64} in sr.netstoich[2]
                    jac_prototype[get_index(comp_j,stoich[1],nS),get_index(comp_j,sub,nS)] = 1.0
                end
            end
        end
    end
    return jac_prototype
end

# Gets the index of a species (or a node's species) in the u array.
get_index(container::Int64,species::Int64,nS) = (container-1)*nS + species
get_indexes(container::Int64,nS) = (container-1)*nS+1:container*nS