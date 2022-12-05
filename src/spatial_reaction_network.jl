### Spatial Reaction Structure. ###
# Describing a spatial reaction that involves species from two neighbouring compartments.
# Currently only permit constant rate.
struct SpatialReaction
    rate::Symbol
    substrates::Tuple{Vector{Pair{Symbol,Int64}},Vector{Pair{Symbol,Int64}}}
    products::Tuple{Vector{Pair{Symbol,Int64}},Vector{Pair{Symbol,Int64}}}
    net_stoich_src::Vector{Pair{Symbol, Int64}}
    net_stoich_dst::Vector{Pair{Symbol, Int64}}

    function SpatialReaction(rate::Symbol, substrates_in::Tuple{Vector,Vector}, products_in::Tuple{Vector,Vector})
        substrates = process_sr_species.(substrates_in)
        products = process_sr_species.(products_in)
        new(rate, substrates, products, net_stoich(substrates,products,1), net_stoich(substrates,products,2))
    end
end
process_sr_species(species::Vector) = map(s -> (s isa Symbol) ? s => 1 : s, species);
function net_stoich(substrates, products, idx)
    sub_dict = Dict(Pair.(first.(substrates[idx]), -last.(substrates[idx])))
    prod_dict = Dict(products[idx])
    return collect(mergewith(-, prod_dict, sub_dict))
end;

function diffusion_reaction(rate, species)
    return SpatialReaction(rate, ([species], []), ([]), [species])
end;

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
    """The graph on which the lattice is defined."""
    lattice::DiGraph
    """Dependent (state) variables representing amount of each species. Must not contain the
    independent variable."""
    
    function LatticeReactionSystem(rs, spatial_reactions, lattice::DiGraph)
        return new(rs, spatial_reactions, lattice)
    end
    function LatticeReactionSystem(rs, spatial_reactions, lattice::SimpleGraph)
        return new(rs, spatial_reactions, DiGraph(lattice))
    end
end

### ODEProblem ###
# Creates an ODEProblem from a LatticeReactionSystem.
function DiffEqBase.ODEProblem(lrs::LatticeReactionSystem, u0, tspan,
                               p = DiffEqBase.NullParameters(), args...;
                               kwargs...)
    @unpack rs, spatial_reactions, lattice = lrs       
                        
    spatial_params = unique(getfield.(spatial_reactions, :rate))
    pV_in, pE_in = split_parameters(p, spatial_params)
    nV, nE = length.([vertices(lattice), edges(lattice)])
    u_idxs = Dict(reverse.(enumerate(Symbolics.getname.(states(rs)))))
    pV_idxes = Dict(reverse.(enumerate(Symbol.(parameters(rs)))))
    pE_idxes = Dict(reverse.(enumerate(unique(getfield.(spatial_reactions, :rate)))))

    u0 = matrix_form(u0, nV, u_idxs)
    pV = matrix_form(pV_in, nV, pV_idxes)
    pE = matrix_form(pE_in, nE, pE_idxes);

    return ODEProblem(build_f(lrs, u_idxs, pE_idxes), u0, tspan, (pV,pE), args...; kwargs...)
end

# Splits parameters into those for the compartments and those for the connections.
split_parameters(parameters::Tuple, spatial_params) = parameters
split_parameters(parameters::Vector, spatial_params) = filter(p -> !in(p[1], spatial_params), parameters), filter(p -> in(p[1], spatial_params), parameters);

# Converts species and parameters to matrices form.
matrix_form(input::Matrix, args...) = input
matrix_form(input::Vector{Pair{Symbol,Vector{Float64}}}, n, index_dict) = mapreduce(permutedims, vcat, last.(sort(input, by=i -> index_dict[i[1]])))
matrix_form(input::Vector, n, index_dict) = matrix_form(map(i -> (i[2] isa Vector) ? i[1] => i[2] : i[1] => fill(i[2], n), input), n, index_dict);

# Creates a function for simulating the spatial ODE with spatial reactions.
function build_f(lrs::LatticeReactionSystem, u_idxs::Dict{Symbol,Int64}, pE_idxes::Dict{Symbol,Int64})
    ofunc = ODEFunction(convert(ODESystem, lrs.rs))

    return function internal___spatial___f(du, u, p, t)
        # Updates for non-spatial reactions.
        for comp_i in 1:size(u,2)
            ofunc((@view du[:,comp_i]), (@view u[:,comp_i]), p[1], t)
        end
        
        # Updates for spatial reactions.
        for comp_i in 1:size(u,2)
            for comp_j::Int64 in (lrs.lattice.fadjlist::Vector{Vector{Int64}})[comp_i], sr::SpatialReaction in lrs.spatial_reactions::Vector{SpatialReaction}
                rate = get_rate(sr,p[2],(@view u[:,comp_i]),(@view u[:,comp_j]), u_idxs, pE_idxes)
                foreach(stoich -> du[u_idxs[stoich[1]],comp_i] += rate*stoich[2], sr.net_stoich_src)
                foreach(stoich -> du[u_idxs[stoich[1]],comp_j] += rate*stoich[2], sr.net_stoich_dst)
            end
        end
    end
end

# Get the rate of a specific reaction.
function get_rate(sr, pE, u_src, u_dst, u_idxs, pE_idxes)
    product = pE[pE_idxes[sr.rate]]
    !isempty(sr.substrates[1]) && (product *= prod(u_src[u_idxs[subs[1]]]^subs[2] / factorial(subs[2]) for subs in sr.substrates[1]))
    !isempty(sr.substrates[2]) && (product *= prod(u_dst[u_idxs[subs[1]]]^subs[2] / factorial(subs[2]) for subs in sr.substrates[2]))
    return product::Float64
end