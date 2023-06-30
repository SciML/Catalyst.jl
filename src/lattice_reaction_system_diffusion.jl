### Spatial Reaction Structure. ###

# Abstract spatial reaction structures.
abstract type AbstractSpatialReaction end

# A diffusion reaction. These are simple to hanlde, and should cover most types of spatial reactions.
# Currently only permit constant rates.
struct DiffusionReaction <: AbstractSpatialReaction
    """The rate function (excluding mass action terms). Currentl only constants supported"""
    rate::Symbol
    """The species that is subject to difusion."""
    species::Symbol
end

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
    """A list of parameters that occur in the spatial reactions."""
    spatial_params::Vector{Symbol}
    """The number of compartments."""
    nC::Int64
    """The number of species."""
    nS::Int64

    function LatticeReactionSystem(rs, spatial_reactions, lattice::DiGraph)
        return new(rs, spatial_reactions, lattice,
                   unique(getfield.(spatial_reactions, :rate)), length(vertices(lattice)),
                   length(species(rs)))
    end
    function LatticeReactionSystem(rs, spatial_reactions, lattice::SimpleGraph)
        return LatticeReactionSystem(rs, spatial_reactions, DiGraph(lattice))
    end
end

### ODEProblem ###
# Creates an ODEProblem from a LatticeReactionSystem.
function DiffEqBase.ODEProblem(lrs::LatticeReactionSystem, u0_in, tspan,
                               p_in = DiffEqBase.NullParameters(), args...;
                               jac = true, sparse = true, kwargs...)
    u0 = resort_values(u0_in, [Symbol(s.f) for s in species(lrs.rs)])
    u0 = [get_component_value(u0, species, comp) for comp in 1:(lrs.nC)
          for species in 1:(lrs.nS)]
    pV, pE = split_parameters(p_in, lrs.spatial_params)
    pV = resort_values(pV, Symbol.(parameters(lrs.rs)))
    pE = resort_values(pE, lrs.spatial_params)

    ofun = build_odefunction(lrs, pV, pE, jac, sparse)
    return ODEProblem(ofun, u0, tspan, pV, args...; kwargs...)
end

# Splits parameters into those for the compartments and those for the connections.
split_parameters(ps::Tuple{<:Any, <:Any}, spatial_params::Vector{Symbol}) = ps
function split_parameters(ps::Vector{<:Pair}, spatial_params::Vector{Symbol})
    (filter(p -> !(p[1] in spatial_params), ps), filter(p -> p[1] in spatial_params, ps))
end
function split_parameters(ps::Vector{<:Number}, spatial_params::Vector{Symbol})
    (ps[1:(length(ps) - length(spatial_params))],
     ps[(length(ps) - length(spatial_params) + 1):end])
end
# Sorts a parameter (or species) vector along parameter (or species) index, and remove the Symbol in the pair.
function resort_values(values::Vector{<:Pair}, symbols::Vector{Symbol})
    last.(sort(values; by = val -> findfirst(val[1] .== symbols)))
end
resort_values(values::Any, symbols::Vector{Symbol}) = values

# Builds an ODEFunction.
function build_odefunction(lrs::LatticeReactionSystem, pV, pE, use_jac::Bool, sparse::Bool)
    ofunc = ODEFunction(convert(ODESystem, lrs.rs); jac = use_jac, sparse = false)
    ofunc_sparse = ODEFunction(convert(ODESystem, lrs.rs); jac = use_jac, sparse = true)
    diffusion_species = Int64[findfirst(s .== [Symbol(s.f) for s in species(lrs.rs)])
                              for s in getfield.(lrs.spatial_reactions, :species)]

    f = build_f(ofunc, pV, pE, diffusion_species, lrs)
    jac_prototype = (sparse ?
                     build_jac_prototype(ofunc_sparse.jac_prototype, pE, diffusion_species,
                                         lrs; set_nonzero = true) : nothing)
    jac = (use_jac ? build_jac(ofunc, pV, pE, diffusion_species, lrs, (isnothing(jac_prototype) ? build_jac_prototype(ofunc_sparse.jac_prototype, pE, diffusion_species, lrs; set_nonzero = true) : jac_prototype); sparse = sparse) : nothing)
    return ODEFunction(f; jac = jac, jac_prototype = jac_prototype)
end

# Creates a function for simulating the spatial ODE with spatial reactions.
function build_f(ofunc::SciMLBase.AbstractODEFunction{true}, pV, pE,
                 diffusion_species::Vector{Int64}, lrs::LatticeReactionSystem)
    leaving_rates = zeros(length(diffusion_species), lrs.nC)
    for (s_idx, species) in enumerate(diffusion_species),
        (e_idx, e) in enumerate(edges(lrs.lattice))

        leaving_rates[s_idx, e.src] += get_component_value(pE, s_idx, e_idx)
    end
    p_base = deepcopy(first.(pV))
    p_update_idx = (p_base isa Vector) ? findall(typeof.(p_base) .== Vector{Float64}) : []
    enumerated_edges = deepcopy(enumerate(edges(lrs.lattice)))

    return function (du, u, p, t)
        # Updates for non-spatial reactions.
        for comp_i::Int64 in 1:(lrs.nC)
            ofunc((@view du[get_indexes(comp_i, lrs.nS)]),
                  (@view u[get_indexes(comp_i, lrs.nS)]),
                  make_p_vector!(p_base, p, p_update_idx, comp_i), t)
        end

        for (s_idx::Int64, species::Int64) in enumerate(diffusion_species)
            for comp_i::Int64 in 1:(lrs.nC)
                du[get_index(comp_i, species, lrs.nS)] -= leaving_rates[s_idx, comp_i] *
                                                          u[get_index(comp_i, species,
                                                                      lrs.nS)]
            end
            for (e_idx::Int64, edge::Graphs.SimpleGraphs.SimpleEdge{Int64}) in enumerated_edges
                du[get_index(edge.dst, species, lrs.nS)] += get_component_value(pE, s_idx,
                                                                                e_idx) *
                                                            u[get_index(edge.src, species,
                                                                        lrs.nS)]
            end
        end
    end
end

function build_jac(ofunc::SciMLBase.AbstractODEFunction{true}, pV, pE,
                   diffusion_species::Vector{Int64}, lrs::LatticeReactionSystem,
                   jac_prototype::SparseMatrixCSC{Float64, Int64}; sparse = true)
    p_base = deepcopy(first.(pV))
    p_update_idx = (p_base isa Vector) ? findall(typeof.(p_base) .== Vector{Float64}) : []

    if sparse
        return function (J, u, p, t)
            # Sets the base according to the spatial reactions.
            J.nzval .= jac_prototype.nzval

            # Updates for non-spatial reactions.
            for comp_i::Int64 in 1:(lrs.nC)
                ofunc.jac((@view J[get_indexes(comp_i, lrs.nS),
                                   get_indexes(comp_i, lrs.nS)]),
                          (@view u[get_indexes(comp_i, lrs.nS)]),
                          make_p_vector!(p_base, p, p_update_idx, comp_i), t)
            end
        end
    else
        jac_diffusion = Matrix(jac_prototype)
        return function (J, u, p, t)
            # Sets the base according to the spatial reactions.
            J .= jac_diffusion

            # Updates for non-spatial reactions.
            for comp_i::Int64 in 1:(lrs.nC)
                ofunc.jac((@view J[get_indexes(comp_i, lrs.nS),
                                   get_indexes(comp_i, lrs.nS)]),
                          (@view u[get_indexes(comp_i, lrs.nS)]),
                          make_p_vector!(p_base, p, p_update_idx, comp_i), t)
            end
        end
    end
end

function build_jac_prototype(ns_jac_prototype::SparseMatrixCSC{Float64, Int64}, pE,
                             diffusion_species::Vector{Int64}, lrs::LatticeReactionSystem;
                             set_nonzero = false)
    only_diff = [(s in diffusion_species) && !Base.isstored(ns_jac_prototype, s, s)
                 for s in 1:(lrs.nS)]

    # Declares sparse array content.
    J_colptr = fill(1, lrs.nC * lrs.nS + 1)
    J_nzval = fill(0.0,
                   lrs.nC * (nnz(ns_jac_prototype) + count(only_diff)) +
                   length(edges(lrs.lattice)) * length(diffusion_species))
    J_rowval = fill(0, length(J_nzval))

    # Finds filled elements.
    for comp in 1:(lrs.nC), s in 1:(lrs.nS)
        col_idx = get_index(comp, s, lrs.nS)

        # Column values.
        local_elements = in(s, diffusion_species) *
                         (length(lrs.lattice.fadjlist[comp]) + only_diff[s])
        diffusion_elements = -(ns_jac_prototype.colptr[(s + 1):-1:s]...)
        J_colptr[col_idx + 1] = J_colptr[col_idx] + local_elements + diffusion_elements

        # Row values.
        rows = ns_jac_prototype.rowval[ns_jac_prototype.colptr[s]:(ns_jac_prototype.colptr[s + 1] - 1)] .+
               (comp - 1) * lrs.nS
        if in(s, diffusion_species)
            diffusion_rows = (lrs.lattice.fadjlist[comp] .- 1) .* lrs.nS .+ s
            split_idx = findfirst(diffusion_rows .> rows[1])
            isnothing(split_idx) && (split_idx = length(diffusion_rows) + 1)
            rows = vcat(diffusion_rows[1:(split_idx - 1)], rows,
                        diffusion_rows[split_idx:end])
            if only_diff[s]
                split_idx = findfirst(rows .> get_index(comp, s, lrs.nS))
                isnothing(split_idx) && (split_idx = length(rows) + 1)
                insert!(rows, split_idx, get_index(comp, s, lrs.nS))
            end
        end
        J_rowval[J_colptr[col_idx]:(J_colptr[col_idx + 1] - 1)] = rows
    end

    # Set element values.
    if set_nonzero
        J_nzval .= 1.0
    else
        for (s_idx, s) in enumerate(diffusion_species),
            (e_idx, edge) in enumerate(edges(lrs.lattice))

            col_start = J_colptr[get_index(edge.src, s, lrs.nS)]
            col_end = J_colptr[get_index(edge.src, s, lrs.nS) + 1] - 1
            column_view = @view J_rowval[col_start:col_end]

            # Updates the source value.
            val_idx_src = col_start +
                          findfirst(column_view .== get_index(edge.src, s_idx, lrs.nS)) - 1
            J_nzval[val_idx_src] -= get_component_value(pE, s_idx, e_idx)

            # Updates the destination value.
            val_idx_dst = col_start +
                          findfirst(column_view .== get_index(edge.dst, s_idx, lrs.nS)) - 1
            J_nzval[val_idx_dst] += get_component_value(pE, s_idx, e_idx)
        end
    end

    return SparseMatrixCSC(lrs.nS * lrs.nC, lrs.nS * lrs.nC, J_colptr, J_rowval, J_nzval)
end

# Gets the index of a species in the u array.
get_index(comp::Int64, species::Int64, nS::Int64) = (comp - 1) * nS + species
# Gets the indexes of a compartment's species in the u array.
get_indexes(comp::Int64, nS::Int64) = ((comp - 1) * nS + 1):(comp * nS)

# For set of values (stoed in a variety of possible forms), a given component (species or parameter), and a place (eitehr a compartment or edge), find that components value at that place.
function get_component_value(vals::Matrix{Float64}, component::Int64, place::Int64)
    vals[component, place]
end
function get_component_value(vals::Vector, component::Int64, place::Int64)
    get_component_value(vals[component], place)
end
get_component_value(vals::Vector{Float64}, place::Int64) = vals[place]
get_component_value(vals::Float64, place::Int64) = vals

# Updated the base parameter vector with the values for a specific compartment.
make_p_vector!(p_base, p::Matrix{Float64}, p_update_idx, comp) = p[:, comp]
function make_p_vector!(p_base, p, p_update_idx, comp_i)
    for idx in p_update_idx
        p_base[idx] = p[idx][comp_i]
    end
    return p_base
end

# Gets all the values in a specific palce.
function get_component_values(vals::Matrix{Float64}, place::Int64, nComponents::Int64)
    (@view vals[:, place])
end
function get_component_values(vals, place::Int64, nComponents::Int64)
    [get_component_value(vals, component, place) for component in 1:nComponents]
end
