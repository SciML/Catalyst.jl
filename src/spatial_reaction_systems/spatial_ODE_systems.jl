

### Spatial ODE Functor Structures ###

# Functor structure containg the information for the forcing function of a spatial ODE with diffusion on a lattice.
struct LatticeDiffusionODEf{R,S,T}
    ofunc::R
    nC::Int64
    nS::Int64
    pC::Vector{Vector{Float64}}
    work_pC::Vector{Float64}    
    enumerated_pC_idx_types::Base.Iterators.Enumerate{BitVector}
    diffusion_rates::Vector{S}
    leaving_rates::Matrix{Float64}
    enumerated_edges::T
    
    function LatticeDiffusionODEf(ofunc::R, pC, diffusion_rates::Vector{S}, lrs::LatticeReactionSystem) where {R, S, T}
        leaving_rates = zeros(length(diffusion_rates), lrs.nC)
        for (s_idx, rates) in enumerate(last.(diffusion_rates)),
            (e_idx, e) in enumerate(edges(lrs.lattice))
    
            leaving_rates[s_idx, e.src] += get_component_value(rates, e_idx)
        end
        work_pC = zeros(lrs.nC)
        enumerated_pC_idx_types = enumerate(length.(pC) .== 1)
        enumerated_edges = deepcopy(enumerate(edges(lrs.lattice)))
        new{R,S,typeof(enumerated_edges)}(ofunc, lrs.nC, lrs.nS, pC, work_pC, enumerated_pC_idx_types, diffusion_rates, leaving_rates, enumerated_edges)
    end
end

# Functor structure containg the information for the forcing function of a spatial ODE with diffusion on a lattice.
struct LatticeDiffusionODEjac{S,T}
    ofunc::S
    nC::Int64
    nS::Int64
    pC::Vector{Vector{Float64}}
    work_pC::Vector{Float64}    
    enumerated_pC_idx_types::Base.Iterators.Enumerate{BitVector}
    sparse::Bool
    jac_values::T

    function LatticeDiffusionODEjac(ofunc::S, pC, lrs::LatticeReactionSystem, jac_prototype::Union{Nothing, SparseMatrixCSC{Float64, Int64}}, sparse::Bool) where {S, T}
        work_pC = zeros(lrs.nC)
        enumerated_pC_idx_types = enumerate(length.(pC) .== 1)
        jac_values = sparse ? jac_prototype.nzval : Matrix(jac_prototype)
        new{S,typeof(jac_values)}(ofunc, lrs.nC, lrs.nS, pC, work_pC, enumerated_pC_idx_types, sparse, jac_values)
    end
end

### ODEProblem ###

# Creates an ODEProblem from a LatticeReactionSystem.
function DiffEqBase.ODEProblem(lrs::LatticeReactionSystem, u0_in, tspan,
                               p_in = DiffEqBase.NullParameters(), args...;
                               jac = true, sparse = jac, kwargs...)
    u0 = lattice_process_u0(u0_in, ModelingToolkit.getname.(species(lrs)), lrs.nC)
    pC, pD = lattice_process_p(p_in, Symbol.(compartment_parameters(lrs)),
                               Symbol.(diffusion_parameters(lrs)), lrs)
    ofun = build_odefunction(lrs, pC, pD, jac, sparse)
    return ODEProblem(ofun, u0, tspan, pC, args...; kwargs...)
end

# Builds an ODEFunction for a spatial ODEProblem.
function build_odefunction(lrs::LatticeReactionSystem, pC::Vector{Vector{Float64}},
                           pD::Vector{Vector{Float64}}, use_jac::Bool, sparse::Bool)
    # Prepeares (non-spatial) ODE functions and list of diffusing species and their rates.
    ofunc = ODEFunction(convert(ODESystem, lrs.rs); jac = use_jac, sparse = false)
    ofunc_sparse = ODEFunction(convert(ODESystem, lrs.rs); jac = use_jac, sparse = true)
    diffusion_rates_speciesmap = compute_all_diffusion_rates(pC, pD, lrs)
    diffusion_rates = [findfirst(isequal(diff_rates[1]), states(lrs.rs)) => diff_rates[2]
                       for diff_rates in diffusion_rates_speciesmap]

    f = LatticeDiffusionODEf(ofunc, pC, diffusion_rates, lrs)
    jac_prototype = (use_jac || sparse) ?
                    build_jac_prototype(ofunc_sparse.jac_prototype, diffusion_rates,
                                        lrs; set_nonzero = use_jac) : nothing
    jac = use_jac ? LatticeDiffusionODEjac(ofunc, pC, lrs, jac_prototype, sparse) : nothing
    return ODEFunction(f; jac = jac, jac_prototype = (sparse ? jac_prototype : nothing))
end

# Builds a jacobian prototype. If requested, populate it with the Jacobian's (constant) values as well.
function build_jac_prototype(ns_jac_prototype::SparseMatrixCSC{Float64, Int64},
                             diffusion_rates, lrs::LatticeReactionSystem;
                             set_nonzero = false)
    diff_species = first.(diffusion_rates)
    # Gets list of indexes for species that diffuse, but are invovled in no other reaction.
    only_diff = [(s in diff_species) && !Base.isstored(ns_jac_prototype, s, s)
                 for s in 1:(lrs.nS)]

    # Declares sparse array content.
    J_colptr = fill(1, lrs.nC * lrs.nS + 1)
    J_nzval = fill(0.0,
                   lrs.nC * (nnz(ns_jac_prototype) + count(only_diff)) +
                   length(edges(lrs.lattice)) * length(diffusion_rates))
    J_rowval = fill(0, length(J_nzval))

    # Finds filled elements.
    for comp in 1:(lrs.nC), s in 1:(lrs.nS)
        col_idx = get_index(comp, s, lrs.nS)

        # Column values.
        local_elements = in(s, diff_species) *
                         (length(lrs.lattice.fadjlist[comp]) + only_diff[s])
        diffusion_elements = -(ns_jac_prototype.colptr[(s + 1):-1:s]...)
        J_colptr[col_idx + 1] = J_colptr[col_idx] + local_elements + diffusion_elements

        # Row values.
        rows = ns_jac_prototype.rowval[ns_jac_prototype.colptr[s]:(ns_jac_prototype.colptr[s + 1] - 1)] .+
               (comp - 1) * lrs.nS
        if in(s, diff_species)
            # Finds the location of the diffusion elements, and inserts the elements from the non-spatial part into this.
            diffusion_rows = (lrs.lattice.fadjlist[comp] .- 1) .* lrs.nS .+ s
            split_idx = isempty(rows) ? 1 : findfirst(diffusion_rows .> rows[1])
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
    if !set_nonzero
        J_nzval .= 1.0
    else
        for (s_idx, (s, rates)) in enumerate(diffusion_rates),
            (e_idx, edge) in enumerate(edges(lrs.lattice))

            col_start = J_colptr[get_index(edge.src, s, lrs.nS)]
            col_end = J_colptr[get_index(edge.src, s, lrs.nS) + 1] - 1
            column_view = @view J_rowval[col_start:col_end]

            # Updates the source value.
            val_idx_src = col_start +
                          findfirst(column_view .== get_index(edge.src, s, lrs.nS)) - 1
            J_nzval[val_idx_src] -= get_component_value(rates, e_idx)

            # Updates the destination value.
            val_idx_dst = col_start +
                          findfirst(column_view .== get_index(edge.dst, s, lrs.nS)) - 1
            J_nzval[val_idx_dst] += get_component_value(rates, e_idx)
        end
    end

    return SparseMatrixCSC(lrs.nS * lrs.nC, lrs.nS * lrs.nC, J_colptr, J_rowval, J_nzval)
end

# Defines the forcing functors effect on the (spatial) ODE system.
function (f_func::LatticeDiffusionODEf)(du, u, p, t)
    # Updates for non-spatial reactions.
    for comp_i::Int64 in 1:(f_func.nC)
        f_func.ofunc((@view du[get_indexes(comp_i, f_func.nS)]),
              (@view u[get_indexes(comp_i, f_func.nS)]),
              view_pC_vector!(f_func.work_pC, p, comp_i, f_func.enumerated_pC_idx_types), t)
    end

    # Updates for spatial diffusion reactions.
    for (s_idx, (s, rates)) in enumerate(f_func.diffusion_rates)
        for comp_i::Int64 in 1:(f_func.nC)
            du[get_index(comp_i, s, f_func.nS)] -= f_func.leaving_rates[s_idx, comp_i] *
                                                u[get_index(comp_i, s,
                                                f_func.nS)]
        end
        for (e_idx::Int64, edge::Graphs.SimpleGraphs.SimpleEdge{Int64}) in f_func.enumerated_edges
            du[get_index(edge.dst, s, f_func.nS)] += get_component_value(rates, e_idx) *
                                                  u[get_index(edge.src, s,
                                                  f_func.nS)]
        end
    end
end

# Defines the jacobian functors effect on the (spatial) ODE system.
function (jac_func::LatticeDiffusionODEjac)(J, u, p, t)
    # Because of weird stuff where the Jacobian is not reset that I don't understand properly.
    reset_J_vals!(J)

    # Updates for non-spatial reactions.
    for comp_i::Int64 in 1:(jac_func.nC)
        jac_func.ofunc.jac((@view J[get_indexes(comp_i, jac_func.nS),
                           get_indexes(comp_i, jac_func.nS)]),
                  (@view u[get_indexes(comp_i, jac_func.nS)]),
                  view_pC_vector!(jac_func.work_pC, p, comp_i, jac_func.enumerated_pC_idx_types), t)
    end

    # Updates for the spatial reactions.
    add_diff_J_vals!(J, jac_func)
end
# Resets the jacobian matrix within a jac call.
reset_J_vals!(J::Matrix) = (J .= 0.0)
reset_J_vals!(J::SparseMatrixCSC) = (J.nzval .= 0.0)
# Updates the jacobian matrix with the difussion values.
add_diff_J_vals!(J::SparseMatrixCSC, jac_func::LatticeDiffusionODEjac) = (J.nzval .+= jac_func.jac_values)
add_diff_J_vals!(J::Matrix, jac_func::LatticeDiffusionODEjac) = (J .+= jac_func.jac_values)