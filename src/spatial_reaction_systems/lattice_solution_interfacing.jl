### Rudimentary Interfacing Function ###
# A single function, `get_lrs_vals`, which contain all interfacing functionality. However, 
# long-term it should be replaced with a sleeker interface. Ideally as MTK-wider support for 
# lattice problems and solutions are introduced. 

"""
    get_lrs_vals(sol, sp, lrs::LatticeReactionSystem; t = nothing)

A function for retrieving the solution of a `LatticeReactionSystem`-based simulation on various
desired forms. Generally, for `LatticeReactionSystem`s, the values in `sol` is ordered in a 
way which is not directly interpretable by the user. Furthermore, the normal Catalyst interface
for solutions (e.g. `sol[:X]`) does not work for these solutions. Hence this function is used instead.

The output is a vector, which in each position contain sp's value (either at a time step of time,
depending on the input `t`). Its shape depends on the lattice (using a similar form as heterogeneous
initial conditions). I.e. for a NxM cartesian grid, the values are NxM matrices. For a masked grid,
the values are sparse matrices. For a graph lattice, the values are vectors (where the value in
the n'th position corresponds to sp's value in the n'th vertex).

Arguments:
- `sol`: The solution from which we wish to retrieve some values.
- `sp`: The species which values we wish to retrieve. Can be either a symbol (e.g. `:X`) or a symbolic
variable (e.g. `X`).
- `lrs`: The `LatticeReactionSystem` which was simulated to generate the solution.
- `t = nothing`: If `nothing`, we simply returns the solution across all saved timesteps. If `t`
instead is a vector (or range of values), returns the solutions interpolated at these timepoints.

Notes:
- The `get_lrs_vals` is not optimised for performance. However, it should still be quite performant,
but there might be some limitations if called a very large number of times.
- Long-term it is likely that this function gets replaced with a sleeker interface.

Example:
```julia
using Catalyst, OrdinaryDiffEq

# Prepare `LatticeReactionSystem`s.
rs = @reaction_network begin
    (k1,k2), X1 <--> X2
end
tr = @transport_reaction D X1 
lrs = LatticeReactionSystem(rs, [tr], CartesianGrid((2,2)))

# Create problems.
u0 = [:X1 => 1, :X2 => 2]
tspan = (0.0, 10.0)
ps = [:k1 => 1, :k2 => 2.0, :D => 0.1]

oprob = ODEProblem(lrs1, u0, tspan, ps)
osol = solve(oprob1, Tsit5())
get_lrs_vals(osol, :X1, lrs) # Returns the value of X1 at each timestep.
get_lrs_vals(osol, :X1, lrs; t = 0.0:10.0) # Returns the value of X1 at times 0.0, 1.0, ..., 10.0
```
"""
function get_lrs_vals(sol, sp, lrs::LatticeReactionSystem; t = nothing)
    # Figures out which species we wish to fetch information about.
    (sp isa Symbol) && (sp = Catalyst._symbol_to_var(lrs, sp))
    sp_idx = findfirst(isequal(sp), species(lrs))
    sp_tot = length(species(lrs))

    # Extracts the lattice and calls the next function. Masked grids (Array of Bools) are converted
    # to sparse array using the same template size as we wish to shape the data to.
    lattice = Catalyst.lattice(lrs)
    if has_masked_lattice(lrs) 
        if grid_dims(lrs) == 3
            error("The `get_lrs_vals` function is not defined for systems based on 3d sparse arrays. Please raise an issue at the Catalyst GitHub site if this is something which would be useful to you.")
        end
        lattice = sparse(lattice)
    end
    get_lrs_vals(sol, lattice, t, sp_idx, sp_tot)
end

# Function which handles the input in the case where `t` is `nothing` (i.e. return `sp`s value
# across all sample points).
function get_lrs_vals(sol, lattice, t::Nothing, sp_idx, sp_tot)
    # ODE simulations contain, in each data point, all values in a single vector. Jump simulations
    # instead in a matrix (NxM, where N is the number of species and M the number of vertices). We
    # must consider each case separately.
    if sol.prob isa ODEProblem
        return [reshape_vals(vals[sp_idx:sp_tot:end], lattice) for vals in sol.u]
    elseif sol.prob isa DiscreteProblem
        return [reshape_vals(vals[sp_idx,:], lattice) for vals in sol.u]
    else
        error("Unknown type of solution provided to `get_lrs_vals`. Only ODE or Jump solutions are supported.")
    end
end

# Function which handles the input in the case where `t` is a range of values (i.e. return `sp`s 
# value at all designated time points.
function get_lrs_vals(sol, lattice, t::AbstractVector{T}, sp_idx, sp_tot) where {T <: Number}
    if (minimum(t) < sol.t[1]) || (maximum(t) > sol.t[end])
        error("The range of the t values provided for sampling, ($(minimum(t)),$(maximum(t))) is not fully within the range of the simulation time span ($(sol.t[1]),$(sol.t[end])).")
    end

    # ODE simulations contain, in each data point, all values in a single vector. Jump simulations
    # instead in a matrix (NxM, where N is the number of species and M the number of vertices). We
    # must consider each case separately.
    if sol.prob isa ODEProblem
        return [reshape_vals(sol(ti)[sp_idx:sp_tot:end], lattice) for ti in t]
    elseif sol.prob isa DiscreteProblem
        return [reshape_vals(sol(ti)[sp_idx,:], lattice) for ti in t]
    else
        error("Unknown type of solution provided to `get_lrs_vals`. Only ODE or Jump solutions are supported.")
    end
end

# Functions which in each sample point reshapes the vector of values to the correct form (depending
# on the type of lattice used).
function reshape_vals(vals, lattice::CartesianGridRej{N, T}) where {N,T}
    return reshape(vals, lattice.dims...)
end
function reshape_vals(vals, lattice::AbstractSparseArray{Bool, Int64, 1})
    return SparseVector(lattice.n, lattice.nzind, vals)
end
function reshape_vals(vals, lattice::AbstractSparseArray{Bool, Int64, 2})
    return SparseMatrixCSC(lattice.m, lattice.n, lattice.colptr, lattice.rowval, vals)
end
function reshape_vals(vals, lattice::DiGraph)
    return vals
end

