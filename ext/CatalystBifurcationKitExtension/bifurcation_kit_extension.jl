### Dispatch for BifurcationKit BifurcationProblems ###

# Creates a BifurcationProblem, using a ReactionSystem as an input.
function BK.BifurcationProblem(rs::ReactionSystem, u0_bif, ps, bif_par, args...;
                               plot_var=nothing, record_from_solution=BK.record_sol_default, jac=true, u0=[], kwargs...)
    if !is_autonomous(rs) 
        error("Attempting to create a `BifurcationProblem` for a non-autonomous system (e.g. where some rate depend on $(rs.iv)). This is not possible.")
    end

    # Converts symbols to symbolics.
    (bif_par isa Symbol) && (bif_par = ModelingToolkit.get_var_to_name(rs)[bif_par])
    (plot_var isa Symbol) && (plot_var = ModelingToolkit.get_var_to_name(rs)[plot_var])
    ((u0_bif isa Vector{<:Pair{Symbol,<:Any}}) || (u0_bif isa Dict{Symbol, <:Any})) && (u0_bif = symmap_to_varmap(rs, u0_bif))
    ((ps isa Vector{<:Pair{Symbol,<:Any}}) || (ps isa Dict{Symbol, <:Any})) && (ps = symmap_to_varmap(rs, ps))
    ((u0 isa Vector{<:Pair{Symbol,<:Any}}) || (u0 isa Dict{Symbol, <:Any})) && (u0 = symmap_to_varmap(rs, u0))

    # Creates NonlinearSystem.
    Catalyst.conservationlaw_errorcheck(rs, vcat(ps, u0))
    nsys = complete(convert(NonlinearSystem, rs; remove_conserved=true, defaults=Dict(u0)))

    # Makes BifurcationProblem (this call goes through the ModelingToolkit-based BifurcationKit extension).
    return BK.BifurcationProblem(nsys, u0_bif, ps, bif_par, args...; plot_var=plot_var,
                                 record_from_solution=record_from_solution, jac=jac, kwargs...)
end