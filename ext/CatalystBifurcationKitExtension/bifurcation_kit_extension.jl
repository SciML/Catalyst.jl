### Dispatch for BifurcationKit BifurcationProblems ###

# Creates a BifurcationProblem, using a ReactionSystem as an input.
function BK.BifurcationProblem(rs::ReactionSystem, u0_bif, ps, bif_par, args...;
        plot_var = nothing, record_from_solution = BK.record_sol_default, jac = true, u0 = [], kwargs...)
    if !isautonomous(rs)
        error("Attempting to create a `BifurcationProblem` for a non-autonomous system (e.g. where some rate depend on $(get_iv(rs))). This is not possible.")
    end

    # Converts symbols to symbolics.
    (bif_par isa Symbol) && (bif_par = ModelingToolkit.get_var_to_name(rs)[bif_par])
    (plot_var isa Symbol) && (plot_var = ModelingToolkit.get_var_to_name(rs)[plot_var])
    if (u0_bif isa Vector{<:Pair{Symbol, <:Any}}) || (u0_bif isa Dict{Symbol, <:Any})
        u0_bif = symmap_to_varmap(rs, u0_bif)
    end
    if (ps isa Vector{<:Pair{Symbol, <:Any}}) || (ps isa Dict{Symbol, <:Any})
        ps = symmap_to_varmap(rs, ps)
    end
    if (u0 isa Vector{<:Pair{Symbol, <:Any}}) || (u0 isa Dict{Symbol, <:Any})
        u0 = symmap_to_varmap(rs, u0)
    end

    # Creates NonlinearSystem.
    Catalyst.conservationlaw_errorcheck(rs, vcat(ps, u0))
    nsys = bkext_make_nsys(rs, u0)

    # Makes BifurcationProblem (this call goes through the ModelingToolkit-based BifurcationKit extension).
    return BK.BifurcationProblem(nsys, u0_bif, ps, bif_par, args...; plot_var,
        record_from_solution, jac, kwargs...)
end

# Creates the NonlinearSystem for the bifurcation problem. Used to be straightforward, but MTK
# updates have made handling of conservation laws more complicated, so we now have to do
# more things here.
function bkext_make_nsys(rs, u0)
    cons_eqs = conservationlaw_constants(rs)
    cons_default = [cons_eq.rhs for cons_eq in cons_eqs]
    cons_default = Catalyst.get_networkproperties(rs).conservedconst => cons_default
    defaults = Dict([u0; cons_default])
    nsys = make_rre_algeqs(rs; defaults, remove_conserved = true, conseqs_remake_warn = false)
    return complete(nsys)
end
