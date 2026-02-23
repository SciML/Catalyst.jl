### Dispatch for BifurcationKit BifurcationProblems ###

# Creates a BifurcationProblem, using a ReactionSystem as an input.
function BK.BifurcationProblem(rs::ReactionSystem, u0_bif, ps, bif_par, args...;
        plot_var = nothing, record_from_solution = BK.record_sol_default, jac = true, u0 = [], kwargs...)
    if !isautonomous(rs)
        error("Attempting to create a `BifurcationProblem` for a non-autonomous system (e.g. where some rate depend on $(get_iv(rs))). This is not possible.")
    end

    # Converts symbols to symbolics.
    (bif_par isa Symbol) && (bif_par = ModelingToolkitBase.get_var_to_name(rs)[bif_par])
    (plot_var isa Symbol) && (plot_var = ModelingToolkitBase.get_var_to_name(rs)[plot_var])
    if (u0_bif isa Vector{<:Pair{Symbol, <:Any}}) || (u0_bif isa Dict{Symbol, <:Any})
        u0_bif = Catalyst.symmap_to_varmap(rs, u0_bif)
    end
    if (ps isa Vector{<:Pair{Symbol, <:Any}}) || (ps isa Dict{Symbol, <:Any})
        ps = Catalyst.symmap_to_varmap(rs, ps)
    end
    if (u0 isa Vector{<:Pair{Symbol, <:Any}}) || (u0 isa Dict{Symbol, <:Any})
        u0 = Catalyst.symmap_to_varmap(rs, u0)
    end

    # Creates nonlinear System. If there are conservation laws, these are manually added as 
    # equations, and the conservationlaw parameter values are manually computed and added to `ps`.
    Catalyst.conservationlaw_errorcheck(rs, vcat(ps, u0))
    nsys = if Catalyst.num_cons_laws(rs) == 0
        complete(ss_ode_model(rs))
    else
        Γ_vals = Catalyst.get_networkproperties(rs).conservedconst => 
            [Symbolics.substitute(ceq.rhs, u0) for ceq in conservationlaw_constants(rs)]
        ps = [ps; Γ_vals]
        complete(ss_ode_model(rs; remove_conserved = true, include_cl_as_eqs = true))
    end

    # Makes BifurcationProblem (this call goes through the ModelingToolkit-based BifurcationKit extension).
    return BK.BifurcationProblem(nsys, u0_bif, ps, bif_par, args...; plot_var,
        record_from_solution, jac, kwargs...)
end