### Dispatch for BifurcationKit BifurcationProblems ###

# Creates a BifurcationProblem, using a ReactionSystem as an input.
function BK.BifurcationProblem(rs::ReactionSystem, u0_bif, ps, bif_par, args...; plot_var=nothing, record_from_solution=BK.record_sol_default, jac=true, u0=[], kwargs...)
    
    println("HERE")

    # Converts symbols to symbolics.
    (bif_par isa Symbol) && (bif_par = rs.var_to_name[bif_par])
    (plot_var isa Symbol) && (plot_var = rs.var_to_name[plot_var])
    (u0_bif isa Vector{Symbol}) && (u0_bif = symmap_to_varmap(rs, u0_bif))
    (ps isa Vector{Symbol}) && (ps = symmap_to_varmap(rs, ps))

    println("HERE")

    # Creates NonlinearSystem.
    conservationlaw_errorcheck(rs, vcat(ps, u0))
    nsys = covnert(NonlinearSystem, rs; defaults=u0)

    # Makes BifurcationProblem.
    return BK.BifurcationProblem(nsys, u0_bif, ps, bif_par, args...; plot_var=plot_var, record_from_solution=record_from_solution, jac=jac, kwargs...)
end

