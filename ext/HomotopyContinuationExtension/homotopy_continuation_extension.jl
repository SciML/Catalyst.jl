### Homotopy Continuation Based Steady State Finding ###

"""
    hc_steady_states(rs::ReactionSystem, ps; filter_negative=true, neg_thres=-1e-20, u0=typeof(ps)(), kwargs...)

Uses homotopy continuation to find the steady states of the ODE system corresponding to the provided reaction system.

Arguments:
- `rs::ReactionSystem`: The reaction system for which we want to find the steady states.
- `ps`: The parameter values for which we want to find the steady states.
- `filter_negative=true`: If set to true, solutions with any species concentration <0 is removed from the output.
- `neg_thres=-1e-20`: Determine the minimum values for which a species concentration is to be considred non-negative. Species conentrations ``> neg_thres`` but `< 0.0` are set to `0.0`.
- `u0=typeof(ps)()`: Initial conditions for which we want to find the steady states. For systems with conservation laws this are required to compute conserved quantities.
- `kwargs...`: any additional arguments (like `show_progress= true`) are passed into HomotopyContinuation.jl's `solve` call. 

Examples:
```@repl
rs = @reaction_network begin
    k1, Y --> 2X
    k2, 2X --> X + Y
    k3, X + Y --> Y
    k4, X --> 0
end
ps = [:k3 => 1.0, :k2 => 2.0, :k4 => 1.5, :k1=>8.0]
hc_sol = hc_steady_states(rs, ps)
```
gives
```
[0.5000000000000002, 2.0000000000000004]
[0.0, 0.0]
[4.499999999999999, 5.999999999999999]
```

Notes:
- This is a wrapper around the `solve` function provided by HomotopyContinuation.jl, all credit for this functionality to that package's authors.
```
  """
function Catalyst.hc_steady_states(rs::ReactionSystem, ps; filter_negative=true, neg_thres=-1e-20, u0=typeof(ps)(), kwargs...)
    ss_poly = steady_state_polynomial(rs, ps, u0)
    sols = HC.real_solutions(HC.solve(ss_poly; kwargs...))
    reorder_sols!(sols, ss_poly, rs)
    return (filter_negative ? filter_negative_f(sols; neg_thres=neg_thres) : sols)
end

# For a given reaction system, paraemter values, and initial conditions, find the polynomial that HC solves to find steady states.
function steady_state_polynomial(rs::ReactionSystem, ps, u0)
    ns = convert(NonlinearSystem, rs; remove_conserved = true)
    pre_varmap = map(pair -> pair::Pair{Num,Float64}, [symmap_to_varmap(rs,u0); symmap_to_varmap(rs,ps)]) # Needed in case ps or u0 are empty (making the combination a Vector{Any}, causing an error).
    conservationlaw_errorcheck(rs, pre_varmap)
    p_vals = MT.varmap_to_vars(pre_varmap, parameters(ns); defaults = MT.defaults(ns))
    p_dict = Dict(parameters(ns) .=> p_vals)
    eqs = vcat(equations(ns), conservedequations(rs))
    eqs = map(eq -> substitute(eq.rhs - eq.lhs, p_dict), eqs)
    return Catalyst.to_multivariate_poly(eqs)
end

# If u0s are not given while conservation laws are present, throws an error.
function conservationlaw_errorcheck(rs, pre_varmap)
    vars_with_vals = union(first.(pre_varmap), keys(ModelingToolkit.defaults(rs)))
    isempty(intersect(species(rs), vars_with_vals)) || return
    isempty(conservedequations(rs)) && return 
    error("The system have conservation laws but no initial conditions were provided. Please provide initial conditions.")
end

# HC orders the solution vector according to the lexiographic values of the variable names. This reorders the output acording to the species index in the reaction system species vector.
function reorder_sols!(sols, ss_poly, rs::ReactionSystem)
    var_names_extended = String.(Symbol.(HC.variables(ss_poly)))
    var_names = [Symbol(s[1:prevind(s, findlast('_', s))]) for s in var_names_extended]
    sort_pattern = indexin(MT.getname.(species(rs)), var_names)
    foreach(sol -> permute!(sol, sort_pattern), sols)
end

# Filters away solutions with negative species concentrations (and for neg_thres < val < 0.0, sets val=0.0).
function filter_negative_f(sols; neg_thres=-1e-20)
    for sol in sols, idx in 1:length(sol)
        (neg_thres < sol[idx] < 0.0) && (sol[idx] = 0.0)
    end
    return filter(sol -> all(>=(0.0), sol), sols)
end