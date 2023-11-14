### Homotopy Continuation Based Steady State Finding ###

"""
    hc_steady_states(rs::ReactionSystem, ps; filter_negative=true, neg_thres=-1e-20, u0=typeof(ps)(), kwargs...)

Uses homotopy continuation via HomotopyContinuation.jl to find the steady states of the ODE system corresponding to the provided reaction system.

Arguments:
- `rs::ReactionSystem`: The reaction system for which we want to find the steady states.
- `ps`: The parameter values for which we want to find the steady states.
- `filter_negative=true`: If set to true, solutions with any species concentration <neg_thres is removed from the output.
- `neg_thres=-1e-20`: Determine the minimum values for which a species concentration is to be considered non-negative. Species concentrations ``> neg_thres`` but `< 0.0` are set to `0.0`.
- `u0=nothing`: Initial conditions for which we want to find the steady states. For systems with conservation laws this are required to compute conserved quantities. Initial conditions are not required for all species, only those involved in conserved quantities (if this set is unknown, it is recommended to provide initial conditions for all species). 
- `kwargs...`: any additional arguments (like `show_progress= true`) are passed into HomotopyContinuation.jl's `solve` call. 

Examples
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
- Homotopy-based steady state finding only works when all rates are rational polynomials (e.g. constant, linear, mm, or hill functions).
```
  """
function Catalyst.hc_steady_states(rs::ReactionSystem, ps; filter_negative=true, neg_thres=-1e-20, u0=[], kwargs...)
    ss_poly = steady_state_polynomial(rs, ps, u0)
    sols = HC.real_solutions(HC.solve(ss_poly; kwargs...))
    reorder_sols!(sols, ss_poly, rs)
    return (filter_negative ? filter_negative_f(sols; neg_thres) : sols)
end

# For a given reaction system, parameter values, and initial conditions, find the polynomial that HC solves to find steady states.
function steady_state_polynomial(rs::ReactionSystem, ps, u0)
    rs = expand_registered_functions(rs)
    ns = convert(NonlinearSystem, rs; remove_conserved = true)
    pre_varmap = [symmap_to_varmap(rs,u0)..., symmap_to_varmap(rs,ps)...]
    conservationlaw_errorcheck(rs, pre_varmap)
    p_vals = ModelingToolkit.varmap_to_vars(pre_varmap, parameters(ns); defaults = ModelingToolkit.defaults(ns))
    p_dict  = Dict(parameters(ns) .=> p_vals)
    eqs_pars_funcs = vcat(equations(ns), conservedequations(rs))
    eqs = map(eq -> substitute(eq.rhs - eq.lhs, p_dict), eqs_pars_funcs)
    eqs_intexp = make_int_exps.(eqs)
    return Catalyst.to_multivariate_poly(remove_denominators.(eqs_intexp))
end

# If u0s are not given while conservation laws are present, throws an error.
function conservationlaw_errorcheck(rs, pre_varmap)
    vars_with_vals = Set(p[1] for p in pre_varmap)
    any(s -> s in vars_with_vals, species(rs)) && return
    isempty(conservedequations(rs)) || 
        error("The system has conservation laws but initial conditions were not provided for some species.")
end

# Parses and expression and return a version where any exponents that are Float64 (but an int, like 2.0) are turned into Int64s.
make_int_exps(expr) = wrap(Rewriters.Postwalk(Rewriters.PassThrough(___make_int_exps))(unwrap(expr))).val
function ___make_int_exps(expr)
    !istree(expr) && return expr
    if (operation(expr) == ^) 
        if isinteger(arguments(expr)[2])
            return arguments(expr)[1] ^ Int64(arguments(expr)[2])
        else
            error("An non integer ($(arguments(expr)[2])) was found as a variable exponent. Non-integer exponents are not supported for homotopy continuation based steady state finding.")
        end
    end
end

# If the input is a fraction, removes the denominator.
function remove_denominators(expr)
    s_expr = simplify_fractions(expr)
    !istree(expr) && return expr
    if operation(s_expr) == /
        return remove_denominators(arguments(s_expr)[1])
    end
    if operation(s_expr) == +
        return sum(remove_denominators(arg) for arg in arguments(s_expr))
    end
    return s_expr
end

# HC orders the solution vector according to the lexicographic values of the variable names. This reorders the output according to the species index in the reaction system species vector.
function reorder_sols!(sols, ss_poly, rs::ReactionSystem)
    var_names_extended = String.(Symbol.(HC.variables(ss_poly)))
    var_names = [Symbol(s[1:prevind(s, findlast('_', s))]) for s in var_names_extended]
    sort_pattern = indexin(MT.getname.(species(rs)), var_names)
    foreach(sol -> permute!(sol, sort_pattern), sols)
end

# Filters away solutions with negative species concentrations (and for neg_thres < val < 0.0, sets val=0.0).
function filter_negative_f(sols; neg_thres=-1e-20)
    for sol in sols, idx in 1:length(sol)
        (neg_thres < sol[idx] < 0) && (sol[idx] = 0)
    end
    return filter(sol -> all(>=(0), sol), sols)
end