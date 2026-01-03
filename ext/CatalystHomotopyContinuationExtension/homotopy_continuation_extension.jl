### Homotopy Continuation Based Steady State Finding ###

"""
    hc_steady_states(rs::ReactionSystem, ps; filter_negative = true, neg_thres = -1e-16, u0 = typeof(ps)(), kwargs...)

Uses homotopy continuation via HomotopyContinuation.jl to find the steady states of the ODE system corresponding to the provided reaction system.

Arguments:
- `rs::ReactionSystem`: The reaction system for which we want to find the steady states.
- `ps`: The parameter values for which we want to find the steady states.
- `filter_negative = true`: If set to true, solutions with any species concentration <neg_thres is removed from the output.
- `neg_thres = -1e-15`: Determine the minimum values for which a species concentration is to be considered non-negative. Species concentrations ``> neg_thres`` but `< 0.0` are set to `0.0`.
- `u0 = nothing`: Initial conditions for which we want to find the steady states. For systems with conservation laws this are required to compute conserved quantities. Initial conditions are not required for all species, only those involved in conserved quantities (if this set is unknown, it is recommended to provide initial conditions for all species).
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
function Catalyst.hc_steady_states(rs::ReactionSystem, ps; filter_negative = true,
        neg_thres = -1e-15, u0 = [], kwargs...)
    if !isautonomous(rs)
        error("Attempting to compute steady state for a non-autonomous system (e.g. where some rate depend on $(get_iv(rs))). This is not possible.")
    end
    ss_poly = steady_state_polynomial(rs, ps, u0)
    sols = HC.real_solutions(HC.solve(ss_poly; kwargs...))
    reorder_sols!(sols, ss_poly, rs)
    return (filter_negative ? filter_negative_f(sols; neg_thres) : sols)
end

# For a given reaction system, parameter values, and initial conditions, find the polynomial that HC solves to find steady states.
function steady_state_polynomial(rs::ReactionSystem, ps, u0)
    # Creates the appropriate nonlinear system, and converts parameters to a form that can
    # be substituted in later.
    rs = Catalyst.expand_registered_functions(rs)
    ns = complete(make_rre_algeqs(rs; remove_conserved = true, conseqs_remake_warn = false))
    pre_varmap = [symmap_to_varmap(rs, u0)..., symmap_to_varmap(rs, ps)...]
    Catalyst.conservationlaw_errorcheck(rs, pre_varmap)
    p_dict = make_p_val_dict(pre_varmap, rs, ns)

    # Step by step convert the equation to something HC can work on (adds conserved equations,
    # inserts parameter values and put everything on a side, converts e.g. 2.0 to 2, remove
    # denominators to avoid rational polynomials).
    eqs = vcat(equations(ns), conservedequations(rs))
    eqs = [substitute(eq.rhs - eq.lhs, p_dict) for eq in eqs]
    eqs = make_int_exps.(eqs)
    eqs = [remove_denominators(common_denominator(eq)) for eq in eqs]
    ss_poly = Catalyst.to_multivariate_poly(eqs)
    return poly_type_convert(ss_poly)
end

# Function for making a map from parameter values (including conservation law constants) to
# their values. Previously a simple call, but made complicated due to various MTK updates so
# not all delegated to a specific function.
function make_p_val_dict(pre_varmap, rs, ns)
    # Creates a parameter vector with conservation parameters expanded (i.e. Γ[1], Γ[2], not `Γ`).
    ps_no_cons = filter(p -> !Catalyst.isconserved(p), parameters(ns))
    ps_cons_expanded = getfield.(conservationlaw_constants(rs), :lhs)
    ps = [ps_no_cons; ps_cons_expanded]

    # Creates a dictionary with the correct default values/equations (again expanding `Γ` to Γ[1], Γ[2]).
    defaults = MT.initial_conditions(ns)
    filter!(p -> !Catalyst.isconserved(p[1]), defaults)
    foreach(conseq -> defaults[conseq.lhs] = conseq.rhs, conservationlaw_constants(rs))

    # Creates and return the full parameter value dictionary.p_vals = MT.varmap_to_vars(pre_varmap, all_ps; defaults = def_dict)
    p_vals = varmap_to_vars_mtkv9(pre_varmap, ps; defaults)
    return Dict(ps .=> p_vals)
end

# Parses and expression and return a version where any exponents that are Float64 (but an int, like 2.0) are turned into Int64s.
function make_int_exps(expr)
    wrap(Rewriters.Postwalk(Rewriters.PassThrough(___make_int_exps))(unwrap(expr))).val
end
function ___make_int_exps(expr)
    !iscall(expr) && return expr
    if (operation(expr) == ^)
        if isinteger(sorted_arguments(expr)[2])
            return sorted_arguments(expr)[1]^Int64(sorted_arguments(expr)[2])
        else
            error("An non integer ($(sorted_arguments(expr)[2])) was found as a variable exponent. Non-integer exponents are not supported for homotopy continuation based steady state finding.")
        end
    end
end

# Converts an expression of the form p(x)/q(x) + r(x)/s(x) to t(x)/u(x) (i.e. puts everything
# above a single denominator, which is what `remove_denominators` is able to simplify away).
function common_denominator(expr)
    iscall(expr) || return expr
    if operation(expr) == /
        num, den = arguments(expr)
        num = common_denominator(num)
        den = common_denominator(den)
        return num / den
    end
    if operation(expr) == +
        num = 0
        den = 1
        for arg in arguments(expr)
            arg = common_denominator(arg)
            if iscall(arg) && operation(arg) == /
                n, d = arguments(arg)
            else
                n = arg
                d = 1
            end
            num = num * d + den * n
            den *= d
        end
        return num / den
    end
    if operation(expr) == ^
        base, pow = arguments(expr)
        base = common_denominator(base)
        if iscall(base) && operation(base) == /
            num, den = arguments(base)
        else
            num, den = base, 1
        end
        num ^= pow
        den ^= pow
        return num / den
    end
    return maketerm(BasicSymbolic, operation(expr), map(common_denominator, arguments(expr)), metadata(expr))
end

# If the input is a fraction, removes the denominator.
function remove_denominators(expr)
    s_expr = simplify_fractions(expr)
    !iscall(expr) && return expr
    if operation(s_expr) == /
        return remove_denominators(sorted_arguments(s_expr)[1])
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
    sort_pattern = indexin(MT.getname.(unknowns(rs)), var_names)
    foreach(sol -> permute!(sol, sort_pattern), sols)
end

# Filters away solutions with negative species concentrations (and for neg_thres < val < 0.0, sets val=0.0).
function filter_negative_f(sols; neg_thres = -1e-15)
    for sol in sols, idx in 1:length(sol)

        (neg_thres < sol[idx] < 0) && (sol[idx] = 0)
    end
    return filter(sol -> all(>=(0), sol), sols)
end

# Sometimes (when polynomials are created from coupled CRN/DAEs), the steady state polynomial have the wrong type.
# This converts it to the correct type, which homotopy continuation can handle.
const WRONG_POLY_TYPE = Vector{DynamicPolynomials.Polynomial{
    DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder},
    DynamicPolynomials.Graded{DynamicPolynomials.LexOrder}}}
const CORRECT_POLY_TYPE = Vector{DynamicPolynomials.Polynomial{
    DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder},
    DynamicPolynomials.Graded{DynamicPolynomials.LexOrder}, Float64}}
function poly_type_convert(ss_poly)
    (typeof(ss_poly) == WRONG_POLY_TYPE) && return convert(CORRECT_POLY_TYPE, ss_poly)
    return ss_poly
end



### SAVED ARCHIVED MTK FUNCTION - REMOVE SOME TIME ###
# pre-v10 version of function
function varmap_to_vars_mtkv9(varmap, varlist; defaults = Dict(), check = true,
        toterm = MT.default_toterm, promotetoconcrete = nothing,
        tofloat = true, use_union = true)
    varlist = collect(map(unwrap, varlist))

    # Edge cases where one of the arguments is effectively empty.
    is_incomplete_initialization = varmap isa DiffEqBase.NullParameters ||
                                   varmap === nothing
    if is_incomplete_initialization || isempty(varmap)
        if isempty(defaults)
            if !is_incomplete_initialization && check
                isempty(varlist) || throw(MT.MissingVariablesError(varlist))
            end
            return nothing
        else
            varmap = Dict()
        end
    end

    # We respect the input type if it's a static array
    # otherwise canonicalize to a normal array
    # container_type = T <: Union{Dict,Tuple} ? Array : T
    if varmap isa MT.StaticArray
        container_type = typeof(varmap)
    else
        container_type = Array
    end

    vals = if eltype(varmap) <: Pair # `varmap` is a dict or an array of pairs
        varmap = MT.todict(varmap)
        _varmap_to_vars_mtkv9(varmap, varlist; defaults = defaults, check = check,
            toterm = toterm)
    else # plain array-like initialization
        varmap
    end

    promotetoconcrete === nothing && (promotetoconcrete = container_type <: AbstractArray)
    if promotetoconcrete
        vals = MT.promote_to_concrete(vals; tofloat = tofloat, use_union = use_union)
    end

    if isempty(vals)
        return nothing
    elseif container_type <: Tuple
        (vals...,)
    else
        SymbolicUtils.Code.create_array(container_type, eltype(vals), Val{1}(),
            Val(length(vals)), vals...)
    end
end

function _varmap_to_vars_mtkv9(varmap::Dict, varlist; defaults = Dict(), check = false,
        toterm = Symbolics.diff2term, initialization_phase = false)
    varmap = canonicalize_varmap_mtkv9(varmap; toterm)
    defaults = canonicalize_varmap_mtkv9(defaults; toterm)
    varmap = merge(defaults, varmap)
    values = Dict()

    T = Union{}
    for var in varlist
        var = MT.unwrap(var)
        val = MT.unwrap(MT.fixpoint_sub(var, varmap; operator = Symbolics.Operator))
        if !isequal(val, var)
            values[var] = val
        end
    end
    missingvars = setdiff(varlist, collect(keys(values)))
    check && (isempty(missingvars) || throw(MT.MissingVariablesError(missingvars)))
    return [values[MT.unwrap(var)] for var in varlist]
end

function canonicalize_varmap_mtkv9(varmap; toterm = Symbolics.diff2term)
    new_varmap = Dict()
    for (k, v) in varmap
        k = MT.unwrap(k)
        v = MT.unwrap(v)
        new_varmap[k] = v
        new_varmap[toterm(k)] = v
        if Symbolics.isarraysymbolic(k) && Symbolics.shape(k) !== Symbolics.Unknown()
            for i in eachindex(k)
                new_varmap[k[i]] = v[i]
                new_varmap[toterm(k[i])] = v[i]
            end
        end
    end
    return new_varmap
end
