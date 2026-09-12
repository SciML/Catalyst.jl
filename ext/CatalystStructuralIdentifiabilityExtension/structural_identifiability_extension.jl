### Structural Identifiability ODE Creation ###

# For a reaction system, list of measured quantities and known parameters, generate a StructuralIdentifiability compatible ODE.
"""
make_si_ode(rs::ReactionSystem; measured_quantities=observed(rs), known_p = [], ignore_no_measured_warn=false, remove_conserved = true, mtkcompile = false)

Creates a reaction rate equation ODE system of the form used within the StructuralIdentifiability.jl package. The output system is compatible with all StructuralIdentifiability functions.

Arguments:
- `rs::ReactionSystem`; The reaction system we wish to convert to an ODE.
- `measured_quantities=[]`: The quantities of the system we can measure. May either be equations (e.g. `x1 + x2`), or single species (e.g. the symbolic `x`, `rs.x`, or the symbol `:x`).
- `known_p = []`: List of parameters for which their values are assumed to be known.
- `ignore_no_measured_warn = false`: If set to `true`, no warning is provided when the `measured_quantities` vector is empty.
- `remove_conserved = true`: Whether to eliminate conservation laws when computing the ode (this can reduce runtime of identifiability analysis significantly).
- `mtkcompile = false`: Whether to apply `ModelingToolkitBase.mtkcompile` to the generated ODE `System`. This is required for `ReactionSystem`s coupled with algebraic equations. Variables eliminated by `mtkcompile` (e.g. those defined by algebraic equations) may still be used in `measured_quantities` and `funcs_to_check`.

Example:
```julia
using Catalyst, StructuralIdentifiability
rs = @reaction_network begin
    (p,d), 0 <--> X
end
make_si_ode(rs; measured_quantities = [:X], known_p = [:p])
```

Notes:
- This function is part of the StructuralIdentifiability.jl extension. StructuralIdentifiability.jl must be imported to access it.
- `measured_quantities` and `known_p` input may also be symbolic (e.g. measured_quantities = [rs.X])
"""
function Catalyst.make_si_ode(rs::ReactionSystem; measured_quantities = [], known_p = [],
        ignore_no_measured_warn = false, remove_conserved = true, mtkcompile::Bool = false)
    # Creates a MTK ODE System, and a list of measured quantities (there are equations).
    # Gives these to SI to create an SI ode model of its preferred form.
    osys, conseqs, _, _ = make_osys(rs; remove_conserved, mtkcompile)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p,
        conseqs; ignore_no_measured_warn)
    return SI.mtk_to_si(osys, measured_quantities)[1]
end

### Structural Identifiability Wrappers ###

"""
assess_local_identifiability(rs::ReactionSystem, args...; measured_quantities = [], known_p = [], remove_conserved = true, ignore_no_measured_warn=false, mtkcompile = false, kwargs...)

Applies StructuralIdentifiability.jl's `assess_local_identifiability` function to the reaction rate equation ODE model for the given Catalyst `ReactionSystem`. Automatically converts the `ReactionSystem` to an appropriate ODE `System` and computes structural identifiability,

Arguments:
- `rs::ReactionSystem`; The reaction system for which we wish to compute structural identifiability of the associated reaction rate equation ODE model.
- `measured_quantities=[]`: The quantities of the system we can measure. May either be equations (e.g. `x1 + x2`), or single species (e.g. the symbolic `x`, `rs.s`, or the symbol `:x`).
- `known_p = []`: List of parameters which values are known.
- `ignore_no_measured_warn = false`: If set to `true`, no warning is provided when the `measured_quantities` vector is empty.
- `remove_conserved = true`: Whether to eliminate conservation laws when computing the ode (this can reduce runtime of identifiability analysis significantly).
- `mtkcompile = false`: Whether to apply `ModelingToolkitBase.mtkcompile` to the generated ODE `System`. This is required for `ReactionSystem`s coupled with algebraic equations. Variables eliminated by `mtkcompile` (e.g. those defined by algebraic equations) may still be used in `measured_quantities` and `funcs_to_check`.

Example:
```julia
using Catalyst, StructuralIdentifiability
rs = @reaction_network begin
    (p,d), 0 <--> X
end
assess_local_identifiability(rs; measured_quantities = [:X], known_p = [:p])
```

Notes:
- This function is part of the StructuralIdentifiability.jl extension. StructuralIdentifiability.jl must be imported to access it.
- `measured_quantities` and `known_p` input may also be symbolic (e.g. measured_quantities = [rs.X])
"""
function SI.assess_local_identifiability(rs::ReactionSystem, args...;
        measured_quantities = [], known_p = [], funcs_to_check = Vector(),
        remove_conserved = true, ignore_no_measured_warn = false,
        mtkcompile::Bool = false, kwargs...)
    # Creates an ODE System, list of measured quantities, and functions to check, of SI's preferred form.
    osys, conseqs, consconsts, vars = make_osys(rs; remove_conserved, mtkcompile)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p,
        conseqs; ignore_no_measured_warn)
    funcs_to_check, si_funcs_to_check = make_ftc(funcs_to_check, conseqs, vars)

    # Computes identifiability and converts it to a easy to read form.
    out = SI.assess_local_identifiability(osys, args...; measured_quantities,
        funcs_to_check = si_funcs_to_check, kwargs...)
    return make_output(out, funcs_to_check, si_funcs_to_check, consconsts)
end

"""
assess_identifiability(rs::ReactionSystem, args...; measured_quantities = [], known_p = [], remove_conserved = true, ignore_no_measured_warn=false, mtkcompile = false, kwargs...)

Applies StructuralIdentifiability.jl's `assess_identifiability` function to a Catalyst `ReactionSystem`. Internally it is converted to an ODE `System`, for which structural identifiability is computed.

Arguments:
- `rs::ReactionSystem`; The reaction system we wish to compute structural identifiability for.
- `measured_quantities=[]`: The quantities of the system we can measure. May either be equations (e.g. `x1 + x2`), or single species (e.g. the symbolic `x`, `rs.s`, or the symbol `:x`).
- `known_p = []`: List of parameters which values are known.
- `ignore_no_measured_warn = false`: If set to `true`, no warning is provided when the `measured_quantities` vector is empty.
- `remove_conserved = true`: Whether to eliminate conservation laws when computing the ode (this can reduce runtime of identifiability analysis significantly).
- `mtkcompile = false`: Whether to apply `ModelingToolkitBase.mtkcompile` to the generated ODE `System`. This is required for `ReactionSystem`s coupled with algebraic equations. Variables eliminated by `mtkcompile` (e.g. those defined by algebraic equations) may still be used in `measured_quantities` and `funcs_to_check`.

Example:
```julia
using Catalyst, StructuralIdentifiability
rs = @reaction_network begin
    (p,d), 0 <--> X
end
assess_identifiability(rs; measured_quantities = [:X], known_p = [:p])
```

Notes:
- This function is part of the StructuralIdentifiability.jl extension. StructuralIdentifiability.jl must be imported to access it.
- `measured_quantities` and `known_p` input may also be symbolic (e.g. measured_quantities = [rs.X])
"""
function SI.assess_identifiability(rs::ReactionSystem, args...;
        measured_quantities = [], known_p = [], funcs_to_check = Vector(),
        remove_conserved = true, ignore_no_measured_warn = false,
        mtkcompile::Bool = false, kwargs...)
    # Creates an ODE System, list of measured quantities, and functions to check, of SI's preferred form.
    osys, conseqs, consconsts, vars = make_osys(rs; remove_conserved, mtkcompile)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p,
        conseqs; ignore_no_measured_warn)
    funcs_to_check, si_funcs_to_check = make_ftc(funcs_to_check, conseqs, vars)

    # Computes identifiability and converts it to a easy to read form.
    # The `::System` designation fixes: https://github.com/SciML/StructuralIdentifiability.jl/issues/360,
    # however, the exact mechanisms of this is still not fully clear.
    out = SI.assess_identifiability(osys::System, args...; measured_quantities,
        funcs_to_check = si_funcs_to_check, kwargs...)
    return make_output(out, funcs_to_check, si_funcs_to_check, consconsts)
end

"""
find_identifiable_functions(rs::ReactionSystem, args...; measured_quantities = [], known_p = [], remove_conserved = true, ignore_no_measured_warn=false, mtkcompile = false, kwargs...)

Applies StructuralIdentifiability.jl's `find_identifiable_functions` function to a Catalyst `ReactionSystem`. Internally it is converted to an ODE `System`, for which structurally identifiable functions are computed.

Arguments:
- `rs::ReactionSystem`; The reaction system we wish to compute structural identifiability for.
- `measured_quantities=[]`: The quantities of the system we can measure. May either be equations (e.g. `x1 + x2`), or single species (e.g. the symbolic `x`, `rs.s`, or the symbol `:x`).
- `known_p = []`: List of parameters which values are known.
- `ignore_no_measured_warn = false`: If set to `true`, no warning is provided when the `measured_quantities` vector is empty.
- `remove_conserved = true`: Whether to eliminate conservation laws when computing the ode (this can reduce runtime of identifiability analysis significantly).
- `mtkcompile = false`: Whether to apply `ModelingToolkitBase.mtkcompile` to the generated ODE `System`. This is required for `ReactionSystem`s coupled with algebraic equations. Variables eliminated by `mtkcompile` (e.g. those defined by algebraic equations) may still be used in `measured_quantities` and `funcs_to_check`.

Example:
```julia
using Catalyst, StructuralIdentifiability
rs = @reaction_network begin
    (p,d), 0 <--> X
end
find_identifiable_functions(rs; measured_quantities = [:X], known_p = [:p])
```

Notes:
- This function is part of the StructuralIdentifiability.jl extension. StructuralIdentifiability.jl must be imported to access it.
- `measured_quantities` and `known_p` input may also be symbolic (e.g. measured_quantities = [rs.X])
"""
function SI.find_identifiable_functions(rs::ReactionSystem, args...;
        measured_quantities = [], known_p = [], remove_conserved = true,
        ignore_no_measured_warn = false, mtkcompile::Bool = false, kwargs...)
    # Creates an ODE System, and list of measured quantities, of SI's preferred form.
    osys, conseqs, consconsts, _ = make_osys(rs; remove_conserved, mtkcompile)
    measured_quantities = make_measured_quantities(rs, measured_quantities, known_p,
        conseqs; ignore_no_measured_warn)

    # Computes identifiable functions and converts it to a easy to read form.
    out = SI.find_identifiable_functions(osys, args...; measured_quantities, kwargs...)
    return vector_subs(out, consconsts)
end

### Helper Functions ###

# From a reaction system, creates the corresponding MTK-style ODE System for SI application
# Also compute the, later needed, conservation law equations and list of system symbols (unknowns and parameters).
function make_osys(rs::ReactionSystem; remove_conserved = true, mtkcompile::Bool = false)
    # Creates the ODE System corresponding to the ReactionSystem (expanding functions and flattening it).
    # Creates a list of the systems all symbols (unknowns and parameters).
    if !ModelingToolkitBase.iscomplete(rs)
        error("Identifiability should only be computed for complete systems. A ReactionSystem can be marked as complete using the `complete` function.")
    end
    rs = complete(Catalyst.expand_registered_functions(flatten(rs)))
    osys = ode_model(rs; remove_conserved)
    if mtkcompile
        osys = ModelingToolkitBase.mtkcompile(osys)
    elseif ModelingToolkitBase.has_alg_equations(rs)
        error("The input ReactionSystem has algebraic equations. This requires setting `mtkcompile = true`.")
    else
        osys = complete(osys)
    end
    vars = [unknowns(rs); parameters(rs)]

    # Computes equations for system conservation laws.
    # `conseqs` and `consconsts` must be `Vector{Pair{Any, Any}}`s also when empty.
    if remove_conserved
        conseqs = Pair{Any, Any}[ceq.lhs => ceq.rhs for ceq in conservedequations(rs)]
        consconsts = Pair{Any, Any}[cconst.lhs => cconst.rhs
                                    for cconst in conservationlaw_constants(rs)]
    else
        conseqs = Pair{Any, Any}[]
        consconsts = Pair{Any, Any}[]
    end

    # Variables eliminated by `mtkcompile` (e.g. algebraic variables) are no longer unknowns of
    # `osys`, so any reference to them in measured quantities or functions to check must be
    # replaced by the corresponding observed expression.
    mtkcompile && append!(conseqs, observed_subs(osys, conseqs))

    return osys, conseqs, consconsts, vars
end

# Creates substitution rules for the observed variables of `osys` (skipping those already
# covered by `conseqs`). Observed equations may depend on each other, so rules are expanded
# until they only contain unknowns and parameters of `osys`.
function observed_subs(osys, conseqs)
    obs = Pair{Any, Any}[eq.lhs => eq.rhs
                         for eq in observed(osys)
                         if !any(ceq -> isequal(ceq[1], eq.lhs), conseqs)]
    for _ in 1:length(obs)
        expanded = Pair{Any, Any}[lhs => substitute(rhs, obs) for (lhs, rhs) in obs]
        all(isequal(e[2], o[2]) for (e, o) in zip(expanded, obs)) && break
        obs = expanded
    end
    return obs
end

# Creates a list of measured quantities of a form that SI can read.
# Each measured quantity must have a form like:
# `obs_var ~ X` # (Here, `obs_var` is a variable, and X is whatever we can measure).
function make_measured_quantities(rs::ReactionSystem, measured_quantities::Vector{T},
        known_p::Vector{S}, conseqs; ignore_no_measured_warn = false) where {T, S}
    # Warning if the user didn't give any measured quantities.
    if !ignore_no_measured_warn && isempty(measured_quantities)
        @warn "No measured quantity provided to the `measured_quantities` argument, any further identifiability analysis will likely fail. You can disable this warning by setting `ignore_no_measured_warn = true`."
    end

    # Appends the known parameters to the measured_quantities vector. Converts any Symbols to symbolics.
    mqiterator = Iterators.flatten((measured_quantities, known_p))
    mqs = [(q isa Symbol) ? Catalyst._symbol_to_var(rs, q) : q for q in mqiterator]
    mqs = vector_subs(mqs, conseqs)

    # Creates one internal observation variable for each measured quantity (`___internal_observables`).
    # Creates a vector of equations, setting each measured quantity equal to one observation variable.
    @variables (___internal_observables(Catalyst.get_iv(rs)))[1:length(mqs)]
    return Equation[(q isa Equation) ? q : (___internal_observables[i] ~ q)
                    for (i, q) in enumerate(mqs)]
end

# Creates the functions that we wish to check for identifiability.
# If no `funcs_to_check` are given, defaults to checking identifiability for all unknowns and parameters.
# Returns both the functions as requested by the user and the versions passed to SI, in which
# conserved quantities (and other eliminated variables) have been substituted away.
# E.g. for `X1 <--> X2`, replaces `X2` with `Γ[1] - X2`.
# Removing conserved quantities makes SI's algorithms much more performant.
function make_ftc(funcs_to_check, conseqs, vars)
    isempty(funcs_to_check) && (funcs_to_check = vars)
    return funcs_to_check, vector_subs(funcs_to_check, conseqs)
end

# Processes the outputs to a better form.
# SI keys its output by the functions it was given (`si_funcs_to_check`). These are mapped back
# to the user-requested functions (`funcs_to_check`), so that e.g. Γ or eliminated variables
# are not displayed, and sorted according to their input order.
function make_output(out, funcs_to_check, si_funcs_to_check, consconsts)
    si_keys = vector_subs(si_funcs_to_check, consconsts)
    out_keys = vector_subs(keys(out), consconsts)
    key_to_idx = Dict(k => i for (i, k) in enumerate(si_keys))
    return OrderedDict(funcs_to_check[key_to_idx[k]] => v
    for (k, v) in zip(out_keys, values(out)))
end

# For a vector of expressions and a conservation law, substitutes the law into every equation.
vector_subs(eqs, subs) = [substitute(eq, subs) for eq in eqs]
