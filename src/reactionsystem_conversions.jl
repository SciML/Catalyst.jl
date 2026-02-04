### ODE & SDE Assembly ###

"""
    oderatelaw(rx; combinatoric_ratelaw=true)

Given a [`Reaction`](@ref), return the symbolic reaction rate law used in
generated ODEs for the reaction. Note, for a reaction defined by

`k*X*Y, X+Z --> 2X + Y`

the expression that is returned will be `k*X(t)^2*Y(t)*Z(t)`. For a reaction of
the form

`k, 2X+3Y --> Z`

the expression that is returned will be `k * (X(t)^2/2) * (Y(t)^3/6)`.

Notes:
- Allocates
- `combinatoric_ratelaw=true` uses factorial scaling factors in calculating the
    rate law, i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. If
    `combinatoric_ratelaw=false` then the ratelaw is `k*S^2`, i.e. the scaling
    factor is ignored.
"""
function oderatelaw(rx; combinatoric_ratelaw = true, expand_catalyst_funs = true)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate
    expand_catalyst_funs && (rl = expand_registered_functions(rl))

    # if the stoichiometric coefficients are not integers error if asking to scale rates
    !all(s -> s isa Union{Integer, SymbolicT}, substoich) &&
        (combinatoric_ratelaw == true) &&
        error("Non-integer stoichiometric coefficients require the combinatoric_ratelaw=false keyword to oderatelaw, or passing combinatoric_ratelaws=false to convert or ODEProblem.")

    if !only_use_rate
        coef = eltype(substoich) <: Number ? one(eltype(substoich)) : 1
        for (i, stoich) in enumerate(substoich)
            combinatoric_ratelaw && (coef *= factorial(stoich))
            rl *= isequal(stoich, one(stoich)) ? substrates[i] : substrates[i]^stoich
        end
        combinatoric_ratelaw && (!isequal(coef, one(coef))) && (rl /= coef)
    end
    rl
end

# Function returning `true` for species which shouldn't change from the reactions,
# including non-species variables.
drop_dynamics(s) = isconstant(s) || isbc(s) || (!isspecies(s))

# Compute signed stoichiometry term: stoich * expr, optimized for common cases.
# Used in both ODE RHS assembly and noise coefficient computation.
function _signed_stoich_term(stoich, expr)
    if stoich isa SymbolicT
        stoich * expr
    else
        signed_expr = (stoich > zero(stoich)) ? expr : -expr
        isone(abs(stoich)) ? signed_expr : stoich * expr
    end
end

function assemble_oderhs(rs, ispcs; combinatoric_ratelaws = true, remove_conserved = false,
        physical_scales = nothing, expand_catalyst_funs = true)
    nps = get_networkproperties(rs)
    species_to_idx = Dict(x => i for (i, x) in enumerate(ispcs))
    rhsvec = Any[0 for _ in ispcs]
    depspec_submap = if remove_conserved
        Dict(eq.lhs => eq.rhs for eq in nps.conservedeqs)
    else
        Dict()
    end

    for (rxidx, rx) in enumerate(get_rxs(rs))
        # check this reaction should be treated as an ODE
        !((physical_scales === nothing) ||
            (physical_scales[rxidx] == PhysicalScale.ODE)) && continue

        rl = oderatelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws,
            expand_catalyst_funs)
        remove_conserved && (rl = substitute(rl, depspec_submap))
        for (spec, stoich) in rx.netstoich
            # dependent species don't get an ODE, so are skipped
            remove_conserved && (spec in nps.depspecs) && continue

            # constant or BC species also do not get equations
            drop_dynamics(spec) && continue

            i = species_to_idx[spec]
            if _iszero(rhsvec[i])
                if stoich isa SymbolicT
                    rhsvec[i] = stoich * rl
                else
                    signedrl = (stoich > zero(stoich)) ? rl : -rl
                    rhsvec[i] = isone(abs(stoich)) ? signedrl : stoich * rl
                end
            else
                if stoich isa SymbolicT
                    rhsvec[i] += stoich * rl
                else
                    Δspec = isone(abs(stoich)) ? rl : abs(stoich) * rl
                    rhsvec[i] = (stoich > zero(stoich)) ? (rhsvec[i] + Δspec) :
                                (rhsvec[i] - Δspec)
                end
            end
        end
    end

    rhsvec
end

function assemble_drift(rs, ispcs; combinatoric_ratelaws = true, as_odes = true,
        include_zero_odes = true, remove_conserved = false, physical_scales = nothing,
        expand_catalyst_funs = true)
    rhsvec = assemble_oderhs(rs, ispcs; combinatoric_ratelaws, remove_conserved,
        physical_scales, expand_catalyst_funs)
    if as_odes
        D = Differential(get_iv(rs))
        eqs = [Equation(D(x), rhs)
               for (x, rhs) in zip(ispcs, rhsvec) if (include_zero_odes || (!_iszero(rhs)))]
    else
        eqs = [Equation(0, rhs) for rhs in rhsvec if (include_zero_odes || (!_iszero(rhs)))]
    end
    eqs
end

"""
    foreach_noise_coeff(f, rx, species_to_idx, nps, depspec_submap; kwargs...)

Iterate over (species_idx, noise_coef) pairs for a reaction and call `f(i, coef)` for each.
The noise coefficient is `stoich * sqrt(|ratelaw|) * [noise_scaling]`.

This is a shared helper used by both `assemble_diffusion` (legacy noise matrix path) and
`add_noise_to_rhs!` (Brownian-based path) to avoid code duplication.
"""
function foreach_noise_coeff(f, rx, species_to_idx, nps, depspec_submap;
        combinatoric_ratelaws = true, remove_conserved = false,
        expand_catalyst_funs = true)
    rl = oderatelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws, expand_catalyst_funs)
    rlsqrt = sqrt(abs(rl))
    hasnoisescaling(rx) && (rlsqrt *= getnoisescaling(rx))
    remove_conserved && (rlsqrt = substitute(rlsqrt, depspec_submap))

    for (spec, stoich) in rx.netstoich
        remove_conserved && (spec in nps.depspecs) && continue
        drop_dynamics(spec) && continue
        if !haskey(species_to_idx, spec)
            error("Species $spec appears in reaction $rx but is not in the independent species list. " *
                  "This indicates a problem with the reaction system structure.")
        end
        i = species_to_idx[spec]
        coef = _signed_stoich_term(stoich, rlsqrt)
        f(i, coef)
    end
end

# this doesn't work with constraint equations currently
function assemble_diffusion(rs, sts, ispcs; combinatoric_ratelaws = true,
        remove_conserved = false, expand_catalyst_funs = true)
    # as BC species should ultimately get an equation, we include them in the noise matrix
    num_bcsts = count(isbc, get_unknowns(rs))

    # we make a matrix sized by the number of reactions
    eqs = Matrix{Any}(undef, length(sts) + num_bcsts, length(get_rxs(rs)))
    eqs .= 0
    species_to_idx = Dict(x => i for (i, x) in enumerate(ispcs))
    nps = get_networkproperties(rs)
    depspec_submap = if remove_conserved
        Dict(eq.lhs => eq.rhs for eq in nps.conservedeqs)
    else
        Dict()
    end

    for (j, rx) in enumerate(get_rxs(rs))
        foreach_noise_coeff(rx, species_to_idx, nps, depspec_submap;
                combinatoric_ratelaws, remove_conserved, expand_catalyst_funs) do i, coef
            eqs[i, j] = coef
        end
    end
    eqs
end

### Brownian Noise Helpers ###

"""
    create_sde_brownians(flatrs, scales)

Create one scalar Brownian variable per SDE-scale reaction. Returns a tuple
`(brownian_vars, brownian_map)` where `brownian_vars` is a `Vector{SymbolicT}` of
the Brownian variables and `brownian_map` is a `Vector{Pair{Int, SymbolicT}}` mapping
reaction indices to their Brownian variable.
"""
function create_sde_brownians(scales)
    sde_indices = [i for (i, s) in enumerate(scales) if s == PhysicalScale.SDE]
    brownian_vars = SymbolicT[]
    brownian_map = Pair{Int, SymbolicT}[]

    for (j, rx_idx) in enumerate(sde_indices)
        B = unwrap(only(@brownians $(Symbol(:B, :_, j))))
        push!(brownian_vars, B)
        push!(brownian_map, rx_idx => B)
    end
    brownian_vars, brownian_map
end

"""
    add_noise_to_rhs!(rhsvec, rs, ispcs, brownian_map; kwargs...)

Mutate the RHS vector `rhsvec` to add Brownian noise terms for each SDE-scale
reaction. For each such reaction, adds `stoich * sqrt(|ratelaw|) * [noise_scaling] * B_j`
to the corresponding species' RHS entry.
"""
function add_noise_to_rhs!(rhsvec, rs, ispcs, brownian_map;
        combinatoric_ratelaws = true, remove_conserved = false,
        expand_catalyst_funs = true)
    nps = get_networkproperties(rs)
    species_to_idx = Dict(x => i for (i, x) in enumerate(ispcs))
    depspec_submap = if remove_conserved
        Dict(eq.lhs => eq.rhs for eq in nps.conservedeqs)
    else
        Dict()
    end

    for (rx_idx, B_j) in brownian_map
        rx = get_rxs(rs)[rx_idx]
        foreach_noise_coeff(rx, species_to_idx, nps, depspec_submap;
                combinatoric_ratelaws, remove_conserved, expand_catalyst_funs) do i, coef
            rhsvec[i] += coef * B_j
        end
    end
end

### Jumps Assembly ###

"""
    jumpratelaw(rx; combinatoric_ratelaw=true)

Given a [`Reaction`](@ref), return the symbolic reaction rate law used in
generated stochastic chemical kinetics model SSAs for the reaction. Note,
for a reaction defined by

`k*X*Y, X+Z --> 2X + Y`

the expression that is returned will be `k*X^2*Y*Z`. For a reaction of
the form

`k, 2X+3Y --> Z`

the expression that is returned will be `k * binomial(X,2) *
binomial(Y,3)`.

Notes:
- Allocates
- `combinatoric_ratelaw=true` uses binomials in calculating the rate law, i.e. for `2S ->
  0` at rate `k` the ratelaw would be `k*S*(S-1)/2`. If `combinatoric_ratelaw=false` then
  the ratelaw is `k*S*(S-1)`, i.e. the rate law is not normalized by the scaling
  factor.
"""
function jumpratelaw(rx; combinatoric_ratelaw = true, expand_catalyst_funs = true)
    @unpack rate, substrates, substoich, only_use_rate = rx

    rl = rate
    expand_catalyst_funs && (rl = expand_registered_functions(rl))

    if !only_use_rate
        coef = eltype(substoich) <: Number ? one(eltype(substoich)) : 1
        for (i, stoich) in enumerate(substoich)
            s = substrates[i]
            if stoich isa SymbolicT
                rl *= combinatoric_ratelaw ? binomial(s, stoich) :
                      factorial(s) / factorial(s - stoich)
            else
                rl *= s
                isone(stoich) && continue
                for i in one(stoich):(stoich - one(stoich))
                    rl *= (s - i)
                end
                combinatoric_ratelaw && (coef *= factorial(stoich))
            end
        end
        combinatoric_ratelaw && !isequal(coef, one(coef)) && (rl /= coef)
    end
    rl
end

# if haveivdep=false then time dependent rates will still be classified as mass action
"""
```julia
ismassaction(rx, rs; rxvars = get_variables(rx.rate),
                              haveivdep = nothing,
                              unknownset = Set(unknowns(rs)),
                              ivset = nothing)
```

True if a given reaction is of mass action form, i.e. `rx.rate` does not depend
on any chemical species that correspond to unknowns of the system, and does not
depend explicitly on the independent variable (usually time).

# Arguments
- `rx`, the [`Reaction`](@ref).
- `rs`, a [`ReactionSystem`](@ref) containing the reaction.
- Optional: `rxvars`, `Variable`s which are not in `rxvars` are ignored as
  possible dependencies.
- Optional: `haveivdep`, `true` if the [`Reaction`](@ref) `rate` field
  explicitly depends on any independent variable (i.e. t or for spatial systems x,y,etc).
  If not set, will be automatically calculated.
- Optional: `unknownset`, set of unknowns which if the rxvars are within mean rx is
  non-mass action.
- Optional: `ivset`, a `Set` of the independent variables of the system. If not provided and
  the system is spatial, i.e. `isspatial(rs) == true`, it will be created with all the
  spatial variables and the time variable. If the rate expression contains any element of
  `ivset`, then `ismassaction(rx,rs) == false`. Pass a custom set to control this behavior.

Notes:
- Non-integer stoichiometry is treated as non-mass action. This includes
  symbolic variables/terms or floating point numbers for stoichiometric
  coefficients.
"""
function ismassaction(rx, rs; rxvars = get_variables(rx.rate),
        haveivdep::Union{Nothing, Bool} = nothing,
        unknownset = Set(get_unknowns(rs)), ivset = nothing)

    # we define non-integer (i.e. float or symbolic) stoich to be non-mass action
    ((eltype(rx.substoich) <: Integer) && (eltype(rx.prodstoich) <: Integer)) ||
        return false

    # if no dependencies must be zero order
    isempty(rxvars) && return true

    if (haveivdep === nothing)
        if isspatial(rs)
            if (ivset === nothing)
                ivs = Set(get_sivs(rs))
                push!(ivs, get_iv(rs))
                ivdep = any(in(ivs), rxvars)
            else
                ivdep = any(in(ivset), rxvars)
            end
        else
            ivdep = any(isequal(get_iv(rs)), rxvars)
        end
        ivdep && return false
    else
        haveivdep && return false
    end
    rx.only_use_rate && return false
    @inbounds for var in rxvars
        # not mass action if have a non-constant, variable species in the rate expression
        (var in unknownset) && return false
    end

    return true
end

@inline function makemajump(rx; combinatoric_ratelaw = true)
    @unpack rate, substrates, substoich, netstoich = rx
    zeroorder = (length(substoich) == 0)
    reactant_stoch = Vector{Pair{Any, eltype(substoich)}}()
    @inbounds for (i, spec) in enumerate(substrates)
        # move constant species into the rate
        if isconstant(spec)
            rate *= spec
            isone(substoich[i]) && continue
            for i in 1:(substoich[i] - 1)
                rate *= spec - i
            end
        else
            push!(reactant_stoch, substrates[i] => substoich[i])
        end
    end

    if (!zeroorder) && combinatoric_ratelaw
        coef = prod(factorial, substoich)
        (!isone(coef)) && (rate /= coef)
    end

    net_stoch = filter(p -> !drop_dynamics(p[1]), netstoich)
    isempty(net_stoch) &&
        error("$rx has no net stoichiometry change once accounting for constant and boundary condition species. This is not supported.")

    MassActionJump(Num(rate), reactant_stoch, net_stoch, scale_rates = false,
        useiszero = false)
end

# recursively visit each neighbor's rooted tree and mark everything in it as vrj
function dfs_mark!(isvrjvec, visited, depgraph, i)
    visited[i] = true
    nhbrs = depgraph[i]
    for nhbr in nhbrs
        if !visited[nhbr]
            isvrjvec[nhbr] = true
            dfs_mark!(isvrjvec, visited, depgraph, nhbr)
        end
    end
    nothing
end

# get_depgraph(rs)[i] is the list of reactions with rates depending on species changed by
# i'th reaction.
function get_depgraph(rs)
    eqs = reactions(rs)
    jdeps = asgraph(rs; eqs)
    vdeps = variable_dependencies(rs; eqs)
    eqeq_dependencies(jdeps, vdeps).fadjlist
end

# note that reactions that are not constant rate are treated as vrjs here
function classify_vrjs(rs, physcales)
    # first we determine vrjs with an explicit time-dependent rate
    rxs = get_rxs(rs)
    isvrjvec = falses(length(rxs))
    havevrjs = false
    rxvars = Set{SymbolicT}()
    for (i, rx) in enumerate(rxs)
        if physcales[i] in NON_CONSTANT_JUMP_SCALES
            isvrjvec[i] = true
            havevrjs = true
            continue
        end

        empty!(rxvars)
        (rx.rate isa SymbolicT) && get_variables!(rxvars, rx.rate)
        @inbounds for rxvar in rxvars
            if isequal(rxvar, get_iv(rs)) || (!MT.isparameter(rxvar) && !isspecies(rxvar))
                isvrjvec[i] = true
                havevrjs = true
                break
            end
        end
    end

    # now we determine vrj's that depend on species modified by a previous vrj
    if havevrjs
        depgraph = get_depgraph(rs)
        visited = falses(length(isvrjvec))
        for (i, isvrj) in enumerate(isvrjvec)
            if isvrj && !visited[i]
                # dfs from the vrj node to propagate vrj classification
                dfs_mark!(isvrjvec, visited, depgraph, i)
            end
        end
    end

    isvrjvec
end

function assemble_jumps(rs; combinatoric_ratelaws = true, physical_scales = nothing,
        expand_catalyst_funs = true, save_positions = (true, true))
    meqs = MassActionJump[]
    ceqs = ConstantRateJump[]
    veqs = VariableRateJump[]
    unknownset = Set(get_unknowns(rs))
    rxs = get_rxs(rs)

    if physical_scales === nothing
        physcales = [PhysicalScale.Jump for _ in enumerate(rxs)]
    else
        physcales = physical_scales
    end

    (isempty(get_rxs(rs)) || !any(in(JUMP_SCALES), physcales)) &&
        error("Must have at least one reaction that will be represented as a jump when constructing a JumpSystem.")

    # note isvrjvec indicates which reactions are not constant rate jumps
    # it may be that a given jump has isvrjvec[i] = true but has a physical
    isvrjvec = classify_vrjs(rs, physcales)

    rxvars = Set{SymbolicT}()
    for (i, rx) in enumerate(rxs)
        # only process reactions that should give jumps
        (physcales[i] in JUMP_SCALES) || continue

        empty!(rxvars)
        (rx.rate isa SymbolicT) && get_variables!(rxvars, rx.rate)

        isvrj = isvrjvec[i]
        if (!isvrj) && ismassaction(rx, rs; rxvars, haveivdep = false, unknownset)
            push!(meqs, makemajump(rx; combinatoric_ratelaw = combinatoric_ratelaws))
        else
            rl = jumpratelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws,
                expand_catalyst_funs)
            affect = Vector{Equation}()
            for (spec, stoich) in rx.netstoich
                # don't change species that are constant or BCs
                !drop_dynamics(spec) && push!(affect, spec ~ Pre(spec) + Pre(stoich))
            end
            if isvrj
                push!(veqs, VariableRateJump(rl, affect; save_positions))
            else
                push!(ceqs, ConstantRateJump(rl, affect))
            end
        end
    end
    reduce(vcat, (meqs, ceqs, veqs); init = JumpType[])
end

### Equation Coupling ###

# merge constraint components with the ReactionSystem components
# also handles removing BC and constant species
function addconstraints!(eqs, rs::ReactionSystem, ists, ispcs; remove_conserved = false,
        treat_conserved_as_eqs = false)
    # if there are BC species, put them after the independent species
    rssts = get_unknowns(rs)
    sts = any(isbc, rssts) ? vcat(ists, filter(isbc, rssts)) : ists
    ps = get_ps(rs)
    initeqs = Equation[]
    defs = MT.initial_conditions(rs)
    obs = MT.observed(rs)

    # make dependent species observables and add conservation constants as parameters
    if remove_conserved && !isempty(conservedequations(rs))
        nps = get_networkproperties(rs)

        # add the conservation constants as parameters and set their values
        ps = copy(ps)
        push!(ps, nps.conservedconst)

        if treat_conserved_as_eqs
            # add back previously removed dependent species
            sts = union(sts, nps.depspecs)

            # treat conserved eqs as normal eqs (lhs must be `0` in case structural simplify is not used)
            append!(eqs, [0 ~ eq.rhs - eq.lhs for eq in conservationlaw_constants(rs)])

            # add initialization equations for conserved parameters
            initialmap = Dict(u => Initial(u) for u in species(rs))
            conseqs = conservationlaw_constants(rs)
            initeqs = [Symbolics.substitute(conseq, initialmap) for conseq in conseqs]
        else
            # add the dependent species as observed
            obs = copy(obs)
            append!(obs, conservedequations(rs))
        end
    end

    ceqs = Equation[eq for eq in get_eqs(rs) if eq isa Equation]
    if !isempty(ceqs)
        if remove_conserved
            @info """
                  Be careful mixing ODEs or algebraic equations and elimination of
                  conservation laws. Catalyst does not check that the conserved equations
                  still hold for the final coupled system of equations. Consider using
                  `remove_conserved = false` and instead calling
                  ModelingToolkitBase.structural_simplify to simplify any generated ODESystem or
                  NonlinearSystem.
                  """
        end
        append!(eqs, ceqs)
    end

    eqs, sts, ps, obs, defs, initeqs
end

### Utility ###

# Throws an error when attempting to convert a spatial system to an unsupported type.
function spatial_convert_err(rs::ReactionSystem, systype)
    isspatial(rs) && error("Conversion to $systype is not supported for spatial networks.")
end

# Finds and differentials in an expression, and sets these to 0.
function remove_diffs(expr)
    if hasnode(Symbolics.is_derivative, expr)
        return replacenode(expr, diff_2_zero)
    else
        return expr
    end
end
diff_2_zero(expr) = (Symbolics.is_derivative(expr) ? 0 : expr)

COMPLETENESS_ERROR = "A ReactionSystem must be complete before it can be converted to other system types. A ReactionSystem can be marked as complete using the `complete` function."

### System Conversions ###

"""
    make_hybrid_model(rs::ReactionSystem; kwargs...)

Convert a [`ReactionSystem`](@ref) to a unified `ModelingToolkitBase.System` that can
contain ODE equations, Brownian noise terms, and/or jump processes depending on each
reaction's assigned [`PhysicalScale`](@ref).

# Keyword Arguments
- `name = nameof(rs)`: name for the returned `System`.
- `physical_scales = nothing`: overrides for per-reaction physical scales. Can be an
  iterable of `index => PhysicalScale` pairs, or a `Vector{PhysicalScale.T}` with one
  entry per reaction in the flattened system.
- `default_scale = PhysicalScale.Auto`: fallback scale for reactions with
  `PhysicalScale.Auto`. If any reaction remains `Auto` after resolution, an error is thrown.
- `combinatoric_ratelaws = get_combinatoric_ratelaws(rs)`: whether to use factorial/binomial
  scaling in rate laws.
- `include_zero_odes = true`: whether to include ODE equations with zero RHS.
- `remove_conserved = false`: whether to apply conservation law elimination. Not compatible
  with jump-scale reactions.
- `expand_catalyst_funs = true`: replace Catalyst functions (e.g. `hill`) with their
  rational form.
- `save_positions = (true, true)`: for `VariableRateJump`s, whether to save the solution
  before and/or after the jump.
- `checks = false`: whether to run `System` constructor checks.
- `initial_conditions = Dict()`: additional initial conditions to merge into the system defaults.

# Scale Resolution Order
1. `physical_scales` kwarg (user override per reaction index)
2. Per-reaction metadata (`get_physical_scale(rx)`)
3. `default_scale` kwarg (fallback for `Auto`)
4. If still `Auto` after all three → error
"""
function make_hybrid_model(rs::ReactionSystem;
        name = nameof(rs),
        physical_scales = nothing,
        default_scale = PhysicalScale.Auto,
        _override_all_scales = nothing,
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        include_zero_odes = true,
        remove_conserved = false,
        expand_catalyst_funs = true,
        save_positions = (true, true),
        checks = false,
        initial_conditions = Dict(),
        kwargs...)

    # Error checks.
    iscomplete(rs) || error(COMPLETENESS_ERROR)
    spatial_convert_err(rs, MT.System)

    flatrs = Catalyst.flatten(rs)

    # Resolve scales: _override_all_scales (internal) takes precedence over everything.
    if _override_all_scales !== nothing
        scales = fill(_override_all_scales, length(reactions(flatrs)))
    else
        scales = merge_physical_scales(reactions(flatrs), physical_scales, default_scale)
    end
    any(==(PhysicalScale.Auto), scales) &&
        error("Unresolved PhysicalScale.Auto scales remain. Provide `default_scale` or per-reaction `physical_scales`.")

    # Get user-provided brownians and jumps from the flattened ReactionSystem.
    user_brownians = MT.brownians(flatrs)
    user_jumps = MT.jumps(flatrs)

    # Layer 1: Reaction-based scale detection.
    # The || with _override_all_scales ensures the empty scales case is handled correctly.
    has_rxn_ode = any(==(PhysicalScale.ODE), scales) || (_override_all_scales == PhysicalScale.ODE)
    has_rxn_sde = any(==(PhysicalScale.SDE), scales) || (_override_all_scales == PhysicalScale.SDE)
    has_rxn_jump = any(in(JUMP_SCALES), scales)

    # Layer 2: User-provided elements (from ReactionSystem brownians/jumps fields).
    has_user_sde = !isempty(user_brownians)
    has_user_jump = !isempty(user_jumps)

    # Layer 3: Combined totals.
    has_ode = has_rxn_ode
    has_sde = has_rxn_sde || has_user_sde
    has_jump = has_rxn_jump || has_user_jump
    has_continuous = has_ode || has_sde || has_nonreactions(flatrs)

    # Conservation law elimination is not compatible with jump reactions.
    remove_conserved && has_jump &&
        throw(ArgumentError("Cannot remove conserved species with Jump-scale reactions."))
    remove_conserved && conservationlaws(flatrs)

    ists, ispcs = get_indep_sts(flatrs, remove_conserved)

    # --- Build drift RHS (ODE + SDE reactions both contribute drift) ---
    eqs = Equation[]
    brownian_vars = SymbolicT[]

    if has_continuous
        # Both ODE and SDE reactions contribute to drift; remap SDE→ODE for filtering
        # so that assemble_oderhs includes SDE reactions in the drift.
        drift_scales = copy(scales)
        for i in eachindex(drift_scales)
            drift_scales[i] == PhysicalScale.SDE && (drift_scales[i] = PhysicalScale.ODE)
        end

        rhsvec = assemble_oderhs(flatrs, ispcs; combinatoric_ratelaws, remove_conserved,
            physical_scales = drift_scales, expand_catalyst_funs)

        # Add Brownian noise terms for SDE-scale reactions.
        if has_rxn_sde
            rxn_brownian_vars, brownian_map = create_sde_brownians(scales)
            add_noise_to_rhs!(rhsvec, flatrs, ispcs, brownian_map;
                combinatoric_ratelaws, remove_conserved, expand_catalyst_funs)
            # Merge reaction-generated brownians with user-provided brownians.
            brownian_vars = unique(vcat(rxn_brownian_vars, user_brownians))
        else
            # Only user-provided brownians (no SDE-scale reactions).
            brownian_vars = user_brownians
        end

        # Convert RHS vector to D(x) ~ rhs equations.
        D = Differential(get_iv(flatrs))
        eqs = [Equation(D(x), rhs)
               for (x, rhs) in zip(ispcs, rhsvec)
               if (include_zero_odes || (!_iszero(rhs)))]
    end

    # --- Build jumps (Jump + VariableRateJump reactions only) ---
    rxn_jumps = JumpType[]
    if has_rxn_jump
        rxn_jumps = assemble_jumps(flatrs; combinatoric_ratelaws, expand_catalyst_funs,
            physical_scales = scales, save_positions)
    end
    # Merge reaction-generated jumps with user-provided jumps.
    jumps = vcat(rxn_jumps, user_jumps)

    # --- Add constraints (BC species, constraint equations, conserved species) ---
    if has_continuous
        eqs, us, ps, obs, defs = addconstraints!(eqs, flatrs, ists, ispcs; remove_conserved)
    else
        # Pure jump case.
        any(isbc, get_unknowns(flatrs)) &&
            (ists = vcat(ists, filter(isbc, get_unknowns(flatrs))))
        us = ists
        ps = get_ps(flatrs)
        obs = MT.observed(flatrs)
        defs = MT.initial_conditions(flatrs)
    end

    # --- Construct unified System ---
    # Note: brownians is a positional arg (5th) in the System constructor.
    MT.System(eqs, get_iv(flatrs), us, ps, brownian_vars;
        jumps,
        observed = obs,
        name,
        initial_conditions = merge(initial_conditions, defs),
        checks,
        continuous_events = MT.get_continuous_events(flatrs),
        discrete_events = MT.get_discrete_events(flatrs),
        metadata = MT.get_metadata(rs),
        kwargs...)
end

"""
```julia
Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
```
Convert a [`ReactionSystem`](@ref) to an `ModelingToolkitBase.ODESystem`.

Keyword args and default values:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate law,
  i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
  `combinatoric_ratelaws=false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
  ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
  itself defaults to true).
- `remove_conserved=false`, if set to `true` will calculate conservation laws of the
  underlying set of reactions (ignoring constraint equations), and then apply them to reduce
  the number of equations.
- `expand_catalyst_funs = true`, replaces Catalyst defined functions like `hill(A,B,C,D)`
  with their rational function representation when converting to another system type. Set to
  `false`` to disable.
"""
function make_rre_ode(rs::ReactionSystem; name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        include_zero_odes = true, remove_conserved = false, checks = false,
        initial_conditions = Dict(), expand_catalyst_funs = true,
        kwargs...)
    # Error if ReactionSystem has coupled brownians or jumps.
    flatrs = Catalyst.flatten(rs)
    if !isempty(MT.brownians(flatrs))
        error("""Cannot convert ReactionSystem with coupled brownian noise to a pure ODE system.
        Found brownians: $(MT.brownians(flatrs))
        Use `SDEProblem` or `HybridProblem` instead.""")
    end
    if !isempty(MT.jumps(flatrs))
        error("""Cannot convert ReactionSystem with coupled jumps to a pure ODE system.
        Found $(length(MT.jumps(flatrs))) jump(s).
        Use `JumpProblem` or `HybridProblem` instead.""")
    end

    make_hybrid_model(rs;
        _override_all_scales = PhysicalScale.ODE,
        name, combinatoric_ratelaws, include_zero_odes,
        remove_conserved, checks, initial_conditions,
        expand_catalyst_funs, kwargs...)
end

const NONLIN_PROB_REMAKE_WARNING = """
    Note, when constructing `NonlinearSystem`s with `remove_conserved = true`, possibly via \
    calling `NonlinearProblem`, `remake(::NonlinearProblem)` has some \
    limitations. If in `remake` the value of the conserved constant is \
    explicitly updated, it is not possible to have one variable's `u0` value \
    auto-update to preserve the conservation law. See the [FAQ \
    entry](https://docs.sciml.ai/Catalyst/stable/faqs/#faq_remake_nonlinprob) for more \
    details. This warning can be disabled by passing `conseqs_remake_warn = false`."""

function is_autonomous_error(iv)
    return """
    Attempting to convert a non-autonomous `ReactionSystem` (e.g. where some rate depends \
    on $(iv)) to a `NonlinearSystem`. This is not possible. if you are intending to \
    compute system steady states, consider creating and solving a `SteadyStateProblem."""
end

"""
```julia
Base.convert(::Type{<:NonlinearSystem},rs::ReactionSystem)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkitBase.NonlinearSystem`.

Keyword args and default values:
- `combinatoric_ratelaws = true` uses factorial scaling factors in calculating the rate law,
  i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
  `combinatoric_ratelaws=false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
  ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
  itself defaults to true).
- `remove_conserved = false`, if set to `true` will calculate conservation laws of the
  underlying set of reactions (ignoring coupled ODE or algebraic equations). For each
  conservation law one steady-state equation is eliminated, and replaced with the
  conservation law. This ensures a non-singular Jacobian.
- `conseqs_remake_warn = true`, set to false to disable warning about `remake` and
  conservation laws. See the [FAQ
  entry](https://docs.sciml.ai/Catalyst/stable/faqs/#faq_remake_nonlinprob) for more
  details.
- `expand_catalyst_funs = true`, replaces Catalyst defined functions like `hill(A,B,C,D)`
  with their rational function representation when converting to another system type. Set to
  `false`` to disable.
"""
function make_rre_algeqs(rs::ReactionSystem; name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        remove_conserved = false, conseqs_remake_warn = true, checks = false,
        initial_conditions = Dict(),
        all_differentials_permitted = false, expand_catalyst_funs = true, kwargs...)
    # Error checks.
    iscomplete(rs) || error(COMPLETENESS_ERROR)
    spatial_convert_err(rs::ReactionSystem, NonlinearSystem)
    remove_conserved && conseqs_remake_warn && (@warn NONLIN_PROB_REMAKE_WARNING)
    isautonomous(rs) || error(is_autonomous_error(get_iv(rs)))

    # Generates system equations.
    fullrs = Catalyst.flatten(rs)
    remove_conserved && conservationlaws(fullrs)
    ists, ispcs = get_indep_sts(fullrs, remove_conserved)
    eqs = assemble_drift(fullrs, ispcs; combinatoric_ratelaws, remove_conserved,
        as_odes = false, include_zero_odes = false, expand_catalyst_funs)
    eqs, us, ps, obs, defs, initeqs = addconstraints!(eqs, fullrs, ists, ispcs;
        remove_conserved, treat_conserved_as_eqs = true)

    # Throws a warning if there are differential equations in non-standard format.
    # Next, sets all differential terms to `0`.
    all_differentials_permitted || nonlinear_convert_differentials_check(rs)
    eqs = [remove_diffs(eq.lhs) ~ remove_diffs(eq.rhs) for eq in eqs]

    NonlinearSystem(eqs, us, ps;
        name,
        observed = obs, initialization_eqs = initeqs,
        initial_conditions = merge(initial_conditions, defs),
        checks,
        metadata = MT.get_metadata(rs),
        kwargs...)
end

# Ideally, when `ReactionSystem`s are converted to `NonlinearSystem`s, any coupled ODEs should be
# on the form D(X) ~ ..., where lhs is the time derivative w.r.t. a single variable, and the rhs
# does not contain any differentials. If this is not the case, we throw a warning to let the user
# know that they should be careful.
function nonlinear_convert_differentials_check(rs::ReactionSystem)
    for eq in filter(is_diff_equation, equations(rs))
        # For each differential equation, checks in order:
        # If there is a differential on the right hand side.
        # If the lhs is not on the form D(...).
        # If the lhs upper level function is not a differential w.r.t. time.
        # If the content of the differential is not a variable (and nothing more).
        # If either of this is a case, throws the warning.
        if hasnode(Symbolics.is_derivative, eq.rhs) ||
           !Symbolics.is_derivative(eq.lhs) ||
           !isequal(Symbolics.operation(eq.lhs), Differential(get_iv(rs))) ||
           (length(arguments(eq.lhs)) != 1) ||
           !any(isequal(arguments(eq.lhs)[1]), nonspecies(rs))
            error("You are attempting to convert a `ReactionSystem` coupled with differential equations to a `NonlinearSystem`. However, some of these differentials are not of the form `D(x) ~ ...` where:
                    (1) The left-hand side is a differential of a single variable with respect to the time independent variable, and
                    (2) The right-hand side does not contain any differentials.
                This is generally not permitted.

                If you still would like to perform this conversion, please use the `all_differentials_permitted = true` option. In this case, all differentials will be set to `0`.
                However, it is recommended to proceed with caution to ensure that the produced nonlinear equation makes sense for your intended application."
            )
        end
    end
end

"""
```julia
Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkitBase.SDESystem`.

Notes:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate law,
  i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
  `combinatoric_ratelaws=false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
  ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
  itself defaults to true).
- `remove_conserved=false`, if set to `true` will calculate conservation laws of the
  underlying set of reactions (ignoring constraint equations), and then apply them to reduce
  the number of equations.
- `expand_catalyst_funs = true`, replaces Catalyst defined functions like `hill(A,B,C,D)`
  with their rational function representation when converting to another system type. Set to
  `false`` to disable.
- `use_legacy_noise = true`, for simple SDE systems without constraints (no algebraic
  equations, no BC species), use the traditional `noise_eqs` matrix approach which avoids
  the need for `mtkcompile`. Set to `false` to use the Brownian-based approach.
"""
function make_cle_sde(rs::ReactionSystem;
        name = nameof(rs), combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        include_zero_odes = true, checks = false, remove_conserved = false,
        initial_conditions = Dict(), expand_catalyst_funs = true,
        use_legacy_noise = true,
        kwargs...)

    # Flatten once upfront and check for constraints.
    flatrs = Catalyst.flatten(rs)

    # Error if ReactionSystem has coupled jumps (SDE + jumps hybrid not supported yet).
    if !isempty(MT.jumps(flatrs))
        error("""Cannot convert ReactionSystem with coupled jumps to a pure SDE system.
        Found $(length(MT.jumps(flatrs))) jump(s).
        Use `HybridProblem` instead. Note: SDE+Jump hybrids require special handling.""")
    end

    has_constraints = has_alg_equations(flatrs) || any(isbc, get_unknowns(flatrs))
    has_user_brownians = !isempty(MT.brownians(flatrs))

    # For simple SDE systems without constraints and without user brownians,
    # use legacy noise_eqs matrix approach (avoids mtkcompile overhead).
    # If user brownians are present, use the new Brownian-based path.
    if use_legacy_noise && !has_constraints && !has_user_brownians
        iscomplete(rs) || error(COMPLETENESS_ERROR)
        spatial_convert_err(rs, MT.System)

        remove_conserved && conservationlaws(flatrs)
        ists, ispcs = get_indep_sts(flatrs, remove_conserved)

        eqs = assemble_drift(flatrs, ispcs; combinatoric_ratelaws, include_zero_odes,
            remove_conserved, expand_catalyst_funs)
        noiseeqs = assemble_diffusion(flatrs, ists, ispcs; combinatoric_ratelaws,
            remove_conserved, expand_catalyst_funs)
        eqs, us, ps, obs, defs = addconstraints!(eqs, flatrs, ists, ispcs; remove_conserved)

        if any(isbc, get_unknowns(flatrs))
            @info "Boundary condition species detected. As constraint equations are not currently supported when converting to SDESystems, the resulting system will be undetermined. Consider using constant species instead."
        end

        return MT.System(eqs, get_iv(flatrs), us, ps;
            noise_eqs = noiseeqs,
            observed = obs,
            name,
            initial_conditions = merge(initial_conditions, defs),
            checks,
            continuous_events = MT.get_continuous_events(flatrs),
            discrete_events = MT.get_discrete_events(flatrs),
            metadata = MT.get_metadata(rs),
            kwargs...)
    else
        # New path: Brownians via make_hybrid_model (requires mtkcompile for SDEProblem).
        return make_hybrid_model(flatrs;
            _override_all_scales = PhysicalScale.SDE,
            name, combinatoric_ratelaws, include_zero_odes,
            remove_conserved, checks, initial_conditions,
            expand_catalyst_funs, kwargs...)
    end
end

"""
    merge_physical_scales(rxs, physical_scales; default = PhysicalScale.Auto)

Merge physical scales for a set of reactions.

# Arguments
- `rxs`, a vector of `Reaction`s.
- `physical_scales`, an iterable of pairs mapping integer reaction indices to
  `PhysicalScale`s.
- `default`, the default physical scale to use for reactions that set PhysicalScale.Auto.
"""
function merge_physical_scales(rxs, physical_scales, default)
    scales = get_physical_scale.(rxs)

    # override metadata attached scales
    if physical_scales !== nothing
        for (key, scale) in physical_scales
            scales[key] = scale
        end
    end

    # transform any "Auto" scales to the default
    for (idx, scale) in enumerate(scales)
        if scale == PhysicalScale.Auto
            scales[idx] = default
        end
    end

    scales
end

# Overload for when physical_scales is already a fully-resolved vector of scales.
function merge_physical_scales(rxs, physical_scales::AbstractVector{<:PhysicalScale.T}, default)
    length(physical_scales) == length(rxs) ||
        error("Length of physical_scales ($(length(physical_scales))) must match number of reactions ($(length(rxs))).")
    scales = copy(physical_scales)
    for (idx, s) in enumerate(scales)
        s == PhysicalScale.Auto && (scales[idx] = default)
    end
    scales
end

# Returns (has_ode, has_sde, has_jump) for the resolved scales.
function detect_scale_types(scales)
    has_ode = any(==(PhysicalScale.ODE), scales)
    has_sde = any(==(PhysicalScale.SDE), scales)
    has_jump = any(s -> s in (PhysicalScale.Jump, PhysicalScale.VariableRateJump), scales)
    (has_ode, has_sde, has_jump)
end

"""
```julia
Base.convert(::Type{<:JumpSystem},rs::ReactionSystem; combinatoric_ratelaws=true)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkitBase.JumpSystem`.

Notes:
- `combinatoric_ratelaws=true` uses binomials in calculating the rate law, i.e. for `2S ->
  0` at rate `k` the ratelaw would be `k*S*(S-1)/2`. If `combinatoric_ratelaws=false` then
  the ratelaw is `k*S*(S-1)`, i.e. the rate law is not normalized by the scaling factor.
  Defaults to the value given when the `ReactionSystem` was constructed (which itself
  defaults to true).
- Does not currently support `ReactionSystem`s that include coupled algebraic or
  differential equations.
- Does not currently support continuous events as these are not supported by
  `ModelingToolkitBase.JumpSystems`.
- `expand_catalyst_funs = true`, replaces Catalyst defined functions like `hill(A,B,C,D)`
  with their rational function representation when converting to another system type. Set to
  `false`` to disable.
- `save_positions = (true, true)`, indicates whether for any reaction classified as a
  `VariableRateJump` to save the solution before and/or after the jump occurs. Defaults to
  true for both.
"""
function make_sck_jump(rs::ReactionSystem; name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        remove_conserved = nothing, checks = false, initial_conditions = Dict(),
        expand_catalyst_funs = true, save_positions = (true, true),
        physical_scales = nothing, kwargs...)
    (remove_conserved !== nothing) &&
        throw(ArgumentError("Catalyst does not support removing conserved species when converting to JumpSystems."))

    # Force all reactions to Jump, only preserving VariableRateJump metadata.
    # ODE/SDE metadata is ignored - use HybridProblem for hybrid systems.
    flatrs = Catalyst.flatten(rs)

    # Error on non-reaction ODE/algebraic/SDE equations (pure Jump only supports reactions).
    # This also catches brownians since they appear in SDE equations.
    non_rxn_eqs = filter(eq -> !(eq isa Reaction), equations(flatrs))
    if !isempty(non_rxn_eqs)
        error("""Cannot convert ReactionSystem with ODE, SDE, or algebraic equations to a pure Jump system.
        Found $(length(non_rxn_eqs)) non-reaction equation(s).
        Use `HybridProblem` instead for mixed ODE+Jump or SDE+Jump systems.""")
    end

    jump_scales = map(reactions(flatrs)) do rx
        get_physical_scale(rx) == PhysicalScale.VariableRateJump ?
            PhysicalScale.VariableRateJump : PhysicalScale.Jump
    end

    # Apply user overrides on top (if provided).
    if physical_scales !== nothing
        for (key, scale) in physical_scales
            jump_scales[key] = scale
        end
    end

    make_hybrid_model(flatrs;
        physical_scales = jump_scales,
        default_scale = PhysicalScale.Jump,
        name, combinatoric_ratelaws, checks, initial_conditions,
        expand_catalyst_funs, save_positions, kwargs...)
end

### Problems ###

# ODEProblem from AbstractReactionNetwork
function DiffEqBase.ODEProblem(rs::ReactionSystem, u0, tspan,
        p = DiffEqBase.NullParameters(), args...;
        check_length = false, name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        include_zero_odes = true, remove_conserved = false, checks = false,
        expand_catalyst_funs = true, structural_simplify = false, kwargs...)
    osys = make_rre_ode(rs; name, combinatoric_ratelaws, include_zero_odes, checks,
        remove_conserved, expand_catalyst_funs)

    # Handles potential differential algebraic equations (which requires `structural_simplify`).
    if structural_simplify
        osys = MT.mtkcompile(osys)
    elseif has_alg_equations(rs)
        error("The input ReactionSystem has algebraic equations. This requires setting `structural_simplify=true` within `ODEProblem` call.")
    else
        osys = complete(osys)
    end

    prob_cond = (p isa DiffEqBase.NullParameters) ? u0 : merge(Dict(u0), Dict(p))
    return ODEProblem(osys, prob_cond, tspan, args...; check_length, kwargs...)
end

"""
```julia
DiffEqBase.NonlinearProblem(rs::ReactionSystem, u0,
        p = DiffEqBase.NullParameters(), args...;
        name = nameof(rs), combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        remove_conserved = false, checks = false, check_length = false,
        structural_simplify = remove_conserved, all_differentials_permitted = false,
        kwargs...)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkitBase.NonlinearSystem`.

Keyword args and default values:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate law,
  i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
  `combinatoric_ratelaws=false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
  ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
  itself defaults to true).
- `remove_conserved=false`, if set to `true` will calculate conservation laws of the
  underlying set of reactions (ignoring coupled ODE or algebraic equations). For each
  conservation law one steady-state equation is eliminated, and replaced with the
  conservation law. This ensures a non-singular Jacobian. When using this option, it is
  recommended to call `ModelingToolkitBase.structural_simplify` on the converted system to then
  eliminate the conservation laws from the system equations.
- `conseqs_remake_warn = true`, set to false to disable warning about `remake` and
  conservation laws. See the [FAQ
  entry](https://docs.sciml.ai/Catalyst/stable/faqs/#faq_remake_nonlinprob) for more
  details.
- `expand_catalyst_funs = true`, replaces Catalyst defined functions like `hill(A,B,C,D)`
  with their rational function representation when converting to another system type. Set to
  `false`` to disable.
"""
function DiffEqBase.NonlinearProblem(rs::ReactionSystem, u0,
        p = DiffEqBase.NullParameters(), args...;
        name = nameof(rs), combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        remove_conserved = false, conseqs_remake_warn = true, checks = false,
        check_length = false, expand_catalyst_funs = true,
        structural_simplify = false, all_differentials_permitted = false, kwargs...)
    nlsys = make_rre_algeqs(rs; name, combinatoric_ratelaws, checks,
        all_differentials_permitted, remove_conserved, conseqs_remake_warn,
        expand_catalyst_funs)
    nlsys = structural_simplify ? MT.mtkcompile(nlsys) : complete(nlsys)
    prob_cond = (p isa DiffEqBase.NullParameters) ? u0 : merge(Dict(u0), Dict(p))
    return NonlinearProblem(nlsys, prob_cond, args...; check_length,
        kwargs...)
end

# SDEProblem from AbstractReactionNetwork
function DiffEqBase.SDEProblem(rs::ReactionSystem, u0, tspan,
        p = DiffEqBase.NullParameters(), args...;
        name = nameof(rs), combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        include_zero_odes = true, checks = false, check_length = false,
        remove_conserved = false, structural_simplify = false,
        expand_catalyst_funs = true, use_legacy_noise = true, kwargs...)

    # Flatten once upfront and pass to make_cle_sde.
    flatrs = Catalyst.flatten(rs)
    sde_sys = make_cle_sde(flatrs; name, combinatoric_ratelaws, expand_catalyst_funs,
        include_zero_odes, checks, remove_conserved, use_legacy_noise)

    # Determine if we need mtkcompile:
    # - If using Brownian-based approach (not legacy), mtkcompile extracts the noise matrix
    # - If there are algebraic equations, mtkcompile handles structural simplification
    # - If structural_simplify is requested explicitly
    has_constraints = has_alg_equations(flatrs) || any(isbc, get_unknowns(flatrs))
    needs_mtkcompile = structural_simplify ||
                       has_alg_equations(flatrs) ||
                       !use_legacy_noise ||
                       has_constraints

    prob_cond = (p isa DiffEqBase.NullParameters) ? u0 : merge(Dict(u0), Dict(p))

    if needs_mtkcompile
        if !structural_simplify && has_alg_equations(flatrs)
            error("The input ReactionSystem has algebraic equations. This requires setting `structural_simplify=true` within `SDEProblem` call.")
        end
        sde_sys = MT.mtkcompile(sde_sys)
        return SDEProblem(sde_sys, prob_cond, tspan, args...; check_length, kwargs...)
    else
        # Legacy path: complete + noise_rate_prototype
        sde_sys = complete(sde_sys)
        p_matrix = zeros(length(get_unknowns(sde_sys)), numreactions(flatrs))
        return SDEProblem(sde_sys, prob_cond, tspan, args...; check_length,
            noise_rate_prototype = p_matrix, kwargs...)
    end
end

# JumpProblem from ReactionSystem
# Note: For hybrid ODE+Jump systems, use HybridProblem instead.
function JumpProcesses.JumpProblem(rs::ReactionSystem, u0, tspan,
        p = DiffEqBase.NullParameters();
        name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        expand_catalyst_funs = true,
        save_positions = (true, true),
        checks = false,
        kwargs...)
    # Pure jump system - use HybridProblem for hybrid ODE+SDE+Jump systems.
    jsys = complete(make_sck_jump(rs; name, combinatoric_ratelaws, checks,
        expand_catalyst_funs, save_positions))
    op = (p isa DiffEqBase.NullParameters) ? u0 : merge(Dict(u0), Dict(p))
    return JumpProblem(jsys, op, tspan; save_positions, kwargs...)
end

"""
    HybridProblem(rs::ReactionSystem, u0, tspan, p = nothing;
                  physical_scales = nothing, default_scale = PhysicalScale.Jump, ...)

Create a problem from a [`ReactionSystem`](@ref) with per-reaction scale control.

This function uses `make_hybrid_model` internally and respects per-reaction
`PhysicalScale` metadata as well as `physical_scales` kwarg overrides.

The return type depends on which reaction scales are present:
- Pure ODE (only ODE-scale reactions) → `ODEProblem`
- Pure SDE or ODE+SDE (no jumps) → `SDEProblem`
- Any jumps present (ODE+Jump, SDE+Jump, ODE+SDE+Jump) → `JumpProblem`

For SDE+Jump combinations, the returned `JumpProblem` wraps an `SDEProblem` internally.

# Arguments
- `rs`: The ReactionSystem to convert.
- `u0`: Initial conditions as a mapping (e.g., `[:S => 100.0, :P => 0.0]`).
- `tspan`: Time span as a tuple (e.g., `(0.0, 10.0)`).
- `p`: Parameters as a mapping (e.g., `[:k1 => 1.0, :k2 => 0.5]`).

# Keyword Arguments
- `physical_scales = nothing`: Per-reaction scale overrides. Can be an iterable of
  `index => PhysicalScale` pairs.
- `default_scale = PhysicalScale.Jump`: Fallback for reactions with `PhysicalScale.Auto`.
  Defaults to `Jump` so that only reactions explicitly tagged as ODE/SDE are treated as continuous.
- `combinatoric_ratelaws = get_combinatoric_ratelaws(rs)`: Use factorial/binomial scaling.
- `save_positions = (true, true)`: For VariableRateJumps, save before/after jump.
- `structural_simplify = false`: Apply structural simplification (required for algebraic equations).
- Other kwargs passed to the underlying problem constructor.

# Returns
- `ODEProblem` if all reactions are ODE-scale
- `SDEProblem` if reactions are ODE/SDE-scale with no jumps
- `JumpProblem` if any jumps are present (wrapping `ODEProblem` for ODE+Jump, or `SDEProblem` for SDE+Jump)

# Example
```julia
# Hybrid ODE+Jump system
rn = @reaction_network begin
    k1, S --> P, [physical_scale = PhysicalScale.ODE]
    k2, P --> S, [physical_scale = PhysicalScale.Jump]
end
prob = HybridProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 10.0), [:k1 => 1.0, :k2 => 0.5])
sol = solve(prob, Tsit5())

# Pure ODE via HybridProblem
prob_ode = HybridProblem(rn, [:S => 100.0, :P => 0.0], (0.0, 10.0), [:k1 => 1.0, :k2 => 0.5];
    default_scale = PhysicalScale.ODE)
sol_ode = solve(prob_ode, Tsit5())

# SDE+Jump hybrid system (requires SDE solver like SRIW1 from StochasticDiffEq)
rn_sde_jump = @reaction_network begin
    k1, S --> P, [physical_scale = PhysicalScale.SDE]
    k2, P --> S, [physical_scale = PhysicalScale.Jump]
end
prob_sde_jump = HybridProblem(rn_sde_jump, [:S => 100.0, :P => 0.0], (0.0, 10.0), [:k1 => 1.0, :k2 => 0.5])
# prob_sde_jump.prob isa SDEProblem  # true - JumpProblem wraps SDEProblem
sol = solve(prob_sde_jump, SRIW1())
```
"""
function HybridProblem(rs::ReactionSystem, u0, tspan,
        p = DiffEqBase.NullParameters();
        name = nameof(rs),
        physical_scales = nothing,
        default_scale = PhysicalScale.Jump,
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        expand_catalyst_funs = true,
        save_positions = (true, true),
        checks = false,
        structural_simplify = false,
        kwargs...)

    # Determine which scale types are present.
    flatrs = Catalyst.flatten(rs)
    resolved_scales = merge_physical_scales(reactions(flatrs), physical_scales, default_scale)
    has_ode, has_sde, has_jump = detect_scale_types(resolved_scales)

    # Also check for user-provided brownians/jumps in the ReactionSystem.
    user_has_sde = !isempty(MT.brownians(flatrs))
    user_has_jump = !isempty(MT.jumps(flatrs))

    # Combine with reaction-detected scales
    has_sde = has_sde || user_has_sde
    has_jump = has_jump || user_has_jump

    # Build the unified System from the flattened ReactionSystem.
    sys = make_hybrid_model(flatrs; name, physical_scales, default_scale,
        combinatoric_ratelaws, expand_catalyst_funs, save_positions, checks)

    # Build problem conditions (u0 + p merged).
    prob_cond = (p isa DiffEqBase.NullParameters) ? u0 : merge(Dict(u0), Dict(p))

    if has_jump
        # Any jumps present → JumpProblem (wrapping ODEProblem or SDEProblem as needed)
        # For SDE+Jump: mtkcompile converts brownians → noise_eqs via extract_brownians_to_noise_eqs
        # For pure Jump: complete is sufficient (avoids unnecessary mtkcompile overhead)
        sys = has_sde ? MT.mtkcompile(sys) : complete(sys)
        return JumpProblem(sys, prob_cond, tspan; save_positions, kwargs...)

    elseif has_sde
        # SDE (with or without ODE) → SDEProblem
        # using Brownian variables for SDEs, so mtkcompile is always needed
        sys = MT.mtkcompile(sys)
        return SDEProblem(sys, prob_cond, tspan; kwargs...)

    else
        # Pure ODE → ODEProblem
        if structural_simplify
            sys = MT.mtkcompile(sys)
        elseif has_alg_equations(flatrs)
            error("The input ReactionSystem has algebraic equations. This requires setting `structural_simplify=true` within `HybridProblem` call.")
        else
            sys = complete(sys)
        end
        return ODEProblem(sys, prob_cond, tspan; kwargs...)
    end
end

# SteadyStateProblem from AbstractReactionNetwork
function DiffEqBase.SteadyStateProblem(rs::ReactionSystem, u0,
        p = DiffEqBase.NullParameters(), args...;
        check_length = false, name = nameof(rs),
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        remove_conserved = false, include_zero_odes = true, checks = false,
        expand_catalyst_funs = true, structural_simplify = false, kwargs...)
    osys = make_rre_ode(rs; name, combinatoric_ratelaws, include_zero_odes, checks,
        remove_conserved, expand_catalyst_funs)

    # Handles potential differential algebraic equations (which requires `structural_simplify`).
    if structural_simplify
        (osys = MT.mtkcompile(osys))
    elseif has_alg_equations(rs)
        error("The input ReactionSystem has algebraic equations. This requires setting `structural_simplify=true` within `ODEProblem` call.")
    else
        osys = complete(osys)
    end

    prob_cond = (p isa DiffEqBase.NullParameters) ? u0 : merge(Dict(u0), Dict(p))
    return SteadyStateProblem(osys, prob_cond, args...; check_length, kwargs...)
end

### Symbolic Variable/Symbol Conversions ###

# convert symbol of the form :sys.a.b.c to a symbolic a.b.c
function _symbol_to_var(sys, sym)
    if hasproperty(sys, sym)
        var = getproperty(sys, sym, namespace = false)
    else
        strs = split(String(sym), MT.NAMESPACE_SEPARATOR)   # need to check if this should be split of not!!!
        if length(strs) > 1
            var = getproperty(sys, Symbol(strs[1]), namespace = false)
            for str in view(strs, 2:length(strs))
                var = getproperty(var, Symbol(str), namespace = true)
            end
        else
            throw(ArgumentError("System $(nameof(sys)): variable $sym does not exist"))
        end
    end
    var
end

"""
    symmap_to_varmap(sys, symmap)

Given a system and map of `Symbol`s to values, generates a map from
corresponding symbolic variables/parameters to the values that can be used to
pass initial conditions and parameter mappings.

For example,
```julia
sir = @reaction_network sir begin
    β, S + I --> 2I
    ν, I --> R
end
subsys = @reaction_network subsys begin
    k, A --> B
end
@named sys = compose(sir, [subsys])
```
gives
```
Model sys with 3 equations
Unknowns (5):
  S(t)
  I(t)
  R(t)
  subsys₊A(t)
  subsys₊B(t)
Parameters (3):
  β
  ν
  subsys₊k
```
to specify initial condition and parameter mappings from *symbols* we can use
```julia
symmap = [:S => 1.0, :I => 1.0, :R => 1.0, :subsys₊A => 1.0, :subsys₊B => 1.0]
u0map  = symmap_to_varmap(sys, symmap)
pmap   = symmap_to_varmap(sys, [:β => 1.0, :ν => 1.0, :subsys₊k => 1.0])
```
`u0map` and `pmap` can then be used as input to various problem types.

Notes:
- Any `Symbol`, `sym`, within `symmap` must be a valid field of `sys`. i.e.
  `sys.sym` must be defined.
"""
function symmap_to_varmap(sys, symmap::Tuple)
    if all(p -> p isa Pair{Symbol}, symmap)
        return ((_symbol_to_var(sys, sym) => val for (sym, val) in symmap)...,)
    else  # if not all entries map a symbol to value pass through
        return symmap
    end
end

function symmap_to_varmap(sys, symmap::AbstractArray{Pair{Symbol, T}}) where {T}
    [_symbol_to_var(sys, sym) => val for (sym, val) in symmap]
end

function symmap_to_varmap(sys, symmap::Dict{Symbol, T}) where {T}
    Dict(_symbol_to_var(sys, sym) => val for (sym, val) in symmap)
end

# don't permute any other types and let varmap_to_vars handle erroring
symmap_to_varmap(sys, symmap) = symmap
#error("symmap_to_varmap requires a Dict, AbstractArray or Tuple to map Symbols to values.")

### Other Conversion-related Functions ###

"""
    to_multivariate_poly(polyeqs::AbstractVector{SymbolicT})

Convert the given system of polynomial equations to multivariate polynomial representation.
For example, this can be used in HomotopyContinuation.jl functions.
"""
function to_multivariate_poly(polyeqs::AbstractVector{Symbolics.SymbolicT})
    @assert length(polyeqs)>=1 "At least one expression must be passed to `multivariate_poly`."

    poly_to_bs = Dict{SymbolicUtils.PolyVarT, Symbolics.SymbolicT}()
    bs_to_poly = Dict{Symbolics.SymbolicT, SymbolicUtils.PolyVarT}()
    ps = map(polyeqs) do x
        if iscall(x) && operation(x) == (/)
            error("We should not be able to get here, please contact the package authors.")
        else
            SymbolicUtils.to_poly!(poly_to_bs, bs_to_poly, x, false)
        end
    end

    ps
end
