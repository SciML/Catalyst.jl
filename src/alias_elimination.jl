### Alias Elimination ###

# An alias `A ~ B` declares that `A` and `B` are the same entity: `A` is eliminated in favor of
# `B`. `flatten` (and so `complete`, which calls it) resolves a system's aliases: it validates
# them, substitutes the surviving symbol for each eliminated one throughout the system, keeps
# eliminated unknowns as observables (`A ~ B`) and eliminated parameters as bindings (`k1 => k2`),
# and stores the eliminated-to-surviving map in the resolved system's metadata.

### Alias Classes ###

"""
    AliasClass

Internal classification of a symbol for alias validation. An alias is only allowed between two
symbols of the same supported class.
"""
@enumx AliasClass begin
    Species          # ordinary species (non-BC, non-constant, non-compound)
    BCSpecies        # boundary condition species
    ConstantSpecies  # constant species (stored as parameters)
    CompoundSpecies  # compound species
    Unknown          # non-species unknowns (variables)
    Parameter        # ordinary parameters (excluding constant species)
    Brownian
    Poissonian
    Bound            # binding key
    Observable       # LHS of an observed equation
    Unsupported      # anything else (e.g. the independent variable, or a symbol not in the system)
end

# The classes that can be aliased (to a symbol of the same class).
const ALIASABLE_CLASSES = (
    AliasClass.Species, AliasClass.BCSpecies, AliasClass.ConstantSpecies,
    AliasClass.Unknown, AliasClass.Parameter,
)

# Classes whose aliases are resolved by the parameter map (the others by the unknown map).
const PARAMETER_ALIAS_CLASSES = (AliasClass.ConstantSpecies, AliasClass.Parameter)

### Classification ###

# Precomputed symbol sets of a flattened system, used to classify alias sides.
struct AliasClassificationContext
    species_set::Set{SymbolicT}
    unknown_set::Set{SymbolicT}
    param_set::Set{SymbolicT}
    brownian_set::Set{SymbolicT}
    poissonian_set::Set{SymbolicT}
    binding_keys::Set{SymbolicT}
    obs_lhs_set::Set{SymbolicT}
end

# Requires a flattened system (uses the non-recursive accessors).
function AliasClassificationContext(rs::ReactionSystem)
    return AliasClassificationContext(
        Set{SymbolicT}(get_species(rs)),
        Set{SymbolicT}(get_unknowns(rs)),
        Set{SymbolicT}(get_ps(rs)),
        Set{SymbolicT}(MT.get_brownians(rs)),
        Set{SymbolicT}(MT.get_poissonians(rs)),
        Set{SymbolicT}(keys(MT.get_bindings(rs))),
        Set{SymbolicT}(eq.lhs for eq in MT.get_observed(rs))
    )
end

# Classifies `x`; more specific classes are checked first.
function alias_class(x, ctx::AliasClassificationContext)
    x in ctx.binding_keys && return AliasClass.Bound
    x in ctx.obs_lhs_set && return AliasClass.Observable
    x in ctx.brownian_set && return AliasClass.Brownian
    x in ctx.poissonian_set && return AliasClass.Poissonian
    if x in ctx.species_set
        isbc(x) && return AliasClass.BCSpecies
        iscompound(x) && return AliasClass.CompoundSpecies
        return AliasClass.Species
    end
    x in ctx.unknown_set && return AliasClass.Unknown
    if x in ctx.param_set
        isconstant(x) && return AliasClass.ConstantSpecies
        return AliasClass.Parameter
    end
    return AliasClass.Unsupported
end

alias_compatible(lhs_class, rhs_class) = (lhs_class == rhs_class) && (lhs_class in ALIASABLE_CLASSES)

function alias_compatibility_error(lhs, rhs, lhs_class, rhs_class)
    for (x, class) in ((lhs, lhs_class), (rhs, rhs_class))
        class == AliasClass.Unsupported &&
            return ArgumentError("Cannot alias $lhs ~ $rhs: $x is not a species, variable or parameter of the system.")
    end
    if lhs_class == rhs_class
        return ArgumentError("Cannot alias $lhs ~ $rhs: aliasing $(lhs_class) symbols is not supported.")
    end
    return ArgumentError("Cannot alias $lhs ($(lhs_class)) to $rhs ($(rhs_class)). Only symbols of the same kind can be aliased: species with species, boundary condition species with boundary condition species, constant species with constant species, variables with variables, and parameters with parameters.")
end

### Validation ###

# Checks that the two sides of an alias have matching units and, for parameters, matching
# time dependence and type.
function check_alias_metadata(lhs, rhs, class)
    lhs_unit = catalyst_get_unit(lhs)
    rhs_unit = catalyst_get_unit(rhs)
    _units_match(lhs_unit, rhs_unit) ||
        throw(ArgumentError("Unit mismatch in alias $lhs ~ $rhs: $lhs has unit $lhs_unit, but $rhs has unit $rhs_unit."))
    class in PARAMETER_ALIAS_CLASSES || return nothing

    # `iscalledparameter` is not public in ModelingToolkitBase.
    if MT.iscalledparameter(lhs) != MT.iscalledparameter(rhs)
        throw(ArgumentError("Cannot alias $lhs ~ $rhs: one is a time-dependent (discrete) parameter and the other is not."))
    end
    lhs_type = MT.symtype(unwrap(lhs))
    rhs_type = MT.symtype(unwrap(rhs))
    isequal(lhs_type, rhs_type) ||
        throw(ArgumentError("Type mismatch in alias $lhs ~ $rhs: $lhs has type $lhs_type, but $rhs has type $rhs_type."))
    return nothing
end

# Cheap structural checks, run by the `ReactionSystem` constructor (full validation happens
# when aliases are resolved, as it needs the flattened system).
function check_aliases(aliases::Vector{Equation})
    seen_lhs = Set{Any}()
    for eq in aliases
        isequal(eq.lhs, eq.rhs) && throw(ArgumentError("Self-alias not allowed: $(eq.lhs) ~ $(eq.rhs)."))
        (eq.lhs in seen_lhs) &&
            throw(ArgumentError("Duplicate alias: $(eq.lhs) appears as the left-hand side of several aliases."))
        push!(seen_lhs, eq.lhs)
    end
    return nothing
end

### Resolution ###

@enum VisitState UNSEEN VISITING DONE

# The resolved aliases of a flattened system.
struct AliasResolution
    # Eliminated symbol => surviving symbol (chains resolved), for unknowns and parameters.
    unknown_map::Dict{SymbolicT, SymbolicT}
    param_map::Dict{SymbolicT, SymbolicT}
    # Eliminated symbol => its alias equation (for error messages).
    alias_of::Dict{SymbolicT, Equation}
end

# Follows each alias chain to its surviving symbol, with cycle detection.
function resolve_chains(graph::Dict{SymbolicT, SymbolicT})
    canonical = Dict{SymbolicT, SymbolicT}()
    visit_state = Dict{SymbolicT, VisitState}()
    path = SymbolicT[]
    for start in keys(graph)
        haskey(canonical, start) && continue
        empty!(path)
        current = start
        canon = start
        while true
            state = get(visit_state, current, UNSEEN)
            if state == DONE
                canon = get(canonical, current, current)
                break
            elseif state == VISITING
                cycle = path[findfirst(isequal(current), path):end]
                throw(ArgumentError("Alias cycle detected involving: $(join(cycle, ", "))."))
            end
            visit_state[current] = VISITING
            push!(path, current)
            if haskey(graph, current)
                current = graph[current]
            else
                canon = current
                break
            end
        end
        for node in path
            visit_state[node] = DONE
            isequal(node, canon) || (canonical[node] = canon)
        end
    end
    return canonical
end

# Validates the aliases of a flattened system and resolves their chains.
function resolve_aliases(alias_eqs::Vector{Equation}, rs::ReactionSystem)
    check_aliases(alias_eqs)
    ctx = AliasClassificationContext(rs)
    unknown_graph = Dict{SymbolicT, SymbolicT}()
    param_graph = Dict{SymbolicT, SymbolicT}()
    alias_of = Dict{SymbolicT, Equation}()
    for eq in alias_eqs
        lhs, rhs = unwrap(eq.lhs), unwrap(eq.rhs)
        lhs_class = alias_class(lhs, ctx)
        rhs_class = alias_class(rhs, ctx)
        alias_compatible(lhs_class, rhs_class) ||
            throw(alias_compatibility_error(lhs, rhs, lhs_class, rhs_class))
        check_alias_metadata(lhs, rhs, lhs_class)
        graph = (lhs_class in PARAMETER_ALIAS_CLASSES) ? param_graph : unknown_graph
        graph[lhs] = rhs
        alias_of[lhs] = eq
    end
    return AliasResolution(resolve_chains(unknown_graph), resolve_chains(param_graph), alias_of)
end

# The alias equations that lead from `x` to its surviving symbol (for error messages).
function alias_chain(x, res::AliasResolution)
    chain = Equation[]
    while haskey(res.alias_of, x)
        eq = res.alias_of[x]
        push!(chain, eq)
        x = unwrap(eq.rhs)
    end
    return chain
end

### Substitution Helpers ###

# Substitutes through an expression, including inside differentials and other operators.
subs_alias(x, submap) = Symbolics.substitute(x, submap; filterer = Returns(true))
subs_alias(x::Number, submap) = x
subs_alias(x::Num, submap) = Num(subs_alias(unwrap(x), submap))
subs_alias(::Nothing, submap) = nothing
subs_alias(eq::Equation, submap) = subs_alias(eq.lhs, submap) ~ subs_alias(eq.rhs, submap)

# Whether a (possibly symbolic) stoichiometric coefficient is provably zero.
function is_zero_stoich(x)
    v = Symbolics.value(x)
    return v isa Number && iszero(v)
end

# Substitutes and merges duplicate species in a substrate or product list (keeping first-seen
# order). Stoichiometric coefficients are substituted too (parametric stoichiometry).
function merge_species_list(species, stoich, submap)
    merged_species = SymbolicT[]
    merged_stoich = Any[]
    seen_index = Dict{SymbolicT, Int}()
    for (sp, st) in zip(species, stoich)
        canonical = get(submap, sp, sp)
        st = subs_alias(st, submap)
        idx = get(seen_index, canonical, 0)
        if idx == 0
            push!(merged_species, canonical)
            push!(merged_stoich, st)
            seen_index[canonical] = length(merged_species)
        else
            merged_stoich[idx] += st
        end
    end
    # Numeric stoichiometries keep a concrete type; symbolic ones are passed as `Num`s, which the
    # `Reaction` constructor recognizes as symbolic stoichiometry.
    isempty(merged_stoich) && return merged_species, eltype(stoich)[]
    all(st -> st isa Number, merged_stoich) && return merged_species, map(identity, merged_stoich)
    return merged_species, Num[st for st in merged_stoich]
end

# Whether a reaction's net change is provably zero for every species.
function has_zero_net_change(subs, substoich, prods, prodstoich)
    net = Dict{SymbolicT, Any}()
    for (sp, st) in zip(subs, substoich)
        net[sp] = get(net, sp, 0) - st
    end
    for (sp, st) in zip(prods, prodstoich)
        net[sp] = get(net, sp, 0) + st
    end
    return all(is_zero_stoich, values(net))
end

function substitute_reaction_metadata(metadata, submap)
    return Pair{Symbol, Any}[key => ((val isa Num || val isa SymbolicT) ? subs_alias(val, submap) : val)
        for (key, val) in metadata]
end

# Substitutes the aliases through reaction `rx` (the `idx`th reaction). Throws an error if this
# turns a reaction with a net change into one without (a self-loop).
function resolve_reaction(rx::Reaction, idx, submap, res::AliasResolution)
    subs, substoich = merge_species_list(rx.substrates, rx.substoich, submap)
    prods, prodstoich = merge_species_list(rx.products, rx.prodstoich, submap)
    if !has_zero_net_change(rx.substrates, rx.substoich, rx.products, rx.prodstoich) &&
            has_zero_net_change(subs, substoich, prods, prodstoich)
        involved = filter(sp -> haskey(submap, sp), [rx.substrates; rx.products])
        chains = unique(reduce(vcat, (alias_chain(sp, res) for sp in involved); init = Equation[]))
        new_rx = Reaction(rx.rate, subs, prods, substoich, prodstoich; only_use_rate = rx.only_use_rate)
        throw(ArgumentError("Alias $(join(("`$eq`" for eq in chains), ", ")) collapses reaction $idx, `$rx`, to `$new_rx`, which has no net change. Remove the redundant reaction or keep the species distinct."))
    end
    rate = subs_alias(rx.rate, submap)
    metadata = substitute_reaction_metadata(rx.metadata, submap)
    return Reaction(rate, subs, prods, substoich, prodstoich; only_use_rate = rx.only_use_rate, metadata)
end

# Substitutes through an event affect; only symbolic affects are supported.
subs_affect(::Nothing, submap) = nothing
function subs_affect(aff::MT.SymbolicAffect, submap)
    eqs = Equation[subs_alias(eq, submap) for eq in aff.affect]
    discretes = SymbolicT[get(submap, p, p) for p in aff.discrete_parameters]
    return MT.SymbolicAffect(eqs, discretes)
end
function subs_affect(aff, submap)
    throw(ArgumentError("Aliases can't yet be resolved in systems with events whose affects are not symbolic equations (got an affect of type $(typeof(aff)))."))
end

# Substitutes through an event's symbolic fields, keeping all its other settings.
function resolve_event(evt::MT.SymbolicContinuousCallback, submap)
    @set! evt.conditions = Equation[subs_alias(c, submap) for c in evt.conditions]
    @set! evt.affect = subs_affect(evt.affect, submap)
    @set! evt.affect_neg = subs_affect(evt.affect_neg, submap)
    @set! evt.initialize = subs_affect(evt.initialize, submap)
    @set! evt.finalize = subs_affect(evt.finalize, submap)
    return evt
end
function resolve_event(evt::MT.SymbolicDiscreteCallback, submap)
    (evt.conditions isa SymbolicT) && (@set! evt.conditions = subs_alias(evt.conditions, submap))
    @set! evt.affect = subs_affect(evt.affect, submap)
    @set! evt.initialize = subs_affect(evt.initialize, submap)
    @set! evt.finalize = subs_affect(evt.finalize, submap)
    return evt
end

# Substitutes through a (symbolic) jump, keeping its other fields.
function resolve_jump(j::Union{ConstantRateJump, VariableRateJump}, submap)
    @set! j.rate = subs_alias(j.rate, submap)
    @set! j.affect! = Equation[subs_alias(eq, submap) for eq in j.affect!]
    return j
end
function resolve_jump(j::MassActionJump, submap)
    throw(ArgumentError("Aliases can't yet be resolved in systems with a `MassActionJump` among their jumps. Use a `ConstantRateJump` or `VariableRateJump` instead."))
end

# Substitutes through a poissonian's rate (stored as metadata on the poissonian;
# `getpoissonianrate` and `topoissonian` are not public in ModelingToolkitBase).
function resolve_poissonian(dN, submap)
    rate = MT.getpoissonianrate(dN)
    rate === nothing && return dN
    new_rate = subs_alias(rate, submap)
    return isequal(new_rate, rate) ? dN : MT.topoissonian(dN, new_rate)
end

# Substitutes through the default values stored in the metadata of the surviving symbols (the
# `ReactionSystem` constructor reads these again).
function resolve_symbol_defaults(syms, submap)
    return map(syms) do sym
        MT.hasdefault(sym) || return sym
        old = MT.getdefault(sym)
        new = subs_alias(old, submap)
        return isequal(old, new) ? sym : MT.setdefault(sym, new)
    end
end

# Removes the default value stored in a symbol's metadata.
clear_default(sym) = MT.hasdefault(sym) ? setmetadata(sym, Symbolics.VariableDefaultValue, nothing) : sym

# Moves initial conditions from eliminated symbols to surviving ones. Values are resolved before
# comparing, so `X => k1` and `Y => k2` agree for `X ~ Y` and `k1 ~ k2`.
function resolve_initial_conditions(ics, submap)
    new_ics = Dict{SymbolicT, Any}()
    for (key, val) in ics
        new_key = get(submap, key, key)
        new_val = subs_alias(val, submap)
        if haskey(new_ics, new_key) && !isequal(new_ics[new_key], new_val)
            throw(ArgumentError("Conflicting initial conditions for aliased symbols: $new_key would get both $(new_ics[new_key]) and $new_val."))
        end
        new_ics[new_key] = new_val
    end
    return new_ics
end

### Alias Map Storage ###

# Metadata key for the eliminated-to-surviving maps of a system whose aliases were resolved.
struct AliasSubstitutionMaps end

# The maps `(; unknown_map, param_map)` stored on `rs` itself, or `nothing`.
get_alias_maps(rs) = MT.getmetadata(rs, AliasSubstitutionMaps, nothing)

# The maps stored on `rs` and (namespaced) on its subsystems, which may have had their aliases
# resolved before being composed.
function collect_alias_maps(rs::ReactionSystem)
    unknown_map = Dict{SymbolicT, SymbolicT}()
    param_map = Dict{SymbolicT, SymbolicT}()
    maps = get_alias_maps(rs)
    if maps !== nothing
        merge!(unknown_map, maps.unknown_map)
        merge!(param_map, maps.param_map)
    end
    for sys in get_systems(rs)
        sys_maps = collect_alias_maps(sys)
        ns(x) = unwrap(MT.renamespace(nameof(sys), x))
        for (k, v) in sys_maps.unknown_map
            unknown_map[ns(k)] = ns(v)
        end
        for (k, v) in sys_maps.param_map
            param_map[ns(k)] = ns(v)
        end
    end
    return (; unknown_map, param_map)
end

# Whether the aliases of the system or of one of its subsystems were already resolved
# (eliminating at least one symbol).
function has_resolved_aliases(rs::ReactionSystem)
    maps = collect_alias_maps(rs)
    return !(isempty(maps.unknown_map) && isempty(maps.param_map))
end

# The metadata for `extend(sys, rs)`: the merged maps of the two systems, if either has any (their
# subsystems' maps are collected from the subsystems). `MetadataT` is not public in
# ModelingToolkitBase.
function extended_alias_metadata(rs::ReactionSystem, sys::ReactionSystem)
    maps = filter(!isnothing, [get_alias_maps(rs), get_alias_maps(sys)])
    isempty(maps) && return MT.MetadataT()
    unknown_map = merge(Dict{SymbolicT, SymbolicT}(), (m.unknown_map for m in maps)...)
    param_map = merge(Dict{SymbolicT, SymbolicT}(), (m.param_map for m in maps)...)
    return Base.ImmutableDict(MT.MetadataT(), AliasSubstitutionMaps => (; unknown_map, param_map))
end

# The metadata of `rs` for its flattened system, which keeps the maps of `rs` and its subsystems.
function flattened_alias_metadata(rs::ReactionSystem)
    metadata = MT.get_metadata(rs)
    has_resolved_aliases(rs) || return metadata
    return Base.ImmutableDict(metadata, AliasSubstitutionMaps => collect_alias_maps(rs))
end

### Main Resolution Function ###

# Resolves the aliases of a flattened system (see the top of this file). Called by `flatten`.
function eliminate_aliases(rs::ReactionSystem; name = nameof(rs))
    isempty(get_systems(rs)) || error("Internal error: aliases can only be resolved on a flattened system.")
    alias_eqs = get_aliases(rs)
    isempty(alias_eqs) && return rs

    res = resolve_aliases(alias_eqs, rs)
    submap = merge(res.unknown_map, res.param_map)

    # Reactions (no reaction is dropped, so reaction indices stay valid).
    rxs = Reaction[resolve_reaction(rx, idx, submap, res) for (idx, rx) in enumerate(get_rxs(rs))]

    # Non-reaction equations; at most one differential equation per variable.
    noneqs = Equation[subs_alias(eq, submap) for eq in get_eqs(rs) if !(eq isa Reaction)]
    defined = Dict{Any, Equation}()
    for (eq, new_eq) in zip((eq for eq in get_eqs(rs) if !(eq isa Reaction)), noneqs)
        lhs = unwrap(new_eq.lhs)
        Symbolics.is_derivative(lhs) || continue
        var = only(Symbolics.arguments(lhs))
        if haskey(defined, var)
            throw(ArgumentError("The aliases give $var two defining equations: `$(defined[var])` and `$eq`. Remove one of them or keep the variables distinct."))
        end
        defined[var] = eq
    end

    # Observables: existing ones are rewritten, eliminated unknowns are added.
    observed = Equation[eq.lhs ~ subs_alias(eq.rhs, submap) for eq in MT.get_observed(rs)]
    for (eliminated, surviving) in res.unknown_map
        push!(observed, eliminated ~ surviving)
    end

    # Bindings: existing ones are rewritten, eliminated parameters are bound to surviving ones.
    bindings = Dict{SymbolicT, Any}(k => subs_alias(v, submap) for (k, v) in MT.get_bindings(rs))
    for (eliminated, surviving) in res.param_map
        bindings[eliminated] = surviving
    end

    # Events, jumps, poissonians, tstops and initial conditions.
    cevents = [resolve_event(evt, submap) for evt in MT.get_continuous_events(rs)]
    devents = [resolve_event(evt, submap) for evt in MT.get_discrete_events(rs)]
    jumps = JumpType[resolve_jump(j, submap) for j in MT.get_jumps(rs)]
    poissonians = SymbolicT[resolve_poissonian(dN, submap) for dN in MT.get_poissonians(rs)]
    tstops = Any[subs_alias(ts, submap) for ts in MT.get_tstops(rs)]
    ics = resolve_initial_conditions(MT.get_initial_conditions(rs), submap)

    # Unknowns lose the eliminated ones; eliminated parameters stay (as bound parameters). The
    # default values of eliminated parameters, moved to the parameters replacing them with the
    # initial conditions above, are removed from their metadata (which the `ReactionSystem`
    # constructor reads again), as bound parameters can't have initial conditions.
    unknowns = resolve_symbol_defaults(
        filter(u -> !haskey(res.unknown_map, u), get_unknowns(rs)), submap)
    ps = resolve_symbol_defaults(
        [haskey(res.param_map, p) ? clear_default(p) : p for p in get_ps(rs)], submap)

    # The maps also keep the symbols eliminated earlier (e.g. in a subsystem whose aliases were
    # resolved before it was composed), with the symbols that replaced them resolved in turn.
    # `MetadataT` is not public in ModelingToolkitBase.
    old_maps = something(get_alias_maps(rs), (; unknown_map = Dict{SymbolicT, SymbolicT}(),
        param_map = Dict{SymbolicT, SymbolicT}()))
    unknown_map = merge(Dict{SymbolicT, SymbolicT}(k => get(submap, v, v) for (k, v) in old_maps.unknown_map),
        res.unknown_map)
    param_map = merge(Dict{SymbolicT, SymbolicT}(k => get(submap, v, v) for (k, v) in old_maps.param_map),
        res.param_map)
    metadata = Base.ImmutableDict(MT.get_metadata(rs), AliasSubstitutionMaps => (; unknown_map, param_map))

    new_rs = ReactionSystem(
        CatalystEqType[rxs; noneqs], get_iv(rs), unknowns, ps, MT.get_brownians(rs);
        poissonians, jumps, observed, name, bindings, initial_conditions = ics,
        checks = false, combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        balanced_bc_check = false, spatial_ivs = get_sivs(rs),
        continuous_events = cevents, discrete_events = devents, tstops, metadata,
        aliases = Equation[])
    return MT.iscomplete(rs) ? complete(new_rs; flatten = false) : new_rs
end

### Querying Aliases ###

# The system with its aliases resolved (itself if they already are).
alias_resolved_system(rs::ReactionSystem) = aliases_present(rs) ? flatten(rs) : rs

"""
    alias_map(rs::ReactionSystem)

Return a `Dict` mapping each symbol eliminated by an alias to the symbol it was replaced by (with
alias chains followed to the end). Empty if `rs` has no aliases.

`remake` does not translate eliminated symbols, so use the surviving symbols with it; see also
[`uneliminated_unknowns`](@ref) and [`uneliminated_parameters`](@ref).
"""
function alias_map(rs::ReactionSystem)
    maps = collect_alias_maps(alias_resolved_system(rs))
    return merge(maps.unknown_map, maps.param_map)
end

"""
    uneliminated_unknowns(rs::ReactionSystem)

Return the unknowns of `rs` that remain after its aliases are eliminated: the ones to give initial
conditions for, and to use with `remake`.
"""
uneliminated_unknowns(rs::ReactionSystem) = unknowns(alias_resolved_system(rs))

"""
    uneliminated_parameters(rs::ReactionSystem)

Return the parameters of `rs` that remain after its aliases are eliminated (eliminated parameters
stay readable, bound to the parameter that replaced them, but can't be given values): the ones to
give values for, and to use with `remake`.
"""
function uneliminated_parameters(rs::ReactionSystem)
    resolved = alias_resolved_system(rs)
    param_map = collect_alias_maps(resolved).param_map
    return filter(p -> !haskey(param_map, unwrap(p)), parameters(resolved))
end

### Problem Inputs ###

# Remaps eliminated symbols in a problem's `u0` or `p` input to the symbols that replaced them,
# keeping the input's form (a `Tuple` stays a `Tuple`, so values keep their types). Errors if
# two entries give an eliminated symbol and its replacement different values.
function remap_alias_inputs(input, rs::ReactionSystem)
    (input isa Union{AbstractDict, Tuple, AbstractVector}) || return input
    isempty(input) && return input
    all(x -> x isa Pair, input) || return input
    submap = alias_map(rs)
    isempty(submap) && return input

    # Each key as a symbolic variable (a `Symbol` that names none in `rs` is kept as is), so that
    # a `Symbol` and a symbolic key for the same symbol are recognized as such.
    pairs = collect(input)
    vars = map(((key, _),) -> (key isa Symbol) ? _try_symbol_to_var(rs, key) : unwrap(key), pairs)
    any(var -> haskey(submap, var), vars) || return input

    # The symbol each entry gives a value for. Entries for the same symbol must agree, and the
    # first of them is kept.
    targets = [get(submap, var, var) for var in vars]
    first_entry = Dict{Any, Int}()
    for (i, target) in enumerate(targets)
        j = get!(first_entry, target, i)
        (j == i || isequal(last(pairs[j]), last(pairs[i]))) && continue
        throw(ArgumentError("Conflicting values for aliased symbols: $(first(pairs[j])) is given $(last(pairs[j])) and $(first(pairs[i])) is given $(last(pairs[i])), but both refer to $target."))
    end
    remapped = [targets[i] => last(pairs[i]) for i in sort!(collect(values(first_entry)))]
    input isa Tuple && return Tuple(remapped)
    input isa AbstractDict && return Dict(remapped)
    return remapped
end

# `Symbol` key → its symbolic variable in `rs`, or the `Symbol` itself if `rs` has none.
function _try_symbol_to_var(rs, sym::Symbol)
    var = try
        _symbol_to_var(rs, sym)
    catch
        return sym
    end
    return unwrap(var)
end
