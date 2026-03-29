### Alias Elimination ###

# This file contains the alias classification, validation, resolution, and elimination
# logic for ReactionSystems. Aliases declare that two symbols represent the same entity
# (e.g. species A is the same as species B). Alias elimination substitutes eliminated
# symbols throughout the system, producing a simplified system.

### AliasClass EnumX ###

"""
    @enumx AliasClass

Classification of a symbolic variable for alias validation purposes.
Used to enforce same-category aliasing and reject unsupported combinations.
"""
@enumx AliasClass begin
    Species          # ordinary species (non-BC, non-constant, non-compound)
    BCSpecies        # boundary condition species (in species/unknowns, has BC metadata)
    ConstantSpecies  # constant species (in parameters, has ParameterConstantSpecies metadata)
    CompoundSpecies  # compound species (components reference other species)
    Unknown          # non-species unknowns (variables)
    Parameter        # ordinary parameters (excludes constant species)
    Brownian
    Poissonian
    Bound            # binding key
    Observable       # LHS of an observed equation
    Unsupported      # anything else (e.g. independent variable)
end

### Classification Context ###

# Precomputed sets for efficient classification. Built once per validation/elimination call.
struct AliasClassificationContext
    species_set::Set{SymbolicT}
    unknown_set::Set{SymbolicT}
    param_set::Set{SymbolicT}
    brownian_set::Set{SymbolicT}
    poissonian_set::Set{SymbolicT}
    binding_keys::Set{SymbolicT}
    obs_lhs_set::Set{SymbolicT}
end

function _make_classification_context(rs::ReactionSystem)
    AliasClassificationContext(
        Set(get_species(rs)),
        Set(get_unknowns(rs)),
        Set(get_ps(rs)),
        Set(MT.get_brownians(rs)),
        Set(MT.get_poissonians(rs)),
        Set(keys(MT.get_bindings(rs))),
        Set(eq.lhs for eq in MT.get_observed(rs)))
end

### Classification Function (flattened systems only) ###

"""
    alias_class(x, ctx::AliasClassificationContext) -> AliasClass.T

Classify a symbolic variable for alias validation. Requires a flattened system
(uses non-recursive `get_*` accessors via the precomputed context).
Order matters — more specific checks come first.
"""
function alias_class(x, ctx::AliasClassificationContext)
    x in ctx.binding_keys && return AliasClass.Bound
    x in ctx.obs_lhs_set && return AliasClass.Observable
    x in ctx.brownian_set && return AliasClass.Brownian
    x in ctx.poissonian_set && return AliasClass.Poissonian
    if x in ctx.species_set
        isbc(x) && return AliasClass.BCSpecies
        hasmetadata(x, CompoundSpecies) && return AliasClass.CompoundSpecies
        return AliasClass.Species
    end
    x in ctx.unknown_set && return AliasClass.Unknown
    if x in ctx.param_set
        isconstant(x) && return AliasClass.ConstantSpecies
        return AliasClass.Parameter
    end
    return AliasClass.Unsupported
end

### Compatibility Predicates ###

# V1 allowed pairings — must be identical class on both sides.
const _ALIAS_COMPATIBLE_CLASSES = Set([
    AliasClass.Species,
    AliasClass.BCSpecies,
    AliasClass.ConstantSpecies,
    AliasClass.Unknown,
    AliasClass.Parameter,
])

"""
    alias_compatible(lhs_class::AliasClass.T, rhs_class::AliasClass.T) -> Bool

Check whether two alias classes are compatible for aliasing in v1.
Both sides must be the same supported class.
"""
function alias_compatible(lhs_class::AliasClass.T, rhs_class::AliasClass.T)
    return lhs_class == rhs_class && lhs_class in _ALIAS_COMPATIBLE_CLASSES
end

"""
    alias_compatibility_error(lhs, rhs, lhs_class, rhs_class) -> ErrorException

Construct an informative error message for incompatible alias classes.
"""
function alias_compatibility_error(lhs, rhs, lhs_class::AliasClass.T, rhs_class::AliasClass.T)
    # Same class but unsupported
    if lhs_class == rhs_class
        return ErrorException(
            "Cannot alias $(lhs_class) variables: $lhs ~ $rhs. " *
            "$(lhs_class) aliasing is not yet supported in Catalyst.")
    end
    # Different classes
    return ErrorException(
        "Cannot alias $lhs ($(lhs_class)) to $rhs ($(rhs_class)). " *
        "Only same-type aliasing is allowed: Species↔Species, BCSpecies↔BCSpecies, " *
        "ConstantSpecies↔ConstantSpecies, Unknown↔Unknown, Parameter↔Parameter.")
end

### Metadata Compatibility ###

"""
    alias_metadata_compatible(lhs, rhs, alias_type::AliasClass.T)

Check metadata compatibility for an alias pair. Throws informative error on mismatch.
V1 uses strict same-category matching — no metadata is ever propagated or mutated.
"""
function alias_metadata_compatible(lhs, rhs, alias_type::AliasClass.T)
    # Unit check — applies to all supported types
    _check_alias_units(lhs, rhs)

    if alias_type == AliasClass.Parameter || alias_type == AliasClass.ConstantSpecies
        _check_alias_parameter_metadata(lhs, rhs)
    end
    nothing
end

# Check unit metadata compatibility.
function _check_alias_units(lhs, rhs)
    lhs_unit = MT.getmetadata(lhs, MT.VariableUnit, nothing)
    rhs_unit = MT.getmetadata(rhs, MT.VariableUnit, nothing)
    if lhs_unit !== nothing && rhs_unit !== nothing
        lhs_unit == rhs_unit || error(
            "Unit mismatch in alias $lhs ~ $rhs: " *
            "$lhs has unit $lhs_unit but $rhs has unit $rhs_unit.")
    elseif lhs_unit !== nothing || rhs_unit !== nothing
        error("Unit mismatch in alias $lhs ~ $rhs: " *
              "one side has units and the other does not. " *
              "Both sides must have matching units or neither should have units.")
    end
end

# Check parameter-specific metadata: time-dependent, discrete, declared type.
function _check_alias_parameter_metadata(lhs, rhs)
    # Time-dependent (called parameter): both must match
    lhs_called = MT.iscalledparameter(lhs)
    rhs_called = MT.iscalledparameter(rhs)
    if lhs_called != rhs_called
        error("Time-dependency mismatch in alias $lhs ~ $rhs: " *
              "$(lhs_called ? "$lhs is time-dependent" : "$lhs is not time-dependent") but " *
              "$(rhs_called ? "$rhs is time-dependent" : "$rhs is not time-dependent"). " *
              "Both must match.")
    end

    # Declared type: if both have type metadata, must match
    lhs_type = MT.symtype(unwrap(lhs))
    rhs_type = MT.symtype(unwrap(rhs))
    if lhs_type != rhs_type
        error("Type mismatch in alias $lhs ~ $rhs: " *
              "$lhs has type $lhs_type but $rhs has type $rhs_type. Both must match.")
    end
end

### Constructor-Time Checks (Cheap, Structural Only) ###

"""
    check_aliases(aliases::Vector{Equation})

Cheap structural checks for aliases, called during constructor when `checks=true`.
Does NOT require a flattened system or symbol lookups. Full validation
(type classification, metadata, cycles) is deferred to `eliminate_aliases`.
"""
function check_aliases(aliases::Vector{Equation})
    seen_lhs = Set{Any}()
    for eq in aliases
        lhs, rhs = eq.lhs, eq.rhs
        # Self-alias check
        isequal(lhs, rhs) && error("Self-alias not allowed: $lhs ~ $rhs")
        # Duplicate LHS check
        lhs in seen_lhs && error("Duplicate alias: $lhs appears as LHS in multiple aliases")
        push!(seen_lhs, lhs)
    end
end

### Alias Resolution ###

# Visitor states for cycle detection.
@enum VisitState UNSEEN VISITING DONE

# Result of alias resolution — substitution maps for unknowns and parameters.
struct AliasResolutionResult
    unknown_submap::Dict{SymbolicT, SymbolicT}
    param_submap::Dict{SymbolicT, SymbolicT}
end

"""
    build_alias_graphs(aliases, rs) -> (unknown_graph, param_graph)

Build directed alias graphs from alias equations. Performs full validation:
type classification, compatibility, and metadata checks. Requires flattened system.
"""
function build_alias_graphs(aliases::Vector{Equation}, rs::ReactionSystem)
    ctx = _make_classification_context(rs)

    unknown_graph = Dict{SymbolicT, SymbolicT}()
    param_graph = Dict{SymbolicT, SymbolicT}()
    seen_lhs = Set{SymbolicT}()

    for eq in aliases
        lhs, rhs = eq.lhs, eq.rhs

        # Self-alias check
        isequal(lhs, rhs) && error("Self-alias not allowed: $lhs ~ $rhs")

        # Duplicate LHS check
        lhs in seen_lhs && error("Duplicate alias: $lhs appears as LHS in multiple aliases")
        push!(seen_lhs, lhs)

        # Classify both sides
        lhs_class = alias_class(lhs, ctx)
        rhs_class = alias_class(rhs, ctx)

        # Compatibility check
        alias_compatible(lhs_class, rhs_class) ||
            throw(alias_compatibility_error(lhs, rhs, lhs_class, rhs_class))

        # Metadata compatibility
        alias_metadata_compatible(lhs, rhs, lhs_class)

        # Insert into appropriate graph
        if lhs_class in (AliasClass.Parameter, AliasClass.ConstantSpecies)
            param_graph[lhs] = rhs
        else
            unknown_graph[lhs] = rhs
        end
    end

    return unknown_graph, param_graph
end

"""
    resolve_canonical!(graph, canonical, visit_state)

Resolve alias chains to canonical representatives with cycle detection.
Uses iterative DFS with path compression. O(n) complexity.
"""
function resolve_canonical!(
    graph::Dict{SymbolicT, SymbolicT},
    canonical::Dict{SymbolicT, SymbolicT},
    visit_state::Dict{SymbolicT, VisitState}
)
    path = SymbolicT[]  # reused across iterations
    for start in keys(graph)
        haskey(canonical, start) && continue

        empty!(path)
        current = start
        local canon

        while true
            state = get(visit_state, current, UNSEEN)

            if state == DONE
                canon = get(canonical, current, current)
                break
            elseif state == VISITING
                cycle_start = findfirst(isequal(current), path)
                cycle_syms = path[cycle_start:end]
                error("Alias cycle detected involving: $cycle_syms")
            else  # UNSEEN
                visit_state[current] = VISITING
                push!(path, current)

                if haskey(graph, current)
                    current = graph[current]
                else
                    canon = current
                    break
                end
            end
        end

        # Path compression — only store eliminated nodes (not the canonical root itself)
        for node in path
            visit_state[node] = DONE
            isequal(node, canon) || (canonical[node] = canon)
        end
    end
end

"""
    validate_and_resolve_aliases(aliases, rs) -> AliasResolutionResult

Full validation and resolution of alias equations on a flattened system.
Returns substitution maps for unknowns and parameters, or `nothing` if empty.
"""
function validate_and_resolve_aliases(aliases::Vector{Equation}, rs::ReactionSystem)
    isempty(aliases) && return nothing

    unknown_graph, param_graph = build_alias_graphs(aliases, rs)

    unknown_canonical = Dict{SymbolicT, SymbolicT}()
    param_canonical = Dict{SymbolicT, SymbolicT}()
    visit_state = Dict{SymbolicT, VisitState}()

    resolve_canonical!(unknown_graph, unknown_canonical, visit_state)
    resolve_canonical!(param_graph, param_canonical, visit_state)

    return AliasResolutionResult(unknown_canonical, param_canonical)
end

### Substitution Helpers ###

"""
    substitute_reaction(rx, unknown_submap, combined, has_unknown_aliases,
        has_param_aliases, eliminated_unknown_set)

Substitute aliases through a reaction. Returns `nothing` if the reaction becomes a no-op.
Uses cheap structural guards to skip work when possible.
"""
function substitute_reaction(rx::Reaction, unknown_submap, combined,
        has_unknown_aliases, has_param_aliases, eliminated_unknown_set)
    # Rate substitution is unconditional (expression-level guards deferred; see plan).
    new_rate = Symbolics.substitute(rx.rate, combined)

    # Skip species list merging if no unknown aliases affect this reaction's species
    if has_unknown_aliases &&
            ((!isnothing(rx.substrates) && any(in(eliminated_unknown_set), rx.substrates)) ||
             (!isnothing(rx.products) && any(in(eliminated_unknown_set), rx.products)))
        new_subs, new_substoich = merge_species_list(rx.substrates, rx.substoich, unknown_submap)
        new_prods, new_prodstoich = merge_species_list(rx.products, rx.prodstoich, unknown_submap)
    else
        new_subs, new_substoich = rx.substrates, rx.substoich
        new_prods, new_prodstoich = rx.products, rx.prodstoich
    end

    # No-op detection
    if is_noop_reaction(new_subs, new_substoich, new_prods, new_prodstoich)
        return nothing
    end

    # Substitute through symbolic reaction metadata (e.g. noise_scaling)
    new_metadata = isempty(combined) ? rx.metadata :
        substitute_reaction_metadata(rx.metadata, combined)

    Reaction(new_rate, new_subs, new_prods, new_substoich, new_prodstoich;
        only_use_rate = rx.only_use_rate, metadata = new_metadata)
end

"""
    merge_species_list(species_list, stoich_list, submap)

Substitute and merge duplicate species in a substrate/product list.
Preserves first-seen order deterministically.
"""
function merge_species_list(species_list, stoich_list, submap)
    isnothing(species_list) && return (nothing, nothing)

    merged_species = SymbolicT[]
    merged_stoich = eltype(stoich_list)[]
    seen_index = Dict{SymbolicT, Int}()

    for (sp, st) in zip(species_list, stoich_list)
        canonical = get(submap, sp, sp)
        idx = get(seen_index, canonical, 0)
        if idx == 0
            push!(merged_species, canonical)
            push!(merged_stoich, st)
            seen_index[canonical] = length(merged_species)
        else
            merged_stoich[idx] += st
        end
    end

    return (merged_species, merged_stoich)
end

"""
    is_noop_reaction(subs, substoich, prods, prodstoich)

Return `true` if a reaction has zero net stoichiometry after merging.
"""
function is_noop_reaction(subs, substoich, prods, prodstoich)
    (isnothing(subs) || isempty(subs)) && (isnothing(prods) || isempty(prods)) && return true

    net = Dict{SymbolicT, Any}()
    if !isnothing(subs)
        for (sp, st) in zip(subs, substoich)
            net[sp] = get(net, sp, 0) - st
        end
    end
    if !isnothing(prods)
        for (sp, st) in zip(prods, prodstoich)
            net[sp] = get(net, sp, 0) + st
        end
    end

    return all(iszero, values(net))
end

"""
    substitute_reaction_metadata(metadata, combined_submap)

Substitute aliases through symbolic reaction metadata values (e.g. noise_scaling).
Non-symbolic values are passed through unchanged.
"""
function substitute_reaction_metadata(metadata, combined_submap)
    return [key => (Symbolics.issym(val) || Symbolics.iscall(val) ?
        Symbolics.substitute(val, combined_submap) : val)
        for (key, val) in metadata]
end

"""
    substitute_equation(eq, combined_submap)

Substitute aliases through both sides of an equation.
"""
function substitute_equation(eq::Equation, combined_submap)
    new_lhs = Symbolics.substitute(eq.lhs, combined_submap)
    new_rhs = Symbolics.substitute(eq.rhs, combined_submap)
    return new_lhs ~ new_rhs
end

# Substitute aliases through a SymbolicContinuousCallback.
function substitute_continuous_event(evt::MT.SymbolicContinuousCallback, combined_submap)
    new_conds = [Symbolics.substitute(c, combined_submap) for c in evt.conditions]
    new_affect = _substitute_affect(evt.affect, combined_submap)
    new_affect_neg = _substitute_affect(evt.affect_neg, combined_submap)
    MT.SymbolicContinuousCallback(new_conds, new_affect; affect_neg = new_affect_neg)
end

# Substitute aliases through a SymbolicDiscreteCallback.
function substitute_discrete_event(evt::MT.SymbolicDiscreteCallback, combined_submap)
    new_conds = Symbolics.substitute(evt.conditions, combined_submap)
    new_affect = _substitute_affect(evt.affect, combined_submap)
    MT.SymbolicDiscreteCallback(new_conds, new_affect)
end

# Substitute through a SymbolicAffect (or nothing).
function _substitute_affect(aff, combined_submap)
    isnothing(aff) && return nothing
    new_eqs = [substitute_equation(eq, combined_submap) for eq in aff.affect]
    # Remap discrete_parameters through alias submap (aliased discrete params stay stale otherwise).
    new_disc_ps = [get(combined_submap, dp, dp) for dp in aff.discrete_parameters]
    MT.SymbolicAffect(new_eqs; discrete_parameters = new_disc_ps)
end

# Substitute through a jump. Per-type handling:
# - ConstantRateJump: substitute rate and affect!
# - VariableRateJump: use @set! to preserve all extra fields (save_positions, bounds, etc.)
# - MassActionJump: error guard (stoichiometry remapping deferred to future version)
function substitute_jump(j::ConstantRateJump, combined_submap)
    new_rate = Symbolics.substitute(j.rate, combined_submap)
    new_affect = Symbolics.substitute(j.affect!, combined_submap)
    ConstantRateJump(new_rate, new_affect)
end

function substitute_jump(j::VariableRateJump, combined_submap)
    j_new = j
    @set! j_new.rate = Symbolics.substitute(j.rate, combined_submap)
    @set! j_new.affect! = Symbolics.substitute(j.affect!, combined_submap)
    return j_new
end

function substitute_jump(j::MassActionJump, combined_submap)
    error("Alias elimination through MassActionJump is not yet supported. " *
          "Please eliminate aliases before adding MassActionJumps, or use " *
          "ConstantRateJump/VariableRateJump instead.")
end

"""
    substitute_variable_defaults(syms, combined_submap)

Substitute aliases through VariableDefaultValue metadata on surviving symbols.
The MTKBase System constructor re-extracts VariableDefaultValue via process_variables!,
so stale references must be fixed. Only creates new symbol objects when the default
actually references an eliminated symbol (isequal short-circuit).
"""
function substitute_variable_defaults(syms::Vector{SymbolicT}, combined_submap)
    return map(syms) do sym
        MT.hasdefault(sym) || return sym
        old_def = MT.getdefault(sym)
        new_def = Symbolics.substitute(old_def, combined_submap)
        isequal(old_def, new_def) && return sym
        return MT.setdefault(sym, new_def)
    end
end

### Initial Conditions Transfer ###

"""
    transfer_initial_conditions(ics, unknown_submap, param_submap, combined)

Transfer initial conditions from eliminated symbols to canonical ones.
Errors on conflicting values. Also substitutes through values.
"""
function transfer_initial_conditions(ics, unknown_submap, param_submap, combined)
    new_ics = copy(ics)

    # Phase 1: Remap keys (eliminated → canonical) with conflict detection
    for submap in (unknown_submap, param_submap)
        for (eliminated, canonical) in submap
            elim_val = get(new_ics, eliminated, nothing)
            canon_val = get(new_ics, canonical, nothing)

            if elim_val !== nothing && canon_val !== nothing
                if !isequal(elim_val, canon_val)
                    error("Initial condition conflict for alias $eliminated ~ $canonical: " *
                          "$eliminated has value $elim_val but $canonical has value $canon_val. " *
                          "These must agree or only one should be specified.")
                end
            elseif elim_val !== nothing
                new_ics[canonical] = elim_val
            end

            delete!(new_ics, eliminated)
        end
    end

    # Phase 2: Substitute through values (they can be symbolic expressions
    # referencing aliased parameters, e.g. [B => 2*k1] with k1 aliased)
    for (key, val) in new_ics
        new_ics[key] = Symbolics.substitute(val, combined)
    end

    return new_ics
end

### Alias Map Metadata Storage ###

"""
    AliasSubstitutionMaps

Metadata key type for storing alias substitution maps on a system.
Following the pattern from `reactionsystem_metadata.jl` (U0Map, ParameterMap).
"""
struct AliasSubstitutionMaps end

function store_alias_maps_in_metadata(metadata, unknown_submap, param_submap)
    # MetadataT is ImmutableDict{DataType, Any} — extend by constructing a new one.
    return Base.ImmutableDict(metadata,
        AliasSubstitutionMaps => (;
            unknown_submap = Dict(unknown_submap),
            param_submap = Dict(param_submap)))
end

"""
    get_alias_maps(sys)

Retrieve alias substitution maps from system metadata, or `nothing` if not present.
"""
function get_alias_maps(sys)
    MT.getmetadata(sys, AliasSubstitutionMaps, nothing)
end

### Main Elimination Function ###

"""
    eliminate_aliases(rs::ReactionSystem; name = nameof(rs))

Return a new `ReactionSystem` with all alias relationships resolved by substitution.
Eliminated unknowns become observables. Eliminated parameters are removed.
Requires a flattened system (no subsystems).

Supported alias types (v1): ordinary species, BC-species, constant species,
non-species unknowns (variables), ordinary parameters.

Not supported: Brownians, Poissonians, compounds, bindings, observables.
"""
function eliminate_aliases(rs::ReactionSystem; name = nameof(rs))
    # Preconditions
    !isempty(get_systems(rs)) &&
        error("System must be flattened before alias elimination. Call `flatten` first.")
    alias_eqs = get_aliases(rs)
    isempty(alias_eqs) && return rs

    # Validate and resolve
    result = validate_and_resolve_aliases(alias_eqs, rs)
    unknown_submap = result.unknown_submap
    param_submap = result.param_submap
    combined = merge(unknown_submap, param_submap)

    # Pre-build sets for cheap membership checks (avoids per-expression allocation).
    has_unknown_aliases = !isempty(unknown_submap)
    has_param_aliases = !isempty(param_submap)
    eliminated_unknown_set = has_unknown_aliases ? Set(keys(unknown_submap)) : Set{SymbolicT}()

    # Substitute through reactions
    new_rxs = Reaction[]
    for rx in get_rxs(rs)
        new_rx = substitute_reaction(rx, unknown_submap, combined,
            has_unknown_aliases, has_param_aliases, eliminated_unknown_set)
        new_rx === nothing && continue
        push!(new_rxs, new_rx)
    end

    # Substitute through non-reaction equations (skip if no aliases touch expressions)
    new_noneq_eqs = Equation[]
    for eq in get_eqs(rs)
        eq isa Reaction && continue
        push!(new_noneq_eqs, substitute_equation(eq, combined))
    end
    new_eqs = CatalystEqType[new_rxs; new_noneq_eqs]

    # Rewrite existing observables (RHS only)
    new_obs = Equation[]
    for obs_eq in MT.get_observed(rs)
        new_rhs = Symbolics.substitute(obs_eq.rhs, combined)
        push!(new_obs, obs_eq.lhs ~ new_rhs)
    end

    # Add eliminated unknowns as new observables (pointing to canonical)
    for (eliminated, canonical) in unknown_submap
        push!(new_obs, eliminated ~ canonical)
    end

    # Substitute through events (skip if empty)
    cevents = MT.get_continuous_events(rs)
    new_cevents = isempty(cevents) ? cevents :
        [substitute_continuous_event(evt, combined) for evt in cevents]
    devents = MT.get_discrete_events(rs)
    new_devents = isempty(devents) ? devents :
        [substitute_discrete_event(evt, combined) for evt in devents]

    # Substitute through jumps (skip if empty)
    old_jumps = MT.get_jumps(rs)
    new_jumps = isempty(old_jumps) ? old_jumps :
        [substitute_jump(j, combined) for j in old_jumps]

    # Substitute through tstops (skip if empty or no param aliases)
    old_tstops = MT.get_tstops(rs)
    new_tstops = (isempty(old_tstops) || !has_param_aliases) ? collect(old_tstops) :
        [Symbolics.substitute(ts, combined) for ts in old_tstops]

    # Handle initial_conditions
    new_ics = transfer_initial_conditions(
        MT.initial_conditions(rs), unknown_submap, param_submap, combined)

    # Rebuild bindings with substituted values (skip if empty)
    old_bindings = MT.get_bindings(rs)
    if isempty(old_bindings)
        new_bindings_dict = SymmapT()
    else
        new_bindings_dict = SymmapT()
        for (key, val) in old_bindings
            new_bindings_dict[key] = Symbolics.substitute(val, combined)
        end
    end

    # Remove eliminated symbols
    new_unknowns = has_unknown_aliases ?
        filter(u -> u ∉ eliminated_unknown_set, get_unknowns(rs)) :
        collect(get_unknowns(rs))
    eliminated_params = has_param_aliases ? Set(keys(param_submap)) : Set{SymbolicT}()
    new_params = has_param_aliases ?
        filter(p -> p ∉ eliminated_params, get_ps(rs)) :
        collect(get_ps(rs))

    # Substitute through VariableDefaultValue on surviving symbols
    new_unknowns = substitute_variable_defaults(new_unknowns, combined)
    new_params = substitute_variable_defaults(new_params, combined)

    # Store alias maps in metadata for downstream input remapping
    new_metadata = store_alias_maps_in_metadata(
        MT.get_metadata(rs), unknown_submap, param_submap)

    # Construct new system with cleared aliases
    new_rs = ReactionSystem(new_eqs, get_iv(rs), new_unknowns, new_params, MT.get_brownians(rs);
        poissonians = MT.get_poissonians(rs),
        jumps = new_jumps,
        observed = new_obs,
        name,
        bindings = new_bindings_dict,
        initial_conditions = new_ics,
        checks = false,
        combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
        balanced_bc_check = false,
        spatial_ivs = get_sivs(rs),
        continuous_events = new_cevents,
        discrete_events = new_devents,
        tstops = new_tstops,
        metadata = new_metadata,
        aliases = Equation[])

    # Preserve completeness status from the input system.
    if MT.iscomplete(rs)
        new_rs = complete(new_rs)
    end
    return new_rs
end

### Input Remapping ###

"""
    remap_alias_inputs(u0_or_p, sys)

Remap eliminated symbol keys in user-provided `u0` or `p` maps to their canonical
equivalents. Works on any system with alias maps in metadata.
Handles both `Dict` and `Vector{Pair}` inputs (preserving format and ordering).
Errors if both eliminated and canonical keys are provided with different values.
"""
function remap_alias_inputs(u0_or_p, sys)
    maps = get_alias_maps(sys)
    isnothing(maps) && return u0_or_p

    combined = merge(maps.unknown_submap, maps.param_submap)
    isempty(combined) && return u0_or_p

    return _remap_inputs(u0_or_p, combined)
end

# No-op for NullParameters (common default for p in problem constructors).
remap_alias_inputs(p::DiffEqBase.NullParameters, sys) = p

# Dict path: uses get/setindex!/delete! directly.
function _remap_inputs(d::AbstractDict, combined)
    result = copy(d)
    for (key, val) in d
        canonical = get(combined, key, nothing)
        isnothing(canonical) && continue
        _check_remap_conflict(result, key, canonical, val)
        result[canonical] = val
        delete!(result, key)
    end
    return result
end

# Vector{Pair} path: rebuild pairs list preserving stable ordering.
# Two-pass: first collect all values by canonical key, then emit deduplicated pairs.
function _remap_inputs(pairs::AbstractVector{<:Pair}, combined)
    # First pass: collect values keyed by canonical target, detect conflicts.
    canonical_vals = Dict{Any, Any}()
    for (key, val) in pairs
        canon = get(combined, key, nothing)
        target = isnothing(canon) ? key : canon
        prev = get(canonical_vals, target, nothing)
        if prev !== nothing && !isequal(prev, val)
            error("Conflict in input map: $key has value $val " *
                  "but $target already has value $prev. " *
                  "These must agree or only one should be specified.")
        end
        canonical_vals[target] = val
    end

    # Second pass: emit pairs in first-seen order, skipping duplicates.
    result = similar(pairs, 0)
    sizehint!(result, length(canonical_vals))
    emitted = Set{Any}()
    for (key, _) in pairs
        canon = get(combined, key, nothing)
        target = isnothing(canon) ? key : canon
        if target ∉ emitted
            push!(result, target => canonical_vals[target])
            push!(emitted, target)
        end
    end
    return result
end

# Helper: check for conflict when remapping a Dict entry.
function _check_remap_conflict(result::AbstractDict, key, canonical, val)
    existing = get(result, canonical, nothing)
    if existing !== nothing && !isequal(existing, val)
        error("Conflict in input map: eliminated symbol $key has value $val " *
              "but canonical symbol $canonical has value $existing. " *
              "These must agree or only one should be specified.")
    end
end

### Conversion Pipeline Helper ###

"""
    prepare_aliases_for_conversion(flatrs; eliminate_aliases=true, allow_constraints=false)

Process aliases on a flattened `ReactionSystem` for conversion. Hot path: when no aliases
are present, returns immediately with zero overhead (single `isempty` check).

When `eliminate_aliases=true` (default): calls `eliminate_aliases(flatrs)`.
When `eliminate_aliases=false` and `allow_constraints=true`: materializes unknown aliases
as algebraic constraints (`lhs ~ rhs`) via `@set!`. Parameter aliases always error.
When `eliminate_aliases=false` and `allow_constraints=false` (jump path): errors.
"""
function prepare_aliases_for_conversion(flatrs::ReactionSystem;
        eliminate_aliases::Bool = true,
        allow_constraints::Bool = false)
    !has_aliases(flatrs) && return flatrs

    if eliminate_aliases
        return Catalyst.eliminate_aliases(flatrs)
    end

    # Not eliminating — check constraints
    if !isempty(get_parameter_aliases(flatrs))
        error("Parameter aliases must be eliminated before conversion. " *
              "Use `eliminate_aliases=true` (default).")
    end

    if !allow_constraints
        error("This conversion path does not support algebraic constraints from " *
              "uneliminated aliases. Use `eliminate_aliases=true` (default).")
    end

    # Materialize unknown aliases as algebraic constraints (lhs ~ rhs kept as-is).
    # Only two fields change, so use @set! instead of full reconstruction.
    result = flatrs
    @set! result.eqs = CatalystEqType[get_eqs(flatrs); get_unknown_aliases(flatrs)]
    @set! result.aliases = Equation[]
    return result
end

# Helper: filter aliases to just unknown (species/variable) aliases.
function get_unknown_aliases(rs::ReactionSystem)
    ctx = _make_classification_context(rs)
    filter(get_aliases(rs)) do eq
        cls = alias_class(eq.lhs, ctx)
        cls in (AliasClass.Species, AliasClass.BCSpecies, AliasClass.Unknown)
    end
end

# Helper: filter aliases to just parameter aliases.
function get_parameter_aliases(rs::ReactionSystem)
    ctx = _make_classification_context(rs)
    filter(get_aliases(rs)) do eq
        cls = alias_class(eq.lhs, ctx)
        cls in (AliasClass.Parameter, AliasClass.ConstantSpecies)
    end
end
