### ReactionComplex Structures ###

"""
$(TYPEDEF)
One reaction complex element

# Fields
$(FIELDS)
"""
struct ReactionComplexElement{T}
    """The integer id of the species representing this element."""
    speciesid::Int
    """The stoichiometric coefficient of this species."""
    speciesstoich::T
end

"""
$(TYPEDEF)
One reaction complex.

# Fields
$(FIELDS)
"""
struct ReactionComplex{V <: Integer} <: AbstractVector{ReactionComplexElement{V}}
    """The integer ids of all species participating in this complex."""
    speciesids::Vector{Int}
    """The stoichiometric coefficients of all species participating in this complex."""
    speciesstoichs::Vector{V}

    function ReactionComplex{V}(speciesids::Vector{Int},
            speciesstoichs::Vector{V}) where {V <: Integer}
        new{V}(speciesids, speciesstoichs)
    end
end

# Special constructor.
function ReactionComplex(speciesids::Vector{Int},
        speciesstoichs::Vector{V}) where {V <: Integer}
    (length(speciesids) == length(speciesstoichs)) ||
        error("Creating a complex with different number of species ids and associated stoichiometries.")
    ReactionComplex{V}(speciesids, speciesstoichs)
end

# Defines base function overloads for `ReactionComplex`.
function (==)(a::ReactionComplex{V}, b::ReactionComplex{V}) where {V <: Integer}
    (a.speciesids == b.speciesids) &&
        (a.speciesstoichs == b.speciesstoichs)
end

function hash(rc::ReactionComplex, h::UInt)
    Base.hash(rc.speciesids, Base.hash(rc.speciesstoichs, h))
end

Base.size(rc::ReactionComplex) = size(rc.speciesids)
Base.length(rc::ReactionComplex) = length(rc.speciesids)

function Base.getindex(rc::ReactionComplex, i...)
    ReactionComplexElement(getindex(rc.speciesids, i...), getindex(rc.speciesstoichs, i...))
end

function Base.setindex!(rc::ReactionComplex, t::ReactionComplexElement, i...)
    setindex!(rc.speciesids, t.speciesid, i...)
    setindex!(rc.speciesstoichs, t.speciesstoich, i...)
    rc
end

function Base.isless(a::ReactionComplexElement, b::ReactionComplexElement)
    isless(a.speciesid, b.speciesid)
end

Base.Sort.defalg(::ReactionComplex) = Base.DEFAULT_UNSTABLE

### NetworkProperties Structure ###
const __UNINITIALIZED_CONSERVED_CONSTS = unwrap(only(@parameters __UNINITIALIZED[1:1]))

#! format: off
# Internal cache for various ReactionSystem calculated properties
Base.@kwdef mutable struct NetworkProperties{I <: Integer, V <: SymbolicT}
    """Indicates if the network properties have been computed yet. true if not yet computed."""
    isempty::Bool = true
    """The network stoichiometric matrix. Rows correspond to species, columns to reactions."""
    netstoichmat::Union{Matrix{Int}, SparseMatrixCSC{Int, Int}} = Matrix{Int}(undef, 0, 0)
    """The conservation matrix. Rows correspond to conservation laws, columns to species."""
    conservationmat::Matrix{I} = Matrix{I}(undef, 0, 0)
    """The cycle matrix."""
    cyclemat::Matrix{I} = Matrix{I}(undef, 0, 0)
    """The column order of the stoichiometric matrix after row echelon form reduction. First
    `rank` entries are the indices of independent species, next `nullity` entries are the
    indices of dependent species."""
    col_order::Vector{Int} = Int[]
    """The rank of the stoichiometric matrix, i.e. number of independent species."""
    rank::Int = 0
    """The nullity of the stoichiometric matrix, i.e. number of dependent species."""
    nullity::Int = 0
    """The independent species."""
    indepspecs::Set{V} = Set{V}()
    """The dependent species."""
    depspecs::Set{V} = Set{V}()
    """The conserved equations in the form dependent_species = conserved_constant - ..."""
    conservedeqs::Vector{Equation} = Equation[]
    """The definitions of the conserved constants in the form conserved_constant = dependent_species + ..."""
    constantdefs::Vector{Equation} = Equation[]
    """The conserved constant symbolic vector, or a default value if not yet initialized."""
    conservedconst::SymbolicT = __UNINITIALIZED_CONSERVED_CONSTS 
    """Map from symbolics for each species to their index in the species vector."""
    speciesmap::Dict{V, Int} = Dict{V, Int}()
    complextorxsmap::OrderedDict{ReactionComplex{Int}, Vector{Pair{Int, Int}}} = OrderedDict{ReactionComplex{Int},Vector{Pair{Int,Int}}}()
    complexes::Vector{ReactionComplex{Int}} = Vector{ReactionComplex{Int}}(undef, 0)
    incidencemat::Union{Matrix{Int}, SparseMatrixCSC{Int, Int}} = Matrix{Int}(undef, 0, 0)
    complexstoichmat::Union{Matrix{Int}, SparseMatrixCSC{Int, Int}} = Matrix{Int}(undef, 0, 0)
    complexoutgoingmat::Union{Matrix{Int}, SparseMatrixCSC{Int, Int}} = Matrix{Int}(undef, 0, 0)
    incidencegraph::Graphs.SimpleDiGraph{Int} = Graphs.DiGraph()
    linkageclasses::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, 0)
    stronglinkageclasses::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, 0)
    terminallinkageclasses::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, 0)
    checkedrobust::Bool = false
    robustspecies::Vector{Int} = Vector{Int}(undef, 0)
    deficiency::Int = -1
end
#! format: on

# Defines base function overloads for `NetworkProperties`.
function Base.show(io::IO, nps::NetworkProperties)
    if (nps.conservationmat !== nothing)
        println(io, "Conserved Equations: ")
        foreach(eq -> println(io, eq), nps.conservedeqs)
        println()
    end
end
Base.isempty(nps::NetworkProperties) = getfield(nps, :isempty)
function Base.setproperty!(nps::NetworkProperties, sym::Symbol, x)
    (sym !== :isempty) && setfield!(nps, :isempty, false)
    setfield!(nps, sym, x)
end

# Resets computed properties.
function reset!(nps::NetworkProperties{I, V}) where {I, V}
    nps.isempty && return
    nps.netstoichmat = Matrix{Int}(undef, 0, 0)
    nps.conservationmat = Matrix{I}(undef, 0, 0)
    nps.cyclemat = Matrix{Int}(undef, 0, 0)
    empty!(nps.col_order)
    nps.rank = 0
    nps.nullity = 0
    nps.conservedconst = __UNINITIALIZED_CONSERVED_CONSTS
    empty!(nps.indepspecs)
    empty!(nps.depspecs)
    empty!(nps.conservedeqs)
    empty!(nps.constantdefs)
    empty!(nps.speciesmap)
    empty!(nps.complextorxsmap)
    empty!(nps.complexes)
    nps.incidencemat = Matrix{Int}(undef, 0, 0)
    nps.complexstoichmat = Matrix{Int}(undef, 0, 0)
    nps.complexoutgoingmat = Matrix{Int}(undef, 0, 0)
    nps.incidencegraph = Graphs.DiGraph()
    empty!(nps.linkageclasses)
    empty!(nps.stronglinkageclasses)
    empty!(nps.terminallinkageclasses)
    nps.deficiency = -1
    empty!(nps.robustspecies)
    nps.checkedrobust = false

    # this needs to be last due to setproperty! setting it to false
    nps.isempty = true
    nothing
end

function initialized_conserved(nps::NetworkProperties)
    nps.conservedconst !== __UNINITIALIZED_CONSERVED_CONSTS
end

### ReactionSystem Constructor Functions ###

# Used to sort the reaction/equation vector as reactions first, equations second.
eqsortby(eq::CatalystEqType) = eq isa Reaction ? 1 : 2

# Figures out a type.
function get_speciestype(iv, unknowns, systems)
    T = Nothing
    !isempty(unknowns) && (T = typeof(first(unknowns)))

    if !isempty(systems)
        for sys in Iterators.filter(s -> s isa ReactionSystem, systems)
            sts = MT.unknowns(sys)
            if !isempty(sts)
                T = typeof(first(sts))
                break
            end
        end
    end

    if T <: Nothing
        @variables A($iv)
        T = typeof(unwrap(A))
    end

    T
end

# search the symbolic expression for parameters or unknowns
# and save in ps and us respectively. vars is used to cache results
function findvars!(ps, us, exprtosearch, ivs, vars)
    MT.get_variables!(vars, exprtosearch; is_atomic = pre_is_not_atomic)
    for var in vars
        (var ∈ ivs) && continue
        if MT.isparameter(var)
            push!(ps, var)
        else
            push!(us, var)
        end
    end
    empty!(vars)
end
# Special dispatch for equations, applied `findvars!` to left-hand and right-hand sides.
function findvars!(ps, us, eq_to_search::Equation, ivs, vars)
    findvars!(ps, us, eq_to_search.lhs, ivs, vars)
    findvars!(ps, us, eq_to_search.rhs, ivs, vars)
end
# Special dispatch for Vectors (applies it to each vector element).
function findvars!(ps, us, exprs_to_search::Vector, ivs, vars)
    foreach(exprtosearch -> findvars!(ps, us, exprtosearch, ivs, vars), exprs_to_search)
end

# Used to ensure that `X`, not `Pre(X)` is sxtracted from expressions like `Pre(x)`.
pre_is_not_atomic(x) = SymbolicUtils.default_is_atomic(x) && !(iscall(x) && operation(x) isa Pre)

# Loops through all events in an supplied event vector, adding all unknowns and parameters found in
# its condition and affect functions to their respective vectors (`ps` and `us`).
function find_event_vars!(ps, us, events::Vector, ivs, vars)
    foreach(event -> find_event_vars!(ps, us, event, ivs, vars), events)
end
# For a single event, adds quantities from its condition and affect expression(s) to `ps` and `us`.
# Applies `findvars!` to the event's condition (`event[1])` and affec (`event[2]`).
# Two dispatches required, and the event can be given as a MTK structure of a Pair of symbolic expressions/equations.
function find_event_vars!(ps, us, event::Pair, ivs, vars)
    findvars!(ps, us, event[1], ivs, vars)
    findvars!(ps, us, event[2], ivs, vars)
end
function find_event_vars!(ps, us, event::MT.AbstractCallback, ivs, vars)
    findvars!(ps, us, event.conditions, ivs, vars)
    findvars!(ps, us, event.affect.affect, ivs, vars)
end

### ReactionSystem Structure ###

"""
WARNING!!!

The following variable is used to check that code that should be updated when the `ReactionSystem`
fields are updated has in fact been updated. Do not just blindly update this without first checking
all such code and updating it appropriately (e.g. serialization). Please use a search for
`reactionsystem_fields` throughout the package to ensure all places which should be updated, are updated.
"""
# Constant storing all reaction system fields (in order). Used to check whether the `ReactionSystem`
# structure have been updated (in the `reactionsystem_uptodate_check` function).
const reactionsystem_fields = (
    :eqs, :rxs, :iv, :sivs, :unknowns, :species, :ps, :var_to_name,
    :observed, :name, :systems, :initial_conditions, :connection_type,
    :networkproperties, :combinatoric_ratelaws, :continuous_events,
    :discrete_events, :metadata, :complete, :parent)

"""
$(TYPEDEF)

A system of chemical reactions.

# Fields
$(FIELDS)

# Example
Continuing from the example in the [`Reaction`](@ref) definition:
```julia
# simple constructor that infers species and parameters
@named rs = ReactionSystem(rxs, t)

# allows specification of species and parameters
@named rs = ReactionSystem(rxs, t, [A,B,C,D], k)
```

Keyword Arguments:
- `observed::Vector{Equation}`, equations specifying observed variables.
- `systems::Vector{AbstractSystems}`, vector of sub-systems. Can be `ReactionSystem`s,
  `ODESystem`s, or `NonlinearSystem`s.
- `name::Symbol`, the name of the system (must be provided, or `@named` must be used).
- `initial_conditions::SymmapT`, a dictionary mapping parameters and species to their initial
  values.
- `checks = true`, boolean for whether to check units.
- `networkproperties = NetworkProperties()`, cache for network properties calculated via API
  functions.
- `combinatoric_ratelaws = true`, sets the default value of `combinatoric_ratelaws` used in
  calls to `convert` or calling various problem types with the `ReactionSystem`.
- `balanced_bc_check = true`, sets whether to check that BC species appearing in reactions
  are balanced (i.e appear as both a substrate and a product with the same stoichiometry).

Notes:
- ReactionSystems currently do rudimentary unit checking, requiring that all species have
  the same units, and all reactions have rate laws with units of (species units) / (time
  units). Unit checking can be disabled by passing the keyword argument `checks=false`.
"""
struct ReactionSystem{V <: NetworkProperties} <: MT.AbstractSystem
    """The equations (reactions and algebraic/differential) defining the system."""
    eqs::Vector{CatalystEqType}
    """The Reactions defining the system. """
    rxs::Vector{Reaction}
    """Independent variable (usually time)."""
    iv::SymbolicT
    """Spatial independent variables"""
    sivs::Vector{SymbolicT}
    """All dependent (unknown) variables, species and non-species. Must not contain the
    independent variable."""
    unknowns::Vector{SymbolicT}
    """Dependent unknown variables representing species"""
    species::Vector{SymbolicT}
    """Parameter variables. Must not contain the independent variable."""
    ps::Vector{SymbolicT}
    """Maps Symbol to corresponding variable."""
    var_to_name::Dict{Symbol, SymbolicT}
    """Equations for observed variables."""
    observed::Vector{Equation}
    """The name of the system"""
    name::Symbol
    """Internal sub-systems"""
    systems::Vector
    """
    The initial values to use when initial conditions and/or
    parameters are not supplied in `ODEProblem`.
    """
    initial_conditions::SymmapT
    """Type of the system"""
    connection_type::Any
    """`NetworkProperties` object that can be filled in by API functions. INTERNAL -- not
    considered part of the public API."""
    networkproperties::V
    """Sets whether to use combinatoric scalings in rate laws. true by default."""
    combinatoric_ratelaws::Bool
    """
    continuous_events: A `Vector{SymbolicContinuousCallback}` that model events.
    The integrator will use root finding to guarantee that it steps at each zero crossing.
    """
    continuous_events::Vector{MT.SymbolicContinuousCallback}
    """
    discrete_events: A `Vector{SymbolicDiscreteCallback}` that models events. Symbolic
    analog to `SciMLBase.DiscreteCallback` that executes an affect when a given condition is
    true at the end of an integration step.
    """
    discrete_events::Vector{MT.SymbolicDiscreteCallback}
    """
    Metadata for the system, to be used by downstream packages.
    """
    metadata::MT.MetadataT
    """
    complete: if a model `sys` is complete, then `sys.x` no longer performs namespacing.
    """
    complete::Bool
    """
    The hierarchical parent system before simplification that MTK now seems to require for
    hierarchical namespacing to work in indexing.
    """
    parent::Any

    # inner constructor is considered private and may change between non-breaking releases.
    function ReactionSystem(eqs, rxs, iv, sivs, unknowns, spcs, ps, var_to_name, observed,
            name, systems, defaults, connection_type, nps, cls, cevs, devs,
            metadata, complete = false, parent = nothing; checks::Bool = true)

        # unit checks are for ODEs and Reactions only currently
        nonrx_eqs = Equation[eq for eq in eqs if eq isa Equation]
        if checks && isempty(sivs)
            check_variables(unknowns, iv)
            check_parameters(ps, iv)
            nonrx_eqs = Equation[eq for eq in eqs if eq isa Equation]
            !isempty(nonrx_eqs) && check_equations(nonrx_eqs, iv)
            !isempty(cevs) && check_equations(equations(cevs), iv)
        end

        if isempty(sivs) && (checks == true || (checks & MT.CheckUnits) > 0)
            if !all(unitless_symvar(sym) for sym in [unknowns; ps; iv])
                for eq in eqs
                    (eq isa Equation) && check_units(eq)
                end
            end
        end

        # Checks that no (non-reaction) equation contains a differential w.r.t. a species.
        for eq in eqs
            (eq isa Reaction) && continue
            (hasnode(is_species_diff, eq.lhs) || hasnode(is_species_diff, eq.rhs)) &&
                error("An equation ($eq) contains a differential with respect to a species. This is currently not supported. If this is a functionality you require, please raise an issue on the Catalyst GitHub page and we can consider the best way to implement it.")
        end
        rs = new{typeof(nps)}(
            eqs, rxs, iv, sivs, unknowns, spcs, ps, var_to_name, observed,
            name, systems, defaults, connection_type, nps, cls, cevs,
            devs, metadata, complete, parent)
        checks && validate(rs)
        rs
    end
end

# Checks if a symbolic expression contains a differential with respect to a species (either directly
# or somewhere within the differential expression).
function is_species_diff(expr)
    Symbolics.is_derivative(expr) || return false
    return hasnode(ex -> (ex isa SymbolicT) && isspecies(ex) && !isbc(ex), expr)
end

# Four-argument constructor. Permits additional inputs as optional arguments.
# Calls the full constructor.
function ReactionSystem(eqs, iv, unknowns, ps;
        observed = Equation[],
        systems = [],
        name = nothing,
        initial_conditions = SymmapT(),
        connection_type = nothing,
        checks = true,
        networkproperties = nothing,
        combinatoric_ratelaws = true,
        balanced_bc_check = true,
        spatial_ivs = nothing,
        continuous_events = nothing,
        discrete_events = nothing,
        metadata = MT.MetadataT())

    # Error checks
    name === nothing &&
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    sysnames = nameof.(systems)
    (length(unique(sysnames)) == length(sysnames)) ||
        throw(ArgumentError("System names must be unique."))

    # Process initial_conditions to unwrap Num wrappers.
    initial_conditions = SymmapT(value(entry[1]) => value(entry[2]) for entry in initial_conditions)

    # handles "bindings". Recently introduced in MTK/Symbolics, explicit Catalyst support need to be
    # implemented. Left empty for now.
    bindings = MT.SymmapT()

    # Extracts independent variables (iv and sivs), dependent variables (species and variables)
    # and parameters. Sorts so that species comes before variables in unknowns vector.
    iv′ = unwrap(iv)
    sivs′ = if spatial_ivs === nothing
        Vector{typeof(iv′)}()
    else
        unwrap.(spatial_ivs)
    end

    unknowns′ = isempty(unknowns) ? SymbolicT[] : sort!(unwrap.(unknowns), by = !isspecies)
    spcs = filter(isspecies, unknowns′)
    ps′ = isempty(ps) ? SymbolicT[] : unwrap.(ps)

    # Checks that no (by Catalyst) forbidden symbols are used.
    allsyms = Iterators.flatten((ps′, unknowns′))
    if !all(sym -> getname(sym) ∉ forbidden_symbols_error, allsyms)
        error("Catalyst reserves the symbols $forbidden_symbols_error for internal use. Please do not use these symbols as parameters or unknowns/species.")
    end

    # Handles reactions and equations. Sorts so that reactions are before equations in the equations vector.
    eqs′ = CatalystEqType[eq for eq in eqs]
    sort!(eqs′; by = eqsortby)
    rxs = Reaction[rx for rx in eqs if rx isa Reaction]

    # Additional error checks.
    if any(MT.isparameter, unknowns′)
        psts = filter(MT.isparameter, unknowns′)
        throw(ArgumentError("Found one or more parameters among the unknowns; this is not allowed. Move: $psts to be parameters."))
    end
    if any(isconstant, unknowns′)
        csts = filter(isconstant, unknowns′)
        throw(ArgumentError("Found one or more constant species among the unknowns; this is not allowed. Move: $csts to be parameters."))
    end
    # If there are BC species, check they are balanced in their reactions.
    if balanced_bc_check && any(isbc, unknowns′)
        for rx in eqs
            if (rx isa Reaction) && !isbcbalanced(rx)
                throw(ErrorException("BC species must be balanced, appearing as a substrate and product with the same stoichiometry. Please fix reaction: $rx"))
            end
        end
    end

    # Adds all unknowns/parameters to the `var_to_name` vector.
    # Adds their (potential) initial values to the initial_conditions dictionary.
    var_to_name = Dict{Symbol, SymbolicT}()
    MT.process_variables!(var_to_name, initial_conditions, bindings, unknowns′)
    MT.process_variables!(var_to_name, initial_conditions, bindings, ps′)
    MT.collect_var_to_name!(var_to_name, convert(Vector{SymbolicT}, [eq.lhs for eq in observed]))

    # Computes network properties.
    nps = if networkproperties === nothing
        NetworkProperties{Int, get_speciestype(iv′, unknowns′, systems)}()
    else
        networkproperties
    end

    # Creates the continuous and discrete events.
    continuous_events = create_symbolic_events(MT.SymbolicContinuousCallback, continuous_events)
    discrete_events = create_symbolic_events(MT.SymbolicDiscreteCallback, discrete_events)

    # handles system metadata.
    metadata = make_metadata(metadata)

    ReactionSystem(
        eqs′, rxs, iv′, sivs′, unknowns′, spcs, ps′, var_to_name, observed, name,
        systems, initial_conditions, connection_type, nps, combinatoric_ratelaws,
        continuous_events, discrete_events, metadata; checks = checks)
end

# Handles that events can be a single event or a vector.
create_symbolic_events(type, events::Vector) = [create_symbolic_event(type, event) for event in events]
create_symbolic_events(type, event) = [create_symbolic_event(type, event)]
create_symbolic_events(type, event::Nothing) = []

# Converts an input event into a form which ModelingToolkit can handle.
create_symbolic_event(type, event::MT.AbstractCallback) = event
create_symbolic_event(type, event) = type(event)

# Reformats the metadata in a format that MTK can handle. Long-term MTK will expose a public function
# that we can (and should) use instead.
function make_metadata(metadata)
    if isempty(metadata)
        metadata = MT.MetadataT()
    elseif metadata isa MT.MetadataT
        metadata = metadata
    else
        meta = MT.MetadataT()
        for kvp in metadata
            meta = Base.ImmutableDict(meta, kvp)
        end
        metadata = meta
    end
    metadata = MT.refreshed_metadata(metadata)
end

# Two-argument constructor (reactions/equations and time variable).
# Calls the `make_ReactionSystem_internal`, which in turn calls the four-argument constructor.
function ReactionSystem(rxs::Vector, iv = Catalyst.DEFAULT_IV; kwargs...)
    make_ReactionSystem_internal(rxs, iv, [], []; kwargs...)
end

# One-argument constructor. Creates an empty `ReactionSystem` from a time independent variable only.
function ReactionSystem(iv; kwargs...)
    ReactionSystem(Reaction[], iv, [], []; kwargs...)
end

# Called internally (whether DSL-based or programmatic model creation is used).
# Creates a sorted reactions + equations vector, also ensuring reaction is first in this vector.
# Extracts potential species, variables, and parameters from the input (if not provided as part of
# the model creation) and creates the corresponding vectors.
# While species are ordered before variables in the unknowns vector, this ordering is not imposed here,
# but carried out at a later stage.
function make_ReactionSystem_internal(rxs_and_eqs::Vector, iv, us_in, ps_in;
        spatial_ivs = nothing, continuous_events = [], discrete_events = [],
        observed = [], kwargs...)

    # Error if any observables have been declared a species or variable
    obs_vars = Set(obs_eq.lhs for obs_eq in observed)
    any(in(obs_vars), us_in) &&
        error("Found an observable in the list of unknowns. This is not allowed.")

    # Creates a combined iv vector (iv and sivs). This is used later in the function (so that
    # independent variables can be excluded when encountered quantities are added to `us` and `ps`).
    t = unwrap(iv)
    ivs = Set([t])
    if (spatial_ivs !== nothing)
        for siv in (spatial_ivs)
            push!(ivs, unwrap(siv))
        end
    end

    # Initialises the new unknowns and parameter vectors.
    # Preallocates the `vars` set, which is used by `findvars!`
    us = OrderedSet{SymbolicT}(us_in)
    ps = OrderedSet{SymbolicT}(ps_in)
    vars = OrderedSet{SymbolicT}()

    # Extracts the reactions and equations from the combined reactions + equations input vector.
    all(eq -> eq isa Union{Reaction, Equation}, rxs_and_eqs)
    rxs = Reaction[eq for eq in rxs_and_eqs if eq isa Reaction]
    eqs = Equation[eq for eq in rxs_and_eqs if eq isa Equation]

    # Loops through all reactions, adding encountered quantities to the unknown and parameter vectors.
    for rx in rxs
        MT.collect_vars!(us, ps, rx, t)
    end

    # Extracts any species, variables, and parameters that occur in (non-reaction) equations.
    # Creates the new reactions + equations vector, `fulleqs` (sorted reactions first, equations next).
    if !isempty(eqs)
        osys = ODESystem(eqs, iv; name = gensym())
        fulleqs = CatalystEqType[rxs; equations(osys)]
        union!(us, unknowns(osys))
        union!(ps, parameters(osys))
    else
        fulleqs = rxs
    end

    # get variables in subsystems with scope at this level
    for ssys in get(kwargs, :systems, [])
        MT.collect_scoped_vars!(us, ps, ssys, t)
    end

    # Loops through all events, adding encountered quantities to the unknown and parameter vectors.
    find_event_vars!(ps, us, continuous_events, ivs, vars)
    find_event_vars!(ps, us, discrete_events, ivs, vars)
    
    # Converts the found unknowns and parameters to vectors.
    usv = collect(us)

    new_ps = OrderedSet{SymbolicT}()
    for p in ps
        if iscall(p) && operation(p) === getindex
            par = arguments(p)[begin]
            if MT.symbolic_has_known_size(par) &&
               all(par[i] in ps for i in eachindex(par))
                push!(new_ps, par)
            else
                push!(new_ps, p)
            end
        else
            push!(new_ps, p)
        end
    end
    psv = collect(new_ps)

    # Passes the processed input into the next `ReactionSystem` call.
    ReactionSystem(fulleqs, t, usv, psv; spatial_ivs, continuous_events,
        discrete_events, observed, kwargs...)
end

### Base Function Dispatches ###

# Generic helper for comparing event vectors by content rather than identity.
# Works around MTK issue #3907 where SymbolicContinuousCallback/SymbolicDiscreteCallback
# don't support proper equality. Takes a matching function to compare individual events.
function events_equal(evts1, evts2, match_fn)
    length(evts1) != length(evts2) && return false
    isempty(evts1) && return true
    matched = falses(length(evts2))
    for evt1 in evts1
        idx = findfirst(j -> !matched[j] && match_fn(evt1, evts2[j]), eachindex(evts2))
        isnothing(idx) && return false
        matched[idx] = true
    end
    return true
end

# Compare two SymbolicAffects by their content.
function symbolic_affect_matches(aff1, aff2)
    issetequal(aff1.affect, aff2.affect) && issetequal(aff1.discrete_parameters, aff2.discrete_parameters)
end

# Compare two SymbolicContinuousCallbacks by their conditions and affects.
function continuous_event_matches(evt1, evt2)
    issetequal(evt1.conditions, evt2.conditions) &&
        symbolic_affect_matches(evt1.affect, evt2.affect) &&
        symbolic_affect_matches(evt1.affect_neg, evt2.affect_neg)
end

# Compare two SymbolicDiscreteCallbacks by their conditions and affects.
function discrete_event_matches(evt1, evt2)
    isequal(evt1.condition, evt2.condition) && symbolic_affect_matches(evt1.affect, evt2.affect)
end

# Convenience wrappers for continuous and discrete event comparison.
continuous_events_equal(evts1, evts2) = events_equal(evts1, evts2, continuous_event_matches)
discrete_events_equal(evts1, evts2) = events_equal(evts1, evts2, discrete_event_matches)

function debug_comparer(fun, prop1, prop2, propname; debug = false)
    if fun(prop1, prop2)
        return true
    else
        debug && println("Comparison was false for property: ", propname,
            "\n    Found: ", prop1, " vs ", prop2)
        return false
    end
end

"""
    isequivalent(rn1::ReactionSystem, rn2::ReactionSystem; ignorenames = true,
        debug = false)

Tests whether the underlying species, parameters and reactions are the same in the two
[`ReactionSystem`](@ref)s. Ignores the names of the systems in testing equality.

Notes:
- This function is primarily intended for testing purposes and is considered internal. It
    may change without a major release and does not guarantee true equality in notions that
    others may want.
- *Does not* currently simplify rates, so a rate of `A^2+2*A+1` would be considered
    different than `(A+1)^2`.
- `ignorenames = false` is used when checking equality of sub and parent systems.
- Does not check that `parent` systems are the same.
- Pass `debug = true` to print out the field that caused the two systems to be considered
    different.
"""
function isequivalent(rn1::ReactionSystem, rn2::ReactionSystem; ignorenames = true,
        debug = false)

    # metadata type fields
    if !ignorenames
        debug_comparer(==, nameof(rn1), nameof(rn2), "name"; debug) || return false
    end
    debug_comparer(==, get_combinatoric_ratelaws(rn1), get_combinatoric_ratelaws(rn2),
        "combinatoric_ratelaws"; debug) || return false
    debug_comparer(==, MT.iscomplete(rn1), MT.iscomplete(rn2), "complete"; debug) ||
        return false

    # symbolic variables and parameters
    debug_comparer(isequal, get_iv(rn1), get_iv(rn2), "ivs"; debug) || return false
    debug_comparer(issetequal, get_sivs(rn1), get_sivs(rn2), "sivs"; debug) || return false
    debug_comparer(issetequal, get_unknowns(rn1), get_unknowns(rn2), "unknowns"; debug) ||
        return false
    debug_comparer(issetequal, get_species(rn1), get_species(rn2), "species"; debug) ||
        return false
    debug_comparer(issetequal, get_ps(rn1), get_ps(rn2), "ps"; debug) || return false
    debug_comparer(
        issetequal, MT.get_initial_conditions(rn1), MT.get_initial_conditions(rn2), "defaults"; debug) ||
        return false

    # equations and reactions
    debug_comparer(
        issetequal, MT.get_observed(rn1), MT.get_observed(rn2), "observed"; debug) ||
        return false
    debug_comparer(issetequal, get_eqs(rn1), get_eqs(rn2), "eqs"; debug) || return false
    # Use custom event comparison functions to work around MTK issue #3907
    debug_comparer(continuous_events_equal, MT.get_continuous_events(rn1),
        MT.get_continuous_events(rn2), "cevents"; debug) || return false
    debug_comparer(discrete_events_equal, MT.get_discrete_events(rn1),
        MT.get_discrete_events(rn2), "devents"; debug) || return false

    # coupled systems
    if (length(get_systems(rn1)) != length(get_systems(rn2)))
        debug && println("Systems have different numbers of subsystems.")
        return false
    end
    # Use isequivalent recursively for subsystems instead of issetequal (which uses ==).
    # This is needed because == on ReactionSystems uses object identity.
    if !systems_are_equivalent(get_systems(rn1), get_systems(rn2); ignorenames, debug)
        debug && println("Comparison was false for property: systems",
            "\n    Found: ", get_systems(rn1), " vs ", get_systems(rn2))
        return false
    end

    true
end

# Helper for isequivalent: checks if two collections of ReactionSystems are equivalent.
# For each system in sys1, checks there's an equivalent one in sys2 (using isequivalent).
function systems_are_equivalent(sys1, sys2; ignorenames = true, debug = false)
    length(sys1) != length(sys2) && return false
    isempty(sys1) && return true
    matched = falses(length(sys2))
    for s1 in sys1
        found = false
        for (j, s2) in enumerate(sys2)
            matched[j] && continue  # already matched to another system
            if isequivalent(s1, s2; ignorenames, debug)
                matched[j] = true
                found = true
                break
            end
        end
        !found && return false
    end
    return true
end

### Basic `ReactionSystem`-specific Accessors ###

"""
    get_species(sys::ReactionSystem)

Return the current dependent variables that represent species in `sys` (toplevel system
only).
"""
get_species(sys::ReactionSystem) = getfield(sys, :species)
has_species(sys::ReactionSystem) = isdefined(sys, :species)

"""
    get_rxs(sys::ReactionSystem)

Return the system's `Reaction` vector (toplevel system only).
"""
get_rxs(sys::ReactionSystem) = getfield(sys, :rxs)
has_rxs(sys::ReactionSystem) = isdefined(sys, :rxs)

"""
    get_sivs(sys::ReactionSystem)

Return the current spatial ivs, if the system is non-spatial returns an empty vector.
"""
get_sivs(sys::ReactionSystem) = getfield(sys, :sivs)
has_sivs(sys::ReactionSystem) = isdefined(sys, :sivs)

"""
    get_networkproperties(sys::ReactionSystem)

Return the current network properties of `sys`.
"""
get_networkproperties(sys::ReactionSystem) = getfield(sys, :networkproperties)

"""
    get_combinatoric_ratelaws(sys::ReactionSystem)

Returns true if the default for the system is to rescale ratelaws, see
https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/introduction_to_catalyst/#Reaction-rate-laws-used-in-simulations
for details. Can be overridden via passing `combinatoric_ratelaws` to `convert` or the
`*Problem` functions.
"""
get_combinatoric_ratelaws(sys::ReactionSystem) = getfield(sys, :combinatoric_ratelaws)

"""
    combinatoric_ratelaws(sys::ReactionSystem)

Returns the effective (default) `combinatoric_ratelaw` value for a compositional system,
calculated by taking the logical or of each component `ReactionSystem`. Can be overridden
during calls to `convert` of problem constructors.
"""
function combinatoric_ratelaws(sys::ReactionSystem)
    crl = get_combinatoric_ratelaws(sys)
    subsys = Iterators.filter(s -> s isa ReactionSystem, get_systems(sys))
    mapreduce(combinatoric_ratelaws, |, subsys; init = crl)
end

# Gets sub systems that are also reaction systems.
# Used by several subsequent API functions.
function filter_nonrxsys(network)
    systems = get_systems(network)
    rxsystems = ReactionSystem[]
    for sys in systems
        (sys isa ReactionSystem) && push!(rxsystems, sys)
    end
    rxsystems
end

# Special species but which take a set of states and names spaces them according to another
# `ReactionSystem`.
# Used by `species(network)`.
function species(network::ReactionSystem, sts)
    [MT.renamespace(network, st) for st in sts]
end

"""
    species(network)

Given a [`ReactionSystem`](@ref), return a vector of all species defined in the system and
any subsystems that are of type `ReactionSystem`. To get the species and non-species
variables in the system and all subsystems, including non-`ReactionSystem` subsystems, uses
`unknowns(network)`.

Notes:
- If `ModelingToolkitBase.get_systems(network)` is non-empty will allocate.
"""
function species(network)
    sts = get_species(network)
    systems = filter_nonrxsys(network)
    isempty(systems) && return sts
    unique([sts; reduce(vcat, map(sys -> species(sys, species(sys)), systems))])
end

"""
    numspecies(network)

Return the total number of species within the given [`ReactionSystem`](@ref) and
subsystems that are `ReactionSystem`s.
"""
function numspecies(network)
    numspcs = length(get_species(network))
    for sys in get_systems(network)
        (sys isa ReactionSystem) && (numspcs += numspecies(sys))
    end
    numspcs
end

"""
    nonspecies(network)

Return the non-species variables within the network, i.e. those unknowns for which `isspecies
== false`.

Notes:
- Allocates a new array to store the non-species variables.
"""
function nonspecies(network)
    unknowns(network)[(numspecies(network) + 1):end]
end

"""
    speciesmap(network)

Given a [`ReactionSystem`](@ref), return a Dictionary mapping species that
participate in `Reaction`s to their index within [`species(network)`](@ref).
"""
function speciesmap(network)
    nps = get_networkproperties(network)
    if isempty(nps.speciesmap)
        nps.speciesmap = Dict(S => i for (i, S) in enumerate(species(network)))
    end
    nps.speciesmap
end

# get the non-bc, independent unknown variables and independent species, preserving their
# relative order in get_unknowns(rs). ASSUMES system has been validated to have no constant
# species as unknowns and is flattened.
function get_indep_sts(rs::ReactionSystem, remove_conserved = false)
    sts = get_unknowns(rs)
    nps = get_networkproperties(rs)
    indepsts = if remove_conserved
        filter(s -> ((s ∈ nps.indepspecs) || (!isspecies(s))) && (!isbc(s)), sts)
    else
        filter(s -> !isbc(s), sts)
    end
    indepsts, filter(isspecies, indepsts)
end

"""
    numparams(network)

Return the total number of parameters within the given system and all subsystems.
"""
function numparams(network)
    nps = length(get_ps(network))
    for sys in get_systems(network)
        nps += numparams(sys)
    end
    nps
end

"""
    paramsmap(network)

Given a [`ReactionSystem`](@ref), return a Dictionary mapping from all
parameters that appear within the system to their index within
`parameters(network)`.
"""
function paramsmap(network)
    Dict(p => i for (i, p) in enumerate(parameters(network)))
end

# used in the next function (`reactions(network)`).
function namespace_reactions(network::ReactionSystem)
    rxs = reactions(network)
    isempty(rxs) && return Reaction[]
    map(rx -> namespace_equation(rx, network), rxs)
end

"""
    reactions(network)

Given a [`ReactionSystem`](@ref), return a vector of all `Reactions` in the system.

Notes:
- If `ModelingToolkitBase.get_systems(network)` is not empty, will allocate.
"""
function reactions(network)
    rxs = get_rxs(network)
    systems = filter_nonrxsys(network)
    isempty(systems) && (return rxs)
    [rxs; reduce(vcat, namespace_reactions.(systems); init = Reaction[])]
end

"""
    numreactions(network)

Return the total number of reactions within the given [`ReactionSystem`](@ref)
and subsystems that are `ReactionSystem`s.
"""
function numreactions(network)
    nr = length(get_rxs(network))
    for sys in get_systems(network)
        (sys isa ReactionSystem) && (nr += numreactions(sys))
    end
    nr
end

"""
    has_nonreactions(network)

Check if the given `network` has any non-reaction equations such as ODEs or algebraic
equations.
"""
function has_nonreactions(network)
    numreactions(network) != length(equations(network))
end

"""
    nonreactions(network)

Return the non-reaction equations within the network (i.e. algebraic and differential equations).

Notes:
- Allocates a new array to store the non-species variables.
"""
function nonreactions(network)
    equations(network)[(numreactions(network) + 1):end]
end

"""
    reactionrates(network)

Given a [`ReactionSystem`](@ref), returns a vector of the symbolic reaction
rates for each reaction.
"""
function reactionrates(rn)
    [r.rate for r in reactions(rn)]
end

"""
    isspatial(rn::ReactionSystem)

Returns whether `rn` has any spatial independent variables (i.e. is a spatial network).
"""
isspatial(rn::ReactionSystem) = !isempty(get_sivs(rn))

### ModelingToolkit Function Dispatches ###

# Retrieves events.
MT.get_continuous_events(sys::ReactionSystem) = getfield(sys, :continuous_events)
# `MT.get_discrete_events(sys::ReactionSystem) = getfield(sys, :get_discrete_events)` should be added here.

# need a custom equations since ReactionSystem.eqs are a mix of Reactions and Equations
function MT.equations(sys::ReactionSystem)
    ivs = independent_variables(sys)
    eqs = get_eqs(sys)
    systems = get_systems(sys)
    if !isempty(systems)
        eqs = CatalystEqType[eqs; 
            reduce(vcat, MT.namespace_equations.(systems); init = Any[])]
        return sort!(eqs; by = eqsortby)
    end
    return eqs
end

function MT.unknowns(sys::ReactionSystem)
    sts = get_unknowns(sys)
    systems = get_systems(sys)
    if !isempty(systems)
        sts = unique!([sts; reduce(vcat, namespace_variables.(systems))])
        sort!(sts; by = !isspecies)
        return sts
    end
    return sts
end

function MT.complete(sys::ReactionSystem; flatten = true, kwargs...)
    newunknowns = OrderedSet{SymbolicT}()
    newparams = OrderedSet{SymbolicT}()
    iv = get_iv(sys)
    MT.collect_scoped_vars!(newunknowns, newparams, sys, iv; depth = -1)
    # don't update unknowns to not disturb `structural_simplify` order
    # `GlobalScope`d unknowns will be picked up and added there
    @set! sys.ps = unique!(vcat(get_ps(sys), collect(newparams)))
    if flatten
        newsys = Catalyst.flatten(sys)
        if MT.has_parent(newsys) && MT.get_parent(sys) === nothing
            @set! newsys.parent = complete(sys; split = false, flatten = false)
        end
        sys = newsys
    end
    isdefined(sys, :complete) ? (@set! sys.complete = true) : sys
end

### Network Matrix Representations ###

"""
    substoichmat(rn; sparse=false)

Returns the substrate stoichiometry matrix, ``S``, with ``S_{i j}`` the stoichiometric
coefficient of the ith substrate within the jth reaction.

Note:
- Set sparse=true for a sparse matrix representation
- Note that constant species are not considered substrates, but just components that modify
  the associated rate law.
"""
function substoichmat(::Type{SparseMatrixCSC{T, Int}},
        rn::ReactionSystem) where {T <: Number}
    Is = Int[]
    Js = Int[]
    Vs = T[]
    smap = speciesmap(rn)
    for (k, rx) in enumerate(reactions(rn))
        stoich = rx.substoich
        for (i, sub) in enumerate(rx.substrates)
            isconstant(sub) && continue
            push!(Js, k)
            push!(Is, smap[sub])
            push!(Vs, stoich[i])
        end
    end
    sparse(Is, Js, Vs, numspecies(rn), numreactions(rn))
end

function substoichmat(::Type{Matrix{T}}, rn::ReactionSystem) where {T <: Number}
    smap = speciesmap(rn)
    smat = zeros(T, numspecies(rn), numreactions(rn))
    for (k, rx) in enumerate(reactions(rn))
        stoich = rx.substoich
        for (i, sub) in enumerate(rx.substrates)
            isconstant(sub) && continue
            smat[smap[sub], k] = stoich[i]
        end
    end
    smat
end

function substoichmat(rn::ReactionSystem; sparse::Bool = false)
    isempty(get_systems(rn)) || error("substoichmat does not currently support subsystems.")
    T = reduce(promote_type, eltype(rx.substoich) for rx in reactions(rn))
    (T == Any) &&
        error("Stoichiometry matrices with symbolic stoichiometry are not supported")
    sparse ? substoichmat(SparseMatrixCSC{T, Int}, rn) : substoichmat(Matrix{T}, rn)
end

"""
    prodstoichmat(rn; sparse=false)

Returns the product stoichiometry matrix, ``P``, with ``P_{i j}`` the stoichiometric
coefficient of the ith product within the jth reaction.

Note:
- Set sparse=true for a sparse matrix representation
- Note that constant species are not treated as products, but just components that modify
  the associated rate law.
"""
function prodstoichmat(::Type{SparseMatrixCSC{T, Int}},
        rn::ReactionSystem) where {T <: Number}
    Is = Int[]
    Js = Int[]
    Vs = T[]
    smap = speciesmap(rn)
    for (k, rx) in enumerate(reactions(rn))
        stoich = rx.prodstoich
        for (i, prod) in enumerate(rx.products)
            isconstant(prod) && continue
            push!(Js, k)
            push!(Is, smap[prod])
            push!(Vs, stoich[i])
        end
    end
    sparse(Is, Js, Vs, numspecies(rn), numreactions(rn))
end

function prodstoichmat(::Type{Matrix{T}}, rn::ReactionSystem) where {T <: Number}
    smap = speciesmap(rn)
    pmat = zeros(T, numspecies(rn), numreactions(rn))
    for (k, rx) in enumerate(reactions(rn))
        stoich = rx.prodstoich
        for (i, prod) in enumerate(rx.products)
            isconstant(prod) && continue
            pmat[smap[prod], k] = stoich[i]
        end
    end
    pmat
end

function prodstoichmat(rn::ReactionSystem; sparse = false)
    isempty(get_systems(rn)) ||
        error("prodstoichmat does not currently support subsystems.")

    T = reduce(promote_type, eltype(rx.prodstoich) for rx in reactions(rn))
    (T == Any) &&
        error("Stoichiometry matrices with symbolic stoichiometry are not supported")
    sparse ? prodstoichmat(SparseMatrixCSC{T, Int}, rn) : prodstoichmat(Matrix{T}, rn)
end

# Used in `netstoichmat` function.
netstoichtype(::Vector{Pair{S, T}}) where {S, T} = T

"""
    netstoichmat(rn, sparse=false)

Returns the net stoichiometry matrix, ``N``, with ``N_{i j}`` the net stoichiometric
coefficient of the ith species within the jth reaction.

Notes:
- Set sparse=true for a sparse matrix representation
- Caches the matrix internally within `rn` so subsequent calls are fast.
- Note that constant species are not treated as reactants, but just components that modify
  the associated rate law. As such they do not contribute to the net stoichiometry matrix.
"""
function netstoichmat(::Type{SparseMatrixCSC{T, Int}},
        rn::ReactionSystem) where {T <: Number}
    Is = Int[]
    Js = Int[]
    Vs = Vector{T}()
    smap = speciesmap(rn)
    for (k, rx) in pairs(reactions(rn))
        for (spec, coef) in rx.netstoich
            isconstant(spec) && continue
            push!(Js, k)
            push!(Is, smap[spec])
            push!(Vs, coef)
        end
    end
    sparse(Is, Js, Vs, numspecies(rn), numreactions(rn))
end

function netstoichmat(::Type{Matrix{T}}, rn::ReactionSystem) where {T <: Number}
    smap = speciesmap(rn)
    nmat = zeros(T, numspecies(rn), numreactions(rn))
    for (k, rx) in pairs(reactions(rn))
        for (spec, coef) in rx.netstoich
            isconstant(spec) && continue
            nmat[smap[spec], k] = coef
        end
    end
    nmat
end

function netstoichmat(rn::ReactionSystem; sparse = false)
    isempty(get_systems(rn)) ||
        error("netstoichmat does not currently support subsystems, please create a flattened system before calling.")

    nps = get_networkproperties(rn)

    # if it is already calculated and has the right type
    !isempty(nps.netstoichmat) && (sparse == issparse(nps.netstoichmat)) &&
        (return nps.netstoichmat)

    # identify a common stoichiometry type
    T = reduce(promote_type, netstoichtype(rx.netstoich) for rx in reactions(rn))
    (T == Any) &&
        error("Stoichiometry matrices are not supported with symbolic stoichiometry.")

    if sparse
        nsmat = netstoichmat(SparseMatrixCSC{T, Int}, rn)
    else
        nsmat = netstoichmat(Matrix{T}, rn)
    end

    # only cache if it is integer
    if T == Int
        nps.netstoichmat = nsmat
    end

    nsmat
end

### General `ReactionSystem`-specific Functions ###

# Checks if the `ReactionSystem` structure have been updated without also updating the
# `reactionsystem_fields` constant. If this is the case, returns `false`. This is used in
# certain functionalities which would break if the `ReactionSystem` structure is updated without
# also updating these functionalities.
function reactionsystem_uptodate_check()
    if fieldnames(ReactionSystem) != reactionsystem_fields
        @warn "The `ReactionSystem` structure have been modified without this being taken into account in the functionality you are attempting to use. Please report this at https://github.com/SciML/Catalyst.jl/issues. Proceed with caution, as there might be errors in whichever functionality you are attempting to use."
    end
end

"""
    reset_networkproperties!(rn::ReactionSystem)

Clears the cache of various properties (like the netstoichiometry matrix). Use if such
properties need to be recalculated for some reason.
"""
function reset_networkproperties!(rn::ReactionSystem)
    reset!(get_networkproperties(rn))
    nothing
end

"""
    dependents(rx, network)

Given a [`Reaction`](@ref) and a [`ReactionSystem`](@ref), return a vector of the
*non-constant* species and variables the reaction rate law depends on. e.g., for

`k*W, 2X + 3Y --> 5Z + W`

the returned vector would be `[W(t),X(t),Y(t)]`.

Notes:
- Allocates
- Does not check for dependents within any subsystems.
- Constant species are not considered dependents since they are internally treated as
  parameters.
- If the rate expression depends on a non-species unknown variable that will be included in
  the dependents, i.e. in
  ```julia
  t = default_t()
  @parameters k
  @variables V(t)
  @species A(t) B(t) C(t)
  rx = Reaction(k*V, [A, B], [C])
  @named rs = ReactionSystem([rx], t)
  issetequal(dependents(rx, rs), [A,B,V]) == true
  ```
"""
function dependents(rx, network)
    if rx.rate isa Number
        return rx.substrates
    else
        rvars = get_variables(rx.rate, unknowns(network))
        return union!(rvars, rx.substrates)
    end
end

"""
    dependents(rx, network)

See documentation for [`dependents`](@ref).
"""
function dependants(rx, network)
    dependents(rx, network)
end

"""
isautonomous(rs::ReactionSystem)

Checks if a system is autonomous (i.e. no rate or equation depend on the independent variable(s)).
Example:
```julia
rs1 = @reaction_system
    (p,d), 0 <--> X
end
isautonomous(rs1) # Returns `true`.

rs2 = @reaction_system
    (p/t,d), 0 <--> X
end
isautonomous(rs2) # Returns `false`.
```
"""
function isautonomous(rs::ReactionSystem)
    # Get all variables occurring in reactions and equations.
    vars = Set{SymbolicT}()
    for eq in equations(rs)
        (eq isa Reaction) ? get_variables!(vars, eq.rate) : get_variables!(vars, eq)
    end

    # Checks for iv and spatial ivs
    (get_iv(rs) in vars) && return false
    any(in(vars), get_sivs(rs)) && return false
    return true
end

### `ReactionSystem` Remaking ###

"""
    remake_ReactionSystem_internal(rs::ReactionSystem;
        default_reaction_metadata::Vector{Pair{Symbol, T}} = Vector{Pair{Symbol, Any}}()) where {T}

Takes a `ReactionSystem` and remakes it, returning a modified `ReactionSystem`. Modifications depend
on which additional arguments are provided. The input `ReactionSystem` is not mutated. Updating
default reaction metadata is currently the only supported feature.

Arguments:
- `rs::ReactionSystem`: The `ReactionSystem` which you wish to remake.
- `default_reaction_metadata::Vector{Pair{Symbol, T}}`: A vector with default `Reaction` metadata values.
    Each metadata in each `Reaction` of the updated `ReactionSystem` will have the value designated in
    `default_reaction_metadata` (however, `Reaction`s that already have that metadata designated will not
    have their value updated).
"""
function remake_ReactionSystem_internal(rs::ReactionSystem; default_reaction_metadata = [])
    rs = set_default_metadata(rs; default_reaction_metadata)
    return rs
end

# For a `ReactionSystem`, updates all `Reaction`'s default metadata.
function set_default_metadata(rs::ReactionSystem; default_reaction_metadata = [])
    # Updates reaction metadata for reactions in this specific system.
    function eqtransform(eq)
        eq isa Reaction ? set_default_metadata(eq, default_reaction_metadata) : eq
    end
    updated_equations = map(eqtransform, get_eqs(rs))
    @set! rs.eqs = updated_equations
    @set! rs.rxs = Reaction[rx for rx in updated_equations if rx isa Reaction]

    # Special routine to handle `Reaction` metadata that can contain new symbolic variables.
    # Currently, `noise_scaling` is the only relevant metadata supported this way.
    drm_dict = Dict(default_reaction_metadata)
    if haskey(drm_dict, :noise_scaling)
        # Finds parameters, species, and variables in the noise scaling term.
        ns_expr = drm_dict[:noise_scaling]
        ns_syms = [unwrap(sym) for sym in get_variables(ns_expr)]
        ns_ps = Iterators.filter(MT.isparameter, ns_syms)
        ns_sps = Iterators.filter(Catalyst.isspecies, ns_syms)
        ns_vs = Iterators.filter(
            sym -> !Catalyst.isspecies(sym) &&
                   !MT.isparameter(sym), ns_syms)
        # Adds parameters, species, and variables to the `ReactionSystem`.
        @set! rs.ps = union(get_ps(rs), ns_ps)
        sps_new = union(get_species(rs), ns_sps)
        @set! rs.species = sps_new
        vs_old = @view get_unknowns(rs)[(length(get_species(rs)) + 1):end]
        @set! rs.unknowns = union(sps_new, vs_old, ns_vs)
    end

    # Updates reaction metadata for all its subsystems.
    new_sub_systems = similar(get_systems(rs))
    for (i, sub_system) in enumerate(get_systems(rs))
        new_sub_systems[i] = set_default_metadata(sub_system; default_reaction_metadata)
    end
    @set! rs.systems = new_sub_systems

    # Returns the updated system.
    return rs
end

# For a `Reaction`, adds missing default metadata values. Equations are passed back unmodified.
function set_default_metadata(rx::Reaction, default_metadata)
    missing_metadata = filter(
        md -> !in(md[1], entry[1] for entry in rx.metadata), default_metadata)
    updated_metadata = vcat(rx.metadata, missing_metadata)
    updated_metadata = convert(Vector{Pair{Symbol, Any}}, updated_metadata)
    return @set rx.metadata = updated_metadata
end
set_default_metadata(eq::Equation, default_metadata) = eq

"""
set_default_noise_scaling(rs::ReactionSystem, noise_scaling)

Creates an updated `ReactionSystem`. This is the old `ReactionSystem`, but each `Reaction` that does
not have a `noise_scaling` metadata have its noise_scaling metadata updated. The input `ReactionSystem`
is not mutated. Any subsystems of `rs` have their `noise_scaling` metadata updated as well.

Arguments:
- `rs::ReactionSystem`: The `ReactionSystem` which you wish to remake.
- `noise_scaling`: The updated noise scaling terms
"""
function set_default_noise_scaling(rs::ReactionSystem, noise_scaling)
    return remake_ReactionSystem_internal(
        rs, default_reaction_metadata = [:noise_scaling => noise_scaling])
end

### ReactionSystem Composing & Hierarchical Modelling ###

"""
    make_empty_network(; iv=DEFAULT_IV, name=gensym(:ReactionSystem))

Construct an empty [`ReactionSystem`](@ref). `iv` is the independent variable,
usually time, and `name` is the name to give the `ReactionSystem`.
"""
function make_empty_network(; iv = DEFAULT_IV, name = gensym(:ReactionSystem))
    ReactionSystem(Reaction[], iv, SymbolicT[], SymbolicT[]; name = name)
end

# Checks if a system is an allowed subsystem (i.e. no SDE parts and no jump).
is_allowed_subsystem(sys::ReactionSystem) = true
function is_allowed_subsystem(sys::System)
    return (isnothing(MT.get_noise_eqs(sys)) || isempty(MT.get_noise_eqs(sys))) &&
        (isnothing(MT.get_jumps(sys)) || isempty(MT.get_jumps(sys)))
end
# If neither a `ReactionSystem` or a `System`, it is something weird we do not know what it is.
is_allowed_subsystem(sys::MT.AbstractSystem) = false

"""
    ModelingToolkitBase.flatten(rs::ReactionSystem)

Merges all subsystems of the given [`ReactionSystem`](@ref) up into `rs`.

Notes:
- Returns a new `ReactionSystem` that represents the flattened system.
- All `Reaction`s within subsystems are namespaced and merged into the list of `Reactions`
  of `rs`. The merged list is then available as `reactions(rs)`.
- All algebraic and differential equations are merged in the equations of `rs`.
- Currently only `ReactionSystem`s, `NonlinearSystem`s and `ODESystem`s are supported as
  sub-systems when flattening.
- `rs.networkproperties` is reset upon flattening.
- The default value of `combinatoric_ratelaws` will be the logical or of all
  `ReactionSystem`s.
"""
function MT.flatten(rs::ReactionSystem; name = nameof(rs))
    isempty(get_systems(rs)) && return rs

    # right now we only guarantee tht certain types of systems work with flatten
    allowed_types = (ReactionSystem, NonlinearSystem, ODESystem)
    isnothing(get_systems(rs)) || all(is_allowed_subsystem, get_systems(rs)) ||
        error("flattening is currently only supported for subsystems mixing ReactionSystems, and Systems withour noise equations and jumps.")

    ReactionSystem(equations(rs), get_iv(rs), unknowns(rs), parameters(rs);
        observed = MT.observed(rs),
        name,
        initial_conditions = MT.initial_conditions(rs),
        checks = false,
        combinatoric_ratelaws = combinatoric_ratelaws(rs),
        balanced_bc_check = false,
        spatial_ivs = get_sivs(rs),
        continuous_events = MT.continuous_events(rs),
        discrete_events = MT.discrete_events(rs),
        metadata = MT.get_metadata(rs))
end

function complete_check(sys, method)
    if MT.iscomplete(sys)
        error("$method with one or more `ReactionSystem`s requires systems to not be marked complete, but system: $(MT.get_name(sys)) is marked complete.")
    end
    nothing
end

"""
    ModelingToolkitBase.compose(sys::ReactionSystem, systems::AbstractArray; name = nameof(sys))

Compose the indicated [`ReactionSystem`](@ref) with one or more `AbstractSystem`s.

Notes:
- The `AbstractSystem` being added in must be an `ODESystem`, `NonlinearSystem`,
  or `ReactionSystem` currently.
- Returns a new `ReactionSystem` and does not modify `rs`.
- By default, the new `ReactionSystem` will have the same name as `sys`.
"""
function MT.compose(sys::ReactionSystem, systems::AbstractArray; name = nameof(sys))
    complete_check(sys, "MT.compose")
    foreach(s -> complete_check(s, "MT.compose"), systems)

    nsys = length(systems)
    nsys == 0 && return sys
    @set! sys.name = name
    @set! sys.systems = [get_systems(sys); systems]
    newunknowns = OrderedSet{SymbolicT}()
    newparams = OrderedSet{SymbolicT}()
    iv = has_iv(sys) ? get_iv(sys) : nothing
    for ssys in systems
        MT.collect_scoped_vars!(newunknowns, newparams, ssys, iv)
    end

    if !isempty(newunknowns)
        @set! sys.unknowns = union(get_unknowns(sys), newunknowns)
        sort!(get_unknowns(sys), by = !isspecies)
        @set! sys.species = filter(isspecies, get_unknowns(sys))
    end

    if !isempty(newparams)
        @set! sys.ps = union(get_ps(sys), newparams)
    end

    return sys
end

"""
    ModelingToolkitBase.extend(sys::AbstractSystem, rs::ReactionSystem; name::Symbol=nameof(sys))

Extends the indicated [`ReactionSystem`](@ref) with another `AbstractSystem`.

Notes:
- The `AbstractSystem` being added in must be an `ODESystem`, `NonlinearSystem`,
  or `ReactionSystem` currently.
- Returns a new `ReactionSystem` and does not modify `rs`.
- By default, the new `ReactionSystem` will have the same name as `sys`.
"""
function MT.extend(sys::MT.AbstractSystem, rs::ReactionSystem;
        name::Symbol = nameof(sys))
    complete_check(sys, "MT.extend")
    complete_check(rs, "MT.extend")

    any(T -> sys isa T, (ReactionSystem, ODESystem, NonlinearSystem)) ||
        error("ReactionSystems can only be extended with ReactionSystems, ODESystems and NonlinearSystems currently. Received a $(typeof(sys)) system.")

    t = get_iv(rs)
    if MT.has_iv(sys)
        isequal(get_iv(sys), t) ||
            error("Extending ReactionSystem with iv, $(get_iv(rs)), with a system with iv, $(get_iv(sys)), this is not supported. Please ensure the `ivs` are the same.")
    end

    # generic system properties
    eqs = union(get_eqs(rs), get_eqs(sys))
    sts = union(get_unknowns(rs), get_unknowns(sys))
    ps = union(get_ps(rs), get_ps(sys))
    obs = union(get_observed(rs), get_observed(sys))
    syss = union(get_systems(rs), get_systems(sys))
    defs = merge(MT.get_initial_conditions(rs), MT.get_initial_conditions(sys)) # prefer `sys`
    continuous_events = union(MT.get_continuous_events(rs), MT.get_continuous_events(sys))
    discrete_events = union(MT.get_discrete_events(rs), MT.get_discrete_events(sys))

    # ReactionSystem specific properties
    if sys isa ReactionSystem
        combinatoric_ratelaws = Catalyst.get_combinatoric_ratelaws(sys) |
                                Catalyst.get_combinatoric_ratelaws(rs)
        sivs = union(get_sivs(sys), get_sivs(rs))
    else
        combinatoric_ratelaws = Catalyst.get_combinatoric_ratelaws(rs)
        sysivs = MT.has_ivs(sys) ? filter(!isequal(t), independent_variables(sys)) :
                 Vector{typeof(t)}()
        sivs = (length(sysivs) > 0) ? union(get_sivs(rs), sysivs) : get_sivs(rs)
    end

    ReactionSystem(eqs, t, sts, ps;
        observed = obs,
        systems = syss,
        name,
        initial_conditions = defs,
        checks = false,
        combinatoric_ratelaws,
        balanced_bc_check = false,
        spatial_ivs = sivs,
        continuous_events,
        discrete_events)
end

### Units Handling ###

"""
    validate(rs::ReactionSystem, info::String="")

Check that all species in the [`ReactionSystem`](@ref) have the same units, and
that the rate laws of all reactions reduce to units of (species units) / (time
units).

Notes:
- Does not check subsystems, constraint equations, or non-species variables.
"""
function validate(rs::ReactionSystem, info::String = "")
    specs = get_species(rs)

    # if there are no species we don't check units on the system
    isempty(specs) && return true

    specunits = get_unit(specs[1])
    validated = true
    for spec in specs
        if get_unit(spec) != specunits
            validated = false
            @warn(string("Species are expected to have units of ", specunits,
                " however, species ", spec, " has units ", get_unit(spec), "."))
        end
    end
    timeunits = get_unit(get_iv(rs))

    # no units for species, time or parameters then assume validated
    unitless = Base.get_extension(ModelingToolkitBase, :MTKDynamicQuantitiesExt).unitless
    if (specunits in (unitless, nothing)) && (timeunits in (unitless, nothing))
        all(unitless_symvar(p) for p in get_ps(rs)) && return true
    end

    rateunits = specunits / timeunits
    for rx in get_rxs(rs)
        rxunits = get_unit(rx.rate)
        for (i, sub) in enumerate(rx.substrates)
            rxunits *= get_unit(sub)^rx.substoich[i]
        end

        # Checks that the reaction's combined units is correct, if not, throws a warning.
        # Needs additional checks because for cases: (1.0^n) and (1.0^n1)*(1.0^n2).
        # These are not considered (be default) considered equal to `1.0` for unitless reactions.
        isequal(rxunits, rateunits) && continue
        if iscall(rxunits)
            unitless_exp(rxunits) && continue
            (operation(rxunits) == *) &&
                all(unitless_exp(arg) for arg in arguments(rxunits)) && continue
        end
        validated = false
        @warn(string(
            "Reaction rate laws are expected to have units of ", rateunits, " however, ",
            rx, " has units of ", rxunits, "."))
    end

    validated
end

# Checks if a unit consist of exponents with base 1 (and is this unitless).
unitless_exp(u) = iscall(u) && (operation(u) == ^) && (arguments(u)[1] == 1)

# Checks if a symbolic variable is unitless. Also accounts for callable parameters (for
# which `get_unit`'s` intended behaviour (or whether it should generate an error) is undefined: https://github.com/SciML/ModelingToolkit.jl/issues/3420).
function unitless_symvar(sym)
    return (sym isa Symbolics.CallAndWrap) || (MT.get_unit(sym) == 1)
end
