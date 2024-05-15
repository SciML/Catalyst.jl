

################################## Reaction Complexes ####################################

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

function ReactionComplex(speciesids::Vector{Int},
                         speciesstoichs::Vector{V}) where {V <: Integer}
    (length(speciesids) == length(speciesstoichs)) ||
        error("Creating a complex with different number of species ids and associated stoichiometries.")
    ReactionComplex{V}(speciesids, speciesstoichs)
end

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
    (setindex!(rc.speciesids, t.speciesid, i...);
     setindex!(rc.speciesstoichs,
               t.speciesstoich, i...);
     rc)
end
function Base.isless(a::ReactionComplexElement, b::ReactionComplexElement)
    isless(a.speciesid, b.speciesid)
end
Base.Sort.defalg(::ReactionComplex) = Base.DEFAULT_UNSTABLE

############################### Network Properties ####################################

#! format: off
# Internal cache for various ReactionSystem calculated properties
Base.@kwdef mutable struct NetworkProperties{I <: Integer, V <: BasicSymbolic{Real}}
    isempty::Bool = true
    netstoichmat::Union{Matrix{Int}, SparseMatrixCSC{Int, Int}} = Matrix{Int}(undef, 0, 0)
    conservationmat::Matrix{I} = Matrix{I}(undef, 0, 0)
    col_order::Vector{Int} = Int[]
    rank::Int = 0
    nullity::Int = 0
    indepspecs::Set{V} = Set{V}()
    depspecs::Set{V} = Set{V}()
    conservedeqs::Vector{Equation} = Equation[]
    constantdefs::Vector{Equation} = Equation[]
    speciesmap::Dict{V, Int} = Dict{V, Int}()
    complextorxsmap::OrderedDict{ReactionComplex{Int}, Vector{Pair{Int, Int}}} = OrderedDict{ReactionComplex{Int},Vector{Pair{Int,Int}}}()
    complexes::Vector{ReactionComplex{Int}} = Vector{ReactionComplex{Int}}(undef, 0)
    incidencemat::Union{Matrix{Int}, SparseMatrixCSC{Int, Int}} = Matrix{Int}(undef, 0, 0)
    complexstoichmat::Union{Matrix{Int}, SparseMatrixCSC{Int, Int}} = Matrix{Int}(undef, 0, 0)
    complexoutgoingmat::Union{Matrix{Int}, SparseMatrixCSC{Int, Int}} = Matrix{Int}(undef, 0, 0)
    incidencegraph::Graphs.SimpleDiGraph{Int} = Graphs.DiGraph()
    linkageclasses::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, 0)
    deficiency::Int = 0
end
#! format: on

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

function reset!(nps::NetworkProperties{I, V}) where {I, V}
    nps.isempty && return
    nps.netstoichmat = Matrix{Int}(undef, 0, 0)
    nps.conservationmat = Matrix{I}(undef, 0, 0)
    empty!(nps.col_order)
    nps.rank = 0
    nps.nullity = 0
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
    nps.deficiency = 0

    # this needs to be last due to setproperty! setting it to false
    nps.isempty = true
    nothing
end

############################### Reaction Systems ####################################

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
- `defaults::Dict`, a dictionary mapping parameters to their default values and species to
  their default initial values.
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
struct ReactionSystem{V <: NetworkProperties} <:
       MT.AbstractTimeDependentSystem
    """The equations (reactions and algebraic/differential) defining the system."""
    eqs::Vector{CatalystEqType}
    """The Reactions defining the system. """
    rxs::Vector{Reaction}
    """Independent variable (usually time)."""
    iv::BasicSymbolic{Real}
    """Spatial independent variables"""
    sivs::Vector{BasicSymbolic{Real}}
    """All dependent (unknown) variables, species and non-species. Must not contain the
    independent variable."""
    unknowns::Vector{BasicSymbolic{Real}}
    """Dependent unknown variables representing species"""
    species::Vector{BasicSymbolic{Real}}
    """Parameter variables. Must not contain the independent variable."""
    ps::Vector{Any}
    """Maps Symbol to corresponding variable."""
    var_to_name::Dict{Symbol, Any}
    """Equations for observed variables."""
    observed::Vector{Equation}
    """The name of the system"""
    name::Symbol
    """Internal sub-systems"""
    systems::Vector
    """
    The default values to use when initial conditions and/or
    parameters are not supplied in `ODEProblem`.
    """
    defaults::Dict
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
    metadata::Any
    """
    complete: if a model `sys` is complete, then `sys.x` no longer performs namespacing.
    """
    complete::Bool

    # inner constructor is considered private and may change between non-breaking releases.
    function ReactionSystem(eqs, rxs, iv, sivs, unknowns, spcs, ps, var_to_name, observed,
                            name, systems, defaults, connection_type, nps, cls, cevs, devs,
                            metadata = nothing, complete = false; checks::Bool = true)
                            
        # Checks that all parameters have the appropriate Symbolics type.
        for p in ps
            (p isa Symbolics.BasicSymbolic) || error("Parameter $p is not a `BasicSymbolic`. This is required.")
        end

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
            if !all(u == 1.0 for u in ModelingToolkit.get_unit([unknowns; ps; iv]))
                for eq in eqs
                    (eq isa Equation) && check_units(eq)
                end
            end
        end

        rs = new{typeof(nps)}(eqs, rxs, iv, sivs, unknowns, spcs, ps, var_to_name, observed,
                              name, systems, defaults, connection_type, nps, cls, cevs,
                              devs, metadata, complete)
        checks && validate(rs)
        rs
    end
end

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
        T = typeof(MT.unwrap(A))
    end

    T
end

eqsortby(eq::CatalystEqType) = eq isa Reaction ? 1 : 2

function ReactionSystem(eqs, iv, unknowns, ps;
                        observed = Equation[],
                        systems = [],
                        name = nothing,
                        default_u0 = Dict(),
                        default_p = Dict(),
                        defaults = _merge(Dict(default_u0), Dict(default_p)),
                        connection_type = nothing,
                        checks = true,
                        networkproperties = nothing,
                        combinatoric_ratelaws = true,
                        balanced_bc_check = true,
                        spatial_ivs = nothing,
                        continuous_events = nothing,
                        discrete_events = nothing,
                        metadata = nothing)
               
    # Error checks
    name === nothing &&
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    sysnames = nameof.(systems)
    (length(unique(sysnames)) == length(sysnames)) || throw(ArgumentError("System names must be unique."))

    # Handle defaults values provided via optional arguments.
    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn("`default_u0` and `default_p` are deprecated. Use `defaults` instead.", :ReactionSystem, force = true)
    end
    defaults = MT.todict(defaults)
    defaults = Dict{Any, Any}(value(k) => value(v) for (k, v) in pairs(defaults))

    # Extracts independent variables (iv and sivs), dependent variables (species and variables)
    # and parameters. Sorts so that species comes before variables in unknowns vector.
    iv′ = value(iv)
    sivs′ = if spatial_ivs === nothing
        Vector{typeof(iv′)}()
    else
        value.(MT.scalarize(spatial_ivs))
    end
    unknowns′ = sort!(value.(MT.scalarize(unknowns)), by = !isspecies) 
    spcs = filter(isspecies, unknowns′)
    ps′ = value.(MT.scalarize(ps))

    # Checks that no (by Catalyst) forbidden symbols are used.
    allsyms = Iterators.flatten((ps′, unknowns′))
    if !all(sym -> getname(sym) ∉ forbidden_symbols_error, allsyms)
        error("Catalyst reserves the symbols $forbidden_symbols_error for internal use. Please do not use these symbols as parameters or unknowns/species.")
    end

    # Handles reactions and equations. Sorts so that reactions are before equaions in the equations vector.
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
    # Adds their (potential) default values to the defaults vector.
    var_to_name = Dict()
    MT.process_variables!(var_to_name, defaults, unknowns′)
    MT.process_variables!(var_to_name, defaults, ps′)
    MT.collect_var_to_name!(var_to_name, eq.lhs for eq in observed)
    #
    # Computes network properties.
    nps = if networkproperties === nothing
        NetworkProperties{Int, get_speciestype(iv′, unknowns′, systems)}()
    else
        networkproperties
    end

    # Creates the continious and discrete callbacks.
    ccallbacks = MT.SymbolicContinuousCallbacks(continuous_events)
    dcallbacks = MT.SymbolicDiscreteCallbacks(discrete_events)

    ReactionSystem(eqs′, rxs, iv′, sivs′, unknowns′, spcs, ps′, var_to_name, observed, name,
                   systems, defaults, connection_type, nps, combinatoric_ratelaws,
                   ccallbacks, dcallbacks, metadata; checks = checks)
end

function ReactionSystem(rxs::Vector, iv = Catalyst.DEFAULT_IV; kwargs...)
    make_ReactionSystem_internal(rxs, iv, Vector{Num}(), Vector{Num}(); kwargs...)
end

# search the symbolic expression for parameters or unknowns
# and save in ps and us respectively. vars is used to cache results
function findvars!(ps, us, exprtosearch, ivs, vars)
    MT.get_variables!(vars, exprtosearch)
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

# Called internally (whether DSL-based or programmtic model creation is used). 
# Creates a sorted reactions + equations vector, also ensuring reaction is first in this vector.
# Extracts potential species, variables, and parameters from the input (if not provided as part of 
# the model creation) and creates the corresponding vectors. 
# While species are ordered before variables in the unknowns vector, this ordering is not imposed here,
# but carried out at a later stage.
function make_ReactionSystem_internal(rxs_and_eqs::Vector, iv, us_in, ps_in; spatial_ivs = nothing, 
                                      continuous_events = [], discrete_events = [], observed = [], kwargs...)

    # Filters away any potential obervables from `states` and `spcs`.
    obs_vars = [obs_eq.lhs for obs_eq in observed]
    us_in = filter(u -> !any(isequal(u, obs_var) for obs_var in obs_vars), us_in)
    
    # Creates a combined iv vector (iv and sivs). This is used later in the function (so that 
    # independent variables can be exluded when encountered quantities are added to `us` and `ps`).
    t = value(iv)
    ivs = Set([t])
    if (spatial_ivs !== nothing)
        for siv in (MT.scalarize(spatial_ivs))
            push!(ivs, value(siv))
        end
    end

    # Initialises the new unknowns and parameter vectors.
    # Preallocates the `vars` set, which is used by `findvars!`
    us = OrderedSet{eltype(us_in)}(us_in)
    ps = OrderedSet{eltype(ps_in)}(ps_in)
    vars = OrderedSet()

    # Extracts the reactions and equations from the combined reactions + equations input vector.
    all(eq -> eq isa Union{Reaction, Equation}, rxs_and_eqs)
    rxs = Reaction[eq for eq in rxs_and_eqs if eq isa Reaction]
    eqs = Equation[eq for eq in rxs_and_eqs if eq isa Equation]

    # Loops through all reactions, adding encountered quantities to the unknown and parameter vectors.
    # Starts by looping through substrates + products only (so these are added to the vector first).
    # Next, the otehr components of reactions (e.g. rates and stoichiometries) are added.
    for rx in rxs
        for reactants in (rx.substrates, rx.products), spec in reactants
            MT.isparameter(spec) ? push!(ps, spec) : push!(us, spec)
        end
    end
    for rx in rxs
        # Adds all quantitites encountered in the reaction's rate.
        findvars!(ps, us, rx.rate, ivs, vars)

        # Extracts all quantitites encountered within stoichiometries.
        for stoichiometry in (rx.substoich, rx.prodstoich), sym in stoichiometry
            (sym isa Symbolic) && findvars!(ps, us, sym, ivs, vars)
        end

        # Will appear here: add stuff from nosie scaling.
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

    # Loops through all events, adding encountered quantities to the unknwon and parameter vectors.
    find_event_vars!(ps, us, continuous_events, ivs, vars)      
    find_event_vars!(ps, us, discrete_events, ivs, vars)  

    # Converts the found unknowns and parameters to vectors.
    usv = collect(us)
    psv = collect(ps)

    # Passes the processed input into the next `ReactionSystem` call.    
    ReactionSystem(fulleqs, t, usv, psv; spatial_ivs, continuous_events, discrete_events, observed, kwargs...)
end

function ReactionSystem(iv; kwargs...)
    ReactionSystem(Reaction[], iv, [], []; kwargs...)
end

# Loops through all events in an supplied event vector, adding all unknowns and parameters found in
# its condition and affect functions to their respective vectors (`ps` and `us`).
function find_event_vars!(ps, us, events::Vector, ivs, vars)
    foreach(event -> find_event_vars!(ps, us, event, ivs, vars), events)
end
# For a single event, adds quantitites from its condition and affect expression(s) to `ps` and `us`.
# Applies `findvars!` to the event's condition (`event[1])` and affec (`event[2]`).
function find_event_vars!(ps, us, event, ivs, vars)
    findvars!(ps, us, event[1], ivs, vars)
    findvars!(ps, us, event[2], ivs, vars)
end

"""
    remake_ReactionSystem_internal(rs::ReactionSystem; 
        default_reaction_metadata::Vector{Pair{Symbol, T}} = Vector{Pair{Symbol, Any}}()) where {T}

Takes a `ReactionSystem` and remakes it, returning a modified `ReactionSystem`. Modifications depend
on which additional arguments are provided. The input `ReactionSystem` is not mutated. Updating
default reaction metadata is currently the only supported feature.

Arguments:
- `rs::ReactionSystem`: The `ReactionSystem` which you wish to remake.
- `default_reaction_metadata::Vector{Pair{Symbol, T}}`: A vector with default `Reaction` metadata values.
    Each metadata in each `Reaction` of the updated `ReactionSystem` will have the value desiganted in
    `default_reaction_metadata` (however, `Reaction`s that already have that metadata designated will not
    have their value updated).
"""
function remake_ReactionSystem_internal(rs::ReactionSystem;  default_reaction_metadata = [])
    rs = set_default_metadata(rs;  default_reaction_metadata)
    return rs
end

# For a `ReactionSystem`, updates all `Reaction`'s default metadata.
function set_default_metadata(rs::ReactionSystem;  default_reaction_metadata = [])
    # Updates reaction metadata for for reactions in this specific system.
    eqtransform(eq) = eq isa Reaction ? set_default_metadata(eq, default_reaction_metadata) : eq
    updated_equations = map(eqtransform, get_eqs(rs))
    @set! rs.eqs = updated_equations
    @set! rs.rxs = Reaction[rx for rx in updated_equations if rx isa Reaction]
    
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
    missing_metadata = filter(md -> !in(md[1], entry[1] for entry in rx.metadata), default_metadata)
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
    return remake_ReactionSystem_internal(rs, default_reaction_metadata = [:noise_scaling => noise_scaling])
end


"""
    isspatial(rn::ReactionSystem)

Returns whether `rn` has any spatial independent variables (i.e. is a spatial network).
"""
isspatial(rn::ReactionSystem) = !isempty(get_sivs(rn))

####################### ModelingToolkit inherited accessors #############################

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
for details. Can be overriden via passing `combinatoric_ratelaws` to `convert` or the
`*Problem` functions.
"""
get_combinatoric_ratelaws(sys::ReactionSystem) = getfield(sys, :combinatoric_ratelaws)

MT.get_continuous_events(sys::ReactionSystem) = getfield(sys, :continuous_events)

# need a custom equations since ReactionSystem.eqs are a mix of Reactions and Equations
function MT.equations(sys::ReactionSystem)
    ivs = independent_variables(sys)
    eqs = get_eqs(sys)
    systems = get_systems(sys)
    if !isempty(systems)
        eqs = CatalystEqType[eqs;
                             reduce(vcat, MT.namespace_equations.(systems, (ivs,));
                                    init = Any[])]
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

"""
    combinatoric_ratelaws(sys::ReactionSystem)

Returns the effective (default) `combinatoric_ratelaw` value for a compositional system,
calculated by taking the logical or of each component `ReactionSystem`. Can be overriden
during calls to `convert` of problem constructors.
"""
function combinatoric_ratelaws(sys::ReactionSystem)
    crl = get_combinatoric_ratelaws(sys)
    subsys = Iterators.filter(s -> s isa ReactionSystem, get_systems(sys))
    mapreduce(combinatoric_ratelaws, |, subsys; init = crl)
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

######################## Other accessors ##############################

"""
has_noise_scaling(reaction::Reaction)

Checks whether a specific reaction has the metadata field `noise_scaing`. If so, returns `true`, else
returns `false`.

Arguments:
- `reaction`: The reaction for which we wish to check.

Example:
```julia
reaction = @reaction k, 0 --> X, [noise_scaling=0.0]
has_noise_scaling(reaction)
"""
function has_noise_scaling(reaction::Reaction)
    return hasmetadata(reaction, :noise_scaling)
end

"""
get_noise_scaling(reaction::Reaction)

Returns the noise_scaling metadata from a specific reaction.

Arguments:
- `reaction`: The reaction for which we wish to retrive all metadata.

Example:
```julia
reaction = @reaction k, 0 --> X, [noise_scaling=0.0]
get_noise_scaling(reaction)
"""
function get_noise_scaling(reaction::Reaction)
    if has_noise_scaling(reaction)
        return getmetadata(reaction, :noise_scaling)
    else
        error("Attempts to access noise_scaling metadata field for a reaction which does not have a value assigned for this metadata.")
    end
end


# These are currently considered internal, but can be used by public accessor functions like get_noise_scaling.

"""
    getmetadata_dict(reaction::Reaction)

Retrives the `ImmutableDict` containing all of the metadata associated with a specific reaction.

Arguments:
- `reaction`: The reaction for which we wish to retrive all metadata.

Example:
```julia
reaction = @reaction k, 0 --> X, [description="Production reaction"]
getmetadata_dict(reaction)
```
"""
function getmetadata_dict(reaction::Reaction)
    return reaction.metadata
end

"""
    hasmetadata(reaction::Reaction, md_key::Symbol)

Checks if a `Reaction` have a certain metadata field. If it does, returns `true` (else returns `false`).

Arguments:
- `reaction`: The reaction for which we wish to check if a specific metadata field exist.
- `md_key`: The metadata for which we wish to check existence of.

Example:
```julia
reaction = @reaction k, 0 --> X, [description="Production reaction"]
hasmetadata(reaction, :description)
```
"""
function hasmetadata(reaction::Reaction, md_key::Symbol)
    return any(isequal(md_key, entry[1]) for entry in getmetadata_dict(reaction))
end

"""
getmetadata(reaction::Reaction, md_key::Symbol)

Retrives a certain metadata value from a `Reaction`. If the metadata does not exists, throws an error.

Arguments:
- `reaction`: The reaction for which we wish to retrive a specific metadata value.
- `md_key`: The metadata for which we wish to retrive.

Example:
```julia
reaction = @reaction k, 0 --> X, [description="Production reaction"]
getmetadata(reaction, :description)
```
"""
function getmetadata(reaction::Reaction, md_key::Symbol)
    if !hasmetadata(reaction, md_key) 
        error("The reaction does not have a metadata field $md_key. It does have the following metadata fields: $(keys(getmetadata_dict(reaction))).")
    end
    metadata = getmetadata_dict(reaction)
    return metadata[findfirst(isequal(md_key, entry[1]) for entry in getmetadata_dict(reaction))][2]
end

####################### dependency graph utilities ########################

# determine which unknowns a reaction depends on
function ModelingToolkit.get_variables!(deps::Set, rx::Reaction, variables)
    (rx.rate isa Symbolic) && get_variables!(deps, rx.rate, variables)
    for s in rx.substrates
        # parametric stoichiometry means may have a parameter as a substrate
        any(isequal(s), variables) && push!(deps, s)
    end
    deps
end

# determine which species a reaction modifies
function ModelingToolkit.modified_unknowns!(munknowns, rx::Reaction, sts::Set)
    for (species, stoich) in rx.netstoich
        (species in sts) && push!(munknowns, species)
    end
    munknowns
end

function ModelingToolkit.modified_unknowns!(munknowns, rx::Reaction, sts::AbstractVector)
    for (species, stoich) in rx.netstoich
        any(isequal(species), sts) && push!(munknowns, species)
    end
    munknowns
end

########################## Compositional Tooling ###########################
function getsubsystypes!(typeset::Set{Type}, sys::T) where {T <: MT.AbstractSystem}
    push!(typeset, T)
    for subsys in get_systems(sys)
        getsubsystypes!(typeset, subsys)
    end
    typeset
end
function getsubsystypes(sys)
    typeset = Set{Type}()
    getsubsystypes!(typeset, sys)
    typeset
end

"""
    Catalyst.flatten(rs::ReactionSystem)

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

    # right now only NonlinearSystems and ODESystems can be handled as subsystems
    subsys_types = getsubsystypes(rs)
    allowed_types = (ReactionSystem, NonlinearSystem, ODESystem)
    all(T -> any(T .<: allowed_types), subsys_types) ||
        error("flattening is currently only supported for subsystems mixing ReactionSystems, NonlinearSystems and ODESystems.")

    ReactionSystem(equations(rs), get_iv(rs), unknowns(rs), parameters(rs);
                   observed = MT.observed(rs),
                   name,
                   defaults = MT.defaults(rs),
                   checks = false,
                   combinatoric_ratelaws = combinatoric_ratelaws(rs),
                   balanced_bc_check = false,
                   spatial_ivs = get_sivs(rs),
                   continuous_events = MT.continuous_events(rs),
                   discrete_events = MT.discrete_events(rs))
end

"""
    ModelingToolkit.extend(sys::AbstractSystem, rs::ReactionSystem; name::Symbol=nameof(sys))

Extends the indicated [`ReactionSystem`](@ref) with another `AbstractSystem`.

Notes:
- The `AbstractSystem` being added in must be an `ODESystem`, `NonlinearSystem`,
  or `ReactionSystem` currently.
- Returns a new `ReactionSystem` and does not modify `rs`.
- By default, the new `ReactionSystem` will have the same name as `sys`.
"""
function ModelingToolkit.extend(sys::MT.AbstractSystem, rs::ReactionSystem;
                                name::Symbol = nameof(sys))
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
    defs = merge(get_defaults(rs), get_defaults(sys)) # prefer `sys`
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
                   defaults = defs,
                   checks = false,
                   combinatoric_ratelaws,
                   balanced_bc_check = false,
                   spatial_ivs = sivs,
                   continuous_events,
                   discrete_events)
end
