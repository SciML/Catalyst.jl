# Catalyst specific symbolics to support SBML
struct ParameterConstantSpecies end
struct VariableBCSpecies end
struct VariableSpecies end
Symbolics.option_to_metadata_type(::Val{:isconstantspecies}) = ParameterConstantSpecies
Symbolics.option_to_metadata_type(::Val{:isbcspecies}) = VariableBCSpecies
Symbolics.option_to_metadata_type(::Val{:isspecies}) = VariableSpecies

"""
    Catalyst.isconstant(s)

Tests if the given symbolic variable corresponds to a constant species.
"""
isconstant(s::Num) = isconstant(MT.value(s))
function isconstant(s)
    MT.getmetadata(s, ParameterConstantSpecies, false)
end

"""
    Catalyst.isbc(s)

Tests if the given symbolic variable corresponds to a boundary condition species.
"""
isbc(s::Num) = isbc(MT.value(s))
function isbc(s)
    MT.getmetadata(s, VariableBCSpecies, false)
end

"""
    isspecies(s)

Tests if the given symbolic variable corresponds to a chemical species.
"""
isspecies(s::Num) = isspecies(MT.value(s))
function isspecies(s)
    MT.getmetadata(s, VariableSpecies, false)
end

"""
    tospecies(s)

Convert the given symbolic variable to be a species by adding the isspecies metadata.

Notes:
- Will error if passed a parameter.
"""
function tospecies(s)
    MT.isparameter(s) &&
        error("Parameters can not be converted to species. Please pass a variable.")
    MT.setmetadata(s, VariableSpecies, true)
end

# true for species which shouldn't change from the reactions, including non-species
# variables
drop_dynamics(s) = isconstant(s) || isbc(s) || (!isspecies(s))

"""
$(TYPEDEF)

One chemical reaction.

# Fields
$(FIELDS)

# Examples

```julia
using Catalyst
@parameters k[1:20]
@variables t
@species A(t) B(t) C(t) D(t)
rxs = [Reaction(k[1], nothing, [A]),            # 0 -> A
       Reaction(k[2], [B], nothing),            # B -> 0
       Reaction(k[3],[A],[C]),                  # A -> C
       Reaction(k[4], [C], [A,B]),              # C -> A + B
       Reaction(k[5], [C], [A], [1], [2]),      # C -> A + A
       Reaction(k[6], [A,B], [C]),              # A + B -> C
       Reaction(k[7], [B], [A], [2], [1]),      # 2B -> A
       Reaction(k[8], [A,B], [A,C]),            # A + B -> A + C
       Reaction(k[9], [A,B], [C,D]),            # A + B -> C + D
       Reaction(k[10], [A], [C,D], [2], [1,1]), # 2A -> C + D
       Reaction(k[11], [A], [A,B], [2], [1,1]), # 2A -> A + B
       Reaction(k[12], [A,B,C], [C,D], [1,3,4], [2, 3]),          # A+3B+4C -> 2C + 3D
       Reaction(k[13], [A,B], nothing, [3,1], nothing),           # 3A+B -> 0
       Reaction(k[14], nothing, [A], nothing, [2]),               # 0 -> 2A
       Reaction(k[15]*A/(2+A), [A], nothing; only_use_rate=true), # A -> 0 with custom rate
       Reaction(k[16], [A], [B]; only_use_rate=true),             # A -> B with custom rate.
       Reaction(k[17]*A*exp(B), [C], [D], [2], [1]),              # 2C -> D with non constant rate.
       Reaction(k[18]*B, nothing, [B], nothing, [2]),             # 0 -> 2B with non constant rate.
       Reaction(k[19]*t, [A], [B]),                                # A -> B with non constant rate.
       Reaction(k[20]*t*A, [B,C], [D],[2,1],[2])                  # 2A +B -> 2C with non constant rate.
  ]
```

Notes:
- `nothing` can be used to indicate a reaction that has no reactants or no products.
  In this case the corresponding stoichiometry vector should also be set to `nothing`.
- The three-argument form assumes all reactant and product stoichiometric coefficients
  are one.
"""
struct Reaction{S, T}
    """The rate function (excluding mass action terms)."""
    rate::Any
    """Reaction substrates."""
    substrates::Vector
    """Reaction products."""
    products::Vector
    """The stoichiometric coefficients of the reactants."""
    substoich::Vector{T}
    """The stoichiometric coefficients of the products."""
    prodstoich::Vector{T}
    """The net stoichiometric coefficients of all species changed by the reaction."""
    netstoich::Vector{Pair{S, T}}
    """
    `false` (default) if `rate` should be multiplied by mass action terms to give the rate law.
    `true` if `rate` represents the full reaction rate law.
    """
    only_use_rate::Bool
end

"""
    isvalidreactant(s)

Test if a species is valid as a reactant (i.e. a species variable or a constant parameter).
"""
isvalidreactant(s) = MT.isparameter(s) ? isconstant(s) : (isspecies(s) && !isconstant(s))

function Reaction(rate, subs, prods, substoich, prodstoich;
                  netstoich = nothing, only_use_rate = false,
                  kwargs...)
    (isnothing(prods) && isnothing(subs)) &&
        throw(ArgumentError("A reaction requires a non-nothing substrate or product vector."))
    (isnothing(prodstoich) && isnothing(substoich)) &&
        throw(ArgumentError("Both substrate and product stochiometry inputs cannot be nothing."))

    if isnothing(subs)
        prodtype = typeof(value(first(prods)))
        subs = Vector{prodtype}()
        !isnothing(substoich) &&
            throw(ArgumentError("If substrates are nothing, substrate stoichiometries have to be so too."))
        substoich = typeof(prodstoich)()
    else
        subs = value.(subs)
    end
    allunique(subs) ||
        throw(ArgumentError("Substrates can not be repeated in the list provided to `Reaction`, please modify the stoichiometry for any repeated substrates instead."))
    S = eltype(substoich)

    if isnothing(prods)
        prods = Vector{eltype(subs)}()
        !isnothing(prodstoich) &&
            throw(ArgumentError("If products are nothing, product stoichiometries have to be so too."))
        prodstoich = typeof(substoich)()
    else
        prods = value.(prods)
    end
    allunique(prods) ||
        throw(ArgumentError("Products can not be repeated in the list provided to `Reaction`, please modify the stoichiometry for any repeated products instead."))
    T = eltype(prodstoich)

    # try to get a common type for stoichiometry, using Any if have Syms
    stoich_type = promote_type(S, T)
    if stoich_type <: Num
        stoich_type = Any
        substoich′ = Any[value(s) for s in substoich]
        prodstoich′ = Any[value(p) for p in prodstoich]
    else
        substoich′ = (S == stoich_type) ? substoich : convert.(stoich_type, substoich)
        prodstoich′ = (T == stoich_type) ? prodstoich : convert.(stoich_type, prodstoich)
    end

    if !(all(isvalidreactant, subs) && all(isvalidreactant, prods))
        badsts = union(filter(!isvalidreactant, subs), filter(!isvalidreactant, prods))
        throw(ArgumentError("""To be a valid substrate or product, non-constant species must be declared via @species, while constant species must be parameters with the isconstantspecies metadata. The following reactants do not follow this convention:\n $badsts"""))
    end

    ns = if netstoich === nothing
        get_netstoich(subs, prods, substoich′, prodstoich′)
    else
        (netstoich_stoichtype(netstoich) != stoich_type) ?
        convert.(stoich_type, netstoich) : netstoich
    end

    Reaction(value(rate), subs, prods, substoich′, prodstoich′, ns, only_use_rate)
end

# three argument constructor assumes stoichiometric coefs are one and integers
function Reaction(rate, subs, prods; kwargs...)
    sstoich = isnothing(subs) ? nothing : ones(Int, length(subs))
    pstoich = isnothing(prods) ? nothing : ones(Int, length(prods))
    Reaction(rate, subs, prods, sstoich, pstoich; kwargs...)
end

function print_rxside(io::IO, specs, stoich)
    # reactants/substrates
    if isempty(specs)
        print(io, "∅")
    else
        for (i, spec) in enumerate(specs)
            prspec = (MT.isparameter(spec) || (MT.operation(spec) == getindex)) ?
                     spec : MT.operation(spec)
            if isequal(stoich[i], one(stoich[i]))
                print(io, prspec)
            elseif istree(stoich[i])
                print(io, "(", stoich[i], ")*", prspec)
            else
                print(io, stoich[i], "*", prspec)
            end

            (i < length(specs)) && print(io, " + ")
        end
    end
    nothing
end

function Base.show(io::IO, rx::Reaction)
    print(io, rx.rate, ", ")
    print_rxside(io, rx.substrates, rx.substoich)
    arrow = rx.only_use_rate ? "⇒" : "-->"
    print(io, " ", arrow, " ")
    print_rxside(io, rx.products, rx.prodstoich)
end

function apply_if_nonempty(f, v)
    isempty(v) && return v
    s = similar(v)
    map!(f, s, v)
    s
end

function ModelingToolkit.namespace_equation(rx::Reaction, name; kw...)
    f = Base.Fix2(namespace_expr, name)
    rate = f(rx.rate)
    subs = apply_if_nonempty(f, rx.substrates)
    prods = apply_if_nonempty(f, rx.products)
    substoich = apply_if_nonempty(f, rx.substoich)
    prodstoich = apply_if_nonempty(f, rx.prodstoich)
    netstoich = if isempty(rx.netstoich)
        rx.netstoich
    else
        ns = similar(rx.netstoich)
        map!(n -> f(n[1]) => f(n[2]), ns, rx.netstoich)
    end
    Reaction(rate, subs, prods, substoich, prodstoich, netstoich, rx.only_use_rate)
end

netstoich_stoichtype(::Vector{Pair{S, T}}) where {S, T} = T

# calculates the net stoichiometry of a reaction as a vector of pairs (sub,substoich)
function get_netstoich(subs, prods, sstoich, pstoich)
    # stoichiometry as a Dictionary
    nsdict = Dict{Any, eltype(sstoich)}(sub => -sstoich[i] for (i, sub) in enumerate(subs))
    for (i, p) in enumerate(prods)
        coef = pstoich[i]
        @inbounds nsdict[p] = haskey(nsdict, p) ? nsdict[p] + coef : coef
    end

    # stoichiometry as a vector
    [el for el in nsdict if !_iszero(el[2])]
end

"""
    isbcbalanced(rx::Reaction)

True if any BC species in `rx` appears as a substrate and product with the same
stoichiometry.
"""
function isbcbalanced(rx::Reaction)
    # any substrate BC must be a product with the same stoichiometry
    for (sidx, sub) in enumerate(rx.substrates)
        if isbc(sub)
            pidx = findfirst(Base.Fix1(isequal, sub), rx.products)
            (pidx === nothing) && return false
            isequal(rx.prodstoich[pidx], rx.substoich[sidx]) || return false
        end
    end

    for prod in rx.products
        if isbc(prod)
            any(Base.Fix1(isequal, prod), rx.substrates) || return false
        end
    end

    true
end

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

const CatalystEqType = Union{Reaction, Equation}

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
    """All dependent (state) variables, species and non-species. Must not contain the
    independent variable."""
    states::Vector{BasicSymbolic{Real}}
    """Dependent state variables representing species"""
    species::Vector{BasicSymbolic{Real}}
    """Parameter variables. Must not contain the independent variable."""
    ps::Vector{BasicSymbolic{Real}}
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
    complete: if a model `sys` is complete, then `sys.x` no longer performs namespacing.
    """
    complete::Bool

    # inner constructor is considered private and may change between non-breaking releases.
    function ReactionSystem(eqs, rxs, iv, sivs, states, spcs, ps, var_to_name, observed,
                            name, systems, defaults, connection_type, nps, cls, cevs, devs,
                            complete::Bool = false; checks::Bool = true)

        # unit checks are for ODEs and Reactions only currently
        nonrx_eqs = Equation[eq for eq in eqs if eq isa Equation]
        if checks && isempty(sivs)
            check_variables(states, iv)
            check_parameters(ps, iv)
            nonrx_eqs = Equation[eq for eq in eqs if eq isa Equation]
            !isempty(nonrx_eqs) && check_equations(nonrx_eqs, iv)
            check_equations(equations(cevs), iv)
        end

        if isempty(sivs) && (checks == true || (checks & MT.CheckUnits) > 0)
            nonrx_eqs = Equation[eq for eq in eqs if eq isa Equation]
            MT.all_dimensionless([states; ps; iv]) || check_units(nonrx_eqs)
        end

        rs = new{typeof(nps)}(eqs, rxs, iv, sivs, states, spcs, ps, var_to_name, observed,
                              name, systems, defaults, connection_type, nps, cls, cevs,
                              devs, complete)
        checks && validate(rs)
        rs
    end

    # Copies a reaction system, but with the option of having some fields replaced
    function ReactionSystem(rs::ReactionSystem; eqs = rs.eqs, rxs = rs.rxs, iv = rs.iv, sivs = rs.sivs, states = rs.states, species = rs.species, ps = rs.ps, var_to_name = rs.var_to_name, observed = rs.observed, name = rs.name, systems = rs.systems, defaults = rs.defaults, connection_type = rs.connection_type, networkproperties = rs.networkproperties, combinatoric_ratelaws = rs.combinatoric_ratelaws, continuous_events = rs.continuous_events, discrete_events = rs.discrete_events, complete = rs.complete)
        new{typeof(networkproperties)}(eqs, rxs, ModelingToolkit.unwrap(iv), ModelingToolkit.unwrap.(sivs), ModelingToolkit.unwrap.(states), ModelingToolkit.unwrap.(species), ModelingToolkit.unwrap.(ps), var_to_name, observed, name, systems, defaults, connection_type, networkproperties, combinatoric_ratelaws, continuous_events, discrete_events, complete)
    end
end

function get_speciestype(iv, states, systems)
    T = Nothing
    !isempty(states) && (T = typeof(first(states)))

    if !isempty(systems)
        for sys in Iterators.filter(s -> s isa ReactionSystem, systems)
            sts = MT.states(sys)
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

function ReactionSystem(eqs, iv, states, ps;
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
                        discrete_events = nothing)
    name === nothing &&
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    sysnames = nameof.(systems)
    (length(unique(sysnames)) == length(sysnames)) ||
        throw(ArgumentError("System names must be unique."))

    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn("`default_u0` and `default_p` are deprecated. Use `defaults` instead.",
                     :ReactionSystem, force = true)
    end
    defaults = MT.todict(defaults)
    defaults = Dict{Any, Any}(value(k) => value(v) for (k, v) in pairs(defaults))

    iv′ = value(iv)
    sivs′ = if spatial_ivs === nothing
        Vector{typeof(iv′)}()
    else
        value.(MT.scalarize(spatial_ivs))
    end
    states′ = sort!(value.(MT.scalarize(states)), by = !isspecies) # species come first
    spcs = filter(isspecies, states′)
    ps′ = value.(MT.scalarize(ps))

    allsyms = Iterators.flatten((ps′, states′))
    all(sym -> getname(sym) ∉ forbidden_symbols_error, allsyms) ||
        error("Catalyst reserves the symbols $forbidden_symbols_error for internal use. Please do not use these symbols as parameters or states/species.")

    # sort Reactions before Equations
    eqs′ = CatalystEqType[eq for eq in eqs]
    sort!(eqs′; by = eqsortby)
    rxs = Reaction[rx for rx in eqs if rx isa Reaction]

    if any(MT.isparameter, states′)
        psts = filter(MT.isparameter, states′)
        throw(ArgumentError("Found one or more parameters among the states; this is not allowed. Move: $psts to be parameters."))
    end

    if any(isconstant, states′)
        csts = filter(isconstant, states′)
        throw(ArgumentError("Found one or more constant species among the states; this is not allowed. Move: $csts to be parameters."))
    end

    # if there are BC species, check they are balanced in their reactions
    if balanced_bc_check && any(isbc, states′)
        for rx in eqs
            if rx isa Reaction
                isbcbalanced(rx) ||
                    throw(ErrorException("BC species must be balanced, appearing as a substrate and product with the same stoichiometry. Please fix reaction: $rx"))
            end
        end
    end

    var_to_name = Dict()
    MT.process_variables!(var_to_name, defaults, states′)
    MT.process_variables!(var_to_name, defaults, ps′)
    MT.collect_var_to_name!(var_to_name, eq.lhs for eq in observed)

    nps = if networkproperties === nothing
        NetworkProperties{Int, get_speciestype(iv′, states′, systems)}()
    else
        networkproperties
    end

    ccallbacks = MT.SymbolicContinuousCallbacks(continuous_events)
    dcallbacks = MT.SymbolicDiscreteCallbacks(discrete_events)

    ReactionSystem(eqs′, rxs, iv′, sivs′, states′, spcs, ps′, var_to_name, observed, name,
                   systems, defaults, connection_type, nps, combinatoric_ratelaws,
                   ccallbacks, dcallbacks; checks = checks)
end

function ReactionSystem(rxs::Vector, iv = Catalyst.DEFAULT_IV; kwargs...)
    make_ReactionSystem_internal(rxs, iv, Vector{Num}(), Vector{Num}(); kwargs...)
end

# search the symbolic expression for parameters or states
# and save in ps and sts respectively. vars is used to cache results
function findvars!(ps, sts, exprtosearch, ivs, vars)
    MT.get_variables!(vars, exprtosearch)
    for var in vars
        (var ∈ ivs) && continue
        if MT.isparameter(var)
            push!(ps, var)
        else
            push!(sts, var)
        end
    end
    empty!(vars)
end

# Only used internally by the @reaction_network macro. Permits giving an initial order to
# the parameters, and then adds additional ones found in the reaction. Name could be
# changed.
function make_ReactionSystem_internal(rxs_and_eqs::Vector, iv, sts_in, ps_in;
                                      spatial_ivs = nothing, kwargs...)
    t = value(iv)
    ivs = Set([t])
    if (spatial_ivs !== nothing)
        for siv in (MT.scalarize(spatial_ivs))
            push!(ivs, value(siv))
        end
    end
    sts = OrderedSet{eltype(sts_in)}(sts_in)
    ps = OrderedSet{eltype(ps_in)}(ps_in)
    vars = OrderedSet()

    all(eq -> eq isa Union{Reaction, Equation}, rxs_and_eqs)
    rxs = Reaction[eq for eq in rxs_and_eqs if eq isa Reaction]
    eqs = Equation[eq for eq in rxs_and_eqs if eq isa Equation]

    # add species / parameters that are substrates / products first
    for rx in rxs, reactants in (rx.substrates, rx.products)
        for spec in reactants
            MT.isparameter(spec) ? push!(ps, spec) : push!(sts, spec)
        end
    end

    for rx in rxs
        findvars!(ps, sts, rx.rate, ivs, vars)
        for s in rx.substoich
            (s isa Symbolic) && findvars!(ps, sts, s, ivs, vars)
        end
        for p in rx.prodstoich
            (p isa Symbolic) && findvars!(ps, sts, p, ivs, vars)
        end
    end

    stsv = collect(sts)
    psv = collect(ps)

    if !isempty(eqs)
        osys = ODESystem(eqs, iv; name = gensym())
        fulleqs = CatalystEqType[rxs; equations(osys)]
        union!(stsv, states(osys))
        union!(psv, parameters(osys))
    else
        fulleqs = rxs
    end

    ReactionSystem(fulleqs, t, stsv, psv; spatial_ivs, kwargs...)
end

function ReactionSystem(iv; kwargs...)
    ReactionSystem(Reaction[], iv, [], []; kwargs...)
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

function MT.states(sys::ReactionSystem)
    sts = get_states(sys)
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

# get the non-bc, independent state variables and independent species, preserving their
# relative order in get_states(rs). ASSUMES system has been validated to have no constant
# species as states and is flattened.
function get_indep_sts(rs::ReactionSystem, remove_conserved = false)
    sts = get_states(rs)
    nps = get_networkproperties(rs)
    indepsts = if remove_conserved
        filter(s -> (s ∈ nps.indepspecs) && (!isbc(s)), sts)
    else
        filter(s -> !isbc(s), sts)
    end
    indepsts, filter(isspecies, indepsts)
end

######################## Conversion to ODEs/SDEs/jump, etc ##############################

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
function oderatelaw(rx; combinatoric_ratelaw = true)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate

    # if the stoichiometric coefficients are not integers error if asking to scale rates
    !all(s -> s isa Union{Integer, Symbolic}, substoich) &&
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

function assemble_oderhs(rs, ispcs; combinatoric_ratelaws = true, remove_conserved = false)
    nps = get_networkproperties(rs)
    species_to_idx = Dict(x => i for (i, x) in enumerate(ispcs))
    rhsvec = Any[0 for _ in ispcs]
    depspec_submap = if remove_conserved
        Dict(eq.lhs => eq.rhs for eq in nps.conservedeqs)
    else
        Dict()
    end

    for rx in get_rxs(rs)
        rl = oderatelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws)
        remove_conserved && (rl = substitute(rl, depspec_submap))
        for (spec, stoich) in rx.netstoich
            # dependent species don't get an ODE, so are skipped
            remove_conserved && (spec in nps.depspecs) && continue

            # constant or BC species also do not get equations
            drop_dynamics(spec) && continue

            i = species_to_idx[spec]
            if _iszero(rhsvec[i])
                if stoich isa Symbolic
                    rhsvec[i] = stoich * rl
                else
                    signedrl = (stoich > zero(stoich)) ? rl : -rl
                    rhsvec[i] = isone(abs(stoich)) ? signedrl : stoich * rl
                end
            else
                if stoich isa Symbolic
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
                        include_zero_odes = true, remove_conserved = false)
    rhsvec = assemble_oderhs(rs, ispcs; combinatoric_ratelaws, remove_conserved)
    if as_odes
        D = Differential(get_iv(rs))
        eqs = [Equation(D(x), rhs)
               for (x, rhs) in zip(ispcs, rhsvec) if (include_zero_odes || (!_iszero(rhs)))]
    else
        eqs = [Equation(0, rhs) for rhs in rhsvec if (include_zero_odes || (!_iszero(rhs)))]
    end
    eqs
end

# this doesn't work with constraint equations currently
function assemble_diffusion(rs, sts, ispcs, noise_scaling; combinatoric_ratelaws = true,
                            remove_conserved = false)
    # as BC species should ultimately get an equation, we include them in the noise matrix
    num_bcsts = count(isbc, get_states(rs))

    # we make a matrix sized by the number of reactions
    eqs = Matrix{Any}(undef, length(sts) + num_bcsts, length(get_rxs(rs)))
    eqs .= 0
    species_to_idx = Dict((x => i for (i, x) in enumerate(ispcs)))
    nps = get_networkproperties(rs)
    depspec_submap = if remove_conserved
        Dict(eq.lhs => eq.rhs for eq in nps.conservedeqs)
    else
        Dict()
    end

    for (j, rx) in enumerate(get_rxs(rs))
        rlsqrt = sqrt(abs(oderatelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws)))
        (noise_scaling !== nothing) && (rlsqrt *= noise_scaling[j])
        remove_conserved && (rlsqrt = substitute(rlsqrt, depspec_submap))

        for (spec, stoich) in rx.netstoich
            # dependent species don't get an equation
            remove_conserved && (spec in nps.depspecs) && continue

            # constant or BC species also do not get equations
            drop_dynamics(spec) && continue

            i = species_to_idx[spec]
            if stoich isa Symbolic
                eqs[i, j] = stoich * rlsqrt
            else
                signedrlsqrt = (stoich > zero(stoich)) ? rlsqrt : -rlsqrt
                eqs[i, j] = isone(abs(stoich)) ? signedrlsqrt : stoich * rlsqrt
            end
        end
    end
    eqs
end

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
function jumpratelaw(rx; combinatoric_ratelaw = true)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate
    if !only_use_rate
        coef = eltype(substoich) <: Number ? one(eltype(substoich)) : 1
        for (i, stoich) in enumerate(substoich)
            s = substrates[i]
            if stoich isa Symbolic
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
                              stateset = Set(states(rs)),
                              ivset = nothing)
```

True if a given reaction is of mass action form, i.e. `rx.rate` does not depend
on any chemical species that correspond to states of the system, and does not
depend explicitly on the independent variable (usually time).

# Arguments
- `rx`, the [`Reaction`](@ref).
- `rs`, a [`ReactionSystem`](@ref) containing the reaction.
- Optional: `rxvars`, `Variable`s which are not in `rxvars` are ignored as
  possible dependencies.
- Optional: `haveivdep`, `true` if the [`Reaction`](@ref) `rate` field
  explicitly depends on any independent variable (i.e. t or for spatial systems x,y,etc).
  If not set, will be automatically calculated.
- Optional: `stateset`, set of states which if the rxvars are within mean rx is
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
                      stateset = Set(get_states(rs)), ivset = nothing)

    # we define non-integer (i.e. float or symbolic) stoich to be non-mass action
    ((eltype(rx.substoich) <: Integer) && (eltype(rx.prodstoich) <: Integer)) ||
        return false

    # if no dependencies must be zero order
    (length(rxvars) == 0) && return true

    if (haveivdep === nothing)
        if isspatial(rs)
            if (ivset === nothing)
                ivs = Set(get_sivs(rs))
                push!(ivs, get_iv(rs))
                ivdep = any(var -> var ∈ ivs, rxvars)
            else
                ivdep = any(var -> var ∈ ivset, rxvars)
            end
        else
            ivdep = any(var -> isequal(get_iv(rs), var), rxvars)
        end
        ivdep && return false
    else
        haveivdep && return false
    end
    rx.only_use_rate && return false
    @inbounds for var in rxvars
        # not mass action if have a non-constant, variable species in the rate expression
        (var in stateset) && return false
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

# get_depgraph(rs)[i] is the list of reactions with rates depending on species changed by
# i'th reaction.
function get_depgraph(rs)
    jdeps = asgraph(rs)
    vdeps = variable_dependencies(rs)
    eqeq_dependencies(jdeps, vdeps).fadjlist
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

function assemble_jumps(rs; combinatoric_ratelaws = true)
    meqs = MassActionJump[]
    ceqs = ConstantRateJump[]
    veqs = VariableRateJump[]
    stateset = Set(get_states(rs))
    all(isspecies, stateset) ||
        error("Conversion to JumpSystem currently requires all states to be species.")
    rxvars = []

    isempty(get_rxs(rs)) &&
        error("Must give at least one reaction before constructing a JumpSystem.")

    # first we determine vrjs with an explicit time-dependent rate
    rxs = get_rxs(rs)
    isvrjvec = falses(length(rxs))
    havevrjs = false
    for (i, rx) in enumerate(rxs)
        empty!(rxvars)
        (rx.rate isa Symbolic) && get_variables!(rxvars, rx.rate)
        @inbounds for rxvar in rxvars
            if isequal(rxvar, get_iv(rs))
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

    for (i, rx) in enumerate(rxs)
        empty!(rxvars)
        (rx.rate isa Symbolic) && get_variables!(rxvars, rx.rate)

        isvrj = isvrjvec[i]
        if (!isvrj) && ismassaction(rx, rs; rxvars = rxvars, haveivdep = false,
                        stateset = stateset)
            push!(meqs, makemajump(rx; combinatoric_ratelaw = combinatoric_ratelaws))
        else
            rl = jumpratelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws)
            affect = Vector{Equation}()
            for (spec, stoich) in rx.netstoich
                # don't change species that are constant or BCs
                (!drop_dynamics(spec)) && push!(affect, spec ~ spec + stoich)
            end
            if isvrj
                push!(veqs, VariableRateJump(rl, affect))
            else
                push!(ceqs, ConstantRateJump(rl, affect))
            end
        end
    end
    vcat(meqs, ceqs, veqs)
end

# merge constraint components with the ReactionSystem components
# also handles removing BC and constant species
function addconstraints!(eqs, rs::ReactionSystem, ists, ispcs; remove_conserved = false)
    # if there are BC species, put them after the independent species
    rssts = get_states(rs)
    sts = any(isbc, rssts) ? vcat(ists, filter(isbc, rssts)) : ists
    ps = get_ps(rs)

    # make dependent species observables and add conservation constants as parameters
    if remove_conserved
        nps = get_networkproperties(rs)

        # add the conservation constants as parameters and set their values
        ps = vcat(ps, collect(eq.lhs for eq in nps.constantdefs))
        defs = copy(MT.defaults(rs))
        for eq in nps.constantdefs
            defs[eq.lhs] = eq.rhs
        end

        # add the dependent species as observed
        obs = copy(MT.observed(rs))
        append!(obs, nps.conservedeqs)
    else
        defs = MT.defaults(rs)
        obs = MT.observed(rs)
    end

    ceqs = Equation[eq for eq in get_eqs(rs) if eq isa Equation]
    if !isempty(ceqs)
        if remove_conserved
            @info """
                  Be careful mixing constraints and elimination of conservation laws.
                  Catalyst does not check that the conserved equations still hold for the
                  final coupled system of equations. Consider using `remove_conserved =
                  false` and instead calling ModelingToolkit.structural_simplify to simplify
                  any generated ODESystem or NonlinearSystem.
                  """
        end
        append!(eqs, ceqs)
    end

    eqs, sts, ps, obs, defs
end

# used by flattened systems that don't support constraint equations currently
function error_if_constraints(::Type{T}, sys::ReactionSystem) where {T <: MT.AbstractSystem}
    any(eq -> eq isa Equation, get_eqs(sys)) &&
        error("Can not convert to a system of type ", T,
              " when there are constraint equations.")
    nothing
end

# used by flattened systems that don't support differential equation constraint eqs
function error_if_constraint_odes(::Type{T},
                                  rs::ReactionSystem) where {T <: MT.AbstractSystem}
    any(eq -> (eq isa Equation) && MT.isdiffeq(eq), get_eqs(rs)) &&
        error("Cannot convert to system type $T when there are ODE constraint equations.")
    nothing
end

function spatial_convert_err(rs::ReactionSystem, systype)
    isspatial(rs) && error("Conversion to $systype is not supported for spatial networks.")
end

"""
```julia
Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
```
Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.ODESystem`.

Keyword args and default values:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate law,
  i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
  `combinatoric_ratelaws=false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
  ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
  itself defaults to true).
- `remove_conserved=false`, if set to `true` will calculate conservation laws of the
  underlying set of reactions (ignoring constraint equations), and then apply them to reduce
  the number of equations.
"""
function Base.convert(::Type{<:ODESystem}, rs::ReactionSystem; name = nameof(rs),
                      combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      include_zero_odes = true, remove_conserved = false, checks = false,
                      default_u0 = Dict(), default_p = Dict(), defaults = _merge(Dict(default_u0), Dict(default_p)),
                      kwargs...)
    spatial_convert_err(rs::ReactionSystem, ODESystem)
    fullrs = Catalyst.flatten(rs)
    remove_conserved && conservationlaws(fullrs)
    ists, ispcs = get_indep_sts(fullrs, remove_conserved)
    eqs = assemble_drift(fullrs, ispcs; combinatoric_ratelaws, remove_conserved,
                         include_zero_odes)
    eqs, sts, ps, obs, defs = addconstraints!(eqs, fullrs, ists, ispcs; remove_conserved)
    
    # Converts expressions like mm(X,v,K) to v*X/(X+K).
    expand_functions  && (eqs = [eq.lhs ~ expand_registered_functions!(eq.rhs) for eq in eqs])

    ODESystem(eqs, get_iv(fullrs), sts, ps;
              observed = obs,
              name,
              defaults = _merge(defaults,defs),
              checks,
              continuous_events = MT.get_continuous_events(fullrs),
              discrete_events = MT.get_discrete_events(fullrs),
              kwargs...)
end

"""
```julia
Base.convert(::Type{<:NonlinearSystem},rs::ReactionSystem)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.NonlinearSystem`.

Keyword args and default values:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate law,
  i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
  `combinatoric_ratelaws=false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
  ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
  itself defaults to true).
- `remove_conserved=false`, if set to `true` will calculate conservation laws of the
  underlying set of reactions (ignoring constraint equations), and then apply them to reduce
  the number of equations.
"""
function Base.convert(::Type{<:NonlinearSystem}, rs::ReactionSystem; name = nameof(rs),
                      combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      include_zero_odes = true, remove_conserved = false, checks = false,
                      kwargs...)
    spatial_convert_err(rs::ReactionSystem, NonlinearSystem)
    fullrs = Catalyst.flatten(rs)
    remove_conserved && conservationlaws(fullrs)
    ists, ispcs = get_indep_sts(fullrs, remove_conserved)
    eqs = assemble_drift(fullrs, ispcs; combinatoric_ratelaws, remove_conserved,
                         as_odes = false, include_zero_odes)
    error_if_constraint_odes(NonlinearSystem, fullrs)
    eqs, sts, ps, obs, defs = addconstraints!(eqs, fullrs, ists, ispcs; remove_conserved)

    NonlinearSystem(eqs, sts, ps;
                    name,
                    observed = obs,
                    defaults = _merge(defaults,defs),
                    checks,
                    kwargs...)
end

"""
```julia
Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.SDESystem`.

Notes:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate law,
  i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. Set
  `combinatoric_ratelaws=false` for a ratelaw of `k*S^2`, i.e. the scaling factor is
  ignored. Defaults to the value given when the `ReactionSystem` was constructed (which
  itself defaults to true).
- `noise_scaling=nothing::Union{Vector{Num},Num,Nothing}` allows for linear scaling of the
  noise in the chemical Langevin equations. If `nothing` is given, the default value as in
  Gillespie 2000 is used. Alternatively, a `Num` can be given, this is added as a parameter
  to the system (at the end of the parameter array). All noise terms are linearly scaled
  with this value. The parameter may be one already declared in the `ReactionSystem`.
  Finally, a `Vector{Num}` can be provided (the length must be equal to the number of
  reactions). Here the noise for each reaction is scaled by the corresponding parameter in
  the input vector. This input may contain repeat parameters.
- `remove_conserved=false`, if set to `true` will calculate conservation laws of the
  underlying set of reactions (ignoring constraint equations), and then apply them to reduce
  the number of equations.
- Does not currently support `ReactionSystem`s that include coupled algebraic or
  differential equations.
"""
function Base.convert(::Type{<:SDESystem}, rs::ReactionSystem;
                      noise_scaling = nothing, name = nameof(rs),
                      combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      include_zero_odes = true, checks = false, remove_conserved = false,
                      kwargs...)
    spatial_convert_err(rs::ReactionSystem, SDESystem)

    flatrs = Catalyst.flatten(rs)
    error_if_constraints(SDESystem, flatrs)

    if noise_scaling isa AbstractArray
        (length(noise_scaling) != numreactions(flatrs)) &&
            error("The number of elements in 'noise_scaling' must be equal " *
                  "to the number of reactions in the flattened reaction system.")
        if !(noise_scaling isa Symbolics.Arr)
            noise_scaling = value.(noise_scaling)
        end
    elseif !isnothing(noise_scaling)
        noise_scaling = fill(value(noise_scaling), numreactions(flatrs))
    end

    remove_conserved && conservationlaws(flatrs)
    ists, ispcs = get_indep_sts(flatrs, remove_conserved)
    eqs = assemble_drift(flatrs, ispcs; combinatoric_ratelaws, include_zero_odes,
                         remove_conserved)
    noiseeqs = assemble_diffusion(flatrs, ists, ispcs, noise_scaling; combinatoric_ratelaws,
                                  remove_conserved)
    eqs, sts, ps, obs, defs = addconstraints!(eqs, flatrs, ists, ispcs; remove_conserved)
    ps = (noise_scaling === nothing) ? ps : vcat(ps, toparam(noise_scaling))

    # Converts expressions like mm(X,v,K) to v*X/(X+K).
    if expand_functions
        eqs = [eq.lhs ~ expand_registered_functions!(eq.rhs) for eq in eqs]
        noiseeqs = [expand_registered_functions!(neq) for neq in noiseeqs]
    end

    if any(isbc, get_states(flatrs))
        @info "Boundary condition species detected. As constraint equations are not currently supported when converting to SDESystems, the resulting system will be undetermined. Consider using constant species instead."
    end

    SDESystem(eqs, noiseeqs, get_iv(flatrs), sts, ps;
              observed = obs,
              name,
              defaults = defs,
              checks,
              continuous_events = MT.get_continuous_events(flatrs),
              discrete_events = MT.get_discrete_events(flatrs),
              kwargs...)
end

"""
```julia
Base.convert(::Type{<:JumpSystem},rs::ReactionSystem; combinatoric_ratelaws=true)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.JumpSystem`.

Notes:
- `combinatoric_ratelaws=true` uses binomials in calculating the rate law, i.e. for `2S ->
  0` at rate `k` the ratelaw would be `k*S*(S-1)/2`. If `combinatoric_ratelaws=false` then
  the ratelaw is `k*S*(S-1)`, i.e. the rate law is not normalized by the scaling factor.
  Defaults to the value given when the `ReactionSystem` was constructed (which itself
  defaults to true).
- Does not currently support `ReactionSystem`s that include coupled algebraic or
  differential equations.
- Does not currently support continuous events as these are not supported by
  `ModelingToolkit.JumpSystems`.
"""
function Base.convert(::Type{<:JumpSystem}, rs::ReactionSystem; name = nameof(rs),
                      combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                      remove_conserved = nothing, checks = false,
                      default_u0 = Dict(), default_p = Dict(), defaults = _merge(Dict(default_u0), Dict(default_p)),
                      kwargs...)
    spatial_convert_err(rs::ReactionSystem, JumpSystem)

    (remove_conserved !== nothing) &&
        error("Catalyst does not support removing conserved species when converting to JumpSystems.")

    flatrs = Catalyst.flatten(rs)
    error_if_constraints(JumpSystem, flatrs)

    (length(MT.continuous_events(flatrs)) > 0) &&
        (@warn "continuous_events will be dropped as they are not currently supported by JumpSystems.")

    eqs = assemble_jumps(flatrs; combinatoric_ratelaws)

    # handle BC species
    sts, ispcs = get_indep_sts(flatrs)
    any(isbc, get_states(flatrs)) && (sts = vcat(sts, filter(isbc, get_states(flatrs))))
    ps = get_ps(flatrs)

    JumpSystem(eqs, get_iv(flatrs), sts, ps;
               observed = MT.observed(flatrs),
               name,
               defaults = _merge(defaults,MT.defaults(flatrs)),
               checks,
               discrete_events = MT.discrete_events(flatrs),
               kwargs...)
end

### Converts a reaction system to ODE or SDE problems ###

# ODEProblem from AbstractReactionNetwork
function DiffEqBase.ODEProblem(rs::ReactionSystem, u0, tspan,
                               p = DiffEqBase.NullParameters(), args...;
                               check_length = false, name = nameof(rs),
                               combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                               include_zero_odes = true, remove_conserved = false,
                               checks = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    osys = convert(ODESystem, rs; name, combinatoric_ratelaws, include_zero_odes, checks,
                   remove_conserved)
    return ODEProblem(osys, u0map, tspan, pmap, args...; check_length, kwargs...)
end

# NonlinearProblem from AbstractReactionNetwork
function DiffEqBase.NonlinearProblem(rs::ReactionSystem, u0,
                                     p = DiffEqBase.NullParameters(), args...;
                                     name = nameof(rs), include_zero_odes = true,
                                     combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                                     remove_conserved = false, checks = false,
                                     check_length = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    nlsys = convert(NonlinearSystem, rs; name, combinatoric_ratelaws, include_zero_odes,
                    checks, remove_conserved)
    return NonlinearProblem(nlsys, u0map, pmap, args...; check_length, kwargs...)
end

# SDEProblem from AbstractReactionNetwork
function DiffEqBase.SDEProblem(rs::ReactionSystem, u0, tspan,
                               p = DiffEqBase.NullParameters(), args...;
                               noise_scaling = nothing, name = nameof(rs),
                               combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                               include_zero_odes = true, checks = false,
                               check_length = false,
                               remove_conserved = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    sde_sys = convert(SDESystem, rs; noise_scaling, name, combinatoric_ratelaws,
                      include_zero_odes, checks, remove_conserved)
    p_matrix = zeros(length(get_states(sde_sys)), numreactions(rs))
    return SDEProblem(sde_sys, u0map, tspan, pmap, args...; check_length,
                      noise_rate_prototype = p_matrix, kwargs...)
end

# DiscreteProblem from AbstractReactionNetwork
function DiffEqBase.DiscreteProblem(rs::ReactionSystem, u0, tspan::Tuple,
                                    p = DiffEqBase.NullParameters(), args...;
                                    name = nameof(rs),
                                    combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                                    checks = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    jsys = convert(JumpSystem, rs; name, combinatoric_ratelaws, checks)
    return DiscreteProblem(jsys, u0map, tspan, pmap, args...; kwargs...)
end

# JumpProblem from AbstractReactionNetwork
function JumpProcesses.JumpProblem(rs::ReactionSystem, prob, aggregator, args...;
                                   name = nameof(rs),
                                   combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                                   checks = false, kwargs...)
    jsys = convert(JumpSystem, rs; name, combinatoric_ratelaws, checks)
    return JumpProblem(jsys, prob, aggregator, args...; kwargs...)
end

# SteadyStateProblem from AbstractReactionNetwork
function DiffEqBase.SteadyStateProblem(rs::ReactionSystem, u0,
                                       p = DiffEqBase.NullParameters(), args...;
                                       check_length = false, name = nameof(rs),
                                       combinatoric_ratelaws = get_combinatoric_ratelaws(rs),
                                       remove_conserved = false, include_zero_odes = true,
                                       checks = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    osys = convert(ODESystem, rs; name, combinatoric_ratelaws, include_zero_odes, checks,
                   remove_conserved)
    return SteadyStateProblem(osys, u0map, pmap, args...; check_length, kwargs...)
end

####################### dependency graph utilities ########################

# determine which states a reaction depends on
function ModelingToolkit.get_variables!(deps::Set, rx::Reaction, variables)
    (rx.rate isa Symbolic) && get_variables!(deps, rx.rate, variables)
    for s in rx.substrates
        # parametric stoichiometry means may have a parameter as a substrate
        any(isequal(s), variables) && push!(deps, s)
    end
    deps
end

# determine which species a reaction modifies
function ModelingToolkit.modified_states!(mstates, rx::Reaction, sts::Set)
    for (species, stoich) in rx.netstoich
        (species in sts) && push!(mstates, species)
    end
    mstates
end

function ModelingToolkit.modified_states!(mstates, rx::Reaction, sts::AbstractVector)
    for (species, stoich) in rx.netstoich
        any(isequal(species), sts) && push!(mstates, species)
    end
    mstates
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

    ReactionSystem(equations(rs), get_iv(rs), states(rs), parameters(rs);
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
    sts = union(get_states(rs), get_states(sys))
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
