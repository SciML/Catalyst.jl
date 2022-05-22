
"""
$(TYPEDEF)

One chemical reaction.

# Fields
$(FIELDS)

# Examples

```julia
using Catalyst
@parameters k[1:20]
@variables t A(t) B(t) C(t) D(t)
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
    rate
    """Reaction substrates."""
    substrates::Vector
    """Reaction products."""
    products::Vector
    """The stoichiometric coefficients of the reactants."""
    substoich::Vector{T}
    """The stoichiometric coefficients of the products."""
    prodstoich::Vector{T}
    """The net stoichiometric coefficients of all species changed by the reaction."""
    netstoich::Vector{Pair{S,T}}
    """
    `false` (default) if `rate` should be multiplied by mass action terms to give the rate law.
    `true` if `rate` represents the full reaction rate law.
    """
    only_use_rate::Bool
end

function Reaction(rate, subs, prods, substoich, prodstoich;
                  netstoich=nothing, only_use_rate=false,
                  kwargs...)

    (isnothing(prods)&&isnothing(subs)) && error("A reaction requires a non-nothing substrate or product vector.")
    (isnothing(prodstoich)&&isnothing(substoich)) && error("Both substrate and product stochiometry inputs cannot be nothing.")
    if isnothing(subs)
        subs = Vector{Term}()
        !isnothing(substoich) && error("If substrates are nothing, substrate stiocihometries have to be so too.")
        substoich = typeof(prodstoich)()
    end
    S = eltype(substoich)

    if isnothing(prods)
        prods = Vector{Term}()
        !isnothing(prodstoich) && error("If products are nothing, product stiocihometries have to be so too.")
        prodstoich = typeof(substoich)()
    end
    T = eltype(prodstoich)

    subs = value.(subs)
    prods = value.(prods)

    # try to get a common type for stoichiometry, using Any if have Syms
    stoich_type = promote_type(S,T)
    if stoich_type <: Num
        stoich_type = Any
        substoich′  = Any[value(s) for s in substoich]
        prodstoich′ = Any[value(p) for p in prodstoich]
    else
        substoich′ = (S == stoich_type) ? substoich : convert.(stoich_type, substoich)
        prodstoich′ = (T == stoich_type) ? prodstoich : convert.(stoich_type, prodstoich)
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
    sstoich = isnothing(subs) ? nothing : ones(Int,length(subs))
    pstoich = isnothing(prods) ? nothing : ones(Int,length(prods))
    Reaction(rate, subs, prods, sstoich, pstoich; kwargs...)
end

function print_rxside(io::IO, specs, stoich)
    # reactants/substrates
    if isempty(specs)
        print(io, "∅")
    else
        for (i,spec) in enumerate(specs)
            if isequal(stoich[i],one(stoich[i]))
                print(io, ModelingToolkit.operation(spec))
            elseif Symbolics.istree(stoich[i])
                print(io, "(", stoich[i], ")*", ModelingToolkit.operation(spec))
            else
                print(io, stoich[i], "*", ModelingToolkit.operation(spec))
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

function ModelingToolkit.namespace_equation(rx::Reaction, name)
    subs  = isempty(rx.substrates) ? rx.substrates : [namespace_expr(sub, name) for sub in rx.substrates]
    prods = isempty(rx.products) ? rx.products : [namespace_expr(prod, name) for prod in rx.products]
    Reaction(namespace_expr(rx.rate, name),
             subs, prods, rx.substoich, rx.prodstoich,
             [namespace_expr(n[1],name) => n[2] for n in rx.netstoich], rx.only_use_rate)
end

netstoich_stoichtype(::Vector{Pair{S,T}}) where {S,T} = T

# calculates the net stoichiometry of a reaction as a vector of pairs (sub,substoich)
function get_netstoich(subs, prods, sstoich, pstoich)
    # stoichiometry as a Dictionary
    nsdict = Dict{Any, eltype(sstoich)}(sub => -sstoich[i] for (i,sub) in enumerate(subs))
    for (i,p) in enumerate(prods)
        coef = pstoich[i]
        @inbounds nsdict[p] = haskey(nsdict, p) ? nsdict[p] + coef : coef
    end

    # stoichiometry as a vector
    [el for el in nsdict if !_iszero(el[2])]
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
struct ReactionComplex{V<:Integer} <: AbstractVector{ReactionComplexElement{V}}
    """The integer ids of all species participating in this complex."""
    speciesids::Vector{Int}
    """The stoichiometric coefficients of all species participating in this complex."""
    speciesstoichs::Vector{V}
end

function (==)(a::ReactionComplex{V},b::ReactionComplex{V}) where {V <: Integer}
    (a.speciesids == b.speciesids) &&
    (a.speciesstoichs == b.speciesstoichs)
end
hash(rc::ReactionComplex,h::UInt) = Base.hash(rc.speciesids,Base.hash(rc.speciesstoichs,h))
Base.size(rc::ReactionComplex) = size(rc.speciesids)
Base.length(rc::ReactionComplex) = length(rc.speciesids)
Base.getindex(rc::ReactionComplex, i...) =
        ReactionComplexElement(getindex(rc.speciesids, i...), getindex(rc.speciesstoichs, i...))
Base.setindex!(rc::ReactionComplex, t::ReactionComplexElement, i...) =
    (setindex!(rc.speciesids, t.speciesid, i...); setindex!(rc.speciesstoichs, t.speciesstoich, i...); rc)
Base.isless(a::ReactionComplexElement, b::ReactionComplexElement) = isless(a.speciesid, b.speciesid)
Base.Sort.defalg(::ReactionComplex{T}) where {T <: Integer} = Base.DEFAULT_UNSTABLE

############################### Network Properties ####################################

# Internal cache for various ReactionSystem calculated properties
Base.@kwdef mutable struct NetworkProperties{I <: Integer, V <: Term}
    isempty::Bool = true
    netstoichmat::Union{Matrix{Int},SparseMatrixCSC{Int,Int}} = Matrix{Int}(undef,0,0)
    conservationmat::Matrix{I} = Matrix{I}(undef,0,0)
    col_order::Vector{Int} = Int[]
    rank::Int = 0
    nullity::Int = 0
    indepspecs::Set{V} = Set{V}()
    depspecs::Set{V} = Set{V}()
    conservedeqs::Vector{Equation} = Equation[]
    constantdefs::Vector{Equation} = Equation[]
    speciesmap::Dict{V,Int} = Dict{V,Int}()
    complextorxsmap::OrderedDict{ReactionComplex{V},Vector{Pair{Int,Int}}} = OrderedDict{ReactionComplex{V},Vector{Pair{Int,Int}}}()
    complexes::Vector{ReactionComplex{V}} = Vector{ReactionComplex{V}}(undef,0)
    incidencemat::Union{Matrix{Int},SparseMatrixCSC{Int,Int}} = Matrix{Int}(undef,0.0)
    complexstoichmat::Union{Matrix{Int},SparseMatrixCSC{Int,Int}} = Matrix{Int}(undef,0.0)
    complexoutgoingmat::Union{Matrix{Int},SparseMatrixCSC{Int,Int}} = Matrix{Int}(undef,0.0)
    incidencegraph::Graphs.SimpleDiGraph{Int} = Graphs.DiGraph()
    linkageclasses::Vector{Vector{Int}} = Vector{Vector{Int}}(undef,0)
    deficiency::Int = 0
    subnetworks::Vector{ReactionSystem} = ReactionSystem[]
end

function Base.show(io::IO, nps::NetworkProperties)
    if (nps.conservationmat !== nothing)
        println(io, "Conserved Equations: ")
        foreach(eq -> println(io,eq), nps.conservedeqs)
        println()
    end
end

Base.isempty(nps::NetworkProperties) = getfield(nps, :isempty)

function Base.setproperty!(nps::NetworkProperties, sym::Symbol, x)
    (sym !== :isempty) && setfield!(nps, :isempty, false)
    setfield!(nps, sym, x)
end

function reset!(nps::NetworkProperties{I,V}) where {I,V}
    nps.isempty && return
    nps.netstoichmat = Matrix{I}(undef,0,0)
    nps.conservationmat = Matrix{I}(undef,0,0)
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
    empty!(nps.incidencemat)
    empty!(nps.complexstoichmat)
    empty!(nps.complexoutgoingmat)
    empty!(nps.incidencegraph)
    empty!(linkageclasses)
    nps.deficiency = 0
    empty!(nps.subnetworks)

    # this needs to be last due to setproperty! setting it to false
    nps.isempty = true
    nothing
end

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
- `systems::Vector{AbstractSystems}`, vector of sub-systems. Can be
  `ReactionSystem`s, `ODESystem`s, or `NonlinearSystem`s.
- `name::Symbol`, the name of the system (must be provided, or `@named` must be
  used).
- `defaults::Dict`, a dictionary mapping parameters to their default values and
  species to their default initial values.
- `checks = true`, boolean for whether to check units.
- `constraints = nothing`, a `NonlinearSystem` or `ODESystem` of coupled
  constraint equations.
- `networkproperties = NetworkProperties()`, cache for network properties calculated via API functions.

Notes:
- ReactionSystems currently do rudimentary unit checking, requiring that all
  species have the same units, and all reactions have rate laws with units of
  (species units) / (time units). Unit checking can be disabled by passing the
  keyword argument `checks=false`.
"""
struct ReactionSystem{U <: Union{Nothing,MT.AbstractSystem}, V <: NetworkProperties} <: MT.AbstractTimeDependentSystem
    """The reactions defining the system."""
    eqs::Vector{Reaction}
    """Independent variable (usually time)."""
    iv::Any
    """Dependent (state) variables representing amount of each species. Must not contain the
    independent variable."""
    states::Vector
    """Parameter variables. Must not contain the independent variable."""
    ps::Vector
    """Maps Symbol to corresponding variable."""
    var_to_name::Dict{Symbol,Any}
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
    """Non-`Reaction` equations that further constrain the system"""
    constraints::U
    """`NetworkProperties` object that can be filled in by API functions. INTERNAL -- not
    considered part of the public API."""
    networkproperties::V

    # inner constructor is considered private and may change between non-breaking releases.
    function ReactionSystem(eqs, iv, states, ps, var_to_name, observed, name, systems,
                            defaults, connection_type, csys, nps;
                            checks::Bool=true)
        if checks
            check_variables(states, iv)
            check_parameters(ps, iv)
        end

        rs = new{typeof(csys), typeof(nps)}(eqs, iv, states, ps, var_to_name, observed, name,
                                            systems, defaults, connection_type, csys, nps)
        checks && validate(rs)
        rs
    end
end

function get_speciestype(iv, states, systems, constraints)
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

function ReactionSystem(eqs, iv, states, ps;
                        observed = [],
                        systems = [],
                        name = nothing,
                        default_u0=Dict(),
                        default_p=Dict(),
                        defaults=_merge(Dict(default_u0), Dict(default_p)),
                        connection_type=nothing,
                        checks = true,
                        constraints = nothing,
                        networkproperties = nothing)
    name === nothing && throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    sysnames = nameof.(systems)
    (length(unique(sysnames)) == length(sysnames)) ||
        throw(ArgumentError("System names must be unique."))

    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn("`default_u0` and `default_p` are deprecated. Use `defaults` instead.", :ReactionSystem, force=true)
    end
    defaults = MT.todict(defaults)
    defaults = Dict{Any,Any}(value(k) => value(v) for (k, v) in pairs(defaults))

    iv′     = value(iv)
    states′ = value.(MT.scalarize(states))
    ps′     = value.(MT.scalarize(ps))
    eqs′    = (eqs isa Vector) ? eqs : collect(eqs)

    var_to_name = Dict()
    MT.process_variables!(var_to_name, defaults, states′)
    MT.process_variables!(var_to_name, defaults, ps′)
    MT.collect_var_to_name!(var_to_name, eq.lhs for eq in observed)

    nps = if networkproperties === nothing
        NetworkProperties{Int, get_speciestype(iv′, states′, systems, constraints)}()
    else
        networkproperties
    end

    ReactionSystem(eqs′, iv′, states′, ps′, var_to_name, observed, name, systems,
                   defaults, connection_type, constraints, nps; checks = checks)
end


# Previous function called by the macro, but still avaiable for general use.
function ReactionSystem(rxs::Vector{<:Reaction}, iv; kwargs...)
    make_ReactionSystem_internal(rxs, iv, nothing, Vector{Num}(); kwargs...)
end

# search the symbolic expression for parameters or states
# and save in ps and sts respectively. vars is used to cache results
function findvars!(ps, sts, exprtosearch, t, vars)
    MT.get_variables!(vars, exprtosearch)
    for var in vars
        isequal(t,var) && continue
        if MT.isparameter(var)
            push!(ps, var)
        else
            push!(sts, var)
        end
    end
    empty!(vars)
end

# Only used internally by the @reaction_network macro. Permits giving an initial order to the parameters,
# and then adds additional ones found in the reaction. Name could be changed.
function make_ReactionSystem_internal(rxs::Vector{<:Reaction}, iv, no_sps::Nothing, ps_in; kwargs...)
    t    = value(iv)
    sts  = OrderedSet(spec for rx in rxs for spec in Iterators.flatten((rx.substrates,rx.products)))
    ps   = OrderedSet{Any}(ps_in)
    vars = OrderedSet()
    for rx in rxs
        findvars!(ps, sts, rx.rate, t, vars)
        for s in rx.substoich
            (s isa Symbolics.Symbolic) && findvars!(ps, sts, s, t, vars)
        end
        for p in rx.prodstoich
            (p isa Symbolics.Symbolic) && findvars!(ps, sts, p, t, vars)
        end
    end
    ReactionSystem(rxs, t, collect(sts), collect(ps); kwargs...)
end


function ReactionSystem(iv; kwargs...)
    ReactionSystem(Reaction[], iv, [], []; kwargs...)
end

####################### ModelingToolkit inherited accessors #############################

"""
    get_constraints(sys::ReactionSystem)

Return the current constraint subsystem, if none is defined will return `nothing`.
"""
get_constraints(sys::ReactionSystem) = getfield(sys, :constraints)
has_constraints(sys::ReactionSystem) = isdefined(sys, :constraints)

"""
    get_networkproperties(sys::ReactionSystem)

Return the current network properties of `sys`.
"""
get_networkproperties(sys::ReactionSystem) = getfield(sys, :networkproperties)

function MT.states(sys::ReactionSystem)
    sts = (get_constraints(sys) === nothing) ? get_states(sys) : vcat(get_states(sys), get_states(get_constraints(sys)))
    systems = get_systems(sys)
    unique(isempty(systems) ? sts : [sts; reduce(vcat,namespace_variables.(systems))])
end

function MT.parameters(sys::ReactionSystem)
    ps = (get_constraints(sys) === nothing) ? get_ps(sys) : vcat(get_ps(sys), get_ps(get_constraints(sys)))
    systems = get_systems(sys)
    unique(isempty(systems) ? ps : [ps; reduce(vcat,namespace_parameters.(systems))])
end

function MT.equations(sys::ReactionSystem)
    eqs = (get_constraints(sys) === nothing) ? get_eqs(sys) : Any[get_eqs(sys); get_eqs(get_constraints(sys))]
    systems = get_systems(sys)
    if !isempty(systems)
        return Any[eqs; reduce(vcat, MT.namespace_equations.(systems); init=Any[])]
    end
    return eqs
end

function MT.defaults(sys::ReactionSystem)
    constraints = get_constraints(sys)
    systems = get_systems(sys)
    defs = get_defaults(sys)

    if constraints !== nothing
        defs = merge(defs, get_defaults(constraints))
    end

    mapfoldr(MT.namespace_defaults, merge, systems, init=defs)
end

function MT.observed(sys::ReactionSystem)
    constraints = get_constraints(sys)
    systems = get_systems(sys)
    obs = get_observed(sys)

    if constraints !== nothing
        obs = vcat(obs, get_observed(constraints))
    end
    reduce(vcat,
           (map(o->MT.namespace_equation(o,s), MT.observed(s)) for s in systems);
           init=obs)
end

function MT.getvar(sys::ReactionSystem, name::Symbol; namespace=false)
    if isdefined(sys, name)
        Base.depwarn("`sys.name` like `sys.$name` is deprecated. Use getters like `get_$name` instead.", "sys.$name")
        return getfield(sys, name)
    end

    systems = get_systems(sys)
    i = findfirst(x->nameof(x)==name, systems)
    i === nothing || return namespace ? rename(systems[i], renamespace(sys, name)) : systems[i]

    constraints = get_constraints(sys)
    if constraints !== nothing && nameof(constraints) == name
        return namespace ? rename(constraints, renamespace(sys, name)) : constraints
    end

    avs = get_var_to_name(sys)
    v = get(avs, name, nothing)
    v === nothing || return namespace ? renamespace(sys, v) : v

    sts = get_states(sys)
    i = findfirst(x->getname(x) == name, sts)
    i === nothing || return namespace ? renamespace(sys, sts[i]) : sts[i]

    ps = get_ps(sys)
    i = findfirst(x->getname(x) == name,ps)
    i === nothing || return namespace ? renamespace(sys, ps[i]) : ps[i]

    obs = get_observed(sys)
    i = findfirst(x->getname(x.lhs)==name,obs)
    i === nothing || return namespace ? renamespace(sys, obs[i]) : obs[i]

    if constraints !== nothing
        return getvar(rename(constraints, nameof(sys)), name; namespace=namespace)
    end

    throw(ArgumentError("System $(nameof(sys)): variable $name does not exist"))
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
function oderatelaw(rx; combinatoric_ratelaw=true)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate

    # if the stoichiometric coefficients are not integers error if asking to scale rates
    !all(s -> s isa Union{Integer,Symbolics.Symbolic}, substoich) && (combinatoric_ratelaw==true) &&
        error("Non-integer stoichiometric coefficients require the combinatoric_ratelaw=false keyword to oderatelaw, or passing combinatoric_ratelaws=false to convert or ODEProblem.")

    if !only_use_rate
        coef = eltype(substoich) <: Number ? one(eltype(substoich)) : 1
        for (i,stoich) in enumerate(substoich)
            combinatoric_ratelaw && (coef *= factorial(stoich))
            rl *= isequal(stoich,one(stoich)) ? substrates[i] : substrates[i]^stoich
        end
        combinatoric_ratelaw && (!isequal(coef,one(coef))) && (rl /= coef)
    end
    rl
end

function assemble_oderhs(rs; combinatoric_ratelaws=true, remove_conserved=false)
    sts = get_states(rs)
    nps = get_networkproperties(rs)
    if remove_conserved
        # get the independent species in an order consistent with sts
        indepsts = Iterators.filter(in(nps.indepspecs), sts)
        species_to_idx = Dict(x => i for (i,x) in enumerate(indepsts))
        rhsvec = Any[0 for i in eachindex(species_to_idx)]
        depspec_submap = Dict(eq.lhs => eq.rhs for eq in nps.conservedeqs)
    else
        species_to_idx = Dict((x => i for (i,x) in enumerate(sts)))
        rhsvec = Any[0 for i in eachindex(sts)]
        depspec_submap = Dict()
    end

    for rx in get_eqs(rs)
        rl = oderatelaw(rx; combinatoric_ratelaw=combinatoric_ratelaws)
        remove_conserved && (rl = substitute(rl, depspec_submap))
        for (spec,stoich) in rx.netstoich
            # dependent species don't get an ODE, so are skipped
            remove_conserved && (spec in nps.depspecs) && continue

            i = species_to_idx[spec]
            if _iszero(rhsvec[i])
                if stoich isa Symbolics.Symbolic
                    rhsvec[i] = stoich * rl
                else
                    signedrl  = (stoich > zero(stoich)) ? rl : -rl
                    rhsvec[i] = isone(abs(stoich)) ? signedrl : stoich * rl
                end
            else
                if stoich isa Symbolics.Symbolic
                    rhsvec[i] += stoich * rl
                else
                    Δspec     = isone(abs(stoich)) ? rl : abs(stoich) * rl
                    rhsvec[i] = (stoich > zero(stoich)) ? (rhsvec[i] + Δspec) : (rhsvec[i] - Δspec)
                end
            end
        end
    end

    rhsvec
end

function assemble_drift(rs; combinatoric_ratelaws=true, as_odes=true, include_zero_odes=true,
                            remove_conserved=false)
    rhsvec = assemble_oderhs(rs; combinatoric_ratelaws, remove_conserved)
    if as_odes
        D = Differential(get_iv(rs))
        if remove_conserved
            indepspecs = get_networkproperties(rs).indepspecs
            sts = Iterators.filter(in(indepspecs), states(rs))
        else
            sts = states(rs)
        end
        eqs = [Equation(D(x),rhs) for (x,rhs) in zip(sts,rhsvec) if (include_zero_odes || (!_iszero(rhs)))]
    else
        eqs = [Equation(0,rhs) for rhs in rhsvec if (include_zero_odes || (!_iszero(rhs)))]
    end
    eqs
end

function assemble_diffusion(rs, noise_scaling; combinatoric_ratelaws=true)
    sts  = get_states(rs)
    eqs  = Matrix{Any}(undef, length(sts), length(get_eqs(rs)))
    eqs .= 0
    species_to_idx = Dict((x => i for (i,x) in enumerate(sts)))

    for (j,rx) in enumerate(get_eqs(rs))
        rlsqrt = sqrt(abs(oderatelaw(rx; combinatoric_ratelaw=combinatoric_ratelaws)))
        (noise_scaling!==nothing) && (rlsqrt *= noise_scaling[j])
        for (spec,stoich) in rx.netstoich
            i = species_to_idx[spec]
            if stoich isa Symbolics.Symbolic
                eqs[i,j] = stoich * rlsqrt
            else
                signedrlsqrt = (stoich > zero(stoich)) ? rlsqrt : -rlsqrt
                eqs[i,j]     = isone(abs(stoich)) ? signedrlsqrt : stoich * rlsqrt
            end
        end
    end
    eqs
end

"""
    jumpratelaw(rx; rxvars=get_variables(rx.rate), combinatoric_ratelaw=true)

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
- `rxvars` should give the `Variable`s, i.e. species and parameters, the rate depends on.
- Allocates
- `combinatoric_ratelaw=true` uses binomials in calculating the rate law, i.e. for `2S ->
  0` at rate `k` the ratelaw would be `k*S*(S-1)/2`. If `combinatoric_ratelaw=false` then
  the ratelaw is `k*S*(S-1)`, i.e. the rate law is not normalized by the scaling
  factor.
"""
function jumpratelaw(rx; rxvars=get_variables(rx.rate), combinatoric_ratelaw=true)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate
    if !only_use_rate
        coef = eltype(substoich) <: Number ? one(eltype(substoich)) : 1
        for (i,stoich) in enumerate(substoich)
            s = substrates[i]
            if stoich isa Symbolics.Symbolic
                rl *= combinatoric_ratelaw ? binomial(s,stoich) : factorial(s) / factorial(s-stoich)
            else
                rl *= s
                isone(stoich) && continue
                for i in one(stoich):(stoich-one(stoich))
                    rl *= (s - i)
                end
                combinatoric_ratelaw && (coef *= factorial(stoich))
            end
        end
        combinatoric_ratelaw && !isequal(coef,one(coef)) && (rl /= coef)
    end
    rl
end

# if haveivdep=false then time dependent rates will still be classified as mass action
"""
```julia
ismassaction(rx, rs; rxvars = get_variables(rx.rate),
                              haveivdep = any(var -> isequal(get_iv(rs),var), rxvars),
                              stateset = Set(states(rs)))
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
  explicitly depends on the independent variable.
- Optional: `stateset`, set of states which if the rxvars are within mean rx is
  non-mass action.

Notes:
- Non-integer stoichiometry is treated as non-mass action. This includes
  symbolic variables/terms or floating point numbers for stoichiometric
  coefficients.
"""
function ismassaction(rx, rs; rxvars = get_variables(rx.rate),
                              haveivdep = any(var -> isequal(get_iv(rs),var), rxvars),
                              stateset = Set(get_states(rs)))

    # we define non-integer (i.e. float or symbolic) stoich to be non-mass action
    ((eltype(rx.substoich) <: Integer) && (eltype(rx.prodstoich) <: Integer)) || return false

    # if no dependencies must be zero order
    (length(rxvars)==0) && return true
    haveivdep && return false
    rx.only_use_rate && return false
    @inbounds for i = 1:length(rxvars)
        (rxvars[i] in stateset) && return false
    end
    return true
end

@inline function makemajump(rx; combinatoric_ratelaw=true)
    @unpack rate, substrates, substoich, netstoich = rx
    zeroorder = (length(substoich) == 0)
    reactant_stoch = Vector{Pair{Any,eltype(substoich)}}(undef, length(substoich))
    @inbounds for i = 1:length(reactant_stoch)
        reactant_stoch[i] = substrates[i] => substoich[i]
    end
    #push!(rstoich, reactant_stoch)
    coef = (zeroorder || (!combinatoric_ratelaw)) ? one(eltype(substoich)) : prod(stoich -> factorial(stoich), substoich)
    (!isone(coef)) && (rate /= coef)
    #push!(rates, rate)
    net_stoch      = [Pair(p[1],p[2]) for p in netstoich]
    #push!(nstoich, net_stoch)
    MassActionJump(Num(rate), reactant_stoch, net_stoch, scale_rates=false, useiszero=false)
end

function assemble_jumps(rs; combinatoric_ratelaws=true)
    meqs = MassActionJump[]; ceqs = ConstantRateJump[]; veqs = VariableRateJump[]
    stateset = Set(get_states(rs))
    #rates = [];  rstoich = []; nstoich = []
    rxvars = []
    ivname = nameof(get_iv(rs))

    isempty(get_eqs(rs)) && error("Must give at least one reaction before constructing a JumpSystem.")
    for rx in get_eqs(rs)
        empty!(rxvars)
        (rx.rate isa Symbolic) && get_variables!(rxvars, rx.rate)
        haveivdep = false
        @inbounds for i = 1:length(rxvars)
            if isequal(rxvars[i], get_iv(rs))
                haveivdep = true
                break
            end
        end
        if ismassaction(rx, rs; rxvars=rxvars, haveivdep=haveivdep, stateset=stateset)
            push!(meqs, makemajump(rx, combinatoric_ratelaw=combinatoric_ratelaws))
        else
            rl     = jumpratelaw(rx, rxvars=rxvars, combinatoric_ratelaw=combinatoric_ratelaws)
            affect = Vector{Equation}()
            for (spec,stoich) in rx.netstoich
                push!(affect, spec ~ spec + stoich)
            end
            if haveivdep
                push!(veqs, VariableRateJump(rl,affect))
            else
                push!(ceqs, ConstantRateJump(rl,affect))
            end
        end
    end
    #eqs[1] = MassActionJump(rates, rstoich, nstoich, scale_rates=false, useiszero=false)
    vcat(meqs,ceqs,veqs)
end

# convert subsystems to the system type T
# function make_systems_with_type!(systems::Vector{T}, rs::ReactionSystem, include_zero_odes=true) where {T <: MT.AbstractSystem}
#     resize!(systems, length(get_systems(rs)))
#     for (i,sys) in enumerate(get_systems(rs))
#         if sys isa ReactionSystem
#             systems[i] = convert(T, sys, include_zero_odes=include_zero_odes)
#         elseif sys isa T
#             systems[i] = sys
#         else
#             try
#                 if (T <: MT.AbstractTimeDependentSystem) && (sys isa MT.AbstractTimeIndependentSystem)
#                     systems[i] = MT.convert_system(T, sys, get_iv(rs))
#                 else
#                     systems[i] = MT.convert_system(T, sys)
#                 end
#             catch e
#                 error("ModelingToolkit does not currently support convert_system($T, $(typeof(sys)))")
#             end
#         end
#     end
#     systems
# end

# merge constraint eqs, states and ps into the top-level eqs, states and ps
function addconstraints!(eqs, rs::ReactionSystem; remove_conserved=false)
    csys = get_constraints(rs)
    nps  = get_networkproperties(rs)

    # make dependent species observables and add conservation constants as parameters
    if remove_conserved
        # preserve the order of the independent species
        sts = filter(in(nps.indepspecs), get_states(rs))

        # add the constants as parameters and set their values
        ps  = Any[p for p in get_ps(rs)]  # as types can be mixed
        append!(ps, eq.lhs for eq in nps.constantdefs)
        defs = copy(MT.defaults(rs))
        for eq in nps.constantdefs
            defs[eq.lhs] = eq.rhs
        end

        # add the dependent species as observed
        obs = copy(MT.observed(rs))
        append!(obs, nps.conservedeqs)
    else
        sts = get_states(rs)
        ps = get_ps(rs)
        defs = MT.defaults(rs)
        obs = MT.observed(rs)
    end

    if csys !== nothing
        sts = vcat(sts, get_states(csys))
        ps = vcat(ps, get_ps(csys))
        ceqs = get_eqs(csys)
        (!isempty(ceqs)) && append!(eqs,ceqs)
        defs = merge(defs, MT.defaults(csys))
        obs = vcat(obs, MT.observed(csys))
    end

    eqs,sts,ps,obs,defs
end

# used by systems that don't support constraint equations currently
function error_if_constraints(::Type{T}, sys::ReactionSystem) where {T <: MT.AbstractSystem}
    (get_constraints(sys) === nothing) ||
            error("Can not convert to a system of type ", T, " when there are constraints.")
end

function error_if_constraint_odes(::Type{T}, rs::ReactionSystem) where {T <: MT.AbstractSystem}
    csys = get_constraints(rs)
    if csys !== nothing
        any(eq -> MT.isdiffeq(eq), equations(csys)) &&
            error("Cannot convert to system type $T when there are ODE constraint equations.")
    end
    nothing
end

"""
```julia
Base.convert(::Type{<:ODESystem},rs::ReactionSystem)
```
Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.ODESystem`.

Notes:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate
law, i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. If
`combinatoric_ratelaws=false` then the ratelaw is `k*S^2`, i.e. the scaling factor is
ignored.
"""
function Base.convert(::Type{<:ODESystem}, rs::ReactionSystem;
                      name=nameof(rs), combinatoric_ratelaws=true, include_zero_odes=true,
                      remove_conserved=false, checks=false, kwargs...)
    fullrs = Catalyst.flatten(rs)
    remove_conserved && conservationlaws(fullrs)
    eqs = assemble_drift(fullrs; combinatoric_ratelaws, remove_conserved, include_zero_odes)
    eqs,sts,ps,obs,defs = addconstraints!(eqs, fullrs; remove_conserved)
    ODESystem(eqs, get_iv(fullrs), sts, ps; name, defaults=defs, observed=obs,
                                            checks, kwargs...)
end

"""
```julia
Base.convert(::Type{<:NonlinearSystem},rs::ReactionSystem)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.NonlinearSystem`.

Notes:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate
law, i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. If
`combinatoric_ratelaws=false` then the ratelaw is `k*S^2`, i.e. the scaling factor is
ignored.
"""
function Base.convert(::Type{<:NonlinearSystem}, rs::ReactionSystem;
                      name=nameof(rs), combinatoric_ratelaws=true, include_zero_odes=true,
                      remove_conserved=false, checks=false, kwargs...)
    fullrs = Catalyst.flatten(rs)
    remove_conserved && conservationlaws(fullrs)
    eqs = assemble_drift(fullrs; combinatoric_ratelaws, remove_conserved, as_odes=false,
                                 include_zero_odes)
    error_if_constraint_odes(NonlinearSystem, fullrs)
    eqs,sts,ps,obs,defs = addconstraints!(eqs, fullrs; remove_conserved)
    NonlinearSystem(eqs, sts, ps; name, defaults=defs, observed=obs, checks, kwargs...)
end

"""
```julia
Base.convert(::Type{<:SDESystem},rs::ReactionSystem)
```

Convert a [`ReactionSystem`](@ref) to an `ModelingToolkit.SDESystem`.

Notes:
- `combinatoric_ratelaws=true` uses factorial scaling factors in calculating the rate
law, i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. If
`combinatoric_ratelaws=false` then the ratelaw is `k*S^2`, i.e. the scaling factor is
ignored.
- `noise_scaling=nothing::Union{Vector{Num},Num,Nothing}` allows for linear
scaling of the noise in the chemical Langevin equations. If `nothing` is given, the default
value as in Gillespie 2000 is used. Alternatively, a `Num` can be given, this is
added as a parameter to the system (at the end of the parameter array). All noise terms
are linearly scaled with this value. The parameter may be one already declared in the `ReactionSystem`.
Finally, a `Vector{Num}` can be provided (the length must be equal to the number of reactions).
Here the noise for each reaction is scaled by the corresponding parameter in the input vector.
This input may contain repeat parameters.
"""
function Base.convert(::Type{<:SDESystem}, rs::ReactionSystem;
                      noise_scaling=nothing, name=nameof(rs), combinatoric_ratelaws=true,
                      include_zero_odes=true, checks = false, kwargs...)

    flatrs = Catalyst.flatten(rs)
    error_if_constraints(SDESystem, flatrs)

    if noise_scaling isa AbstractArray
        (length(noise_scaling)!=numreactions(flatrs)) &&
        error("The number of elements in 'noise_scaling' must be equal " *
              "to the number of reactions in the flattened reaction system.")
        if !(noise_scaling isa Symbolics.Arr)
            noise_scaling = value.(noise_scaling)
        end
    elseif !isnothing(noise_scaling)
        noise_scaling = fill(value(noise_scaling),numreactions(flatrs))
    end

    eqs = assemble_drift(flatrs; combinatoric_ratelaws=combinatoric_ratelaws,
                                 include_zero_odes=include_zero_odes)
    noiseeqs = assemble_diffusion(flatrs, noise_scaling;
                                  combinatoric_ratelaws=combinatoric_ratelaws)
    SDESystem(eqs, noiseeqs, get_iv(flatrs), get_states(flatrs),
              (noise_scaling===nothing) ? get_ps(flatrs) : union(get_ps(flatrs), toparam(noise_scaling));
              name=name,
              defaults=MT.defaults(flatrs),
              observed=MT.observed(flatrs),
              checks = checks,
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
  the ratelaw is `k*S*(S-1)`, i.e. the rate law is not normalized by the scaling
  factor.
"""
function Base.convert(::Type{<:JumpSystem},rs::ReactionSystem;
                      name=nameof(rs), combinatoric_ratelaws=true, checks = false, kwargs...)

    flatrs = Catalyst.flatten(rs)
    error_if_constraints(JumpSystem, flatrs)

    eqs = assemble_jumps(flatrs; combinatoric_ratelaws=combinatoric_ratelaws)
    JumpSystem(eqs, get_iv(flatrs), get_states(flatrs), get_ps(flatrs); name=name,
               defaults=MT.defaults(flatrs), observed=MT.observed(flatrs),
               checks = checks, kwargs...)
end


### Converts a reaction system to ODE or SDE problems ###


# ODEProblem from AbstractReactionNetwork
function DiffEqBase.ODEProblem(rs::ReactionSystem, u0, tspan, p=DiffEqBase.NullParameters(), args...;
                               check_length=false, name=nameof(rs), combinatoric_ratelaws=true,
                               include_zero_odes=true, remove_conserved=false, checks=false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap  = symmap_to_varmap(rs, p)
    osys  = convert(ODESystem, rs; name, combinatoric_ratelaws, include_zero_odes, checks,
                                   remove_conserved)
    return ODEProblem(osys, u0map, tspan, pmap, args...; check_length, kwargs...)
end

# NonlinearProblem from AbstractReactionNetwork
function DiffEqBase.NonlinearProblem(rs::ReactionSystem, u0, p=DiffEqBase.NullParameters(), args...;
                                     name=nameof(rs), combinatoric_ratelaws=true, include_zero_odes=true,
                                     remove_conserved=false, checks = false, check_length=false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap  = symmap_to_varmap(rs, p)
    nlsys = convert(NonlinearSystem, rs; name, combinatoric_ratelaws, include_zero_odes,
                                         checks, remove_conserved)
    return NonlinearProblem(nlsys, u0map, pmap, args...; check_length, kwargs...)
end


# SDEProblem from AbstractReactionNetwork
function DiffEqBase.SDEProblem(rs::ReactionSystem, u0, tspan, p=DiffEqBase.NullParameters(), args...;
                               noise_scaling=nothing, name=nameof(rs), combinatoric_ratelaws=true,
                               include_zero_odes=true, checks = false, check_length=false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap  = symmap_to_varmap(rs, p)
    p_matrix = zeros(length(get_states(rs)), length(get_eqs(rs)))
    sde_sys  = convert(SDESystem, rs; noise_scaling, name, combinatoric_ratelaws,
                                      include_zero_odes, checks)
    return SDEProblem(sde_sys, u0map, tspan, pmap, args...; check_length,
                                                            noise_rate_prototype=p_matrix, kwargs...)
end

# DiscreteProblem from AbstractReactionNetwork
function DiffEqBase.DiscreteProblem(rs::ReactionSystem, u0, tspan::Tuple, p=DiffEqBase.NullParameters(),
                                    args...; name=nameof(rs), combinatoric_ratelaws=true, checks = false, kwargs...)
    u0map = symmap_to_varmap(rs, u0)
    pmap  = symmap_to_varmap(rs, p)
    jsys  = convert(JumpSystem, rs; name, combinatoric_ratelaws, checks)
    return DiscreteProblem(jsys, u0map, tspan, pmap, args...; kwargs...)
end

# JumpProblem from AbstractReactionNetwork
function DiffEqJump.JumpProblem(rs::ReactionSystem, prob, aggregator, args...;
                                name=nameof(rs), combinatoric_ratelaws=true, checks = false, kwargs...)
    jsys = convert(JumpSystem, rs; name, combinatoric_ratelaws, checks)
    return JumpProblem(jsys, prob, aggregator, args...; kwargs...)
end

# SteadyStateProblem from AbstractReactionNetwork
function DiffEqBase.SteadyStateProblem(rs::ReactionSystem, u0, p=DiffEqBase.NullParameters(), args...;
                                       check_length=false, name=nameof(rs), combinatoric_ratelaws=true,
                                       remove_conserved=false, include_zero_odes=true, checks=false, kwargs...)

    u0map = symmap_to_varmap(rs, u0)
    pmap = symmap_to_varmap(rs, p)
    osys = convert(ODESystem, rs; name, combinatoric_ratelaws, include_zero_odes, checks,
                                  remove_conserved)
    return SteadyStateProblem(osys, u0map, pmap, args...; check_length, kwargs...)
end

####################### dependency graph utilities ########################

# determine which species a reaction depends on
function ModelingToolkit.get_variables!(deps::Set, rx::Reaction, variables)
    (rx.rate isa Symbolic) && get_variables!(deps, rx.rate, variables)
    for s in rx.substrates
        push!(deps, s)
    end
    deps
end

# determine which species a reaction modifies
function ModelingToolkit.modified_states!(mstates, rx::Reaction, sts::Set)
    for (species,stoich) in rx.netstoich
        (species in sts) && push!(mstates, species)
    end
    mstates
end


########################## Compositional Tooling ###########################
function getsubsyseqs!(eqs::Vector{Equation}, sys)
    if sys isa ReactionSystem
        push!(eqs, get_eqs(get_constraints(sys)))
    else
        push!(eqs, get_eqs(sys))
    end

    for subsys in get_systems(sys)
        getsubsyseqs!(eqs, subsys)
    end
    eqs
end
function getsubsyseqs(sys)
    getsubsyseqs!(Vector{Equation}(), sys)
end

function getsubsystypes!(typeset::Set{Type}, sys::T) where {T <: MT.AbstractSystem}
    push!(typeset, T)
    csys = (sys isa ReactionSystem) ? get_constraints(sys) : nothing
    csys !== nothing && getsubsystypes!(typeset, csys)
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
- All `Reaction`s within subsystems are namespaced and merged into the list of
  `Reactions` of `rs`. The merged list is then available as `reactions(rs)` or
  `get_eqs(rs)`.
- All algebraic equations are merged into a `NonlinearSystem` or `ODESystem`
  stored as `get_constraints(rs)`. If `get_constraints !== nothing` then the
  algebraic equations are merged with the current constraints in a system of the
  same type as the current constraints, otherwise the new constraint system is
  an `ODESystem`.
- Currently only `ReactionSystem`s, `NonlinearSystem`s and `ODESystem`s are
  supported as sub-systems when flattening.
- `rs.networkproperties` is reset upon flattening.
"""
function MT.flatten(rs::ReactionSystem; name=nameof(rs))

    isempty(get_systems(rs)) && return rs

    # right now only NonlinearSystems and ODESystems can be handled as subsystems
    subsys_types  = getsubsystypes(rs)
    allowed_types = (ReactionSystem, NonlinearSystem, ODESystem)
    all(T -> any(T .<: allowed_types), subsys_types) ||
        error("flattening is currently only supported for subsystems mixing ReactionSystems, NonlinearSystems and ODESystems.")

    specs      = species(rs)
    sts        = states(rs)
    reactionps = reactionparams(rs)
    ps         = parameters(rs)
    alleqs     = equations(rs)
    rxs        = Reaction[rx for rx in alleqs if rx isa Reaction]

    # constraints = states, parameters and equations that do not appear in Reactions
    csts = setdiff(sts, specs)
    cps  = setdiff(ps, reactionps)
    ceqs = Equation[eq for eq in alleqs if eq isa Equation]
    if isempty(csts) && isempty(cps) && isempty(ceqs)
        newcsys = nothing
    else
        if ODESystem in subsys_types
            newcsys = ODESystem(ceqs, get_iv(rs), csts, cps; name=name)
        else # must be a NonlinearSystem
            any(T -> T <: NonlinearSystem, subsys_types) ||
                error("Error, found constraint Equations but no ODESystem or NonlinearSystem associated with them.")
            newcsys = NonlinearSystem(ceqs, csts, cps; name=name)
        end
    end

    ReactionSystem(rxs, get_iv(rs), specs, reactionps; observed = MT.observed(rs),
                                                       name = name,
                                                       defaults = MT.defaults(rs),
                                                       checks = false,
                                                       constraints = newcsys)
end

"""
    ModelingToolkit.extend(sys::Union{NonlinearSystem,ODESystem}, rs::ReactionSystem; name::Symbol=nameof(sys))

Extends the indicated [`ReactionSystem`](@ref) with a
`ModelingToolkit.NonlinearSystem` or `ModelingToolkit.ODESystem`, which will be
stored internally as constraint equations.

Notes:
- Returns a new `ReactionSystem` and does not modify `rs`.
- By default, the new `ReactionSystem` will have the same name as `sys`.
"""
function ModelingToolkit.extend(sys::Union{NonlinearSystem,ODESystem}, rs::ReactionSystem; name::Symbol=nameof(sys))
    csys = (get_constraints(rs) === nothing) ? sys : extend(sys, get_constraints(rs))
    ReactionSystem(get_eqs(rs), get_iv(rs), get_states(rs), get_ps(rs);
                    observed = get_observed(rs), name = name,
                    systems = get_systems(rs), defaults = get_defaults(rs),
                    checks=false, constraints=csys)
end

"""
    ModelingToolkit.extend(sys::ReactionSystem, rs::ReactionSystem; name::Symbol=nameof(sys))

Extends the indicated [`ReactionSystem`](@ref) with another `ReactionSystem`.
Similar to calling `merge!` except constraint systems are allowed (and will also
be merged together).

Notes:
- Returns a new `ReactionSystem` and does not modify `rs`.
- By default, the new `ReactionSystem` will have the same name as `sys`.
"""
function ModelingToolkit.extend(sys::ReactionSystem, rs::ReactionSystem; name::Symbol=nameof(sys))
    eqs  = union(get_eqs(rs), get_eqs(sys))
    sts  = union(get_states(rs), get_states(sys))
    ps   = union(get_ps(rs), get_ps(sys))
    obs  = union(get_observed(rs), get_observed(sys))
    defs = merge(get_defaults(rs), get_defaults(sys)) # prefer `sys`
    syss = union(get_systems(rs), get_systems(sys))

    csys = get_constraints(rs)
    csys2 = get_constraints(sys)
    if csys === nothing
        newcsys = csys2
    else
        newcsys = (csys2 === nothing) ? csys : extend(csys2, csys, name)
    end

    ReactionSystem(eqs, get_iv(rs), sts, ps; observed = obs, name = name,
                    systems = syss, defaults = defs, checks=false, constraints=newcsys)
end
