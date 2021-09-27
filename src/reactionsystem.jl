
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
struct Reaction{S, T <: Number}
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
    if isnothing(prods)
        prods = Vector{Term}()
        !isnothing(prodstoich) && error("If products are nothing, product stiocihometries have to be so too.")
        prodstoich = typeof(substoich)()
    end
    subs = value.(subs)
    prods = value.(prods)
    ns = isnothing(netstoich) ? get_netstoich(subs, prods, substoich, prodstoich) : netstoich
    Reaction(value(rate), subs, prods, substoich, prodstoich, ns, only_use_rate)
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
            if stoich[1] == 1
                print(io, ModelingToolkit.operation(spec))
            else
                print(io, stoich[i], ModelingToolkit.operation(spec))
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
    Reaction(namespace_expr(rx.rate, name), 
             namespace_expr(rx.substrates, name),
             namespace_expr(rx.products, name),
             rx.substoich, rx.prodstoich,            
             [namespace_expr(n[1],name) => n[2] for n in rx.netstoich], rx.only_use_rate)
end

# calculates the net stoichiometry of a reaction as a vector of pairs (sub,substoich)
function get_netstoich(subs, prods, sstoich, pstoich)
    # stoichiometry as a Dictionary
    nsdict = Dict{Any, eltype(sstoich)}(sub => -sstoich[i] for (i,sub) in enumerate(subs))
    for (i,p) in enumerate(prods)
        coef = pstoich[i]
        @inbounds nsdict[p] = haskey(nsdict, p) ? nsdict[p] + coef : coef
    end

    # stoichiometry as a vector
    ns = [el for el in nsdict if el[2] != zero(el[2])]

    ns
end

"""
$(TYPEDEF)

A system of chemical reactions.

# Fields
$(FIELDS)

# Example
Continuing from the example in the [`Reaction`](@ref) definition:
```julia
@named rs = ReactionSystem(rxs, t, [A,B,C,D], k)
```

Notes:
- ReactionSystems currently do rudimentary unit checking, requiring that
  all species have the same units, and all reactions have rate laws with 
  units of (species units) / (time units). Unit checking can be disabled
  by passing the keyword argument `checks=false`.
"""
struct ReactionSystem{U <: Union{Nothing,MT.AbstractSystem}} <: MT.AbstractTimeDependentSystem
    """The reactions defining the system."""
    eqs::Vector{Reaction}
    """Independent variable (usually time)."""
    iv::Any
    """Dependent (state) variables representing amount of each species. Must not contain the independent variable."""
    states::Vector
    """Parameter variables. Must not contain the independent variable."""
    ps::Vector
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

    function ReactionSystem(eqs, iv, states, ps, observed, name, systems, defaults, connection_type, csys; checks::Bool = true)

        iv′ = value(iv)
        states′ = value.(states)
        ps′ = value.(ps)

        if checks
            check_variables(states′, iv′)
            check_parameters(ps′, iv′)
            # check_units(eqs)    # disable as check the newly generated system below
        end
        rs = new{typeof(csys)}(collect(eqs), iv′, states′, ps′, observed, name, systems, defaults, connection_type, csys)
        checks && validate(rs)

        rs
    end
end

function ReactionSystem(eqs, iv, species, ps;
                        observed = [],
                        systems = [],
                        name = nothing,
                        default_u0=Dict(),
                        default_p=Dict(),
                        defaults=_merge(Dict(default_u0), Dict(default_p)),
                        connection_type=nothing,
                        checks = true, 
                        constraints = nothing)
    name === nothing && throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))

    ReactionSystem(eqs, iv, species, ps, observed, name, systems, 
                   defaults, connection_type, constraints; checks = checks)
end

function ReactionSystem(iv; kwargs...)
    ReactionSystem(Reaction[], iv, [], []; kwargs...)
end


####################### ModelingToolkit inherited accessors #############################

get_constraints(sys::ReactionSystem) = getfield(sys, :constraints)
has_constraints(sys::ReactionSystem) = isdefined(sys, :constraints)

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

######################## Conversion to ODEs/SDEs/jump, etc ##############################

"""
    oderatelaw(rx; combinatoric_ratelaw=true)

Given a [`Reaction`](@ref), return the symbolic reaction rate law used in
generated ODEs for the reaction. Note, for a reaction defined by

`k*X*Y, X+Z --> 2X + Y`

the expression that is returned will be `k*X(t)^2*Y(t)*Z(t)`. For a reaction
of the form

`k, 2X+3Y --> Z`

the expression that is returned will be `k * (X(t)^2/2) * (Y(t)^3/6)`.

Notes:
- Allocates
- `combinatoric_ratelaw=true` uses factorial scaling factors in calculating the rate
    law, i.e. for `2S -> 0` at rate `k` the ratelaw would be `k*S^2/2!`. If
    `combinatoric_ratelaw=false` then the ratelaw is `k*S^2`, i.e. the scaling factor is
    ignored.
"""
function oderatelaw(rx; combinatoric_ratelaw=true)
    @unpack rate, substrates, substoich, only_use_rate = rx
    rl = rate
    if !only_use_rate
        coef = one(eltype(substoich))
        for (i,stoich) in enumerate(substoich)
            coef *= factorial(stoich)
            rl   *= isone(stoich) ? substrates[i] : substrates[i]^stoich
        end
        combinatoric_ratelaw && (!isone(coef)) && (rl /= coef)
    end
    rl
end

function assemble_oderhs(rs; combinatoric_ratelaws=true)
    sts = get_states(rs)
    species_to_idx = Dict((x => i for (i,x) in enumerate(sts)))
    rhsvec         = Any[0 for i in eachindex(sts)]

    for rx in get_eqs(rs)
        rl = oderatelaw(rx; combinatoric_ratelaw=combinatoric_ratelaws)
        for (spec,stoich) in rx.netstoich
            i = species_to_idx[spec]
            if _iszero(rhsvec[i])
                signedrl  = (stoich > zero(stoich)) ? rl : -rl
                rhsvec[i] = isone(abs(stoich)) ? signedrl : stoich * rl
            else
                Δspec     = isone(abs(stoich)) ? rl : abs(stoich) * rl
                rhsvec[i] = (stoich > zero(stoich)) ? (rhsvec[i] + Δspec) : (rhsvec[i] - Δspec)
            end
        end
    end

    rhsvec
end

function assemble_drift(rs; combinatoric_ratelaws=true, as_odes=true, include_zero_odes=true)
    rhsvec = assemble_oderhs(rs; combinatoric_ratelaws=combinatoric_ratelaws)
    if as_odes
        D   = Differential(get_iv(rs))
        eqs = [Equation(D(x),rhs) for (x,rhs) in zip(get_states(rs),rhsvec) if (include_zero_odes || (!_iszero(rhs)))]
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
            i            = species_to_idx[spec]
            signedrlsqrt = (stoich > zero(stoich)) ? rlsqrt : -rlsqrt
            eqs[i,j]     = isone(abs(stoich)) ? signedrlsqrt : stoich * rlsqrt
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
        coef = one(eltype(substoich))
        for (i,stoich) in enumerate(substoich)
            s   = substrates[i]
            rl *= s
            isone(stoich) && continue
            for i in one(stoich):(stoich-one(stoich))
                rl *= (s - i)
            end
            combinatoric_ratelaw && (coef *= factorial(stoich))
        end
        !isone(coef) && (rl /= coef)
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
on any chemical species that correspond to states of the system, and does not depend
explicitly on the independent variable (usually time).

# Arguments
- `rx`, the [`Reaction`](@ref).
- `rs`, a [`ReactionSystem`](@ref) containing the reaction.
- Optional: `rxvars`, `Variable`s which are not in `rxvars` are ignored as possible dependencies.
- Optional: `haveivdep`, `true` if the [`Reaction`](@ref) `rate` field explicitly depends on the independent variable.
- Optional: `stateset`, set of states which if the rxvars are within mean rx is non-mass action.
"""
function ismassaction(rx, rs; rxvars = get_variables(rx.rate),
                              haveivdep = any(var -> isequal(get_iv(rs),var), rxvars),
                              stateset = Set(get_states(rs)))
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
function make_systems_with_type!(systems::Vector{T}, rs::ReactionSystem, include_zero_odes=true) where {T <: MT.AbstractSystem}
    resize!(systems, length(get_systems(rs)))
    for (i,sys) in enumerate(get_systems(rs))
        if sys isa ReactionSystem
            systems[i] = convert(T, sys, include_zero_odes=include_zero_odes)
        elseif sys isa T
            systems[i] = sys
        else
            try 
                if (T <: MT.AbstractTimeDependentSystem) && (sys isa MT.AbstractTimeIndependentSystem)
                    systems[i] = MT.convert_system(T, sys, get_iv(rs))
                else
                    systems[i] = MT.convert_system(T, sys)
                end
            catch e
                error("ModelingToolkit does not currently support convert_system($T, $(typeof(sys)))")
            end
        end
    end    
    systems
end

# merge constraint eqs, states and ps into the top-level eqs, states and ps
function addconstraints!(eqs, rs::ReactionSystem)   
    csys = get_constraints(rs)     
    sts  = get_states(rs); ps = get_ps(rs)

    if csys !== nothing
        csts = get_states(csys); cps = get_ps(csys); ceqs = get_eqs(csys)
        sts  = isempty(csts) ? sts : [sts; csts]        
        ps   = isempty(cps) ? ps : [ps; cps]
        (!isempty(ceqs)) && append!(eqs,ceqs)
    end

    eqs,sts,ps
end

# used by systems that don't support constraint equations currently
function error_if_constraints(::Type{T}, sys::ReactionSystem) where {T <: MT.AbstractSystem}
    (get_constraints(sys) === nothing) || 
            error("Can not convert to a system of type ", T, " when there are constraints.")
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
                      checks=false, kwargs...)
    eqs        = assemble_drift(rs; combinatoric_ratelaws=combinatoric_ratelaws, 
                                    include_zero_odes=include_zero_odes)                                 
    eqs,sts,ps = addconstraints!(eqs, rs)
    systems    = make_systems_with_type!(Vector{ODESystem}(), rs, include_zero_odes)
    ODESystem(eqs, get_iv(rs), sts, ps; name=name, systems=systems, 
              defaults=get_defaults(rs), checks=checks, kwargs...)
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
function Base.convert(::Type{<:NonlinearSystem},rs::ReactionSystem;
                      name=nameof(rs), combinatoric_ratelaws=true, include_zero_odes=true, 
                      checks = false, kwargs...)
    eqs        = assemble_drift(rs; combinatoric_ratelaws=combinatoric_ratelaws, as_odes=false, 
                                    include_zero_odes=include_zero_odes)
    eqs,sts,ps = addconstraints!(eqs, rs)
    systems    = make_systems_with_type!(Vector{NonlinearSystem}(), rs, include_zero_odes)
    NonlinearSystem(eqs, sts, ps; name=name, systems=systems, 
                    defaults=get_defaults(rs), checks = checks, kwargs...)
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

    error_if_constraints(SDESystem, rs)

    if noise_scaling isa AbstractArray
        (length(noise_scaling)!=numreactions(rs)) &&
        error("The number of elements in 'noise_scaling' must be equal " *
              "to the number of reactions in the reaction system.")
        if !(noise_scaling isa Symbolics.Arr)
            noise_scaling = value.(noise_scaling)
        end
    elseif !isnothing(noise_scaling)
        noise_scaling = fill(value(noise_scaling),numreactions(rs))
    end

    eqs      = assemble_drift(rs; combinatoric_ratelaws=combinatoric_ratelaws, 
                                  include_zero_odes=include_zero_odes)
    noiseeqs = assemble_diffusion(rs,noise_scaling;
                                  combinatoric_ratelaws=combinatoric_ratelaws)
    systems = make_systems_with_type!(Vector{SDESystem}(), rs, include_zero_odes)
    SDESystem(eqs, noiseeqs, get_iv(rs), get_states(rs),
              (noise_scaling===nothing) ? get_ps(rs) : union(get_ps(rs), toparam(noise_scaling));
              name=name, 
              systems=systems,
              defaults=get_defaults(rs),
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
    
    error_if_constraints(JumpSystem, rs)
    isempty(get_systems(rs)) || error("Conversion to JumpSystems with subsystems is not currently supported.")

    eqs     = assemble_jumps(rs; combinatoric_ratelaws=combinatoric_ratelaws)
    systems = convert.(JumpSystem, get_systems(rs))    
    JumpSystem(eqs, get_iv(rs), get_states(rs), get_ps(rs); name=name, systems=systems, 
               defaults=get_defaults(rs), checks = checks, kwargs...)
end


### Converts a reaction system to ODE or SDE problems ###


# ODEProblem from AbstractReactionNetwork
function DiffEqBase.ODEProblem(rs::ReactionSystem, u0, tspan, p=DiffEqBase.NullParameters(), args...; 
                               check_length=false, kwargs...)
    return ODEProblem(convert(ODESystem,rs; kwargs...),u0,tspan,p, args...; check_length, kwargs...)
end

# NonlinearProblem from AbstractReactionNetwork
function DiffEqBase.NonlinearProblem(rs::ReactionSystem, u0, p=DiffEqBase.NullParameters(), args...; 
                                     check_length=false, kwargs...)
    return NonlinearProblem(convert(NonlinearSystem,rs; kwargs...), u0, p, args...; check_length, kwargs...)
end


# SDEProblem from AbstractReactionNetwork
function DiffEqBase.SDEProblem(rs::ReactionSystem, u0, tspan, p=DiffEqBase.NullParameters(), args...; 
                               noise_scaling=nothing, kwargs...)
    sde_sys  = convert(SDESystem,rs;noise_scaling=noise_scaling, kwargs...)
    p_matrix = zeros(length(get_states(rs)), length(get_eqs(rs)))
    return SDEProblem(sde_sys,u0,tspan,p,args...; noise_rate_prototype=p_matrix,kwargs...)
end

# DiscreteProblem from AbstractReactionNetwork
function DiffEqBase.DiscreteProblem(rs::ReactionSystem, u0, tspan::Tuple, p=DiffEqBase.NullParameters(), 
                                    args...; kwargs...)
    return DiscreteProblem(convert(JumpSystem,rs; kwargs...), u0,tspan,p, args...; kwargs...)
end

# JumpProblem from AbstractReactionNetwork
function DiffEqJump.JumpProblem(rs::ReactionSystem, prob, aggregator, args...; kwargs...)
    return JumpProblem(convert(JumpSystem,rs; kwargs...), prob, aggregator, args...; kwargs...)
end

# SteadyStateProblem from AbstractReactionNetwork
function DiffEqBase.SteadyStateProblem(rs::ReactionSystem, u0, p=DiffEqBase.NullParameters(), args...; 
                                       kwargs...)
    return SteadyStateProblem(ODEFunction(convert(ODESystem,rs; kwargs...)),u0,p, args...; kwargs...)
end

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

function flatten(rs::ReactionSystem)
    systems = get_systems(rs)
    isempty(systems) && return rs
    
    all(T -> any(T .<: (ReactionSystem,NonlinearSystem)), getsubsystypes(rs)) || 
        error("flattening is currently only supported for subsystems mixing ReactionSystems and NonlinearSystems.")
    
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
    ReactionSystem(rxs, get_iv(rs), specs, reactionps;
                   observed = MT.observed(rs),                    
                   name = nameof(rs),
                   defaults = MT.defaults(rs),
                   checks = false,
                   constraints = NonlinearSystem(ceqs,csts,cps,name=nameof(rs)))
end