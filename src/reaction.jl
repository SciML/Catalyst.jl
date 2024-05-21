### Species Symbolic Variable Definition ###

# Defines species-related metadata.
struct ParameterConstantSpecies end
struct VariableBCSpecies end
struct VariableSpecies end
Symbolics.option_to_metadata_type(::Val{:isconstantspecies}) = ParameterConstantSpecies
Symbolics.option_to_metadata_type(::Val{:isbcspecies}) = VariableBCSpecies
Symbolics.option_to_metadata_type(::Val{:isspecies}) = VariableSpecies

# Defines species-related metadata getters.
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

# Other species functions.

"""
    isvalidreactant(s)

Test if a species is valid as a reactant (i.e. a species variable or a constant parameter).
"""
isvalidreactant(s) = MT.isparameter(s) ? isconstant(s) : (isspecies(s) && !isconstant(s))


### Reaction Constructor Functions ###

# Checks if a metadata input has an entry :only_use_rate => true
function metadata_only_use_rate_check(metadata)
    only_use_rate_idx = findfirst(:only_use_rate == entry[1] for entry in metadata)
    isnothing(only_use_rate_idx) && return false
    return Bool(metadata[only_use_rate_idx][2])
end

# Used to promote a vector to the appropriate type. Takes the `vec == Any[]` case into account by
# returning an empty vector of the appropriate type.
function promote_reaction_vector(vec, type)
    isempty(vec) && (return type[])
    type[value(v) for v in vec]
end

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

### Reaction Structure ###

"""
$(TYPEDEF)

One chemical reaction.

# Fields
$(FIELDS)

# Examples

```julia
using Catalyst
t = default_t()
@parameters k[1:20]
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
struct Reaction{T}
    """The rate function (excluding mass action terms)."""
    rate::Any
    """Reaction substrates."""
    substrates::Vector{BasicSymbolic{Real}}
    """Reaction products."""
    products::Vector{BasicSymbolic{Real}}
    """The stoichiometric coefficients of the reactants."""
    substoich::Vector{T}
    """The stoichiometric coefficients of the products."""
    prodstoich::Vector{T}
    """The net stoichiometric coefficients of all species changed by the reaction."""
    netstoich::Vector{Pair{BasicSymbolic{Real}, T}}
    """
    `false` (default) if `rate` should be multiplied by mass action terms to give the rate law.
    `true` if `rate` represents the full reaction rate law.
    """
    only_use_rate::Bool
    """
    Contain additional data, such whenever the reaction have a specific noise-scaling expression for
    the chemical Langevin equation.
    """
    metadata::Vector{Pair{Symbol, Any}}
end

# Five-argument constructor accepting rate, substrates, and products, and their stoichiometries.
function Reaction(rate, subs::Vector, prods::Vector, substoich::Vector{S}, prodstoich::Vector{T};
                   netstoich = [], metadata = Pair{Symbol, Any}[], 
                   only_use_rate = metadata_only_use_rate_check(metadata), kwargs...) where {S,T}
    # Error checks.
    isempty(subs) && isempty(prods) &&
        throw(ArgumentError("A reaction requires either a non-empty substrate or product vector."))
    length(subs) != length(substoich) &&
        throw(ArgumentError("The substrate vector ($(subs)) and the substrate stoichiometry vector ($(substoich)) must have equal length."))
    length(prods) != length(prodstoich) &&
        throw(ArgumentError("The product vector ($(prods)) and the product stoichiometry vector ($(prodstoich)) must have equal length."))
    allunique(subs) ||
        throw(ArgumentError("Substrates can not be repeated in the list provided to `Reaction`, please modify the stoichiometry for any repeated substrates instead."))
    allunique(prods) ||
        throw(ArgumentError("Products can not be repeated in the list provided to `Reaction`, please modify the stoichiometry for any repeated products instead."))

    # Ensures everything have uniform and correct types.
    subs = promote_reaction_vector(subs, BasicSymbolic{Real})
    prods =  promote_reaction_vector(prods, BasicSymbolic{Real})
    stoich_type = promote_type(S, T)
    (stoich_type <: Num) && (stoich_type = Any)
    substoich = promote_reaction_vector(substoich, stoich_type)
    prodstoich = promote_reaction_vector(prodstoich, stoich_type)

    # Checks that all reactants are valid.
    if !(all(isvalidreactant, subs) && all(isvalidreactant, prods))
        badsts = union(filter(!isvalidreactant, subs), filter(!isvalidreactant, prods))
        throw(ArgumentError("To be a valid substrate or product, non-constant species must be declared via @species, while constant species must be parameters with the isconstantspecies metadata. The following reactants do not follow this convention:\n $badsts"))
    end
    
    # Computes the net stoichiometries.
    netstoich = if isnothing(netstoich)
        get_netstoich(subs, prods, substoich, prodstoich)
    else
        if typeof(netstoich) != Vector{Pair{BasicSymbolic{Real}, stoich_type}}
            netstoich = Pair{BasicSymbolic{Real}, stoich_type}[
                            value(ns[1]) => convert(stoich_type, ns[2]) for ns in netstoich]
        end
        netstoich
    end

    # Handles metadata (check that all entries are unique, remove potential `only_use_rate`
    # entries, and converts to the required type.
    if !allunique(entry[1] for entry in metadata)
        error("Repeated entries for the same metadata encountered in the following metadata set: $([entry[1] for entry in metadata]).")
    end
    if any(:only_use_rate == entry[1] for entry in metadata) 
        deleteat!(metadata, findfirst(:only_use_rate == entry[1] for entry in metadata))
    end
    metadata = convert(Vector{Pair{Symbol, Any}}, metadata)

    Reaction{stoich_type}(value(rate), subs, prods, substoich, prodstoich, netstoich, only_use_rate, metadata)
end

# Three-argument constructor. Handles the case where no stoichiometries is given
# (by assuming that all stoichiometries are `1`).
function Reaction(rate, subs::Vector, prods::Vector; kwargs...)
    Reaction(rate, subs, prods, ones(Int, length(subs)), ones(Int, length(prods)); kwargs...)
end

# Handles cases where `nothing` is given (instead of an empty vector). 
function Reaction(rate, subs::Vector, prods::Nothing, substoich::Vector, prodstoich::Nothing; kwargs...)
    Reaction(rate, subs, Int64[], substoich, Int64[]; kwargs...)
end
function Reaction(rate, subs::Nothing, prods::Vector, substoich::Nothing, prodstoich::Vector; kwargs...)
    Reaction(rate, Int64[], prods, Int64[], prodstoich; kwargs...)
end
function Reaction(rate, subs::Vector, prods::Nothing; kwargs...)
    Reaction(rate, subs, Int64[]; kwargs...)
end
function Reaction(rate, subs::Nothing, prods::Vector; kwargs...)
    Reaction(rate, Int64[], prods; kwargs...)
end

# Union type for `Reaction`s and `Equation`s.
const CatalystEqType = Union{Reaction, Equation}

### Base Function Dispatches ###

# Used by `Base.show`.
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

# Show function for `Reaction`s.
function Base.show(io::IO, rx::Reaction)
    print(io, rx.rate, ", ")
    print_rxside(io, rx.substrates, rx.substoich)
    arrow = rx.only_use_rate ? "⇒" : "-->"
    print(io, " ", arrow, " ")
    print_rxside(io, rx.products, rx.prodstoich)
end

"""
    ==(rx1::Reaction, rx2::Reaction)

Tests whether two [`Reaction`](@ref)s are identical.

Notes:
- Ignores the order in which stoichiometry components are listed.
- *Does not* currently simplify rates, so a rate of `A^2+2*A+1` would be
    considered different than `(A+1)^2`.
"""
function (==)(rx1::Reaction, rx2::Reaction)
    isequal(rx1.rate, rx2.rate) || return false
    issetequal(zip(rx1.substrates, rx1.substoich), zip(rx2.substrates, rx2.substoich)) ||
        return false
    issetequal(zip(rx1.products, rx1.prodstoich), zip(rx2.products, rx2.prodstoich)) ||
        return false
    issetequal(rx1.netstoich, rx2.netstoich) || return false
    rx1.only_use_rate == rx2.only_use_rate
end

# Hash function.
function hash(rx::Reaction, h::UInt)
    h = Base.hash(rx.rate, h)
    for s in Iterators.flatten((rx.substrates, rx.products))
        h ⊻= hash(s)
    end
    for s in Iterators.flatten((rx.substoich, rx.prodstoich))
        h ⊻= hash(s)
    end
    for s in rx.netstoich
        h ⊻= hash(s)
    end
    Base.hash(rx.only_use_rate, h)
end


### ModelingToolkit Function Dispatches ###

# Used by ModelingToolkit.namespace_equation.
function apply_if_nonempty(f, v)
    isempty(v) && return v
    s = similar(v)
    map!(f, s, v)
    s
end

# Returns a name-spaced version of a reaction.
function MT.namespace_equation(rx::Reaction, name; kw...)
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
    Reaction(rate, subs, prods, substoich, prodstoich, netstoich, rx.only_use_rate, rx.metadata)
end

# Overwrites equation-type functions to give the correct input for `Reaction`s.
MT.is_diff_equation(rx::Reaction) = false
MT.is_alg_equation(rx::Reaction) = false

"""
    get_symbolics(set, rx::Reaction)

Returns all symbolic variables that are part of a reaction. This includes all variables 
encountered in:
    - Rates.
    - Among substrates and products.
    - Among stoichiometries.
    - Among potential noise scaling metadata.
"""
function get_symbolics(rx::Reaction)
    return ModelingToolkit.get_variables!([], rx)
end

"""
    get_variables!(set, rx::Reaction)

Adds all symbolic variables that are part of a reaction to set. This includes all variables 
encountered in:
    - Rates.
    - Among substrates and products.
    - Among stoichiometries.
    - Among potential noise scaling metadata.
"""
function ModelingToolkit.get_variables!(set, rx::Reaction)
    if isdefined(Main, :Infiltrator)
        Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
      end
    get_variables!(set, rx.rate)
    foreach(sub -> push!(set, sub), rx.substrates)
    foreach(prod -> push!(set, prod), rx.products)
    for stoichs in (rx.substoich, rx.prodstoich), stoich in stoichs
        (stoich isa BasicSymbolic) && get_variables!(set, stoich)
    end
    if hasnoisescaling(rx)
        get_variables!(set, getnoisescaling(rx))
    end
    return (set isa AbstractVector) ? unique!(set) : set
end

### Dependency-related Functions ###

# determine which unknowns a reaction depends on
function MT.get_variables!(deps::Set, rx::Reaction, variables)
    (rx.rate isa Symbolic) && get_variables!(deps, rx.rate, variables)
    for s in rx.substrates
        # parametric stoichiometry means may have a parameter as a substrate
        any(isequal(s), variables) && push!(deps, s)
    end
    deps
end

# determine which species a reaction modifies
function MT.modified_unknowns!(munknowns, rx::Reaction, sts::Set)
    for (species, stoich) in rx.netstoich
        (species in sts) && push!(munknowns, species)
    end
    munknowns
end

function MT.modified_unknowns!(munknowns, rx::Reaction, sts::AbstractVector)
    for (species, stoich) in rx.netstoich
        any(isequal(species), sts) && push!(munknowns, species)
    end
    munknowns
end


### `Reaction`-specific Functions ### 

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


### Reaction Metadata Implementation ###
# These are currently considered internal, but can be used by public accessor functions like getnoisescaling.

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


### Implemented Reaction Metadata ###

# Noise scaling.
"""
hasnoisescaling(reaction::Reaction)

Returns `true` if the input reaction has the `noise_scaing` metadata field assigned, else `false`.

Arguments:
- `reaction`: The reaction we wish to check for the `noise_scaing` metadata field.

Example:
```julia
reaction = @reaction k, 0 --> X, [noise_scaling=0.0]
hasnoisescaling(reaction)
```
"""
function hasnoisescaling(reaction::Reaction)
    return hasmetadata(reaction, :noise_scaling)
end

"""
getnoisescaling(reaction::Reaction)

Returns `noise_scaing` metadata field for the input reaction.

Arguments:
- `reaction`: The reaction we wish to retrieve the `noise_scaing` metadata field.

Example:
```julia
reaction = @reaction k, 0 --> X, [noise_scaling=0.0]
getnoisescaling(reaction)
```
"""
function getnoisescaling(reaction::Reaction)
    if hasnoisescaling(reaction)
        return getmetadata(reaction, :noise_scaling)
    else
        error("Attempts to access noise_scaling metadata field for a reaction which does not have a value assigned for this metadata.")
    end
end

# Description.
"""
hasdescription(reaction::Reaction)

Returns `true` if the input reaction has the `description` metadata field assigned, else `false`.

Arguments:
- `reaction`: The reaction we wish to check for the `description` metadata field.

Example:
```julia
reaction = @reaction k, 0 --> X, [description="A reaction"]
hasdescription(reaction)
```
"""
function hasdescription(reaction::Reaction)
    return hasmetadata(reaction, :description)
end

"""
getdescription(reaction::Reaction)

Returns `description` metadata field for the input reaction.

Arguments:
- `reaction`: The reaction we wish to retrieve the `description` metadata field.

Example:
```julia
reaction = @reaction k, 0 --> X, [description="A reaction"]
getdescription(reaction)
```
"""
function getdescription(reaction::Reaction)
    if hasdescription(reaction)
        return getmetadata(reaction, :description)
    else
        error("Attempts to access `description` metadata field for a reaction which does not have a value assigned for this metadata.")
    end
end

# Misc.
"""
hasmisc(reaction::Reaction)

Returns `true` if the input reaction has the `misc` metadata field assigned, else `false`.

Arguments:
- `reaction`: The reaction we wish to check for the `misc` metadata field.

Example:
```julia
reaction = @reaction k, 0 --> X, [misc="A reaction"]
hasmisc(reaction)
```
"""
function hasmisc(reaction::Reaction)
    return hasmetadata(reaction, :misc)
end

"""
getmisc(reaction::Reaction)

Returns `misc` metadata field for the input reaction.

Arguments:
- `reaction`: The reaction we wish to retrieve the `misc` metadata field.

Example:
```julia
reaction = @reaction k, 0 --> X, [misc="A reaction"]
getmisc(reaction)
```

Notes:
- The `misc` field can contain any valid Julia structure. This mean that Catalyst cannot check it
for symbolci variables that are added here. This means that symbolic variables (e.g. parameters of 
species) that are stored here are not accessible to Catalyst. This can cause troubles when e.g. 
creating a `ReactionSystem` programmatically (in which case any symbolic variables stored in the
`misc` metadata field should also be explicitly provided to the `ReactionSystem` constructor). 

"""
function getmisc(reaction::Reaction)
    if hasmisc(reaction)
        return getmetadata(reaction, :misc)
    else
        error("Attempts to access `misc` metadata field for a reaction which does not have a value assigned for this metadata.")
    end
end


### Units Handling ###

"""
    validate(rx::Reaction; info::String = "")

Check that all substrates and products within the given [`Reaction`](@ref) have
the same units, and that the units of the reaction's rate expression are
internally consistent (i.e. if the rate involves sums, each term in the sum has
the same units).

"""
function validate(rx::Reaction; info::String = "")
    validated = MT._validate([rx.rate], [string(rx, ": rate")], info = info)

    subunits = isempty(rx.substrates) ? nothing : get_unit(rx.substrates[1])
    for i in 2:length(rx.substrates)
        if get_unit(rx.substrates[i]) != subunits
            validated = false
            @warn(string("In ", rx, " the substrates have differing units."))
        end
    end

    produnits = isempty(rx.products) ? nothing : get_unit(rx.products[1])
    for i in 2:length(rx.products)
        if get_unit(rx.products[i]) != produnits
            validated = false
            @warn(string("In ", rx, " the products have differing units."))
        end
    end

    if (subunits !== nothing) && (produnits !== nothing) && (subunits != produnits)
        validated = false
        @warn(string("in ", rx,
                     " the substrate units are not consistent with the product units."))
    end

    validated
end