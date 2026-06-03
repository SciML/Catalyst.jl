### Spatial Reaction Structures ###

# Abstract spatial reaction structures.
abstract type AbstractSpatialReaction end

### Edge Parameter Metadata ###

# Implements the edgeparameter metadata field.
struct EdgeParameter end
Symbolics.option_to_metadata_type(::Val{:edgeparameter}) = EdgeParameter

# Implements the isedgeparameter check function.
"""
    isedgeparameter(p)

Returns `true` if the parameter `p` is an edge parameter (else `false`).
"""
isedgeparameter(x::Num, args...) = isedgeparameter(unwrap(x), args...)
function isedgeparameter(x, default = false)
    if iscall(x) && operation(x) === getindex
        x = first(arguments(x))
    end
    Symbolics.getmetadata(x, EdgeParameter, default)
end

### Transport Reaction Structures ###

# A transport reaction. These are simple to handle, and should cover most types of spatial reactions.
# Only permit constant rates (possibly consisting of several parameters).
struct TransportReaction <: AbstractSpatialReaction
    """The rate function (excluding mass action terms). Currently, only constants supported"""
    rate::Any
    """The species that is subject to diffusion."""
    species::SymbolicT

    # Creates a diffusion reaction.
    function TransportReaction(rate, species)
        if any(!MT.isparameter(var) for var in MT.get_variables(rate))
            error("TransportReaction rate contains variables: $(filter(var -> !MT.isparameter(var), MT.get_variables(rate))). The rate must consist of parameters only.")
        end
        new(rate, species.val)
    end
end

# Macro for creating a `TransportReaction`.
macro transport_reaction(rateex::ExprValues, species::ExprValues)
    make_transport_reaction(striplines(rateex), species)
end
function make_transport_reaction(rateex, species)
    # Handle interpolation of variables
    rateex = esc_dollars!(rateex)
    species = esc_dollars!(species)

    # Parses input expression.
    parameters = ExprValues[]
    find_parameters_in_rate!(parameters, rateex)

    # Checks for input errors.
    forbidden_symbol_check(union([species], parameters))

    # Creates expressions corresponding to actual code from the internal DSL representation.
    sexprs = get_usexpr([species], Dict{Symbol, Expr}())
    pexprs = get_psexpr(parameters, [], Dict{Symbol, Expr}())
    iv = :($(DEFAULT_IV_SYM) = default_t())
    trxexpr = :(TransportReaction($rateex, $species))

    # Appends `edgeparameter` metadata to all declared parameters.
    for idx in 4:2:(2 + 2 * length(parameters))
        insert!(pexprs.args, idx, :([edgeparameter = true]))
    end

    quote
        $pexprs
        $iv
        $sexprs
        $trxexpr
    end
end

# Gets the parameters in a `TransportReaction`.
MT.parameters(tr::TransportReaction) = collect(Symbolics.get_variables(tr.rate))

# Gets the species in a `TransportReaction`.
spatial_species(tr::TransportReaction) = [tr.species]

# Checks that a `TransportReaction` is valid for a given reaction system.
function check_spatial_reaction_validity(rs::ReactionSystem, tr::TransportReaction;
        edge_parameters = [])
    # Checks that the species exist in the reaction system.
    # (ODE simulation code becomes difficult if this is not required,
    # as non-spatial jacobian and f function generated from rs are of the wrong size).
    if !any(isequal(tr.species), species(rs))
        error("Currently, species used in `TransportReaction`s must have previously been declared within the non-spatial ReactionSystem. This is not the case for $(tr.species).")
    end

    # Checks that the rate does not depend on species.
    rate_vars = MT.getname.(Symbolics.get_variables(tr.rate))
    if !isempty(intersect(MT.getname.(species(rs)), rate_vars))
        error("The following species were used in rates of a transport reactions: $(setdiff(MT.getname.(species(rs)), rate_vars)).")
    end

    # Checks that the species does not exist in the system with different metadata.
    if any(isequal(tr.species, s) && !isequivalent(tr.species, s) for s in species(rs))
        error("A transport reaction used a species, $(tr.species), with metadata not matching its discrete space reaction system. Please fetch this species from the reaction system and use it during transport reaction creation.")
    end
    if any(isequal(rs_p, tr_p) && !isequivalent(rs_p, tr_p)
    for rs_p in parameters(rs), tr_p in Symbolics.get_variables(tr.rate))
        error("A transport reaction used a parameter with metadata not matching its discrete space reaction system. Please fetch this parameter from the reaction system and use it during transport reaction creation.")
    end

    # Checks that no edge parameter occurs among rates of non-spatial reactions.
    if any(!isempty(intersect(Symbolics.get_variables(r.rate), edge_parameters))
    for r in reactions(rs))
        error("Edge parameter(s) were found as a rate of a non-spatial reaction.")
    end
end

# Since MTK's "isequal" ignores metadata, we have to use a special function that accounts for this.
# This is important because whether something is an edge parameter is defined in metadata.
# MTK potentially start adding there own metadata to system symbolic variables, to prevent breakage
# for this we now also have al list of `ignored_metadata` (which MTK will add to `ReactionSystem`)
# metadata. Long-term the solution is to permit the creation of spatial reactions within
# a `ReactionSystem` when it is created.
const ep_metadata = Catalyst.EdgeParameter => true
function isequivalent(sym1, sym2; ignored_metadata = [MT.SymScope])
    isequal(sym1, sym2) || (return false)
    if any((md1 != ep_metadata) && (md1[1] ∉ ignored_metadata) && (md1 ∉ sym2.metadata)
    for md1 in sym1.metadata)
        return false
    elseif any((md2 != ep_metadata) && (md2[1] ∉ ignored_metadata) && (md2 ∉ sym1.metadata)
    for md2 in sym2.metadata)
        return false
    elseif typeof(sym1) != typeof(sym2)
        return false
    end
    return true
end

# Implements custom `==`.
"""
    ==(rx1::TransportReaction, rx2::TransportReaction)

Tests whether two [`TransportReaction`](@ref)s are identical.
"""
function (==)(tr1::TransportReaction, tr2::TransportReaction)
    isequal(tr1.rate, tr2.rate) || return false
    isequal(tr1.species, tr2.species) || return false
    return true
end

# Implements custom `hash`.
function hash(tr::TransportReaction, h::UInt)
    h = Base.hash(tr.rate, h)
    Base.hash(tr.species, h)
end

### Utility ###

# Loops through a rate and extracts all parameters.
function find_parameters_in_rate!(parameters, rateex::ExprValues)
    if rateex isa Symbol
        if rateex in [:t, :∅, :Ø, :im, :nothing, CONSERVED_CONSTANT_SYMBOL]
            error("Forbidden term $(rateex) used in transport reaction rate.")
        elseif !(rateex in [:ℯ, :pi, :π])
            push!(parameters, rateex)
        end
    elseif rateex isa Expr
        # Note, this (correctly) skips $(...) expressions.
        for i in 2:length(rateex.args)
            find_parameters_in_rate!(parameters, rateex.args[i])
        end
    end
    return nothing
end
