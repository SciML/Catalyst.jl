### Spatial Reaction Structures ###

# Abstract spatial reaction structures.
abstract type AbstractSpatialReaction end

### Transport Reaction Structures ###

# Implements the transportparameter metadata field.
struct TransportParameter end
Symbolics.option_to_metadata_type(::Val{:transportparameter}) = TransportParameter

istransportparameter(x::Num, args...) = istransportparameter(Symbolics.unwrap(x), args...)
function istransportparameter(x, default = false)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, TransportParameter, default)
end

# A transport reaction. These are simple to hanlde, and should cover most types of spatial reactions.
# Currently only permit constant rates.
struct TransportReaction <: AbstractSpatialReaction
    """The rate function (excluding mass action terms). Currently only constants supported"""
    rate::Num
    """The species that is subject to difusion."""
    species::Num
    """A symbol representation of the species that is subject to difusion."""
    species_sym::Symbol    # Required for identification in certain cases.

    # Creates a diffusion reaction.
    function TransportReaction(rate::Num, species::Num)
        new(rate, species, ModelingToolkit.getname(species))
    end
    function TransportReaction(rate::Number, species::Num)
        new(Num(rate), species, ModelingToolkit.getname(species))
    end
    function TransportReaction(rate::Symbol, species::Num)
        new(Symbolics.variable(rate), species, ModelingToolkit.getname(species))
    end
    function TransportReaction(rate::Num, species::Symbol)
        new(rate, Symbolics.variable(species), species)
    end
    function TransportReaction(rate::Number, species::Symbol)
        new(Num(rate), Symbolics.variable(species), species)
    end
    function TransportReaction(rate::Symbol, species::Symbol)
        new(Symbolics.variable(rate), Symbolics.variable(species), species)
    end
end
# Creates a vector of TransportReactions.
function transport_reactions(transport_reactions)
    [TransportReaction(dr[1], dr[2]) for dr in transport_reactions]
end
# Gets the parameters in a transport reaction.
ModelingToolkit.parameters(dr::TransportReaction) = Symbolics.get_variables(dr.rate)