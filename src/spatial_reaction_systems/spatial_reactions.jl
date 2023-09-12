### Spatial Reaction Structures ###

# Abstract spatial reaction structures.
abstract type AbstractSpatialReaction end

# Implements the edgeparameter metadata field.
struct EdgeParameter end
Symbolics.option_to_metadata_type(::Val{:edgeparameter}) = EdgeParameter

isedgeparameter(x::Num, args...) = isedgeparameter(Symbolics.unwrap(x), args...)
function isedgeparameter(x, default = false)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, EdgeParameter, default)
end

### Transport Reaction Structures ###

# A transport reaction. These are simple to hanlde, and should cover most types of spatial reactions.
# Currently only permit constant rates.
struct TransportReaction <: AbstractSpatialReaction
    """The rate function (excluding mass action terms). Currently only constants supported"""
    rate::Num
    """The species that is subject to difusion."""
    species::BasicSymbolic{Real}
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

# Macro for creating a transport reaction.
macro transport_reaction(species::Symbol, rate::Expr)
    make_transport_reaction(species, MacroTools.striplines(rate))
end
function make_transport_reaction(species, rate)
    quote
        @parameters 
        @variable t
        @species $(species)(t)
        return TransportReaction(species, rate)
    end
end

# Gets the parameters in a transport reaction.
ModelingToolkit.parameters(tr::TransportReaction) = Symbolics.get_variables(tr.rate)

# Gets the species in a transport reaction.
spatial_species(tr::TransportReaction) = [tr.species]

# Checks that a trnasport reaction is valid for a given reaction system.
function check_spatial_reaction_validity(rs::ReactionSystem, tr::TransportReaction)
    # Checks that the rate does not depend on species.    
    isempty(setdiff(ModelingToolkit.getname(species(rs)), ModelingToolkit.getname(Symbolics.get_variables(tr.rate)))) || error("The following species were used in rates of a transport reactions: $(setdiff(ModelingToolkit.getname(species(rs)), ModelingToolkit.getname(Symbolics.get_variables(tr.rate)))).")
    # Checks that the species does not exist in the system with different metadata.
    any([isequal(tr.species, s) && !isequal(tr.species.metadata, s.metadata)])  || error("A transport reaction used a species ($(tr.species)) with metadata not matching its lattice reaction system. Please fetch this species from the reaction system and used in transport reaction creation.")
end
