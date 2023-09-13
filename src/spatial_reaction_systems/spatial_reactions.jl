### Spatial Reaction Structures ###

# Abstract spatial reaction structures.
abstract type AbstractSpatialReaction end

### EdgeParameter Metadata ###

# Implements the edgeparameter metadata field.
struct EdgeParameter end
Symbolics.option_to_metadata_type(::Val{:edgeparameter}) = EdgeParameter

# Implements the isedgeparameter check function.
isedgeparameter(x::Num, args...) = isedgeparameter(Symbolics.unwrap(x), args...)
function isedgeparameter(x, default = false)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, EdgeParameter, default)
end

### Transport Reaction Structures ###

# A transport reaction. These are simple to hanlde, and should cover most types of spatial reactions.
# Only permit constant rates (possibly consisting of several parameters).
struct TransportReaction <: AbstractSpatialReaction
    """The rate function (excluding mass action terms). Currently only constants supported"""
    rate::Num
    """The species that is subject to difusion."""
    species::BasicSymbolic{Real}

    # Creates a diffusion reaction.
    function TransportReaction(rate::Num, species::Num)
        new(rate, species.val)
    end
end
# If the rate is a value, covert that to a numemric expression.
function TransportReaction(rate::Number, args...)
    TransportReaction(Num(rate), args...)
end
# Creates a vector of TransportReactions.
function transport_reactions(transport_reactions)
    [TransportReaction(tr[1], tr[2]) for tr in transport_reactions]
end

# Macro for creating a transport reaction.
macro transport_reaction(rateex::ExprValues, species::Symbol)
    make_transport_reaction(MacroTools.striplines(rateex), species)
end
function make_transport_reaction(rateex, species)
    parameters = [];
    find_parameters_in_rate!(parameters, rateex)
    quote
        @parameters $(parameters...)
        @variables t
        @species $(species)(t)
        TransportReaction($rateex, $species)
    end
end

# Gets the parameters in a transport reaction.
ModelingToolkit.parameters(tr::TransportReaction) = convert(Vector{BasicSymbolic{Real}}, Symbolics.get_variables(tr.rate))

# Gets the species in a transport reaction.
spatial_species(tr::TransportReaction) = [tr.species]

# Checks that a transport reaction is valid for a given reaction system.
function check_spatial_reaction_validity(rs::ReactionSystem, tr::TransportReaction; edge_parameters=[])
    # Checks that the species exist in the reaction system (ODE simulation code becomes difficult if this is not required, as non-spatial jacobian and f function generated from rs is of wrong size).  
    any(isequal(tr.species), species(rs)) || error("Currently, species used in TransportReactions must also be in non-spatial ReactionSystem. This is not the case for $(tr.species).")
    # Checks that the rate does not depend on species.    
    isempty(intersect(ModelingToolkit.getname.(species(rs)), ModelingToolkit.getname.(Symbolics.get_variables(tr.rate)))) || error("The following species were used in rates of a transport reactions: $(setdiff(ModelingToolkit.getname.(species(rs)), ModelingToolkit.getname.(Symbolics.get_variables(tr.rate)))).")
    # Checks that the species does not exist in the system with different metadata.
    any([isequal(tr.species, s) && !isequal(tr.species.metadata, s.metadata) for s in species(rs)]) && error("A transport reaction used a species, $(tr.species), with metadata not matching its lattice reaction system. Please fetch this species from the reaction system and used in transport reaction creation.")
    any([isequal(rs_p, tr_p) && !equivalent_metadata(rs_p, tr_p) for rs_p in parameters(rs), tr_p in Symbolics.get_variables(tr.rate)]) && error("A transport reaction used a parameter with metadata not matching its lattice reaction system. Please fetch this parameter from the reaction system and used in transport reaction creation.")
    # Checks that no edge parameter occur among rates of non-spatial reactions.
    any([!isempty(intersect(Symbolics.get_variables(r.rate), edge_parameters)) for r in reactions(rs)]) && error("Edge paramter(s) were found as a rate of a non-spatial reaction.")
end
equivalent_metadata(p1, p2) = isempty(setdiff(p1.metadata, p2.metadata, [Catalyst.EdgeParameter => true]))

### Utility ###
# Loops through a rate and extract all parameters.
function find_parameters_in_rate!(parameters, rateex::ExprValues)
    if rateex isa Symbol
        if !(rateex in [:ℯ, :pi, :π])
            push!(parameters, rateex)
        elseif rateex in [:t, :∅, forbidden_symbols_error...]
            error("Forbidden term $(rateex) used in transport reaction rate.")
        end
    elseif rateex isa Expr
        # note, this (correctly) skips $(...) expressions
        for i in 2:length(rateex.args)
            find_parameters_in_rate!(parameters, rateex.args[i])
        end
    end
    nothing
end