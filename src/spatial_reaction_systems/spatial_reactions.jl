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

# A transport reaction. These are simple to handle, and should cover most types of spatial reactions.
# Only permit constant rates (possibly consisting of several parameters).
struct TransportReaction <: AbstractSpatialReaction
    """The rate function (excluding mass action terms). Currently only constants supported"""
    rate::Any
    """The species that is subject to diffusion."""
    species::BasicSymbolic{Real}

    # Creates a diffusion reaction.
    function TransportReaction(rate, species)
        if any(!ModelingToolkit.isparameter(var) for var in ModelingToolkit.get_variables(rate)) 
            error("TransportReaction rate contains variables: $(filter(var -> !ModelingToolkit.isparameter(var), ModelingToolkit.get_variables(rate))). The rate must consist of parameters only.")
        end
        new(rate, species.val)
    end
end
# Creates a vector of TransportReactions.
function TransportReactions(transport_reactions)
    [TransportReaction(tr[1], tr[2]) for tr in transport_reactions]
end

# Macro for creating a transport reaction.
macro transport_reaction(rateex::ExprValues, species::ExprValues)
    make_transport_reaction(MacroTools.striplines(rateex), species)
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
    sexprs = get_sexpr([species], Dict{Symbol, Expr}())
    pexprs = get_pexpr(parameters, Dict{Symbol, Expr}())
    iv = :(@variables $(DEFAULT_IV_SYM))
    trxexpr = :(TransportReaction($rateex, $species))

    quote
        $pexprs
        $iv
        $sexprs
        $trxexpr
    end
end

# Gets the parameters in a transport reaction.
ModelingToolkit.parameters(tr::TransportReaction) = Symbolics.get_variables(tr.rate)

# Gets the species in a transport reaction.
spatial_species(tr::TransportReaction) = [tr.species]

# Checks that a transport reaction is valid for a given reaction system.
function check_spatial_reaction_validity(rs::ReactionSystem, tr::TransportReaction; edge_parameters=[])
    # Checks that the species exist in the reaction system.
    # (ODE simulation code becomes difficult if this is not required,
    # as non-spatial jacobian and f function generated from rs is of wrong size).  
    if !any(isequal(tr.species), species(rs)) 
        error("Currently, species used in TransportReactions must have previously been declared within the non-spatial ReactionSystem. This is not the case for $(tr.species).")
    end

    # Checks that the rate does not depend on species.    
    rate_vars = ModelingToolkit.getname.(Symbolics.get_variables(tr.rate))
    if !isempty(intersect(ModelingToolkit.getname.(species(rs)), rate_vars)) 
        error("The following species were used in rates of a transport reactions: $(setdiff(ModelingToolkit.getname.(species(rs)), rate_vars)).")
    end

    # Checks that the species does not exist in the system with different metadata.
    if any([isequal(tr.species, s) && !isequivalent(tr.species, s) for s in species(rs)]) 
        error("A transport reaction used a species, $(tr.species), with metadata not matching its lattice reaction system. Please fetch this species from the reaction system and used in transport reaction creation.")
    end
    if any([isequal(rs_p, tr_p) && !equivalent_metadata(rs_p, tr_p) for rs_p in parameters(rs), tr_p in Symbolics.get_variables(tr.rate)]) 
        error("A transport reaction used a parameter with metadata not matching its lattice reaction system. Please fetch this parameter from the reaction system and used in transport reaction creation.")
    end

    # Checks that no edge parameter occur among rates of non-spatial reactions.
    if any([!isempty(intersect(Symbolics.get_variables(r.rate), edge_parameters)) for r in reactions(rs)])
        error("Edge paramter(s) were found as a rate of a non-spatial reaction.")
    end
end
equivalent_metadata(p1, p2) = isempty(setdiff(p1.metadata, p2.metadata, [Catalyst.EdgeParameter => true]))

# Since MTK's "isequal" ignores metadata, we have to use a special function that accounts for this.
# This is important because whether something is an edge parameter is defined in metadata.
function isequivalent(sym1, sym2)
    !isequal(sym1, sym2) && (return false)
    (sym1.metadata != sym2.metadata) && (return false)
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
# Loops through a rate and extract all parameters.
function find_parameters_in_rate!(parameters, rateex::ExprValues)
    if rateex isa Symbol
        if rateex in [:t, :∅, :im, :nothing, CONSERVED_CONSTANT_SYMBOL]
            error("Forbidden term $(rateex) used in transport reaction rate.")
        elseif !(rateex in [:ℯ, :pi, :π])
            push!(parameters, rateex)
        end
    elseif rateex isa Expr
        # Note, this (correctly) skips $(...) expressions
        for i in 2:length(rateex.args)
            find_parameters_in_rate!(parameters, rateex.args[i])
        end
    end
    nothing
end