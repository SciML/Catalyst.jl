# Functions for querying network properties.

######### Accessors: #########

"""
    species(network)

Given an `ReactionSystem`, return a vector of species `Variable`s.
"""
function species(network)
    states(network)
end

"""
    params(network)

Given an `ReactionSystem`, return a vector of parameter `Variable`s.
"""
function params(network)
    parameters(network)
end

"""
    speciesmap(network)

Given an `ReactionSystem`, return a Dictionary mapping from species
`Variable`s to species indices. (Allocates)
"""
function speciesmap(network)
    Dict(S => i for (i,S) in enumerate(species(network)))
end

"""
    paramsmap(network)

Given an `ReactionSystem`, return a Dictionary mapping from parameter
`Variable`s to parameter indices. (Allocates)
"""
function paramsmap(network)
    Dict(p => i for (i,p) in enumerate(params(network)))
end

"""
    numspecies(network)

Return the number of species within the given `ReactionSystem`.
"""
function numspecies(network)
    length(species(network))
end

"""
    numreactions(network)

Return the number of reactions within the given `ReactionSystem`.
"""
function numreactions(network)
    length(equations(network))
end

"""
    numparams(network)

Return the number of parameters within the given `ReactionSystem`.
"""
function numparams(network)
    length(params(network))
end


### Macros for using DSL notation to modify an already created network. ###

# # Creates a new network by making addition to an already existing one.
# macro add_reactions(network, ex::Expr, parameters...)
#     #To be written
# end

# # modifies an already existing reaction network by adding new reactions.
# macro add_reactions!(network, ex::Expr, parameters...)
#     #To be written
# end


######################## reaction network operators #######################

"""
    ==(rn1::ModelingToolki.ReactionSystem, rn2::ModelingToolkit.ReactionSystem)

Tests whether the underlying species `Variables`s, parameter `Variables`s and reactions
are the same in the two networks. Ignores order network components were defined. 
*Does not* currently simplify rates, so a rate of `A^2+2*A+1` would be considered
different than `(A+1)^2`.
"""
function (==)(rn1::ReactionSystem, rn2::ReactionSystem)
    issetequal(species(rn1), species(rn2)) || return false
    issetequal(params(rn1), params(rn2)) || return false
    isequal(rn1.iv, rn2.iv) || return false
    (numreactions(rn1) == numreactions(rn2)) || return false
    for eq1 in equations(rn1), eq2 in equations(rn2)
        isequal(eq1, eq2) || return false
    end
    for sys1 in rn1.systems, sys2 in rn2.systems
        (sy1 == sys2) || return false
    end
end

######################## functions to extend a network ####################

"""
    addspecies!(network, s::Variable)

Given a `ReactionSystem`, add the species corresponding to the variable `s`
to the network (if it is not already defined). Returns the integer id 
of the species within the system.
"""
function addspecies!(rn::ReactionSystem, s::Variable)
    curidx = findfirst(S -> isequal(S, s), species(rn))
    if curidx === nothing
        push!(rn.states, s)
        return length(species(rn))
    else
        return curidx
    end    
end

"""
    addspecies!(network, speciesop::Operation)

Given a `ReactionSystem`, add the species corresponding to the variable `s`
to the network (if it is not already defined). Returns the integer id 
of the species within the system.
"""
function addspecies!(rn::ReactionSystem, s::Operation) 
    !(s.op isa Variable) && error("If the passed in species is an Operation, it must correspond to an underlying Variable.")        
    addspecies!(rn, convert(Variable,s))    
end

"""
    addparam!(network, p::Variable)

Given a `ReactionSystem`, add the parameter corresponding to the variable `p`
to the network (if it is not already defined). Returns the integer id 
of the parameter within the system.
"""
function addparam!(rn::ReactionSystem, p::Variable)
    curidx = findfirst(S -> isequal(S, p), params(rn))
    if curidx === nothing
        push!(rn.ps, p)
        return length(params(rn))
    else
        return curidx
    end    
end

"""
    addparams!(network, p::Operation)

Given a `ReactionSystem`, add the parameter corresponding to the variable `p`
to the network (if it is not already defined). Returns the integer id 
of the parameter within the system.
"""
function addparam!(rn::ReactionSystem, p::Operation) 
    !(p.op isa Variable) && error("If the passed in parameter is an Operation, it must correspond to an underlying Variable.")        
    addparam!(rn, convert(Variable,p))    
end

# """
#     addreaction!(network, rateex::Union{Expr,Symbol,Int,Float64}, rxexpr::Expr)
# Given an AbstractReaction network, add a reaction with the passed in rate and
# reaction expressions. i.e. a reaction of the form
# ```julia
# k*X, 2X + Y --> 2W
# ```
# would have `rateex=:(k*X)` and `rxexpr=:(2X + Y --> W)`,
# ```julia
# 10.5, 0 --> X
# ```
# would have `rateex=10.5` and `rxexpr=:(0 --> X)`, and
# ```julia
# k, X+X --> Z
# ```
# would have `rateex=:k` and `rxexpr=:(X+X --> Z)`.
# All normal DSL reaction definition notation should be supported.
# """
# function addreaction!(rn::DiffEqBase.ReactionSystem, rateex::ExprValues, rxexpr::Expr)
#     ex = Expr(:block, :(($rateex, $rxexpr)))
#     newrxs = get_reactions(ex)
#     foreach(rx -> push!(rn.reactions,ReactionStruct(rx, species(rn))), newrxs)
#     nothing
# end

# """
#     addreaction!(network, rateex::Union{Expr,Symbol,Int,Float64}, substrates, products)
# Given an AbstractReaction network, add a reaction with the passed in rate,
# `rateex`, substrate stoichiometry, and product stoichiometry. Stoichiometries
# are represented as tuples of `Pair{Symbol,Int}`. i.e. a reaction of the form
# ```julia
# k*X, 2X + Y --> 2W
# ```
# would have `rateex=:(k*X)`, `substrates=(:X=>2, :Y=>2)`` and
# `products=(W=>2,)`,
# ```julia
# 10.5, 0 --> X
# ```
# would have `rateex=10.5`, `substrates=()` and `products=(:X=>1,)`, and
# ```julia
# k, X+X --> Z
# ```
# would have `rateex=:k`, `substrates=(:X=>2,)` and `products=(:Z=>2,)`.
# All normal DSL reaction definition notation should be supported for the
# `rateex`.
# """
# function addreaction!(rn::DiffEqBase.ReactionSystem, rateex::ExprValues,
#                                         subs::Tuple{Vararg{Pair{Symbol,Int}}},
#                                         prods::Tuple{Vararg{Pair{Symbol,Int}}}) where {T <: Number}

#     substrates = ReactantStruct[ReactantStruct(p[1],p[2]) for p in subs]
#     dependents = Symbol[p[1] for p in subs]
#     products = ReactantStruct[ReactantStruct(p[1],p[2]) for p in prods]
#     ns = netstoich(substrates, products)
#     rate_DE = isempty(subs) ? rateex : mass_rate_DE(substrates, true, rateex)
#     rate_SSA = isempty(subs) ? rateex : mass_rate_SSA(substrates, true, rateex)

#     # resolve dependents from rateex
#     if rateex isa Number
#         ismassaction = true
#     elseif rateex isa Symbol
#         if haskey(speciesmap(rn), rateex)
#             ismassaction = false
#             if rateex ∉ dependents
#                 push!(dependents, rateex)
#             end
#         elseif haskey(paramsmap(rn), rateex)
#             ismassaction = true
#         else
#             error("rateex is a symbol that is neither a species or parameter.")
#         end
#     else # isa Expr

#         # mimicing ReactionStruct constructor for now, but this should be optimized...
#         newdeps = unique!(recursive_content(rate_DE, speciesmap(rn), Vector{Symbol}()))
#         ismassaction = issetequal(dependents,newdeps)
#         dependents = newdeps
#     end

#     push!(rn.reactions, ReactionStruct(substrates, products, ns, rateex, rate_DE, rate_SSA, dependents, ismassaction))
#     nothing
# end

