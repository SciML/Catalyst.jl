# Functions for querying network properties.

######### Accessors: #########

"""
    species(network)
Given an `AbstractReactionNetwork`, return a vector of species symbols.
"""
function species(network)
    network.syms
end

"""
    params(network)
Given an `AbstractReactionNetwork`, return a vector of parameter symbols.
"""
function params(network)
    network.params
end

"""
    speciesmap(network)

Given an `AbstractReactionNetwork`, return a Dictionary mapping from species
symbol to species index. 
"""
function speciesmap(network)
    network.syms_to_ints
end

"""
    paramsmap(network)

Given an `AbstractReactionNetwork`, return a Dictionary mapping from parameter
symbol to parameter index.
"""
function paramsmap(network)
    network.params_to_ints
end

"""
    numspecies(network)

Return the number of species within the given `AbstractReactionNetwork`.
"""
function numspecies(network)
    length(speciesmap(network))
end

"""
    numreactions(network)

Return the number of reactions within the given `AbstractReactionNetwork`.
"""
function numreactions(network)
    length(network.reactions)
end

"""
    numparams(network)

Return the number of parameters within the given `AbstractReactionNetwork`.
"""
function numparams(network)
    length(paramsmap(network))
end

######### Generated Functions: #########

"""
    oderhsfun(network)

Given an `AbstractReactionNetwork`, return a function, `f!(du,u,p,t)`, that
evaluates the current value of the ODE model derivative functions, ``du/dt = f(u,t)``, 
within `du`.

*Note,* for a network generated with the `@min_reaction_network` macro `addodes!`
must be called first.
"""
function oderhsfun(network)
    isnothing(network.f) && error("Error, call addodes! first.")
    network.f
end

"""
    jacfun(network)

Given an `AbstractReactionNetwork`, return a function, `jac!(J,u,p,t)`, that
evaluates the current Jacobian matrix, `J`, of the ODE model, ``du/dt = f(u,t)``. 
The Jacobian matrix has entries 

``J_{i,j} = \\partial f_i(u,t) / \\partial u_j``.

*Note,* for a network generated with the `@min_reaction_network` macro `addodes!` 
must be called first.
"""
function jacfun(network)
    isnothing(network.jac) && error("Error, call addodes! first.")
    network.jac
end

"""
    paramjacfun(network)

Given an `AbstractReactionNetwork`, return a function, `pjac(pJ,u,p,t)`, that
evaluates the current parameter Jacobian matrix, `pJ`, of the ODE model, ``du/dt = f(u,t)``. 
The parameter Jacobian matrix has entries 

``pJ_{i,j} = \\partial f_i(u,t) / \\partial p_j``.

*Note,* for a network generated with the `@min_reaction_network` macro `addodes!` 
must be called first.
"""
function paramjacfun(network)
    isnothing(network.paramjac) && error("Error, call addodes! first.")
    network.paramjac
end

"""
    odefun(network)

Given an `AbstractReactionNetwork`, return a `DiffEqBase.ODEFunction` encoding
an ODE model for the reaction network. 

*Note,* for a network generated with the `@min_reaction_network` macro `addodes!` 
must be called first.
"""
function odefun(network)
    isnothing(network.odefun) && error("Error, call addodes! first.")
    network.odefun
end

"""
    noisefun(network)

Given an `AbstractReactionNetwork`, return a function, `g(η,u,p,t)`, that
evaluates the current noise coefficients for each reaction in the Chemical
Langevin Equation representation within `η`.

*Note,* for a network generated with the `@min_reaction_network` macro `addsdes!` 
must be called first.
"""
function noisefun(network)
    isnothing(network.g) && error("Error, call addsdes! first.")
    network.g 
end

"""
    sdefun(network)

Given an `AbstractReactionNetwork`, return a `DiffEqBase.SDEFunction` encoding
a Chemical Langevin Equation SDE model for the reaction network. 

*Note,* for a network generated with the `@min_reaction_network` macro `addsdes!` 
must be called first.
"""
function sdefun(network)
    isnothing(network.sdefun) && error("Error, call addsdes! first.")
    network.sdefun
end

"""
    jumps(network)

Given an `AbstractReactionNetwork`, return a tuple of `AbstractJumps` encoding
a stochastical chemical kinetics representation for the reaction network.

*Note,* for a network generated with the `@min_reaction_network` macro `addjumps!` 
must be called first.
"""
function jumps(network)
    isnothing(network.jumps) && error("Error, call addjumps! first.")
    network.jumps
end

"""
    regularjumps(network)

Given an `AbstractReactionNetwork`, return a `RegularJump` encoding a
stochastical chemical kinetics representation of the reaction network for use in
``\\tau``-leaping approximations.

*Note,* for a network generated with the `@min_reaction_network` macro `addjumps!` 
must be called first.
"""
function regularjumps(network)
    isnothing(network.regular_jumps) && error("Error, call addjumps! first.")
    network.regular_jumps
end


######### Generated Expressions: #########

"""
    odeexprs(network)

Given an `AbstractReactionNetwork`, return a vector of the ODE expressions.

*Note,* for a network generated with the `@min_reaction_network` macro `addodes!` 
must be called first.
"""
function odeexprs(network)
    isnothing(network.f_func) && error("Error, call addodes! first.")
    network.f_func
end

"""
    jacobianexprs(network)

Given an `AbstractReactionNetwork`, return a matrix with the ODE Jacobian
expressions.

*Note,* for a network generated with the `@min_reaction_network` macro `addodes!` 
must be called first.
"""
function jacobianexprs(network)
    isnothing(network.symjac) && error("Error, call addodes! first.")
    network.symjac
end

"""
    noiseexprs(network)

Given an `AbstractReactionNetwork`, return a vector of the SDE noise expressions
for each reaction.

*Note,* for a network generated with the `@min_reaction_network` macro `addsdes!` 
must be called first.
"""
function noiseexprs(network)
    isnothing(network.g_func) && error("Error, call addsdes! first.")
    network.g_func
end

"""
    jumpexprs(network)

Given an `AbstractReactionNetwork`, return a tuple of the jump rates and affects
expressions.

*Note,* for a network generated with the `@min_reaction_network` macro `addjumps!` 
must be called first.
"""
function jumpexprs(network)
    (isnothing(network.jump_rate_expr) || isnothing(network.jump_affect_expr)) && error("Error, call addjumps! first.")
    network.jump_rate_expr,network.jump_affect_expr
end

"""
    rateexpr(network, rxidx)

Given an `AbstractReactionNetwork`, return the reaction rate expression for the
reaction with index `rxidx`. Note, for a reaction defined by

`k*X*Y, X+Z --> 2X + Y`

the expression that is returned will be `:(k*X*Y)`, while the *rate law* used in
ODEs and SDEs would be `k*X^2*Y*Z`.
"""
function rateexpr(network, rxidx)
    network.reactions[rxidx].rate_org
end

"""
    oderatelawexpr(network, rxidx)

Given an `AbstractReactionNetwork`, return the reaction rate law expression used
in generated ODEs for the reaction with index `rxidx`. Note, for a reaction
defined by

`k*X*Y, X+Z --> 2X + Y`

the expression that is returned will be `:(k*X^2*Y*Z)`. For a reaction of the
form 

`k, 2X+3Y --> Z`

the expression that is returned will be `:(k * (X^2/2) * (Y^3/6))`.
""" 
function oderatelawexpr(network, rxidx)
    network.reactions[rxidx].rate_DE
end

"""
    ssaratelawexpr(network, rxidx)

Given an `AbstractReactionNetwork`, return the reaction rate law expression used
in generated stochastic chemical kinetic model SSAs for the reaction with index
`rxidx`. Note, for a reaction defined by

`k*X*Y, X+Z --> 2X + Y`

the expression that is returned will be `:(k*X^2*Y*Z)`. For a reaction of the
form 

`k, 2X+3Y --> Z`

the expression that is returned will be `:(k * binomial(X,2) * binomial(Y,3))`.
""" 
function ssaratelawexpr(network, rxidx)
    network.reactions[rxidx].rate_SSA
end

######### Reaction Properties: #########

function substratestoich(rs::ReactionStruct, specmap)
    sort!( [specmap[s.reactant] => s.stoichiometry for s in rs.substrates] )
end

"""
    substratestoich(network, rxidx)

Given an `AbstractReactionNetwork` and a reaction index, `rxidx`, return a
vector of pairs, mapping ids of species that serve as substrates in the reaction
to the corresponding stoichiometric coefficient as a substrate. 

Allocates a new vector to store the pairs.
""" 
function substratestoich(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    substratestoich(rn.reactions[rxidx], speciesmap(rn))
end

"""
    substratesymstoich(network, rxidx)

Given an `AbstractReactionNetwork` and a reaction index, `rxidx`, return a
Vector of `ReactantStruct`s, mapping the symbols of species that serve as
substrates in the reaction to the corresponding stoichiometric coefficient as a
substrate. 

Non-allocating, returns underlying field within the reaction_network.
""" 
function substratesymstoich(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    rn.reactions[rxidx].substrates
end

function productstoich(rs::ReactionStruct, specmap)
    sort( [specmap[p.reactant] => p.stoichiometry for p in rs.products] )    
end

"""
    productstoich(network, rxidx)

Given an `AbstractReactionNetwork` and a reaction index, `rxidx`, return a
vector of pairs, mapping ids of species that are products in the reaction to the
corresponding stoichiometric coefficient as a product.

Allocates a new vector to store the pairs.
"""
function productstoich(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    productstoich(rn.reactions[rxidx], speciesmap(rn))
end

"""
    productsymstoich(network, rxidx)

Given an `AbstractReactionNetwork` and a reaction index, `rxidx`, return a
Vector of `ReactantStruct`s, mapping the symbols of species that are products in
the reaction to the corresponding stoichiometric coefficient as a product. 

Non-allocating, returns underlying field within the reaction_network.
""" 
function productsymstoich(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    rn.reactions[rxidx].products
end

function netstoich(rs::ReactionStruct, specmap)
    sort!([specmap[ns.reactant] => ns.stoichiometry for ns in rs.netstoich])
end

"""
    netstoich(network, rxidx)

Given an `AbstractReactionNetwork` and a reaction index, `rxidx`, return a
vector of pairs, mapping ids of species that change numbers due to the reaction
to the net stoichiometric coefficient of the species (i.e. net change in the
species due to the reaction).

Allocates a new vector to store the pairs.
"""
function netstoich(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    netstoich(rn.reactions[rxidx], speciesmap(rn))
end

"""
    netsymstoich(network, rxidx)

Given an `AbstractReactionNetwork` and a reaction index, `rxidx`, return a
Vector of `ReactantStruct`s, mapping the symbols of species that change numbers
due to the reaction to the net stoichiometric coefficient of the species (i.e.
net change in the species due to the reaction).

Non-allocating, returns underlying field within the reaction_network.
"""
function netsymstoich(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    rn.reactions[rxidx].netstoich
end

"""
    ismassaction(network, rxidx)

Given an `AbstractReactionNetwork` and a reaction index, `rxidx`, return a
boolean indicating whether the given reaction is of mass action form. For
example, the reaction

`2*k, 2X + 3Y --> 5Z + W`

would return true, while reactions with state-dependent rates like

`k*X, X + Y --> Z`

would return false.

Non-allocating, returns underlying field within the reaction_network.
"""
function ismassaction(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    rn.reactions[rxidx].is_pure_mass_action
end

"""
    dependents(network, rxidx)

Given an `AbstractReactionNetwork` and a reaction index, `rxidx`, return a
vector of symbols of species the *reaction rate law* depends on. i.e. for

`k*W, 2X + 3Y --> 5Z + W`

the returned vector would be `[:W,:X,:Y]`.

Non-allocating, returns underlying field within the reaction_network.
"""
function dependents(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    rn.reactions[rxidx].dependants
end

"""
    dependants(network, rxidx)

See documentation for [`dependents(network, rxidx)`](@ref).
"""
function dependants(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    dependents(rn, rxidx)
end

"""
    substrates(network, rxidx)

Given an `AbstractReactionNetwork` and a reaction index, `rxidx`, return a
vector of symbols of species that correspond to substrates in the reaction. 
i.e. for

`k*W, X + 3Y --> X + W`

the returned vector would be `[:X,:Y]`.

Allocates a new vector to store the symbols.
"""
function substrates(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    [sub.reactant for sub in rn.reactions[rxidx].substrates]
end

"""
    products(network, rxidx)

Given an `AbstractReactionNetwork` and a reaction index, `rxidx`, return a
vector of symbols of species that correspond to products in the reaction. 
i.e. for

`k*W, X + 3Y --> X + W`

the returned vector would be `[:X,:W]`.

Allocates a new vector to store the symbols.
"""
function products(rn::DiffEqBase.AbstractReactionNetwork, rxidx)
    [prod.reactant for prod in rn.reactions[rxidx].products]
end

######### Network Properties: #########

"""
    rxtospecies_depgraph(network)

Given an `AbstractReactionNetwork`, returns a `Vector{Vector{Int}}` mapping each
reaction index to the indices of species that depend on it.
"""
function rxtospecies_depgraph(network)
    specmap = speciesmap(network)

    # map from a reaction to vector of species that depend on it
    [sort!( [ns.first for ns in netstoich(rx, specmap)] ) for rx in network.reactions]
end

"""
    speciestorx_depgraph(network)

Given an `AbstractReactionNetwork`, returns a `Vector{Vector{Int}}` mapping each
species index to the indices of reactions that depend on it.
"""
function speciestorx_depgraph(network)
    numrxs  = numreactions(network)    
    specmap = speciesmap(network)

    # map from a species to reactions that depend on it
    spectorxs = [Vector{Int}() for i=1:numspecies(network)]
    for rx in 1:numrxs
        for specsym in dependents(network, rx) 
            push!(spectorxs[specmap[specsym]], rx)
        end
    end
    foreach(stor -> sort!(stor), spectorxs)

    spectorxs
end

"""
    rxtorx_depgraph(network)

Given an `AbstractReactionNetwork`, returns a `Vector{Vector{Int}}` mapping each
reaction index to the indices of reactions that depend on it.
"""
function rxtorx_depgraph(network, sptorxs=speciestorx_depgraph(network))
    numrxs = numreactions(network)
    dep_sets = [SortedSet{Int}() for n = 1:numrxs]

    # dg as vector of sets
    for rx in 1:numrxs
        net_stoich = netstoich(network, rx)

        for (spec,stoich) in net_stoich
            for dependent_rx in sptorxs[spec]
                push!(dep_sets[rx], dependent_rx)
            end
        end

        # every reaction should depend on itself
        push!(dep_sets[rx], rx)
    end

    # dg as vector of vectors
    dep_graph = Vector{Vector{Int}}(undef, numrxs)
    for rx in 1:numrxs
        dep_graph[rx] = [dep for dep in dep_sets[rx]]
    end

    dep_graph
end