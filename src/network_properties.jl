""" 
Functions for querying network properties.
"""

######### Accessors: #########

"""
Return a Dictionary mapping from species symbol to species index.
"""
function speciesmap(network)
    network.syms_to_ints
end

"""
Return a Dictionary mapping from parameter symbol to parameter index.
"""
function paramsmap(network)
    network.params_to_ints
end

function numspecies(network)
    length(speciesmap(network))
end

function numreactions(network)
    length(network.reactions)
end

function numparams(network)
    length(paramsmap(network))
end

######### Generated Expressions: #########

"""
Return a vector of the ODE derivative expressions.
"""
function odeexprs(network)
    isnothing(network.f_func) && error("Error, call addodes! first.")
    network.f_func
end

"""
Return a matrix with the ODE Jacobian expressions.
"""
function jacobianexprs(network)
    isnothing(network.symjac) && error("Error, call addodes! first.")
    network.symjac
end

"""
Return a vector of the SDE noise expressions for each reaction.
"""
function noiseexprs(network)
    isnothing(network.g_func) && error("Error, call addsdes! first.")
    network.g_func
end

"""
Return a tuple of the jump rates and affects expressions
"""
function jumpexprs(network)
    (isnothing(network.jump_rate_expr) || isnothing(network.jump_affect_expr)) && error("Error, call addjumps! first.")
    network.jump_rate_expr,network.jump_affect_expr
end


######### Reaction Properties: #########

"""
Return a Vector of pairs, mapping ids of species that serve as substrates 
in the reaction to the corresponding stoichiometric coefficient as a substrate.
"""
function get_substrate_stoich(rs::ReactionStruct, specmap)
    reactant_stoich = [specmap[s.reactant] => s.stoichiometry for s in rs.substrates]
    sort!(reactant_stoich)
    reactant_stoich
end

"""
Return a Vector of pairs, mapping ids of species that participate in the 
reaction to the net stoichiometric coefficient of the species (i.e. net change
in the species due to the reaction).
"""
function get_net_stoich(rs::ReactionStruct, specmap)
    nsdict = Dict{Int,Int}(specmap[s.reactant] => -s.stoichiometry for s in rs.substrates)

    for prod in rs.products
        pidx = specmap[prod.reactant]
        nsdict[pidx] = haskey(nsdict, pidx) ? nsdict[pidx] + prod.stoichiometry : prod.stoichiometry
    end

    net_stoich = Vector{Pair{Int,Int}}()
    for stoich_map in sort(collect(nsdict))
        (stoich_map[2] != zero(Int)) && push!(net_stoich, stoich_map)
    end

    net_stoich
end

######### Network Properties: #########

"""
Returns a Vector{Vector{Int}} mapping a reaction index 
to the indices of species that depend on it.
"""
function rxtospecies_depgraph(network)
    specmap = speciesmap(network)

    # map from a reaction to vector of species that depend on it
    [sort!( [ns.first for ns in get_net_stoich(rx, specmap)] ) for rx in network.reactions]
end

"""
Returns a Vector{Vector{Int}} mapping a species index 
to the indices of reactions that depend on it.
"""
function speciestorx_depgraph(network)
    numrxs  = numreactions(network)    
    specmap = speciesmap(network)

    # map from a species to reactions that depend on it
    spectorxs = [Vector{Int}() for i=1:numspecies(network)]
    for rx in 1:numrxs
        for specsym in network.reactions[rx].dependants
            push!(spectorxs[specmap[specsym]], rx)
        end
    end
    foreach(stor -> sort!(stor), spectorxs)

    spectorxs
end

"""
Returns a Vector{Vector{Int}} mapping a reaction index
to the indices of reactions that depend on it.
"""
function rxtorx_depgraph(network, sptorxsmap=nothing)
    sptorxs = isnothing(sptorxsmap) ? speciestorx_depgraph(network) : sptorxsmap
    numrxs = numreactions(network)
    dep_sets = [SortedSet{Int}() for n = 1:numrxs]

    # dg as vector of sets
    for rx in 1:numrxs
        net_stoich = get_net_stoich(network.reactions[rx], speciesmap(network))

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