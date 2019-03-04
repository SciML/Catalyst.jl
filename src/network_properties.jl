""" 
Functions for querying network properties.
"""

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

"""
Return a Vector of pairs, mapping ids of species that serve as substrates 
in the reaction to the corresponding stoichiometric coefficiant as a substrate.
"""
function get_substrate_stoich(rs::ReactionStruct, specmap)
    reactant_stoich = [specmap[s.reactant] => s.stoichiometry for s in rs.substrates]
    sort!(reactant_stoich)
    reactant_stoich
end

"""
Return a Vector of pairs, mapping ids of species that participate in the 
reaction to the net stoichiometric coefficiant of the species (i.e. net change
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

