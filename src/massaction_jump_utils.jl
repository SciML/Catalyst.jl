# collection of utility functions for getting mass action jumps out of 
# DiffEqBiological reaction network structures

# return a map from species symbol to species index
function species_to_indices(network)
    specs = network.syms
    Dict( specs[i] => i for i in eachindex(specs) )    
end

# return a map from reaction param symbol to rate index
function rate_to_indices(network)
    rates = network.params
    Dict( rates[i] => i for i in eachindex(rates) )
end

# given a ReactionStruct and a species map construct a MassActionJump
function make_majump(rs, specmap, ratemap, params)
    reactant_stoich = Vector{Pair{Int,Int}}()
    nsdict = Dict{Int,Int}()
    for substrate in rs.substrates
        specpair = specmap[substrate.reactant] => substrate.stoichiometry
        push!(reactant_stoich, specpair)
        specpair = specmap[substrate.reactant] => -substrate.stoichiometry
        push!(nsdict, specpair)
    end
    sort!(reactant_stoich)

    for product in rs.products
        prodidx = specmap[product.reactant]
        if haskey(nsdict, prodidx)
            nsdict[prodidx] += product.stoichiometry
        else
            push!(nsdict, prodidx => product.stoichiometry)
        end
    end

    net_stoich = Vector{Pair{Int,Int}}()
    for stoich_map in sort(collect(nsdict))
        push!(net_stoich, stoich_map)
    end

    if isempty(net_stoich)
        error("Empty net stoichiometry vectors for mass action reactions are not allowed.")
    end

    if typeof(rs.rate_org) == Symbol
        rateconst = params[ratemap[rs.rate_org]]
    elseif typeof(rs.rate_org) == Expr
        rateconst = eval(rs.rate_org)
    elseif typeof(rs.rate_org) <: Number
        rateconst = Float64(rs.rate_org)
    else
        error("reaction_network.reactions.rate_org must have a type of Symbol, Expr, or Number.")
    end

    MassActionJump(rateconst, reactant_stoich, net_stoich)
end

# given a reaction network and species map, split the ConstantRateJumps and MassActionJumps
function network_to_jumpset(rn, specmap, ratemap, params)

    empty_majump = MassActionJump(0.0, [0=>1], [1=>1])
    majumpvec    = Vector{typeof(empty_majump)}()
    cjumpvec     = Vector{ConstantRateJump}()

    for (i,rs) in enumerate(rn.reactions)        
        if rs.is_pure_mass_action
            push!(majumpvec, make_majump(rs, specmap, ratemap, params))
        else
            push!(cjumpvec, rn.jumps[i])
        end
    end

    if isempty(majumpvec)
        return JumpSet((), cjumpvec, nothing, nothing)
    else
        return JumpSet((), cjumpvec, nothing, majumpvec)
    end
end