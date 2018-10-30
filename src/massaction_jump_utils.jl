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

# get substrate stoichiometry for a reaction
function get_substrate_stoich(rs, specmap)
    reactant_stoich = Vector{Pair{Int,Int}}()
    for substrate in rs.substrates
        specpair = specmap[substrate.reactant] => substrate.stoichiometry
        push!(reactant_stoich, specpair)
    end
    sort!(reactant_stoich)
    reactant_stoich
end

# get the net stoichiometry for a reaction
function get_net_stoich(rs, specmap)
    nsdict = Dict{Int,Int}()
    for substrate in rs.substrates
        specpair = specmap[substrate.reactant] => -substrate.stoichiometry
        push!(nsdict, specpair)
    end

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
        if stoich_map[2] != zero(Int)
            push!(net_stoich, stoich_map)
        end
    end

    net_stoich
end

# given a ReactionStruct and a species map construct a MassActionJump
function make_majump(rs, specmap, ratemap, params)
    reactant_stoich = get_substrate_stoich(rs, specmap)
    net_stoich      = get_net_stoich(rs, specmap)
    if isempty(net_stoich)
        error("Empty net stoichiometry vectors for mass action reactions are not allowed.")
    end

    if typeof(rs.rate_org) == Symbol
        rateconst = params[ratemap[rs.rate_org]]
    elseif typeof(rs.rate_org) == Expr
        rateconst = eval(rs.rate_org)
    elseif typeof(rs.rate_org) <: Number
        rateconst = rs.rate_org
    else
        error("reaction_network.reactions.rate_org must have a type of Symbol, Expr, or Number.")
    end

    MassActionJump(Float64(rateconst), reactant_stoich, net_stoich)
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


######################### dependency graph utilities #########################

# map from a reaction index to the corresponding jump index
# assumes all mass action jumps are ordered before constant rate jumps
function rxidxs_to_jidxs_map(rn, num_majumps)
    majidx = 1
    cjidx  = num_majumps + 1
    numrxs = length(rn.reactions)
    rxi_to_ji = zeros(Int, numrxs)
    for i = 1:numrxs
        if rn.reactions[i].is_pure_mass_action
            rxi_to_ji[i] = majidx
            majidx      += 1
        else
            rxi_to_ji[i] = cjidx
            cjidx       += 1
        end
    end
    rxi_to_ji
end

# map from jump to vector of species it changes
function jump_to_dep_specs_map(rn, specmap, rxidxs_jidxs)
    numrxs  = length(rn.reactions)
    numspec = length(specmap)

    # map from a jump to species that depend on it
    jtos_vec = Vector{Vector{valtype(specmap)}}(undef, numrxs)
    for rx in 1:numrxs
        jidx = rxidxs_jidxs[rx]
        jtos_vec[jidx] = [specmap[specsym] for specsym in rn.reactions[rx].dependants]
        sort!(jtos_vec[jidx])
    end

    jtos_vec
end

# map from species to Set of jumps depending on that species
function spec_to_dep_jumps_map(rn, specmap, rxidxs_to_jidxs)
    numrxs  = length(rn.reactions)
    numspec = length(specmap)

    # map from a species to jumps that depend on it
    spec_to_dep_jumps = [Set{Int}() for n = 1:numspec]
    for rx in 1:numrxs
        for specsym in rn.reactions[rx].dependants
            push!(spec_to_dep_jumps[specmap[specsym]], rxidxs_to_jidxs[rx])
        end
    end

    spec_to_dep_jumps
end

# given a reaction network and species map, construct a jump dependency graph
function depgraph_from_network(rn, specmap, jset, rxidxs_to_jidxs, spec_to_dep_jumps)

    # create map from a jump to jumps depending on it
    numrxs   = length(rn.reactions)
    dep_sets = [SortedSet{Int}() for n = 1:numrxs]
    for rx in 1:numrxs
        jidx = rxidxs_to_jidxs[rx]

        # get the net reaction stoichiometry
        net_stoich = get_net_stoich(rn.reactions[rx], specmap)

        # rx changes spec, hence rxs depending on spec depend on rx
        for (spec,stoch) in net_stoich
            for dependent_jump in spec_to_dep_jumps[spec]
                push!(dep_sets[jidx], dependent_jump)
            end
        end
    end

    # convert to Vectors of Vectors
    dep_graph = Vector{Vector{Int}}(undef,numrxs)
    for jidx in 1:numrxs
        dep_graph[jidx] = [dep for dep in dep_sets[jidx]]
    end

    dep_graph
end
