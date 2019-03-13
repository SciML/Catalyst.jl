# collection of utility functions for getting mass action jumps out of
# DiffEqBiological reaction network structures

# given a ReactionStruct and a species map construct a MassActionJump
function make_majump(rs, specmap, ratemap, params, param_context)
    reactant_stoich = substratestoich(rs, specmap)
    net_stoich      = netstoich(rs, specmap)
    if isempty(net_stoich)
        error("Empty net stoichiometry vectors for mass action reactions are not allowed.")
    end

    if typeof(rs.rate_org) == Symbol
        rateconst = params[ratemap[rs.rate_org]]
    elseif typeof(rs.rate_org) == Expr
        # eval in param_context, in case Expr depends on params
        rateconst = Base.eval(param_context, rs.rate_org)
    elseif typeof(rs.rate_org) <: Number
        rateconst = rs.rate_org
    else
        error("reaction_network.reactions.rate_org must have a type of Symbol, Expr, or Number.")
    end

    MassActionJump(Float64(rateconst), reactant_stoich, net_stoich)
end

# given a reaction network and species map, split the ConstantRateJumps and MassActionJumps
function network_to_jumpset(rn, params)
    jumps        = rn.jumps
    empty_majump = MassActionJump(0.0, [0=>1], [1=>1])
    majumpvec    = Vector{typeof(empty_majump)}()
    cjumpvec     = Vector{ConstantRateJump}()
    specmap      = speciesmap(rn)
    ratemap      = paramsmap(rn)

    # populate dummy module with params as local variables
    # (for eval-ing parameter expressions)
    param_context = Module()
    for (param, index) in ratemap
        Base.eval(param_context, :($param = $(params[index])))
    end

    # check if a (non-mass action) jump is defined for all reactions:
    alljumps = (length(rn.reactions) == length(jumps))
    
    idx = 1
    for rs in rn.reactions
        if rs.is_pure_mass_action
            push!(majumpvec, make_majump(rs, specmap, ratemap, params, param_context))
            alljumps && (idx += 1)
        else
            push!(cjumpvec, jumps[idx])
            idx += 1
        end
    end

    if isempty(majumpvec)
        return JumpSet((), cjumpvec, nothing, nothing)
    else
        return JumpSet((), cjumpvec, nothing, majumpvec)
    end
end


######################### jump dependency graph utilities #########################

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
function jump_to_dep_specs_map(rn, rxidxs_jidxs)
    numrxs  = length(rn.reactions)
    specmap = speciesmap(rn)
    numspec = length(specmap)

    # map from a jump to vector of species that depend on it
    jtos_vec = Vector{Vector{valtype(specmap)}}(undef, numrxs)
    for rx in 1:numrxs
        jidx           = rxidxs_jidxs[rx]
        jtos_vec[jidx] = sort!( [ns.first for ns in netstoich(rn.reactions[rx], specmap)] )
    end

    jtos_vec
end

# map from species to Set of jumps depending on that species
function spec_to_dep_jumps_map(rn, rxidxs_to_jidxs)
    numrxs  = length(rn.reactions)
    specmap = speciesmap(rn)
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
function depgraph_from_network(rn, jset, rxidxs_to_jidxs, spec_to_dep_jumps)

    # create map from a jump to jumps depending on it
    numrxs   = length(rn.reactions)
    dep_sets = [SortedSet{Int}() for n = 1:numrxs]
    for rx in 1:numrxs
        jidx = rxidxs_to_jidxs[rx]

        # get the net reaction stoichiometry
        net_stoich = netstoich(rn.reactions[rx], speciesmap(rn))

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
