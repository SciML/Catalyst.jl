### SDEProblem ###
DiffEqBase.SDEProblem(rn::DiffEqBase.AbstractReactionNetwork, u0::Union{AbstractArray, Number}, args...; kwargs...) =
    SDEProblem(rn, rn.g::Function, u0, args...;noise_rate_prototype=rn.p_matrix, kwargs...)

### JumpProblem ###
function DiffEqJump.JumpProblem(prob,aggregator,rn::DiffEqBase.AbstractReactionNetwork; kwargs...)
    if typeof(prob)<:DiscreteProblem && any(x->typeof(x) <: VariableRateJump, rn.jumps)
        error("When using time dependant reaction rates a DiscreteProblem should not be used (try an ODEProblem). Also, use a continious solver.")
    end

    # map from species symbol to index of species
    spec_to_idx = species_to_indices(rn)

    # map from parameter symbol to index of parameter in prob.p
    param_to_idx = rate_to_indices(rn)

    # get a JumpSet of the possible jumps
    jset = network_to_jumpset(rn, spec_to_idx, param_to_idx, prob.p)

    # # construct dependency graph if needed by aggregator
    if needs_depgraph(aggregator)
        dep_graph = depgraph_from_network(rn, spec_to_idx, jset)
    else
        dep_graph = nothing
    end

    JumpProblem(prob, aggregator, jset; dep_graph=dep_graph, kwargs...)
end

### SteadyStateProblem ###
DiffEqBase.SteadyStateProblem(rn::DiffEqBase.AbstractReactionNetwork, args...; kwargs...) =
    SteadyStateProblem(rn.f, args...; kwargs...)

function DiffEqBase.SteadyStateProblem{isinplace}(rn::DiffEqBase.AbstractReactionNetwork, args...; kwargs...) where isinplace
    SteadyStateProblem{isinplace}(rn.f, args...; kwargs...)
end
