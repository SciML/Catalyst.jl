### SDEProblem ###
DiffEqBase.SDEProblem(rn::AbstractReactionNetwork, u0::Union{AbstractArray, Number}, args...; kwargs...) =
    SDEProblem(rn, rn.g, u0, args...;noise_rate_prototype=rn.p_matrix, kwargs...)

### JumpProblem ###
function DiffEqJump.JumpProblem(prob,aggregator,rn::AbstractReactionNetwork; kwargs...)
    if typeof(prob)<:DiscreteProblem && any(issubtype.(typeof.(rn.jumps),VariableRateJump))
        error("When using time dependant reaction rates a DiscreteProblem should not be used (try an ODEProblem). Also, use a continious solver.")
    end

    # map from species symbol to index of species
    spec_to_idx = species_to_indices(rn)

    # map from parameter symbol to index of parameter in prob.p
    param_to_idx = rate_to_indices(rn)

    # get a JumpSet of the possible jumps
    jset = network_to_jumpset(rn, spec_to_idx, param_to_idx, prob.p)

    JumpProblem(prob, aggregator, jset; kwargs...)
end

### SteadyStateProblem ###
DiffEqBase.SteadyStateProblem(rn::AbstractReactionNetwork, args...; kwargs...) =
    SteadyStateProblem(rn.f, args...; kwargs...)

function DiffEqBase.SteadyStateProblem{isinplace}(rn::AbstractReactionNetwork, args...; kwargs...) where isinplace
    SteadyStateProblem{isinplace}(rn.f, args...; kwargs...)
end
