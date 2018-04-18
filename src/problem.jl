### SDEProblem ###
DiffEqBase.SDEProblem(rn::AbstractReactionNetwork, u0::Union{AbstractArray, Number}, args...; kwargs...) =
    SDEProblem(rn, rn.g, u0, args...;noise_rate_prototype=rn.p_matrix, kwargs...)

### JumpProblem ###
function DiffEqJump.JumpProblem(prob,aggregator::Direct,rn::AbstractReactionNetwork; kwargs...)
    if typeof(prob)<:DiscreteProblem && any(issubtype.(typeof.(rn.jumps),VariableRateJump))
        error("When using time dependant reaction rates a DiscreteProblem should not be used (try an ODEProblem). Also, use a continious solver.")
    end
    JumpProblem(prob,aggregator::Direct,rn.jumps...;kwargs...)
end

### SteadyStateProblem ###
DiffEqBase.SteadyStateProblem(rn::AbstractReactionNetwork, args...; kwargs...) =
    SteadyStateProblem(rn.f, args...; kwargs...)

function DiffEqBase.SteadyStateProblem{isinplace}(rn::AbstractReactionNetwork, args...; kwargs...) where isinplace
    SteadyStateProblem{isinplace}(rn.f, args...; kwargs...)
end
