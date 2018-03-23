### ODEProblem ###
DiffEqBase.ODEProblem(rn::AbstractReactionNetwork, args...; kwargs...) =
    ODEProblem(rn.f, args...; kwargs...)

### SDEProblem ###
DiffEqBase.SDEProblem(rn::AbstractReactionNetwork, args...; kwargs...) =
    SDEProblem(rn.f,rn.g, args...;noise_rate_prototype=rn.p_matrix, kwargs...)

### JumpProblem ###
DiffEqJump.JumpProblem(prob,aggregator::Direct,rn::AbstractReactionNetwork; kwargs...) =
    if typeof(prob)<:DiscreteProblem && any(issubtype.(typeof.(equi_model.jumps),VariableRateJump))
        error("When using time dependant reaction rates a DiscreteProblem should not be used (try an ODEProblem). Also, use a continious solver.")
    end
    JumpProblem(prob,aggregator::Direct,rn.jumps...;kwargs...)

### SteadyStateProblem ###
DiffEqBase.SteadyStateProblem(rn::AbstractReactionNetwork, args...; kwargs...) =
    SteadyStateProblem(rn.f, args...; kwargs...)

function DiffEqBase.SteadyStateProblem{isinplace}(rn::AbstractReactionNetwork, args...; kwargs...) where isinplace
    SteadyStateProblem{isinplace}(rn.f, args...; kwargs...)
end
