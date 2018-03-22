### ODEProblem ###
DiffEqBase.ODEProblem(rn::AbstractReactionNetwork, args...; kwargs...) =
    ODEProblem(rn.f, args...; kwargs...)

### SDEProblem ###
DiffEqBase.SDEProblem(rn::AbstractReactionNetwork, args...; kwargs...) =
    SDEProblem(rn.f,rn.g, args...;noise_rate_prototype=rn.p_matrix, kwargs...)

### JumpProblem ###
DiffEqJump.JumpProblem(prob,aggregator::Direct,rn::AbstractReactionNetwork; kwargs...) =
    JumpProblem(prob,aggregator::Direct,rn.jumps...;kwargs...)

### SteadyStateProblem ###
DiffEqBase.SteadyStateProblem(rn::AbstractReactionNetwork, args...; kwargs...) =
    SteadyStateProblem(rn.f, args...; kwargs...)

function DiffEqBase.SteadyStateProblem{isinplace}(rn::AbstractReactionNetwork, args...; kwargs...) where isinplace
    SteadyStateProblem{isinplace}(rn.f, args...; kwargs...)
end
