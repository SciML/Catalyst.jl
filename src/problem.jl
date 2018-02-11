### ODEProblem ###
DiffEqBase.ODEProblem(rn::AbstractReactionNetwork, args...; kwargs...) =
                 ODEProblem(rn.f, args...; kwargs...)

### SDEProblem ###
DiffEqBase.SDEProblem(rn::AbstractReactionNetwork, args...; kwargs...) =
                 SDEProblem(rn.f,rn.g, args...;noise_rate_prototype=rn.p_matrix, kwargs...)

### JumpProblem ###
DiffEqJump.JumpProblem(prob,aggregator::Direct,rn::AbstractReactionNetwork) =
                 JumpProblem(prob,aggregator::Direct,rn.jumps...)
