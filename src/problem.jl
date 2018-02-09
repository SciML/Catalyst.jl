### ODEProblem ###
DiffEqBase.ODEProblem(ar::AbstractReaction, args...; kwargs...) =
                 ODEProblem(ar.f, args...; kwargs...)

### SDEProblem ###
DiffEqBase.SDEProblem(ar::AbstractReaction, args...; kwargs...) =
                 SDEProblem(ar.f,ar.g, args...;noise_rate_prototype=rn.p_matrix, kwargs...)

### JumpProblem ###
DiffEqJump.JumpProblem(prob,aggregator::Direct,ar::AbstractReaction) =
                 JumpProblem(prob,aggregator::Direct,ar.jumps...)
