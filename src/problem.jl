GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,
                 rs::AbstractReaction...;kwargs...) =
                 JumpProblem(prob,aggregator,build_jumps_from_reaction(rs...);kwargs...)

GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,
                 r::ReactionSet;kwargs...) =
                 JumpProblem(prob,aggregator,build_jumps_from_reaction(r);kwargs...)

### Added by Me ###
### ODEProblem ###
newODEProblem(ar::AbstractReaction, args...; kwargs...) =
                 ODEProblem(ar.f, args...; kwargs...)

### SDEProblem ###
newSDEProblem(ar::AbstractReaction, args...; kwargs...) =
                 SDEProblem(ar.f,ar.g, args...;noise_rate_prototype=rn.p_matrix, kwargs...)

### JumpProblem ###
newJumpProblem(prob,aggregator::Direct,ar::AbstractReaction) =
                 JumpProblem(prob,aggregator::Direct,ar.jumps...)
