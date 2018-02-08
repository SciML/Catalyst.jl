GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,
                 rs::AbstractReaction...;kwargs...) =
                 JumpProblem(prob,aggregator,build_jumps_from_reaction(rs...);kwargs...)

GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,
                 r::ReactionSet;kwargs...) =
                 JumpProblem(prob,aggregator,build_jumps_from_reaction(r);kwargs...)

### Added by Me ###
### ODEProblem ###
ODEProblem(ar::AbstractReaction, args...; kwargs...) =
                 ODEProblem(ar.f, args...; kwargs...)

### SDEProblem ###
SDEProblem(ar::AbstractReaction, args...; kwargs...) =
                 SDEProblem(ar.f,ar.g, args...;noise_rate_prototype=rn.p_matrix, kwargs...)

### JumpProblem ###
JumpProblem(prob,aggregator::Direct,ar::AbstractReaction) =
                 JumpProblem(prob,aggregator::Direct,ar.jumps...)
