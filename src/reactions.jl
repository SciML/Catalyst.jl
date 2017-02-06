immutable Reaction{T,T2,N,N2}
  rate_constant::T
  reactants::NTuple{N,Int}
  stoichiometry::NTuple{N2,NTuple{2,T2}}
end

function Reaction(rate_constant,reactants::AbstractArray,stoichiometry)
  Reaction(rate_constant,tuple(reactants...),stoichiometry)
end

function Reaction(rate_constant,reactants,stoichiometry::AbstractArray)
  Reaction(rate_constant,reactants,tuple(stoichiometry...))
end

function Reaction(rate_constant,reactants::AbstractArray,stoichiometry::AbstractArray)
  Reaction(rate_constant,tuple(reactants...),tuple(stoichiometry...))
end

immutable ReactionSet{T}
  reactions::T
end

ReactionSet(args...) = ReactionSet(tuple(args...))

function build_jumps_from_reaction(r::Reaction;save_positions=(false,true))
  rate = function (t,u)
    val = r.rate_constant
    @inbounds for i in eachindex(r.reactants)
      val *= u[r.reactants[i]]
    end
    val
  end
  affect! = function (integrator)
    @inbounds for i in eachindex(r.stoichiometry)
      integrator.u[r.stoichiometry[i][1]] += r.stoichiometry[i][2]
    end
  end
  ConstantRateJump(rate,affect!,save_positions)
end

function build_jumps_from_reaction(rs::Reaction...;save_positions=(false,true))
  build_jumps_from_reaction(ReactionSet(rs...))
end

function build_jumps_from_reaction(r::ReactionSet;save_positions=(false,true))
  jumps = []
  @inbounds for i in eachindex(r.reactions)
    push!(jumps,build_jumps_from_reaction(r.reactions[i]))
  end
  JumpSet(jumps...)
end
