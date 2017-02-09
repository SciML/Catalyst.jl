immutable Reaction{T,T2,N,N2} <: AbstractReaction
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

immutable VariableRateReaction{T,T2,I,N,N2} <: AbstractReaction
  rate_constant::T
  reactants::NTuple{N,Int}
  stoichiometry::NTuple{N2,NTuple{2,T2}}
  idxs::I
  rootfind::Bool
  interp_points::Int
  abstol::T
  reltol::T2
end

function VariableRateReaction(rate_constant,reactants::AbstractArray,stoichiometry;kwargs...)
  VariableRateReaction(rate_constant,tuple(reactants...),stoichiometry;kwargs...)
end

function VariableRateReaction(rate_constant,reactants,stoichiometry::AbstractArray;kwargs...)
  VariableRateReaction(rate_constant,reactants,tuple(stoichiometry...);kwargs...)
end

function VariableRateReaction(rate_constant,reactants::AbstractArray,stoichiometry::AbstractArray;kwargs...)
  VariableRateReaction(rate_constant,tuple(reactants...),tuple(stoichiometry...);kwargs...)
end

function VariableRateReaction(rate_constant,reactants::Tuple,stoichiometry::Tuple;
                              idxs = nothing,
                              rootfind=true,
                              interp_points=10,
                              abstol=1e-12,reltol=0)
  VariableRateReaction(rate_constant,reactants,
                       stoichiometry,idxs,
                       rootfind,interp_points,
                       abstol,reltol)
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

function build_jumps_from_reaction(r::VariableRateReaction;save_positions=(false,true))
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
  VariableRateJump(rate,affect!,
             r.idxs,r.rootfind,r.interp_points,
             save_positions,r.abstol,r.reltol)
end

function build_jumps_from_reaction(rs::AbstractReaction...;save_positions=(false,true))
  build_jumps_from_reaction(ReactionSet(rs...))
end

function build_jumps_from_reaction(r::ReactionSet;save_positions=(false,true))
  jumps = []
  @inbounds for i in eachindex(r.reactions)
    push!(jumps,build_jumps_from_reaction(r.reactions[i]))
  end
  JumpSet(jumps...)
end
