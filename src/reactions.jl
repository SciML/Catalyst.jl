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
    # for higher order interactions i.e. reactants [1,1] vaule = k*u[1]*(u[1]-1)
    # Since a single entitiy cannot interact with it self! 
    k = 0
    prev_reactant = -1
    @fastmath @inbounds for i in eachindex(r.reactants)
      if prev_reactant == r.reactants[i]
        k += 1
      end
      val *= (u[r.reactants[i]] - k)
      prev_reactant = r.reactants[i]    
    end
    val
  end
  affect! = function (integrator)
    @fastmath @inbounds for i in eachindex(r.stoichiometry)
      integrator.u[r.stoichiometry[i][1]] += r.stoichiometry[i][2]
    end
  end
  ConstantRateJump(rate,affect!,save_positions)
end

function build_jumps_from_reaction(r::VariableRateReaction;save_positions=(false,true))
  rate = function (t,u)
    val = r.rate_constant
    # for higher order interactions i.e. reactants [1,1] vaule = k*u[1]*(u[1]-1)
    # Since a single entitiy cannot interact with it self! 
    k = 0
    prev_reactant = -1
    @fastmath @inbounds for i in eachindex(r.reactants)
      if prev_reactant == r.reactants[i]
        k += 1
      end
      val *= (u[r.reactants[i]] - k)
      prev_reactant = r.reactants[i]    
    end
    val
  end
  affect! = function (integrator)
    @fastmath @inbounds for i in eachindex(r.stoichiometry)
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


macro reaction_network(ex)
  def_reaction_network(ex)
end

function def_reaction_network(ex::Expr)
  reactions = :([])
  mapping = OrderedDict{Symbol,Int}()

  for line in ex.args
    if line.head == :tuple
      rate = line.args[1]
      eqn  = line.args[2]

      reactants = :([])     # indices
      stoichiometry = :([]) # net stoichiometry

      reacs, prods = parse_reaction(eqn, mapping)

      for (id, ix) in mapping
        if haskey(reacs, id)
          push!(reactants.args, ix)
        end

        net = get(prods, id, 0) - get(reacs, id, 0)

        if net != 0
          push!(stoichiometry.args, :(($ix, $net)))
        end
      end

      push!(reactions.args, :(Reaction($rate, $reactants, $stoichiometry)))
    end
  end

  return :(ReactionSet($reactions...))
end

function parse_reaction(ex::Expr, mapping)
  reactants = Dict{Symbol,Int}()
  products  = Dict{Symbol,Int}()

  if ex.head == :-->
    exr = ex.args[1] # LHS
    exp = ex.args[2] # RHS

    add_participants!(reactants, exr, mapping)
    add_participants!(products,  exp, mapping)
  else
    throw("malformed reaction")
  end

  return reactants, products
end

function add_participants!(dict, ex, mapping)
  # found a species symbol
  if isa(ex, Symbol)
    if !haskey(mapping, ex)
      mapping[ex] = length(mapping) + 1
    end

    id = ex
    val = get(dict, id, 0)

    dict[id] = val + 1
  # species symbol has a coefficient
  elseif isa(ex, Expr) && ex.args[1] == :*
    id    = ex.args[3]
    coeff = ex.args[2]

    if !haskey(mapping, id)
      mapping[id] = length(mapping) + 1
    end

    val = get(dict, id, 0)
    dict[id] = val + coeff
  # found something else, probably of the form (a A + b B)
  elseif isa(ex, Expr)
    for i in 2:length(ex.args)
      add_participants!(dict, ex.args[i], mapping)
    end
  end
end
