# Functions for querying network properties.

######### Accessors: #########

"""
    species(network)

Given a [`ReactionSystem`](@ref), return a vector of species `Variable`s.

Notes:
- If `ModelingToolkit.get_systems(network)` is not empty, may allocate. Otherwise returns
  `ModelingToolkit.get_states(network)`.
"""
function species(network)
    isempty(get_systems(network)) ? get_states(network) : states(network)
end

"""
    params(network)

Given a [`ReactionSystem`](@ref), return a vector of parameter `Variable`s.

Notes:
- If `ModelingToolkit.get_systems(network)` is not empty, may allocate. Otherwise returns
  `ModelingToolkit.get_ps(network)`.
"""
function params(network)
    isempty(get_systems(network)) ? get_ps(network) : parameters(network)
end

"""
    reactions(network)

Given a [`ReactionSystem`](@ref), return a vector of all `Reactions` in the system.

Notes:
- If `ModelingToolkit.get_systems(network)` is not empty, may allocate. Otherwise returns
  `ModelingToolkit.get_eqs(network)`.
"""
function reactions(network)
    isempty(get_systems(network)) ? get_eqs(network) : equations(network)
end

"""
    speciesmap(network)

Given a [`ReactionSystem`](@ref), return a Dictionary mapping from species to
species indices. (Allocates)
"""
function speciesmap(network)
    Dict(S => i for (i,S) in enumerate(species(network)))
end

"""
    paramsmap(network)

Given a [`ReactionSystem`](@ref), return a Dictionary mapping from parameters to
parameter indices. (Allocates)
"""
function paramsmap(network)
    Dict(p => i for (i,p) in enumerate(params(network)))
end

"""
    numspecies(network)

Return the number of species within the given [`ReactionSystem`](@ref).
"""
function numspecies(network)
    ns = length(get_states(network))
    for sys in get_systems(network)
        ns += numspecies(sys)
    end
    ns
end

"""
    numreactions(network)

Return the number of reactions within the given [`ReactionSystem`](@ref).
"""
function numreactions(network)
    nr = length(get_eqs(network))
    for sys in get_systems(network)
        nr += numreactions(sys)
    end
    nr
end

"""
    numparams(network)

Return the number of parameters within the given [`ReactionSystem`](@ref).
"""
function numparams(network)
    np = length(get_ps(network))
    for sys in get_systems(network)
        np += numparams(sys)
    end
    np
end

"""
    dependents(rx, network)

Given a [`Reaction`](@ref) and a [`ReactionSystem`](@ref), return a vector of
`ModelingToolkit.Num`s corresponding to species the *reaction rate
law* depends on. E.g., for

`k*W, 2X + 3Y --> 5Z + W`

the returned vector would be `[W(t),X(t),Y(t)]`.

Notes:
- Allocates
- Does not check for dependents within any subsystems.
"""
function dependents(rx, network)
    if rx.rate isa Number
        return rx.substrates
    else
        rvars = get_variables(rx.rate, species(network))
        return union!(rvars, rx.substrates)
    end
end

"""
    dependents(rx, network)

See documentation for [`dependents`](@ref).
"""
function dependants(network, rxidx)
    dependents(network, rxidx)
end


"""
    substoichmat(rn; sparse=false, smap=speciesmap(rn))

Returns the substrate stoichiometry matrix, ``S``, with ``S_{i j}`` the
stoichiometric coefficient of the ith substrate within the jth reaction.

Note:
- Set sparse=true for a sparse matrix representation
"""
function substoichmat(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem; smap=speciesmap(rn))
    Is=Int[];  Js=Int[];  Vs=Int[];
    for (k,rx) in enumerate(reactions(rn))
        stoich = rx.substoich
        for (i,sub) in enumerate(rx.substrates)
            push!(Js, k)
            push!(Is, smap[sub])
            push!(Vs, stoich[i])
        end
    end
    sparse(Is,Js,Vs,numspecies(rn),numreactions(rn))
end
function substoichmat(::Type{Matrix{Int}},rn::ReactionSystem; smap=speciesmap(rn))
    smat = zeros(Int, numspecies(rn), numreactions(rn))
    for (k,rx) in enumerate(reactions(rn))
        stoich = rx.substoich
        for (i,sub) in enumerate(rx.substrates)
            smat[smap[sub],k] = stoich[i]
        end
    end
    smat
end
function substoichmat(rn::ReactionSystem; sparse::Bool=false, smap=speciesmap(rn))
	sparse ? substoichmat(SparseMatrixCSC{Int,Int}, rn; smap=smap) : substoichmat(Matrix{Int}, rn; smap=smap)
end


"""
    prodstoichmat(rn; sparse=false, smap=speciesmap(rn))

Returns the product stoichiometry matrix, ``P``, with ``P_{i j}`` the
stoichiometric coefficient of the ith product within the jth reaction.

Note:
- Set sparse=true for a sparse matrix representation
"""
function prodstoichmat(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem; smap=speciesmap(rn))
    Is=Int[];  Js=Int[];  Vs=Int[];
    for (k,rx) in enumerate(reactions(rn))
        stoich = rx.prodstoich
        for (i,prod) in enumerate(rx.products)
			push!(Js, k)
			push!(Is, smap[prod])
			push!(Vs, stoich[i])
        end
    end
    sparse(Is,Js,Vs,numspecies(rn),numreactions(rn))
end
function prodstoichmat(::Type{Matrix{Int}},rn::ReactionSystem; smap=speciesmap(rn))
    pmat = zeros(Int, numspecies(rn), numreactions(rn))
    for (k,rx) in enumerate(reactions(rn))
        stoich = rx.prodstoich
        for (i,prod) in enumerate(rx.products)
            pmat[smap[prod],k] = stoich[i]
        end
    end
    pmat
end
function prodstoichmat(rn::ReactionSystem; sparse=false, smap=speciesmap(rn))
	sparse ? prodstoichmat(SparseMatrixCSC{Int,Int}, rn; smap=smap) : prodstoichmat(Matrix{Int}, rn; smap=smap)
end


"""
    netstoichmat(rn,sparsity=false; smap=speciesmap(rn))

Returns the net stoichiometry matrix, ``N``, with ``N_{i j}`` the net
stoichiometric coefficient of the ith species within the jth reaction.

Note:
- Set sparse=true for a sparse matrix representation
"""
function netstoichmat(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem; smap=speciesmap(rn))
    Is=Int[];  Js=Int[];  Vs=Int[];
    for (k,rx) in pairs(reactions(rn))
        for (spec,coef) in rx.netstoich
			push!(Js, k)
			push!(Is, smap[spec])
			push!(Vs, coef)
        end
    end
    sparse(Is,Js,Vs,numspecies(rn),numreactions(rn))
end
function netstoichmat(::Type{Matrix{Int}},rn::ReactionSystem; smap=speciesmap(rn))
    nmat = zeros(Int,numspecies(rn),numreactions(rn))
    for (k,rx) in pairs(reactions(rn))
        for (spec,coef) in rx.netstoich
            nmat[smap[spec],k] = coef
        end
    end
    nmat
end
function netstoichmat(rn::ReactionSystem; sparse=false, smap=speciesmap(rn))
	sparse ? netstoichmat(SparseMatrixCSC{Int,Int}, rn; smap=smap) : netstoichmat(Matrix{Int}, rn; smap=smap)
end


######################## reaction complexes and reaction rates ###############################
"""
$(TYPEDEF)
One reaction complex element

# Fields
$(FIELDS)
"""
struct ReactionComplexElement{T}
    """The integer id of the species representing this element."""
    speciesid::Int
    """The stoichiometric coefficient of this species."""
    speciesstoich::T
end

"""
$(TYPEDEF)
One reaction complex.

# Fields
$(FIELDS)
"""
struct ReactionComplex{V<:Integer} <: AbstractVector{ReactionComplexElement{V}}
    """The integer ids of all species participating in this complex."""
    speciesids::Vector{Int}
    """The stoichiometric coefficients of all species participating in this complex."""
    speciesstoichs::Vector{V}
end

function (==)(a::ReactionComplex{V},b::ReactionComplex{V}) where {V <: Integer} 
    (a.speciesids == b.speciesids) &&
    (a.speciesstoichs == b.speciesstoichs) 
end
hash(rc::ReactionComplex,h::UInt) = Base.hash(rc.speciesids,Base.hash(rc.speciesstoichs,h))
Base.size(rc::ReactionComplex) = size(rc.speciesids)
Base.length(rc::ReactionComplex) = length(rc.speciesids)
Base.getindex(rc::ReactionComplex, i...) = 
        ReactionComplexElement(getindex(rc.speciesids, i...), getindex(rc.speciesstoichs, i...))
Base.setindex!(rc::ReactionComplex, t::ReactionComplexElement, i...) = 
    (setindex!(rc.speciesids, t.speciesid, i...); setindex!(rc.speciesstoichs, t.speciesstoich, i...); rc) 
Base.isless(a::ReactionComplexElement, b::ReactionComplexElement) = isless(a.speciesid, b.speciesid)
Base.Sort.defalg(::ReactionComplex{T}) where {T <: Integer} = Base.DEFAULT_UNSTABLE

"""
    reactioncomplexmap(rn::ReactionSystem; smap=speciesmap(rn))

Find each [`ReactionComplex`](@ref) within the specified system, constructing a
mapping from the complex to vectors that indicate which reactions it appears in
as substrates and products.

Notes:
- Each [`ReactionComplex`](@ref) is mapped to a vector of pairs, with each pair
  having the form `reactionidx => ± 1`, where `-1` indicates the complex appears
  as a substrate and `+1` as a product in the reaction with integer label
  `reactionidx`.
"""
function reactioncomplexmap(rn::ReactionSystem; smap=speciesmap(rn))
    rxs = reactions(rn)
    numreactions(rn) > 0 || error("There must be at least one reaction to find reaction complexes.")
    complextorxsmap = OrderedDict{ReactionComplex{eltype(rxs[1].substoich)},Vector{Pair{Int,Int}}}()
    for (i,rx) in enumerate(rxs)
        reactantids = isempty(rx.substrates) ? Vector{Int}() : [smap[sub] for sub in rx.substrates]
        subrc = sort!(ReactionComplex(reactantids, copy(rx.substoich)))
        if haskey(complextorxsmap, subrc)
            push!(complextorxsmap[subrc], i => -1)
        else
            complextorxsmap[subrc] = [i => -1]
        end

        productids = isempty(rx.products) ? Vector{Int}() : [smap[prod] for prod in rx.products]
        prodrc = sort!(ReactionComplex(productids, copy(rx.prodstoich)))
        if haskey(complextorxsmap, prodrc)
            push!(complextorxsmap[prodrc], i => 1)
        else
            complextorxsmap[prodrc] = [i => 1]
        end
    end
    complextorxsmap
end


@doc raw"""
    reactioncomplexes(network; sparse=false, smap=speciesmap(rn))

Calculate the reaction complexes and complex incidence matrix for the given [`ReactionSystem`](@ref). 

Notes:
- returns a pair of a vector of [`ReactionComplex`](@ref)s and the complex incidence matrix.
- An empty [`ReactionComplex`](@ref) denotes the null (∅) state (from reactions like ∅ -> A or A -> ∅).
- The complex incidence matrix, ``B``, is number of complexes by number of reactions with
```math
B_{i j} = \begin{cases}
-1, &\text{if the i'th complex is the substrate of the j'th reaction},\\
1, &\text{if the i'th complex is the product of the j'th reaction},\\
0, &\text{otherwise.}
\end{cases}
```
- Set sparse=true for a sparse matrix representation of the incidence matrix
"""
function reactioncomplexes(::Type{SparseMatrixCSC{Int,Int}},rn::ReactionSystem; 
                            smap=speciesmap(rn), complextorxsmap=reactioncomplexmap(rn; smap=smap))
    complexes = collect(keys(complextorxsmap))
    Is=Int[];  Js=Int[];  Vs=Int[];
	for (i,c) in enumerate(complexes)
        for (j,σ) in complextorxsmap[c]
			push!(Is, i)
			push!(Js, j)
			push!(Vs, σ)
        end
    end
	B = sparse(Is,Js,Vs,length(complexes),numreactions(rn))
    complexes,B
end
function reactioncomplexes(::Type{Matrix{Int}},rn::ReactionSystem; 
                            smap=speciesmap(rn), complextorxsmap=reactioncomplexmap(rn; smap=smap))
    complexes = collect(keys(complextorxsmap))
    B = zeros(Int, length(complexes), numreactions(rn));
    for (i,c) in enumerate(complexes)
        for (j,σ) in complextorxsmap[c]
            B[i,j] = σ
        end
    end
    complexes,B
end
function reactioncomplexes(rn::ReactionSystem; sparse=false, smap=speciesmap(rn))
	sparse ? reactioncomplexes(SparseMatrixCSC{Int,Int}, rn; smap=smap) : reactioncomplexes(Matrix{Int}, rn; smap=smap)
end


"""
    reaction_rates(network)

Given a [`ReactionSystem`](@ref), returns a vector of the symbolic reaction
rates for each reaction.
"""
function reactionrates(rn)
    [r.rate for r in reactions(rn)]
end


"""
    complexstoichmat(network; sparse=false, rcs=keys(reactioncomplexmap(rn)))

Given a [`ReactionSystem`](@ref) and vector of reaction complexes, return a
matrix with positive entries of size num_of_species x num_of_complexes, where
the non-zero positive entries in the kth column denote stoichiometric
coefficients of the species participating in the kth reaction complex.

Notes:
- `rcs` correspond to an iterable of the `ReactionComplexes`, i.e.
  `rcs=keys(reactioncomplexmap(rn))` or `reactioncomplexes(rn)[1]`.
- Set sparse=true for a sparse matrix representation
"""
function complexstoichmat(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem; rcs=keys(reactioncomplexmap(rn)))
    Is=Int[];  Js=Int[];  Vs=Int[];
    for (i,rc) in enumerate(rcs)
        for rcel in rc
			push!(Is, rcel.speciesid)
			push!(Js, i)
			push!(Vs, rcel.speciesstoich)
        end
    end
    Z = sparse(Is,Js,Vs, numspecies(rn), length(rcs))
end
function complexstoichmat(::Type{Matrix{Int}}, rn::ReactionSystem; rcs=keys(reactioncomplexmap(rn)))
    Z=zeros(Int, numspecies(rn), length(rcs))
    for (i,rc) in enumerate(rcs)
        for rcel in rc
            Z[rcel.speciesid,i] = rcel.speciesstoich
        end
    end
    Z
end
function complexstoichmat(rn::ReactionSystem; sparse=false, rcs=keys(reactioncomplexmap(rn)))
	sparse ? complexstoichmat(SparseMatrixCSC{Int,Int}, rn; rcs=rcs) : complexstoichmat(Matrix{Int}, rn; rcs=rcs)
end

@doc raw"""
    complexoutgoingmat(network; sparse=false, B=reactioncomplexes(rn)[2])

Given a [`ReactionSystem`](@ref) and complex incidence matrix, ``B``, return a
matrix of size num of complexes by num of reactions that identifies substrate
complexes.

Notes:
- The complex outgoing matrix, ``\Delta``, is defined by 
```math
\Delta_{i j} = \begin{cases}
    = 0,    &\text{if } B_{i j} = 1, \\
    = B_{i j}, &\text{otherwise.}
\end{cases}
```
- Set sparse=true for a sparse matrix representation
"""
function complexoutgoingmat(::Type{SparseMatrixCSC{Int,Int}}, rn::ReactionSystem; B=reactioncomplexes(rn,sparse=true)[2])    
    n = size(B,2)
	rows = rowvals(B)
	vals = nonzeros(B)
    Is = Int[]; Js = Int[]; Vs = Int[]
    sizehint!(Is, div(length(vals),2))
    sizehint!(Js, div(length(vals),2))
    sizehint!(Vs, div(length(vals),2))
	for j = 1:n
	   for i in nzrange(B, j)
	      if vals[i] != one(eltype(vals)) 
            push!(Is, rows[i])
            push!(Js, j) 
            push!(Vs, vals[i])
          end
	   end
	end
    sparse(Is,Js,Vs,size(B,1),size(B,2))
end
function complexoutgoingmat(::Type{Matrix{Int}}, rn::ReactionSystem; B=reactioncomplexes(rn,sparse=false)[2])
    Δ = copy(B)
    for (I,b) in pairs(Δ)
        (b == 1) && (Δ[I] = 0)
    end
    Δ
end
function complexoutgoingmat(rn::ReactionSystem; sparse=false, B=reactioncomplexes(rn,sparse=sparse)[2])
	sparse ? complexoutgoingmat(SparseMatrixCSC{Int,Int}, rn; B=B) : complexoutgoingmat(Matrix{Int}, rn; B=B)
end


"""
   incidencematgraph(incidencemat::Matrix)
   incidencematgraph(incidencemat::SparseMatrixCSC{Int,Int})

Given an incidence matrix of reaction-network in dense or sparse matrix form,
writtens a {n ,m} directed simple graph,
where n is number of nodes(reaction-complexes)
and m is number of edges(reactions)
"""
function incidencematgraph(incidencemat::Matrix{Int})
   @assert all(∈([-1,0,1]) ,incidencemat)
   n = size(incidencemat,1)  # no. of nodes/complexes
   graph = LightGraphs.DiGraph(n)
   for col in eachcol(incidencemat)
       src = 0; dst = 0;
       for i in eachindex(col)
              (col[i] == -1) && (src = i)
              (col[i] == 1) && (dst = i)
              (src != 0) && (dst != 0) && break
       end
       add_edge!(graph, src, dst)
    end
   return graph
end
function incidencematgraph(incidencemat::SparseMatrixCSC{Int,Int})
   @assert all(∈([-1,0,1]) ,incidencemat)
   m,n = size(incidencemat)  
   graph = LightGraphs.DiGraph(m)
   rows = rowvals(incidencemat)
   vals = nonzeros(incidencemat)
   for j = 1:n
      inds=nzrange(incidencemat, j)
      row = rows[inds];
      val = vals[inds];
      if val[1] == -1
         add_edge!(graph, row[1], row[2])
      else
         add_edge!(graph, row[2], row[1])
      end
   end
   return graph
end


"""
      linkageclasses(incidencegraph)

Given a directed simple graph created from the incidence matrix
of reaction network, returns Vector of indices of directly connected
reaction-complexes(more precisely indices of reactioncomplexes(network)[1])
in reactions
"""
function linkageclasses(incidencegraph::SimpleDiGraph{Int64})
   LightGraphs.connected_components(incidencegraph)
end


"""
   deficiency(network ; sparse=false ,ns ,ig ,lc )

Given a ReactionSystem, returns deficiency of the reaction network
defined as,
deficiency = number of reaction complexes - number of linkage classes - dimension of stochiometric subspace

Note: The kwargs are as follows:
- ns , netstoichmat of the network
- ig, incidencegraph from incidence matrix of the network
- lc , linkage classes from incidencegraph of the network
- sparse, a Bool,defaulted as false, for constructing a dense or sparse netstoichmat ,incidence matrix, etc.
"""
function deficiency(rs::ReactionSystem ; sparse::Bool=false, ns::AbstractMatrix = netstoichmat(rs;sparse=sparse),
                        ig::SimpleDiGraph{Int64} = incidencematgraph(reactioncomplexes(rs;sparse=sparse)[2]),
                        lc::Vector{Vector{Int64}} = linkageclasses(ig))

   LightGraphs.nv(ig) - length(lc) - AbstractAlgebra.rank(AbstractAlgebra.matrix(AbstractAlgebra.ZZ,ns))
end


################################################################################################
######################## conservation laws ###############################

""" 
    conservationlaws(netstoichmat::AbstractMatrix)::Matrix

Given the net stoichiometry matrix of a reaction system, computes a matrix of
conservation laws, each represented as a row in the output. 
"""
function conservationlaws(nsm::AbstractMatrix)::Matrix
    n_spec,n_reac = size(nsm)
    
    # We basically have to compute the left null space of the matrix
    # over the integers; this is best done using its Smith Normal Form.
    # note, we transpose as this was written when netstoichmat was reac by spec
    nsm_conv = AbstractAlgebra.matrix(AbstractAlgebra.ZZ, nsm')
    S, T, U = AbstractAlgebra.snf_with_transform(nsm_conv)
    
    # Zero columns of S (which occur after nonzero columns in SNF)
    # correspond to conserved quantities
    n = findfirst(i -> all(S[:,i] .== 0), 1:n_spec)
    if n === nothing
        return zeros(Int, 0, n_spec)
    end
    
    ret = Matrix(U[:,n:end]')
    
    # If all coefficients for a conservation law are negative
    # we might as well flip them to become positive
    for i in 1:size(ret,1)
        all(ret[i,:] .<= 0) && (ret[i,:] .*= -1)
    end
    
    ret
end

"""
    conservedquantities(state, cons_laws)

Compute conserved quantities for a system with the given conservation laws.
"""
conservedquantities(state, cons_laws) = cons_laws * state

######################## reaction network operators #######################

"""
    ==(rn1::Reaction, rn2::Reaction)

Tests whether two [`Reaction`](@ref)s are identical.

Notes:
- Ignores the order in which stoichiometry components are listed.
- *Does not* currently simplify rates, so a rate of `A^2+2*A+1` would be
    considered different than `(A+1)^2`.
"""
function (==)(rn1::Reaction, rn2::Reaction)
    isequal(rn1.rate, rn2.rate) || return false
    issetequal(zip(rn1.substrates,rn1.substoich), zip(rn2.substrates,rn2.substoich)) || return false
    issetequal(zip(rn1.products,rn1.prodstoich), zip(rn2.products,rn2.prodstoich)) || return false
    issetequal(rn1.netstoich, rn2.netstoich)
end


"""
    ==(rn1::ReactionSystem, rn2::ReactionSystem)

Tests whether the underlying species, parameters and reactions are the same in
the two [`ReactionSystem`](@ref)s. Ignores order network components were
defined.

Notes:
- *Does not* currently simplify rates, so a rate of `A^2+2*A+1` would be
    considered different than `(A+1)^2`.
- Flattens subsystems, and hence may allocate, when checking equality.
"""
function (==)(rn1::ReactionSystem, rn2::ReactionSystem)
    issetequal(species(rn1), species(rn2)) || return false
    issetequal(params(rn1), params(rn2)) || return false
    isequal(get_iv(rn1), get_iv(rn2)) || return false
    (numreactions(rn1) == numreactions(rn2)) || return false

    # the following fails for some reason, so need to use issubset
    #issetequal(equations(rn1), equations(rn2)) || return false
    (issubset(reactions(rn1),reactions(rn2)) && issubset(reactions(rn2),reactions(rn1))) || return false

    # BELOW SHOULD NOT BE NEEDED as species, params and equations flatten
    #issetequal(rn1.systems, rn2.systems) || return false
    # sys1 = rn1.systems; sys2 = rn2.systems
    # (issubset(sys1,sys2) && issubset(sys2,sys1)) || return false
    true
end


######################## functions to extend a network ####################

"""
    make_empty_network(; iv=DEFAULT_IV, name=gensym(:ReactionSystem))

Construct an empty [`ReactionSystem`](@ref). `iv` is the independent variable,
usually time, and `name` is the name to give the `ReactionSystem`.
"""
function make_empty_network(; iv=DEFAULT_IV, name=gensym(:ReactionSystem))
    ReactionSystem(Reaction[], iv, [], []; name=name)
end

"""
    addspecies!(network::ReactionSystem, s::Symbolic; disablechecks=false)

Given a [`ReactionSystem`](@ref), add the species corresponding to the variable
`s` to the network (if it is not already defined). Returns the integer id of the
species within the system.

Notes:
- `disablechecks` will disable checking for whether the passed in variable is
  already defined, which is useful when adding many new variables to the system.
  *Do not disable checks* unless you are sure the passed in variable is a new
  variable, as this will potentially leave the system in an undefined state.
"""
function addspecies!(network::ReactionSystem, s::Symbolic; disablechecks=false)

    # we don't check subsystems since we will add it to the top-level system...
    curidx = disablechecks ? nothing : findfirst(S -> isequal(S, s), get_states(network))
    if curidx === nothing
        push!(get_states(network), s)
        return length(get_states(network))
    else
        return curidx
    end
end


"""
    addspecies!(network::ReactionSystem, s::Num; disablechecks=false)

Given a [`ReactionSystem`](@ref), add the species corresponding to the
variable `s` to the network (if it is not already defined). Returns the
integer id of the species within the system.

- `disablechecks` will disable checking for whether the passed in variable is
  already defined, which is useful when adding many new variables to the system.
  *Do not disable checks* unless you are sure the passed in variable is a new
  variable, as this will potentially leave the system in an undefined state.
"""
function addspecies!(network::ReactionSystem, s::Num; disablechecks=false)
    addspecies!(network, value(s), disablechecks=disablechecks)
end

"""
    addparam!(network::ReactionSystem, p::Symbolic; disablechecks=false)

Given a [`ReactionSystem`](@ref), add the parameter corresponding to the
variable `p` to the network (if it is not already defined). Returns the integer
id of the parameter within the system.

- `disablechecks` will disable checking for whether the passed in variable is
  already defined, which is useful when adding many new variables to the system.
  *Do not disable checks* unless you are sure the passed in variable is a new
  variable, as this will potentially leave the system in an undefined state.
"""
function addparam!(network::ReactionSystem, p::Symbolic; disablechecks=false)
    # we don't check subsystems since we will add it to the top-level system...
    if istree(p) && !(operation(p) isa Sym)
        error("If the passed in parameter is an expression, it must correspond to an underlying Variable.")
    end
    curidx = disablechecks ? nothing : findfirst(S -> isequal(S, p), get_ps(network))
    if curidx === nothing
        push!(get_ps(network), p)
        return length(get_ps(network))
    else
        return curidx
    end
end

"""
    addparam!(network::ReactionSystem, p::Num; disablechecks=false)

Given a [`ReactionSystem`](@ref), add the parameter corresponding to the
variable `p` to the network (if it is not already defined). Returns the
integer id of the parameter within the system.

- `disablechecks` will disable checking for whether the passed in variable is
  already defined, which is useful when adding many new variables to the system.
  *Do not disable checks* unless you are sure the passed in variable is a new
  variable, as this will potentially leave the system in an undefined state.
"""
function addparam!(network::ReactionSystem, p::Num; disablechecks=false)
    addparam!(network, value(p); disablechecks=disablechecks)
end

"""
    addreaction!(network::ReactionSystem, rx::Reaction)

Add the passed in reaction to the [`ReactionSystem`](@ref). Returns the
integer id of `rx` in the list of `Reaction`s within `network`.

Notes:
- Any new species or parameters used in `rx` should be separately added to
    `network` using [`addspecies!`](@ref) and [`addparam!`](@ref).
"""
function addreaction!(network::ReactionSystem, rx::Reaction)
    push!(get_eqs(network), rx)
    length(get_eqs(network))
end


"""
    merge!(network1::ReactionSystem, network2::ReactionSystem)

Merge `network2` into `network1`.

Notes:
- Duplicate reactions between the two networks are not filtered out.
- [`Reaction`](@ref)s are not deepcopied to minimize allocations, so both
  networks will share underlying data arrays.
- Subsystems are not deepcopied between the two networks and will hence be
  shared.
- Returns `network1`.
- Does not currently handle pins.
"""
function Base.merge!(network1::ReactionSystem, network2::ReactionSystem)
    isequal(get_iv(network1), get_iv(network2)) || 
        error("Reaction networks must have the same independent variable to be mergable.")
    union!(get_states(network1), get_states(network2))
    union!(get_ps(network1), get_ps(network2))
    append!(get_eqs(network1), get_eqs(network2))
    append!(get_systems(network1), get_systems(network2))
    network1
end

"""
    merge(network1::ReactionSystem, network2::ReactionSystem)

Create a new network merging `network1` and `network2`.

Notes:
- Duplicate reactions between the two networks are not filtered out.
- [`Reaction`](@ref)s are not deepcopied to minimize allocations, so the new
  network will share underlying data arrays.
- Subsystems are not deepcopied between the two networks and will hence be
  shared.
- Returns the merged network.
- Does not currently handle pins.
"""
function Base.merge(network1::ReactionSystem, network2::ReactionSystem)
    network = make_empty_network()
    merge!(network, network1)
    merge!(network, network2)
    network
end


###############################   units   #####################################

"""
    validate(rx::Reaction; info::String = "")     

Check that all substrates and products within the given [`Reaction`](@ref) have
the same units, and that the units of the reaction's rate expression are
internally consistent (i.e. if the rate involves sums, each term in the sum has
the same units).

"""
function validate(rx::Reaction; info::String = "")     
    validated = ModelingToolkit._validate([rx.rate], [string(rx, ": rate")], info = info)
    
    subunits = isempty(rx.substrates) ? nothing : get_unit(rx.substrates[1])
    for i in 2:length(rx.substrates)
        if get_unit(rx.substrates[i]) != subunits
            validated = false
            @warn(string("In ", rx, " the substrates have differing units."))
        end
    end

    produnits = isempty(rx.products) ? nothing : get_unit(rx.products[1])
    for i in 2:length(rx.products)
        if get_unit(rx.products[i]) != produnits
            validated = false
            @warn(string("In ", rx, " the products have differing units."))
        end
    end

    if (subunits !== nothing) && (produnits !== nothing) && (subunits != produnits)
        validated = false
        @warn(string("in ", rx, " the substrate units are not consistent with the product units."))
    end

    validated
end

"""
    validate(rs::ReactionSystem, info::String="")

Check that all species in the [`ReactionSystem`](@ref) have the same units, and
that the rate laws of all reactions reduce to units of (species units) / (time
units).
"""
function validate(rs::ReactionSystem, info::String="")
    specs = get_states(rs)

    # if there are no species we don't check units on the system
    isempty(specs) && return true   

    specunits = get_unit(specs[1])
    validated = true
    for spec in specs
        if get_unit(spec) != specunits
            validated = false 
            @warn(string("Species are expected to have units of ", specunits, " however, species ", spec, " has units ", get_unit(spec), "."))
        end
    end
    timeunits = get_unit(get_iv(rs))
    rateunits = specunits / timeunits

    for rx in get_eqs(rs)
        rxunits = get_unit(rx.rate)
        for (i,sub) in enumerate(rx.substrates)
            rxunits *= get_unit(sub) ^ rx.substoich[i]
        end

        if rxunits != rateunits
            validated = false
            @warn(string("Reaction rate laws are expected to have units of ", rateunits, " however, ", rx, " has units of ", rxunits, "."))
        end
    end

    validated
end
