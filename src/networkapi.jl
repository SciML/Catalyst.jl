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
    substoichmat(rn; smap=speciesmap(rn))

Returns the substrate stoichiometry matrix
"""
function substoichmat(rn; smap=speciesmap(rn))
    smat = zeros(Int,(numreactions(rn),numspecies(rn)))
    for (k,rx) in enumerate(reactions(rn))
        stoich = rx.substoich
        for (i,sub) in enumerate(rx.substrates)
            smat[k,smap[sub]] = stoich[i]
        end
    end
    smat
end

"""
    prodstoichmat(rn; smap=speciesmap(rn))

Returns the product stoichiometry matrix
"""
function prodstoichmat(rn; smap=speciesmap(rn))
    pmat = zeros(Int,(numreactions(rn),numspecies(rn)))
    for (k,rx) in enumerate(reactions(rn))
        stoich = rx.prodstoich
        for (i,prod) in enumerate(rx.products)
            pmat[k,smap[prod]] = stoich[i]
        end
    end
    pmat
end


"""
    netstoichmat(rn; smap=speciesmap(rn))

Returns the net stoichiometry matrix
"""
function netstoichmat(rn; smap=speciesmap(rn))
    nmat = zeros(Int,(numreactions(rn),numspecies(rn)))
    for (k,rx) in pairs(reactions(rn))
        for (spec,coef) in rx.netstoich
            nmat[k,smap[spec]] = coef
        end
    end
    nmat
end


######################## some reacion complexes matrices and reaction rates ###############################
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
    reactioncomplexes(network, smap=speciesmap(rn))

Calculate the reaction complexes and complex incidence matrix for the given [`ReactionSystem`](@ref). 

Notes:
- returns a pair of a vector of [`ReactionComplex`](@ref)s and the complex incidence matrix.
- An empty [`ReactionComplex`](@ref) denotes the null (∅) state (from reactions like ∅ -> A or A -> ∅).
- The complex incidence matrix, B, is number of complexes by number of reactions with
    Bᵢⱼ = -1, if the i'th complex is the substrate of the j'th reaction,
           1, if the i'th complex is the product of the j'th reaction,
           0, otherwise
"""
function reactioncomplexes(rn; smap=speciesmap(rn))
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
    
    complexes = collect(keys(complextorxsmap))
    B = zeros(Int64, length(complexes), numreactions(rn));
    for (i,c) in enumerate(complexes)
        for (j,σ) in complextorxsmap[c]
            B[i,j] = σ
        end
    end
    complexes,B
end

"""
    reaction_rates(network)

Given a [`ReactionSystem`](@ref), returns a vector of the symbolic reaction rates for each reaction.
"""
function reactionrates(rn)
    return [r.rate for r in reactions(rn)]
end


"""
    complexstoichmat(network; rcs=reactioncomplexes(rn)[1]))

Given a [`ReactionSystem`](@ref) and vector of reaction complexes, return a
matrix with positive entries of size num_of_species x num_of_complexes, where
the non-zero positive entries in the kth column denote stoichiometric
coefficients of the species participating in the kth reaction complex.
"""
function complexstoichmat(rn; rcs=reactioncomplexes(rn)[1])
    Z = zeros(Int64, numspecies(rn), length(rcs));
    for (i,rc) in enumerate(rcs)
        for rcel in rc
            Z[rcel.speciesid,i] = rcel.speciesstoich
        end
    end
    return Z
end

"""
    complexoutgoingmat(network; B=reactioncomplexes(rn)[2])

Given a [`ReactionSystem`](@ref) and complex incidence matrix, B, return a matrix
of size num_of_complexes x num_of_reactions.

Notes
- The complex outgoing matrix, Δ, is defined by 
    Δᵢⱼ = 0,    if Bᵢⱼ = 1 
    Δᵢⱼ = Bᵢⱼ,  otherwise
"""
function complexoutgoingmat(rn; B=reactioncomplexes(rn)[2])
    Δ = copy(B)
    Δ[Δ .== 1] .= 0
    return Δ
end


################################################################################################
######################## conservation laws ###############################

""" 
    conservationlaws(netstoichmat::AbstractMatrix)::Matrix

Given the net stoichiometry matrix of a reaction system, computes a matrix
of conservation laws, each represented as a row in the output. 
"""
function conservationlaws(nsm::AbstractMatrix)::Matrix
    n_reac, n_spec = size(nsm)
    
    # We basically have to compute the left null space of the matrix
    # over the integers; this is best done using its Smith Normal Form.
    nsm_conv = AbstractAlgebra.matrix(AbstractAlgebra.ZZ, nsm)
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
