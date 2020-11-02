# Functions for querying network properties.

######### Accessors: #########

"""
    species(network)

Given a [`ReactionSystem`](@ref), return a vector of species `Variable`s.

Notes:
- If `network.systems` is not empty, may allocate. Otherwise returns
  `network.states`.
"""
function species(network)
    isempty(network.systems) ? network.states : states(network)
end

"""
    params(network)

Given a [`ReactionSystem`](@ref), return a vector of parameter `Variable`s.

Notes:
- If `network.systems` is not empty, may allocate. Otherwise returns
  `network.ps`.
"""
function params(network)
    isempty(network.systems) ? network.ps : parameters(network)
end

"""
    reactions(network)

Given a [`ReactionSystem`](@ref), return a vector of all `Reactions` in the system.

Notes:
- If `network.systems` is not empty, may allocate. Otherwise returns
  `network.eqs`.
"""
function reactions(network)
    isempty(network.systems) ? network.eqs : equations(network)
end

"""
    speciesmap(network)

Given a [`ReactionSystem`](@ref), return a Dictionary mapping from species
`Variable`s to species indices. (Allocates)
"""
function speciesmap(network)
    Dict(S => i for (i,S) in enumerate(species(network)))
end

"""
    paramsmap(network)

Given a [`ReactionSystem`](@ref), return a Dictionary mapping from parameter
`Variable`s to parameter indices. (Allocates)
"""
function paramsmap(network)
    Dict(p => i for (i,p) in enumerate(params(network)))
end

"""
    numspecies(network)

Return the number of species within the given [`ReactionSystem`](@ref).
"""
function numspecies(network)
    ns = length(network.states)
    for sys in network.systems
        ns += numspecies(ns)
    end
    ns
end

"""
    numreactions(network)

Return the number of reactions within the given [`ReactionSystem`](@ref).
"""
function numreactions(network)
    nr = length(network.eqs)
    for sys in network.systems
        nr += numreactions(sys)
    end
    nr
end

"""
    numparams(network)

Return the number of parameters within the given [`ReactionSystem`](@ref).
"""
function numparams(network)
    np = length(network.ps)
    for sys in network.systems
        np += numparams(ns)
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
        rvars = ModelingToolkit.get_variables(rx.rate, species(network))
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

Tests whether the underlying species `Variables`s, parameter `Variables`s
and reactions are the same in the two [`ReactionSystem`](@ref)s. Ignores
order network components were defined.

Notes:
- *Does not* currently simplify rates, so a rate of `A^2+2*A+1` would be
    considered different than `(A+1)^2`.
- Flattens subsystems, and hence may allocate, when checking equality.
"""
function (==)(rn1::ReactionSystem, rn2::ReactionSystem)
    issetequal(species(rn1), species(rn2)) || return false
    issetequal(params(rn1), params(rn2)) || return false
    isequal(rn1.iv, rn2.iv) || return false
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
    make_empty_network(; iv=Sym{ModelingToolkit.Parameter{Real}}(:t))

Construct an empty [`ReactionSystem`](@ref). `iv` is the independent variable, usually time.
"""
function make_empty_network(; iv=Sym{ModelingToolkit.Parameter{Real}}(:t))
    ReactionSystem(Reaction[], iv, Num[], Num[], Sym[], Equation[], gensym(:ReactionSystem), ReactionSystem[])
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
    curidx = disablechecks ? nothing : findfirst(S -> isequal(S, s), network.states)
    if curidx === nothing
        push!(network.states, s)
        return length(network.states)
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
    if p isa Term && !(p.op isa Sym)
        error("If the passed in parameter is an expression, it must correspond to an underlying Variable.")
    end
    curidx = disablechecks ? nothing : findfirst(S -> isequal(S, p), network.ps)
    if curidx === nothing
        push!(network.ps, p)
        return length(network.ps)
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
    push!(network.eqs, rx)
    length(network.eqs)
end


"""
    merge!(network1::ReactionSystem, network2::ReactionSystem)

Merge `network2` into `network1`.

Notes:
- Duplicate reactions between the two networks are not filtered out.
- [`Reaction`](@ref)s are not deepcopied to minimize allocations, so both networks will share underlying data arrays.
- Subsystems are not deepcopied between the two networks and will hence be shared.
- Returns `network1`.
"""
function merge!(network1::ReactionSystem, network2::ReactionSystem)
    isequal(network1.iv, network2.iv) || error("Reaction networks must have the same independent variable to be mergable.")
    union!(network1.states, network2.states)
    union!(network1.ps, network2.ps)

    append!(network1.eqs, network2.eqs)
    append!(network1.systems, network2.systems)

    network1
end

"""
    merge(network1::ReactionSystem, network2::ReactionSystem)

Create a new network merging `network1` and `network2`.

Notes:
- Duplicate reactions between the two networks are not filtered out.
- [`Reaction`](@ref)s are not deepcopied to minimize allocations, so the new network will share underlying data arrays.
- Subsystems are not deepcopied between the two networks and will hence be shared.
- Returns the merged network.
"""
function merge(network1::ReactionSystem, network2::ReactionSystem)
    network = make_empty_network()
    merge!(network, network1)
    merge!(network, network2)
    network
end
