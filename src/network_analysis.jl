



######### Accessors: #########

"""
    dependents(rx, network)

Given a [`Reaction`](@ref) and a [`ReactionSystem`](@ref), return a vector of the
*non-constant* species and variables the reaction rate law depends on. e.g., for

`k*W, 2X + 3Y --> 5Z + W`

the returned vector would be `[W(t),X(t),Y(t)]`.

Notes:
- Allocates
- Does not check for dependents within any subsystems.
- Constant species are not considered dependents since they are internally treated as
  parameters.
- If the rate expression depends on a non-species unknown variable that will be included in
  the dependents, i.e. in
  ```julia
  t = default_t()
  @parameters k
  @variables V(t)
  @species A(t) B(t) C(t)
  rx = Reaction(k*V, [A, B], [C])
  @named rs = ReactionSystem([rx], t)
  issetequal(dependents(rx, rs), [A,B,V]) == true
  ```
"""
function dependents(rx, network)
    if rx.rate isa Number
        return rx.substrates
    else
        rvars = get_variables(rx.rate, unknowns(network))
        return union!(rvars, rx.substrates)
    end
end

"""
    dependents(rx, network)

See documentation for [`dependents`](@ref).
"""
function dependants(rx, network)
    dependents(rx, network)
end

"""
    reactionrates(network)

Given a [`ReactionSystem`](@ref), returns a vector of the symbolic reaction
rates for each reaction.
"""
function reactionrates(rn)
    [r.rate for r in reactions(rn)]
end

"""
    substoichmat(rn; sparse=false)

Returns the substrate stoichiometry matrix, ``S``, with ``S_{i j}`` the stoichiometric
coefficient of the ith substrate within the jth reaction.

Note:
- Set sparse=true for a sparse matrix representation
- Note that constant species are not considered substrates, but just components that modify
  the associated rate law.
"""
function substoichmat(::Type{SparseMatrixCSC{T, Int}},
                      rn::ReactionSystem) where {T <: Number}
    Is = Int[]
    Js = Int[]
    Vs = T[]
    smap = speciesmap(rn)
    for (k, rx) in enumerate(reactions(rn))
        stoich = rx.substoich
        for (i, sub) in enumerate(rx.substrates)
            isconstant(sub) && continue
            push!(Js, k)
            push!(Is, smap[sub])
            push!(Vs, stoich[i])
        end
    end
    sparse(Is, Js, Vs, numspecies(rn), numreactions(rn))
end
function substoichmat(::Type{Matrix{T}}, rn::ReactionSystem) where {T <: Number}
    smap = speciesmap(rn)
    smat = zeros(T, numspecies(rn), numreactions(rn))
    for (k, rx) in enumerate(reactions(rn))
        stoich = rx.substoich
        for (i, sub) in enumerate(rx.substrates)
            isconstant(sub) && continue
            smat[smap[sub], k] = stoich[i]
        end
    end
    smat
end
function substoichmat(rn::ReactionSystem; sparse::Bool = false)
    isempty(get_systems(rn)) || error("substoichmat does not currently support subsystems.")
    T = reduce(promote_type, eltype(rx.substoich) for rx in reactions(rn))
    (T == Any) &&
        error("Stoichiometry matrices with symbolic stoichiometry are not supported")
    sparse ? substoichmat(SparseMatrixCSC{T, Int}, rn) : substoichmat(Matrix{T}, rn)
end

"""
    prodstoichmat(rn; sparse=false)

Returns the product stoichiometry matrix, ``P``, with ``P_{i j}`` the stoichiometric
coefficient of the ith product within the jth reaction.

Note:
- Set sparse=true for a sparse matrix representation
- Note that constant species are not treated as products, but just components that modify
  the associated rate law.
"""
function prodstoichmat(::Type{SparseMatrixCSC{T, Int}},
                       rn::ReactionSystem) where {T <: Number}
    Is = Int[]
    Js = Int[]
    Vs = T[]
    smap = speciesmap(rn)
    for (k, rx) in enumerate(reactions(rn))
        stoich = rx.prodstoich
        for (i, prod) in enumerate(rx.products)
            isconstant(prod) && continue
            push!(Js, k)
            push!(Is, smap[prod])
            push!(Vs, stoich[i])
        end
    end
    sparse(Is, Js, Vs, numspecies(rn), numreactions(rn))
end
function prodstoichmat(::Type{Matrix{T}}, rn::ReactionSystem) where {T <: Number}
    smap = speciesmap(rn)
    pmat = zeros(T, numspecies(rn), numreactions(rn))
    for (k, rx) in enumerate(reactions(rn))
        stoich = rx.prodstoich
        for (i, prod) in enumerate(rx.products)
            isconstant(prod) && continue
            pmat[smap[prod], k] = stoich[i]
        end
    end
    pmat
end
function prodstoichmat(rn::ReactionSystem; sparse = false)
    isempty(get_systems(rn)) ||
        error("prodstoichmat does not currently support subsystems.")

    T = reduce(promote_type, eltype(rx.prodstoich) for rx in reactions(rn))
    (T == Any) &&
        error("Stoichiometry matrices with symbolic stoichiometry are not supported")
    sparse ? prodstoichmat(SparseMatrixCSC{T, Int}, rn) : prodstoichmat(Matrix{T}, rn)
end

"""
    netstoichmat(rn, sparse=false)

Returns the net stoichiometry matrix, ``N``, with ``N_{i j}`` the net stoichiometric
coefficient of the ith species within the jth reaction.

Notes:
- Set sparse=true for a sparse matrix representation
- Caches the matrix internally within `rn` so subsequent calls are fast.
- Note that constant species are not treated as reactants, but just components that modify
  the associated rate law. As such they do not contribute to the net stoichiometry matrix.
"""
function netstoichmat(::Type{SparseMatrixCSC{T, Int}},
                      rn::ReactionSystem) where {T <: Number}
    Is = Int[]
    Js = Int[]
    Vs = Vector{T}()
    smap = speciesmap(rn)
    for (k, rx) in pairs(reactions(rn))
        for (spec, coef) in rx.netstoich
            isconstant(spec) && continue
            push!(Js, k)
            push!(Is, smap[spec])
            push!(Vs, coef)
        end
    end
    sparse(Is, Js, Vs, numspecies(rn), numreactions(rn))
end
function netstoichmat(::Type{Matrix{T}}, rn::ReactionSystem) where {T <: Number}
    smap = speciesmap(rn)
    nmat = zeros(T, numspecies(rn), numreactions(rn))
    for (k, rx) in pairs(reactions(rn))
        for (spec, coef) in rx.netstoich
            isconstant(spec) && continue
            nmat[smap[spec], k] = coef
        end
    end
    nmat
end

netstoichtype(::Vector{Pair{S, T}}) where {S, T} = T

function netstoichmat(rn::ReactionSystem; sparse = false)
    isempty(get_systems(rn)) ||
        error("netstoichmat does not currently support subsystems, please create a flattened system before calling.")

    nps = get_networkproperties(rn)

    # if it is already calculated and has the right type
    !isempty(nps.netstoichmat) && (sparse == issparse(nps.netstoichmat)) &&
        (return nps.netstoichmat)

    # identify a common stoichiometry type
    T = reduce(promote_type, netstoichtype(rx.netstoich) for rx in reactions(rn))
    (T == Any) &&
        error("Stoichiometry matrices are not supported with symbolic stoichiometry.")

    if sparse
        nsmat = netstoichmat(SparseMatrixCSC{T, Int}, rn)
    else
        nsmat = netstoichmat(Matrix{T}, rn)
    end

    # only cache if it is integer
    if T == Int
        nps.netstoichmat = nsmat
    end

    nsmat
end

# the following function is adapted from SymbolicUtils.jl v.19
# later on (Spetember 2023) modified by Torkel and Shashi (now assumes input not on polynomial form, which is handled elsewhere, previous version failed in these cases anyway).
# Copyright (c) 2020: Shashi Gowda, Yingbo Ma, Mason Protter, Julia Computing.
# MIT license
"""
    to_multivariate_poly(polyeqs::AbstractVector{BasicSymbolic{Real}})

Convert the given system of polynomial equations to multivariate polynomial representation.
For example, this can be used in HomotopyContinuation.jl functions.
"""
function to_multivariate_poly(polyeqs::AbstractVector{Symbolics.BasicSymbolic{Real}})
    @assert length(polyeqs)>=1 "At least one expression must be passed to `multivariate_poly`."

    pvar2sym, sym2term = SymbolicUtils.get_pvar2sym(), SymbolicUtils.get_sym2term()
    ps = map(polyeqs) do x
        if istree(x) && operation(x) == (/)
            error("We should not be able to get here, please contact the package authors.")
        else
            PolyForm(x, pvar2sym, sym2term).p
        end
    end

    ps
end
"""
    setdefaults!(rn, newdefs)

Sets the default (initial) values of parameters and species in the
`ReactionSystem`, `rn`.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
setdefaults!(sir, [:S => 999.0, :I => 1.0, :R => 1.0, :β => 1e-4, :ν => .01])

# or
t = default_t()
@parameter β ν
@species S(t) I(t) R(t)
setdefaults!(sir, [S => 999.0, I => 1.0, R => 0.0, β => 1e-4, ν => .01])
```
gives initial/default values to each of `S`, `I` and `β`

Notes:
- Can not be used to set default values for species, variables or parameters of
  subsystems or constraint systems. Either set defaults for those systems
  directly, or [`flatten`](@ref) to collate them into one system before setting
  defaults.
- Defaults can be specified in any iterable container of symbols to value pairs
  or symbolics to value pairs.
"""
function setdefaults!(rn, newdefs)
    defs = eltype(newdefs) <: Pair{Symbol} ? symmap_to_varmap(rn, newdefs) : newdefs
    rndefs = MT.get_defaults(rn)
    for (var, val) in defs
        rndefs[value(var)] = value(val)
    end
    nothing
end

function __unpacksys(rn)
    ex = :(begin end)
    for key in keys(get_var_to_name(rn))
        var = MT.getproperty(rn, key, namespace = false)
        push!(ex.args, :($key = $var))
    end
    ex
end

"""
    @unpacksys sys::ModelingToolkit.AbstractSystem

Loads all species, variables, parameters, and observables defined in `sys` as
variables within the calling module.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
@unpacksys sir
```
will load the symbolic variables, `S`, `I`, `R`, `ν` and `β`.

Notes:
- Can not be used to load species, variables, or parameters of subsystems or
  constraints. Either call `@unpacksys` on those systems directly, or
  [`flatten`](@ref) to collate them into one system before calling.
- Note that this places symbolic variables within the calling module's scope, so
  calling from a function defined in a script or the REPL will still result in
  the symbolic variables being defined in the `Main` module.
"""
macro unpacksys(rn)
    quote
        ex = Catalyst.__unpacksys($(esc(rn)))
        Base.eval($(__module__), ex)
    end
end

# convert symbol of the form :sys.a.b.c to a symbolic a.b.c
function _symbol_to_var(sys, sym)
    if hasproperty(sys, sym)
        var = getproperty(sys, sym, namespace = false)
    else
        strs = split(String(sym), "₊")   # need to check if this should be split of not!!!
        if length(strs) > 1
            var = getproperty(sys, Symbol(strs[1]), namespace = false)
            for str in view(strs, 2:length(strs))
                var = getproperty(var, Symbol(str), namespace = true)
            end
        else
            throw(ArgumentError("System $(nameof(sys)): variable $sym does not exist"))
        end
    end
    var
end

"""
    symmap_to_varmap(sys, symmap)

Given a system and map of `Symbol`s to values, generates a map from
corresponding symbolic variables/parameters to the values that can be used to
pass initial conditions and parameter mappings.

For example,
```julia
sir = @reaction_network sir begin
    β, S + I --> 2I
    ν, I --> R
end
subsys = @reaction_network subsys begin
    k, A --> B
end
@named sys = compose(sir, [subsys])
```
gives
```
Model sys with 3 equations
Unknowns (5):
  S(t)
  I(t)
  R(t)
  subsys₊A(t)
  subsys₊B(t)
Parameters (3):
  β
  ν
  subsys₊k
```
to specify initial condition and parameter mappings from *symbols* we can use
```julia
symmap = [:S => 1.0, :I => 1.0, :R => 1.0, :subsys₊A => 1.0, :subsys₊B => 1.0]
u0map  = symmap_to_varmap(sys, symmap)
pmap   = symmap_to_varmap(sys, [:β => 1.0, :ν => 1.0, :subsys₊k => 1.0])
```
`u0map` and `pmap` can then be used as input to various problem types.

Notes:
- Any `Symbol`, `sym`, within `symmap` must be a valid field of `sys`. i.e.
  `sys.sym` must be defined.
"""
function symmap_to_varmap(sys, symmap::Tuple)
    if all(p -> p isa Pair{Symbol}, symmap)
        return ((_symbol_to_var(sys, sym) => val for (sym, val) in symmap)...,)
    else  # if not all entries map a symbol to value pass through
        return symmap
    end
end

function symmap_to_varmap(sys, symmap::AbstractArray{Pair{Symbol, T}}) where {T}
    [_symbol_to_var(sys, sym) => val for (sym, val) in symmap]
end

function symmap_to_varmap(sys, symmap::Dict{Symbol, T}) where {T}
    Dict(_symbol_to_var(sys, sym) => val for (sym, val) in symmap)
end

# don't permute any other types and let varmap_to_vars handle erroring
symmap_to_varmap(sys, symmap) = symmap
#error("symmap_to_varmap requires a Dict, AbstractArray or Tuple to map Symbols to values.")

######################## reaction complexes and reaction rates ###############################

"""
    reset_networkproperties!(rn::ReactionSystem)

Clears the cache of various properties (like the netstoichiometry matrix). Use if such
properties need to be recalculated for some reason.
"""
function reset_networkproperties!(rn::ReactionSystem)
    reset!(get_networkproperties(rn))
    nothing
end

# get the species indices and stoichiometry while filtering out constant species.
function filter_constspecs(specs, stoich::AbstractVector{V}, smap) where {V <: Integer}
    isempty(specs) && (return Vector{Int}(), Vector{V}())

    if any(isconstant, specs)
        ids = Vector{Int}()
        filtered_stoich = Vector{V}()
        for (i, s) in enumerate(specs)
            if !isconstant(s)
                push!(ids, smap[s])
                push!(filtered_stoich, stoich[i])
            end
        end
    else
        ids = map(Base.Fix1(getindex, smap), specs)
        filtered_stoich = copy(stoich)
    end
    ids, filtered_stoich
end

"""
    reactioncomplexmap(rn::ReactionSystem)

Find each [`ReactionComplex`](@ref) within the specified system, constructing a mapping from
the complex to vectors that indicate which reactions it appears in as substrates and
products.

Notes:
- Each [`ReactionComplex`](@ref) is mapped to a vector of pairs, with each pair having the
  form `reactionidx => ± 1`, where `-1` indicates the complex appears as a substrate and
  `+1` as a product in the reaction with integer label `reactionidx`.
- Constant species are ignored as part of a complex. i.e. if species `A` is constant then
  the reaction `A + B --> C + D` is considered to consist of the complexes `B` and `C + D`.
  Likewise `A --> B` would be treated as the same as `0 --> B`.
"""
function reactioncomplexmap(rn::ReactionSystem)
    isempty(get_systems(rn)) ||
        error("reactioncomplexmap does not currently support subsystems.")

    # check if previously calculated and hence cached
    nps = get_networkproperties(rn)
    !isempty(nps.complextorxsmap) && return nps.complextorxsmap
    complextorxsmap = nps.complextorxsmap

    rxs = reactions(rn)
    smap = speciesmap(rn)
    numreactions(rn) > 0 ||
        error("There must be at least one reaction to find reaction complexes.")
    for (i, rx) in enumerate(rxs)
        subids, substoich = filter_constspecs(rx.substrates, rx.substoich, smap)
        subrc = sort!(ReactionComplex(subids, substoich))
        if haskey(complextorxsmap, subrc)
            push!(complextorxsmap[subrc], i => -1)
        else
            complextorxsmap[subrc] = [i => -1]
        end

        prodids, prodstoich = filter_constspecs(rx.products, rx.prodstoich, smap)
        prodrc = sort!(ReactionComplex(prodids, prodstoich))
        if haskey(complextorxsmap, prodrc)
            push!(complextorxsmap[prodrc], i => 1)
        else
            complextorxsmap[prodrc] = [i => 1]
        end
    end
    complextorxsmap
end

function reactioncomplexes(::Type{SparseMatrixCSC{Int, Int}}, rn::ReactionSystem,
                           complextorxsmap)
    complexes = collect(keys(complextorxsmap))
    Is = Int[]
    Js = Int[]
    Vs = Int[]
    for (i, c) in enumerate(complexes)
        for (j, σ) in complextorxsmap[c]
            push!(Is, i)
            push!(Js, j)
            push!(Vs, σ)
        end
    end
    B = sparse(Is, Js, Vs, length(complexes), numreactions(rn))
    complexes, B
end
function reactioncomplexes(::Type{Matrix{Int}}, rn::ReactionSystem, complextorxsmap)
    complexes = collect(keys(complextorxsmap))
    B = zeros(Int, length(complexes), numreactions(rn))
    for (i, c) in enumerate(complexes)
        for (j, σ) in complextorxsmap[c]
            B[i, j] = σ
        end
    end
    complexes, B
end

@doc raw"""
    reactioncomplexes(network::ReactionSystem; sparse=false)

Calculate the reaction complexes and complex incidence matrix for the given
[`ReactionSystem`](@ref).

Notes:
- returns a pair of a vector of [`ReactionComplex`](@ref)s and the complex incidence matrix.
- An empty [`ReactionComplex`](@ref) denotes the null (∅) state (from reactions like ∅ -> A
  or A -> ∅).
- Constant species are ignored in generating a reaction complex. i.e. if A is constant then
  A --> B consists of the complexes ∅ and B.
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
function reactioncomplexes(rn::ReactionSystem; sparse = false)
    isempty(get_systems(rn)) ||
        error("reactioncomplexes does not currently support subsystems.")
    nps = get_networkproperties(rn)
    if isempty(nps.complexes) || (sparse != issparse(nps.complexes))
        complextorxsmap = reactioncomplexmap(rn)
        nps.complexes, nps.incidencemat = if sparse
            reactioncomplexes(SparseMatrixCSC{Int, Int}, rn, complextorxsmap)
        else
            reactioncomplexes(Matrix{Int}, rn, complextorxsmap)
        end
    end
    nps.complexes, nps.incidencemat
end

"""
    incidencemat(rn::ReactionSystem; sparse=false)

Calculate the incidence matrix of `rn`, see [`reactioncomplexes`](@ref).

Notes:
- Is cached in `rn` so that future calls, assuming the same sparsity, will also be fast.
"""
incidencemat(rn::ReactionSystem; sparse = false) = reactioncomplexes(rn; sparse)[2]

function complexstoichmat(::Type{SparseMatrixCSC{Int, Int}}, rn::ReactionSystem, rcs)
    Is = Int[]
    Js = Int[]
    Vs = Int[]
    for (i, rc) in enumerate(rcs)
        for rcel in rc
            push!(Is, rcel.speciesid)
            push!(Js, i)
            push!(Vs, rcel.speciesstoich)
        end
    end
    Z = sparse(Is, Js, Vs, numspecies(rn), length(rcs))
end
function complexstoichmat(::Type{Matrix{Int}}, rn::ReactionSystem, rcs)
    Z = zeros(Int, numspecies(rn), length(rcs))
    for (i, rc) in enumerate(rcs)
        for rcel in rc
            Z[rcel.speciesid, i] = rcel.speciesstoich
        end
    end
    Z
end

"""
    complexstoichmat(network::ReactionSystem; sparse=false)

Given a [`ReactionSystem`](@ref) and vector of reaction complexes, return a
matrix with positive entries of size number of species by number of complexes,
where the non-zero positive entries in the kth column denote stoichiometric
coefficients of the species participating in the kth reaction complex.

Notes:
- Set sparse=true for a sparse matrix representation
"""
function complexstoichmat(rn::ReactionSystem; sparse = false)
    isempty(get_systems(rn)) ||
        error("complexstoichmat does not currently support subsystems.")
    nps = get_networkproperties(rn)
    if isempty(nps.complexstoichmat) || (sparse != issparse(nps.complexstoichmat))
        nps.complexstoichmat = if sparse
            complexstoichmat(SparseMatrixCSC{Int, Int}, rn, keys(reactioncomplexmap(rn)))
        else
            complexstoichmat(Matrix{Int}, rn, keys(reactioncomplexmap(rn)))
        end
    end
    nps.complexstoichmat
end

function complexoutgoingmat(::Type{SparseMatrixCSC{Int, Int}}, rn::ReactionSystem, B)
    n = size(B, 2)
    rows = rowvals(B)
    vals = nonzeros(B)
    Is = Int[]
    Js = Int[]
    Vs = Int[]
    sizehint!(Is, div(length(vals), 2))
    sizehint!(Js, div(length(vals), 2))
    sizehint!(Vs, div(length(vals), 2))
    for j in 1:n
        for i in nzrange(B, j)
            if vals[i] != one(eltype(vals))
                push!(Is, rows[i])
                push!(Js, j)
                push!(Vs, vals[i])
            end
        end
    end
    sparse(Is, Js, Vs, size(B, 1), size(B, 2))
end
function complexoutgoingmat(::Type{Matrix{Int}}, rn::ReactionSystem, B)
    Δ = copy(B)
    for (I, b) in pairs(Δ)
        (b == 1) && (Δ[I] = 0)
    end
    Δ
end

@doc raw"""
    complexoutgoingmat(network::ReactionSystem; sparse=false)

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
function complexoutgoingmat(rn::ReactionSystem; sparse = false)
    isempty(get_systems(rn)) ||
        error("complexoutgoingmat does not currently support subsystems.")
    nps = get_networkproperties(rn)
    if isempty(nps.complexoutgoingmat) || (sparse != issparse(nps.complexoutgoingmat))
        B = reactioncomplexes(rn, sparse = sparse)[2]
        nps.complexoutgoingmat = if sparse
            complexoutgoingmat(SparseMatrixCSC{Int, Int}, rn, B)
        else
            complexoutgoingmat(Matrix{Int}, rn, B)
        end
    end
    nps.complexoutgoingmat
end

function incidencematgraph(incidencemat::Matrix{Int})
    @assert all(∈([-1, 0, 1]), incidencemat)
    n = size(incidencemat, 1)  # no. of nodes/complexes
    graph = Graphs.DiGraph(n)
    for col in eachcol(incidencemat)
        src = 0
        dst = 0
        for i in eachindex(col)
            (col[i] == -1) && (src = i)
            (col[i] == 1) && (dst = i)
            (src != 0) && (dst != 0) && break
        end
        Graphs.add_edge!(graph, src, dst)
    end
    return graph
end
function incidencematgraph(incidencemat::SparseMatrixCSC{Int, Int})
    @assert all(∈([-1, 0, 1]), incidencemat)
    m, n = size(incidencemat)
    graph = Graphs.DiGraph(m)
    rows = rowvals(incidencemat)
    vals = nonzeros(incidencemat)
    for j in 1:n
        inds = nzrange(incidencemat, j)
        row = rows[inds]
        val = vals[inds]
        if val[1] == -1
            Graphs.add_edge!(graph, row[1], row[2])
        else
            Graphs.add_edge!(graph, row[2], row[1])
        end
    end
    return graph
end

"""
    incidencematgraph(rn::ReactionSystem)

Construct a directed simple graph where nodes correspond to reaction complexes and directed
edges to reactions converting between two complexes.

Notes:
- Requires the `incidencemat` to already be cached in `rn` by a previous call to
  `reactioncomplexes`.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
complexes,incidencemat = reactioncomplexes(sir)
incidencematgraph(sir)
```
"""
function incidencematgraph(rn::ReactionSystem)
    nps = get_networkproperties(rn)
    if Graphs.nv(nps.incidencegraph) == 0
        isempty(nps.incidencemat) &&
            error("Please call reactioncomplexes(rn) first to construct the incidence matrix.")
        nps.incidencegraph = incidencematgraph(nps.incidencemat)
    end
    nps.incidencegraph
end

linkageclasses(incidencegraph) = Graphs.connected_components(incidencegraph)

"""
    linkageclasses(rn::ReactionSystem)

Given the incidence graph of a reaction network, return a vector of the
connected components of the graph (i.e. sub-groups of reaction complexes that
are connected in the incidence graph).

Notes:
- Requires the `incidencemat` to already be cached in `rn` by a previous call to
  `reactioncomplexes`.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
complexes,incidencemat = reactioncomplexes(sir)
linkageclasses(sir)
```
gives
```julia
2-element Vector{Vector{Int64}}:
 [1, 2]
 [3, 4]
```
"""
function linkageclasses(rn::ReactionSystem)
    nps = get_networkproperties(rn)
    if isempty(nps.linkageclasses)
        nps.linkageclasses = linkageclasses(incidencematgraph(rn))
    end
    nps.linkageclasses
end

@doc raw"""
    deficiency(rn::ReactionSystem)

Calculate the deficiency of a reaction network.

Here the deficiency, ``\delta``, of a network with ``n`` reaction complexes,
``\ell`` linkage classes and a rank ``s`` stoichiometric matrix is

```math
\delta = n - \ell - s
```

Notes:
- Requires the `incidencemat` to already be cached in `rn` by a previous call to
  `reactioncomplexes`.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
rcs,incidencemat = reactioncomplexes(sir)
δ = deficiency(sir)
```
"""
function deficiency(rn::ReactionSystem)
    nps = get_networkproperties(rn)
    conservationlaws(rn)
    r = nps.rank
    ig = incidencematgraph(rn)
    lc = linkageclasses(rn)
    nps.deficiency = Graphs.nv(ig) - length(lc) - r
    nps.deficiency
end

function subnetworkmapping(linkageclass, allrxs, complextorxsmap, p)
    rxinds = sort!(collect(Set(rxidx for rcidx in linkageclass
                               for rxidx in complextorxsmap[rcidx])))
    rxs = allrxs[rxinds]
    specset = Set(s for rx in rxs for s in rx.substrates if !isconstant(s))
    for rx in rxs
        for product in rx.products
            !isconstant(product) && push!(specset, product)
        end
    end
    specs = collect(specset)
    newps = Vector{eltype(p)}()
    for rx in rxs
        Symbolics.get_variables!(newps, rx.rate, p)
    end
    rxs, specs, newps   # reactions and species involved in reactions of subnetwork
end

"""
    subnetworks(rn::ReactionSystem)

Find subnetworks corresponding to each linkage class of the reaction network.

Notes:
- Requires the `incidencemat` to already be cached in `rn` by a previous call to
  `reactioncomplexes`.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
complexes,incidencemat = reactioncomplexes(sir)
subnetworks(sir)
```
"""
function subnetworks(rs::ReactionSystem)
    isempty(get_systems(rs)) || error("subnetworks does not currently support subsystems.")
    lcs = linkageclasses(rs)
    rxs = reactions(rs)
    p = parameters(rs)
    t = get_iv(rs)
    spatial_ivs = get_sivs(rs)
    complextorxsmap = [map(first, rcmap) for rcmap in values(reactioncomplexmap(rs))]
    subnetworks = Vector{ReactionSystem}()
    for i in 1:length(lcs)
        reacs, specs, newps = subnetworkmapping(lcs[i], rxs, complextorxsmap, p)
        newname = Symbol(nameof(rs), "_", i)
        push!(subnetworks,
              ReactionSystem(reacs, t, specs, newps; name = newname, spatial_ivs))
    end
    subnetworks
end

"""
    linkagedeficiencies(network::ReactionSystem)

Calculates the deficiency of each sub-reaction network within `network`.

Notes:
- Requires the `incidencemat` to already be cached in `rn` by a previous call to
  `reactioncomplexes`.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
rcs,incidencemat = reactioncomplexes(sir)
linkage_deficiencies = linkagedeficiencies(sir)
```
"""
function linkagedeficiencies(rs::ReactionSystem)
    lcs = linkageclasses(rs)
    subnets = subnetworks(rs)
    δ = zeros(Int, length(lcs))
    for (i, subnet) in enumerate(subnets)
        conservationlaws(subnet)
        nps = get_networkproperties(subnet)
        δ[i] = length(lcs[i]) - 1 - nps.rank
    end
    δ
end

"""
    isreversible(rn::ReactionSystem)

Given a reaction network, returns if the network is reversible or not.

Notes:
- Requires the `incidencemat` to already be cached in `rn` by a previous call to
  `reactioncomplexes`.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
rcs,incidencemat = reactioncomplexes(sir)
isreversible(sir)
```
"""
function isreversible(rn::ReactionSystem)
    ig = incidencematgraph(rn)
    Graphs.reverse(ig) == ig
end

"""
    isweaklyreversible(rn::ReactionSystem, subnetworks)

Determine if the reaction network with the given subnetworks is weakly reversible or not.

Notes:
- Requires the `incidencemat` to already be cached in `rn` by a previous call to
  `reactioncomplexes`.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
rcs,incidencemat = reactioncomplexes(sir)
subnets = subnetworks(rn)
isweaklyreversible(rn, subnets)
```
"""
function isweaklyreversible(rn::ReactionSystem, subnets)
    im = get_networkproperties(rn).incidencemat
    isempty(im) &&
        error("Error, please call reactioncomplexes(rn::ReactionSystem) to ensure the incidence matrix has been cached.")
    sparseig = issparse(im)
    for subnet in subnets
        nps = get_networkproperties(subnet)
        isempty(nps.incidencemat) && reactioncomplexes(subnet; sparse = sparseig)
    end
    all(Graphs.is_strongly_connected ∘ incidencematgraph, subnets)
end

############################################################################################
######################## conservation laws ###############################

"""
    conservedequations(rn::ReactionSystem)

Calculate symbolic equations from conservation laws, writing dependent variables as
functions of independent variables and the conservation law constants.

Notes:
- Caches the resulting equations in `rn`, so will be fast on subsequent calls.

Examples:
```@repl
rn = @reaction_network begin
    k, A + B --> C
    k2, C --> A + B
    end
conservedequations(rn)
```
gives
```
2-element Vector{Equation}:
 B(t) ~ A(t) + Γ[1]
 C(t) ~ Γ[2] - A(t)
```
"""
function conservedequations(rn::ReactionSystem)
    conservationlaws(rn)
    nps = get_networkproperties(rn)
    nps.conservedeqs
end

"""
    conservationlaw_constants(rn::ReactionSystem)

Calculate symbolic equations from conservation laws, writing the conservation law constants
in terms of the dependent and independent variables.

Notes:
- Caches the resulting equations in `rn`, so will be fast on subsequent calls.

Examples:
```@julia
rn = @reaction_network begin
    k, A + B --> C
    k2, C --> A + B
    end
conservationlaw_constants(rn)
```
gives
```
2-element Vector{Equation}:
 Γ[1] ~ B(t) - A(t)
 Γ[2] ~ A(t) + C(t)
```
"""
function conservationlaw_constants(rn::ReactionSystem)
    conservationlaws(rn)
    nps = get_networkproperties(rn)
    nps.constantdefs
end

"""
    conservationlaws(netstoichmat::AbstractMatrix)::Matrix

Given the net stoichiometry matrix of a reaction system, computes a matrix of
conservation laws, each represented as a row in the output.
"""
function conservationlaws(nsm::T; col_order = nothing) where {T <: AbstractMatrix}

    # compute the left nullspace over the integers
    N = MT.nullspace(nsm'; col_order)

    # if all coefficients for a conservation law are negative, make positive
    for Nrow in eachcol(N)
        all(r -> r <= 0, Nrow) && (Nrow .*= -1)
    end

    # check we haven't overflowed
    iszero(N' * nsm) || error("Calculation of the conservation law matrix was inaccurate, "
          * "likely due to numerical overflow. Please use a larger integer "
          * "type like Int128 or BigInt for the net stoichiometry matrix.")

    T(N')
end

function cache_conservationlaw_eqs!(rn::ReactionSystem, N::AbstractMatrix, col_order)
    nullity = size(N, 1)
    r = numspecies(rn) - nullity     # rank of the netstoichmat
    sts = species(rn)
    indepidxs = col_order[begin:r]
    indepspecs = sts[indepidxs]
    depidxs = col_order[(r + 1):end]
    depspecs = sts[depidxs]
    constants = MT.unwrap.(MT.scalarize((@parameters Γ[1:nullity])[1]))

    conservedeqs = Equation[]
    constantdefs = Equation[]
    for (i, depidx) in enumerate(depidxs)
        scaleby = (N[i, depidx] != 1) ? N[i, depidx] : one(eltype(N))
        (scaleby != 0) || error("Error, found a zero in the conservation law matrix where "
              *
              "one was not expected.")
        coefs = @view N[i, indepidxs]
        terms = sum(p -> p[1] / scaleby * p[2], zip(coefs, indepspecs))
        eq = depspecs[i] ~ constants[i] - terms
        push!(conservedeqs, eq)
        eq = constants[i] ~ depspecs[i] + terms
        push!(constantdefs, eq)
    end

    # cache in the system
    nps = get_networkproperties(rn)
    nps.rank = r
    nps.nullity = nullity
    nps.indepspecs = Set(indepspecs)
    nps.depspecs = Set(depspecs)
    nps.conservedeqs = conservedeqs
    nps.constantdefs = constantdefs

    nothing
end

"""
    conservationlaws(rs::ReactionSystem)

Return the conservation law matrix of the given `ReactionSystem`, calculating it if it is
not already stored within the system, or returning an alias to it.

Notes:
- The first time being called it is calculated and cached in `rn`, subsequent calls should
  be fast.
"""
function conservationlaws(rs::ReactionSystem)
    nps = get_networkproperties(rs)
    !isempty(nps.conservationmat) && (return nps.conservationmat)
    nsm = netstoichmat(rs)
    nps.conservationmat = conservationlaws(nsm; col_order = nps.col_order)
    cache_conservationlaw_eqs!(rs, nps.conservationmat, nps.col_order)
    nps.conservationmat
end

"""
    conservedquantities(state, cons_laws)

Compute conserved quantities for a system with the given conservation laws.
"""
conservedquantities(state, cons_laws) = cons_laws * state


# If u0s are not given while conservation laws are present, throws an error.
# Used in HomotopyContinuation and BifurcationKit extensions.
# Currently only checks if any u0s are given
# (not whether these are enough for computing conserved quantitites, this will yield a less informative error).
function conservationlaw_errorcheck(rs, pre_varmap)
    vars_with_vals = Set(p[1] for p in pre_varmap)
    any(s -> s in vars_with_vals, species(rs)) && return
    isempty(conservedequations(Catalyst.flatten(rs))) ||
        error("The system has conservation laws but initial conditions were not provided for some species.")
end
