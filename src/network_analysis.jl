### Reaction Complex Handling ###

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

"""
    incidencemat(rn::ReactionSystem; sparse=false)

Calculate the incidence matrix of `rn`, see [`reactioncomplexes`](@ref).

Notes:
- Is cached in `rn` so that future calls, assuming the same sparsity, will also be fast.
"""
incidencemat(rn::ReactionSystem; sparse = false) = reactioncomplexes(rn; sparse)[2]

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

### Linkage, Deficiency, Reversibility ###

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

linkageclasses(incidencegraph) = Graphs.connected_components(incidencegraph)

"""
    stronglinkageclasses(rn::ReactionSystem)

    Return the strongly connected components of a reaction network's incidence graph (i.e. sub-groups of reaction complexes such that every complex is reachable from every other one in the sub-group).
"""

function stronglinkageclasses(rn::ReactionSystem)
    nps = get_networkproperties(rn)
    if isempty(nps.stronglinkageclasses)
        nps.stronglinkageclasses = stronglinkageclasses(incidencematgraph(rn))
    end
    nps.stronglinkageclasses
end

stronglinkageclasses(incidencegraph) = Graphs.strongly_connected_components(incidencegraph)

"""
    terminallinkageclasses(rn::ReactionSystem)

    Return the terminal strongly connected components of a reaction network's incidence graph (i.e. sub-groups of reaction complexes that are 1) strongly connected and 2) every outgoing reaction from a complex in the component produces a complex also in the component).
"""

function terminallinkageclasses(rn::ReactionSystem)
    nps = get_networkproperties(rn)
    if isempty(nps.terminallinkageclasses)
        slcs = stronglinkageclasses(rn)
        tslcs = filter(lc -> isterminal(lc, rn), slcs)
        nps.terminallinkageclasses = tslcs
    end
    nps.terminallinkageclasses
end

# Check whether a given linkage class in a reaction network is terminal, i.e. all outgoing reactions from complexes in the linkage class produce a complex also in hte linkage class
function isterminal(lc::Vector, rn::ReactionSystem)
    imat = incidencemat(rn)

    for r in 1:size(imat, 2)
        # Find the index of the reactant complex for a given reaction
        s = findfirst(==(-1), @view imat[:, r])

        # If the reactant complex is in the linkage class, check whether the product complex is also in the linkage class. If any of them are not, return false. 
        if s in Set(lc)
            p = findfirst(==(1), @view imat[:, r])
            p in Set(lc) ? continue : return false
        end
    end
    true
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

    # Check if deficiency has been computed already (initialized to -1)
    if nps.deficiency == -1
        conservationlaws(rn)
        r = nps.rank
        ig = incidencematgraph(rn)
        lc = linkageclasses(rn)
        nps.deficiency = Graphs.nv(ig) - length(lc) - r
    end
    nps.deficiency
end

# Used in the subsequent function.
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

### Conservation Laws ###

# Implements the `conserved` parameter metadata.
struct ConservedParameter end
Symbolics.option_to_metadata_type(::Val{:conserved}) = ConservedParameter

"""
isconserved(p)

Checks if the input parameter (`p`) is a conserved quantity (i.e. have the `conserved`)
metadata.
"""
isconserved(x::Num, args...) = isconserved(Symbolics.unwrap(x), args...)
function isconserved(x, default = false)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, ConservedParameter, default)
end

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

# Used in the subsequent function.
function cache_conservationlaw_eqs!(rn::ReactionSystem, N::AbstractMatrix, col_order)
    nullity = size(N, 1)
    r = numspecies(rn) - nullity     # rank of the netstoichmat
    sts = species(rn)
    indepidxs = col_order[begin:r]
    indepspecs = sts[indepidxs]
    depidxs = col_order[(r + 1):end]
    depspecs = sts[depidxs]
    constants = MT.unwrap.(MT.scalarize(only(
        @parameters $(CONSERVED_CONSTANT_SYMBOL)[1:nullity] [conserved = true])))

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

"""
    iscomplexbalanced(rs::ReactionSystem, parametermap)

Constructively compute whether a network will have complex-balanced equilibrium
solutions, following the method in van der Schaft et al., [2015](https://link.springer.com/article/10.1007/s10910-015-0498-2#Sec3). Accepts a dictionary, vector, or tuple of variable-to-value mappings, e.g. [k1 => 1.0, k2 => 2.0,...]. 
"""

function iscomplexbalanced(rs::ReactionSystem, parametermap::Dict)
    if length(parametermap) != numparams(rs)
        error("Incorrect number of parameters specified.")
    end

    pmap = symmap_to_varmap(rs, parametermap)
    pmap = Dict(ModelingToolkit.value(k) => v for (k, v) in pmap)

    sm = speciesmap(rs)
    cm = reactioncomplexmap(rs)
    complexes, D = reactioncomplexes(rs)
    rxns = reactions(rs)
    nc = length(complexes)
    nr = numreactions(rs)
    nm = numspecies(rs)

    if !all(r -> ismassaction(r, rs), rxns)
        error("The supplied ReactionSystem has reactions that are not ismassaction. Testing for being complex balanced is currently only supported for pure mass action networks.")
    end

    rates = [substitute(rate, pmap) for rate in reactionrates(rs)]

    # Construct kinetic matrix, K
    K = zeros(nr, nc)
    for c in 1:nc
        complex = complexes[c]
        for (r, dir) in cm[complex]
            rxn = rxns[r]
            if dir == -1
                K[r, c] = rates[r]
            end
        end
    end

    L = -D * K
    S = netstoichmat(rs)

    # Compute ρ using the matrix-tree theorem
    g = incidencematgraph(rs)
    R = ratematrix(rs, rates)
    ρ = matrixtree(g, R)

    # Determine if 1) ρ is positive and 2) D^T Ln ρ lies in the image of S^T
    if all(>(0), ρ)
        img = D' * log.(ρ)
        if rank(S') == rank(hcat(S', img))
            return true
        else
            return false
        end
    else
        return false
    end
end

function iscomplexbalanced(rs::ReactionSystem, parametermap::Vector{Pair{Symbol, Float64}})
    pdict = Dict(parametermap)
    iscomplexbalanced(rs, pdict)
end

function iscomplexbalanced(rs::ReactionSystem, parametermap::Tuple{Pair{Symbol, Float64}})
    pdict = Dict(parametermap)
    iscomplexbalanced(rs, pdict)
end

function iscomplexbalanced(rs::ReactionSystem, parametermap)
    error("Parameter map must be a dictionary, tuple, or vector of symbol/value pairs.")
end

"""
    ratematrix(rs::ReactionSystem, parametermap)

    Given a reaction system with n complexes, outputs an n-by-n matrix where R_{ij} is the rate constant of the reaction between complex i and complex j. Accepts a dictionary, vector, or tuple of variable-to-value mappings, e.g. [k1 => 1.0, k2 => 2.0,...]. 
"""

function ratematrix(rs::ReactionSystem, rates::Vector{Float64})
    complexes, D = reactioncomplexes(rs)
    n = length(complexes)
    rxns = reactions(rs)
    ratematrix = zeros(n, n)

    for r in 1:length(rxns)
        rxn = rxns[r]
        s = findfirst(==(-1), @view D[:, r])
        p = findfirst(==(1), @view D[:, r])
        ratematrix[s, p] = rates[r]
    end
    ratematrix
end

function ratematrix(rs::ReactionSystem, parametermap::Dict)
    if length(parametermap) != numparams(rs)
        error("Incorrect number of parameters specified.")
    end

    pmap = symmap_to_varmap(rs, parametermap)
    pmap = Dict(ModelingToolkit.value(k) => v for (k, v) in pmap)

    rates = [substitute(rate, pmap) for rate in reactionrates(rs)]
    ratematrix(rs, rates)
end

function ratematrix(rs::ReactionSystem, parametermap::Vector{Pair{Symbol, Float64}})
    pdict = Dict(parametermap)
    ratematrix(rs, pdict)
end

function ratematrix(rs::ReactionSystem, parametermap::Tuple{Pair{Symbol, Float64}})
    pdict = Dict(parametermap)
    ratematrix(rs, pdict)
end

function ratematrix(rs::ReactionSystem, parametermap)
    error("Parameter map must be a dictionary, tuple, or vector of symbol/value pairs.")
end

### BELOW: Helper functions for iscomplexbalanced

function matrixtree(g::SimpleDiGraph, distmx::Matrix)
    n = nv(g)
    if size(distmx) != (n, n)
        error("Size of distance matrix is incorrect")
    end

    π = zeros(n)

    if !Graphs.is_connected(g)
        ccs = Graphs.connected_components(g)
        for cc in ccs
            sg, vmap = Graphs.induced_subgraph(g, cc)
            distmx_s = distmx[cc, cc]
            π_j = matrixtree(sg, distmx_s)
            π[cc] = π_j
        end
        return π
    end

    # generate all spanning trees
    ug = SimpleGraph(SimpleDiGraph(g))
    trees = collect(Combinatorics.combinations(collect(edges(ug)), n - 1))
    trees = SimpleGraph.(trees)
    trees = filter!(t -> isempty(Graphs.cycle_basis(t)), trees)

    # constructed rooted trees for every vertex, compute sum
    for v in 1:n
        rootedTrees = [reverse(Graphs.bfs_tree(t, v, dir = :in)) for t in trees]
        π[v] = sum([treeweight(t, g, distmx) for t in rootedTrees])
    end

    # sum the contributions
    return π
end

function treeweight(t::SimpleDiGraph, g::SimpleDiGraph, distmx::Matrix)
    prod = 1
    for e in edges(t)
        s = Graphs.src(e)
        t = Graphs.dst(e)
        prod *= distmx[s, t]
    end
    prod
end

### Deficiency one

"""
    satisfiesdeficiencyone(rn::ReactionSystem)

    Check if a reaction network obeys the conditions of the deficiency one theorem, which ensures that there is only one equilibrium for every positive stoichiometric compatibility class.
"""

function satisfiesdeficiencyone(rn::ReactionSystem)
    complexes, D = reactioncomplexes(rn)
    δ = deficiency(rn)
    δ_l = linkagedeficiencies(rn)

    lcs = linkageclasses(rn)
    tslcs = terminallinkageclasses(rn)

    # Check the conditions for the deficiency one theorem: 
    #   1) the deficiency of each individual linkage class is at most 1; 
    #   2) the sum of the linkage deficiencies is the total deficiency, and 
    #   3) there is only one terminal linkage class per linkage class. 
    all(<=(1), δ_l) && (sum(δ_l) == δ) && (length(lcs) == length(tslcs))
end

"""
    robustspecies(rn::ReactionSystem)

    Return a vector of indices corresponding to species that are concentration robust, i.e. for every positive equilbrium, the concentration of species s will be the same. 
"""

function robustspecies(rn::ReactionSystem)
    complexes, D = reactioncomplexes(rn)
    nps = get_networkproperties(rn)

    if deficiency(rn) != 1
        error("This algorithm currently only checks for robust species in networks with deficiency one.")
    end

    # A species is concnetration robust in a deficiency one network if there are two non-terminal complexes (i.e. complexes 
    # belonging to a linkage class that is not terminal) that differ only in species s (i.e. their difference is some 
    # multiple of s. (A + B, A) differ only in B. (A + 2B, B) differ in both A and B, since A + 2B - B = A + B). 
    
    if !nps.checkedrobust
        tslcs = terminallinkageclasses(rn)
        Z = complexstoichmat(rn)

        # Find the complexes that do not belong to a terminal linkage class 
        nonterminal_complexes = deleteat!([1:length(complexes);], vcat(tslcs...))
        robust_species = Int64[]

        for (c_s, c_p) in collect(Combinatorics.combinations(nonterminal_complexes, 2))
            # Check the difference of all the combinations of complexes. The support is the set of indices that are non-zero 
            supp = findall(!=(0), Z[:, c_s] - Z[:, c_p])

            # If the support has length one, then they differ in only one species, and that species is concentration robust. 
            length(supp) == 1 && supp[1] ∉ robust_species && push!(robust_species, supp...)
        end
        nps.checkedrobust = true
        nps.robustspecies = robust_species
    end

    nps.robustspecies
end

"""
    isconcentrationrobust(rn::ReactionSystem, species::Int)

    Given a reaction network and an index of a species, check if that species is concentration robust. 
"""

function isconcentrationrobust(rn::ReactionSystem, species::Int)
    robust_species = robustspecies(rn)
    return species in robust_species
end
