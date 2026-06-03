# This file contain the code required for MTK's `nullspace` function. The code is fetched from a
# MTKv9.80.1 (MIT-licensed). Later copies of MTK have moved this functionality to the AGPL license.
# A single line have been modified (to remove a `ModelingToolkit.`). This is marked with a
# comment at the appropriate line. The code consists of the full "bareiss.jl" file, and a single function
# from the "alias_elimination.jl" file. The code here needed for various network analysis
# functionality, and was removed from MTK around the v10/v11 changes.

### The "structural_transformation/bareiss.jl" file ###

using LinearAlgebra
using SparseArrays
using SparseArrays: AbstractSparseMatrixCSC, getcolptr

macro swap(a, b)
    return esc(:(($a, $b) = ($b, $a)))
end

import Base: swaprows!

function bareiss_update!(
        zero!, M::StridedMatrix, k, swapto, pivot,
        prev_pivot::Base.BitInteger
    )
    flag = zero(prev_pivot)
    prev_pivot = Base.MultiplicativeInverses.SignedMultiplicativeInverse(prev_pivot)
    @inbounds for i in (k + 1):size(M, 2)
        Mki = M[k, i]
        @simd ivdep for j in (k + 1):size(M, 1)
            M[j, i], r = divrem(M[j, i] * pivot - M[j, k] * Mki, prev_pivot)
            flag = flag | r
        end
    end
    iszero(flag) || error("Overflow occurred")
    return zero!(M, (k + 1):size(M, 1), k)
end

function bareiss_update!(zero!, M::StridedMatrix, k, swapto, pivot, prev_pivot)
    @inbounds for i in (k + 1):size(M, 2), j in (k + 1):size(M, 1)
        M[j, i] = exactdiv(M[j, i] * pivot - M[j, k] * M[k, i], prev_pivot)
    end
    return zero!(M, (k + 1):size(M, 1), k)
end

@views function bareiss_update!(zero!, M::AbstractMatrix, k, swapto, pivot, prev_pivot)
    if prev_pivot isa Base.BitInteger
        prev_pivot = Base.MultiplicativeInverses.SignedMultiplicativeInverse(prev_pivot)
    end
    V = M[(k + 1):end, (k + 1):end]
    V .= exactdiv.(V .* pivot .- M[(k + 1):end, k] * M[k, (k + 1):end]', prev_pivot)
    zero!(M, (k + 1):size(M, 1), k)
    if M isa AbstractSparseMatrixCSC
        dropzeros!(M)
    end
end

function bareiss_update_virtual_colswap!(
        zero!, M::AbstractMatrix, k, swapto, pivot,
        prev_pivot
    )
    if prev_pivot isa Base.BitInteger
        prev_pivot = Base.MultiplicativeInverses.SignedMultiplicativeInverse(prev_pivot)
    end
    V = @view M[(k + 1):end, :]
    V .= @views exactdiv.(V .* pivot .- M[(k + 1):end, swapto[2]] * M[k, :]', prev_pivot)
    return zero!(M, (k + 1):size(M, 1), swapto[2])
end

bareiss_zero!(M, i, j) = M[i, j] .= zero(eltype(M))

function find_pivot_col(M, i)
    p = findfirst(!iszero, @view M[i, i:end])
    p === nothing && return nothing
    idx = CartesianIndex(i, p + i - 1)
    return (idx, M[idx])
end

function find_pivot_any(M, i)
    p = findfirst(!iszero, @view M[i:end, i:end])
    p === nothing && return nothing
    idx = p + CartesianIndex(i - 1, i - 1)
    return (idx, M[idx])
end

const bareiss_colswap = (Base.swapcols!, swaprows!, bareiss_update!, bareiss_zero!)
const bareiss_virtcolswap = (
    (M, i, j) -> nothing, swaprows!,
    bareiss_update_virtual_colswap!, bareiss_zero!,
)

"""
    bareiss!(M, [swap_strategy])

Perform Bareiss's fraction-free row-reduction algorithm on the matrix `M`.
Optionally, a specific pivoting method may be specified.

swap_strategy is an optional argument that determines how the swapping of rows and columns is performed.
bareiss_colswap (the default) swaps the columns and rows normally.
bareiss_virtcolswap pretends to swap the columns which can be faster for sparse matrices.
"""
function bareiss!(
        M::AbstractMatrix{T}, swap_strategy = bareiss_colswap;
        find_pivot = find_pivot_any, column_pivots = nothing
    ) where {T}
    swapcols!, swaprows!, update!, zero! = swap_strategy
    prev = one(eltype(M))
    n = size(M, 1)
    pivot = one(T)
    column_permuted = false
    for k in 1:n
        r = find_pivot(M, k)
        r === nothing && return (k - 1, pivot, column_permuted)
        (swapto, pivot) = r
        if column_pivots !== nothing && k != swapto[2]
            column_pivots[k] = swapto[2]
            column_permuted |= true
        end
        if CartesianIndex(k, k) != swapto
            swapcols!(M, k, swapto[2])
            swaprows!(M, k, swapto[1])
        end
        update!(zero!, M, k, swapto, pivot, prev)
        prev = pivot
    end
    return (n, pivot, column_permuted)
end

function nullspace(A; col_order = nothing)
    n = size(A, 2)
    workspace = zeros(Int, 2 * n)
    column_pivots = @view workspace[1:n]
    pivots_cache = @view workspace[(n + 1):(2n)]
    @inbounds for i in 1:n
        column_pivots[i] = i
    end
    B = copy(A)
    (rank, d, column_permuted) = bareiss!(B; column_pivots)
    reduce_echelon!(B, rank, d, pivots_cache)

    # The first rank entries in col_order are columns that give a basis
    # for the column space. The remainder give the free variables.
    if col_order !== nothing
        resize!(col_order, size(A, 2))
        col_order .= 1:size(A, 2)
        for (i, cp) in enumerate(column_pivots)
            @swap(col_order[i], col_order[cp])
        end
    end

    fill!(pivots_cache, 0)

    # Modified: Was `N = ModelingToolkit.reduced_echelon_nullspace(rank, B, pivots_cache)`.
    N = reduced_echelon_nullspace(rank, B, pivots_cache)
    return apply_inv_pivot_rows!(N, column_pivots)
end

function apply_inv_pivot_rows!(M, ipiv)
    for i in size(M, 1):-1:1
        swaprows!(M, i, ipiv[i])
    end
    return M
end

###
### Modified from AbstractAlgebra.jl
###
### https://github.com/Nemocas/AbstractAlgebra.jl/blob/4803548c7a945f3f7bd8c63f8bb7c79fac92b11a/LICENSE.md
function reduce_echelon!(
        A::AbstractMatrix{T}, rank, d,
        pivots_cache = zeros(Int, size(A, 2))
    ) where {T}
    m, n = size(A)
    isreduced = true
    @inbounds for i in 1:rank
        for j in 1:(i - 1)
            if A[j, i] != zero(T)
                isreduced = false
                @goto out
            end
        end
        if A[i, i] != one(T)
            isreduced = false
            @goto out
        end
    end
    @label out
    @inbounds for i in (rank + 1):m, j in 1:n
        A[i, j] = zero(T)
    end
    isreduced && return A

    @inbounds if rank > 1
        t = zero(T)
        q = zero(T)
        d = -d
        pivots = pivots_cache
        np = rank
        j = k = 1
        for i in 1:rank
            while iszero(A[i, j])
                pivots[np + k] = j
                j += 1
                k += 1
            end
            pivots[i] = j
            j += 1
        end
        while k <= n - rank
            pivots[np + k] = j
            j += 1
            k += 1
        end
        for k in 1:(n - rank)
            for i in (rank - 1):-1:1
                t = A[i, pivots[np + k]] * d
                for j in (i + 1):rank
                    t += A[i, pivots[j]] * A[j, pivots[np + k]] + q
                end
                A[i, pivots[np + k]] = exactdiv(-t, A[i, pivots[i]])
            end
        end
        d = -d
        for i in 1:rank
            for j in 1:rank
                if i == j
                    A[j, pivots[i]] = d
                else
                    A[j, pivots[i]] = zero(T)
                end
            end
        end
    end
    return A
end

function reduced_echelon_nullspace(
        rank, A::AbstractMatrix{T},
        pivots_cache = zeros(Int, size(A, 2))
    ) where {T}
    n = size(A, 2)
    nullity = n - rank
    U = zeros(T, n, nullity)
    @inbounds if rank == 0
        for i in 1:nullity
            U[i, i] = one(T)
        end
    elseif nullity != 0
        pivots = @view pivots_cache[1:rank]
        nonpivots = @view pivots_cache[(rank + 1):n]
        j = k = 1
        for i in 1:rank
            while iszero(A[i, j])
                nonpivots[k] = j
                j += 1
                k += 1
            end
            pivots[i] = j
            j += 1
        end
        while k <= nullity
            nonpivots[k] = j
            j += 1
            k += 1
        end
        d = -A[1, pivots[1]]
        for i in 1:nullity
            for j in 1:rank
                U[pivots[j], i] = A[j, nonpivots[i]]
            end
            U[nonpivots[i], i] = d
        end
    end
    return U
end


### From the "systems/alias_elimination.jl file ###

function exactdiv(a::Integer, b)
    d, r = divrem(a, b)
    @assert r == 0
    return d
end
