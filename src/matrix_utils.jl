check_square(mtx::AbstractMatrix, label::AbstractString = "Matrix") =
    (size(mtx, 1) == size(mtx, 2)) ||
        throw(DimensionMismatch("$label needs to be square, $(size(mtx)) found"))

# check if the value `a` should be indexed by sortedvalues()/indexvalues()
isindexed(a, test::EdgeTest{T}) where T =
    (isnothing(skipval(test)) || (a != skipval(test))) &&
    (isnothing(test.threshold) || !isstronger(a, test.threshold, rev=isreverse(test)))

function sortedvalues!(res::AbstractVector{T}, A::AbstractArray{T},
                       test::EdgeTest{T} = EdgeTest{T}();
                       alg=Base.Sort.defalg(A)) where T
    empty!(res)
    isempty(A) && return res
    sizehint!(res, isnothing(skipval(test)) ? length(A) : sum(w -> w != skipval(test), A))
    @inbounds for a in A
        isindexed(a, test) && push!(res, a)
    end
    unique!(sort!(res, alg=alg, rev=isreverse(test)))
    return res
end

sortedvalues(A::AbstractArray{T}, test::EdgeTest{T} = EdgeTest{T}();
             kwargs...) where T =
    sortedvalues!(Vector{T}(), A, test; kwargs...)

function indexvalues!(iA::AbstractMatrix{I}, weights::AbstractVector{T},
                      A::AbstractArray{T}, test::EdgeTest{T} = EdgeTest{T}();
                      kwargs...) where {T, I <: Integer}
    sortedvalues!(weights, A, test; kwargs...)
    weightdict = Dict{T, I}(val => i for (i, val) in enumerate(weights))
    # convert adjmtx to weight indices. higher index=stronger edge
    iA = fill!(reshape(resize!(vec(iA), length(A)), size(A)), 0)
    @inbounds for (i, a) in enumerate(A)
        isindexed(a, test) && (iA[i] = weightdict[a]) #searchsortedfirst(weights, a, rev=rev))
    end
    return iA, weights
end

indexvalues(::Type{I}, A::AbstractArray{T}, test::EdgeTest{T} = EdgeTest{T}();
            kwargs...) where {T, I <: Integer} =
    indexvalues!(similar(A, I), Vector{T}(), A, test; kwargs...)

empty_adjacencymatrix(::Type{M}, ::Type{T} = eltype(M)) where {M <: AbstractMatrix, T} =
    Matrix{T}(undef, (0, 0))

subgraph_adjacencymatrix(adjmtx::AbstractMatrix, comp_indices::AbstractVector{<:Integer},
    pool::Union{ArrayPool, Nothing} = nothing) =
    view(adjmtx, comp_indices, comp_indices)

"""
"Condenses" the matrix `A` by aggregating the values in the blocks of its
elements defined by `row_groups` and `col_groups`.

# Arguments
* `rev::Bool`, defaults to `false`: if false, the aggregated value is the
  maximal value of the block (i.e. the edge with the highest weight). Otherwise,
  it's the minimal value (the edge with the smallest weight).
"""
condense(A::AbstractMatrix{T}, node_groups::AbstractPartition,
         test::EdgeTest = EdgeTest{T}()) where T =
    condense!(similar(A, length(node_groups), length(node_groups)), A,
              node_groups, test)

condense(A::AbstractMatrix{T},
         row_groups::AbstractPartition, col_groups::AbstractPartition,
         test::EdgeTest{T} = EdgeTest{T}();
         zerodiag::Bool = false) where T =
  condense!(similar(A, length(row_groups), length(col_groups)), A,
            row_groups, col_groups, test, zerodiag=zerodiag)

"""
"Condenses" the matrix `A` by aggregating the values in the blocks of its
elements defined by `row_groups` and `col_groups`.

# Arguments
* `rev::Bool`, defaults to `false`: if false, the aggregated value is the
  maximal value of the block (i.e. the edge with the highest weight). Otherwise,
  it's the minimal value (the edge with the smallest weight).
"""
function condense!(B::AbstractMatrix{T}, A::AbstractMatrix{T},
                   row_groups::AbstractPartition, col_groups::AbstractPartition,
                   test::EdgeTest{T} = EdgeTest{T}();
                   zerodiag::Bool = false) where T
    nrows = nelems(row_groups)
    ncols = nelems(col_groups)
    size(A) == (nrows, ncols) ||
        throw(DimensionMismatch("A size ($(size(A))) and row/col labels sizes ($nrows, $ncols) do not match."))
    size(B) == (length(row_groups), length(col_groups)) ||
        throw(DimensionMismatch("B size ($(size(B))) and row/col group number ($(length(row_groups)), $(length(col_groups))) do not match."))
    @inbounds for (jj, cols) in enumerate(col_groups)
        Bjj = view(B, :, jj)
        for j in cols
            Aj = view(A, :, j)
            for (ii, rows) in enumerate(row_groups)
                firstcol = j == first(cols)
                if zerodiag && (ii == jj)
                    firstcol && (Bjj[ii] = zero(T))
                elseif firstcol && (length(rows) == 1)
                    Bjj[ii] = Aj[first(rows)]
                else
                    Bij = ifelse(firstcol, defaultweight(test), Bjj[ii])
                    @simd for i in rows
                        w = Aj[i]
                        Bij = ifelse((isnothing(skipval(test)) || (w != skipval(test))) &&
                                     ((!isnothing(skipval(test)) && (Bij == skipval(test))) ||
                                      isweaker(Bij, w, rev=isreverse(test))),
                                     w, Bij)
                    end
                    Bjj[ii] = Bij
                end
            end
        end
    end
    return B
end

condense!(B::AbstractMatrix, A::AbstractMatrix{T},
          node_groups::AbstractPartition,
          test::EdgeTest{T} = EdgeTest{T}();
          zerodiag::Bool = false) where T =
    condense!(B, A, node_groups, node_groups, test, zerodiag=zerodiag)
