function sortedvalues!(res::AbstractVector{T}, A::AbstractArray{T};
                       skipval::Union{Number, Nothing} = zero(eltype(A)),
                       threshold::Union{Number, Nothing}=nothing,
                       alg=Base.Sort.defalg(A),
                       rev::Bool=false) where T
    empty!(res)
    isempty(A) && return res
    sizehint!(res, isnothing(skipval) ? length(A) : sum(w -> w != skipval, A))
    for a in A
        if (isnothing(skipval) || (a != skipval)) &&
           (isnothing(threshold) || !isstronger(a, threshold, rev=rev))
            push!(res, a)
        end
    end
    unique!(sort!(res, alg=alg, rev=rev))
    return res
end

sortedvalues(A::AbstractArray{T}; kwargs...) where T =
    sortedvalues!(Vector{T}(), A; kwargs...)

function indexvalues!(iA::AbstractMatrix{I}, weights::AbstractVector{T},
                      A::AbstractArray{T};
                      skipval::Union{Number, Nothing} = zero(eltype(A)),
                      kwargs...) where {T, I <: Integer}
    sortedvalues!(weights, A; skipval=skipval, kwargs...)
    weightdict = Dict{T, I}(val => i for (i, val) in enumerate(weights))
    # convert adjmtx to weight indices. higher index=stronger edge
    iA = fill!(reshape(resize!(vec(iA), length(A)), size(A)), 0)
    @inbounds for (i, a) in enumerate(A)
        if isnothing(skipval) || a != skipval
            iA[i] = weightdict[a]#searchsortedfirst(weights, a, rev=rev)
        end
    end
    return iA, weights
end

indexvalues(::Type{I}, A::AbstractArray{T}; kwargs...) where {T, I <: Integer} =
    indexvalues!(similar(A, I), Vector{T}(), A; kwargs...)

subgraph_adjacencymatrix(adjmtx::AbstractMatrix, comp_indices::AbstractVector{<:Integer}) =
    view(adjmtx, comp_indices, comp_indices)

"""
"Condenses" the matrix `A` by aggregating the values in the blocks of its
elements defined by `row_groups` and `col_groups`.

# Arguments
* `rev::Bool`, defaults to `false`: if false, the aggregated value is the
  maximal value of the block (i.e. the edge with the highest weight). Otherwise,
  it's the minimal value (the edge with the smallest weight).
"""
condense(A::AbstractMatrix, node_groups::AbstractPartition; kwargs...) =
    condense!(similar(A, length(node_groups), length(node_groups)), A,
              node_groups; kwargs...)

condense(A::AbstractMatrix,
         row_groups::AbstractPartition, col_groups::AbstractPartition;
         kwargs...) =
  condense!(similar(A, length(row_groups), length(col_groups)), A,
            row_groups, col_groups; kwargs...)

"""
"Condenses" the matrix `A` by aggregating the values in the blocks of its
elements defined by `row_groups` and `col_groups`.

# Arguments
* `rev::Bool`, defaults to `false`: if false, the aggregated value is the
  maximal value of the block (i.e. the edge with the highest weight). Otherwise,
  it's the minimal value (the edge with the smallest weight).
"""
function condense!(B::AbstractMatrix,
                   A::AbstractMatrix,
                   row_groups::AbstractPartition,
                   col_groups::AbstractPartition = row_groups;
                   skipval::Union{Number, Nothing} = zero(eltype(A)),
                   rev::Bool=false)
    nrows = nelems(row_groups)
    ncols = nelems(col_groups)
    size(A) == (nrows, ncols) ||
        throw(DimensionMismatch("A size ($(size(A))) and row/col labels sizes ($nrows, $ncols) do not match."))
    size(B) == (length(row_groups), length(col_groups)) ||
        throw(DimensionMismatch("B size ($(size(B))) and row/col group number ($(length(row_groups)), $(length(col_groups))) do not match."))
    fill!(B, defaultweight(eltype(A), skipval=skipval, rev=rev))
    @inbounds for (jj, cols) in enumerate(col_groups)
        B_j = view(B, :, jj)
        for j in cols
            Aj = view(A, :, j)
            for (ii, rows) in enumerate(row_groups)
                Bij = B_j[ii]
                for i in rows
                    w = Aj[i]
                    !isnothing(skipval) && (w == skipval) && continue
                    if (!isnothing(skipval) && (Bij == skipval)) || isweaker(Bij, w, rev=rev)
                        Bij = w
                    end
                end
                B_j[ii] = Bij
            end
        end
    end
    return B
end
