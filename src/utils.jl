@inline isweaker(a::Number, b::Number; rev::Bool=false) =
    (!rev && (a < b)) || (rev && (a > b))
@inline isstronger(a::Number, b::Number; rev::Bool=false) =
    (!rev && (a > b)) || (rev && (a < b))

@inline isvalidedge(w::Number;
            skipval::Union{Number, Nothing} = zero(w),
            threshold::Union{Number, Nothing} = nothing,
            rev::Bool = false) =
    (isnothing(skipval) || (w != skipval)) &&
    (isnothing(threshold) || !isweaker(w, threshold, rev=rev))

@inline defaultweight(::Type{T};
              skipval::Union{Number, Nothing} = zero(T),
              rev::Bool = false) where T =
    isnothing(skipval) ? (rev ? typemax(T) : typemin(T)) : skipval


function valuegroupsizes(maxvalue::Int, values::AbstractVector{Int})
    sizes = fill(0, maxvalue)
    for (i, val) in enumerate(values)
        sizes[val] += 1
    end
    return sizes
end

valuegroupsizes(values::AbstractVector{Int}) = valuegroupsizes(maximum(values), values)

function valuegroups(maxvalue::Int, values::AbstractVector{Int})
    groupstarts = push!(valuegroupsizes(maxvalue, values), 0)
    curstart = 1
    @inbounds for i in eachindex(groupstarts)
        groupsize = groupstarts[i]
        groupstarts[i] = curstart
        curstart += groupsize
    end
    poses = copy(groupstarts)
    order = Vector{Int}(undef, length(values))
    @inbounds for (i, val) in enumerate(values)
        order[poses[val]] = i
        poses[val] += 1
    end
    return IndicesPartition(order, groupstarts)
end

valuegroups(values::AbstractVector{Int}) = valuegroups(maximum(values), values)

function insertsorted!(arr::AbstractArray, item::Any;
                       rev::Bool=false, unique::Bool=false)
    ix = searchsortedfirst(arr, item, rev=rev)
    if ix > length(arr)
        push!(arr, item)
    elseif !unique || (@inbounds(arr[ix]) != item)
        insert!(arr, ix, item)
    end
    return arr
end

function sortedvalues(A::AbstractArray{T};
                      skipval::Union{Number, Nothing} = zero(eltype(A)),
                      threshold::Union{Number, Nothing}=nothing,
                      alg=Base.Sort.defalg(A),
                      rev::Bool=false) where T
    isempty(A) && return Vector{T}()
    res = Vector{T}()
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

"""
"Condenses" the matrix `A` by aggregating the values in the blocks of its
elements defined by `row_groups` and `col_groups`.

# Arguments
* `rev::Bool`, defaults to `false`: if false, the aggregated value is the
  maximal value of the block (i.e. the edge with the highest weight). Otherwise,
  it's the minimal value (the edge with the smallest weight).
"""
function condense(A::AbstractMatrix,
                  row_groups::AbstractPartition,
                  col_groups::AbstractPartition = row_groups;
                  skipval::Union{Number, Nothing} = zero(eltype(A)),
                  rev::Bool=false)
    nrows = nelems(row_groups)
    ncols = nelems(col_groups)
    size(A) == (nrows, ncols) ||
        throw(DimensionMismatch("A size ($(size(A))) and row/col labels sizes ($nrows, $ncols) do not match."))
    B = fill!(similar(A, length(row_groups), length(col_groups)),
              defaultweight(eltype(A), skipval=skipval, rev=rev))
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
