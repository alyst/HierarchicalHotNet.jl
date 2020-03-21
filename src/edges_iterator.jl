abstract type AbstractOutedgesIterator{W <: Number} end

Base.IteratorEltype(::Type{<:AbstractOutedgesIterator}) = Base.HasEltype()
Base.eltype(::Type{AbstractOutedgesIterator{W}}) where W = Pair{Int, W}
Base.eltype(it::AbstractOutedgesIterator) = eltype(typeof(it))

Base.IteratorSize(::Type{<:AbstractOutedgesIterator}) = Base.SizeUnknown()

isweightreverse(it::AbstractOutedgesIterator) = isweightreverse(typeof(it))
isvalidedge(w::Number, it::AbstractOutedgesIterator) =
    isvalidedge(w, skipval=skipval(it), threshold=it.threshold, rev=isweightreverse(it))

struct MatrixOutedgesIterator{W <: Number, V <: AbstractVector, S, T, R} <: AbstractOutedgesIterator{W}
    col::V
    skipval::S
    threshold::T

    @inline function MatrixOutedgesIterator(mtx::AbstractMatrix, v::Integer;
        skipval::Union{Number, Nothing}=zero(eltype(mtx)),
        threshold::Union{Number, Nothing}=nothing,
        rev::Bool=false
    )
        W = eltype(mtx)
        S = isnothing(skipval) ? Nothing : W
        T = isnothing(threshold) ? Nothing : W
        col = view(mtx, :, v)
        new{W, typeof(col), S, T, rev}(col, skipval, threshold)
    end
end

outedges(mtx::AbstractMatrix, v::Integer; kwargs...) =
    MatrixOutedgesIterator(mtx, v; kwargs...)

skipval(it::MatrixOutedgesIterator) = it.skipval
isweightreverse(::Type{<:MatrixOutedgesIterator{<:Any, <:Any, <:Any, <:Any, R}}) where R = R
isweightreverse(it::MatrixOutedgesIterator) = isweightreverse(typeof(it))
isvalidedge(w::Number, it::MatrixOutedgesIterator) =
    isvalidedge(w, skipval=it.skipval, threshold=it.threshold, rev=isweightreverse(it))

@inline function Base.iterate(it::MatrixOutedgesIterator, i::Integer = 0)
    while i < length(it.col)
        i += 1
        w = @inbounds(it.col[i])
        isvalidedge(w, it) && return (i => w, i)
    end
    return nothing
end
