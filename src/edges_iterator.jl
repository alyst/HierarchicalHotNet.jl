abstract type AbstractOutedgesIterator{W <: Number} end

Base.IteratorEltype(::Type{<:AbstractOutedgesIterator}) = Base.HasEltype()
Base.eltype(::Type{AbstractOutedgesIterator{W}}) where W = Pair{Int, W}
Base.eltype(it::AbstractOutedgesIterator) = eltype(typeof(it))

Base.IteratorSize(::Type{<:AbstractOutedgesIterator}) = Base.SizeUnknown()

isweightreverse(it::AbstractOutedgesIterator) = isweightreverse(typeof(it))
isvalidedge(w::Number, it::AbstractOutedgesIterator) =
    isvalidedge(w, skipval=skipval(it), threshold=it.threshold, rev=isweightreverse(it))

struct MatrixOutedgesIterator{W <: Number, V <: AbstractVector, E <: EdgeTest} <: AbstractOutedgesIterator{W}
    col::V
    test::E

    @inline function MatrixOutedgesIterator(mtx::AbstractMatrix{W},
                                            v::Integer, test::EdgeTest{W}) where W
        col = view(mtx, :, v)
        new{W, typeof(col), typeof(test)}(col, test)
    end
end

@inline outedges(mtx::AbstractMatrix{W}, v::Integer,
                 test::EdgeTest = EdgeTest{W}()) where W =
    MatrixOutedgesIterator(mtx, v, test)

@inline function Base.iterate(it::MatrixOutedgesIterator, i::Integer = 0)
    l = length(it.col)
    while i < l
        i += 1
        w = @inbounds(it.col[i])
        isvalidedge(w, it.test) && return (i => w, i)
    end
    return nothing
end
