using SparseArrays

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

struct SparseMatrixOutedgesIterator{W <: Number, R <: AbstractVector, E <: EdgeTest} <: AbstractOutedgesIterator{W}
    mtx::SparseMatrixCSC{W}
    outrange::R
    test::E

    @inline function SparseMatrixOutedgesIterator(mtx::SparseMatrixCSC{W},
                                                  v::Integer, test::EdgeTest{W}) where W
        !isnothing(test.threshold) && (test.threshold <= 0) &&
            throw(ArgumentError("EdgeTest with nonpositive thresholds are not supported for SparseMatrices"))
        outrange = nzrange(mtx, v)
        new{W, typeof(outrange), typeof(test)}(mtx, outrange, test)
    end
end

@inline outedges(mtx::SparseMatrixCSC{W}, v::Integer,
                 test::EdgeTest = EdgeTest{W}()) where W =
    SparseMatrixOutedgesIterator(mtx, v, test)

@inline function Base.iterate(it::SparseMatrixOutedgesIterator, i::Integer = 0)
    if i != 0
        @assert i >= first(it.outrange)-1
    else
        i = first(it.outrange)-1
    end
    l = last(it.outrange)
    while i < l
        i += 1
        w = @inbounds(nonzeros(it.mtx)[i])
        isvalidedge(w, it.test) && return (@inbounds(rowvals(it.mtx)[i]) => w, i)
    end
    return nothing
end
