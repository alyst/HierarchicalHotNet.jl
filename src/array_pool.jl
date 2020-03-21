"""
Helps to maintain the pool of reusable arrays of different sizes
and reduce the burden on garbage collection.
"""
mutable struct ArrayPool{T}
    arrays::Vector{Vector{T}}
    nborrowed::Int
    borrow_limit::Int

    ArrayPool{T}(borrow_limit::Integer = 0) where T =
        new{T}(Vector{Vector{T}}(), 0, borrow_limit)
end

Base.eltype(::Type{ArrayPool{T}}) where T = T
Base.eltype(pool::ArrayPool) = eltype(typeof(pool))

"""
Gets an array of specific size from the pool.
The returned array should be returned back to the pool using `release!()`.
"""
function borrow!(pool::ArrayPool{T}, len::Integer = 0) where T
    (pool.borrow_limit > 0) && (pool.nborrowed >= pool.borrow_limit) &&
        error("Pool has reached the borrow limit ($(pool.nborrowed) array(s)), check your code")
    pool.nborrowed += 1
    return isempty(pool.arrays) ?
            Vector{T}(undef, len) :
            resize!(pop!(pool.arrays), len)
end

borrow!(pool::ArrayPool{T}, size::NTuple{N, Int}) where {T, N} =
    reshape(borrow!(pool, prod(size)), size)

borrow!(::Type{T}, pool::ArrayPool{T}, size=0) where T =
    borrow!(pool, size)

"""
Releases an array returned by `acquire!()` back into the pool.
"""
function release!(pool::ArrayPool{T}, arr::Array{T}) where T
    (pool.nborrowed > 0) ||
        error("Returning an array that was never borrowed from the pool")
    push!(pool.arrays, vec(arr))
    pool.nborrowed -= 1
    return pool
end

# no-pool versions of borrow and release to make
# the ArrayPool usage optional and transparent
borrow!(::Type{T}, ::Nothing, len::Integer = 0) where T =
    Vector{T}(undef, len)
borrow!(::Type{T}, ::Nothing, size::NTuple{T, N}) where {T, N} =
    Array{T, B}(undef, size)

# do nothing
release!(::Nothing, arr::Array) = ()
