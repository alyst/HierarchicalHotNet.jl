abstract type AbstractObjectPool{T} end

"""
Helps to maintain the pool of reusable objects
and reduce the burden on garbage collection.
"""
mutable struct ObjectPool{T} <: AbstractObjectPool{T}
    objects::Vector{T}
    nborrowed::Int
    borrow_limit::Int

    ObjectPool{T}(borrow_limit::Integer = 0) where T =
        new{T}(Vector{T}(), 0, borrow_limit)
end

"""
Releases an object acquired by `borrow!()` back into the pool.
"""
function release!(pool::ObjectPool{T}, obj::T) where T
    (pool.nborrowed > 0) ||
        error("Returning an object that was never borrowed from the pool")
    push!(pool.objects, obj)
    pool.nborrowed -= 1
    return pool
end

"""
Gets an object from the pool.
The returned object should be returned back to the pool using `release!()`.
"""
function borrow!(pool::ObjectPool{T}) where T
    (pool.borrow_limit > 0) && (pool.nborrowed >= pool.borrow_limit) &&
        error("Pool has reached the borrow limit ($(pool.nborrowed) object(s)), check your code")
    pool.nborrowed += 1
    return isempty(pool.objects) ? T() : pop!(pool.objects)
end

objtype(::Type{<:AbstractObjectPool{T}}) where T = T
objtype(pool::AbstractObjectPool) = objtype(typeof(pool))

const AbstractArrayPool{T} = AbstractObjectPool{Vector{T}}

const ArrayPool{T} = ObjectPool{Vector{T}}

Base.eltype(::Type{<:AbstractArrayPool{T}}) where T = T
Base.eltype(pool::AbstractArrayPool) = eltype(typeof(pool))

"""
Gets an array of specific size from the pool.
The returned array should be returned back to the pool using `release!()`.
"""
function borrow!(pool::ArrayPool{T}, len::Integer = 0) where T
    (pool.borrow_limit > 0) && (pool.nborrowed >= pool.borrow_limit) &&
        error("Pool has reached the borrow limit ($(pool.nborrowed) array(s)), check your code")
    pool.nborrowed += 1
    return isempty(pool.objects) ?
            Vector{T}(undef, len) :
            resize!(pop!(pool.objects), len)
end

borrow!(pool::AbstractArrayPool{T}, size::NTuple{N}) where {T, N} =
    reshape(borrow!(pool, prod(size)), size)

borrow!(::Type{T}, pool::AbstractArrayPool{T}, size=0) where T =
    borrow!(pool, size)

"""
Releases an array returned by `borrow!()` back into the pool.
"""
function release!(pool::ArrayPool{T}, arr::Array{T}) where T
    (pool.nborrowed > 0) ||
        error("Returning an array that was never borrowed from the pool")
    push!(pool.objects, vec(arr))
    pool.nborrowed -= 1
    return pool
end

# no-pool versions of borrow and release to make
# the ArrayPool usage optional and transparent
borrow!(::Type{T}, ::Nothing) where T = T()
borrow!(::Type{Vector{T}}, ::Nothing, len::Integer) where T =
    Vector{T}(undef, len)
borrow!(::Type{Vector{T}}, ::Nothing, size::NTuple{N}) where {T, N} =
    Array{T, N}(undef, size)

# do nothing
release!(::Nothing, obj::Any) = ()
