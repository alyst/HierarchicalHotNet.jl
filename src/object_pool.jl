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
    release!(pool::ObjectPool{T}, obj::T) where T

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
    borrow!(pool::ObjectPool{T}) where T -> T

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

"""
Alias for [`ObjectPool`](@ref HierarchicalHotNet.ObjectPool) for arrays.
"""
const ArrayPool{T} = ObjectPool{Vector{T}}

Base.eltype(::Type{<:AbstractArrayPool{T}}) where T = T
Base.eltype(pool::AbstractArrayPool) = eltype(typeof(pool))

"""
    borrow!(pool::ArrayPool{T}, len::Integer = 0) where T -> T

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
    release!(pool::ArrayPool{T}, arr::Array{T}) where T

Releases an array returned by `borrow!()` back into the pool.
"""
function release!(pool::ArrayPool{T}, arr::Array{T}) where T
    (pool.nborrowed > 0) ||
        error("Returning an array that was never borrowed from the pool")
    push!(pool.objects, vec(arr))
    pool.nborrowed -= 1
    return pool
end

"""
Fake array pool that just constructs arrays.
"""
struct NoopObjectPool{T} <: AbstractObjectPool{T}
    NoopObjectPool{T}(borrow_limit::Integer = 0) where T =
        new{T}()
end

borrow!(::NoopObjectPool{T}) where T = T()

# do nothing
release!(::NoopObjectPool{T}, obj::T) where T = ()

const NoopArrayPool{T} = NoopObjectPool{Vector{T}}

borrow!(::NoopArrayPool{T}, len::Integer) where T =
    Vector{T}(undef, len)
borrow!(::NoopArrayPool{T}, size::NTuple{N}) where {T, N} =
    Array{T, N}(undef, size)

release!(::NoopArrayPool{T}, arr::Array{T}) where T = ()

# no-pool versions of borrow and release to make
# the ArrayPool usage optional and transparent
borrow!(::Type{T}, ::Nothing) where T = T()
borrow!(::Type{Vector{T}}, ::Nothing, len::Integer) where T =
    Vector{T}(undef, len)
borrow!(::Type{Vector{T}}, ::Nothing, size::NTuple{N}) where {T, N} =
    Array{T, N}(undef, size)

# do nothing
release!(::Nothing, obj::Any) = ()

"""
Pool of object pools.
Manages `ObjectPool` objects for different types.
"""
struct ObjectPools
    typepools::Dict{Type, ObjectPool}
    default_borrow_limit::Int;

    ObjectPools(default_borrow_limit::Integer = 0) =
        new(Dict{Type, ObjectPool}(), default_borrow_limit)
end

objpool(pools::ObjectPools, ::Type{T};
        borrow_limit::Union{Integer, Nothing} = nothing) where T =
    get!(() -> ObjectPool{T}(isnothing(borrow_limit) ?
                             pools.default_borrow_limit : borrow_limit),
        pools.typepools, T)::ObjectPool{T}

objpool(::Nothing, ::Type{T};
        borrow_limit::Union{Integer, Nothing} = nothing) where T =
    NoopObjectPool{T}()

arraypool(pools::Union{ObjectPools, Nothing}, type::Type{T};
          kwargs...) where T =
    objpool(pools, Vector{T}, kwargs...)

Base.getindex(pools::ObjectPools, type::Type; kwargs...) =
    objpool(pools, type; kwargs...)

Base.haskey(pools::ObjectPools, type::Type) = haskey(pools.typepools, type)
