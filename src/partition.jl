"""
Base class for partitions of elements of type `T`.
"""
const AbstractPartition{T} = AbstractVector{<:AbstractVector{T}}

"""
    nelems(ptn::AbstractPartition)

Total number of elements in the partition.
"""
nelems(ptn::AbstractPartition) = mapreduce(length, +, ptn, init=0)

"""
The type of parts returned by `Partition{T}`
"""
const PartitionPart{T} = SubArray{T,1,Vector{T},Tuple{UnitRange{Int64}},true}

"""
Efficiently stores the partition of a vector of
elements of type `T` into disjoint sets.

Supports iterator interface.
"""
struct Partition{T} <: AbstractVector{PartitionPart{T}}
    elems::Vector{T}
    starts::Vector{Int} # starts of parts + the last element is length(element)+1

    Partition{T}() where T = new{T}(Vector{T}(), [1])
    # FIXME check starts?
    Partition{T}(elems::Vector{T},
                 starts::Vector{Int} = isempty(elems) ? [1] : [1, length(elems)+1]) where T =
        new{T}(elems, starts)
    Partition(elems::Vector{T},
              starts::Vector{Int} = isempty(elems) ? [1] : [1, length(elems)+1]) where T = Partition{T}(elems, starts)
end

"""
    nelems(ptn::AbstractPartition) -> Int

Number of parts in the partition.

Same as `length(ptn)`.
"""
nparts(ptn::Partition) = length(ptn.starts) - 1
nelems(ptn::Partition) = length(ptn.elems)

"""
    elems(ptn::Partition{T}) where T -> Vector{T}

Vector of elements in the partition.

Rerturns the vectors as they are stored in `ptn`:
the elements are ordered by their parts (elements from the first part come first etc)
and maintain the same order as within each part.
"""
elems(ptn::Partition) = ptn.elems

"""
    partrange(ptn::Partition, i::Integer) -> UnitRange{Int}

Range of indices associated with `i`-th part.
"""
Base.@propagate_inbounds partrange(ptn::Partition, i::Integer) = ptn.starts[i]:ptn.starts[i+1]-1
"""
    partlength(ptn::Partition, i::Integer) -> Int

Number of elements in `i`-th part.
"""
Base.@propagate_inbounds partlength(ptn::Partition, i::Integer) = (ptn.starts[i+1]-ptn.starts[i])
"""
    ispartempty(ptn::Partition, i::Integer) -> Bool

Whether `i`-th part is empty.
"""
Base.@propagate_inbounds ispartempty(ptn::Partition, i::Integer) = partlength(ptn, i) == 0

# fallback
Base.@propagate_inbounds partlength(ptn::AbstractVector{<:AbstractVector}, i::Integer) = length(ptn[i])
Base.@propagate_inbounds ispartempty(ptn::AbstractVector{<:AbstractVector}, i::Integer) = isempty(ptn[i])

"""
    pushelem!(ptn::Partition, el) -> ptn

Push the element into the last part of partition.
"""
pushelem!(ptn::Partition, el) = push!(ptn.elems, el)

"""
    closepart!(ptn::Partition) -> ptn

Finalize the last part of the partition and append
the new empty part.

All subsequent calls to [`pushelem!`](@ref) will add elements
into the new part.
"""
function closepart!(ptn::Partition)
    push!(ptn.starts, length(ptn.elems) + 1)
    return ptn
end

Base.IteratorEltype(::Type{<:Partition}) = Base.HasEltype()
Base.eltype(::Partition{T}) where T = PartitionPart{T}

Base.IteratorSize(::Type{<:Partition}) = Base.HasLength()
Base.length(ptn::Partition) = nparts(ptn)
Base.size(ptn::Partition) = (length(ptn),)

Base.@propagate_inbounds Base.iterate(ptn::Partition, i::Integer = 0) =
    i < length(ptn) ? (@inbounds(ptn[i+1]), i+1) : nothing

Base.@propagate_inbounds Base.getindex(ptn::Partition, i::Integer) =
    view(ptn.elems, partrange(ptn, i))

function Base.empty!(ptn::Partition)
    empty!(ptn.elems)
    empty!(ptn.starts)
    push!(ptn.starts, 1)
    return ptn
end

function Base.copy!(a::Partition{T}, b::Partition{T}) where T
    resize!(a.elems, length(b.elems))
    copyto!(a.elems, b.elems)
    resize!(a.starts, length(b.starts))
    copyto!(a.starts, b.starts)
    return a
end

function Base.filter!(f, dest::Partition{T}, src::Partition{T}) where T
    empty!(dest)
    for part in src
        f(part) || continue
        append!(dest.elems, part)
        closepart!(dest)
    end
    return dest
end

Base.filter(f, src::Partition{T}) where T = filter!(f, Partition{T}(), src)

"""
    repeat!(ptn::Partition, n::Integer) -> ptn

Repeat the parts of `ptn` `n` times.

```jldoctest
julia> using HierarchicalHotNet

julia> ptn = HierarchicalHotNet.IndicesPartition(5, nparts=1)
1-element HierarchicalHotNet.Partition{Int64}:
 [1, 2, 3, 4, 5]

julia> HierarchicalHotNet.repeat!(ptn, 3)
3-element HierarchicalHotNet.Partition{Int64}:
 [1, 2, 3, 4, 5]
 [1, 2, 3, 4, 5]
 [1, 2, 3, 4, 5]
```
"""
function repeat!(ptn::Partition, n::Integer)
    (n == 1) && return ptn
    (n == 0) && return empty!(ptn)
    (n < 0) && throw(ArgumentError("n must be non-negative"))
    if !isempty(ptn.elems)
        l = length(ptn.elems)
        resize!(ptn.elems, l*n)
        @inbounds for i in 1:(n-1)
            copyto!(ptn.elems, 1 + l*i, ptn.elems, 1, l)
        end
        np = length(ptn.starts)
        laststart = last(ptn.starts)
        @inbounds for _ in 1:(n-1), j in 2:np
            laststart += ptn.starts[j] - ptn.starts[j-1]
            push!(ptn.starts, laststart)
        end
    else
        fill!(resize!(ptn.starts, (length(ptn.starts)-1)*n + 1), 1)
    end
    return ptn
end

"""
Partition of an integer vector
"""
const IndicesPartition = Partition{Int}
const IndicesPartitionPart = eltype(IndicesPartition)

IndicesPartition(n::Integer=0; nparts::Integer=n) =
    reset!(IndicesPartition(Vector{Int}(undef, n), Vector{Int}()), nparts=nparts)

"""
    reset!(p::IndicesPartition, n::Integer=nelems(p); nparts::Integer=n) -> p

Resets the integer partition `p`.

If `nparts == 1`, the partition is reset into `[[1, 2, 3, ..., n]]`.
If `n == nparts`, sets the partition to `[[1], [2], [3], ..., [n]]`.
"""
function reset!(p::IndicesPartition, n::Integer=nelems(p); nparts::Integer=n)
    if (nparts != n) && (nparts != 1)
        throw(ArgumentError("No default method to divide $n element(s) into $nparts group(s)"))
    end
    resize!(p.elems, n) .= 1:n
    resize!(p.starts, nparts+1)
    if nparts == 1
        p.starts[1] = 1
        p.starts[2] = n+1
    elseif nparts == n
        p.starts .= 1:n+1
    end
    return p
end
