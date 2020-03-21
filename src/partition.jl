const AbstractPartition{T} = AbstractVector{<:AbstractVector{T}}
nelems(ptn::AbstractPartition) = mapreduce(length, +, ptn, init=0)

# type of parts returned by `Partition{T}`
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

nparts(ptn::Partition) = length(ptn.starts) - 1
nelems(ptn::Partition) = length(ptn.elems)

elems(ptn::Partition) = ptn.elems

# range of indices associated with i-th part
partrange(ptn::Partition, i::Integer) = ptn.starts[i]:ptn.starts[i+1]-1
# push element to the last part of partition
pushelem!(ptn::Partition, el) = push!(ptn.elems, el)
# finalize the current part
# all subsequent pushelem!() will add elements to the next part
function closepart!(ptn::Partition)
    push!(ptn.starts, length(ptn.elems) + 1)
    return ptn
end

Base.IteratorEltype(::Type{<:Partition}) = Base.HasEltype()
Base.eltype(::Partition{T}) where T = PartitionPart{T}

Base.IteratorSize(::Type{<:Partition}) = Base.HasLength()
Base.length(ptn::Partition) = nparts(ptn)
Base.size(ptn::Partition) = (length(ptn),)

Base.iterate(ptn::Partition, i::Integer = 0) =
    i < length(ptn) ? (ptn[i+1], i+1) : nothing

Base.getindex(ptn::Partition, i::Integer) =
    view(ptn.elems, partrange(ptn, i))

function Base.empty!(ptn::Partition)
    empty!(ptn.elems)
    empty!(ptn.starts)
    push!(ptn.starts, 1)
    return ptn
end

function Base.copy!(a::Partition{T}, b::Partition{T}) where T
    copy!(a.elems, b.elems)
    copy!(a.starts, b.starts)
    return a
end

# partition of an integer vector
const IndicesPartition = Partition{Int}

IndicesPartition(n::Integer=0; ngroups::Integer=n) =
    reset!(IndicesPartition(Vector{Int}(undef, n), Vector{Int}()), ngroups=ngroups)

"""
Resets `IntegerPartition`.

If `ngroups == 1`, the partition is reset into `[[1, 2, 3, ..., n]]`.
If `n == ngroups`, sets the partition to `[[1], [2], [3], ..., [n]]`.
"""
function reset!(p::IndicesPartition, n::Integer=length(p.elems); ngroups::Integer=n)
    if (ngroups != n) && (ngroups != 1)
        throw(ArgumentError("No default method to divide $n element(s) into $ngroups group(s)"))
    end
    resize!(p.elems, n) .= 1:n
    resize!(p.starts, ngroups+1)
    if ngroups == 1
        p.starts[1] = 1
        p.starts[2] = n+1
    elseif ngroups == n
        p.starts .= 1:n+1
    end
    return p
end
