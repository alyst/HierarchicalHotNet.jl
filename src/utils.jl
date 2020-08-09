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

# position on the Hilbert curve
"""
    hilbertorder(i::Integer, j::Integer, n::Integer)

Get the index of the point `(i, j)` on the *Hilbert curve* that passes through
the points of `n×n` grid, where `n` is some power of 2.
"""
function hilbertorder(i::Integer, j::Integer, n::Integer)
    (n > 0) || throw(ArgumentError("n must be positive, n=$n given"))
    (i > 0) || throw(ArgumentError("i must be positive, i=$i given"))
    (j > 0) || throw(ArgumentError("j must be positive, j=$j given"))
    d = 1
    s = n
    i -= 1
    j -= 1
    n -= 1
    while s > 1
        s >>= 1
        ri = (i & s) > 0
        rj = (j & s) > 0
        if ri && rj
            d += 2s*s
        elseif ri && !rj
            d += s*s
        elseif !ri && rj
            d += 3s*s
            i, j = n - j, n - i
        else
            i, j = j, i
        end
    end
    return d
end

## Empirical CDF

struct EmCDF{T <: Real, W <: Real, I}
    sorted_values::Vector{Tuple{T, W, W, W}}
end

weighttype(::Type{<:EmCDF{<:Any, W}}) where W = W
weighttype(ecdf::EmCDF) = weighttype(typeof(ecdf))
isinterpolating(::Type{<:EmCDF{<:Any, <:Any, I}}) where I = I
isinterpolating(ecdf::EmCDF) = isinterpolating(typeof(ecdf))

function (ecdf::EmCDF)(x::Real)
    isnan(x) && return NaN
    pos = searchsortedlast(ecdf.sorted_values, x, by=first)
    (pos == 0) && return zero(weighttype(ecdf))
    @inbounds val, cdf_val, val_invdelta, cdf_delta = ecdf.sorted_values[pos]
    if !isinterpolating(ecdf) || (val == x) || (pos == length(ecdf.sorted_values))
        return cdf_val
    else
        return muladd(cdf_delta, min((x - val) * val_invdelta, one(weighttype(ecdf))), cdf_val)
    end
end

# broadcasts ecdf() over an array
# caching the previous calculated value
function Base.Broadcast.broadcasted(ecdf::EmCDF, v::AbstractArray)
    res = similar(v, weighttype(ecdf))
    @inbounds for i in eachindex(v)
        res[i] = i == 1 || v[i] != v[i-1] ? ecdf(v[i]) : res[i-1]
    end
    return res
end

function ecdf(X::AbstractVector,
              weights::Union{Nothing, AbstractVector}=nothing,
              interpolate::Bool=false)
    any(isnan, X) && throw(ArgumentError("ecdf can not include NaN values"))
    evenweights = isnothing(weights) || isempty(weights)
    evenweights || (length(X) == length(weights)) ||
        throw(ArgumentError("data and weight vectors must be the same size, " *
                            "got $(length(X)) and $(length(weights))"))
    T = eltype(X)
    W0 = evenweights ? Int : eltype(weights)
    W = isnothing(weights) ? Float64 : eltype(one(W0)/sum(weights))

    wsum = evenweights ? length(X) : sum(weights)
    ord = sortperm(X)

    sorted_vals = sizehint!(Vector{Tuple{T, W, W, W}}(), length(X))
    isempty(X) && return EmCDF{T, W, interpolate}(sorted_vals)

    valprev = val = X[first(ord)]
    wsumprev = zero(W0)
    valw = zero(W0)

    push_valprev!() = push!(sorted_vals, (valprev, min(wsumprev/wsum, one(W)),
                                          inv(val - valprev), valw/wsum))

    @inbounds for i in ord
        valnew = X[i]
        if (val != valnew) || (i == last(ord))
            (wsumprev > 0) && push_valprev!()
            valprev = val
            val = valnew
            wsumprev += valw
            valw = zero(W0)
        end
        valw += evenweights ? one(W0) : weights[i]
    end
    #@assert valw + wsumprev == wsum # may fail due to fp-arithmetic
    (wsumprev > 0) && push_valprev!()
    # last value
    push!(sorted_vals, (val, one(W), zero(W), zero(W)))
    return EmCDF{T,W,interpolate}(sorted_vals)
end

Base.minimum(ecdf::EmCDF) = first(ecdf.sorted_values)

Base.maximum(ecdf::EmCDF) = last(ecdf.sorted_values)

Base.extrema(ecdf::EmCDF) = (minimum(ecdf), maximum(ecdf))
