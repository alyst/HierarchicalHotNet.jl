"""
Distribute the vertices of the weighted graph `g` into `nbins` bins,
so that the vertices that have similar sum of outgoing/incoming edges are
put into the same bin.
"""
function vertexbins(g::AbstractSimpleWeightedGraph;
                    by::Symbol=:out, nbins::Integer=10)
    if (by == :in) || (by == :out)
        wvtxs = vec(sum(weights(g), dims=by == :out ? 1 : 2))
        vorder = sortperm(wvtxs, rev=true)
    else
        throw(ArgumentError("Unsupported by=$by"))
    end
    binstarts = Vector{Int}(undef, nbins+1)
    maxbinsize = fld1(length(vorder), nbins)
    curpos = 1
    for i in 1:nbins
        binstarts[i] = curpos
        curpos += maxbinsize
        if curpos > length(vorder)
            curpos = length(vorder) + 1
        end
    end
    binstarts[end] = length(vorder) + 1
    return IndicesPartition(vorder, binstarts)
end

"""
Generates a random permutation of elements within each specified group.
"""
function randpermgroups!(v::AbstractVector{Int}, groups::AbstractPartition)
    length(v) == nelems(groups) ||
        throw(DimensionMismatch("Length of output vector ($(length(v))) does not match the number of grouped indices ($(nindices(groups)))"))
    group_perm = Vector{Int}()
    for group in groups
        randperm!(resize!(group_perm, length(group)))
        @inbounds for (i, j) in zip(group, group_perm)
            v[i] = group[j]
        end
    end
    return v
end

randpermgroups(groups::AbstractPartition) =
    randpermgroups!(Vector{Int}(undef, nelems(groups)), groups)
