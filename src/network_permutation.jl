"""
Distribute the vertices of the weighted graph `g` into `nbins` bins,
so that the vertices that have similar sum of outgoing edges are put into
the same bin.
"""
function vertexbins(g::AbstractSimpleWeightedGraph;
                    nbins::Integer=10)
    wvtxs = vec(sum(weights(g), dims=1))
    vdescorder = sortperm(wvtxs, rev=true)
    binstarts = Vector{Int}(undef, nbins+1)
    maxbinsize = fld1(length(wvtxs), nbins)
    curpos = 1
    for i in 1:nbins
        binstarts[i] = curpos
        curpos += maxbinsize
        if curpos > length(vdescorder)
            curpos = length(vdescorder) + 1
        end
    end
    binstarts[end] = length(vdescorder) + 1
    return IndicesPartition(vdescorder, binstarts)
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
