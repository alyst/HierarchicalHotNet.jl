"""
Distribute the vertices of the weighted graph `g` into `nbins` bins,
so that the vertices that have similar sum of outgoing/incoming edges are
put into the same bin.
"""
function vertexbins(g::AbstractSimpleWeightedGraph;
                    by::Symbol=:out, method=:tree, kwargs...)
    if method == :sort
        return vertexbins_sort(g; by=by, kwargs...)
    elseif method == :tree
        return vertexbins_tree(g; by=by, kwargs...)
    else
        throw(ArgumentError("Unsupported vertex binning method=:$method"))
    end
end

function vertexbins_sort(g::AbstractSimpleWeightedGraph;
                         by::Symbol=:out, nbins::Integer=10)
    if (by == :in) || (by == :out)
        wvtxs = vec(sum(LightGraphs.weights(g), dims=by == :out ? 1 : 2))
        vorder = sortperm(wvtxs, rev=true)
    elseif by == :outXin
        nbins_inner = 2^ceil(Int, log2(sqrt(nbins)))
        # bin vertices by in and out edges
        bin_out = fill(0, nv(g))
        for (i, els) in enumerate(vertexbins(g, by=:out, nbins=nbins_inner))
            @inbounds bin_out[els] .= i
        end
        bin_in = fill(0, nv(g))
        for (i, els) in enumerate(vertexbins(g, by=:in, nbins=nbins_inner))
            @inbounds bin_in[els] .= i
        end
        # calculate vertex positions on the hibert curve going through the in and out bin ids
        vtx_hilberpos = hilbertorder.(bin_out, bin_in, nbins_inner)
        # sort vertices by hilber position
        vorder = sortperm(vtx_hilberpos)
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

function vertexbins_tree(g::AbstractSimpleWeightedGraph;
                         by::Symbol=:out, nbins::Integer=10,
                         weightf=sqrt, distf=sqrt)
    # vertex distances based on the in- and/or outcoming weights
    if (by == :in) || (by == :out)
        wvtxs = reshape(sum(weightf, LightGraphs.weights(g), dims=by == :out ? 1 : 2), (1, nv(g)))
    elseif by == :outXin
        wvtxs = vcat(reshape(sum(weightf, LightGraphs.weights(g), dims=1), (1, nv(g))),
                     reshape(sum(weightf, LightGraphs.weights(g), dims=2), (1, nv(g))))
    else
        throw(ArgumentError("Unsupported by=$by"))
    end
    wdists = pairwise(Euclidean(), wvtxs, dims=2)
    wtree = hclust(distf.(wdists), linkage=:ward, branchorder=:optimal)
    return valuegroups(cutree(wtree, k=nbins))
end

"""
Generates a random permutation of elements within each specified group.
"""
function randpermgroups!(v::AbstractVector{Int}, groups::AbstractPartition)
    length(v) == nelems(groups) ||
        throw(DimensionMismatch("Length of output vector ($(length(v))) does not match the number of grouped indices ($(nelems(groups)))"))
    group_perm = Vector{Int}() # temporary array of permuted indices within the group
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
