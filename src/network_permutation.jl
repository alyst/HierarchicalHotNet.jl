"""
    vertexbins(g::AbstractSimpleWeightedGraph,
               vertices::AbstractArray{Int}=1:nv(g);
               nbins::Integer=10, by::Symbol=:out, method=:tree) ->
        IndicesPartition

Partition the `vertices` of the weighted graph `g` into `nbins` bins,
so that the vertices of the same bin have similar sum of edge weights.

# Keyword arguments
* `nbins::Integer`: the number of vertex bins. Use *nbins* &ge; 1 for
   bigger (*nv* &ge; 1) networks
* `by::Symbol`: what edges to consider for grouping vertices.
  The supported options are:
  - `:out` (*default*, as in the original paper) --
    use only the outgoing edges
  - `:in` -- use only the incoming edges
  - `:outXin` (*recommended*) -- use all edges, so that the vertices of
    the same bin have similar sums of both incoming and outgoing edges
* `method::Symbol`: the method for binning vertices
  - `:sort` (as in the original paper) -- order the vertices by the
    sum of weights, then split the sorted array into equal bins.
    This method doesn't handle the `by=:outXin` well.
  - `:tree` (*default*) -- build the *hierarchical tree* of vertices
    based on their sum of weights, then cut the tree to get the
    requested number of bins.
"""
function vertexbins(g::AbstractSimpleWeightedGraph,
                    vertices::AbstractArray{Int}=1:nv(g);
                    by::Symbol=:out, method=:tree, kwargs...)
    if method == :sort
        return vertexbins_sort(g, vertices; by=by, kwargs...)
    elseif method == :tree
        return vertexbins_tree(g, vertices; by=by, kwargs...)
    else
        throw(ArgumentError("Unsupported vertex binning method=:$method"))
    end
end

function vertexbins_sort(g::AbstractSimpleWeightedGraph,
                         vertices::AbstractArray{Int};
                         by::Symbol=:out, nbins::Integer=10,
                         output::Symbol=:vertex)
    if (by == :in) || (by == :out)
        wvtxs = vec(sum(LightGraphs.weights(g), dims=by == :out ? 1 : 2))[vertices]
        vorder = sortperm(wvtxs, rev=true)
    elseif by == :outXin
        nbins_inner = 2^ceil(Int, log2(sqrt(nbins)))
        # bin vertices by in and out edges
        bin_out = fill(0, length(vertices))
        for (i, els) in enumerate(vertexbins_sort(g, vertices, by=:out, nbins=nbins_inner, output=:index))
            @inbounds bin_out[els] .= i
        end
        bin_in = fill(0, length(vertices))
        for (i, els) in enumerate(vertexbins_sort(g, vertices, by=:in, nbins=nbins_inner, output=:index))
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
    return IndicesPartition(output == :vertex ? vertices[vorder] : vorder, binstarts)
end

function vertexbins_tree(g::AbstractSimpleWeightedGraph,
                         vertices::AbstractArray{Int};
                         by::Symbol=:out, nbins::Integer=10,
                         weightf=sqrt, distf=sqrt)
    # vertex distances based on the in- and/or outcoming weights
    if (by == :in) || (by == :out)
        wvtxs = reshape(sum(weightf, LightGraphs.weights(g), dims=by == :out ? 1 : 2), (1, nv(g)))[:, vertices]
    elseif by == :outXin
        wvtxs = vcat(reshape(sum(weightf, LightGraphs.weights(g), dims=1), (1, nv(g))),
                     reshape(sum(weightf, LightGraphs.weights(g), dims=2), (1, nv(g))))[:, vertices]
    else
        throw(ArgumentError("Unsupported by=$by"))
    end
    wdists = pairwise(Euclidean(), wvtxs, dims=2)
    wtree = hclust(distf.(wdists), linkage=:ward, branchorder=:optimal)
    res = valuegroups(cutree(wtree, k=nbins))
    # convert from indices in vertices array to vertices
    @inbounds for (i, el) in enumerate(res.elems)
        res.elems[i] = vertices[el]
    end
    return res
end

"""
    randpermgroups!(v::AbstractVector{Int}, groups::AbstractPartition) -> v

Randomly reshuffle the elements of `groups` within each group and put
the result into `v`, group by group.
"""
function randpermgroups!(v::AbstractVector{Int}, groups::AbstractPartition)
    (length(v) >= nelems(groups)) ||
        throw(DimensionMismatch("Length of output vector ($(length(v))) is smaller than the number of grouped indices ($(nelems(groups)))"))
    group_perm = Vector{Int}() # temporary array of permuted indices within the group
    for group in groups
        randperm!(resize!(group_perm, length(group)))
        @inbounds for (i, j) in zip(group, group_perm)
            v[i] = group[j]
        end
    end
    return v
end

"""
    randpermgroups(groups::AbstractPartition) -> Vector{Int}

Randomly reshuffle the elements of `groups` within each group.
"""
randpermgroups(groups::AbstractPartition) =
    randpermgroups!(Vector{Int}(undef, nelems(groups)), groups)
