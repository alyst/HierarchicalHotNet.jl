"""
Node of the [`SCCTree`](@ref).
"""
struct SCCTreeNode{T}
    parent::Int     # index of the parent node, 0 if no parent (= root node)
    threshold::T    # min/max threshold at which the node is a connected component
    nvertices::Int  # total number of vertices in a component represented by a node

    children::Vector{Int}   # children node indices
    vertices::Vector{Int}   # indices of graph vertices directly attached to the node
end

nvertices(node::SCCTreeNode) = node.nvertices

Base.:(==)(a::SCCTreeNode{T}, b::SCCTreeNode{T}) where T =
    (a.parent == b.parent) &&
    isequal(a.threshold, b.threshold) && # there could be NaNs
    (a.nvertices == b.nvertices) &&
    (a.children == b.children) &&
    (a.vertices == b.vertices)

"""
The tree of *strongly connected components*.
Organizes the strongly connected components of the weighted directed graph into
a tree.

Cutting the tree at specific threshold with [`cut`](@ref) gives the corresponding strongly
connected components. The root of the tree corresponds to the weakest threshold.
"""
struct SCCTree{T}
    rev::Bool                           # if greater threshold is weaker than lower
    nodes::Vector{SCCTreeNode{T}}       # tree nodes, sorted root to bottom
    vertexnodes::Vector{Int}            # vertex index -> index of the node it is attached to
    thresholds::Vector{T}               # all thresholds producing different SCCs, from weak to strong

    function SCCTree(tree::SCCSeedling{T}) where T
        # calculate parents and node sizes
        parents = fill(0, length(tree.nodes))
        @inbounds for (nodeix, node) in enumerate(tree.nodes)
            @inbounds parents[node.children] .= nodeix
        end
        root = findfirst(iszero, parents)
        @assert isempty(tree.nodes) && isnothing(root) || root > 0
        @assert isnothing(root) || (findnext(iszero, parents, root+1) === nothing)

        # sort nodes from root to vertices by
        # a) threshold
        # b) the index of the first vertex
        nodesorder = sortperm(tree.nodes, lt=(a, b) ->
            # weaker threshold comes first (including root with undefined=0 threshold)
            (a.threshold < b.threshold) ||
            # component with lower minimal-index vertex comes first
            ((a.threshold == b.threshold) && (a.firstvertex <= b.firstvertex)))
        # map between old and new indices
        old2new = invperm(nodesorder)
        @assert isnothing(root) || (old2new[root] == 1) # root is always the first

        # create SCCTree nodes
        thresholds = Vector{T}()
        newnodes = sizehint!(Vector{SCCTreeNode{T}}(), length(tree.nodes))
        for (newix, oldix) in enumerate(nodesorder)
            @assert newix == length(newnodes) + 1
            node = tree.nodes[oldix]
            if (node.threshold != 0) &&
               (isempty(thresholds) || last(thresholds) != tree.weights[node.threshold])
                push!(thresholds, tree.weights[node.threshold])
            end
            # update children node indices
            newchildren = node.children
            @inbounds for (i, oldchild) in enumerate(newchildren)
                newchildren[i] = old2new[oldchild]
            end
            sort!(newchildren)
            push!(newnodes, SCCTreeNode{T}(oldix != root ? @inbounds(old2new[parents[oldix]]) : 0,
                                           node.threshold == 0 ? NaN : tree.weights[node.threshold],
                                           node.nvertices,
                                           newchildren, node.vertices))
        end
        @assert isempty(newnodes) || (newnodes[1].parent == 0) # root
        # update indices of nodes that bind directly to graph vertices
        vertexnodes = tree.vertexnodes
        @inbounds for (i, oldix) in enumerate(vertexnodes)
            vertexnodes[i] = old2new[oldix]
        end
        return new{T}(tree.rev, newnodes, vertexnodes, thresholds)
    end
end

weighttype(::Type{SCCTree{T}}) where T = T
weighttype(tree::SCCTree) = weighttype(typeof(tree))

nvertices(tree::SCCTree) = length(tree.vertexnodes)
vertices(tree::SCCTree, parent::Integer) = appendvertices!(Vector{Int}(), tree, parent)

Base.:(==)(a::SCCTree{T}, b::SCCTree{T}) where T =
    (a.rev == b.rev) &&
    (a.nodes == b.nodes) &&
    (a.vertexnodes == b.vertexnodes) &&
    (a.thresholds == b.thresholds) # FIXME there could be NaNs?

# Appends vertices from the subtree with the root at `nodeindex`
# to the `vertices` vector/indices partition
function appendvertices!(vertices::AbstractVector, tree::SCCTree,
                         nodeindex::Integer)
    (nodeindex >= 1) || throw(ArgumentError("Incorrect node index ($nodeindex)"))
    @inbounds node = tree.nodes[nodeindex]
    @inbounds for vertex in node.vertices
        if vertices isa Partition # FIXME ugly way to append to a partition
            pushelem!(vertices, vertex)
        else
            push!(vertices, vertex)
        end
    end
    for child in node.children
        appendvertices!(vertices, tree, child)
    end
    nothing
end

"""
    cut(tree, threshold; minsize=1) -> IndicesPartition

Cuts the tree at the given `threshold` to get the corresponding
strongly connected components of the original graph.
`minsize` optional parameter specifies whether components smaller than that
would be ignored.
"""
cut(tree::SCCTree, threshold::Number; kwargs...) =
    cut!(IndicesPartition(0), tree, threshold; kwargs...)

function cut!(comps::IndicesPartition, tree::SCCTree, threshold::Number;
              minsize::Integer=1)
    reset!(comps, 0)
    if !isempty(tree.nodes)
        cut_subtree!(comps, tree, threshold, minsize, 1,
                     first(tree.nodes).threshold)
    end
    return comps
end

# helper method for cut!()
# cuts the subtree specified by the nodeindex at given threshold
# the cutting (generating new components) depends on parent_threshold (threshold of the parents)
function cut_subtree!(comps::IndicesPartition, tree::SCCTree,
                      threshold::Number, minsize::Integer,
                      nodeindex::Int, parent_threshold::Number)
    node = tree.nodes[nodeindex]
    (nvertices(node) < minsize) && return
    notweaker = !isnan(node.threshold) && !isweaker(node.threshold, threshold, rev=tree.rev)
    if notweaker
        #@info "node #$(nodeindex) threshold=$(node.threshold), marking as component #$nextcompindex"
        appendvertices!(comps, tree, nodeindex)
    else
        for child in node.children
            cut_subtree!(comps, tree, threshold, minsize, child, node.threshold)
        end
        # vertices are in a component of its own
        if minsize <= 1
            @inbounds for vertex in node.vertices
                pushelem!(comps, vertex)
                closepart!(comps, sort=true)
            end
        end
        return
    end
    if (nodeindex == 1) || isnan(parent_threshold) ||
       isweaker(parent_threshold, threshold, rev=tree.rev)
        #@info "node #$(nodeindex) (size=$(node.nvertices), threshold=$(node.threshold)): parent threshold is $(parent_threshold), so advancing component $(nextcompindex+1)"
        closepart!(comps, sort=true)
    end
    nothing
end
