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
        empty!(tree.nodes) # avoid reusing the nodes as their arrays belong to SCCTree now
        @assert isempty(newnodes) || (newnodes[1].parent == 0) # root
        # update indices of nodes that bind directly to graph vertices
        vertexnodes = similar(tree.vertexnodes)
        @inbounds for (i, oldix) in enumerate(tree.vertexnodes)
            vertexnodes[i] = old2new[oldix]
        end
        return new{T}(tree.rev, newnodes, vertexnodes, thresholds)
    end
end

weighttype(::Type{SCCTree{T}}) where T = T
weighttype(tree::SCCTree) = weighttype(typeof(tree))

nthresholds(tree::SCCTree) = length(tree.thresholds)
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
    nodestack = [nodeindex]
    while !isempty(nodestack)
        node = tree.nodes[pop!(nodestack)]

        @inbounds for vertex in node.vertices
            if vertices isa Partition # FIXME ugly way to append to a partition
                pushelem!(vertices, vertex)
            else
                push!(vertices, vertex)
            end
        end
        # descend to children (reverse to keep the processing order from the first to the last child)
        append!(nodestack, Iterators.reverse(node.children))
    end
    return vertices
end

"""
    cut(tree::SCCTree, threshold; minsize=1) -> IndicesPartition

Cuts the tree at the given edge weight `threshold` to get the corresponding
strongly connected components of the original graph.

# Keyword arguments
* `minsize`: the minimal number of vertices in the component.
  Smaller connected components are skipped. By default returns components of any size.
"""
cut(tree::SCCTree, threshold::Number; kwargs...) =
    cut!(IndicesPartition(Vector{Int}(undef, nvertices(tree)), Vector{Int}()),
         tree, threshold; kwargs...)

function cut!(comps::IndicesPartition, tree::SCCTree, threshold::Number;
              minsize::Integer=1)
    empty!(comps)
    isempty(tree.nodes) && return comps

    nodestack = [(1, first(tree.nodes).threshold)]
    while !isempty(nodestack)
        nodeindex, parent_threshold = pop!(nodestack)
        if nodeindex < 0 # second processing of the node
            # just add the node vertices (to its own component)
            # without descending to its children as it was already done before
            node = tree.nodes[-nodeindex]
            @inbounds for vertex in node.vertices
                pushelem!(comps, vertex)
                closepart!(comps)
            end
            continue
        end
        node = tree.nodes[nodeindex]
        (nvertices(node) < minsize) && continue # skip small nodes
        notweaker = !isnan(node.threshold) && !isweaker(node.threshold, threshold, rev=tree.rev)
        if notweaker
            #@info "node #$(nodeindex) threshold=$(node.threshold), marking as component #$nextcompindex"
            appendvertices!(comps, tree, nodeindex)
            if (nodeindex == 1) || isnan(parent_threshold) ||
                    isweaker(parent_threshold, threshold, rev=tree.rev)
                #@info "node #$(nodeindex) (size=$(node.nvertices), threshold=$(node.threshold)): parent threshold is $(parent_threshold), so advancing component $(nextcompindex+1)"
                closepart!(comps)
                sort!(comps[end]) # sort the vertices in a component
            end
        else
            if !isempty(node.children)
                if minsize <= 1 && !isempty(node.vertices)
                    # append the vertices of the node after its children
                    push!(nodestack, (-nodeindex, node.threshold))
                end
                # descend to children, so that the children would be processed from first to last
                for child in Iterators.reverse(node.children)
                    push!(nodestack, (child, node.threshold))
                end
            elseif minsize <= 1
                # each vertex is in a component of its own
                @inbounds for vertex in node.vertices
                    pushelem!(comps, vertex)
                    closepart!(comps)
                end
            end
        end
    end
    return comps
end
