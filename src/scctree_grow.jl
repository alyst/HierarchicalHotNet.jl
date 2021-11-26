#  Construction of SCCTree

# Intermediate `SCCTreeNode` node representation used by `SCCSeedling`.
struct SCCSeedlingNode{I}
    threshold::I        # weakest threshold the node is a part of some component
    nvertices::Int      # number of vertices in a component represented by a node
    firstvertex::Int    # graph vertex with the smallest index associated with the node

    children::Vector{Int}   # indices of children nodes
    vertices::Vector{Int}   # indices of graph vertices directly attached to the node
end

# Helper structure for building the `SCCTree`.
# Using indices of adjmtx weights instead of the weights itelsf allows
# simplifying the logic of direct/reversed weights ordering and
# collecting unique weights for each subtree adjacency matrix much more
# efficiently.
mutable struct SCCSeedling{T, I <: Integer, M <: AbstractMatrix}
    rev::Bool

    iadjmtx::M          # original adjmtx with weights replaced by their indices in `weights` vector
    weights::Vector{T}  # unique weights of `adjmtx` sorted from weakest to strongest
    nodes::Vector{SCCSeedlingNode{I}} # nodes of the growing tree
    vertexnodes::Vector{Int} # maps graph verticex to the index of the node it is directly attached to

    stamp::Int # next stamp id to use in sortediweights()
    weight_stamps::Vector{Int} # helper array to stamp values observed in adjmtx

    pools::ObjectPools

    function SCCSeedling{T,I}(adjmtx::AbstractMatrix{T}; rev::Bool = false) where {T, I}
        empty_adjmtx = empty_adjacencymatrix(typeof(adjmtx), I)
        new{T,I,typeof(empty_adjmtx)}(rev, empty_adjmtx,
            Vector{T}(), Vector{SCCSeedlingNode{I}}(),
            Vector{Int}(),
            1, Vector{Int}(),
            ObjectPools())
    end
end

# FIXME better/dynamic than I=Int32?
# a) Int32 is enough for most usecases (matrices up to 46000x46000)
# b) Int32 is different from Int (on most machines),
#    which is good since sizes for indices and weight vectors are different
SCCSeedling(adjmtx::AbstractMatrix{T}, ::Type{I} = Int32; kwargs...) where {T, I} =
    reset!(SCCSeedling{T, I}(adjmtx), adjmtx; kwargs...)

iweighttype(::Type{<:SCCSeedling{<:Any, I}}) where I = I
iweighttype(tree::SCCSeedling) = iweighttype(typeof(tree))

function reset!(tree::SCCSeedling{T,I},
                adjmtx::AbstractMatrix{T};
                skipval::Union{Number, Nothing}=zero(eltype(adjmtx)),
                rev::Bool=false) where {T,I}
    test = EdgeTest{T}(skipval=skipval, rev=rev, threshold=nothing)
    iadjmtx, _ = indexvalues!(tree.iadjmtx, tree.weights, adjmtx, test)
    tree.iadjmtx = iadjmtx
    nv = size(adjmtx, 1)
    tree.rev = rev
    empty!(tree.nodes)
    fill!(resize!(tree.vertexnodes, nv), 0)

    tree.stamp = 1
    fill!(resize!(tree.weight_stamps, length(tree.weights)), 0)

    arraypool(tree.pools, I).borrow_limit = max(10, nv*10)
    arraypool(tree.pools, Int).borrow_limit = max(10, nv*10)
    arraypool(tree.pools, DFSState).borrow_limit = 1
    objpool(tree.pools, IndicesPartition).borrow_limit = max(10, length(tree.weights))

    return tree
end

nvertices(tree::SCCSeedling) = length(tree.vertexnodes)

# using pool seems slower
weightpool(tree::SCCSeedling{T, I}) where {T, I} =
    NoopArrayPool{I}()#arraypool(tree.pools, I)

matrixpool(tree::SCCSeedling{T, I}) where {T, I} =
    NoopArrayPool{I}()#arraypool(tree.pools, I)

# Gets the sorted unique values observed in `arr`.
# To improve the performance, the method uses `superset`, a sorted superset
# of unique values present in `arr`.
function sortediweights(tree::SCCSeedling{T, I}, arr::AbstractArray{I},
                        superset::AbstractVector) where {T, I}
    iweights = empty!(borrow!(weightpool(tree), length(superset)))
    # generate the unique stamp
    stamp = (tree.stamp += 1)
    # stamp all observed values
    @inbounds for iw in arr
        if iw != 0
            tree.weight_stamps[iw] = stamp
        end
    end
    # only look at the superset values and collect those that have the current stamp
    @inbounds for iw in superset
        if tree.weight_stamps[iw] == stamp
            push!(iweights, iw)
        end
    end
    return iweights
end

partition(tree::SCCSeedling, n::Integer; nparts::Integer=n) =
    reset!(borrow!(objpool(tree.pools, IndicesPartition)), n, nparts=nparts)

release!(tree::SCCSeedling, ptn::IndicesPartition) =
    release!(objpool(tree.pools, IndicesPartition), ptn)

# appends one node to the tree for each part `comps`
# `comps` parts represent strongly connected components at `threshold` level
# the indices in `comps` refer to the elements of `subtree`,
# whereas the elements of `subtree` are either existing `tree` nodes (positive indices)
# or yet-unattached vertices of the original graph (negative indices)
function append_components!(tree::SCCSeedling,
                            comps::AbstractPartition,
                            subtree::AbstractVector{Int},
                            threshold::Number)
    for (i, comp) in enumerate(comps)
        if (length(comps) == 1) || (length(comp) > 1)
            # component contains several subcomponents
            # or the nodes have to be attached to the root
            comp_node = length(tree.nodes) + 1
            # distribute subtree references to nodes and vertices
            vertices = Vector{Int}()
            nodes = Vector{Int}()
            @inbounds for j in comp
                ref = subtree[j]
                if ref > 0 # ref is a node
                    insertsorted!(nodes, ref, unique=true)
                elseif ref < 0 # ref is an unattached vertex
                    insertsorted!(vertices, -ref) # vertices are unique
                else
                    error("Ref could not be zero")
                end
                # update subtree to point to the current component node
                subtree[j] = comp_node
            end
            ntotalvertices = length(vertices)
            firstvertex = isempty(vertices) ? nvertices(tree) + 1 : minimum(vertices)
            if !isempty(nodes)
                firstvertex = min(firstvertex, minimum(i -> tree.nodes[i].firstvertex, nodes))
                ntotalvertices += sum(i -> tree.nodes[i].nvertices, nodes)
            end
            push!(tree.nodes, eltype(tree.nodes)(threshold, ntotalvertices, firstvertex,
                                                 nodes, vertices))
            #tree.verbose && @info("New component node (#$comp_node)->[threshold=$threshold]->($(comp_subnodes))")
            @inbounds tree.vertexnodes[vertices] .= comp_node
        end
    end
    return tree
end

"""
    scctree(g::Union{AbstractSimpleWeightedGraph, AbstractMatrix};
            method=:bisect, skipval=0, rev=false) -> SCCTree

Computes the hierarchical decomposition of the weighted directed graph `g` into
*strongly connected components*.

Supports weighted graph objects as well as their adjacency matrix representations.

# Keyword arguments
* `method::Symbol` (defaults to `:bisect`): the method for partitioning.
   The supported methods are `:bisect` (the fastest) and `:bottomup` (slow,
   but simpler).
* `skipval::Number` (defaults to zero): what value of the adjacency matrix
   should be treated as "no edge".
* `rev::Bool` (defaults to `false`): if `true`, bigger edges weights are considered
   weaker than smaller ones
"""
scctree(g::AbstractSimpleWeightedGraph; kwargs...) =
    scctree(Graphs.weights(g), skipval=zero(eltype(Graphs.weights(g))); kwargs...)

function scctree(adjmtx::AbstractMatrix; method::Symbol=:bisect,
                 skipval::Union{Number, Nothing} = zero(eltype(adjmtx)),
                 rev::Bool=false,
                 seedling::Union{SCCSeedling, Nothing} = nothing,
                 kwargs...)
    nvertices = size(adjmtx, 1)
    check_square(adjmtx, "Adjacency matrix")
    tree = isnothing(seedling) ? SCCSeedling(adjmtx, skipval=skipval, rev=rev) :
           reset!(seedling, adjmtx, skipval=skipval, rev=rev)
    if method == :bisect
        scctree_bisect!(tree; kwargs...)
    elseif method == :bottomup
        scctree_bottomup!(tree; kwargs...)
    else
        throw(ArgumentError("Unsupported graph hierarchical clustering method: $method"))
    end
    return SCCTree(tree)
end

# Constructs SCCTree using the "bottom-up" algorithm:
# Detect strongly connected components of the graph starting from the strongest
# threshold and gradually merge the components into bigger ones as the threshold
# gets weaker.
function scctree_bottomup!(tree::SCCSeedling; verbose::Bool=false)
    @assert isempty(tree.nodes)
    vertex_roots = collect(-1:-1:-nvertices(tree)) # start with each vertex being a root of a subtree
    comps = IndicesPartition(nvertices(tree), nparts=nvertices(tree)) # reusable partition of vertices into components
    ncomps = length(comps)
    I = iweighttype(tree)
    for threshold in length(tree.weights):-1:1
        threshold_cut = EdgeTest{I}(skipval=0, threshold=threshold, rev=false)
        strongly_connected_components!(comps, tree.iadjmtx, threshold_cut, tree.pools)
        (length(comps) == ncomps) && continue # same components as for the stronger threshold
        ncomps = length(comps) # new partition
        verbose && @info("$(ncomps) components detected at threshold=$threshold")
        append_components!(tree, comps, vertex_roots, threshold)
        (ncomps == 1) && break # all nodes are already in the single connected component
    end
    if (ncomps > 1) || (nvertices(tree) == 1) # multiple connected components or single vertex
        # add the no-threshold root to join disconnected components
        append_components!(tree, reset!(comps, nparts=1), vertex_roots, 0)
    end

    return tree
end

# Constructs SCCTree using the "bisection" algorithm:
# Detect strongly connected components using the median threshold, then
# recursively repeat the procedure for the resulting subgraphs and for
# the "condensed" version of the original graph, where the connected components
# are replaced by single nodes.
function scctree_bisect!(tree::SCCSeedling; verbose::Bool=false)
    @assert isempty(tree.nodes)
    scctree_bisect_subtree!(tree, tree.iadjmtx, collect(-1:-1:-nvertices(tree)), 0, 0,
                            eachindex(tree.weights), verbose=verbose)
    return tree
end

# The "worker" function of "bisection" algorithm.
# Recursively calls itself for the subgraphs and the "condensed supergraph".
function scctree_bisect_subtree!(tree::SCCSeedling, adjmtx::AbstractMatrix{<:Integer},
                                 subtree::AbstractVector{Int},
                                 subtree_threshold::Integer, nodes_threshold::Integer,
                                 parent_weights::AbstractVector;
                                 verbose::Bool=false) where T
    verbose && @info "scctree_bisect_range!($(subtree) nodes, subtree_threshold=$subtree_threshold, nodes_threshold=$nodes_threshold)"
    weights = sortediweights(tree, adjmtx, parent_weights)
    intpool = arraypool(tree.pools, Int)
    I = iweighttype(tree)
    wpool = weightpool(tree)
    mtxpool = matrixpool(tree)
    # release!(tree.iweights_pool, weights) should be called before returning from this function
    if isempty(weights) # (condensed) graph with no edges
        # i.e. the original graph has multple connected components
        @assert subtree_threshold == 0 # no edge weights could be only in the root (condensed) graph
        verbose && @info("No edges, creating the root node to join connected components")
        onecomp = partition(tree, length(subtree), nparts=1)
        append_components!(tree, onecomp, subtree, subtree_threshold)
        release!(wpool, weights)
        release!(tree, onecomp)
        return length(tree.nodes)
    end

    # find the weight indices corresponding to subtree and nodes thresholds
    if subtree_threshold != 0
        subtree_lev = searchsortedfirst(weights, subtree_threshold)
        # correct the subtree threshold for this component
        # it might be stronger than estimated on the upper level
        if subtree_threshold < weights[subtree_lev]
            subtree_threshold = weights[subtree_lev]
        end
    else
        subtree_lev = 1
    end
    @assert 1 <= subtree_lev <= length(weights) "subtree_threshold=$subtree_threshold outside of ($(first(weights))..$(last(weights)))"

    if nodes_threshold != 0
        nodes_lev = searchsortedfirst(weights, nodes_threshold)
        (nodes_lev > length(weights)) && (nodes_lev = length(weights))
    else
        nodes_lev = length(weights)
    end
    @assert 1 <= nodes_lev <= length(weights) "nodes_threshold=$nodes_threshold outside of ($(first(weights))..$(last(weights)))"
    verbose && @info("subtree_threshold=weights[$subtree_lev]=$(weights[subtree_lev]), nodes_threshold=weights[$nodes_lev]=$(weights[nodes_lev]) of $weights")

    nnodes = length(subtree)
    # initialize with trivial components
    comps = partition(tree, nnodes, nparts=nnodes)
    # search for non-trivial dissection
    bisect_threshold = NaN
    while true
        if nodes_lev <= subtree_lev + 1 # all bisections checked
            if subtree_threshold == 0 # it's the root subtree (all others should have subtree threshold)
                # check if graph has multiple connected components
                subtree_cut = EdgeTest{I}(skipval=0, threshold=weights[subtree_lev], rev=false)
                strongly_connected_components!(comps, adjmtx, subtree_cut, tree.pools)
                if length(comps) == 1
                    # it's a connected graph, so assign the threshold to the root
                    subtree_threshold = weights[subtree_lev]
                elseif length(comps) < nnodes # nontrivial SCCs!
                    bisect_lev = subtree_lev
                    bisect_threshold = weights[bisect_lev]
                    break
                end
            end
            # merge components into a single node and return its index
            verbose && @info("Creating a node for subtree $subtree at subtree_lev=$subtree_lev, nodes_lev=$nodes_lev")
            # reset to a single component
            append_components!(tree, reset!(comps, nparts=1), subtree, subtree_threshold)
            release!(wpool, weights)
            release!(tree, comps)
            return length(tree.nodes)
        end

        # bisect current sub-scctree (use the weight from the middle as SCC threshold)
        bisect_lev = fld1(subtree_lev + nodes_lev, 2)
        @assert subtree_lev < bisect_lev < nodes_lev
        bisect_threshold = weights[bisect_lev]
        bisect = EdgeTest{I}(skipval=0, threshold=bisect_threshold, rev=false)
        strongly_connected_components!(comps, adjmtx, bisect, tree.pools)
        verbose && @info("Identified $ncomps strongly connected component(s) at threshold[$bisect_lev]=$bisect_threshold: $index2comp")
        if length(comps) == 1 # single component
            verbose && @info("Single SCC")
            subtree_lev = bisect_lev # +1
            subtree_threshold = bisect_threshold
        elseif length(comps) == nnodes # no non-trivial components
            verbose && @info("Trivial SCCs")
            nodes_lev = bisect_lev #-1
            nodes_threshold = bisect_threshold
        else # multiple non-trivial components
            break
        end
    end

    # multiple components
    verbose && @info("Building subtrees of each of $ncomps SCCs")
    # build a graph of components relationships
    comp_roots = borrow!(intpool, 0)  # nodes of the roots of subtrees for each component
    comp_subtree = borrow!(intpool, 0) # reusable vector for component subtree refs
    # recurse into each of multiple components
    for (i, comp_indices) in enumerate(comps)
        # build the subtree for the current component
        if length(comp_indices) > 1 # recursively cluster i-th component
            verbose && @info("scctree_scc!(component #$i)")
            comp_adjmtx = subgraph_adjacencymatrix(adjmtx, comp_indices, intpool)
            # recode indices in the current subtree into node/vertex refs
            sizehint!(empty!(comp_subtree), length(comp_indices))
            @inbounds for i in comp_indices
                push!(comp_subtree, subtree[i])
            end
            comp_root = scctree_bisect_subtree!(tree, comp_adjmtx, comp_subtree,
                                                bisect_threshold, nodes_threshold,
                                                weights, verbose=verbose)
        else # single node/vertex, don't recurse
            comp_root = subtree[comp_indices[1]]
        end
        push!(comp_roots, comp_root)
    end
    release!(intpool, comp_subtree)

    # build a graph of components relationships
    comps_adjmtx = condense!(borrow!(mtxpool, (length(comps), length(comps))),
                             adjmtx, comps, zerodiag=true)
    # cluster it and attach the resulting subtree to the current root node
    verbose && @info("scctree_scc!(condensed components graph)")
    res = scctree_bisect_subtree!(tree, comps_adjmtx, comp_roots,
                                  subtree_threshold, bisect_threshold,
                                  weights, verbose=verbose)
    release!(mtxpool, comps_adjmtx)
    release!(wpool, weights)
    release!(intpool, comp_roots)
    release!(tree, comps)
    return res
end
