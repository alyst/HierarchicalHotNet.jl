using SparseArrays

const VertexId = Int # ID of the original graph vertex
const NodeId = Int   # ID of the SCC graph vertex (component id)
const NodeIx = Int   # index of the SCC graph vertex in adjmtx (its column/row)
const FlowId = Int   # ID of the source flow
const PtId = Int     # ID of the source flow subgraph vertex
const IWeight = Int
const FlowAdjMatrix = SparseMatrixCSC{IWeight, VertexId}

# vertex of the SourceFlowsPeel graph
struct FlowVertexInfo
    nodeid::NodeId  # id of the node (in SCCTree)
    issink::Bool    # whether the node contains the sinks
end

# SCCTree node that participates in the source -> sink flow(s)
mutable struct FlowsPeelNode
    id::NodeId                  # node id, same as SCCNode.id
    index::NodeIx               # node index in SCCTreeFlowsPeelingIterator.iadjmtx
    vertices::Vector{VertexId}  # vertices of the node
    sources::Vector{VertexId}   # source vertices contained in the node
    sinks::Vector{VertexId}     # sink vertices contained in the node

    flows::Set{FlowId}          # flows that go through this node
end

# source/intermediate node -> sink path (in SourceFlowsPeel graph) info
struct FlowSinkPath
    sinkid::NodeId      # sink node id (0 = undefined)
    len::Int            # path length
    minweight::IWeight  # minimal weight of the edge along the path

    FlowSinkPath(sinkid::NodeId = 0, len::Integer = 0, minweight::Number = 0) =
        new(sinkid, len, minweight)
end

isdefined(state::FlowSinkPath) = state.sinkid > 0

# merges the two pathes to the same destination:
# select the shortest, and with maxweight if both have similar length
function best(ainfo::FlowSinkPath, binfo::FlowSinkPath)
    @assert isdefined(binfo)
    if isdefined(ainfo)
        ainfo.sinkid == binfo.sinkid || throw(ArgumentError("Different destination sinks for merged paths"))
        return FlowSinkPath(ainfo.sinkid, min(ainfo.len, binfo.len),
                        ainfo.len == binfo.len ? max(ainfo.minweight, binfo.minweight) :
                        ainfo.len < binfo.len ? ainfo.minweight : binfo.minweight)
    else
        return binfo
    end
end

addedge(info::FlowSinkPath, w::IWeight) =
    FlowSinkPath(info.sinkid, info.len + 1, info.len > 0 ? min(info.minweight, w) : w)

# flow from the given source to sinks
mutable struct SourceFlowsPeel
    id::FlowId                  # id of the flow
    srcid::NodeId               # id of the source node

    dirty::Bool                 # whether an update is required
    minweight::IWeight          # minimal indexed weight of the edge along the path from source to sink
    iadjmtx::FlowAdjMatrix      # indexed adjacent matrix for a Direct Acyclic Graph induced by the flow
    vtxinfo::Vector{FlowVertexInfo}   # graph vertex to node mapping
    paths::Vector{FlowSinkPath} # info about all source -> sink paths within the flow
end

# hash of the paths
function pathshash(peel::SourceFlowsPeel, nodes::AbstractDict{NodeId, FlowsPeelNode})
    res = hash(length(peel.paths))
    for path in peel.paths
        res = hash(nodes[path.sinkid].sinks, hash(path.minweight, hash(path.len, res)))
    end
    return res
end

Base.copy(flow::SourceFlowsPeel) =
    SourceFlowsPeel(flow.id, flow.srcid, flow.dirty, flow.minweight, flow.iadjmtx,
                    flow.vtxinfo, similar(flow.paths, 0))

# immutable part of SCCTreeFlowsPeelingIterator
struct SCCTreeFlowsPeeling{T}
    tree::SCCTree{T}

    vertex_iadjmtx::FlowAdjMatrix   # original adjmtx with weights replaced by their indices in `weights` vector
    weights::Vector{T}              # unique weights of `adjmtx` sorted from weakest to strongest
    ithresholds::Vector{IWeight}    # indexed version of SCCTree.thresholds (indices in weights array that correspond to thresholds that result in different tree cuts)

    sources::Set{VertexId}          # ids of the source vertices in the original graph
    sinks::Set{VertexId}            # ids of the sink vertices in the original graph

    node_indices::Dict{NodeId, Int}    # maps ids of SCCTree node to the index of the corresponding part in node_vertices_parts
    node_vertices::IndicesPartition    # each part are the indices of the vertices of some SCCTree node

    function SCCTreeFlowsPeeling(tree::SCCTree{T}, adjmtx::AbstractMatrix,
                                sources::AbstractVector{VertexId},
                                sinks::AbstractVector{VertexId};
                                sortvertices::Bool = false,
                                verbose::Bool = false) where T
        iadjmtx, weights = indexvalues(IWeight, adjmtx, EdgeTest{eltype(adjmtx)}(rev=tree.rev))
        ithresholds = Vector{IWeight}(undef, length(tree.thresholds))
        # index the thresholds
        for (i, thresh) in enumerate(tree.thresholds)
            ithresholds[i] = ithresh = searchsortedfirst(weights, thresh)
            if i > 1
                lastithresh = ithresholds[i - 1]
                (lastithresh >= ithresh) && @warn "threshold[$i]=$thresh indexed the same ($ithresh) as threshold[$(i-1)]=$(tree.thresholds[i-1]) ($lastithresh)"
            end
        end
        spiadjmtx = sparse(iadjmtx)
        verbose && @info "SCCTreeFlowsPeels: $(size(spiadjmtx, 1)) vertice(s), $(nnz(spiadjmtx)) edge(s) with $(length(ithresholds)) threshold(s) imported"

        return new{T}(tree, spiadjmtx, weights, ithresholds, Set(sources), Set(sinks),
                      cache_node_vertices(tree, sortvertices=sortvertices)...)
    end
end

nthresholds(peeling::SCCTreeFlowsPeeling) = length(peeling.ithresholds)

function cache_node_vertices(tree::SCCTree; sortvertices::Bool=false)
    vtxparts = IndicesPartition()
    nodeindices = Dict{NodeId, Int}()

    nodestack = Vector{Tuple{NodeId, Int}}()
    nvertices(tree) > 0 && push!(nodestack, (1, 0))
    while !isempty(nodestack)
        id, childix = last(nodestack)
        node = tree.nodes[id]
        if childix < length(node.children)
            # process the next child node
            childix += 1
            nodestack[end] = (id, childix)
            push!(nodestack, (node.children[childix], 0))
        else
            # all children processed, fill the node
            append!(elems(vtxparts), node.vertices)
            for childid in node.children
                childix = nodeindices[childid]
                childrange = partrange(vtxparts, childix)
                nelms = nelems(vtxparts)
                resize!(elems(vtxparts), nelms + length(childrange))
                unsafe_copyto!(elems(vtxparts), nelms + 1, elems(vtxparts), first(childrange), length(childrange))
            end
            closepart!(vtxparts)
            sortvertices && sort!(vtxparts[end]) # sort vertices
            # add the node
            nodeindices[id] = length(vtxparts)
            pop!(nodestack)
        end
    end
    @assert length(nodeindices) == length(vtxparts) == length(tree.nodes)
    return nodeindices, vtxparts
end

node_vertices(peeling::SCCTreeFlowsPeeling, id::NodeId) =
    peeling.node_vertices[peeling.node_indices[id]]

# Implements efficient enumeration of source->sink flows
# for each SCCTree cutting threshold.
# Iterates through the cutting thresholds (from weakest to strongest),
# remembers the flows at given threshold and uses them to more efficiently
# enumerate the flows at the next threshold.
mutable struct SCCTreeFlowsPeelingIterator{T <: Number}
    peeling::SCCTreeFlowsPeeling{T} # immutable part

    ithresholdpos::Int      # current position at peeling.ithresholds
    ithreshold::IWeight     # current indexed threshold

    iadjmtx::FlowAdjMatrix  # indexed adjacency matrix for a graph between the current nodes
    nodeix2id::Vector{NodeId} # map from the node graph vertex index (col/row in iadjmtx) to the node id in SCCTree

    nodes::Dict{NodeId, FlowsPeelNode}  # current active (ones that participate in any flow) SCCTree nodes
    flows::Dict{FlowId, SourceFlowsPeel} # current flows (from some active source node to any sink it can reach)
    lastflowid::FlowId      # the id of the last generated flow

    verbose::Bool           # whether to generate diagnostic messagees

    SCCTreeFlowsPeelingIterator(peeling::SCCTreeFlowsPeeling{T}; verbose::Bool=false) where T =
        new{T}(peeling, 0, 0, spdiagm(0 => IWeight[]), NodeId[],
               Dict{NodeId, FlowsPeelNode}(), Dict{FlowId, SourceFlowsPeel}(), 0, verbose)
end

function eachflowspeel(tree::SCCTree, adjmtx::AbstractMatrix,
                       sources::AbstractVector{VertexId},
                       sinks::AbstractVector{VertexId};
                       sortvertices::Bool=false,
                       verbose::Bool = false
)
    peeling = SCCTreeFlowsPeeling(tree, adjmtx, sources, sinks, sortvertices=sortvertices, verbose=verbose)
    return SCCTreeFlowsPeelingIterator(peeling, verbose=verbose)
end

# the value returned by SCCTreeFlowsPeelingIterator
struct SCCTreeFlowsPeels{T <: Number}
    ithreshold::IWeight     # index of the current SCCTree cutting threshold
    weights::Vector{T}      # reference to SCCTreeFlowsPeeling.weights
    nodes::Dict{NodeId, FlowsPeelNode}  # active SCCTree nodes
    flows::Dict{FlowId, SourceFlowsPeel} # active source->sink flows
    modified::Bool          # the flows modified in comparison to the previous iteration

    SCCTreeFlowsPeels(it::SCCTreeFlowsPeelingIterator{T}, modified::Bool) where T =
        new{T}(it.ithreshold, it.peeling.weights, it.nodes, it.flows, modified)
end

threshold(peels::SCCTreeFlowsPeels) = peels.ithreshold > 0 ? peels.weights[peels.ithreshold] : NaN
thresholdindex(peels::SCCTreeFlowsPeels) = peels.ithreshold

threshold(it::SCCTreeFlowsPeelingIterator) = it.ithreshold > 0 ? it.peeling.weights[it.ithreshold] : NaN
thresholdindex(it::SCCTreeFlowsPeelingIterator) = it.ithreshold

Base.IteratorEltype(::Type{<:SCCTreeFlowsPeelingIterator}) = Base.HasEltype()
Base.eltype(::Type{SCCTreeFlowsPeelingIterator{T}}) where T = SCCTreeFlowsPeels{T}
Base.eltype(it::SCCTreeFlowsPeelingIterator) = eltype(typeof(it))

Base.IteratorSize(::Type{<:SCCTreeFlowsPeelingIterator}) = Base.HasLength()
Base.length(it::SCCTreeFlowsPeelingIterator) = nthresholds(it.peeling)

edgetest(it::SCCTreeFlowsPeelingIterator) =
    EdgeTest{IWeight}(rev=false, threshold=it.ithreshold > 0 ? it.ithreshold : nothing)

# filter flow graph by the edgemask
function induce!(flow::SourceFlowsPeel, edgemask::AbstractVector{Bool},
                 it::SCCTreeFlowsPeelingIterator,
                 pools::Union{ObjectPools, Nothing} = nothing)
    it.verbose && @info "Flow #$(flow.id): filtering unused edges and vertices"
    isempty(flow.paths) && return flow # skip induction if no paths (flow gets deleted)

    nv = size(flow.iadjmtx, 1)
    boolpool = arraypool(pools, Bool)
    # which vertices are in the new graph
    vused = fill!(borrow!(boolpool, nv), false)
    @inbounds for v in 1:nv
        for eix in nzrange(flow.iadjmtx, v)
            if edgemask[eix]
                vused[v] = true
                vused[rowvals(flow.iadjmtx)[eix]] = true
            end
        end
    end
    # use the source self-loop, which might not have an edge
    for path in flow.paths
        if path.sinkid == flow.srcid
            vused[1] = true
        end
    end
    @assert vused[1] # source vertex is always used
    it.verbose && @info "Flow #$(flow.id): $(count(vused))/$(length(vused)) vertices used"

    # update the column ranges
    colptrs = sizehint!(Vector{Int}(), count(vused) + 1)
    usedvtx_indices = fill(0, length(vused))
    push!(colptrs, 1)
    @inbounds for (v, isused) in enumerate(vused)
        if isused
            vlast = last(colptrs)
            for eix in nzrange(flow.iadjmtx, v)
                edgemask[eix] && (vlast += 1)
            end
            usedvtx_indices[v] = length(colptrs) # reindex the used node
            push!(colptrs, vlast)
        else # unlink unused node from the flow
            nodeid = flow.vtxinfo[v].nodeid
            it.verbose && @info "Node #$nodeid: detaching flow #$(flow.id) (unused)"
            node = it.nodes[nodeid]
            @assert node.id == nodeid
            pop!(node.flows, flow.id)
        end
    end

    # replace the graph matrix
    @assert all(>(0), nonzeros(flow.iadjmtx)[edgemask])
    flow.iadjmtx = FlowAdjMatrix(length(colptrs) - 1, length(colptrs) - 1, colptrs,
                                usedvtx_indices[rowvals(flow.iadjmtx)[edgemask]],
                                nonzeros(flow.iadjmtx)[edgemask])
    if nnz(flow.iadjmtx) > 0
        flow.minweight = minimum(nonzeros(flow.iadjmtx))
        @assert flow.minweight > 0
    else
        flow.minweight = 0
    end
    # update vertex info
    flow.vtxinfo = flow.vtxinfo[vused]
    flow.dirty = false
    release!(boolpool, vused)
    return flow
end

# replace the nodes of SourceFlowsPeel with the graphs of their internal structure
function local_nodes_expand!(
        flow::SourceFlowsPeel, expanded::AbstractDict{NodeId, IndicesPartitionPart},
        it::SCCTreeFlowsPeelingIterator, reattach::Bool,
        etest = edgetest(it),
        pools::Union{ObjectPools, Nothing} = nothing;
        verbose::Bool = it.verbose
)
    # build the mapping from the new local node index to the old local index
    # the new local node index to the global index
    newinfo = sizehint!(similar(flow.vtxinfo, 0), length(flow.vtxinfo))
    newix2loldix = sizehint!(similar(Vector{PtId}(), 0), length(flow.vtxinfo))
    newix2gix = sizehint!(similar(Vector{PtId}(), 0), length(flow.vtxinfo))

    for (i, vinfo) in enumerate(flow.vtxinfo)
        newids = get(expanded, vinfo.nodeid, nothing)
        oldnode = it.nodes[vinfo.nodeid]
        if newids !== nothing # the old node was expanded
            if reattach
                verbose && @info "Expanded Node #$(oldnode.id): detaching flow #$(flow.id)"
                pop!(oldnode.flows, flow.id) # detach the old flow from the old node
            end
            for newid in newids
                node = it.nodes[newid]
                push!(newinfo, FlowVertexInfo(newid, !isempty(node.sinks)))
                push!(newix2gix, node.index)
                push!(newix2loldix, -i)
                if reattach
                    verbose && @info "Subnode #$(node.id): reattaching flow #$(flow.id)"
                    push!(node.flows, flow.id) # attach flow to the new node
                end
            end
        else # the old node is kept
            push!(newinfo, vinfo)
            push!(newix2gix, oldnode.index)
            push!(newix2loldix, i)
        end
    end
    # the first vertex is the source (although the node id may change)
    if flow.srcid != newinfo[1].nodeid
        verbose && @info "Flow #$(flow.id): changing the source from $(flow.srcid) to $(newinfo[1].nodeid)"
        flow.srcid = newinfo[1].nodeid
        @assert !reattach || (flow.id in it.nodes[flow.srcid].flows) "Flow #$(flow.id) is not attached to its new source node #$(flow.srcid)"
        @assert !isempty(it.nodes[flow.srcid].sources)
    end

    # from the global index to new local index
    gix2newix = fill(0, length(it.nodeix2id))
    @inbounds for (newix, gix) in enumerate(newix2gix)
        gix2newix[gix] = newix
    end

    # generate the new adjmtx
    newcolptrs = sizehint!(Vector{Int}(), length(newix2loldix)+1)
    push!(newcolptrs, 1)
    newrowvals = Vector{PtId}()
    newweights = Vector{IWeight}()
    colrows = Vector{Pair{PtId, IWeight}}() # to collect nonzero rows for a current column
    voldedges = Vector{IWeight}(undef, length(flow.vtxinfo))
    lastvold = 0

    for (vnew, vold_) in enumerate(newix2loldix)
        vexpanded = vold_ < 0
        vold = ifelse(vexpanded, -vold_, vold_)
        # build the vold edge mask
        @inbounds if lastvold != vold
            fill!(voldedges, 0)
            for (uold, iw) in outedges(flow.iadjmtx, vold, etest)
                voldedges[uold] = iw
            end
            lastvold = vold
        end
        vgedges = nzrange(it.iadjmtx, newix2gix[vnew])
        empty!(colrows)
        for eix in vgedges
            iw = nonzeros(it.iadjmtx)[eix]
            isvalidedge(iw, etest) || continue
            gix = rowvals(it.iadjmtx)[eix]
            unew = gix2newix[gix]
            (unew > 0) || continue
            uold_ = newix2loldix[unew]
            uold = ifelse(uold_ < 0, -uold_, uold_)
            if (voldedges[uold] > 0) || # non-expanded preserved
               (vexpanded && vold == uold) # expanded diagonal
               push!(colrows, unew => iw)
            end
        end
        sort!(colrows, by=first) # sort rows by their index
        for (row, w) in colrows
            push!(newrowvals, row)
            push!(newweights, w)
        end
        push!(newcolptrs, length(newrowvals) + 1)
    end
    @assert length(newcolptrs) == length(newinfo) + 1

    flow.iadjmtx = FlowAdjMatrix(length(newinfo), length(newinfo), newcolptrs, newrowvals, newweights)
    flow.vtxinfo = newinfo
    flow.dirty = true # the updated flow graph has to be processed
    return flow
end

# expand the flows that go through the expanded nodes
function expand_flows!(
        it::SCCTreeFlowsPeelingIterator,
        expanded::AbstractDict{NodeId, IndicesPartitionPart},
        pools::Union{ObjectPools, Nothing} = nothing;
        verbose::Bool = it.verbose
)
    flowids_to_expand = Set{FlowId}()
    for nodeid in keys(expanded)
        # expand the local nodes of every flow that is incident with nodeid
        union!(flowids_to_expand, it.nodes[nodeid].flows)
    end

    swapfirstelem!(::Nothing, pos::Integer) = nothing

    function swapfirstelem!(v::AbstractVector, pos::Integer)
        if pos != 1
            v[1], v[pos] = v[pos], v[1]
        end
        return v
    end

    for flowid in flowids_to_expand
        flow = it.flows[flowid]
        newsources = get(expanded, flow.srcid, nothing)
        # we fork the flow if the sources of the expanded node are distributed into more than one subnodes
        newsourcepos = !isnothing(newsources) ?
                findall(subid -> !isempty(it.nodes[subid].sources), newsources) : nothing
        if !isnothing(newsourcepos) && length(newsourcepos) > 1
            verbose && @info "Flow #$(flow.id): forking"
            # remove the old flow from all its nodes
            for vinfo in flow.vtxinfo
                verbose && @info "Node #$(vinfo.nodeid): detaching flow #$(flow.id) (old)"
                pop!(it.nodes[vinfo.nodeid].flows, flow.id)
            end
            pop!(it.flows, flow.id) # and from the registry
            # fork the flow for each subnode that still contains the source vertices
            for srcpos in newsourcepos
                subflow = copy(flow)
                subflow.id = (it.lastflowid += 1)
                it.flows[subflow.id] = subflow
                swapfirstelem!(newsources, srcpos) # make newnode the first node in newnodes
                local_nodes_expand!(subflow, expanded, it, false, edgetest(it), pools, verbose=verbose)
                swapfirstelem!(newsources, srcpos) # revert the order
                # reference the new flow from all its nodes
                for vinfo in subflow.vtxinfo
                    verbose && @info "Node #$(vinfo.nodeid): attaching flow #$(subflow.id) (new)"
                    push!(it.nodes[vinfo.nodeid].flows, subflow.id)
                end
            end
        else
            # make newnode the first node in newnodes
            verbose && @info "Expanding flow #$(flow.id)"
            srcpos = !isnothing(newsourcepos) ? newsourcepos[1] : 1
            swapfirstelem!(newsources, srcpos)
            local_nodes_expand!(flow, expanded, it, true, edgetest(it), pools, verbose=verbose)
            swapfirstelem!(newsources, srcpos)
        end
    end
    return it
end

# traverse the graph removing flow edges that are below minweight
# and updating source -> sink paths
function update!(flow::SourceFlowsPeel,
                 it::SCCTreeFlowsPeelingIterator,
                 etest = edgetest(it),
                 pools::Union{ObjectPools, Nothing} = nothing)
    nv = size(flow.iadjmtx, 1)
    nv == size(flow.iadjmtx, 2) ||
        throw(DimensionMismatch("Adjacency matrix must be square ($(size(flow.iadjmtx)) given)"))
    ne = length(nonzeros(flow.iadjmtx))

    boolpool = arraypool(pools, Bool)
    visited = fill!(borrow!(boolpool, nv), false) # DFS-visited vertices
    haspathtosink = fill!(borrow!(boolpool, nv), false)
    eused = fill!(borrow!(boolpool, ne), false) # flags edges used to reach from the source to sinks
    sinkpathpool = arraypool(pools, FlowSinkPath)
    sinkpath = fill!(borrow!(sinkpathpool, (nv, nv)), FlowSinkPath()) # whether the sink u (row) is reachable from v (col)
    dfspool = arraypool(pools, DFSState)
    dfs_stack = borrow!(dfspool)
    empty!(flow.paths)

    # updates paths to sinks through vertex v using the paths from u and edge u->v weight (uv)
    function updatesinkpaths(v::VertexId, u::VertexId, w::IWeight)
        haspathtosink[u] || return nothing
        vpaths = view(sinkpath, :, v) # paths from v
        upaths = view(sinkpath, :, u) # paths from u
        for (i, uipath) in enumerate(upaths)
            if isdefined(uipath)
                @assert haspathtosink[i]
                vuipath = addedge(uipath, w)
                vpaths[i] = best(vpaths[i], vuipath)
                #it.verbose && @info "adding: $v -> $u ($w) to $u ~> $i ($(uistate.minweight)): $(vstates[i].minweight)"
                haspathtosink[v] = true
            end
        end
        return nothing
    end

    src = 1 # the first vertex is always source
    @assert flow.vtxinfo[src].nodeid == flow.srcid
    push!(dfs_stack, (src, 0))
    while !isempty(dfs_stack)
        v, edgeix = dfs_stack[end]
        @assert edgeix >= 0
        vinfo = flow.vtxinfo[v]
        if !visited[v]
            @assert edgeix == 0 # not started iterating if not visited
            visited[v] = true
            if vinfo.issink
                sinkpath[v, v] = FlowSinkPath(vinfo.nodeid, 0, 0) # min weight is updated later when iterating
                haspathtosink[v] = true
            end
        end
        edgeit = outedges(flow.iadjmtx, v, etest)
        u = 0
        while true
            iter = iterate(edgeit, edgeix)
            isnothing(iter) && break
            (i, iw), edgeix = iter
            if i == v # update self-loop path FlowInfo
                vinfo.issink && (sinkpath[v, v] = FlowSinkPath(vinfo.nodeid, 0, iw))
            elseif visited[i] && haspathtosink[i]
                # since we require that flowgraph is DAG
                # vertex i should not be in the current dfs stack
                # all vertices reachable from i should also be reachable from v via edgeix
                updatesinkpaths(v, i, iw)
                eused[edgeix] = true
            elseif !visited[i]
                u = i # select v -> i edge
                # update v edges iterator
                dfs_stack[end] = v, edgeix
                break
            end
        end
        if u > 0 # follow v->u Diedge
            push!(dfs_stack, (u, 0))
            @assert !haspathtosink[u] # never traveled from u before
        else # no more outedges to travel from v, backtrack from v
            pop!(dfs_stack)
            if !isempty(dfs_stack)
                # add the prev -> v edge to the graph
                # and update the paths from prev
                prev, prevedgeix = dfs_stack[end]
                if haspathtosink[v]
                    eused[prevedgeix] = true
                    updatesinkpaths(prev, v, nonzeros(flow.iadjmtx)[prevedgeix])
                end
            end
        end
    end
    # collect all detected source -> sinks paths
    if haspathtosink[src]
        for (i, destinfo) in enumerate(view(sinkpath, :, src))
            if isdefined(destinfo)
                @assert flow.vtxinfo[i].issink
                push!(flow.paths, destinfo)
            end
        end
    end
    it.verbose && @info "Flow #$(flow.id): found $(length(flow.paths)) path(s) using $(count(eused))/$(length(eused)) edge(s)"

    induce!(flow, eused, it, pools)

    release!(sinkpathpool, sinkpath)
    release!(boolpool, visited)
    release!(boolpool, eused)
    release!(dfspool, dfs_stack)

    return flow
end

function collect_subnodes(node::FlowsPeelNode, thresh::Number, peeling::SCCTreeFlowsPeeling)
    tree = peeling.tree
    (node.id > length(tree.nodes)) && return nothing # it's a vertex
    # skip nodes that are stronger than the threshold
    node_thresh = tree.nodes[node.id].threshold
    (node_thresh >= thresh && !tree.rev ||
     node_thresh <= thresh && tree.rev) && return nothing

    sccnode = tree.nodes[node.id]
    subnode2vertices = IndicesPartition()
    sizehint!(elems(subnode2vertices), nvertices(sccnode))
    subnodeids = Vector{NodeId}()
    for vertex in sccnode.vertices
        pushelem!(subnode2vertices, vertex)
        closepart!(subnode2vertices)
        push!(subnodeids, length(tree.nodes) + vertex) # vertex without component
    end
    @inbounds for childid in sccnode.children
        nelms = nelems(subnode2vertices)
        childvtxs = node_vertices(peeling, childid)
        childrange = parentindices(childvtxs)[1]
        resize!(elems(subnode2vertices), nelms + length(childrange))
        unsafe_copyto!(elems(subnode2vertices), nelms + 1, parent(childvtxs), first(childrange), length(childrange))
        closepart!(subnode2vertices)
        push!(subnodeids, childid)
    end
    return subnodeids => subnode2vertices
end

function collect_expanded_nodes!(it::SCCTreeFlowsPeelingIterator)
    # cut the nodes by the threshold
    expanded = Dict{NodeId, IndicesPartitionPart}()
    nodes_subnodes = IndicesPartition()
    thresh = threshold(it)
    for node in values(it.nodes)
        subnodes = collect_subnodes(node, thresh, it.peeling)
        if subnodes !== nothing
            subids, subnode2vertices = subnodes
            @assert length(subids) > 1
            # create the new nodes
            for (subid, subvertices) in zip(subids, subnode2vertices)
                @assert subid > 0
                subsources = Vector{VertexId}()
                subsinks = Vector{VertexId}()
                for v in subvertices
                    (v in it.peeling.sources) && push!(subsources, v)
                    (v in it.peeling.sinks) && push!(subsinks, v)
                end
                it.nodes[subid] = FlowsPeelNode(subid, 0, subvertices, subsources, subsinks, Set{FlowId}())
                pushelem!(nodes_subnodes, subid)
            end
            closepart!(nodes_subnodes)
            expanded[node.id] = nodes_subnodes[end]
        end
    end
    return expanded
end

# expand the nodes adjacency matrix
# given how the nodes are being expanded
function expand_adjmatrix!(it::SCCTreeFlowsPeelingIterator, expanded::AbstractDict{NodeId, IndicesPartitionPart})
    # reindex the nodes
    newix2id = sizehint!(Vector{NodeId}(), length(it.nodeix2id))
    for id in it.nodeix2id
        subnodeids = get(expanded, id, nothing)
        if subnodeids === nothing
            haskey(it.nodes, id) || continue # skip removed nodes
            push!(newix2id, id)
        else
            append!(newix2id, subnodeids)
        end
    end
    # update the node and vertex indices
    vtx_adjmtx = it.peeling.vertex_iadjmtx
    vtx2newix = fill(0, size(vtx_adjmtx, 1))
    @inbounds for (ix, id) in enumerate(newix2id)
        node = it.nodes[id]
        node.index = ix
        vtx2newix[node.vertices] .= ix
    end

    ithresh = it.ithreshold
    newcolptr = push!(sizehint!(Vector{NodeIx}(), length(newix2id) + 1), 1)
    newrowvals = sizehint!(Vector{NodeIx}(), length(rowvals(it.iadjmtx)))
    newweights = sizehint!(Vector{IWeight}(), length(nonzeros(it.iadjmtx)))

    # add column to the new adjmtx
    function appendcolumn(colw::AbstractVector{IWeight})
        @inbounds for (i, w) in enumerate(colw)
            if w >= ithresh
                push!(newrowvals, i)
                push!(newweights, w)
            end
        end
        push!(newcolptr, length(newrowvals) + 1)
    end

    tmpcolw = fill(0, length(newix2id)) # dense version of the newadjmtx column
    vtxw = nonzeros(vtx_adjmtx)
    vtxrowvals = rowvals(vtx_adjmtx)
    vtxorder = sortperm(vtx2newix) # order the vertices by the indices in newadjmtx
    @inbounds for i in vtxorder
        newix = vtx2newix[i]
        if newix == 0
            continue # skip unused vertices
        elseif newix > length(newcolptr) # the new newadjmtx column started, add the column for the previous one
            appendcolumn(tmpcolw)
            fill!(tmpcolw, 0) # reset the column
        else
            @assert newix == length(newcolptr)
        end
        # update the tmpcolw with the edges of the i-th vertex
        vtxnzrange = nzrange(vtx_adjmtx, i)
        for eix in vtxnzrange
            w = vtxw[eix]
            if w >= ithresh
                j = vtxrowvals[eix]
                newjx = vtx2newix[j]
                if newjx > 0
                    tmpcolw[newjx] = max(tmpcolw[newjx], w)
                end
            end
        end
    end
    # add the edges from the last newix to the new adjmtx
    if vtx2newix[vtxorder[end]] > 0
        appendcolumn(tmpcolw)
    end
    #@show newcolptr newix2id vtx2newix
    @assert length(newcolptr) == length(newix2id) + 1
    @assert all(>=(ithresh), newweights)
    it.iadjmtx = FlowAdjMatrix(length(newix2id), length(newix2id), newcolptr, newrowvals, newweights)
    it.nodeix2id = newix2id
    return it
end

# reset the iterator to the start where
# there's a single root node
function reset!(it::SCCTreeFlowsPeelingIterator)
    it.verbose && @info "Resetting SCCTreeFlowsPeelingIterator"
    empty!(it.nodes)
    empty!(it.flows)
    it.lastflowid = 0
    it.ithresholdpos = 0
    it.ithreshold = 0

    peeling = it.peeling
    if !isempty(peeling.sources)
        # initialize the single root node
        rootid = 1

        it.nodes[rootid] = FlowsPeelNode(rootid, 1, 1:nvertices(peeling.tree),
                                         sort!(collect(it.peeling.sources)), sort!(collect(it.peeling.sinks)), Set(1))
        it.nodeix2id = [rootid]
        it.iadjmtx = spdiagm(0 => [it.ithreshold])
        it.flows[it.lastflowid += 1] = SourceFlowsPeel(1, rootid, false, it.ithreshold,
                                                       spdiagm(0 => [it.ithreshold]), [FlowVertexInfo(rootid, !isempty(it.peeling.sinks))],
                                                       FlowSinkPath[])
    end
    return it
end

# advance to the next cutting threshold
# and update the adjmtx and source flows accordingly
function peel!(it::SCCTreeFlowsPeelingIterator; verbose::Bool = it.verbose)
    if it.ithresholdpos >= length(it.peeling.ithresholds)
        throw(BoundsError(it.peeling, it.ithresholdpos))
    end
    # select the next threshold
    it.ithresholdpos += 1
    it.ithreshold = it.peeling.ithresholds[it.ithresholdpos]
    verbose && @info "Threshold #$(it.ithresholdpos)/$(length(it)): w[$(thresholdindex(it))]=$(threshold(it))"

    expanded_nodes = collect_expanded_nodes!(it)
    #@show expanded_nodes
    expand_adjmatrix!(it, expanded_nodes)

    prevlastflowid = it.lastflowid
    expand_flows!(it, expanded_nodes, verbose=verbose)
    modified = prevlastflowid != it.lastflowid # flows modified
    for nodeid in keys(expanded_nodes)
        @assert isempty(it.nodes[nodeid].flows) "Expanded node $nodeid still has some flows attached"
    end
    for flow in values(it.flows)
        !flow.dirty && (flow.minweight >= it.ithreshold) && continue # skip stronger flows
        verbose && @info "Flow #$(flow.id): updating"
        prevpathshash = modified ? hash(0) : pathshash(flow, it.nodes)
        update!(flow, it, edgetest(it))
        modified = modified || (prevpathshash != pathshash(flow, it.nodes))
        if isempty(flow.paths) # disconnect empty flows from their nodes
            for vtxinfo in flow.vtxinfo
                node = get(it.nodes, vtxinfo.nodeid, nothing)
                if node !== nothing
                    verbose && @info "Node #$(node.id): detaching flow #$(flow.id) (empty)"
                    pop!(node.flows, flow.id)
                end
            end
        end
    end
    # remove empty flows
    filter!(pairs(it.flows)) do kv
        if isempty(kv[2].paths)
            verbose && @info "Deleting flow #$(kv[1]): empty paths"
            return false
        else
            return true
        end
    end
    # remove nodes not attached to any flows
    filter!(pairs(it.nodes)) do kv
        if isempty(kv[2].flows)
            verbose && @info "Deleting node #$(kv[1]): no flows"
            return false
        else
            return true
        end
    end
    verbose && @info "$(length(it.flows)) flow(s) of $(length(it.nodes)) node(s)"
    # check that all expanded nodes are removed
    for nodeid in keys(expanded_nodes)
        @assert !haskey(it.nodes, nodeid) "Expanded node $nodeid not removed"
    end
    return modified
end

function Base.iterate(it::SCCTreeFlowsPeelingIterator, ithreshpos::Integer = 0; verbose::Bool = it.verbose)
    if it.ithresholdpos != ithreshpos
        throw(ArgumentError("Iterator ithresholdpos ($(it.ithresholdpos)) and state ($(ithreshpos)) differ"))
    end
    # the last threshold is when each vertex is in its own SCC
    if ithreshpos >= length(it)
        return nothing
    end
    modified = true
    if it.ithreshold == 0
        reset!(it)
        if it.ithreshold == 0
            it.verbose && @info "Initial peeling"
            peel!(it, verbose=verbose)
        end
    else
        modified = peel!(it, verbose=verbose)
    end

    return SCCTreeFlowsPeels(it, modified), it.ithresholdpos
end

function flowstats(
        peels::SCCTreeFlowsPeels;
        sourcesinkweights::Union{AbstractMatrix, Nothing} = nothing,
        maxweight::Union{Number, Nothing} = nothing,
        used_sources::Union{AbstractVector{Int}, Nothing} = nothing,
        used_sinks::Union{AbstractVector{Int}, Nothing} = nothing,
        verbose::Bool = false
)
    nvtxflows = 0
    flowlen_sum = 0
    flowinvlen_sum = 0.0
    floweight_sum = 0.0
    flow_minedgeweight_sum = isempty(peels.weights) ? NaN : 0.0
    flowavghopweight_sum = 0.0
    compflowlen_sum = 0
    compflowinvlen_sum = 0.0
    compfloweight_sum = 0.0
    compflow_minedgeweight_sum = isempty(peels.weights) ? NaN : 0.0
    compflowlen_max = 0
    ncompflows = 0
    isnothing(used_sources) || empty!(used_sources)
    isnothing(used_sinks) || empty!(used_sinks)

    sinknodes = Set{NodeId}()
    for flow in values(peels.flows)
        verbose && @info "Processing flow #$(flow.id)"
        @assert !flow.dirty
        @assert !isempty(flow.paths)
        pathsources = peels.nodes[flow.srcid].sources
        ncompflows += length(flow.paths)
        if !isnothing(used_sources)
            append!(used_sources, pathsources)
        end
        for path in flow.paths
            pathsinks = peels.nodes[path.sinkid].sinks
            npairs = length(pathsources) * length(pathsinks)
            nvtxflows += npairs

            if !isnothing(sourcesinkweights)
                cur_floweight_sum = 0.0 # all pairwise flows between the vertices of current components
                @inbounds for src in pathsources
                    src_flows = view(sourcesinkweights, :, src)
                    for snk in pathsinks
                        floweight = src_flows[snk]
                        if !isnothing(maxweight) && (floweight > maxweight)
                            floweight = maxweight
                        end
                        cur_floweight_sum += floweight
                    end
                end
            else
                cur_floweight_sum = NaN
            end
            if !isnothing(used_sinks) && !(path.sinkid in sinknodes)
                append!(used_sinks, pathsinks)
            end
            push!(sinknodes, path.sinkid)
            floweight_sum += cur_floweight_sum
            compfloweight_sum += cur_floweight_sum / npairs

            flowlen_sum += npairs * path.len
            flowavghopweight_sum += cur_floweight_sum / (path.len + 1)
            flowinvlen_sum += npairs / (path.len + 1)

            compflowlen_sum += path.len
            compflowinvlen_sum += inv(path.len + 1)
            compflowlen_max = max(compflowlen_max, path.len)
            if !isempty(peels.weights)
                comp_minedgeweight = path.minweight > 0 ? peels.weights[path.minweight] : peels.weights[peels.ithreshold]
                flow_minedgeweight_sum += comp_minedgeweight * npairs
                compflow_minedgeweight_sum += comp_minedgeweight
            end
        end
    end

    return (nflows = nvtxflows, ncompflows=ncompflows,
            flowlen_sum = flowlen_sum, compflowlen_sum = compflowlen_sum,
            flowinvlen_sum = flowinvlen_sum, compflowinvlen_sum = compflowinvlen_sum,
            compflowlen_max = compflowlen_max,
            floweight_sum = floweight_sum, compfloweight_sum = compfloweight_sum,
            flowavghopweight_sum = flowavghopweight_sum,
            flow_minedgeweight_sum = flow_minedgeweight_sum,
            compflow_minedgeweight_sum = compflow_minedgeweight_sum,
            ncompsources = length(peels.flows),
            ncompsinks = length(sinknodes))
end
