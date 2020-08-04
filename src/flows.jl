const Diedge = Pair{Int, Int}               # directed edge
const FlowInfoWIP = Int                     # in-progress FlowInfo
struct FlowInfo
    len::Int        # shortest path
end
const Flow = Tuple{Diedge, FlowInfo}        # flow: source=>target, data
const CompDiedge = Tuple{Diedge, Diedge}    # directed edge from one component to the other (source comp => target comp, source vertex => target vertex)
const CompFlow = Tuple{Diedge, Diedge, FlowInfo} # directed flow from one component to the other

function componentsflowgraph!(
    subgraph::Union{AbstractVector{Diedge}, Nothing},
    flows::AbstractVector{Flow},
    adjmtx::AbstractMatrix{T},
    compsources::AbstractVector{<:AbstractVector},
    compsinks::AbstractVector{<:AbstractVector},
    test::EdgeTest{T} = EdgeTest{T}(),
    pools::Union{ObjectPools, Nothing} = nothing
) where T
    nv = size(adjmtx, 1)
    nv == size(adjmtx, 2) ||
        throw(DimensionMismatch("Adjacency matrix must be square ($(size(adjmtx)) given)"))
    nv == length(compsources) || throw(DimensionMismatch("compsources length ($(length(compsources))) doesn't match vertex count ($(size(adjmtx, 1)))"))
    nv == length(compsinks) || throw(DimensionMismatch("compsinks length ($(length(compsinks))) doesn't match vertex count ($(size(adjmtx, 1)))"))

    setpool = objpool(pools, Dict{Int, FlowInfoWIP})
    boolpool = arraypool(pools, Bool)
    visited = fill!(borrow!(boolpool, nv), false)
    reachablepool = arraypool(pools, Union{Dict{Int, FlowInfoWIP}, Nothing})
    reachable = fill!(borrow!(reachablepool, nv), nothing) # which sinks are reachable from this vertex
    dfspool = arraypool(pools, DFSState)
    dfs_stack = borrow!(dfspool)
    empty!(flows)
    !isnothing(subgraph) && empty!(subgraph)

    # updates flows through vertex v using the flows from u and edge u->v weight (uv)
    function updatereachable!(v::Integer, u::Integer)
        usinks = reachable[u]
        if isnothing(reachable[v])
            reachable[v] = vsinks = sizehint!(empty!(borrow!(setpool)), length(usinks))
            for (i, ui_len) in usinks
                vsinks[i] = ui_len+1
            end
        else
            vsinks = reachable[v]
            for (i, ui_len) in usinks
                vsinks[i] = min(get(vsinks, i, ui_len+1), ui_len+1)
            end
        end
        return reachable
    end

    @inbounds for v in axes(adjmtx, 2)
        (visited[v] || ispartempty(compsources, v)) && continue
        @assert isempty(dfs_stack)
        push!(dfs_stack, (v, 0))
        @assert isnothing(reachable[v])
        if !ispartempty(compsinks, v) # initialize reachable list of not-yet-visited vertices
            vsinks = empty!(borrow!(setpool))
            vsinks[v] = 0
            reachable[v] = vsinks
        end
        while !isempty(dfs_stack)
            v, edgeitstate = dfs_stack[end]
            @assert edgeitstate >= 0
            edgeit = outedges(adjmtx, v, test)
            visited[v] = true
            u = 0
            while true
                iter = iterate(edgeit, edgeitstate)
                isnothing(iter) && break
                (i, _), edgeitstate = iter
                if (i != v) && visited[i] &&
                   !isnothing(reachable[i]) && !isempty(reachable[i])
                    # since we require that flowgraph is DAG
                    # i should not be in the current dfs stack
                    # append flow to the reachable vertex to the used_edges list
                    updatereachable!(v, i)
                    !isnothing(subgraph) && push!(subgraph, (v => i))
                elseif !visited[i]
                    u = i
                    # update v edges iterator
                    dfs_stack[end] = v, edgeitstate
                    break
                end
            end
            if u > 0 # follow v->u Diedge
                push!(dfs_stack, (u, 0))
                @assert isnothing(reachable[u])
                if !ispartempty(compsinks, u)
                    usinks = empty!(borrow!(setpool))
                    usinks[u] = 0
                    reachable[u] = usinks
                end
            else # backtrack from v
                pop!(dfs_stack)
                if !ispartempty(compsources, v)
                    if !isnothing(reachable[v])
                        for (dest, vinfo) in reachable[v]
                            push!(flows, (v => dest, FlowInfo(vinfo)))
                        end
                        # self-loop
                        !ispartempty(compsinks, v) && !isnothing(subgraph) && push!(subgraph, v => v)
                    end
                end
                if !isempty(dfs_stack)
                    prev, _ = dfs_stack[end]
                    if !isnothing(reachable[v])
                        !isnothing(subgraph) && push!(subgraph, prev => v)
                        updatereachable!(prev, v)
                    end
                end
            end
        end
    end
    for set in reachable
        isnothing(set) || release!(setpool, set)
    end
    release!(reachablepool, reachable)
    release!(boolpool, visited)
    release!(dfspool, dfs_stack)
    return subgraph, flows
end

componentsflowgraph(adjmtx::AbstractMatrix{T},
        compsources::AbstractVector{<:AbstractVector},
        compsinks::AbstractVector{<:AbstractVector},
        test::EdgeTest{T} = EdgeTest{T}(),
        pools::Union{ObjectPools, Nothing} = nothing
) where T =
    componentsflowgraph!(Vector{Diedge}(), Vector{Flow}(),
            adjmtx, compsources, compsinks, test, pools)

function componentsflowgraph!(
    compgraph::Union{AbstractVector{Diedge}, Nothing},
    compflows::AbstractVector{Flow},
    compsources::IndicesPartition, compsinks::IndicesPartition,
    comps::IndicesPartition,
    adjmtx::AbstractMatrix,
    sources::AbstractVector{Int}, sinks::AbstractVector{Int},
    test::EdgeTest,
    pools::Union{ObjectPools, Nothing} = nothing
)
    # flow graph between connected components
    # using pool actually slows it down
    mtxpool = arraypool(pools, eltype(adjmtx))#NoopArrayPool{eltype(adjmtx)}()#
    compmtx = condense!(borrow!(mtxpool, (length(comps), length(comps))),
                        adjmtx, comps, test, zerodiag=true)

    # distribute sources and sinks into connected components
    ptnpool = objpool(pools, IndicesPartition)
    empty!(compsources)
    empty!(compsinks)
    @inbounds for comp in comps
        for i in comp
            isrc = searchsortedfirst(sources, i)
            if (isrc <= length(sources)) && (sources[isrc] == i)
                pushelem!(compsources, i)
            end
            isnk = searchsortedfirst(sinks, i)
            if (isnk <= length(sinks)) && (sinks[isnk] == i)
                pushelem!(compsinks, i)
            end
        end
        closepart!(compsources)
        closepart!(compsinks)
    end

    componentsflowgraph!(compgraph, compflows, compmtx,
                         compsources, compsinks, test, pools)
    release!(mtxpool, compmtx)

    return compgraph, compflows
end

# expand components flows and subgraph into vertex flows and subgraph
function expand_componentsflowgraph!(
    subgraph::Union{AbstractVector{CompDiedge}, Nothing},
    flows::AbstractVector{CompFlow},
    compgraph::Union{AbstractVector{Diedge}, Nothing},
    compflows::AbstractVector{Flow},
    compsources::IndicesPartition, compsinks::IndicesPartition,
    comps::IndicesPartition,
    adjmtx::AbstractMatrix,
    sources::AbstractVector{Int}, sinks::AbstractVector{Int},
    test::EdgeTest,
    pools::Union{ObjectPools, Nothing} = nothing
)
    # expand component flows into source/sink node flows
    empty!(flows)
    @inbounds for ((i, j), info) in compflows
        for src in compsources[i], snk in compsinks[j]
            push!(flows, (src => snk, i => j, info))
        end
    end

    # expand components subgraph into vertex subgraph
    if !isnothing(subgraph)
        empty!(subgraph)
        boolpool = arraypool(pools, Bool)
        used_comps = fill!(borrow!(boolpool, length(comps)), false)
        @inbounds for (i, j) in compgraph
            used_comps[i] = true
            used_comps[j] = true
            (i == j) && continue # skip components self-loops
            for src in comps[i]
                srcedges = view(adjmtx, :, src)
                for trg in comps[j]
                    if isvalidedge(srcedges[trg], test)
                        push!(subgraph, (src => trg, i => j))
                    end
                end
            end
        end

        # expand all used components
        # (both those that have self-loops in components graph and those that don't)
        @inbounds for (i, used) in enumerate(used_comps)
            used || continue
            for src in comps[i]
                srcedges = view(adjmtx, :, src)
                for trg in comps[i]
                    if isvalidedge(srcedges[trg], test)
                        push!(subgraph, (src => trg, i => i))
                    end
                end
            end
        end
        release!(boolpool, used_comps)
    end

    return subgraph, flows
end

function flowgraph!(
    subgraph::Union{AbstractVector{CompDiedge}, Nothing},
    flows::AbstractVector{CompFlow},
    comps::IndicesPartition,
    adjmtx::AbstractMatrix,
    sources::AbstractVector{Int}, sinks::AbstractVector{Int},
    test::EdgeTest,
    pools::Union{ObjectPools, Nothing} = nothing
)
    flowpool = arraypool(pools, Flow)
    compflows = borrow!(flowpool)
    edgepool = arraypool(pools, Diedge)
    compgraph = !isnothing(subgraph) ? borrow!(edgepool) : nothing

    ptnpool = objpool(pools, IndicesPartition)
    compsources = borrow!(ptnpool)
    compsinks = borrow!(ptnpool)

    componentsflowgraph!(compgraph, compflows, compsources, compsinks, comps,
                         adjmtx, sources, sinks, test, pools)
    expand_componentsflowgraph!(subgraph, flows, compgraph, compflows,
                                compsources, compsinks, comps,
                                adjmtx, sources, sinks, test, pools)

    release!(ptnpool, compsinks)
    release!(ptnpool, compsources)

    !isnothing(compgraph) && release!(edgepool, compgraph)
    release!(flowpool, compflows)

    return subgraph, flows, comps
end

flowgraph!(
    subgraph::Union{AbstractVector{CompDiedge}, Nothing},
    flows::AbstractVector{CompFlow},
    comps::IndicesPartition,
    tree::SCCTree, adjmtx::AbstractMatrix,
    sources::AbstractVector{Int}, sinks::AbstractVector{Int},
    test::EdgeTest,
    pools::Union{ObjectPools, Nothing} = nothing
) = flowgraph!(subgraph, flows, cut!(empty!(comps), tree, test.threshold),
               adjmtx, sources, sinks, test, pools)

flowgraph(tree::SCCTree, adjmtx::AbstractMatrix,
          sources::AbstractVector{Int}, sinks::AbstractVector{Int},
          test::EdgeTest,
          pools::Union{ObjectPools, Nothing} = nothing
) = flowgraph!(Vector{CompDiedge}(), Vector{CompFlow}(), IndicesPartition(),
               tree, adjmtx, sources, sinks, test, pools)

function nflows(
    comps::IndicesPartition,
    adjmtx::AbstractMatrix,
    sources::AbstractVector{Int}, sinks::AbstractVector{Int},
    test::EdgeTest,
    pools::Union{ObjectPools, Nothing} = nothing
)
    flowpool = arraypool(pools, Flow)
    compflows = borrow!(flowpool)

    ptnpool = objpool(pools, IndicesPartition)
    compsources = borrow!(ptnpool)
    compsinks = borrow!(ptnpool)
    componentsflowgraph!(nothing, compflows, compsources, compsinks, comps,
                         adjmtx, sources, sinks, test, pools)
    nvtxflows = 0
    flowlen_sum = 0
    compflowlen_sum = 0
    compflowlen_max = 0
    @inbounds for ((compi, compj), info) in compflows
        npairs = length(compsources[compi])*length(compsinks[compj])
        nvtxflows += npairs
        flowlen_sum += npairs * info.len
        compflowlen_sum += info.len
        compflowlen_max = max(compflowlen_max, info.len)
    end

    release!(flowpool, compflows)
    release!(ptnpool, compsinks)
    release!(ptnpool, compsources)

    return (nflows = nvtxflows, ncompflows=length(compflows),
            flowlen_sum = flowlen_sum, compflowlen_sum = compflowlen_sum,
            compflowlen_max = compflowlen_max,
            ncompsources = sum(!isempty, compsources),
            ncompsinks = sum(!isempty, compsinks))
end

function nflows(
    tree::SCCTree, adjmtx::AbstractMatrix,
    sources::AbstractVector{Int}, sinks::AbstractVector{Int},
    test::EdgeTest,
    pools::Union{ObjectPools, Nothing} = nothing
)
    ptnpool = objpool(pools, IndicesPartition)
    comps = cut!(empty!(borrow!(ptnpool)), tree, test.threshold)
    res = nflows(comps, adjmtx, sources, sinks, test, pools)
    release!(ptnpool, comps)
    return res
end
