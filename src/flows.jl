const Diedge = Pair{Int, Int}               # directed edge
const FlowInfoWIP = Int                     # in-progress FlowInfo
struct FlowInfo
    len::Int        # shortest path
    weight::Float64 # start => end walk weight
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
                            push!(flows, (v => dest, FlowInfo(vinfo, adjmtx[dest, v])))
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
                        adjmtx, comps, test, zerodiag=false)

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
    pools::Union{ObjectPools, Nothing} = nothing;
    mincompsize::Union{Integer, Nothing} = nothing
)
    # expand component flows into source/sink node flows
    empty!(flows)
    @inbounds for ((i, j), info) in compflows
        if !isnothing(mincompsize) # skip small components
            ((partlength(comps, i) < mincompsize) || (partlength(comps, j) < mincompsize)) && continue
        end
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
            if !isnothing(mincompsize) # skip small components
                ((partlength(comps, i) < mincompsize) || (partlength(comps, j) < mincompsize)) && continue
            end
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
    pools::Union{ObjectPools, Nothing} = nothing;
    mincompsize::Union{Integer, Nothing} = nothing
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
                                adjmtx, sources, sinks, test, pools;
                                mincompsize=mincompsize)

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
    pools::Union{ObjectPools, Nothing} = nothing;
    kwargs...
) = flowgraph!(subgraph, flows, cut!(empty!(comps), tree, test.threshold),
               adjmtx, sources, sinks, test, pools; kwargs...)

flowgraph(tree::SCCTree, adjmtx::AbstractMatrix,
          sources::AbstractVector{Int}, sinks::AbstractVector{Int},
          test::EdgeTest,
          pools::Union{ObjectPools, Nothing} = nothing;
          kwargs...
) = flowgraph!(Vector{CompDiedge}(), Vector{CompFlow}(), IndicesPartition(),
               tree, adjmtx, sources, sinks, test, pools; kwargs...)

"""
    nflows(comps::IndicesPartition, adjmtx::AbstractMatrix,
           sources::AbstractVector{Int}, sinks::AbstractVector{Int},
           test::EdgeTest;
           maxweight::Union{Number, Nothing} = nothing,
           used_sources::Union{AbstractVector{Int}, Nothing}=nothing,
           used_sinks::Union{AbstractVector{Int}, Nothing}=nothing) -> NamedTupe

Calculates statistics from the flows from *sources* to *sinks* vertices
in the weighted directed graph defined by *adjmtx* and *test* and the *comps* vertex components.

Returns the *NamedTuple* with the following fields
 * `nflows`: the number of *source* → *sink* flows
 * `ncompoflows`: the number of unique *component(source)* → *component(sink)* pairs for each *source* → *sink* flow
 * `flowlen_sum`: total length of all flows (flow length = the number of SCCs it crosses)
 * `compflowlen_sum`: total length of flows (every unique pairs of source/sink components is counted once)
 * `flowinvlen_sum`: the sum of ``1/\\mathrm{flowlength}``
 * `compflowinvlen_sum`: the sum of ``1/\\mathrm{flowlength}`` (each unique source/sink component pair is counted once)
 * `compflowlen_max`:
 * `floweight_sum`:
 * `floweight_sum`:
 * `compfloweight_sum`:
 * `flowavghopweight_sum`:
 * `ncompsources`: the number of distinct SCCs that have source nodes
 * `ncompsinks`: the number of distinct SCCs that have sink nodes

# Keyword arguments
* `maxweight`: if specified, limits the flow weight by `maxweight`
* `used_sources::AbstractVector{Int}`: if specified, acts as an output parameter that contains the
   sorted list of sources that have outgoing flows
* `used_sinks::AbstractVector{Int}`: if specified, acts as an output parameter that contains the
   sorted list of sinks that have incoming flows
"""
function nflows(
    comps::IndicesPartition,
    adjmtx::AbstractMatrix,
    sources::AbstractVector{Int}, sinks::AbstractVector{Int},
    test::EdgeTest,
    pools::Union{ObjectPools, Nothing} = nothing;
    maxweight::Union{Number, Nothing} = nothing,
    used_sources::Union{AbstractVector{Int}, Nothing}=nothing,
    used_sinks::Union{AbstractVector{Int}, Nothing}=nothing
)
    nvtxflows = 0
    flowlen_sum = 0
    flowinvlen_sum = 0.0
    floweight_sum = 0.0
    flowavghopweight_sum = 0.0
    compflowlen_sum = 0
    compflowinvlen_sum = 0
    compfloweight_sum = 0.0
    compflowlen_max = 0
    ncompflows = 0
    ncompsources = 0
    ncompsinks = 0
    isnothing(used_sources) || empty!(used_sources)
    isnothing(used_sinks) || empty!(used_sinks)
    if length(sources) > 0 || length(sinks) > 0 # || just to compute compsources/compsinks, for flows should be &&
        flowpool = arraypool(pools, Flow)
        compflows = borrow!(flowpool)
        ptnpool = objpool(pools, IndicesPartition)
        compsources = borrow!(ptnpool)
        compsinks = borrow!(ptnpool)
        componentsflowgraph!(nothing, compflows, compsources, compsinks, comps,
                             adjmtx, sources, sinks, test, pools)
        ncompsources = count(!isempty, compsources)
        ncompsinks = count(!isempty, compsinks)
        used_compsources = falses(length(compsources))
        used_compsinks = falses(length(compsinks))

        @inbounds for ((compi, compj), info) in compflows
            npairs = length(compsources[compi])*length(compsinks[compj])
            nvtxflows += npairs

            cur_floweight_sum = 0.0 # all pairwise flows between the vertices of current components
            for src in compsources[compi]
                src_flows = view(adjmtx, :, src)
                for snk in compsinks[compj]
                    floweight = src_flows[snk]
                    if !isnothing(maxweight) && (floweight > maxweight)
                        floweight = maxweight
                    end
                    cur_floweight_sum += floweight
                end
            end
            if !isnothing(used_sources) && !used_compsources[compi]
                used_compsources[compi] = true
                for src in compsources[compi]
                    isrc = searchsortedfirst(used_sources, src)
                    if (isrc > length(used_sources)) || (used_sources[isrc] != src)
                        insert!(used_sources, isrc, src)
                    end
                end
            end
            if !isnothing(used_sinks) && !used_compsinks[compj]
                used_compsinks[compj] = true
                for snk in compsinks[compj]
                    isnk = searchsortedfirst(used_sinks, snk)
                    if (isnk > length(used_sinks)) || (used_sinks[isnk] != snk)
                        insert!(used_sinks, isnk, snk)
                    end
                end
            end
            floweight_sum += cur_floweight_sum
            flowlen_sum += npairs * info.len
            flowavghopweight_sum += cur_floweight_sum / (info.len + 1)
            flowinvlen_sum += npairs / (info.len + 1)

            compflowlen_sum += info.len
            compflowinvlen_sum += inv(info.len + 1)
            compfloweight_sum += info.weight
            compflowlen_max = max(compflowlen_max, info.len)
        end
        ncompflows = length(compflows)
        release!(flowpool, compflows)
        release!(ptnpool, compsinks)
        release!(ptnpool, compsources)
    end

    return (nflows = nvtxflows, ncompflows=ncompflows,
            flowlen_sum = flowlen_sum, compflowlen_sum = compflowlen_sum,
            flowinvlen_sum = flowinvlen_sum, compflowinvlen_sum = compflowinvlen_sum,
            compflowlen_max = compflowlen_max,
            floweight_sum = floweight_sum, compfloweight_sum = compfloweight_sum,
            flowavghopweight_sum = flowavghopweight_sum,
            ncompsources = ncompsources,
            ncompsinks = ncompsinks)
end

function nflows(
    tree::SCCTree, adjmtx::AbstractMatrix,
    sources::AbstractVector{Int}, sinks::AbstractVector{Int},
    test::EdgeTest,
    pools::Union{ObjectPools, Nothing} = nothing;
    kwargs...
)
    ptnpool = objpool(pools, IndicesPartition)
    comps = cut!(empty!(borrow!(ptnpool)), tree, test.threshold)
    res = nflows(comps, adjmtx, sources, sinks, test, pools; kwargs...)
    release!(ptnpool, comps)
    return res
end

"""
    traceflows!(flow2paths::AbstractDict{Diedge, Partition{Int}},
                step_adjmtx::AbstractMatrix,
                steptest::EdgeTest,
                walk_adjmtx::AbstractMatrix,
                walktest::EdgeTest;
                sources::Union{AbstractVector, AbstractSet, Nothing} = nothing,
                sinks::Union{AbstractVector, AbstractSet, Nothing} = nothing,
                maxsteps::Integer=2) -> Dict{Diedge, Partition{Int}}

Trace the random walk (specified by *walk_adjmtx* and *walktest*) steps in the original graph
(given by *step_adjmtx* and *steptest*).

### Keyword arguments
* `sources` (optional): the indices of vertices to use as path starts.
  If not specified, all vertices are used as path starts.
* `sinks` (optional): the indices of vertices to use as path ends.
  If not specified, all vertices are used as path ends.
* `maxsteps=2`: maximal number of steps in a traced path, longer paths are discarded.

Returns the mapping from the flow diedges to the [`HierarchicalHotNet.Parition`](@ref)
object. Each part corresponds to the path, from diedge start to diedge end, in the original network,
and the part elements are the indices of the intermedidate vertices along the path (start and end not included).
"""
function traceflows!(
    flow2paths::AbstractDict{Diedge, Partition{Int}},
    step_adjmtx::AbstractMatrix,
    steptest::EdgeTest,
    walk_adjmtx::AbstractMatrix,
    walktest::EdgeTest,
    pools::Union{ObjectPools, Nothing} = nothing;
    sources::Union{AbstractVector, AbstractSet, Nothing} = nothing,
    sinks::Union{AbstractVector, AbstractSet, Nothing} = nothing,
    maxsteps::Integer=2
)
    nv = size(step_adjmtx, 1)
    nv == size(step_adjmtx, 2) ||
        throw(DimensionMismatch("Steps adjacency matrix must be square ($(size(step_adjmtx)) given)"))
    size(walk_adjmtx) == size(step_adjmtx) ||
        throw(DimensionMismatch("Steps and walk adjacency matrices must have equal sizes ($(size(step_adjmtx)) and $(size(walk_adjmtx)) given)"))
    (maxsteps > 0) || throw(ArgumentError("maxsteps must be positive ($(maxsteps) given)"))

    vtxsetpool = objpool(pools, Set{Int})
    visited = empty!(borrow!(vtxsetpool)) # local visited vertex
    dfspool = arraypool(pools, DFSState)
    dfs_stack = borrow!(dfspool) # local depth-first search stack
    @inbounds for v in axes(step_adjmtx, 2)
        isnothing(sources) || in(v, sources) || continue
        empty!(dfs_stack)
        push!(dfs_stack, (v, 0))
        empty!(visited)
        push!(visited, v)
        walkedges = view(walk_adjmtx, :, v)

        while !isempty(dfs_stack)
            w, edgeitstate = dfs_stack[end]
            @assert edgeitstate >= 0
            edgeit = outedges(step_adjmtx, w, steptest)

            # pick the edge to follow
            u = 0
            while true
                iter = iterate(edgeit, edgeitstate)
                isnothing(iter) && break
                (i, _), edgeitstate = iter
                if !(i in visited)
                    # update w edges iterator
                    dfs_stack[end] = w, edgeitstate
                    u = i
                    break
                end
            end
            if u > 0 # follow v->u Diedge
                # v -...-> u path is present in the walk graph, and u is a sink,
                # add its intermediate vertices to the result
                if isvalidedge(walkedges[u], walktest) && (isnothing(sinks) || in(u, sinks))
                    flowpaths = get!(() -> Partition{Int}(), flow2paths, v => u)
                    @inbounds for i in 2:length(dfs_stack)
                        pushelem!(flowpaths, dfs_stack[i][1])
                    end
                    closepart!(flowpaths) # finish given path
                end
                # if v ->..-> u path has not reached maximum length, continue DFS
                if length(dfs_stack) < maxsteps
                    push!(dfs_stack, (u, 0))
                    push!(visited, u)
                end
            else # backtrack from v
                pop!(dfs_stack)
            end
        end
    end
    # fill the stepgraph with traced steps
    release!(dfspool, dfs_stack)
    release!(vtxsetpool, visited)
    return flow2paths
end

"""
    traceflows(step_adjmtx::AbstractMatrix, steptest::EdgeTest,
               walk_adjmtx::AbstractMatrix, walktest::EdgeTest;
               kwargs...) -> Dict{Diedge, Partition{Int}}

Trace the random walk (specified by *walk_adjmtx* and *walktest*) steps in the original graph
(given by *step_adjmtx* and *steptest*).

See [`HierarchicalHotNet.traceflows!`](@ref) for the detailed description.
"""
traceflows(step_adjmtx::AbstractMatrix, steptest::EdgeTest,
           walk_adjmtx::AbstractMatrix, walktest::EdgeTest,
           pools::Union{ObjectPools, Nothing} = nothing; kwargs...) =
    traceflows!(Dict{Diedge, Partition{Int}}(),
                step_adjmtx, steptest, walk_adjmtx, walktest, pools; kwargs...)
