# state of Depth-First-Search: vertex index + edge iterator state
const DFSState = Tuple{Int, Int}

strongly_connected_components(adjmtx::AbstractMatrix{T},
                              test::EdgeTest{T} = EdgeTest{T}(),
                              pools::Union{ObjectPools, Nothing} = nothing) where T =
    strongly_connected_components!(IndicesPartition(Vector{Int}(undef, size(adjmtx, 1)), Vector{Int}()),
                                   adjmtx, test, pools)

# Adapted/adopted from LightGraphs.jl
function strongly_connected_components!(components::IndicesPartition,
                                        adjmtx::AbstractMatrix{T},
                                        test::EdgeTest{T},
                                        pools::Union{ObjectPools, Nothing} = nothing) where T
    nnodes = size(adjmtx, 1)
    nnodes == size(adjmtx, 2) ||
        throw(DimensionMismatch("adjmtx has to be square, $(size(adjmtx)) found"))
    empty!(components)
    intpool = arraypool(pools, Int)
    index = fill!(borrow!(intpool, nnodes), 0)    # first time the vertex was visited, 0=unseen
    stack = borrow!(intpool, 0)   # stores vertices which have been discovered and not yet assigned to any component
    boolpool = arraypool(pools, Bool)
    onstack = fill!(borrow!(boolpool, nnodes), false)  # if a vertex is waiting in the stack to receive a component assignment
    lowlink = fill!(borrow!(intpool, nnodes), 0)  # lowest index vertex that it can reach through back edge (index array not vertex id number)
    parent = fill!(borrow!(intpool, nnodes), 0)   # parent of every vertex in dfs

    count = 1
    dfspool = arraypool(pools, DFSState)
    dfs_stack = borrow!(dfspool) # depth-first search stack
    @inbounds for s in axes(adjmtx, 2)
        (index[s] == 0) || continue # skip visited
        index[s] = count
        lowlink[s] = count
        onstack[s] = true
        parent[s] = s
        push!(stack, s)
        count += 1

        # start dfs from 's'
        push!(dfs_stack, (s, 0))
        while !isempty(dfs_stack)
            v, edgeitstate = dfs_stack[end] #end is the most recently added item
            u = 0 # index of the first unvisited neighbour
            @inbounds edgeit = outedges(adjmtx, v, test)
            @inbounds while true
                iter = iterate(edgeit, edgeitstate)
                isnothing(iter) && break
                (vout, _), edgeitstate = iter
                vout_index = index[vout]
                if vout_index == 0
                    # unvisited neighbor found
                    u = vout
                    break
                    #GOTO A push u onto DFS stack and continue DFS
                elseif onstack[vout] && (lowlink[v] > vout_index)
                    # we have already seen vout, but can update the lowlink of v
                    # which has the effect of possibly keeping v on the stack until n is ready to pop.
                    # update lowest index 'v' can reach through out neighbors
                    lowlink[v] = vout_index
                end
            end
            if u == 0
                # All out neighbors already visited or no out neighbors
                # we have fully explored the DFS tree from v.
                # time to start popping.
                popped, _ = pop!(dfs_stack)
                lowlink[parent[popped]] = min(lowlink[parent[popped]], lowlink[popped])

                if index[v] == lowlink[v]
                    # found a cycle in a completed dfs tree.
                    while !isempty(stack) #break when popped == v
                        # drain stack until we see v.
                        # everything on the stack until we see v is in the SCC rooted at v.
                        popped = pop!(stack)
                        pushelem!(components, popped)
                        onstack[popped] = false
                        # popped has been assigned a component, so we will never see it again.
                        if popped == v
                            # we have drained the stack of an entire component.
                            break
                        end
                    end
                    # finish the component
                    closepart!(components)
                end
            else #LABEL A
                # update edge iterator state in DFS stack
                dfs_stack[end] = v, edgeitstate
                # add unvisited neighbor to dfs
                index[u] = count
                lowlink[u] = count
                onstack[u] = true
                parent[u] = v
                count += 1
                push!(stack, u)
                push!(dfs_stack, (u, 0))
                # next iteration of while loop will expand the DFS tree from ui
            end
        end
    end
    @assert nelems(components) == nnodes
    release!(dfspool, dfs_stack)
    release!(intpool, parent)
    release!(intpool, lowlink)
    release!(intpool, stack)
    release!(intpool, index)
    release!(boolpool, onstack)
    return components
end
