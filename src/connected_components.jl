strongly_connected_components(adjmtx::AbstractMatrix,
                              pool::Union{ArrayPool{Int}, Nothing} = nothing;
                              kwargs...) =
    strongly_connected_components!(IndicesPartition(Vector{Int}(undef, size(adjmtx, 1)), Vector{Int}()),
                                   adjmtx, pool; kwargs...)

# Adapted/adopted from LightGraphs.jl
function strongly_connected_components!(components::IndicesPartition,
                                        adjmtx::AbstractMatrix{T},
                                        pool::Union{ArrayPool{Int}, Nothing} = nothing;
                                        skipval::Union{Number, Nothing}=zero(eltype(adjmtx)),
                                        threshold::Union{Number, Nothing}=nothing,
                                        rev::Bool=false) where T
    nnodes = size(adjmtx, 1)
    nnodes == size(adjmtx, 2) ||
        throw(DimensionMismatch("adjmtx has to be square, $(size(adjmtx)) found"))
    reset!(components, 0)
    index = fill!(borrow!(Int, pool, nnodes), 0)    # first time the vertex was visited, 0=unseen
    stack = borrow!(Int, pool, 0)   # stores vertices which have been discovered and not yet assigned to any component
    onstack = fill!(borrow!(Int, pool, nnodes), 0)  # 1 if a vertex is waiting in the stack to receive a component assignment
    lowlink = fill!(borrow!(Int, pool, nnodes), 0)  # lowest index vertex that it can reach through back edge (index array not vertex id number)
    parent = fill!(borrow!(Int, pool, nnodes), 0)   # parent of every vertex in dfs

    count = 1
    dfs_stack = borrow!(Int, pool, 0)   # depth-first search stack
    @inbounds for s in axes(adjmtx, 2)
        (index[s] == 0) || continue # skip visited
        index[s] = count
        lowlink[s] = count
        onstack[s] = 1
        parent[s] = s
        push!(stack, s)
        count += 1

        # start dfs from 's'
        push!(dfs_stack, s)

        while !isempty(dfs_stack)
            v = dfs_stack[end] #end is the most recently added item
            v_outedges = view(adjmtx, :, v)
            u = 0 # index of the first unvisited neighbour
            @inbounds for (vout, w) in enumerate(v_outedges)
                isvalidedge(w, skipval=skipval, threshold=threshold, rev=rev) || continue
                vout_index = index[vout]
                if vout_index == 0
                    # unvisited neighbor found
                    u = vout
                    break
                    #GOTO A push u onto DFS stack and continue DFS
                elseif (onstack[vout] == 1) && (lowlink[v] > vout_index)
                    # we have already seen n, but can update the lowlink of v
                    # which has the effect of possibly keeping v on the stack until n is ready to pop.
                    # update lowest index 'v' can reach through out neighbors
                    lowlink[v] = vout_index
                end
            end
            if u == 0
                # All out neighbors already visited or no out neighbors
                # we have fully explored the DFS tree from v.
                # time to start popping.
                popped = pop!(dfs_stack)
                lowlink[parent[popped]] = min(lowlink[parent[popped]], lowlink[popped])

                if index[v] == lowlink[v]
                    # found a cycle in a completed dfs tree.
                    while !isempty(stack) #break when popped == v
                        # drain stack until we see v.
                        # everything on the stack until we see v is in the SCC rooted at v.
                        popped = pop!(stack)
                        pushelem!(components, popped)
                        onstack[popped] = 0
                        # popped has been assigned a component, so we will never see it again.
                        if popped == v
                            # we have drained the stack of an entire component.
                            break
                        end
                    end
                    # finish the component
                    closepart!(components, sort=true)
                end
            else #LABEL A
                # add unvisited neighbor to dfs
                index[u] = count
                lowlink[u] = count
                onstack[u] = 1
                parent[u] = v
                count += 1
                push!(stack, u)
                push!(dfs_stack, u)
                # next iteration of while loop will expand the DFS tree from ui
            end
        end
    end
    @assert nelems(components) == nnodes
    release!(pool, dfs_stack)
    release!(pool, stack)
    release!(pool, onstack)
    release!(pool, index)
    release!(pool, lowlink)
    release!(pool, parent)
    return components
end
