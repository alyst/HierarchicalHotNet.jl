"""
Adjacency matrix of a modified directed weighted graph.
For each of ``N_{te}`` *tunnel entry* vertices it creates a new mirror vertex
that has the same incoming edges as the original one and outgoing edges
(*tunnels*) that connect it with the designated mirror *tunnel exit* nodes.
There's a single unique *tunnel exit* mirror mode for each tunnel. It has the
single incoming edge -- from the mirror *tunnel entry* node, while the outgoing
the edges are the same as for the original *exit* vertex (except self-loops).

If there's a directed path from the original *tunel exit* to *tunel entry* node,
it transforms into a directed cycle going throw the corresponding *tunnel edge*
between the mirror nodes. This directed cycle, therefore, is a part of some
strongly connected component.
"""
mutable struct TunnelsMatrix{W, M <: AbstractMatrix} <: AbstractMatrix{W}
    parent::M               # original adjacency matrix
    nparent::Int            # number of nodes (1:nparent) from the original graph exposed
    entries::Vector{Int}    # indices of parent "entry" rows
    tunnels::IndicesPartition # one part for each entry with indices of the parent "exit" columns
    tunnel_weight::W        # weight of the tunnel edge connecting entry with exit

    function TunnelsMatrix(parent::AbstractMatrix{W},
        entries::AbstractVector{Int}, tunnels::AbstractPartition{Int},
        nparent::Integer = min(size(parent, 1), size(parent, 2));
        tunnel_weight::Number = !isempty(parent) ? maximum(parent) : zero(W)
    ) where W
        minside = min(size(parent, 1), size(parent, 2))
        (0 <= nparent <= minside) ||
            throw(ArgumentError("The number of exposed parent nodes should be in 0:$(minside) range, $nparent given"))
        length(entries) == length(tunnels) ||
            throw(DimensionMismatch("The number of entry points ($(length(entries))) and the number of tunnel parts ($(length(tunnels))) does not match"))
        new{W, typeof(parent)}(parent, nparent, entries, tunnels, tunnel_weight)
    end

    # constructor for tunnels between all entries and exits
    TunnelsMatrix(parent::AbstractMatrix,
                  entries::AbstractVector{Int}, exits::AbstractVector{Int};
                  kwargs...) =
        TunnelsMatrix(parent, entries,
                      repeat!(IndicesPartition(exits), length(entries)),
                      min(size(parent, 1), size(parent, 2)); kwargs...)
end

parenttype(::Type{TunnelsMatrix{W, M}}) where {W, M} = M
parenttype(mtx::TunnelsMatrix) = parenttype(typeof(mtx))

nentries(mtx::TunnelsMatrix) = length(mtx.entries)
ntunnels(mtx::TunnelsMatrix) = nelems(mtx.tunnels)

@inline Base.size(mtx::TunnelsMatrix, i::Integer) =
    mtx.nparent + nentries(mtx) + ntunnels(mtx)

@inline Base.size(mtx::TunnelsMatrix) = (size(mtx, 1), size(mtx, 2))

Base.@propagate_inbounds function Base.getindex(mtx::TunnelsMatrix,
                                                i::Integer, j::Integer)
    s1 = size(mtx.parent, 1)
    s2 = size(mtx.parent, 2)
    if j <= mtx.nparent # outgoing edges of the original node
        if i <= mtx.nparent # parent matrix
            return getindex(mtx.parent, i, j)
        elseif i <= mtx.nparent + nentries(mtx) # j -> entry(i)=ie
            # entries preserve all incoming edges of their original nodes, except self-loops
            ie = mtx.entries[i - mtx.nparent]
            return j != ie ? getindex(mtx.parent, ie, j) : zero(eltype(mtx))
        else # original nodes cannot enter tunnels
            return zero(eltype(mtx))
        end
    elseif j <= mtx.nparent + nentries(mtx) # outgoing edges of tunnel entries
        if i <= mtx.nparent + nentries(mtx)
            # tunnel entry can only lead into tunnel exit
            return zero(eltype(mtx))
        else
            jem = j - mtx.nparent
            ixm = i - (mtx.nparent + nentries(mtx))
            # there are tunnels for js-th entry to each exit,
            # they have indices ((js-1)*nf) + 1:nf
            return in(ixm, partrange(mtx.tunnels, jem)) ?
                mtx.tunnel_weight : zero(eltype(mtx))
        end
    else # outgoing edges of tunnel exits
        if i <= mtx.nparent
            # exits preserve the outgoing edges of their original nodes, except self-loops
            jxm = j - (mtx.nparent + nentries(mtx))
            @inbounds jx = mtx.tunnels.elems[jxm]
            return i != jx ? getindex(mtx.parent, i, jx) : zero(eltype(mtx))
        else
            return zero(eltype(mtx))
        end
    end
end

struct TunnelsMatrixOutedgesIterator{W <: Number, V <: AbstractVector, M <: AbstractMatrix, T} <:
            AbstractOutedgesIterator{W}
    col::V
    threshold::T
    rev::Bool

    mtx::TunnelsMatrix{W, M}
    colv::Int
    exitv::Int
    entryv::Int

    function TunnelsMatrixOutedgesIterator(mtx::TunnelsMatrix, v::Integer;
        threshold::Union{Number, Nothing}=nothing,
        rev::Bool=false
    )
        W = eltype(mtx)
        T = isnothing(threshold) ? Nothing : W

        checkindex(Bool, axes(mtx, 2), v) ||
            throw(BoundsError("attempt to access column $v of $(size(mtx.parent, 2))-column TunnelsMatrix"))
        lastentry = mtx.nparent + nentries(mtx)
        if v <= mtx.nparent # original vertex
            colv = v
            exitv = 0
            entryv = 0
        elseif v <= lastentry # entry mirror vertex
            colv = 1 # no original outgoing edges for entry, so it doesn't matter
            exitv = 0
            entryv = v - mtx.nparent
        else # exit mirror vertex
            # original exit vertex
            exitv = colv = mtx.tunnels.elems[v - mtx.nparent - nentries(mtx)]
            entryv = 0
        end
        col = view(mtx.parent, :, colv)
        new{W, typeof(col), typeof(mtx.parent), T}(
            col, threshold, rev,
            mtx, colv, exitv, entryv)
    end
end

function outedges(mtx::TunnelsMatrix, v::Integer;
                  skipval::Union{Number, Nothing} = zero(eltype(mtx)), kwargs...)
    skipval != zero(eltype(mtx)) && throw(ArgumentError("TunnelsMatrix outedges iterator only supports zero skipvals"))
    return TunnelsMatrixOutedgesIterator(mtx, v; kwargs...)
end

skipval(it::TunnelsMatrixOutedgesIterator{W}) where W = zero(W)
isweightreverse(it::TunnelsMatrixOutedgesIterator) = it.rev

function Base.iterate(it::TunnelsMatrixOutedgesIterator, i::Integer = 0)
    lastentry = it.mtx.nparent + nentries(it.mtx)
    if it.entryv > 0 # v is entry, so it only connects to tunnel exits
        # assuming that tunnel_weight is always above threshold
        if i == 0
            i = lastentry + @inbounds(it.mtx.tunnels.starts[it.entryv])
        else
            i += 1
        end
        if i < lastentry + @inbounds(it.mtx.tunnels.starts[it.entryv+1])
            return (i => it.mtx.tunnel_weight), i
        else
            return nothing
        end
    end

    # exit or the original graph vertex
    while i < it.mtx.nparent
        i += 1
        (i == it.exitv) && continue
        w = @inbounds(it.col[i])
        isvalidedge(w, it) && return (i => w, i)
    end
    (it.exitv > 0) && return nothing # it's exit vertex, so no further edges

    # outgoing edges from the original vertices into entry mirrors
    l = it.mtx.nparent
    while i < lastentry
        i += 1
        ie = it.mtx.entries[i - l]
        (ie == it.colv) && continue # no outgoing edge to own mirror
        w = @inbounds(it.col[ie])
        isvalidedge(w, it) && return (i => w, i)
    end
    return nothing
end

function indexvalues!(iA::TunnelsMatrix, weights::AbstractVector,
                      A::TunnelsMatrix; kwargs...)
    iparent, weights = indexvalues!(iA.parent, weights, A.parent; kwargs...)
    @assert last(weights) == A.tunnel_weight
    iA.parent = iparent
    iA.nparent = A.nparent
    iA.tunnel_weight = length(weights)
    copy!(iA.entries, A.entries)
    copy!(iA.tunnels, A.tunnels)
    return iA, weights
end

indexvalues(::Type{I}, A::TunnelsMatrix{T}; kwargs...) where {T, I <: Integer} =
    indexvalues!(TunnelsMatrix(Matrix{I}(undef, (0, 0)), Vector{Int}(),
                               IndicesPartition()),
                 Vector{T}(), A; kwargs...)

function subgraph_adjacencymatrix(adjmtx::TunnelsMatrix,
                                  comp_indices::AbstractVector{<:Integer})
    parent_cols = parent_rows = Vector{Int}()
    newentries = Vector{Int}()
    newexits = Vector{Int}()
    newtunnelstarts = Vector{Int}()
    nnewparents = 0
    lastentry = adjmtx.nparent + nentries(adjmtx)
    old2new_iem = fill(0, nentries(adjmtx))
    old2new_ix = fill(0, maximum(adjmtx.tunnels.elems))
    last_new_iem = 0
    for i in comp_indices
        if i <= adjmtx.nparent
            @assert nnewparents == 0
            push!(parent_cols, i)
        elseif i <= lastentry
            if nnewparents == 0
                # mark the bound of visible parents and remember rows
                nnewparents = length(parent_cols)
                parent_rows = copy(parent_cols)
            end
            old_iem = i - adjmtx.nparent
            old_ie = adjmtx.entries[old_iem]
            new_ie = searchsortedfirst(view(parent_rows, 1:nnewparents), old_ie)
            if (new_ie > length(parent_rows)) || (parent_rows[new_ie] != old_ie)
                push!(parent_rows, old_ie)
                new_ie = length(parent_rows)
            end
            push!(newentries, new_ie)
            old2new_iem[old_iem] = length(newentries) # remember
        else # tunnels exits
            ixm = i - lastentry
            # find which new entry it belongs to:
            # find its old entry
            old_iem = searchsortedlast(adjmtx.tunnels.starts, ixm)
            @assert 1 <= old_iem <= nentries(adjmtx)
            # find the matching new entry
            new_iem = old2new_iem[old_iem]
            (new_iem > 0) || throw(ArgumentError("The view doesn't keep the #$(old_iem) entry mirror vertex, but keeps its exit #$i"))
            # fix the number of tunnels of the all the entries before
            @assert new_iem >= last_new_iem
            while last_new_iem < new_iem
                last_new_iem += 1
                push!(newtunnelstarts, length(newexits)+1)
            end
            # find the new exit
            old_ix = adjmtx.tunnels.elems[ixm]
            new_ix = old2new_ix[old_ix]
            if new_ix == 0 # not yet encountered
                # search within old columns
                new_ix = searchsortedfirst(view(parent_cols, 1:nnewparents), old_ix)
                if (new_ix > length(parent_cols)) || (parent_cols[new_ix] != old_ix)
                    push!(parent_cols, old_ix)
                    new_ix = length(parent_cols)
                end
                old2new_ix[old_ix] = new_ix # remember
            end
            push!(newexits, new_ix)
        end
    end
    if nnewparents == 0 # no nodes beyond parents
        # mark the bound of visible parents and remember rows
        nnewparents = length(parent_cols)
        # note that parent_rows point to parent_cols
    end
    # close the last tunnel part(s)
    while length(newtunnelstarts) <= length(newentries)
        push!(newtunnelstarts, length(newexits) + 1)
    end
    parentmtx = view(adjmtx.parent, parent_rows, parent_cols)
    return TunnelsMatrix(parentmtx, newentries,
                         IndicesPartition(newexits, newtunnelstarts),
                         nnewparents, tunnel_weight=adjmtx.tunnel_weight)
end

function condense!(B::AbstractMatrix,
                   A::TunnelsMatrix,
                   node_groups::AbstractPartition;
                   skipval::Union{Number, Nothing} = zero(eltype(A)),
                   rev::Bool=false)
    nnodes = nelems(node_groups)
    size(A) == (nnodes, nnodes) ||
        throw(DimensionMismatch("A size ($(size(A))) and row/col labels sizes ($nnodes, $nnodes) do not match."))
    size(B) == (length(node_groups), length(node_groups)) ||
        throw(DimensionMismatch("B size ($(size(B))) and row/col group number ($(length(node_groups)), $(length(node_groups))) do not match."))
    fill!(B, defaultweight(eltype(A), skipval=skipval, rev=rev))

    lastentry = A.nparent + nentries(A)
    @inbounds for (jj, cols) in enumerate(node_groups)
        B_j = view(B, :, jj)
        for j in cols
            if j <= A.nparent # original vertex
                colj = j
                exitj = 0
                entryj = 0
            elseif j <= lastentry # entry mirror vertex
                colj = 1 # no original outgoing edges for entry, so it doesn't matter
                exitj = 0
                entryj = j - A.nparent
            else # exit mirror vertex
                # original exit vertex
                exitj = colj = A.tunnels.elems[j - lastentry]
                entryj = 0
            end

            Aj = view(A.parent, :, colj)
            for (ii, rows) in enumerate(node_groups)
                (B_j[ii] == A.tunnel_weight) && continue # already at maximum
                lastentry = A.nparent + nentries(A)
                if entryj > 0 # j is entry, so it only connects to tunnel exits
                    starti = lastentry + A.tunnels.starts[entryj]
                    ix = searchsortedfirst(rows, starti)
                    if ix <= length(rows) # some indices after the first tunnel or match it
                        endi = lastentry + A.tunnels.starts[entryj+1]-1
                        if starti <= rows[ix] <= endi # there's at least one tunnel
                            B_j[ii] = A.tunnel_weight
                        end
                    end
                else
                    Bij = B_j[ii]
                    for i in rows
                        if i <= A.nparent
                            (i == exitj) && continue
                            w = Aj[i]
                        elseif (exitj != 0)
                            break # it's an exit vertex, no further outgoing edges
                        elseif i <= lastentry
                            # check for outgoing edges to the entries
                            ie = A.entries[i - A.nparent]
                            (ie == j) && continue
                            w = Aj[ie]
                        else
                            break # no further non-zero entries
                        end
                        !isnothing(skipval) && (w == skipval) && continue
                        if (!isnothing(skipval) && (Bij == skipval)) || isweaker(Bij, w, rev=rev)
                            Bij = w
                        end
                    end
                    B_j[ii] = Bij
                end
            end
        end
    end
    return B
end
