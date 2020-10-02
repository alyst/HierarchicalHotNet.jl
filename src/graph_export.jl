function export_flowgraph(
    tree::SCCTree{T}, threshold::Number,
    walkmatrix::AbstractMatrix,
    sources::AbstractVector{<:Integer}, sinks::AbstractVector{<:Integer};
    orig_diedges::Union{AbstractDataFrame, Nothing} = nothing,
    vertices_stats::Union{AbstractDataFrame, Nothing} = nothing,
    diedges_stats::Union{AbstractDataFrame, Nothing} = nothing,
    stepmatrix::Union{AbstractMatrix, Nothing} = nothing,
    step_threshold::Number = 0.75 * threshold, maxsteps::Integer = 2,
    flow_edges::Bool=false,
    pvalue_mw_max::Number=0.05,
    pvalue_fisher_max::Number=0.05,
    verbose::Bool=false,
    pools::Union{ObjectPools, Nothing}=nothing
) where T
    nvertices(tree) == size(walkmatrix, 1) == size(walkmatrix, 2) ||
        throw(DimensionMismatch("Number of tree vertices ($(nvertices(tree))) doesn't match the walk matrix dimensions ($(size(walkmatrix)))"))
    W = eltype(walkmatrix)
    subgraph, flows, conncomps = flowgraph(tree, walkmatrix, sources, sinks,
                                           EdgeTest{T}(threshold=threshold),
                                           pools)
    components_df = conncomponents_stats(conncomps,
                                         average_weights=false,
                                         mannwhitney_tests=false)
    components_df[!, :is_used] .= false
    for ((_, _), (a, b)) in subgraph
        components_df.is_used[a] = true
        components_df.is_used[b] = true
    end
    used_comps = IndicesPartition()
    for (i, comp) in enumerate(conncomps)
        if components_df[i, :is_used]
            append!(used_comps.elems, comp)
            closepart!(used_comps)
        end
    end
    components_df[!, :component0] = components_df[!, :component]
    comp2new = fill(0, nrow(components_df))
    # FIXME what kind of filtering?, should it be before flowgraph()?
    #filter!(r -> r.is_used, components_df)
    components_df[!, :component] = 1:nrow(components_df)
    comp2new[components_df.component0] .= components_df.component

    verbose && @info "$(nelems(used_comps)) vertices in $(length(used_comps)) connected components"
    pos2vertex = copy(used_comps.elems)
    vertex2pos = fill(0, nvertices(tree))
    vertex2pos[pos2vertex] = eachindex(pos2vertex)
    if !isnothing(stepmatrix)
        # trace the steps between the vertices of the flowgraph
        flow_walkmatrix = fill(zero(W), size(walkmatrix))
        @inbounds flow_walkmatrix[pos2vertex, pos2vertex] .= view(walkmatrix, pos2vertex, pos2vertex)
        stepedges = tracesteps(Matrix(stepmatrix), EdgeTest{W}(threshold=step_threshold),
                               flow_walkmatrix, EdgeTest{W}(threshold=threshold),
                               pools, maxsteps=maxsteps)
        verbose && @info("$(length(stepedges)) step edge(s) traced (step_threshold=$(step_threshold), maxsteps=$maxsteps)")
        # add vertices from steps that are not yet in the flowgraph
        extra_vtxs_set = Set{Int}()
        @inbounds for (src, trg) in stepedges
            (vertex2pos[src] == 0) && push!(extra_vtxs_set, src)
            (vertex2pos[trg] == 0) && push!(extra_vtxs_set, trg)
        end
        extra_vtxs = sort!(collect(extra_vtxs_set))
        verbose && @info("$(length(extra_vtxs)) extra vertice(s) for step edges")
        # update the mapping
        vertex2pos[extra_vtxs] = length(pos2vertex) .+ eachindex(extra_vtxs)
        append!(pos2vertex, extra_vtxs)
    else
        stepedges = nothing
        extra_vtxs = Vector{Int}()
    end

    # create vertices dataframe
    vertices_df = DataFrame(vertex = pos2vertex)
    vertices_df[!, :component] .= 0
    for i in eachindex(used_comps)
        vertices_df[partrange(used_comps, i), :component] .= i
    end
    #@assert all(>(0), vertices_df.component) doesn't hold if there are intermediate step edges
    vertices_df[!, :is_source] .= false
    @inbounds for v in sources
        pos = vertex2pos[v]
        (pos > 0) && (vertices_df[pos, :is_source] = true)
    end
    vertices_df[!, :is_sink] .= false
    @inbounds for v in sinks
        pos = vertex2pos[v]
        (pos > 0) && (vertices_df[pos, :is_sink] = true)
    end
    vertices_df[!, :is_steptrace] .= false
    @inbounds for v in extra_vtxs
        pos = vertex2pos[v]
        (pos > 0) && (vertices_df[pos, :is_steptrace] = true)
    end
    if !isnothing(vertices_stats)
        vertices_df = leftjoin(vertices_df, vertices_stats, on=:vertex)
    end

    filter!(e -> comp2new[e[2][1]] > 0 && comp2new[e[2][2]] > 0, subgraph)
    filter!(e -> comp2new[e[2][1]] > 0 && comp2new[e[2][2]] > 0, flows)
    diedges_df = DataFrame(source = Vector{Int}(),
                           target = Vector{Int}(),
                           walkweight = Vector{Union{W, Missing}}(),
                           walkweight_rev = Vector{Union{W, Missing}}())
    flows_df = copy(diedges_df)
    for ((src, trg), (srccomp, trgcomp)) in subgraph
        (comp2new[srccomp] != 0) && (comp2new[trgcomp] != 0) || continue
        wj2i = walkmatrix[trg, src]
        @assert wj2i >= threshold
        push!(diedges_df, (source = src,
                           target = trg,
                           walkweight = wj2i,
                           walkweight_rev = walkmatrix[src, trg]))
    end
    if !isnothing(diedges_stats)
        diedges_df = leftjoin(diedges_df, diedges_stats, on=[:source, :target])
    end

    flows_df.flow = Vector{String}()
    flows_df.flowlen = Vector{Int}()
    flows_df.floweight = Vector{W}()
    for ((src, trg), (srccomp, trgcomp), info) in flows
        (comp2new[srccomp] != 0) && (comp2new[trgcomp] != 0) || continue
        push!(flows_df, (source = src,
                         target = trg,
                         walkweight = threshold,
                         walkweight_rev = threshold,
                         flow = srccomp == trgcomp ? "loop" : "flow",
                         flowlen = info.len,
                         floweight = info.weight))
    end

    source_stats_df = combine(groupby(flows_df, :source)) do outedges_df
        sinks = sort!(unique(collect(zip(outedges_df.flowlen,
                                         outedges_df.target))))
        sourcesinks = sort!(unique!(outedges_df[outedges_df.flow .== "loop", :target]))
        DataFrame(flows_to = isempty(sinks) ? missing :
                              join(string.(last.(sinks), '(', first.(sinks),')'), ' '),
                  nflows_to = length(sinks),
                  loops_through = isempty(sourcesinks) ? missing : join(sourcesinks, ' '),
                  nloops_through = length(sourcesinks))
    end
    target_stats_df = combine(groupby(flows_df, :target)) do inedges_df
        sources = sort!(unique(collect(zip(inedges_df.flowlen,
                                           inedges_df.source))))
        DataFrame(flows_from = isempty(sources) ? missing :
                               join(string.(last.(sources), '(', first.(sources),')'), ' '),
                  nflows_from = length(sources))
    end

    if !isnothing(stepedges)
        steps_df = DataFrame(source = first.(stepedges),
                             target = last.(stepedges))
        extra_diedges_df = antijoin(steps_df, diedges_df, on=[:source, :target])
        extra_diedges_df[!, :flow] .= "trace"
        append!(diedges_df, extra_diedges_df, cols=:union)
        verbose && @info("$(nrow(extra_diedges_df)) traced step edge(s) added")
    end
    append!(diedges_df, flows_df, cols=:union)

    if !isnothing(orig_diedges)
        diedges_df = leftjoin(diedges_df, orig_diedges, on=[:source, :target])
        diedges_df.is_original = .!ismissing.(coalesce.(diedges_df[!, :weight]))
    end
    used_vertices = union!(Set(diedges_df.source), Set(diedges_df.target))
    filter!(r -> r.vertex ∈ used_vertices, vertices_df)
    vertices_df = leftjoin(vertices_df, rename!(source_stats_df, :source => :vertex),
                           on=:vertex)
    vertices_df = leftjoin(vertices_df, rename!(target_stats_df, :target => :vertex),
                           on=:vertex)
    outedges_df = filter(r -> coalesce(r.walkweight, zero(W)) >= coalesce(r.walkweight_rev, zero(W)), diedges_df)
    outedges_df[!, :is_reverse] .= false
    inedges_df = filter(r -> coalesce(r.walkweight, zero(W)) < coalesce(r.walkweight_rev, zero(W)), diedges_df)
    inedges_df[!, :is_reverse] .= true
    rename!(inedges_df, :source => :target, :target => :source)

    function combine_diedges(edge_df::AbstractDataFrame)
        res = DataFrame(
            has_flow = !all(ismissing, edge_df.flowlen),
            has_trace = all(x -> coalesce(x, "") == "trace", edge_df.flow),
            has_walk = any(ismissing, edge_df.flowlen),
        )
        edge_mask = ismissing.(edge_df.flowlen) .& .!edge_df.is_reverse
        edgerev_mask = ismissing.(edge_df.flowlen) .& edge_df.is_reverse
        for col in [:walkweight, :prob_perm_walkweight_greater,
                    :walkpermweight_median, :walkpermweight_mad,
                    :walkpermweight_mean, :walkpermweight_std]
            hasproperty(edge_df, col) || continue
            res[!, col] .= any(edge_mask) ? maximum(edge_df[edge_mask, col]) : missing
            res[!, Symbol(col, "_rev")] .= any(edgerev_mask) ? maximum(edge_df[edgerev_mask, col]) : missing
        end
        if hasproperty(edge_df, :is_original)
            orig1st = findfirst(r -> r.is_original && !r.is_reverse, eachrow(edge_df))
            orig1st_rev = findfirst(r -> r.is_original && r.is_reverse, eachrow(edge_df))
            res.has_original = !isnothing(orig1st)
            res.has_original_rev = !isnothing(orig1st_rev)
        else
            orig1st = orig1st_rev = nothing
        end
        if hasproperty(edge_df, :diedge_type)
            res.target_type = isnothing(orig1st) ? missing : edge_df.diedge_type[orig1st]
            res.source_type = isnothing(orig1st_rev) ? missing : edge_df.diedge_type[orig1st_rev]
        end
        if hasproperty(edge_df, :interaction_type)
            res.interaction_type = isnothing(orig1st) ? missing : edge_df.interaction_type[orig1st]
        end
        if !flow_edges && !any(res.has_walk)
            return filter!(r -> false, res) # pure flows are not exported
        end
        return res
    end

    if nrow(outedges_df) + nrow(inedges_df) > 0
        edges_df = combine(combine_diedges, groupby(vcat(outedges_df, inedges_df),
                      [:source, :target]))
    else
        # workaround: add missing columns to the empty frame
        edges_df = combine_diedges(outedges_df)
        edges_df.source = Vector{Int}()
        edges_df.target = Vector{Int}()
    end
    return (components = components_df,
            vertices = vertices_df,
            diedges = diedges_df,
            edges = edges_df)
end
