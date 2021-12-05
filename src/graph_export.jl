"""
    export_flowgraph(tree::SCCTree{T}, threshold::Number,
                     walkmatrix::AbstractMatrix,
                     sources::AbstractVector{<:Integer}, sinks::AbstractVector{<:Integer};
                     orig_diedges::Union{AbstractDataFrame, Nothing} = nothing,
                     vertices_stats::Union{AbstractDataFrame, Nothing} = nothing,
                     diedges_stats::Union{AbstractDataFrame, Nothing} = nothing,
                     flowpaths::Symbol = :skip,
                     stepmatrix::Union{AbstractMatrix, Nothing} = nothing,
                     step_threshold::Number = 0.75 * threshold, maxsteps::Integer = 2,
                     step_sinks::Union{AbstractVector, AbstractSet, Nothing} = nothing,
                     step_sources::Union{AbstractVector, AbstractSet, Nothing} = nothing,
                     flow_edges::Bool=false,
                     pvalue_mw_max::Number=0.05,
                     pvalue_fisher_max::Number=0.05,
                     verbose::Bool=false,
                     pools::Union{ObjectPools, Nothing}=nothing,
                     mincompsize::Union{Integer, Nothing}=nothing,
                     exported_sinks::AbstractVector{<:Integer}=sinks
    ) -> NamedTuple

Cuts the `tree` at `threshold` and exports the resulting SCC network
as the collection of dataframes.

If specified, calculates the flows from `sources` to `sinks` and returns
them as additional edges.

### Keyword arguments
  * `orig_diedges::AbstractDataFrame`: optional collection of the original edges.
    The metadata from this frame is added to the overlapping diedges of the SCC network.
  * `vertices_stats::AbstractDataFrame`: optional vertices statistics
  * `diedges_stats::AbstarctDataFrame`: optional directed edges statistics
  * `flowpaths::Symbol`: how the flows should be traced
     * `skip` (default): no tracing
     * `flowattr`: trace the flows (see [`HierarchicalHotNet.traceflows`](@ref)) and add as `flowpaths` column to *diedges* data frame
     * `steps`: trace the flows (see [`HierarchicalHotNet.traceflows`](@ref)) and add as extra diedges of type `step` to *diedges* data frame

### Returns

Named tuple with fields
  * `components::DataFrame`: statistics for *Strongly Connected Components*
  * `vertices::DataFrame`: network vertices
  * `diedges::DataFrame`: directed edges
  * `edges::DataFrame`: undirected edges

### See also
[`HierarchicalHotNet.cut`](@ref)
"""
function export_flowgraph(
    tree::SCCTree{T}, threshold::Number,
    walkmatrix::AbstractMatrix,
    sources::AbstractVector{<:Integer}, sinks::AbstractVector{<:Integer};
    orig_diedges::Union{AbstractDataFrame, Nothing} = nothing,
    vertices_stats::Union{AbstractDataFrame, Nothing} = nothing,
    diedges_stats::Union{AbstractDataFrame, Nothing} = nothing,
    flowpaths::Symbol = :skip,
    stepmatrix::Union{AbstractMatrix, Nothing} = nothing,
    step_threshold::Number = 0.75 * threshold, maxsteps::Integer = 2,
    step_sinks::Union{AbstractVector, AbstractSet, Nothing} = nothing,
    step_sources::Union{AbstractVector, AbstractSet, Nothing} = nothing,
    flow_edges::Bool=false,
    pvalue_mw_max::Number=0.05,
    pvalue_fisher_max::Number=0.05,
    verbose::Bool=false,
    pools::Union{ObjectPools, Nothing}=nothing,
    mincompsize::Union{Integer, Nothing}=nothing,
    exported_sinks::Union{AbstractVector{<:Integer}, Nothing}=sinks
) where T
    nvertices(tree) == size(walkmatrix, 1) == size(walkmatrix, 2) ||
        throw(DimensionMismatch("Number of tree vertices ($(nvertices(tree))) doesn't match the walk matrix dimensions ($(size(walkmatrix)))"))
    W = eltype(walkmatrix)
    verbose && isempty(sinks) && @warn "No sinks provided, flowgraph will be empty"
    subgraph, flows, conncomps = flowgraph(tree, walkmatrix, sources, sinks,
                                           EdgeTest{T}(threshold=threshold),
                                           pools, mincompsize=mincompsize)
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
    extra_vtxs = Vector{Int}()
    if flowpaths != :skip
        verbose && @info("Tracing flows paths (step_threshold=$(step_threshold), maxsteps=$maxsteps)...")
        isnothing(stepmatrix) && throw(ArgumentError("stepmatrix is required for tracing flows, but it was not provided"))

        (flowpaths ∈ [:stepedges, :flowattr]) ||
            throw(ArgumentError("Unsupported flowpaths=:$(flowpaths) mode. Only :skip, :stepedges and :flowattr are supported"))

        # trace the steps between the vertices of the flowgraph
        flow_walkmatrix = fill(zero(W), size(walkmatrix))
        @inbounds flow_walkmatrix[pos2vertex, pos2vertex] .= view(walkmatrix, pos2vertex, pos2vertex)
        flow2paths = traceflows(Matrix(stepmatrix), EdgeTest{W}(threshold=step_threshold),
                                flow_walkmatrix, EdgeTest{W}(threshold=threshold),
                                pools, maxsteps=maxsteps, sources=step_sources, sinks=step_sinks)
        verbose && @info("$(isempty(flow2paths) ? 0 : sum(ptn -> nelems(ptn) + nparts(ptn), values(flow2paths))) step edge(s) of $(length(flow2paths)) flow(s) traced")
        if flowpaths ==:stepedges
            # collect extra vertices and trace edges that are not yet in the flowgraph
            flowpath_steps = Set{Diedge}()
            extra_vtxs_set = Set{Int}()
            for ((src, trg), paths) in pairs(flow2paths)
                prev = src
                @inbounds for v in elems(paths)
                    push!(flowpath_steps, prev => v)
                    (vertex2pos[v] == 0) && push!(extra_vtxs_set, v)
                    prev = v
                end
                push!(flowpath_steps, prev => trg)
            end
            extra_vtxs = append!(extra_vtxs, extra_vtxs_set)
            verbose && @info("$(length(extra_vtxs)) extra vertice(s) for $(length(flowpath_steps)) step edges")
            flowpath_vertices = nothing
            flow2paths = nothing # don't use it later to create .flowpaths column
            # update the mapping
            vertex2pos[extra_vtxs] = length(pos2vertex) .+ eachindex(extra_vtxs)
            append!(pos2vertex, extra_vtxs)
        elseif flowpaths == :flowattr
            flowpath_steps = nothing
        else
            error("Unsupported flowpaths=:$(flowpaths) mode. Developers note: should have been checked earlier")
        end
    else
        flowpath_steps = nothing
        flow2paths = nothing
    end
    sort!(extra_vtxs)

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
    if !isnothing(exported_sinks)
        @inbounds for v in exported_sinks
            pos = vertex2pos[v]
            (pos > 0) && (vertices_df[pos, :is_sink] = true)
        end
    end
    vertices_df[!, :is_steptrace] .= false
    @inbounds for v in extra_vtxs
        pos = vertex2pos[v]
        (pos > 0) && (vertices_df[pos, :is_steptrace] = true)
    end
    if !isnothing(vertices_stats)
        if verbose
            missed_vertices_df = antijoin(vertices_df, vertices_stats, on=:vertex)
            (nrow(missed_vertices_df)>0) && @warn("vertices_stats missing information on $(nrow(missed_vertices_df)) of $(nrow(vertices_df)) vertice(s)")
        end
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
    if !isnothing(flow2paths)
        diedges_df.flowpaths = [get(flow2paths, r.source => r.target, missing) for r in eachrow(diedges_df)]
    end
    if !isnothing(diedges_stats)
        if verbose
            missed_diedges_df = antijoin(diedges_df, diedges_stats, on=[:source, :target])
            (nrow(missed_diedges_df)>0) && @warn("diedges_stats missing information on $(nrow(missed_diedges_df)) of $(nrow(diedges_df)) diedge(s)")
        end
        diedges_df = leftjoin(diedges_df, diedges_stats, on=[:source, :target])
    end

    flows_df.flow = Vector{String}()
    flows_df.flowlen = Vector{Int}()
    flows_df.floweight = Vector{W}()
    for ((src, trg), (srccomp, trgcomp), info) in flows
        (vertex2pos[trg] > 0 && vertices_df.is_sink[vertex2pos[trg]]) || continue
        (comp2new[srccomp] != 0) && (comp2new[trgcomp] != 0) || continue
        push!(flows_df, (source = src,
                         target = trg,
                         walkweight = threshold,
                         walkweight_rev = threshold,
                         flow = srccomp == trgcomp ? "loop" : "flow",
                         flowlen = info.len,
                         floweight = info.minweight))
    end

    if nrow(flows_df) > 0
        source_stats_df = combine(groupby(flows_df, :source)) do outedges_df
            sinks = sort!(unique(collect(zip(outedges_df.target, outedges_df.flowlen))),
                        by=x -> (x[2], x[1]))
            sourcesinks = sort!(unique!(outedges_df[outedges_df.flow .== "loop", :target]))
            DataFrame(flows_to = isempty(sinks) ? missing : [sinks],
                      nflows_to = length(sinks),
                      loops_through = isempty(sourcesinks) ? missing : [sourcesinks],
                      nloops_through = length(sourcesinks))
        end
        target_stats_df = combine(groupby(flows_df, :target)) do inedges_df
            sources = sort!(unique(collect(zip(inedges_df.source, inedges_df.flowlen))),
                            by=x -> (x[2], x[1]))
            DataFrame(flows_from = isempty(sources) ? missing : [sources],
                      nflows_from = length(sources))
        end
    else
        source_stats_df = DataFrame(source = Int[],
                                    flows_to = Vector{Int}[],
                                    nflows_to = Int[],
                                    loops_through = Vector{Int}[],
                                    nloops_through = Int[],)
        target_stats_df = DataFrame(target = Int[],
                                    flows_from = Vector{Int}[],
                                    nflows_from = Int[],)
    end

    if !isnothing(flowpath_steps)
        steps_df = DataFrame(source = first.(flowpath_steps),
                             target = last.(flowpath_steps))
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
            res[!, :has_original] .= !isnothing(orig1st)
            res[!, :has_original_rev] .= !isnothing(orig1st_rev)
        else
            orig1st = orig1st_rev = nothing
        end
        if hasproperty(edge_df, :diedge_type)
            res[!, :original_type] .= isnothing(orig1st) ? missing : edge_df.diedge_type[orig1st]
            res[!, :original_rev_type] .= isnothing(orig1st_rev) ? missing : edge_df.diedge_type[orig1st_rev]
        end
        if hasproperty(edge_df, :flowpaths)
            paths1st = findfirst(r -> !ismissing(r.flowpaths) && !r.is_reverse, eachrow(edge_df))
            paths1st_rev = findfirst(r -> !ismissing(r.flowpaths) && r.is_reverse, eachrow(edge_df))

            res[!, :flowpaths] .= isnothing(paths1st) ? missing : Ref(edge_df.flowpaths[paths1st])
            res[!, :flowpaths_rev] .= isnothing(paths1st_rev) ? missing : Ref(edge_df.flowpaths[paths1st_rev])
        end
        # include data of the original diedges
        if orig_diedges !== nothing
            for col in names(orig_diedges)
                (col == "source" || col == "target" || col == "diedge_type") && continue
                res[!, col] .= isnothing(orig1st) ? missing : edge_df[orig1st, col]
                res[!, col * "_rev"] .= isnothing(orig1st_rev) ? missing : edge_df[orig1st_rev, col]
            end
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
