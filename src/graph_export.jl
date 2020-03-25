function export_flowgraph(
    tree::SCCTree{T}, threshold::Number,
    walkmatrix::AbstractMatrix,
    sources::AbstractVector{<:Integer}, sinks::AbstractVector{<:Integer};
    orig_diedges::Union{AbstractDataFrame, Nothing} = nothing,
    vertices_stats::Union{AbstractDataFrame, Nothing} = nothing,
    flow_edges::Bool=false,
    pvalue_mw_max::Number=0.05,
    pvalue_fisher_max::Number=0.05,
    verbose::Bool=false,
    pools::Union{ObjectPools, Nothing}=nothing
) where T
    nvertices(tree) == size(walkmatrix, 1) == size(walkmatrix, 2) ||
        throw(DimensionMismatch("Number of tree vertices ($(nvertices(tree))) doesn't match the walk matrix dimensions ($(size(walkmatrix)))"))
    subgraph, flows, conncomps = flowgraph(tree, walkmatrix, sources, sinks,
                                           EdgeTest{T}(threshold=threshold),
                                           pools)
    components_df = conncomponents_stats(conncomps, vertices_stats,
                                         average_weights=true,
                                         mannwhitney_tests=true)
    components_df[!, :is_used] .= true
    used_comps = IndicesPartition()
    for (i, comp) in enumerate(conncomps)
        if (length(comp) >= 1 #= FIXME =#) &&
            (!hasproperty(components_df, :pvalue_walkweight_mw) ||
                (components_df.pvalue_walkweight_mw[i] <= pvalue_mw_max)) #&&
           #(comp_stats_df.pvalue_walkweights_mw[i] <= pvalue_fisher_max)
            append!(used_comps.elems, comp)
            closepart!(used_comps)
        else
            components_df[i, :is_used] = false
        end
    end
    components_df[!, :component0] = components_df[!, :component]
    comp2new = fill(0, nrow(components_df))
    filter!(r -> r.is_used, components_df)
    components_df[!, :component] = 1:nrow(components_df)
    comp2new[components_df.component0] .= components_df.component

    verbose && @info "$(nelems(used_comps)) vertices in $(length(used_comps)) connected components"
    vertex2pos = fill(0, nvertices(tree))
    vertex2pos[used_comps.elems] = 1:nelems(used_comps)
    vertices_df = DataFrame(vertex = used_comps.elems)
    vertices_df[!, :component] .= 0
    for i in eachindex(used_comps)
        vertices_df[partrange(used_comps, i), :component] .= i
    end
    @assert all(>(0), vertices_df.component)
    vertices_df[!, :is_source] .= false
    vertices_df[!, :is_sink] .= false
    @inbounds for src in sources
        pos = vertex2pos[src]
        (pos > 0) && (vertices_df[pos, :is_source] = true)
    end
    @inbounds for snk in sinks
        pos = vertex2pos[snk]
        (pos > 0) && (vertices_df[pos, :is_sink] = true)
    end
    if !isnothing(vertices_stats)
        vertices_df = join(vertices_df, vertex_info, on=:vertex, kind=:left)
    end

    filter!(e -> comp2new[e[2][1]] > 0 && comp2new[e[2][2]] > 0, subgraph)
    filter!(e -> comp2new[e[2][1]] > 0 && comp2new[e[2][2]] > 0, flows)
    diedges_df = DataFrame(source = Vector{Int}(),
                           target = Vector{Int}(),
                           walkweight = Vector{Float64}(),
                           walkweight_rev = Vector{Float64}())
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

    flows_df.flow = Vector{String}()
    for ((src, trg), (srccomp, trgcomp)) in flows
        (comp2new[srccomp] != 0) && (comp2new[trgcomp] != 0) || continue
        push!(flows_df, (source = src,
                         target = trg,
                         walkweight = threshold,
                         walkweight_rev = threshold,
                         flow = srccomp == trgcomp ? "loop" : "linear"))
    end
    source_stats_df = by(flows_df, :source) do outedges_df
        sinks = sort!(unique(outedges_df.target))
        sourcesinks = sort!(unique!(outedges_df[outedges_df.flow .== "circular", :target]))
        DataFrame(flows_to = isempty(sinks) ? missing : join(sinks, ' '),
                  nflows_to = length(sinks),
                  loops_through = isempty(sourcesinks) ? missing : join(sourcesinks, ' '),
                  nloops_through = length(sourcesinks))
    end
    target_stats_df = by(flows_df, :target) do inedges_df
        sources = sort!(unique(inedges_df.source))
        DataFrame(flows_from = isempty(sources) ? missing : join(sources, ' '),
                  nflows_from = length(sources))
    end
    diedges_df.flow = missings(String, nrow(diedges_df))
    append!(diedges_df, flows_df)
    if !isnothing(orig_diedges)
        diedges_df = join(diedges_df, orig_diedges, on=[:source, :target], kind=:left)
        diedges_df[!, :is_original] .= .!ismissing.(coalesce.(diedges_df[!, :weight]))
    end
    used_vertices = union!(Set(diedges_df.source), Set(diedges_df.target))
    filter!(r -> r.vertex ∈ used_vertices, vertices_df)
    vertices_df = join(vertices_df, rename!(source_stats_df, :source => :vertex),
                       on=:vertex, kind=:left)
    vertices_df = join(vertices_df, rename!(target_stats_df, :target => :vertex),
                       on=:vertex, kind=:left)
    outedges_df = filter(r -> r.walkweight >= r.walkweight_rev, diedges_df)
    outedges_df[!, :is_reverse] .= false
    inedges_df = filter(r -> r.walkweight < r.walkweight_rev, diedges_df)
    inedges_df[!, :is_reverse] .= true
    rename!(inedges_df, :source => :target, :target => :source)
    edges_df = by(vcat(outedges_df, inedges_df),
                  [:source, :target]) do edge_df
        res = DataFrame(
            has_flow = !all(ismissing, edge_df.flow),
            has_walk = any(ismissing, edge_df.flow),
            walkweight = any(r -> ismissing(r.flow) && !r.is_reverse, eachrow(edge_df)) ?
                    maximum(edge_df[ismissing.(edge_df.flow) .& .!edge_df.is_reverse, :walkweight]) : missing,
            walkweight_rev = any(r -> ismissing(r.flow) && r.is_reverse, eachrow(edge_df)) ?
                    maximum(edge_df[ismissing.(edge_df.flow) .& edge_df.is_reverse, :walkweight]) : missing
        )
        if hasproperty(edge_df, :is_original)
            res[!, :has_original] .= any(r -> r.is_original && !r.is_reverse, eachrow(edge_df))
            res[!, :has_original_rev] .= any(r -> r.is_original && r.is_reverse, eachrow(edge_df))
        end
        if hasproperty(edge_df, :diedge_type)
            orig1st = findfirst(r -> r.is_original && !r.is_reverse, eachrow(edge_df))
            res[!, :has_original] .= !isnothing(orig1st)
            res[!, :target_type] .= isnothing(orig1st) ? missing : edge_df.diedge_type[orig1st]
            if hasproperty(edge_df, :interaction_type)
                res[!, :interaction_type] .= isnothing(orig1st) ? missing : edge_df.interaction_type[orig1st]
            end
            orig1st_rev = findfirst(r -> r.is_original && r.is_reverse, eachrow(edge_df))
            res[!, :has_original_rev] .= !isnothing(orig1st_rev)
            res[!, :source_type] .= isnothing(orig1st_rev) ? missing : edge_df.diedge_type[orig1st_rev]
        end
        if !flow_edges && !res.has_walk[1]
            return filter!(r -> false, res) # pure flows are not exported
        end
        return res
    end
    return (components = components_df,
            vertices = vertices_df,
            diedges = diedges_df,
            edges = edges_df)
end
