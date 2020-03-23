function conncomponents_graph(
    tree::SCCTree, threshold::Number,
    walkmatrix::AbstractMatrix;
    orig_diedges::Union{AbstractDataFrame, Nothing} = nothing,
    vertices_stats::Union{AbstractDataFrame, Nothing} = nothing,
    minsize::Integer=2,
    pvalue_mw_max::Number=0.05,
    pvalue_fisher_max::Number=0.05,
    verbose::Bool=false
)
    nvertices(tree) == size(walkmatrix, 1) == size(walkmatrix, 2) ||
        throw(DimensionMismatch("Number of tree vertices ($(nvertices(tree))) doesn't match the walk matrix dimensions ($(size(walkmatrix)))"))
    conncomps = cut(tree, threshold, minsize=1)
    components_df = conncomponents_stats(conncomps, vertices_stats,
                                         average_weights=true,
                                         mannwhitney_tests=true)
    components_df[!, :is_used] .= true
    used_comps = IndicesPartition()
    for (i, comp) in enumerate(conncomps)
        if (length(comp) >= minsize) &&
            (!hasproperty(components_df, :pvalue_walkweight_mw) ||
                (components_df.pvalue_walkweight_mw[i] <= pvalue_mw_max)) #&&
           #(comp_stats_df.pvalue_walkweights_mw[i] <= pvalue_fisher_max)
            append!(used_comps.elems, comp)
            closepart!(used_comps)
        else
            components_df[i, :is_used] = false
        end
    end
    filter!(r -> r.is_used, components_df)
    components_df[!, :component] = 1:nrow(components_df)
    verbose && @info "$(nelems(used_comps)) vertices in $(length(used_comps)) connected components"
    vertices_df = DataFrame(vertex = used_comps.elems)
    vertices_df[!, :component] .= 0
    for i in eachindex(used_comps)
        vertices_df[partrange(used_comps, i), :component] .= i
    end
    @assert all(>(0), vertices_df.component)
    if !isnothing(vertices_stats)
        vertices_df = join(vertices_df, vertex_info, on=:vertex, kind=:left)
    end
    vertices_df[!, :covertex] .= vertices_df[!, :vertex]
    if walkmatrix isa TunnelsMatrix
        vertices_df[!, :is_entry_ref] .= false
        vertices_df[!, :is_exit_ref] .= false
        vertices_df[!, :is_entry] .= false
        vertices_df[!, :is_exit] .= false
        lastentry = walkmatrix.nparent + nentries(walkmatrix)
        for (i, v) in enumerate(vertices_df.vertex)
            if v > walkmatrix.nparent
                if v <= lastentry
                    vref = walkmatrix.entries[v - walkmatrix.nparent]
                    vertices_df[i, :is_entry] = true
                    vrefpos = searchsortedfirst(vertices_df.vertex, vref)
                    if (vrefpos <= nrow(vertices_df)) &&
                       (vertices_df.vertex[vrefpos] == vref)
                        vertices_df[vrefpos, :is_entry_ref] = true
                    end
                else
                    vref = walkmatrix.tunnels.elems[v - lastentry]
                    vertices_df[i, :is_exit] = true
                    vrefpos = searchsortedfirst(vertices_df.vertex, vref)
                    if (vrefpos <= nrow(vertices_df)) &&
                       (vertices_df.vertex[vrefpos] == vref)
                        vertices_df[vrefpos, :is_exit_ref] = true
                    end
                end
                vertices_df[i, :covertex] = vref
            end
        end
    end
    # merge together gene vertices and their sinks/sources mirrors
    covertices_df = by(vertices_df, [:component, :covertex]) do covertex_df
        vertices = sort!(covertex_df.vertex)
        res = DataFrame(
            vertices = join(vertices, ' '),
            nvertices = length(vertices),
            is_entry = any(covertex_df.is_entry),
            is_exit = any(covertex_df.is_exit),
            is_entry_ref = any(covertex_df.is_entry_ref),
            is_exit_ref = any(covertex_df.is_exit_ref),
            n_entries = sum(covertex_df.is_entry),
            n_exits = sum(covertex_df.is_exit)
        )
        if hasproperty(covertex_df, :weight)
            res.weight = maximum(covertex_df.weight)
        end
        if hasproperty(covertex_df, :walkweight)
            res.walkweight = maximum(covertex_df.walkweight)
        end
        return res
    end

    edges_matrix = walkmatrix[used_comps.elems, used_comps.elems]
    diedges_df = DataFrame(source = Vector{Int}(),
                           target = Vector{Int}(),
                           component = Vector{Int}(),
                           is_tunnel = Vector{Bool}(),
                           walkweight = Vector{Float64}(),
                           walkweight_rev = Vector{Float64}())
    for i in axes(edges_matrix, 1), j in axes(edges_matrix, 2)
        (vertices_df.component[i] == vertices_df.component[j]) || continue
        wj2i = edges_matrix[i, j]
        (wj2i >= threshold) || continue
        push!(diedges_df, (source = used_comps.elems[j],
                           target = used_comps.elems[i],
                           component = vertices_df.component[j],
                           is_tunnel = vertices_df.is_entry[j] && vertices_df.is_exit[i],
                           walkweight = wj2i,
                           walkweight_rev = edges_matrix[j, i]))
    end
    if !isnothing(orig_diedges)
        diedges_df = join(diedges_df, orig_diedges, on=[:source, :target], kind=:left)
        diedges_df[!, :is_original] .= .!ismissing.(coalesce.(diedges_df[!, :weight]))
    end
    diedges_df = join(diedges_df, rename!(vertices_df[!, [:vertex, :covertex]],
                                          :vertex=>:source, :covertex=>:cosource),
                      on=:source, kind=:left)
    diedges_df = join(diedges_df, rename!(vertices_df[!, [:vertex, :covertex]],
                                          :vertex=>:target, :covertex=>:cotarget),
                      on=:target, kind=:left)
    cosource_stats_df = by(diedges_df, :cosource) do codiedges_df
        cotargets = sort!(unique!(codiedges_df[codiedges_df.is_tunnel, :cotarget]))
        DataFrame(tunnels_to = isempty(cotargets) ? missing : join(cotargets, ' '),
                  ntunnels_to = length(cotargets))
    end
    cotarget_stats_df = by(diedges_df, :cotarget) do codiedges_df
        cosources = sort!(unique!(codiedges_df[codiedges_df.is_tunnel, :cosource]))
        DataFrame(tunnels_from = isempty(cosources) ? missing : join(cosources, ' '),
                  ntunnels_from = length(cosources))
    end
    covertices_df = join(covertices_df, rename!(cosource_stats_df, :cosource => :covertex), on=:covertex, kind=:left)
    covertices_df = join(covertices_df, rename!(cotarget_stats_df, :cotarget => :covertex), on=:covertex, kind=:left)
    diedges_fwd_df = filter(r -> r.walkweight >= r.walkweight_rev, diedges_df)
    diedges_fwd_df[!, :is_reverse] .= false
    diedges_rev_df = filter(r -> r.walkweight < r.walkweight_rev, diedges_df)
    diedges_rev_df[!, :is_reverse] .= true
    rename!(diedges_rev_df, :source => :target, :cosource => :cotarget,
                            :target => :source, :cotarget => :cosource)
    coedges_df = by(vcat(diedges_fwd_df, diedges_rev_df),
                    [:component, :cosource, :cotarget]) do coedge_df
        res = DataFrame(
            has_tunnel = any(coedge_df.is_tunnel),
            has_walk = !all(coedge_df.is_tunnel),
            walkweight = any(r -> !r.is_tunnel && !r.is_reverse, eachrow(coedge_df)) ?
                    maximum(coedge_df[.!coedge_df.is_tunnel .& .!coedge_df.is_reverse, :walkweight]) : missing,
            walkweight_rev = any(r -> !r.is_tunnel && r.is_reverse, eachrow(coedge_df)) ?
                  maximum(coedge_df[.!coedge_df.is_tunnel .& coedge_df.is_reverse, :walkweight]) : missing
        )
        if hasproperty(coedge_df, :is_original)
            res[!, :has_original] .= any(r -> r.is_original && !r.is_reverse, eachrow(coedge_df))
            res[!, :has_original_rev] .= any(r -> r.is_original && r.is_reverse, eachrow(coedge_df))
        end
        if hasproperty(coedge_df, :diedge_type)
            orig1st = findfirst(r -> r.is_original && !r.is_reverse, eachrow(coedge_df))
            res[!, :has_original] .= !isnothing(orig1st)
            res[!, :target_type] .= isnothing(orig1st) ? missing : coedge_df.diedge_type[orig1st]
            if hasproperty(coedge_df, :interaction_type)
                res[!, :interaction_type] .= isnothing(orig1st) ? missing : coedge_df.interaction_type[orig1st]
            end
            orig1st_rev = findfirst(r -> r.is_original && r.is_reverse, eachrow(coedge_df))
            res[!, :has_original_rev] .= !isnothing(orig1st_rev)
            res[!, :source_type] .= isnothing(orig1st_rev) ? missing : coedge_df.diedge_type[orig1st_rev]
        end
        return res
    end
    return (components = components_df,
            vertices = vertices_df,
            covertices = covertices_df,
            diedges = diedges_df,
            coedges = coedges_df)
end
