function _graph_stats!(res::AbstractDataFrame,
    weights::AbstractArray{<:Number, N},
    walkweights::AbstractArray{<:Number, N},
    permweights::AbstractArray{<:Number},
    walkpermweights::AbstractArray{<:Number}
) where N
    permw_mtx = reshape(permweights, (prod(size(permweights)[1:N]), size(permweights, N+1)))
    permwalkw_mtx = reshape(walkpermweights, (prod(size(walkpermweights)[1:N]), size(permweights, N+1)))
    res.weight = copy(weights) |> vec
    res.permweight_mean = mean(permw_mtx, dims=2) |> vec
    res.permweight_std = std(permw_mtx, dims=2) |> vec
    res.permweight_median = median(permw_mtx, dims=2) |> vec
    res.permweight_mad = [mad(weights) for weights in eachrow(permw_mtx)]
    res.walkweight = copy(walkweights) |> vec
    res.walkpermweight_mean = mean(permwalkw_mtx, dims=2) |> vec
    res.walkpermweight_std = std(permwalkw_mtx, dims=2) |> vec
    res.walkpermweight_median = median(permwalkw_mtx, dims=2) |> vec
    res.walkpermweight_mad = [mad(walkweights) for walkweights in eachrow(permwalkw_mtx)]
    res.walkweight_delta = res.walkweight - res.walkpermweight_median
    res.nless_weight = [sum(>(w), permw)
                        for (w, permw) in zip(res.weight, eachrow(permw_mtx))]
    res.ngreater_weight = [sum(<(w), permw)
                           for (w, permw) in zip(res.weight, eachrow(permw_mtx))]
    res.nless_walkweight = [sum(>(w), permw)
                            for (w, permw) in zip(res.walkweight, eachrow(permwalkw_mtx))]
    res.ngreater_walkweight = [sum(<(w), permw)
                               for (w, permw) in zip(res.walkweight, eachrow(permwalkw_mtx))]
    return res
end

function vertex_stats(weights::AbstractVector{<:Number},
                      walkweights::AbstractVector{<:Number},
                      permweights::AbstractMatrix{<:Number},
                      walkpermweights::AbstractMatrix{<:Number})
    (length(weights) == size(permweights, 1)) ||
        throw(DimensionMismatch("Original weights length ($(length(weights))) doesn't match the permuted weights length $(size(permweights, 1))"))
    (length(weights) == length(walkweights)) ||
        throw(DimensionMismatch("Weights length ($(length(weights))) doesn't match the random walk weights length $(length(walkweights))"))
    (length(weights) == size(walkpermweights, 1)) ||
        throw(DimensionMismatch("Weights length ($(length(weights))) doesn't match the permuted random walk weights length $(size(walkpermweights, 1))"))
    (size(permweights, 2) == size(walkpermweights, 2)) ||
        throw(DimensionMismatch("Weights permutations count ($(size(permweights, 2))) doesn't match the random walk permutations count $(size(walkpermweights, 2))"))

    res = _graph_stats!(DataFrame(vertex = eachindex(weights)),
                        weights, walkweights, permweights, walkpermweights)
    return res
end

function conncomponents_stats(
    conncomps::IndicesPartition,
    vertex_stats::Union{AbstractDataFrame, Nothing} = nothing;
    average_weights::Bool = true,
    weight_tests::Bool = false,
    mannwhitney_tests::Bool = false
)
    isnothing(vertex_stats) || (vertex_stats.vertex == 1:nrow(vertex_stats)) ||
        throw(ArgumentError("vertex_stats should be ordered by vertex IDs from 1 to N without gaps"))
    res = DataFrame(
        component = eachindex(conncomps),
        nvertices = length.(conncomps)
    )
    isnothing(vertex_stats) && return res

    if mannwhitney_tests
        if weight_tests
            pvalues_weight_mw = sizehint!(Vector{Float64}(), length(conncomps))
        end
        pvalues_walkweight_mw = sizehint!(Vector{Float64}(), length(conncomps))
    end
    if average_weights
        weight_medians = sizehint!(Vector{Float64}(), length(conncomps))
        walkweight_medians = sizehint!(Vector{Float64}(), length(conncomps))
        weight_means = sizehint!(Vector{Float64}(), length(conncomps))
        walkweight_means = sizehint!(Vector{Float64}(), length(conncomps))
        permweight_medians = sizehint!(Vector{Float64}(), length(conncomps))
        walkpermweight_medians = sizehint!(Vector{Float64}(), length(conncomps))
        permweight_means = sizehint!(Vector{Float64}(), length(conncomps))
        walkpermweight_means = sizehint!(Vector{Float64}(), length(conncomps))
    end
    if average_weights || mannwhitney_tests
        comp_weights = Vector{Float64}()
        comp_walkweights = Vector{Float64}()
        comp_permweights = Vector{Float64}()
        comp_walkpermweights = Vector{Float64}()
        for (i, comp) in enumerate(conncomps)
            empty!(comp_weights)
            empty!(comp_walkweights)
            empty!(comp_permweights)
            empty!(comp_walkpermweights)
            @inbounds for v in comp
                push!(comp_weights, vertex_stats.weight[v])
                push!(comp_walkweights, vertex_stats.walkweight[v])
                push!(comp_permweights, vertex_stats.permweight_median[v])
                push!(comp_walkpermweights, vertex_stats.walkpermweight_median[v])
            end
            if average_weights
                push!(weight_medians, median(comp_weights))
                push!(weight_means, mean(comp_weights))
                push!(permweight_medians, median(comp_permweights))
                push!(permweight_means, mean(comp_permweights))
                push!(walkweight_medians, median(comp_walkweights))
                push!(walkweight_means, mean(comp_walkweights))
                push!(walkpermweight_medians, median(comp_walkpermweights))
                push!(walkpermweight_means, mean(comp_walkpermweights))
            end
            if mannwhitney_tests
                if weight_tests
                    push!(pvalues_weight_mw,
                          pvalue(MannWhitneyUTest(comp_weights, comp_permweights), tail=:right))
                end
                push!(pvalues_walkweight_mw,
                      pvalue(MannWhitneyUTest(comp_walkweights, comp_walkpermweights), tail=:right))
            end
        end
    end

    ntotalhits = sum(vertex_stats.is_hit)
    res[!, :nhits] = sum.(v -> @inbounds(vertex_stats.is_hit[v]), conncomps)
    res[!, :pvalue_walkweight_fisher] = pvalue.(
        Hypergeometric.(ntotalhits, nrow(vertex_stats) - ntotalhits, res.nvertices),
        res.nhits, tail=:right)
    if mannwhitney_tests
        if weight_tests
            res.pvalue_weight_mw = pvalues_weight_mw
        end
        res.pvalue_walkweight_mw = pvalues_walkweight_mw
    end
    if average_weights
        res.weight_median = weight_medians
        res.permweight_median = permweight_medians
        res.weight_mean = weight_means
        res.permweight_mean = permweight_means
        res.walkweight_median = walkweight_medians
        res.walkpermweight_median = walkpermweight_medians
        res.walkweight_mean = walkweight_means
        res.walkpermweight_mean = walkpermweight_means
    end
    return res
end

function treecut_stats(tree::SCCTree,
                       vertex_stats::Union{AbstractDataFrame, Nothing}=nothing;
                       walkmatrix::Union{AbstractMatrix, Nothing}=nothing,
                       sources::Union{AbstractVector, Nothing}=nothing,
                       sinks::Union{AbstractVector, Nothing}=nothing,
                       top_count::Integer=5,
                       mannwhitney_tests::Bool = false,
                       pvalue_mw_max::Number=0.05,
                       pvalue_fisher_max::Number=0.05,
                       pools::Union{ObjectPools, Nothing} = nothing,
                       nflows_ratio::Number=0.9)
    ptnpool = objpool(pools, IndicesPartition)
    comps = borrow!(ptnpool)
    comps_nontriv = borrow!(ptnpool)
    ncomps = Vector{Int}()
    ncomps_nontriv = Vector{Int}()
    maxcomp_sizes = Vector{Int}()
    topn_sizesum = Vector{Int}()
    intpool = arraypool(pools, Int)
    comps_size = borrow!(intpool)
    comps_perm = borrow!(intpool)
    if vertex_stats !== nothing
        nsignif_fisher = Vector{Int}()
        signif_sizesum_fisher = Vector{Int}()
        if mannwhitney_tests
            nsignif_mw = Vector{Int}()
            signif_sizesum_mw = Vector{Int}()
        end
    end
    if !isnothing(sources)
        topn_nsources = Vector{Int}()
        sourceset = Set(sources)
    end
    if !isnothing(sinks)
        topn_nsinks = Vector{Int}()
        sinkset = Set(sinks)
    end
    if !isnothing(sources) && !isnothing(sinks)
        if isnothing(walkmatrix)
            @warn "To count source-sink flows, walkmatrix must be specified"
            nflows_v = nothing
        else
            #iwalkmatrix, weights = indexvalues!(borrow!(arraypool(pools, Int32), size(walkmatrix)),
            #                                borrow!(arraypool(pools, eltype(walkmatrix))),
            #                                walkmatrix, EdgeTest{eltype(walkmatrix)}(rev=tree.rev))
            nflows_v = Vector{Int}()
            ncompflows_v = Vector{Int}()
            ncompsources_v = Vector{Int}()
            ncompsinks_v = Vector{Int}()
            flow_avglen_v = Vector{Float64}()
            flow_avgweight_v = Vector{Float64}()
            compflow_avglen_v = Vector{Float64}()
            compflow_avgweight_v = Vector{Float64}()
            flow_dist_v = Vector{Float64}()
            compflow_dist_v = Vector{Float64}()
            #if !isempty(tree.thresholds)
            #    iminthresh = searchsortedfirst(weights, first(tree.thresholds))
            #    imaxthresh = searchsortedfirst(weights, last(tree.thresholds))
            #else
            #    iminthresh = imaxthresh = 0
            #end
        end
    else
        nflows_v = nothing
    end
    lastcompsquares = 0
    #for i in 100:120
    #    thresh = tree.thresholds[i]
    for (i, thresh) in enumerate(tree.thresholds)
        cut!(comps, tree, thresh)
        push!(maxcomp_sizes, maximum(length, comps))
        push!(ncomps, length(comps))
        push!(ncomps_nontriv, sum(comp -> length(comp)>1, comps))
        map!(length, resize!(comps_size, length(comps)), comps)
        sortperm!(resize!(comps_perm, length(comps)), comps_size, rev=true)
        topn_pos = min(top_count, length(comps))
        topn_comp_ixs = view(comps_perm, 1:topn_pos)
        topn_comps = view(comps, topn_comp_ixs)
        push!(topn_sizesum, sum(view(comps_size, topn_comp_ixs)))
        if !isnothing(sources)
            push!(topn_nsources, sum(comp -> sum(i -> in(i, sourceset), comp), topn_comps))
        end
        if !isnothing(sinks)
            push!(topn_nsinks, sum(comp -> sum(i -> in(i, sinkset), comp), topn_comps))
        end
        if !isnothing(nflows_v)
            compsquares = sum(l -> ifelse(l > 1, abs2(float(l)), 0.2), comps_size)
            # check if the components changed significantly enough
            if (lastcompsquares == 0) || (i == length(tree.thresholds)) ||
               (compsquares/lastcompsquares <= nflows_ratio)
                lastcompsquares = compsquares
                #ithresh = searchsortedfirst(weights, thresh)
                #@assert (ithresh <= length(weights)) && (weights[ithresh] == thresh)
                foreach(sort!, comps) # sorting improves condense!(iwalkmatrix) performace
                flowstats =
                    nflows(comps, walkmatrix, sources, sinks, EdgeTest{eltype(walkmatrix)}(threshold=thresh), pools)
                nvtxflows_max = length(sources)*length(sinks)
                ncompflows_max = flowstats.ncompsources*flowstats.ncompsinks

                push!(ncompsources_v, flowstats.ncompsources)
                push!(ncompsinks_v, flowstats.ncompsinks)
                push!(nflows_v, flowstats.nflows)
                push!(ncompflows_v, flowstats.ncompflows)
                push!(flow_avglen_v, flowstats.flowlen_sum/flowstats.nflows)
                push!(compflow_avglen_v, flowstats.compflowlen_sum/flowstats.ncompflows)
                push!(flow_avgweight_v, flowstats.floweight_sum / nvtxflows_max)
                push!(compflow_avgweight_v, flowstats.compfloweight_sum / ncompflows_max)
                push!(flow_dist_v, ((nvtxflows_max - flowstats.nflows) * (length(comps) + 1) + flowstats.flowlen_sum) / nvtxflows_max)
                push!(compflow_dist_v, ((ncompflows_max - flowstats.ncompflows) * (length(comps) + 1) + flowstats.compflowlen_sum) / ncompflows_max)
            else # duplicate the last nflows
                push!(ncompsources_v, last(ncompsources_v))
                push!(ncompsinks_v, last(ncompsinks_v))
                push!(nflows_v, last(nflows_v))
                push!(ncompflows_v, last(ncompflows_v))
                push!(flow_avglen_v, last(flow_avglen_v))
                push!(compflow_avglen_v, last(compflow_avglen_v))
                push!(flow_avgweight_v, last(flow_avgweight_v))
                push!(compflow_avgweight_v, last(compflow_avgweight_v))
                push!(flow_dist_v, last(flow_dist_v))
                push!(compflow_dist_v, last(compflow_dist_v))
            end
        end
        if vertex_stats !== nothing
            filter!(comp -> length(comp) > 1, comps_nontriv, comps)
            compstats_df = conncomponents_stats(comps_nontriv, vertex_stats, average_weights=false, mannwhitney_tests=mannwhitney_tests)
            push!(nsignif_fisher, sum(<=(pvalue_fisher_max), compstats_df.pvalue_walkweight_fisher))
            push!(signif_sizesum_fisher, sum(r -> ifelse(r.pvalue_walkweight_fisher <= pvalue_fisher_max, r.nvertices, 0),
                                             eachrow(compstats_df)))
            if mannwhitney_tests
                push!(nsignif_mw, sum(<=(pvalue_mw_max), compstats_df.pvalue_walkweight_mw))
                push!(signif_sizesum_mw, sum(r -> ifelse(r.pvalue_walkweight_mw <= pvalue_mw_max, r.nvertices, 0),
                                             eachrow(compstats_df)))
            end
        end
    end
    res = DataFrame(
        threshold = tree.thresholds[1:length(ncomps)],
        ncomponents = ncomps,
        ncomponents_nontrivial = ncomps_nontriv,
        maxcomponent_size = maxcomp_sizes,
        log10_maxcomponent_size = log10.(maxcomp_sizes),
        topn_components_sizesum = topn_sizesum,
        log10_topn_components_sizesum = log10.(topn_sizesum)
    )
    if vertex_stats !== nothing
        res.ncomponents_signif_fisher = nsignif_fisher
        res.components_signif_sizesum_fisher = signif_sizesum_fisher
        if mannwhitney_tests
            res.ncomponents_signif_mw = nsignif_mw
            res.components_signif_sizesum_mw = signif_sizesum_mw
        end
    end
    if !isnothing(sources)
        res.topn_nsources = topn_nsources
    end
    if !isnothing(sinks)
        res.topn_nsinks = topn_nsinks
    end
    if !isnothing(nflows_v)
        res.ncompsources = ncompsources_v
        res.ncompsinks = ncompsinks_v
        res.nflows = nflows_v
        res.ncompflows = ncompflows_v
        res.flow_avglen = flow_avglen_v
        res.compflow_avglen = compflow_avglen_v
        res.flow_avgweight = flow_avgweight_v
        res.compflow_avgweight = compflow_avgweight_v
        res.flow_distance = flow_dist_v
        res.compflow_distance = compflow_dist_v
        # calculate avgweight quantiles
        # prepare vertex of weights for ECDF
        # the matrix has a lot of zeros, so to help sorting, we compress them into single entry
        mtxvals = sizehint!(Vector{eltype(walkmatrix)}(), length(walkmatrix))
        nzeros = 0
        for w in walkmatrix
            if w != 0
                push!(mtxvals, w)
            else
                nzeros += 1
            end
        end
        if nzeros > 0
            mtxweights = fill(1, length(mtxvals))
            push!(mtxvals, 0)
            push!(mtxweights, nzeros)
            walkmatrix_cdf = HierarchicalHotNet.ecdf(mtxvals, mtxweights, true)
        else
            walkmatrix_cdf = HierarchicalHotNet.ecdf(mtxvals, nothing, true)
        end
        res.flow_avgweight_qtl = walkmatrix_cdf.(res.flow_avgweight)
        res.compflow_avgweight_qtl = walkmatrix_cdf.(res.compflow_avgweight)

        #release!(arraypool(pools, Int32), iwalkmatrix)
        #release!(arraypool(pools, eltype(walkmatrix)), weights)
    end

    release!(ptnpool, comps_nontriv)
    release!(ptnpool, comps)
    release!(intpool, comps_perm)
    release!(intpool, comps_size)
    return res
end

const treecut_metrics = [
    :ncomponents, :ncomponents_nontrivial,
    :ncomponents_signif_mw, :ncomponents_signif_fisher,
    :components_signif_sizesum_mw, :components_signif_sizesum_fisher,
    :maxcomponent_size, :log10_maxcomponent_size,
    :topn_components_sizesum, :log10_topn_components_sizesum,
    :topn_nsources, :topn_nsinks,
    :ncompsources, :ncompsinks,
    :nflows, :ncompflows, :flow_avglen, :compflow_avglen,
    :flow_avgweight, :flow_avgweight_qtl,
    :compflow_avgweight, :compflow_avgweight_qtl,
    :flow_distance, :compflow_distance]

function bin_treecut_stats(
    cutstats_df::AbstractDataFrame;
    by_cols::Union{AbstractVector{Symbol}, Symbol, Nothing} = nothing,
    threshold_range::Union{NTuple{2, Float64}, Nothing} = nothing,
    threshold_nbins::Integer = 100,
    stat_cols::AbstractVector{Symbol} = intersect(treecut_metrics, propertynames(cutstats_df)),
)
    used_thresholds = threshold_range === nothing ? cutstats_df.threshold :
        filter(t -> threshold_range[1] <= t <= threshold_range[2],
               cutstats_df.threshold)
    used_thresholds_range = extrema(used_thresholds)
    threshold_bins = quantile(used_thresholds, 0.0:(1/threshold_nbins):1.0)
    threshold_bin_centers = 0.5 .* (threshold_bins[1:(length(threshold_bins)-1)] .+ threshold_bins[2:end])
    cutstats_df = isa(cutstats_df, DataFrame) ? copy(cutstats_df, copycols=false) : copy(cutstats_df)
    cutstats_df.threshold_bin = searchsortedlast.(Ref(threshold_bins), cutstats_df.threshold)

    if by_cols isa AbstractVector
        used_by_cols = push!(copy(by_cols), :threshold_bin)
    elseif by_cols isa Symbol
        used_by_cols = [by_cols, :threshold_bin]
    else
        used_by_cols = [:threshold_bin]
    end

    binstats_df = combine(groupby(cutstats_df, used_by_cols),
                          [col => median => col for col in stat_cols]...)
    binstats_df.threshold = [0 < bin <= threshold_nbins ? threshold_bin_centers[bin] : missing
                             for bin in binstats_df.threshold_bin]
    return binstats_df
end

quantile_suffix(q::Real) =
    reduce(replace, init=@sprintf("%02.1f", 100*q),
           [r"\.0$" => "",
            r"^(\d)\." => s"0\1",
            r"\.(\d)" => s"\1"])

function aggregate_treecut_binstats(
    binstats_df::AbstractDataFrame;
    by_cols::Union{AbstractVector{Symbol}, Symbol, Nothing} = nothing,
    quantiles::AbstractVector{<:Number} = [0.025, 0.25, 0.5, 0.75, 0.975],
    stat_cols::AbstractVector{Symbol} = intersect(treecut_metrics, propertynames(binstats_df))
)
    used_quantiles = sort!(copy(quantiles))
    med_pos = searchsortedfirst(used_quantiles, 0.5)
    if used_quantiles[med_pos] != 0.5
        insert!(used_quantiles, med_pos, 0.5)
    end

    thresholds_df = binstats_df[.!nonunique(binstats_df, [:threshold_bin, :threshold]),
                                   [:threshold_bin, :threshold]]
    @assert nrow(thresholds_df) == length(unique(thresholds_df.threshold_bin))

    if by_cols isa AbstractVector
        used_by_cols = push!(copy(by_cols), :threshold_bin)
    elseif by_cols isa Symbol
        used_by_cols = [by_cols, :threshold_bin]
    else
        used_by_cols = [:threshold_bin]
    end
    unique!(used_by_cols)
    deleteat!(used_by_cols, findall(==(:threshold), used_by_cols))

    binstats_df_grouped = groupby(binstats_df, used_by_cols)
    aggstats_df = reduce(vcat, [begin
        res = combine(binstats_df_grouped,
                      [col => (x -> any(x -> ismissing(x) || isnan(x), x) ? NaN : quantile(x, qtl)) => col
                       for col in stat_cols]...)
        res[!, :quantile] .= qtl
        res
    end for qtl in used_quantiles])
    leftjoin(aggstats_df, thresholds_df, on=:threshold_bin)
end

function extreme_treecut_binstats(
    binstats_df::AbstractDataFrame,
    perm_aggstats_df::AbstractDataFrame;
    extra_join_cols::Union{Nothing, AbstractVector{Symbol}} = nothing,
    stat_cols::AbstractVector{Symbol} = intersect(treecut_metrics, propertynames(binstats_df))
)
    perm_aggstats_cols = intersect(treecut_metrics, propertynames(perm_aggstats_df))
    join_cols = [:threshold, :threshold_bin]
    isnothing(extra_join_cols) || unique!(append!(join_cols, extra_join_cols))
    joinstats_df = leftjoin(select(binstats_df, [join_cols; stat_cols]),
                            select(perm_aggstats_df, [join_cols; :quantile; perm_aggstats_cols]),
                            on=join_cols, makeunique=true)
    by_cols = [:quantile]
    isnothing(extra_join_cols) || unique!(append!(by_cols, extra_join_cols))
    combine(groupby(joinstats_df, by_cols)) do df
        reduce(vcat, [begin
            perm_col = Symbol(col, "_1")
            res = reduce(vcat, [begin
                deltas = [ismissing(r[col]) || ismissing(r[perm_col]) ? 0.0 : r[col] - r[perm_col]
                          for r in eachrow(df)]
                delta_val, delta_pos = aggfun(deltas)
                DataFrame(
                    :metric => col,
                    :value => df[delta_pos, col],
                    :value_type => aggtype,
                    :permuted_value => df[delta_pos, perm_col],
                    :delta => delta_val,
                    :threshold_bin => df.threshold_bin[delta_pos],
                    :threshold => df.threshold[delta_pos]
                )
            end for (aggtype, aggfun) in ["min_delta" => findmin, "max_delta" => findmax]])
            res
        end for col in stat_cols])
    end
end
