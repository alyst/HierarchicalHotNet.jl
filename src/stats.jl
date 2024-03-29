# Statistics for network diffusion results and SCCTrees:
# comparison of the results based on real data and random permutations

# compute statistics for permutation-based weights
# and compare them with the real data-based ones
# the columns in `permweights` correspond to different permutations
function _permweight_stats(weights::AbstractVector, permweights::AbstractMatrix)
    @assert length(weights) == size(permweights, 1)
    means = similar(weights)
    stds = similar(weights)
    medians = similar(weights)
    mads = similar(weights)
    nless = similar(weights, Int)
    ngreater = similar(nless)
    Threads.@threads for i in axes(permweights, 1)
        permw = permweights[i, :]
        means[i] = mean(permw)
        stds[i] = std(permw)
        medians[i] = mean(permw)
        mads[i] = mad(permw)
        nless[i] = sum(>(weights[i]), permw)
        ngreater[i] = sum(<(weights[i]), permw)
    end
    return means, stds, medians, mads, nless, ngreater
end

function _graph_stats!(res::AbstractDataFrame,
    weights::AbstractArray{<:Number, N},
    walkweights::AbstractArray{<:Number, N},
    permweights::Union{AbstractArray{<:Number}, Nothing} = nothing,
    walkpermweights::Union{AbstractArray{<:Number}, Nothing} = nothing
) where N
    res.weight = copy(weights) |> vec
    if !isnothing(permweights)
        permw_mtx = reshape(permweights, (prod(size(permweights)[1:N]), size(permweights, N+1)))
        res.permweight_mean, res.permweight_std,
        res.permweight_median, res.permweight_mad,
        res.nless_weight, res.ngreater_weight = _permweight_stats(res.weight, permw_mtx)
        res.weight_delta = res.weight - res.permweight_median
    end

    res.walkweight = copy(walkweights) |> vec
    if !isnothing(walkpermweights)
        permwalkw_mtx = reshape(walkpermweights, (prod(size(walkpermweights)[1:N]), size(walkpermweights, N+1)))
        res.walkpermweight_mean, res.walkpermweight_std,
        res.walkpermweight_median, res.walkpermweight_mad,
        res.nless_walkweight, res.ngreater_walkweight = _permweight_stats(res.walkweight, permwalkw_mtx)
        res.walkweight_delta = res.walkweight - res.walkpermweight_median
    end
    return res
end

"""
    vertex_stats(weights::AbstractVector{<:Number},
                 walkweights::AbstractVector{<:Number},
                 [permweights::AbstractMatrix{<:Number}],
                 [walkpermweights::AbstractMatrix{<:Number}]) -> DataFrame

Calculates statistics for the permuted vertex weights distribution and how it
is different from the actual weights.

Returns the data frame with per-vertex mean, standard deviation,
median, MAD and the probability that permuted value is greater/lower than the corresponding real value
for the weights of the original matrix (`weights`) as well as random walk matrix weights (`walkweights`).

### Parameters
* `weights`: weights of the vertices in the original network
* `walkweights`: weights of the vertices after network diffusion analysis
  (stationary random walk distribution)
* `permweights`: matrix of permuted weights; rows correspond to vertices,
  columns -- to permutations
* `walkpermweights`: matrix of vertex weights based on network diffusion analysis
  using `permweights` as input; rows correspond to vertices, columns to permutations
"""
function vertex_stats(weights::AbstractVector{<:Number},
                      walkweights::AbstractVector{<:Number},
                      permweights::Union{AbstractMatrix{<:Number}, Nothing} = nothing,
                      walkpermweights::Union{AbstractMatrix{<:Number}, Nothing} = nothing
)
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


"""
    diedge_stats(weights::AbstractVector{<:Number},
                 walkweights::AbstractVector{<:Number},
                 [permweights::AbstractMatrix{<:Number}],
                 [walkpermweights::AbstractMatrix{<:Number}]) -> DataFrame

Calculates statistics for the directed edges permuted weights distribution and how it
is different from the actual weights of directed edges.

The output is similar to [`HierarchicalHotNet.vertex_stats`](@ref) for vertices.
"""
function diedge_stats(weights::AbstractMatrix{<:Number},
                      walkweights::AbstractMatrix{<:Number},
                      permweights::Union{AbstractArray{<:Number, 3}, Nothing} = nothing,
                      walkpermweights::Union{AbstractArray{<:Number, 3}, Nothing} = nothing
)
    check_square(weights, "weight matrix")
    (size(weights) == size(walkweights)) ||
        throw(DimensionMismatch("Original weights size ($(size(weights))) doesn't match the random walk weights size $(size(walkweights))"))
    (size(weights) == size(permweights, (1, 2))) ||
        throw(DimensionMismatch("Original weights size ($(size(weights))) doesn't match the permuted weights size $(size(permweights, (1, 2)))"))
    (size(weights) == size(walkpermweights, (1, 2))) ||
        throw(DimensionMismatch("Original weights size ($(size(weights))) doesn't match the permuted random walk weights size $(size(walkpermweights, (1, 2)))"))

    res = DataFrame(
            src = getindex.(Ref(CartesianIndices(weights)), 2),
            dest = getindex.(Ref(CartesianIndices(weights)), 1))
    _graph_stats!(res, weights, walkweights, permweights, walkpermweights)
    return res
end

function diedge_stats(nvertices::Integer, diedge_indices::AbstractVector{<:Integer},
                      weights::AbstractVector{<:Number},
                      walkweights::AbstractVector{<:Number},
                      permweights::Union{AbstractMatrix{<:Number}, Nothing} = nothing,
                      walkpermweights::Union{AbstractMatrix{<:Number}, Nothing} = nothing
)
    (length(diedge_indices) == length(weights)) ||
        throw(DimensionMismatch("Diedge indices count ($(length(diedge_indices))) doesn't match the weights count $(length(weights))"))
    (length(weights) == length(walkweights)) ||
        throw(DimensionMismatch("Weights length ($(length(weights))) doesn't match the random walk weights length $(length(walkweights))"))
    isnothing(permweights) || (length(weights) == size(permweights, 1)) ||
        throw(DimensionMismatch("Original weights length ($(length(weights))) doesn't match the permuted weights length $(size(permweights, 1))"))
    isnothing(walkpermweights) || (length(weights) == size(walkpermweights, 1)) ||
        throw(DimensionMismatch("Original weights length ($(length(weights))) doesn't match the permuted random walk weights length $(size(walkpermweights, 1))"))

    all_diedges = CartesianIndices((nvertices, nvertices))
    res = DataFrame(
        src = [all_diedges[i][2] for i in diedge_indices],
        dest = [all_diedges[i][1] for i in diedge_indices])
    _graph_stats!(res, weights, walkweights, permweights, walkpermweights)
    return res
end

"""
Calculate per-connected component statistics.
"""
function conncomponents_stats(
    conncomps::IndicesPartition,
    vertex_weights::Union{AbstractVector, Nothing} = nothing,
    vertex_walkweights::Union{AbstractVector, Nothing} = nothing,
    perm_vertex_weights::Union{AbstractMatrix, Nothing} = nothing,
    perm_vertex_walkweights::Union{AbstractMatrix, Nothing} = nothing;
    average_weights::Bool = true,
    weight_tests::Bool = false,
    mannwhitney_tests::Bool = false
)
    res = DataFrame(
        component = eachindex(conncomps),
        nvertices = length.(conncomps)
    )
    isnothing(vertex_weights) && return res

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
                isnothing(vertex_weights) || push!(comp_weights, view(vertex_weights, v))
                isnothing(comp_permweights) || append!(comp_permweights, view(perm_vertex_weights, v, :) |> vec)
                isnothing(vertex_walkweights) || push!(comp_walkweights, view(vertex_walkweights, v))
                isnothing(comp_walkpermweights) || append!(comp_walkpermweights, view(perm_vertex_walkweights, v, :) |> vec)
            end
            if average_weights
                push!(weight_medians, median(comp_weights))
                push!(weight_means, mean(comp_weights))
                push!(walkweight_medians, median(comp_walkweights))
                push!(walkweight_means, mean(comp_walkweights))
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

"""
    treecut_stats(tree::SCCTree;
                  [walkmatrix::AbstractMatrix],
                  [maxweight::Number],
                  [sources], [sinks], [sourcesinkweights], [top_count],
                  [pools]) -> DataFrame

Calculate SCC network statistic for each cutting threshold of `tree`.
"""
function treecut_stats(tree::SCCTree;
                       walkmatrix::Union{AbstractMatrix, Nothing}=nothing,
                       maxweight::Union{Number, Nothing}=nothing,
                       sources::Union{AbstractVector, Nothing}=nothing,
                       sinks::Union{AbstractVector, Nothing}=nothing,
                       sourcesinkweights::Union{AbstractMatrix, Nothing}=nothing,
                       top_count::Integer=5,
                       pools::Union{ObjectPools, Nothing} = nothing)
    ptnpool = objpool(pools, IndicesPartition)
    comps = borrow!(ptnpool)
    comps_nontriv = borrow!(ptnpool)
    intpool = arraypool(pools, Int)
    comps_size = borrow!(intpool)
    comps_perm = borrow!(intpool)

    res = DataFrame(
        threshold = similar(tree.thresholds, 0),
        ncomponents = Int[],
        ncomponents_nontrivial = Int[],
        maxcomponent_size = Int[],
        topn_components_sizesum = Int[],
    )
    if !isnothing(walkmatrix) && !isnothing(sources) && !isnothing(sinks)
        res.nflows = Int[]
        res.ncompflows = Int[]
        res.flow_avglen = Float64[]
        res.flow_avginvlen = Float64[]
        res.compflow_avglen = Float64[]
        res.compflow_avginvlen = Float64[]
        res.flow_avgweight = Float64[]
        res.flow_avghopweight = Float64[]
        res.compflow_avgweight = Float64[]
        res.flow_avgminedgeweight = Float64[]
        res.compflow_avgminedgeweight = Float64[]
        res.flow_distance = Float64[]
        res.compflow_distance = Float64[]

        peelit = eachflowspeel(tree, walkmatrix, sources, sinks, sortvertices=true, verbose=false)
        #iwalkmatrix, weights = indexvalues!(borrow!(arraypool(pools, Int32), size(walkmatrix)),
        #                                borrow!(arraypool(pools, eltype(walkmatrix))),
        #                                walkmatrix, EdgeTest{eltype(walkmatrix)}(rev=tree.rev))

        #if !isempty(tree.thresholds)
        #    iminthresh = searchsortedfirst(weights, first(tree.thresholds))
        #    imaxthresh = searchsortedfirst(weights, last(tree.thresholds))
        #else
        #    iminthresh = imaxthresh = 0
        #end
    elseif !isnothing(walkmatrix) || !isnothing(sources) || !isnothing(sinks)
        @warn "To count source-sink flows, walkmatrix, sources and sinks must be specified"
        peelit = nothing
    end
    if !isnothing(sources)
        insources = in(Set(sources))
        res.ncompsources = Int[]
        res.topn_nsources = Int[]
    end
    if !isnothing(sinks)
        insinks = in(Set(sinks))
        res.ncompsinks = Int[]
        res.topn_nsinks = Int[]
    end

    newrow = Dict{Symbol, Any}()
    ithresholdpos = 0
    #for i in 100:120
    #    thresh = tree.thresholds[i]
    for (i, thresh) in enumerate(tree.thresholds)
        cut!(comps, tree, thresh)
        map!(length, resize!(comps_size, length(comps)), comps)
        sortperm!(resize!(comps_perm, length(comps)), comps_size, rev=true)
        topn_pos = min(top_count, length(comps))
        topn_comp_ixs = view(comps_perm, 1:topn_pos)
        topn_comps = view(comps, topn_comp_ixs)
        push!(newrow, :threshold => thresh,
                      :maxcomponent_size => maximum(length, comps),
                      :ncomponents => length(comps),
                      :ncomponents_nontrivial => count(comp -> length(comp)>1, comps),
                      :topn_components_sizesum => sum(view(comps_size, topn_comp_ixs)))
        if hasproperty(res, :ncompsources)
            push!(newrow, :topn_nsources => sum(comp -> count(insources, comp), topn_comps),
                          :ncompsources => count(comp -> any(insources, comp), comps))
        end
        if hasproperty(res, :ncompsinks)
            push!(newrow, :topn_nsinks => sum(comp -> count(insinks, comp), topn_comps),
                          :ncompsinks => count(comp -> any(insinks, comp), comps))
        end
        if hasproperty(res, :nflows)
            flows, ithresholdpos = iterate(peelit, ithresholdpos)
            @assert threshold(flows) == thresh
            if flows.modified # if flows unchanged, newrow already contains uptodate information
                flstats = flowstats(flows, maxweight=maxweight, sourcesinkweights=sourcesinkweights)
                nvtxflows_max = length(sources)*length(sinks)
                ncompflows_max = (newrow[:ncompsources]::Int)*(newrow[:ncompsinks]::Int)

                push!(newrow, :nflows => flstats.nflows,
                    :ncompflows => flstats.ncompflows,
                    :flow_avglen => flstats.flowlen_sum/flstats.nflows,
                    :flow_avginvlen => flstats.flowinvlen_sum/nvtxflows_max,
                    :compflow_avglen => flstats.compflowlen_sum/flstats.ncompflows,
                    :compflow_avginvlen => flstats.compflowinvlen_sum/ncompflows_max,
                    :flow_avgweight => flstats.floweight_sum / nvtxflows_max,
                    :flow_avghopweight => flstats.flowavghopweight_sum / nvtxflows_max,
                    :compflow_avgweight => flstats.compfloweight_sum / ncompflows_max,
                    :flow_avgminedgeweight => flstats.flow_minedgeweight_sum / nvtxflows_max,
                    :compflow_avgminedgeweight => flstats.compflow_minedgeweight_sum / ncompflows_max,
                    :flow_distance => ((nvtxflows_max - flstats.nflows) * (length(comps) + 1) + flstats.flowlen_sum) / nvtxflows_max,
                    :compflow_distance => ((ncompflows_max - flstats.ncompflows) * (length(comps) + 1) + flstats.compflowlen_sum) / ncompflows_max)
            end
        end
        push!(res, newrow)
    end
    res.log10_maxcomponent_size = log10.(res.maxcomponent_size)
    res.log10_topn_components_sizesum = log10.(res.topn_components_sizesum)
    if hasproperty(res, :nflows)
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
        #release!(arraypool(pools, Int32), iwalkmatrix)
        #release!(arraypool(pools, eltype(walkmatrix)), weights)
    end

    release!(ptnpool, comps_nontriv)
    release!(ptnpool, comps)
    release!(intpool, comps_perm)
    release!(intpool, comps_size)
    return res
end

"""
    treecut_compstats(tree::SCCTree,
                      vertex_weights::AbstractVector,
                      vertex_walkweights::AbstractVector,
                      perm_vertex_weights::AbstractMatrix,
                      perm_vertex_walkweights::AbstractMatrix;
                      [mannwhitney_tests::Bool],
                      [pvalue_mw_max::Number],
                      [pvalue_fisher_max::Number],
                      [pools]) -> DataFrame

Calculate SCC network statistic for each cutting threshold of `tree`.
"""
function treecut_compstats(tree::SCCTree,
    vertex_weights::AbstractVector,
    vertex_walkweights::AbstractVector,
    perm_vertex_weights::AbstractMatrix,
    perm_vertex_walkweights::AbstractMatrix;
    mannwhitney_tests::Bool = false,
    pvalue_mw_max::Number=0.05,
    pvalue_fisher_max::Number=0.05,
    pools::Union{ObjectPools, Nothing} = nothing
)
    ptnpool = objpool(pools, IndicesPartition)
    comps = borrow!(ptnpool)
    comps_nontriv = borrow!(ptnpool)
    nsignif_fisher = Vector{Int}()
    signif_sizesum_fisher = Vector{Int}()
    if mannwhitney_tests
        nsignif_mw = Vector{Int}()
        signif_sizesum_mw = Vector{Int}()
    end
    for (i, thresh) in enumerate(tree.thresholds)
        cut!(comps, tree, thresh)
        filter!(comp -> length(comp) > 1, comps_nontriv, comps)
        compstats_df = conncomponents_stats(comps_nontriv, vertex_weights, vertex_walkweights,
                perm_vertex_weights, perm_vertex_walkweights,
                average_weights=false, mannwhitney_tests=mannwhitney_tests)
        push!(nsignif_fisher, sum(<=(pvalue_fisher_max), compstats_df.pvalue_walkweight_fisher))
        push!(signif_sizesum_fisher, sum(r -> ifelse(r.pvalue_walkweight_fisher <= pvalue_fisher_max, r.nvertices, 0),
                          eachrow(compstats_df)))
        if mannwhitney_tests
            push!(nsignif_mw, sum(<=(pvalue_mw_max), compstats_df.pvalue_walkweight_mw))
            push!(signif_sizesum_mw, sum(r -> ifelse(r.pvalue_walkweight_mw <= pvalue_mw_max, r.nvertices, 0),
                            eachrow(compstats_df)))
        end
    end
    res = DataFrame(
            threshold = tree.thresholds[1:length(ncomps)],
            ncomponents_signif_fisher = nsignif_fisher,
            components_signif_sizesum_fisher = signif_sizesum_fisher
    )
    if mannwhitney_tests
        res.ncomponents_signif_mw = nsignif_mw
        res.components_signif_sizesum_mw = signif_sizesum_mw
    end
    release!(ptnpool, comps_nontriv)
    release!(ptnpool, comps)
    return res
end

"""
[`treecut_stats()`](@ref HierarchicalHotNet.treecut_stats) metrics (dataframe columns) to consider for
[`bin_treecut_stats()`](@ref HierarchicalHotNet.bin_treecut_stats) and
[`extreme_treecut_stats()`](@ref HierarchicalHotNet.extreme_treecut_stats).
"""
const TreecutMetrics = [
    :ncomponents, :ncomponents_nontrivial,
    :ncomponents_signif_mw, :ncomponents_signif_fisher,
    :components_signif_sizesum_mw, :components_signif_sizesum_fisher,
    :maxcomponent_size, :log10_maxcomponent_size,
    :topn_components_sizesum, :log10_topn_components_sizesum,
    :topn_nsources, :topn_nsinks,
    :ncompsources, :ncompsinks,
    :nflows, :ncompflows,
    :flow_avglen, :compflow_avglen,
    :flow_avginvlen, :compflow_avginvlen,
    :flow_avgweight, :compflow_avgweight,
    :flow_avgminedgeweight, :compflow_avgminedgeweight,
    :flow_avghopweight,
    :flow_distance, :compflow_distance]

function add_bins!(df::AbstractDataFrame, col::Symbol, bin_bounds::AbstractVector{<:Number})
    df[:, string(col, "_bin")] .= searchsortedlast.(Ref(bin_bounds), df[:, col])
    return df
end

function add_bins!(df::DataFrame, col::Symbol, bin_bounds::AbstractVector{<:Number})
    df[!, string(col, "_bin")] = searchsortedlast.(Ref(bin_bounds), df[!, col])
    return df
end

"""
    bin_treecut_stats(cutstats_df::AbstractDataFrame) -> DataFrame

Bin treecut thresholds and calculate average statistics in each bin.

Takes the output of [`HierarchicalHotNet.treecut_stats`](@ref) from
multiple SCC trees (discriminated by `by_cols`), identifies the bind for
treecut thresholds and calculates the average metric values (`stat_cols`)
within each bin.
"""
function bin_treecut_stats(
    cutstats_df::AbstractDataFrame;
    by_cols::Union{AbstractVector{Symbol}, Symbol, Nothing} = nothing,
    threshold_range::Union{NTuple{2, Float64}, Nothing} = nothing,
    threshold_nbins::Integer = 100,
    threshold_maxbinwidth::Union{Nothing, Number} = nothing,
    stat_cols::AbstractVector{Symbol} = intersect(TreecutMetrics, propertynames(cutstats_df)),
)
    used_thresholds = threshold_range === nothing ? cutstats_df.threshold :
        filter(t -> threshold_range[1] <= t <= threshold_range[2],
               cutstats_df.threshold)
    sort!(used_thresholds)
    thresh_bins = sizehint!(similar(used_thresholds, 0), threshold_nbins)
    lastbinindex = 0
    for (i, thresh) in enumerate(used_thresholds)
        if isempty(thresh_bins) || # first threshold
        # enough thresholds in the bin or last threshold
        (last(thresh_bins) < thresh) && (((i + 1 - lastbinindex) * threshold_nbins >= length(used_thresholds)) || (i == length(used_thresholds))) ||
        # too wide threshold
        !isnothing(threshold_maxbinwidth) && ((thresh - last(thresh_bins)) >= threshold_maxbinwidth)
            push!(thresh_bins, thresh)
            lastbinindex = i
        end
    end

    # copy to avoid modifying original frame
    cutstats_df = isa(cutstats_df, DataFrame) ?
            copy(cutstats_df, copycols=false) :
            select(cutstats_df, [stat_cols; by_cols; [:threshold]])
    add_bins!(cutstats_df, :threshold, thresh_bins)

    if by_cols isa AbstractVector
        used_by_cols = push!(copy(by_cols), :threshold_bin)
    elseif by_cols isa Symbol
        used_by_cols = [by_cols, :threshold_bin]
    else
        used_by_cols = [:threshold_bin]
    end

    binstats_df = combine(groupby(cutstats_df, used_by_cols),
                          [col => median => col for col in stat_cols]...)
    binstats_df.threshold_binmin = missings(Float64, nrow(binstats_df))
    binstats_df.threshold_binmid = missings(Float64, nrow(binstats_df))
    binstats_df.threshold_binmax = missings(Float64, nrow(binstats_df))
    @inbounds for (i, bin) in enumerate(binstats_df.threshold_bin)
        (0 < bin < length(thresh_bins)) || continue
        binstats_df.threshold_binmin[i] = binmin = thresh_bins[bin]
        binstats_df.threshold_binmax[i] = binmax = thresh_bins[bin+1]
        binstats_df.threshold_binmid[i] = 0.5*(binmin + binmax)
    end
    return binstats_df
end

quantile_suffix(q::Real) =
    reduce(replace, init=@sprintf("%02.1f", 100*q),
           [r"\.0$" => "",
            r"^(\d)\." => s"0\1",
            r"\.(\d)" => s"\1"])

"""
    aggregate_treecut_binstats(binstats_df::AbstractDataFrame) -> DataFrame

Aggregate the binned treecut statistics across multiple trees.

Takes `binstats_df`, the output of [`HierarchicalHotNet.bin_treecut_stats`](@ref), and
calculates the metric values for the specified quantiles.
"""
function aggregate_treecut_binstats(
    binstats_df::AbstractDataFrame;
    by_cols::Union{AbstractVector{Symbol}, Symbol, Nothing} = nothing,
    quantiles::AbstractVector{<:Number} = [0.025, 0.25, 0.5, 0.75, 0.975],
    metric_cols::AbstractVector{Symbol} = intersect(TreecutMetrics, propertynames(binstats_df))
)
    used_quantiles = sort!(copy(quantiles))
    med_pos = searchsortedfirst(used_quantiles, 0.5)
    if used_quantiles[med_pos] != 0.5
        insert!(used_quantiles, med_pos, 0.5)
    end

    thresholds_df = select(binstats_df, :threshold_bin, :threshold_binmid) |> unique!
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
                       for col in metric_cols]...)
        res[!, :quantile] .= qtl
        res
    end for qtl in used_quantiles])
    leftjoin(aggstats_df, thresholds_df, on=:threshold_bin)
end

"""
    extreme_treecut_stats(stats_df::AbstractDataFrame) -> DataFrame

Calculate the cut threshold and corresponding metric value, where the difference
between real (taken from `stats_df`) and permutation metrics (taken from `perm_aggstats_df`)
are maximal/minimal (depending on the metric).

## Arguments
  * `stats_df`: tree statistics calculated by [`treecut_stats`](@ref)
  * `perm_aggstats_df`: aggregated binned permutated tree statistics calculated by [`aggregate_treecut_binstats`](@ref)
  * `extra_join_cols`: optional columns, in addition to `:threshold_bin` to use for joining `stats_df` and `perm_aggstats_df`
  * `metric_cols`: columns of `stats_df` and `perm_aggstats_df` containing treecut metrics
     to consider for threshold calculation (see [`TreecutMetrics`](@ref HierarchicalHotNet.TreecutMetrics))
  * `start_maxquantile`: if specified, calculates (in addition to minimal and maximal metric)
     the metric corresponding to the given quantile as well as ``1 - quantile``
  * `threshold_range`: if given, contrains metric statistic calculation to given min/max thresholds
  * `threshold_weight`: optional function that takes stats_df row and returns
    the prior weight of the corresponding cut threshold
"""
function extreme_treecut_stats(
    stats_df::AbstractDataFrame,
    perm_aggstats_df::AbstractDataFrame;
    extra_join_cols::Union{Nothing, AbstractVector{Symbol}} = nothing,
    metric_cols::AbstractVector{Symbol} = intersect(TreecutMetrics, propertynames(stats_df)),
    threshold_weight = nothing,
    relative::Bool = false,
    stat_maxquantile::Union{Nothing, Number} = 0.25,
    threshold_range::Union{Tuple{<:Number, <:Number}, Nothing} = nothing
)
    join_cols = [:threshold_bin]
    isnothing(extra_join_cols) || unique!(append!(join_cols, extra_join_cols))

    # bring the values from the 3 consequitive bins in one row
    perm_aggstats_mid_df, perm_aggstats_prev_df, perm_aggstats_next_df = [begin
        df = select(perm_aggstats_df, [[:quantile, :threshold_binmid]; join_cols; metric_cols])
        df.threshold_bin .+= shift
        rename!(df, [col => Symbol(col, suffix) for col in metric_cols])
        rename!(df, :threshold_binmid => Symbol("threshold_bin", suffix))
    end for (suffix, shift) in [("_mid", 0), ("_prev", +1), ("_next", -1)]]
    perm_aggstats3_df = leftjoin(leftjoin(perm_aggstats_mid_df, perm_aggstats_next_df, on=[join_cols; :quantile]),
                                 perm_aggstats_prev_df, on=[join_cols; :quantile])
    # fill missing values on the edges with the closest non-missing ones
    for col in [metric_cols; :threshold_bin]
        perm_aggstats3_df[!, Symbol(col, "_prev")] .= coalesce.(
            perm_aggstats3_df[!, Symbol(col, "_prev")],
            perm_aggstats3_df[!, Symbol(col, "_mid")])
        perm_aggstats3_df[!, Symbol(col, "_next")] .= coalesce.(
            perm_aggstats3_df[!, Symbol(col, "_next")],
            perm_aggstats3_df[!, Symbol(col, "_mid")])
    end
    joinstats_df = leftjoin(select(stats_df, [join_cols; metric_cols; [:threshold]]),
                            perm_aggstats3_df, on=join_cols)
    # linearly interpolate perm metrics at the real threshold
    for col in metric_cols
        midcol = Symbol(col, "_mid")
        prevcol = Symbol(col, "_prev")
        nextcol = Symbol(col, "_next")
        joinstats_df[!, Symbol(col, "_perm")] = [begin
            t = r.threshold
            tmid = r.threshold_bin_mid
            if r.threshold < tmid
                vend = r[prevcol]
                tend = r.threshold_bin_prev
            else
                vend = r[nextcol]
                tend = r.threshold_bin_next
            end
            k = (tmid - t)/(tmid - tend)
            r[midcol] * (1 - k) + vend * k
        end for r in eachrow(joinstats_df)]
    end
    if !isnothing(threshold_range)
        filter!(r -> threshold_range[1] <= r.threshold <= threshold_range[2], joinstats_df)
    end
    by_cols = [:quantile]
    isnothing(extra_join_cols) || unique!(append!(by_cols, extra_join_cols))
    combine(groupby(joinstats_df, by_cols)) do df
        thres_weights = threshold_weight !== nothing ? threshold_weight.(eachrow(df)) : nothing
        reduce(vcat, [begin
            perm_col = Symbol(col, "_perm")
            res = reduce(vcat, [begin
                deltas = [ismissing(r[col]) || ismissing(r[perm_col]) || isnan(r[col]) || isnan(r[perm_col]) || (relative && r[perm_col] == 0) ? naval :
                          (r[col] - r[perm_col]) / (relative ? abs(r[perm_col]) : 1.0)
                          for r in eachrow(df)]
                if thres_weights !== nothing
                    deltas .*= thres_weights
                end
                vals = [("opt", aggfun(deltas)[2])]
                if !isnothing(stat_maxquantile)
                    qrange = (aggfun((0.0, 1.0))[1],
                              aggfun((stat_maxquantile, 1 - stat_maxquantile))[1])
                    if qrange[1] > qrange[2]
                        qrange = (qrange[2], qrange[1])
                    end
                    drange = quantile(deltas, qrange)
                    smallpos = bigpos = 0
                    @inbounds for i in 1:nrow(df)
                        (drange[1] <= deltas[i] <= drange[2]) || continue
                        if (smallpos == 0) || (df.threshold[smallpos] < df.threshold[i])
                            smallpos = i
                        end
                        if (bigpos == 0) || (df.threshold[bigpos] > df.threshold[i])
                            bigpos = i
                        end
                    end
                    (smallpos > 0) && push!(vals, ("small", smallpos))
                    (bigpos > 0) && push!(vals, ("big", bigpos))
                end
                reduce(vcat, [begin
                    DataFrame(
                        :metric => col,
                        :stat => aggtype,
                        :value_type => valtype,
                        :value => df[valpos, col],
                        :permuted_value => df[valpos, perm_col],
                        :delta => df[valpos, col] - df[valpos, perm_col],
                        :threshold => df.threshold[valpos],
                        :threshold_bin => df.threshold_bin[valpos]
                    )
                end for (valtype, valpos) in vals])
            end for (aggtype, (aggfun, naval)) in ["min" => (findmin, Inf), "max" => (findmax, -Inf)]])
            res
        end for col in metric_cols])
    end
end
