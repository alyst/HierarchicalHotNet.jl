var documenterSearchIndex = {"docs":
[{"location":"export.html#export","page":"Export","title":"Export HotNet Networks","text":"","category":"section"},{"location":"export.html","page":"Export","title":"Export","text":"HierarchicalHotNet.export_flowgraph","category":"page"},{"location":"export.html#HierarchicalHotNet.export_flowgraph","page":"Export","title":"HierarchicalHotNet.export_flowgraph","text":"export_flowgraph(tree::SCCTree{T}, threshold::Number,\n                 walkmatrix::AbstractMatrix,\n                 sources::AbstractVector{<:Integer}, sinks::AbstractVector{<:Integer};\n                 orig_diedges::Union{AbstractDataFrame, Nothing} = nothing,\n                 vertices_stats::Union{AbstractDataFrame, Nothing} = nothing,\n                 diedges_stats::Union{AbstractDataFrame, Nothing} = nothing,\n                 flowpaths::Symbol = :skip,\n                 stepmatrix::Union{AbstractMatrix, Nothing} = nothing,\n                 step_threshold::Number = 0.75 * threshold, maxsteps::Integer = 2,\n                 step_sinks::Union{AbstractVector, AbstractSet, Nothing} = nothing,\n                 step_sources::Union{AbstractVector, AbstractSet, Nothing} = nothing,\n                 flow_edges::Bool=false,\n                 pvalue_mw_max::Number=0.05,\n                 pvalue_fisher_max::Number=0.05,\n                 verbose::Bool=false,\n                 pools::Union{ObjectPools, Nothing}=nothing,\n                 mincompsize::Union{Integer, Nothing}=nothing,\n                 exported_sinks::AbstractVector{<:Integer}=sinks\n) -> NamedTuple\n\nCuts the tree at threshold and exports the resulting SCC network as the collection of dataframes.\n\nIf specified, calculates the flows from sources to sinks and returns them as additional edges.\n\nKeyword arguments\n\norig_diedges::AbstractDataFrame: optional collection of the original edges. The metadata from this frame is added to the overlapping diedges of the SCC network.\nvertices_stats::AbstractDataFrame: optional vertices statistics\ndiedges_stats::AbstarctDataFrame: optional directed edges statistics\nflowpaths::Symbol: how the flows should be traced\nskip (default): no tracing\nflowattr: trace the flows (see HierarchicalHotNet.traceflows) and add as flowpaths column to diedges data frame\nsteps: trace the flows (see HierarchicalHotNet.traceflows) and add as extra diedges of type step to diedges data frame\n\nReturns\n\nNamed tuple with fields\n\ncomponents::DataFrame: statistics for Strongly Connected Components\nvertices::DataFrame: network vertices\ndiedges::DataFrame: directed edges\nedges::DataFrame: undirected edges\n\nSee also\n\nHierarchicalHotNet.cut\n\n\n\n\n\n","category":"function"},{"location":"source_sink.html#source_sink","page":"Source → Sink Flows","title":"Source-Sink Analysis","text":"","category":"section"},{"location":"source_sink.html","page":"Source → Sink Flows","title":"Source → Sink Flows","text":"For any edge weight threshold t, cutting the strongly connected components (SCC) tree provides the set of SCC of the original directed weighted network. Within each SCC subnetwork, there is a path from any node to any other node along the edges with weights not smaller than t. The original network may still contain the other edges with weights ≥t, which can connect one SCC subnetwork to the other, but it is not possible to reenter the same SCC again by traveling along these edges (otherwise by definition there is a bigger SCC). In other words, for each t we have a direct acyclic graph of connections between SCCs.","category":"page"},{"location":"source_sink.html","page":"Source → Sink Flows","title":"Source → Sink Flows","text":"This property makes it very convenient to enumerate all paths from one set of nodes (sources) to the other (sinks). It is implemented by the HierarchicalHotNet.flowgraph() method. One can use it to identify signaling subnetworks that connect one biological data (e.g. interactors of a particular protein) to another (e.g. downstream changes resulting from knock out or overexpressing this protein).","category":"page"},{"location":"source_sink.html","page":"Source → Sink Flows","title":"Source → Sink Flows","text":"The output of HierarchicalHotNet.flowgraph() is the subnetwork that consists of selected SCCs and the edges that connect these SCCs, plus the list of paths within this subnetwork from source to sink nodes. It could be shown that the path lengths tend to be smaller for the diffusion networks based on the real data than the ones based on the reshuffled node weights.","category":"page"},{"location":"source_sink.html","page":"Source → Sink Flows","title":"Source → Sink Flows","text":"HierarchicalHotNet.flowgraph\nHierarchicalHotNet.flowgraph!\nHierarchicalHotNet.nflows\nHierarchicalHotNet.traceflows\nHierarchicalHotNet.traceflows!","category":"page"},{"location":"source_sink.html#HierarchicalHotNet.nflows","page":"Source → Sink Flows","title":"HierarchicalHotNet.nflows","text":"nflows(comps::IndicesPartition, adjmtx::AbstractMatrix,\n       sources::AbstractVector{Int}, sinks::AbstractVector{Int},\n       test::EdgeTest;\n       maxweight::Union{Number, Nothing} = nothing,\n       used_sources::Union{AbstractVector{Int}, Nothing}=nothing,\n       used_sinks::Union{AbstractVector{Int}, Nothing}=nothing) -> NamedTupe\n\nCalculates statistics from the flows from sources to sinks vertices in the weighted directed graph defined by adjmtx and test and the comps vertex components.\n\nReturns the NamedTuple with the following fields\n\nnflows: the number of source → sink flows\nncompoflows: the number of unique component(source) → component(sink) pairs for each source → sink flow\nflowlen_sum: total length of all flows (flow length = the number of SCCs it crosses)\ncompflowlen_sum: total length of flows (every unique pairs of source/sink components is counted once)\nflowinvlen_sum: the sum of 1mathrmflowlength\ncompflowinvlen_sum: the sum of 1mathrmflowlength (each unique source/sink component pair is counted once)\ncompflowlen_max:\nfloweight_sum:\nfloweight_sum:\ncompfloweight_sum:\nflowavghopweight_sum:\nncompsources: the number of distinct SCCs that have source nodes\nncompsinks: the number of distinct SCCs that have sink nodes\n\nKeyword arguments\n\nmaxweight: if specified, limits the flow weight by maxweight\nused_sources::AbstractVector{Int}: if specified, acts as an output parameter that contains the  sorted list of sources that have outgoing flows\nused_sinks::AbstractVector{Int}: if specified, acts as an output parameter that contains the  sorted list of sinks that have incoming flows\n\n\n\n\n\n","category":"function"},{"location":"source_sink.html#HierarchicalHotNet.traceflows","page":"Source → Sink Flows","title":"HierarchicalHotNet.traceflows","text":"traceflows(step_adjmtx::AbstractMatrix, steptest::EdgeTest,\n           walk_adjmtx::AbstractMatrix, walktest::EdgeTest;\n           kwargs...) -> Dict{Diedge, Partition{Int}}\n\nTrace the random walk (specified by walk_adjmtx and walktest) steps in the original graph (given by step_adjmtx and steptest).\n\nSee HierarchicalHotNet.traceflows! for the detailed description.\n\n\n\n\n\n","category":"function"},{"location":"source_sink.html#HierarchicalHotNet.traceflows!","page":"Source → Sink Flows","title":"HierarchicalHotNet.traceflows!","text":"traceflows!(flow2paths::AbstractDict{Diedge, Partition{Int}},\n            step_adjmtx::AbstractMatrix,\n            steptest::EdgeTest,\n            walk_adjmtx::AbstractMatrix,\n            walktest::EdgeTest;\n            sources::Union{AbstractVector, AbstractSet, Nothing} = nothing,\n            sinks::Union{AbstractVector, AbstractSet, Nothing} = nothing,\n            maxsteps::Integer=2) -> Dict{Diedge, Partition{Int}}\n\nTrace the random walk (specified by walk_adjmtx and walktest) steps in the original graph (given by step_adjmtx and steptest).\n\nKeyword arguments\n\nsources (optional): the indices of vertices to use as path starts. If not specified, all vertices are used as path starts.\nsinks (optional): the indices of vertices to use as path ends. If not specified, all vertices are used as path ends.\nmaxsteps=2: maximal number of steps in a traced path, longer paths are discarded.\n\nReturns the mapping from the flow diedges to the HierarchicalHotNet.Parition object. Each part corresponds to the path, from diedge start to diedge end, in the original network, and the part elements are the indices of the intermedidate vertices along the path (start and end not included).\n\n\n\n\n\n","category":"function"},{"location":"network_diffusion.html#netdiff","page":"Network Diffusion","title":"Network Diffusion","text":"","category":"section"},{"location":"network_diffusion.html","page":"Network Diffusion","title":"Network Diffusion","text":"Network diffusion is an efficient way to model signal propagation in molecular networks. While it could not be regarded as an accurate simulation of biological system, network diffusion allows assessing complex molecular networks and revealing various topological structures, such as hubs, communities etc.","category":"page"},{"location":"network_diffusion.html","page":"Network Diffusion","title":"Network Diffusion","text":"Currently, HierarchicalHotNet.jl package implements random walk with restart diffusion method. The input weighted graph is processed by `stepmatrix function, which prepares for random walk, and then submitted to `randomwalkwith_restart.","category":"page"},{"location":"network_diffusion.html#netweights","page":"Network Diffusion","title":"Node and edge weights","text":"","category":"section"},{"location":"network_diffusion.html","page":"Network Diffusion","title":"Network Diffusion","text":"HierarchicalHotNet.stepmatrix\nHierarchicalHotNet.random_walk_matrix\nHierarchicalHotNet.similarity_matrix\nHierarchicalHotNet.neighborhood_weights","category":"page"},{"location":"network_diffusion.html#HierarchicalHotNet.stepmatrix","page":"Network Diffusion","title":"HierarchicalHotNet.stepmatrix","text":"stepmatrix(g::Union{AbstractSimpleWeightedGraph, AbstractMatrix}\n           inedge_weights::Union{AbstractVector, Nothing} = nothing,\n           normalize_weights::Bool = true) -> AbstractMatrix\n\nConstruct the step matrix for the random walk.\n\n\n\n\n\n","category":"function"},{"location":"network_diffusion.html#HierarchicalHotNet.random_walk_matrix","page":"Network Diffusion","title":"HierarchicalHotNet.random_walk_matrix","text":"Find the matrix that transforms given initial node probabilities into stationary distribution of visiting probabilities of a random walk with restart.\n\n\n\n\n\n","category":"function"},{"location":"network_diffusion.html#permweights","page":"Network Diffusion","title":"Node weights permutation","text":"","category":"section"},{"location":"network_diffusion.html","page":"Network Diffusion","title":"Network Diffusion","text":"One way to confirm the relevance of network diffusion-based predictions is to show that the results based on real the data differ significantly from the ones based on randomized data, where no meaningful patterns are expected.","category":"page"},{"location":"network_diffusion.html","page":"Network Diffusion","title":"Network Diffusion","text":"In the original HierarchicalHotNet paper by Reina et al (2018) it is proposed to group the vertices with similar in- and out-degrees into bins (see vertexbins) and randomly shuffle the weights of the vertices within each bin (see randpermgroups).","category":"page"},{"location":"network_diffusion.html","page":"Network Diffusion","title":"Network Diffusion","text":"This scheme would reassign the weights of hub vertices to other hubs, and low-degree vertices – to other low-degree vertices. So, in addition to preserving the overall distribution of node weights, the correlation of node degree and its weight would be kept too. However, such reshuffling should eliminate the correlation of paths between the specific nodes and the node weights along these paths.","category":"page"},{"location":"network_diffusion.html","page":"Network Diffusion","title":"Network Diffusion","text":"HierarchicalHotNet.vertexbins\nHierarchicalHotNet.randpermgroups!\nHierarchicalHotNet.randpermgroups","category":"page"},{"location":"network_diffusion.html#HierarchicalHotNet.vertexbins","page":"Network Diffusion","title":"HierarchicalHotNet.vertexbins","text":"vertexbins(g::AbstractSimpleWeightedGraph,\n           vertices::AbstractArray{Int}=1:nv(g);\n           nbins::Integer=10, by::Symbol=:out, method=:tree) ->\n    IndicesPartition\n\nPartition the vertices of the weighted graph g into nbins bins, so that the vertices of the same bin have similar sum of edge weights.\n\nKeyword arguments\n\nnbins::Integer: the number of vertex bins. Use nbins &ge; 1 for  bigger (nv &ge; 1) networks\nby::Symbol: what edges to consider for grouping vertices. The supported options are:\n:out (default, as in the original paper) – use only the outgoing edges\n:in – use only the incoming edges\n:outXin (recommended) – use all edges, so that the vertices of the same bin have similar sums of both incoming and outgoing edges\nmethod::Symbol: the method for binning vertices\n:sort (as in the original paper) – order the vertices by the sum of weights, then split the sorted array into equal bins. This method doesn't handle the by=:outXin well.\n:tree (default) – build the hierarchical tree of vertices based on their sum of weights, then cut the tree to get the requested number of bins.\n\n\n\n\n\n","category":"function"},{"location":"network_diffusion.html#HierarchicalHotNet.randpermgroups!","page":"Network Diffusion","title":"HierarchicalHotNet.randpermgroups!","text":"randpermgroups!(v::AbstractVector{Int}, groups::AbstractPartition) -> v\n\nRandomly reshuffle the elements of groups within each group and put the result into v, group by group.\n\n\n\n\n\n","category":"function"},{"location":"network_diffusion.html#HierarchicalHotNet.randpermgroups","page":"Network Diffusion","title":"HierarchicalHotNet.randpermgroups","text":"randpermgroups(groups::AbstractPartition) -> Vector{Int}\n\nRandomly reshuffle the elements of groups within each group.\n\n\n\n\n\n","category":"function"},{"location":"stats.html#netstats","page":"Statistics","title":"Network Statistics","text":"","category":"section"},{"location":"stats.html","page":"Statistics","title":"Statistics","text":"Permutation allows generating randomized data that mimcs the key properties of the original vertex weight distribution. With enough permutations, it's possible to analyze how different are the results of network diffusion based on real weights in comparison to permuted weights.","category":"page"},{"location":"stats.html","page":"Statistics","title":"Statistics","text":"The package allows doing this analysis at the level of individual vertices (HierarchicalHotNet.vertex_stats), directed edges (HierarchicalHotNet.diedge_stats), connected components (HierarchicalHotNet.conncomponents_stats) etc.","category":"page"},{"location":"stats.html","page":"Statistics","title":"Statistics","text":"The statistcs from multiple permutation and cutting thresholds could be binned  (HierarchicalHotNet.bin_treecut_stats) and then aggregated for calculating the quantiles of resulting distributions (HierarchicalHotNet.aggregate_treecut_binstats). Finally, HierarchicalHotNet.extreme_treecut_stats can find the edge cutting threshold with the maximal difference between the real and permuted weights.","category":"page"},{"location":"stats.html","page":"Statistics","title":"Statistics","text":"HierarchicalHotNet.vertex_stats\nHierarchicalHotNet.diedge_stats\nHierarchicalHotNet.conncomponents_stats\nHierarchicalHotNet.treecut_stats\nHierarchicalHotNet.treecut_compstats\nHierarchicalHotNet.bin_treecut_stats\nHierarchicalHotNet.aggregate_treecut_binstats\nHierarchicalHotNet.extreme_treecut_stats","category":"page"},{"location":"stats.html#HierarchicalHotNet.vertex_stats","page":"Statistics","title":"HierarchicalHotNet.vertex_stats","text":"vertex_stats(weights::AbstractVector{<:Number},\n             walkweights::AbstractVector{<:Number},\n             permweights::AbstractMatrix{<:Number},\n             walkpermweights::AbstractMatrix{<:Number}) -> DataFrame\n\nCalculates statistics for the permuted vertex weights distribution and how it is different from the actual weights.\n\nweights: weights of the vertices in the original network\nwalkweights: weights of the vertices after network diffusion analysis (stationary random walk distribution)\npermweights: matrix of permuted weights; rows correspond to vertices, columns – to permutations\nwalkpermweights: matrix of vertex weights based on network diffusion analysis using permweights as input; rows correspond to vertices, columns to permutations\n\n\n\n\n\n","category":"function"},{"location":"stats.html#HierarchicalHotNet.diedge_stats","page":"Statistics","title":"HierarchicalHotNet.diedge_stats","text":"diedge_stats(weights::AbstractVector{<:Number},\n             walkweights::AbstractVector{<:Number},\n             permweights::AbstractMatrix{<:Number},\n             walkpermweights::AbstractMatrix{<:Number}) -> DataFrame\n\nCalculates statistics for the directed edges permuted weights distribution and how it is different from the actual weights of directed edges.\n\n\n\n\n\n","category":"function"},{"location":"index.html#HierarchicalHotNet.jl-package","page":"Introduction","title":"HierarchicalHotNet.jl package","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"HierarchicalHotNet.jl implements Hierarchical HotNet algorithm (see the original paper by M.A. Reyna et al (2018)) with a few permormance optimization and additional features.","category":"page"},{"location":"index.html#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Network-based analysis is a powerful bioinformatic tool to put high-throughput experimental data in the context of current knowledge of cellular interactions and identify potential connections between the perturbed genes.","category":"page"},{"location":"index.html#Workflow","page":"Introduction","title":"Workflow","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"load the network of gene/protein functional interactions, e.g. ReactomeFI or STRING protein links.\nuse experimental data to assign weights to the genes/proteins in the network\ngenerate randomized data by reshuffling node weights\napply random walk with restart to get the stationary distribution of vertex visiting probabilities and edge transition probabilities\nuse random walk-based weighted graphs to generate trees of Strongly Connected Components\nanalyze the distribution of SCC metrics at each cutting threshold for read data-based SCC tree and the randomized ones\ncut the SCC tree at the optimal edge threshold and export the result as a collection of data frames","category":"page"},{"location":"scctree.html#scctree","page":"SCC Tree","title":"Strongly Connected Components Tree","text":"","category":"section"},{"location":"scctree.html","page":"SCC Tree","title":"SCC Tree","text":"HierarchicalHotNet.SCCTree\nHierarchicalHotNet.scctree\nHierarchicalHotNet.cut","category":"page"},{"location":"scctree.html#HierarchicalHotNet.SCCTree","page":"SCC Tree","title":"HierarchicalHotNet.SCCTree","text":"The tree of strongly connected components. Organizes the strongly connected components of the weighted directed graph into a tree.\n\nCutting the tree at specific threshold with cut gives the corresponding strongly connected components. The root of the tree corresponds to the weakest threshold.\n\n\n\n\n\n","category":"type"},{"location":"scctree.html#HierarchicalHotNet.scctree","page":"SCC Tree","title":"HierarchicalHotNet.scctree","text":"scctree(g::Union{AbstractSimpleWeightedGraph, AbstractMatrix};\n        method=:bisect, skipval=0, rev=false) -> SCCTree\n\nComputes the hierarchical decomposition of the weighted directed graph g into strongly connected components.\n\nSupports weighted graph objects as well as their adjacency matrix representations.\n\nKeyword arguments\n\nmethod::Symbol (defaults to :bisect): the method for partitioning.  The supported methods are :bisect (the fastest) and :bottomup (slow,  but simpler).\nskipval::Number (defaults to zero): what value of the adjacency matrix  should be treated as \"no edge\".\nrev::Bool (defaults to false): if true, bigger edges weights are considered  weaker than smaller ones\n\n\n\n\n\n","category":"function"},{"location":"scctree.html#HierarchicalHotNet.cut","page":"SCC Tree","title":"HierarchicalHotNet.cut","text":"cut(tree::SCCTree, threshold; minsize=1) -> IndicesPartition\n\nCuts the tree at the given edge weight threshold to get the corresponding strongly connected components of the original graph.\n\nKeyword arguments\n\nminsize: the minimal number of vertices in the component. Smaller connected components are skipped. By default returns components of any size.\n\n\n\n\n\n","category":"function"},{"location":"utils.html#partition","page":"Utilities","title":"Partitions","text":"","category":"section"},{"location":"utils.html","page":"Utilities","title":"Utilities","text":"Partition type provides an efficient container for storing multiple vectors of elements. The advantage over Vector{Vector{T}} is that internally it uses single Vector{T} to store the elements of all of its parts. The disadvantage is that it doesn't support adding or removing elements to/from arbitrary part, only the last part could be modified. Partition supports iterator interface for iterating over its parts as well as parts indexing and filter! for elements filtering.","category":"page"},{"location":"utils.html","page":"Utilities","title":"Utilities","text":"HierarchicalHotNet.AbstractPartition\nHierarchicalHotNet.Partition\nHierarchicalHotNet.PartitionPart","category":"page"},{"location":"utils.html#HierarchicalHotNet.AbstractPartition","page":"Utilities","title":"HierarchicalHotNet.AbstractPartition","text":"Base class for partitions of elements of type T.\n\n\n\n\n\n","category":"type"},{"location":"utils.html#HierarchicalHotNet.Partition","page":"Utilities","title":"HierarchicalHotNet.Partition","text":"Efficiently stores the partition of a vector of elements of type T into disjoint sets.\n\nSupports iterator interface.\n\n\n\n\n\n","category":"type"},{"location":"utils.html#HierarchicalHotNet.PartitionPart","page":"Utilities","title":"HierarchicalHotNet.PartitionPart","text":"The type of parts returned by Partition{T}\n\n\n\n\n\n","category":"type"},{"location":"utils.html","page":"Utilities","title":"Utilities","text":"HierarchicalHotNet.nelems\nHierarchicalHotNet.elems\n\nHierarchicalHotNet.nparts\nHierarchicalHotNet.partrange\nHierarchicalHotNet.partlength\nHierarchicalHotNet.ispartempty\n\nHierarchicalHotNet.pushelem!\nHierarchicalHotNet.closepart!\n\nHierarchicalHotNet.repeat!\nHierarchicalHotNet.IndicesPartition\nHierarchicalHotNet.reset!","category":"page"},{"location":"utils.html#HierarchicalHotNet.nelems","page":"Utilities","title":"HierarchicalHotNet.nelems","text":"nelems(ptn::AbstractPartition)\n\nTotal number of elements in the partition.\n\n\n\n\n\n","category":"function"},{"location":"utils.html#HierarchicalHotNet.elems","page":"Utilities","title":"HierarchicalHotNet.elems","text":"elems(ptn::Partition{T}) where T -> Vector{T}\n\nVector of elements in the partition.\n\nRerturns the vectors as they are stored in ptn: the elements are ordered by their parts (elements from the first part come first etc) and maintain the same order as within each part.\n\n\n\n\n\n","category":"function"},{"location":"utils.html#HierarchicalHotNet.nparts","page":"Utilities","title":"HierarchicalHotNet.nparts","text":"nelems(ptn::AbstractPartition) -> Int\n\nNumber of parts in the partition.\n\nSame as length(ptn).\n\n\n\n\n\n","category":"function"},{"location":"utils.html#HierarchicalHotNet.pushelem!","page":"Utilities","title":"HierarchicalHotNet.pushelem!","text":"pushelem!(ptn::Partition, el) -> ptn\n\nPush the element into the last part of partition.\n\n\n\n\n\n","category":"function"},{"location":"utils.html#HierarchicalHotNet.closepart!","page":"Utilities","title":"HierarchicalHotNet.closepart!","text":"closepart!(ptn::Partition) -> ptn\n\nFinalize the last part of the partition and append the new empty part.\n\nAll subsequent calls to pushelem! will add elements into the new part.\n\n\n\n\n\n","category":"function"},{"location":"utils.html#HierarchicalHotNet.repeat!","page":"Utilities","title":"HierarchicalHotNet.repeat!","text":"repeat!(ptn::Partition, n::Integer) -> ptn\n\nRepeat the parts of ptn n times.\n\njulia> using HierarchicalHotNet\n\njulia> ptn = HierarchicalHotNet.IndicesPartition(5, nparts=1)\n1-element HierarchicalHotNet.Partition{Int64}:\n [1, 2, 3, 4, 5]\n\njulia> HierarchicalHotNet.repeat!(ptn, 3)\n3-element HierarchicalHotNet.Partition{Int64}:\n [1, 2, 3, 4, 5]\n [1, 2, 3, 4, 5]\n [1, 2, 3, 4, 5]\n\n\n\n\n\n","category":"function"},{"location":"utils.html#HierarchicalHotNet.IndicesPartition","page":"Utilities","title":"HierarchicalHotNet.IndicesPartition","text":"Partition of an integer vector\n\n\n\n\n\n","category":"type"},{"location":"utils.html#HierarchicalHotNet.reset!","page":"Utilities","title":"HierarchicalHotNet.reset!","text":"reset!(p::IndicesPartition, n::Integer=nelems(p); nparts::Integer=n) -> p\n\nResets the integer partition p.\n\nIf nparts == 1, the partition is reset into [[1, 2, 3, ..., n]]. If n == nparts, sets the partition to [[1], [2], [3], ..., [n]].\n\n\n\n\n\n","category":"function"}]
}
