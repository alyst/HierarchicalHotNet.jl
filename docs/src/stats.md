# [Network Statistics](@id netstats)

[Permutation](@ref permweights) allows generating randomized data that mimcs
the key properties of the original vertex weight distribution.
With enough permutations, it's possible to analyze how different are the results
of network diffusion based on real weights in comparison to permuted weights.

The package allows doing this analysis at the level of individual vertices
([`HierarchicalHotNet.vertex_stats`](@ref)), directed edges ([`HierarchicalHotNet.diedge_stats`](@ref)),
connected components ([`HierarchicalHotNet.conncomponents_stats`](@ref)) etc.

The statistcs from multiple permutation and cutting thresholds could be binned
([`HierarchicalHotNet.bin_treecut_stats`](@ref)) and then aggregated for calculating
the quantiles of resulting distributions ([`HierarchicalHotNet.aggregate_treecut_binstats`](@ref)).
Finally, [`HierarchicalHotNet.extreme_treecut_stats`](@ref) can find the edge cutting threshold
with the maximal difference between the real and permuted weights.

```@docs
HierarchicalHotNet.vertex_stats
HierarchicalHotNet.diedge_stats
HierarchicalHotNet.conncomponents_stats
HierarchicalHotNet.treecut_stats
HierarchicalHotNet.treecut_compstats
HierarchicalHotNet.bin_treecut_stats
HierarchicalHotNet.aggregate_treecut_binstats
HierarchicalHotNet.extreme_treecut_stats
HierarchicalHotNet.TreecutMetrics
```
