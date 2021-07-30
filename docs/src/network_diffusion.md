# [Network Diffusion](@id netdiff)

*Network diffusion* is an efficient way to model signal propagation in molecular networks.
While it could not be regarded as an accurate simulation of biological system, *network
diffusion* allows assessing complex molecular networks and revealing various topological
structures, such as hubs, communities etc.

Currently, *HierarchicalHotNet.jl* package implements [*random walk with restart*](https://en.wikipedia.org/wiki/Random_walk)
diffusion method. The input weighted graph is processed by [`stepmatrix()`](@ref HierarchicalHotNet.stepmatrix) function,
which prepares the *adjacency matrix* for the random walk.
The resulting matrix is then submitted to [`random_walk_with_restart()`](@ref HierarchicalHotNet.random_walk_with_restart) method that generates
the matrix of random walk transition probabilities.

Alternatively, one can use [`stabilized_stepmatrix()`](@ref HierarchicalHotNet.stabilized_stepmatrix) method to weight the
edges of the original network according to how frequently their are being utilized during the random walk.

## [Node and edge weights](@id netweights)

TODO

```@docs
HierarchicalHotNet.stepmatrix
HierarchicalHotNet.random_walk_matrix
HierarchicalHotNet.similarity_matrix
HierarchicalHotNet.stabilized_stepmatrix
HierarchicalHotNet.neighborhood_weights
```

## [Node weights permutation](@id permweights)

One way to confirm the relevance of network diffusion-based predictions is to show that the
results based on real the data differ significantly from the ones based on randomized data, where no
meaningful patterns are expected.

In the original HierarchicalHotNet paper by [_Reina et al (2018)_](https://academic.oup.com/bioinformatics/article/34/17/i972/5093236)
it is proposed to group the vertices with similar in- and out-degrees into bins (see [`vertexbins()`](@ref HierarchicalHotNet.vertexbins))
and randomly shuffle the weights of the vertices within each bin (see [`randpermgroups()`](@ref HierarchicalHotNet.randpermgroups)).

This scheme would reassign the weights of hub vertices to other hubs, and low-degree vertices --
to other low-degree vertices. So, in addition to preserving the overall distribution of node weights,
the correlation of node degree and its weight would be kept too.
However, such reshuffling should eliminate the correlation of paths between the specific nodes and the node weights along these paths.

```@docs
HierarchicalHotNet.vertexbins
HierarchicalHotNet.randpermgroups!
HierarchicalHotNet.randpermgroups
```

