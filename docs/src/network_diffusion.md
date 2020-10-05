# [Network diffusion](@id netdiff)

```@docs
HierarchicalHotNet.stepmatrix
HierarchicalHotNet.random_walk_matrix
HierarchicalHotNet.similarity_matrix
HierarchicalHotNet.neighborhood_weights
```

## Network Permutation

One way to confirm the relevance of network diffusion-based predictions is to show that the
results based on real the data differ significantly from the ones based on randomized data, where no
meaningful patterns are expected.
One scheme for generating realistic randomized dataset, proposed in the original paper,
is to group the vertices with similar in- and out-degrees into bins and randomly shuffle
the weights of the vertices within each bin. In this scheme, the weights of hubs
are randomly reassigned to other hubs, and low-degree vertices to other low-degree vertices,
so that overall the distribution of weights across the nodes looks very similar to real
experiments, except that any correlations between the vertex weights and biological pathways are distrupted.

```@docs
HierarchicalHotNet.vertexbins
HierarchicalHotNet.randpermgroups!
HierarchicalHotNet.randpermgroups
```
