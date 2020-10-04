# [Network diffusion](@id netdiff)

```@docs
HierarchicalHotNet.stepmatrix
HierarchicalHotNet.random_walk_matrix
HierarchicalHotNet.similarity_matrix
HierarchicalHotNet.neighborhood_weights
```

## Network Permutation

One way to confirm the relevance of predictions based on the network diffusion
is to compare them with the predictions based on the randomized data.
One scheme of generating realistic randomized dataset, proposed in the original paper,
is to reshuffle the weights of the vertices of similar degrees, so that hubs weights
are randomly assigned to other hubs, and low-degree vertices to other low-degree vertices.

```@docs
HierarchicalHotNet.vertexbins
HierarchicalHotNet.randpermgroups!
HierarchicalHotNet.randpermgroups
```
