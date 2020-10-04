# [Partitions](@id partition)

`Partition` is an efficient way to store multiple vectors of elements.
The advantage over `Vector{Vector{T}}` is that internally it uses the single
container to store the elements of all its parts. The disadvantage is that it
doesn't support adding or removing elements to/from arbitrary part,
only the last part could be modified.
`Partition` supports iterator interface for iterating over its parts as
well as parts indexing and `filter!` for elements filtering.

```@docs
HierarchicalHotNet.AbstractPartition
HierarchicalHotNet.PartitionPart
HierarchicalHotNet.Partition
```

```@docs
HierarchicalHotNet.nelems
HierarchicalHotNet.elems

HierarchicalHotNet.nparts
HierarchicalHotNet.partrange
HierarchicalHotNet.partlength
HierarchicalHotNet.ispartempty

HierarchicalHotNet.pushelem!
HierarchicalHotNet.closepart!

HierarchicalHotNet.repeat!
HierarchicalHotNet.IndicesPartition
HierarchicalHotNet.reset!
```