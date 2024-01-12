# [Partitions](@id partition)

`Partition` type provides an efficient container for storing multiple vectors of elements.
The advantage over `Vector{Vector{T}}` is that internally it uses single `Vector{T}`
to store the elements of all of its parts. The disadvantage is that it
doesn't support adding or removing elements to/from arbitrary part,
only the last part could be modified.
`Partition` supports iterator interface for iterating over its parts as
well as parts indexing and `filter!` for elements filtering.

```@docs
HierarchicalHotNet.AbstractPartition
HierarchicalHotNet.Partition
HierarchicalHotNet.PartitionPart
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

# Matrix utilities

```@docs
HierarchicalHotNet.condense
HierarchicalHotNet.condense!
```

# Object pools

Set of utilities for managing object pools to reduce the GC stress.

```@docs
HierarchicalHotNet.ObjectPool
HierarchicalHotNet.ArrayPool
HierarchicalHotNet.NoopObjectPool
HierarchicalHotNet.borrow!
HierarchicalHotNet.release!
HierarchicalHotNet.ObjectPools
```

# Other

```@docs
HierarchicalHotNet.hilbertorder
```
