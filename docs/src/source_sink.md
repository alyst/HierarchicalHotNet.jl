# [Source-Sink Analysis](@id source_sink)

For any edge weight threshold *t*, cutting the [*strongly connected components (SCC) tree*](@ref scctree)
provides the set of SCC of the original directed weighted network.
Within each SCC subnetwork, there is a path from any node to any other node along
the edges with weights not smaller than *t*. The original network may still contain
the other edges with weights *≥t*, which can connect one SCC subnetwork to the
other, but it is not possible to reenter the same SCC again by traveling along
these edges (otherwise by definition there is a bigger SCC). In other words,
for each *t* we have a *direct acyclic graph* of connections between SCCs.

This property makes it very convenient to enumerate all paths from one set of nodes
(*sources*) to the other (*sinks*). It is implemented by the *flowgraph()* method.
One can use it to identify signaling subnetworks that connect one biological
data (e.g. interactors of a particular protein) to another (e.g. downstream
changes resulting from knock out or overexpressing this protein).

The output of *flowgraph()* is the subnetwork that consists of selected SCCs
and the edges that connect these SCCs, plus the list of paths within this subnetwork
from *source* to *sink* nodes. It could be shown that the path lengths tend
to be smaller for the diffusion networks based on the real data than the ones
based on the reshuffled node weights.

```@docs
HierarchicalHotNet.flowgraph
HierarchicalHotNet.flowgraph!
HierarchicalHotNet.nflows
HierarchicalHotNet.traceflows
HierarchicalHotNet.traceflows!
```
