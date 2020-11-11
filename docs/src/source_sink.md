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
(*sources*) to the other (*sinks*). It is implemented by the [`HierarchicalHotNet.flowgraph()`](@ref) method.
One can use it to identify signaling subnetworks that connect one biological
data (e.g. interactors of a particular protein) to another (e.g. downstream
changes resulting from knock out or overexpressing this protein).

The output of [`HierarchicalHotNet.flowgraph()`](@ref) is the subnetwork that consists of selected SCCs
and the edges that connect these SCCs, plus the list of paths within this subnetwork
from *source* to *sink* nodes. It could be shown that the path lengths tend
to be smaller for the diffusion networks based on the real data than the ones
based on the reshuffled node weights.

## [Network Metrics](@id network_metrics)

To show that, at specific threshold ``t_*``, *HHotNet* results based
on real data, ``H_{\mathrm{real}}(t_*)``, are significantly
different than random ones, ``H_{\mathrm{perm}, i}(t_*)``, ``i=1,2,\ldots,N_{\mathrm{perm}}``,
the approach of the original publication was to calculate some metric ``m(H)`` and show that
the probability that ``m(H_{\mathrm{real}}(t_*))`` is a sample from ``M_{\mathrm{perm}}(t_*)`` is low:
```math
\min \Big(P\big(m(H_{\mathrm{real}}(t_*)) \leq M_{\mathrm{perm}}(t_*) \big), P\big(m(H_{\mathrm{real}}(t_*)) \geq M_{\mathrm{perm}}(t_*) \big) \Big) \leq p_*,
```
where ``M_{\mathrm{perm}}(t)`` is a random variable derived from the empirical distribution
of ``m(H_{\mathrm{perm}, i}(t))``, ``i = 1, 2, \ldots, N_{\mathrm{perm}}``, and ``p_*`` is the significance level.

In [_M.A. Reyna et al_ (2018)](https://academic.oup.com/bioinformatics/article/34/17/i972/5093236),
the ``m_{\max\mathrm{size}}(H)`` metric -- the size of the maximal strongly connected component -- was used,
and it was shown that ``m_{\max\mathrm{size}}(H_{\mathrm{real}})`` is statistically significantly larger
than at random for some ``t_*``. [HierarchicalHotNet.treecomp_stats](@ref) provides ``m_{\max\mathrm{size}}(t)``
in `maxcomponent_size` column.
It is also possible to use the total size of top *N* components (`topn_components_sizesum` column) or the number of
non-trivial (larger than single node) SCCs (`ncomponents_nontrivial` column) for the same purpose.

### [Flow Metrics](@id sourcesink_metrics)

The analysis of source-to-sink flows requires a different metric, which takes into account the flows that exist
at given threshold *t*. The flow metric that is showing good is *the average inverse of flow length*:
```math
L_{\mathrm{avg}}^{-1}(H) = \frac{1}{N_{\mathrm{source}} \cdot N_{\mathrm{sink}}} \sum_{f \in \mathrm{flows}(H)} \frac{1}{N_{\mathrm{SCC}}(f)},
```
where ``N_{\mathrm{source}}`` and ``N_{\mathrm{sink}}`` are the numbers of sources and sinks, respectively
(so ``N_{\mathrm{source}} \cdot N_{\mathrm{sink}}`` is the number of all possible distinct flows),
``\mathrm{flows}(H)`` are all the source-to-sink flows of ``H``, and ``N_{\mathrm{SCC}}(f)`` is the
number of distinct strongly connected components that contain the nodes of the ``f`` flow.
This metric changes from 1 (all sources and sinks in the same stongly connected component) to 0
(no or infinitely long flows).
This metric is available as `flow_avginvlen` column in the [HierarchicalHotNet.treecomp_stats](@ref) output.

The alternative metrics for the flow analysis provided by [HierarchicalHotNet.treecomp_stats](@ref) are:
 * the average transition probability of the flow (`flow_avgweight` column):
```math
w_{\mathrm{avg. flow}}(H) = \frac{1}{N_{\mathrm{src}} \cdot N_{\mathrm{sink}}} \sum_{f \in \mathrm{flows}(H)} \min \big( w(\mathrm{source}(f), \mathrm{sink}(f)), w_{\max} \big),
```
where ``w(i, j)`` is the transition probability from ``i``-th to ``j``-th node in the random walk with restart.
The distribution of ``w(i, j)`` is essentially non-normal, and to limit the strong influence of the few "outliers" on the
total ``w_{\mathrm{avg. flow}}(H)``, some ``w_{\max}`` constant is used.
* the average flow transition probability per SCC (`flow_avghopweight` column):
```math
w_{\mathrm{avg. hop}}(H) = \frac{1}{N_{\mathrm{src}} \cdot N_{\mathrm{sink}}} \sum_{f \in \mathrm{flows}(H)} \frac{\min \big( w(\mathrm{source}(f), \mathrm{sink}(f)), w_{\max} \big)}{N_{\mathrm{SCC}}(f)},
```

```@docs
HierarchicalHotNet.flowgraph
HierarchicalHotNet.flowgraph!
HierarchicalHotNet.nflows
HierarchicalHotNet.traceflows
HierarchicalHotNet.traceflows!
```
