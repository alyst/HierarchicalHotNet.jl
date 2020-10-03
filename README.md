# HierarchicalHotNet

Julia implementation of *Hierarchical HotNet* algorithm for finding hotspots in
the altered networks.

| Docs | Build | Test |
|:-----|:------|:-----|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://alyst.github.io/HierarchicalHotNet.jl/dev) | [![Build Status](https://travis-ci.org/alyst/HierarchicalHotNet.jl.svg)](https://travis-ci.org/alyst/HierarchicalHotNet.jl) | [![codecov](http://codecov.io/github/alyst/HierarchicalHotNet.jl/branch/master/graph/badge.svg)](http://codecov.io/github/alyst/HierarchicalHotNet.jl) |

The original *Python* implementation by [RaphaelLab](http://compbio.cs.brown.edu/) is available at [hierarchical-hotnet](https://github.com/raphael-group/hierarchical-hotnet) GitHub repository.

In comparison to the original implementation, this Julia package includes a
few optimizations and additional functions:
* reduced memory footprint
* edge weights indexing for faster building of SCC tree
* taking into account node *in-* and *out-degree* when reshuffling the node weights
  for randomized input data
* *source-sink analysis* for identifying paths from the set of *source* nodes to
  the *sink* nodes (see below)

## Source-Sink Analysis

For any edge weight threshold *t*, the *strongly connected components tree*
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

## References
The *Hierarchical HotNet* paper:

> M.A. Reyna, M.D.M. Leiserson, B.J. Raphael. Hierarchical HotNet: identifying hierarchies of altered subnetworks. [_ECCB/Bioinformatics_ **34**(17):i972-980](https://academic.oup.com/bioinformatics/article/34/17/i972/5093236), 2018.
