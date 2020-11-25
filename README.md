# HierarchicalHotNet

Julia implementation of *Hierarchical HotNet* algorithm for finding hotspots in
the altered networks.

| Docs | Build | Test |
|:-----|:------|:-----|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://alyst.github.io/HierarchicalHotNet.jl/dev) | [![Build Status](https://travis-ci.org/alyst/HierarchicalHotNet.jl.svg)](https://travis-ci.org/alyst/HierarchicalHotNet.jl) | [![codecov](http://codecov.io/github/alyst/HierarchicalHotNet.jl/branch/master/graph/badge.svg)](http://codecov.io/github/alyst/HierarchicalHotNet.jl) |

The original *Python* implementation by [RaphaelLab](http://compbio.cs.brown.edu/) is available at [hierarchical-hotnet](https://github.com/raphael-group/hierarchical-hotnet) GitHub repository.

In comparison to the original implementation, this Julia package includes a
few enhancements and additional functionality:
* optimizations and enhancements
  * reduced memory footprint
  * edge weights indexing for faster building of SCC tree
  * taking into account node *in-* and *out-degree* when reshuffling the node weights
    for randomized input data
* [*source-sink analysis*](https://alyst.github.io/HierarchicalHotNet.jl/dev/sourcesink.html) for identifying paths from the set of *source* nodes to
  the *sink* nodes
* calculation of edge *p-*values using randomized networks
* export the network as data frames of vertices, edges and connected components

## References
The *Hierarchical HotNet* paper:

> M.A. Reyna, M.D.M. Leiserson, B.J. Raphael. Hierarchical HotNet: identifying hierarchies of altered subnetworks. [_ECCB/Bioinformatics_ **34**(17):i972-980](https://academic.oup.com/bioinformatics/article/34/17/i972/5093236), 2018.
