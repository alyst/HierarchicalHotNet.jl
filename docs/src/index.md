# HierarchicalHotNet.jl package

*HierarchicalHotNet.jl* implements *Hierarchical HotNet* algorithm (see the original paper by [_M.A. Reyna et al_ (2018)](https://academic.oup.com/bioinformatics/article/34/17/i972/5093236))
with a few permormance optimization and additional features.

## Introduction

Network-based analysis is a powerful bioinformatic tool to put high-throughput experimental data in the context of
current knowledge of cellular interactions and identify potential connections between the perturbed genes.

## Workflow

 * load the network of gene/protein functional interactions, e.g. [ReactomeFI](https://reactome.org/tools/reactome-fiviz) or [STRING](https://string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Homo+sapiens) protein links.
 * use experimental data to assign weights to the genes/proteins in the network
 * generate randomized data by [reshuffling node weights](@ref permweights)
 * apply [*random walk with restart*](@ref netdiff) to get the *stationary distribution* of vertex
   visiting probabilities and *edge transition* probabilities
 * use random walk-based weighted graphs to generate [trees of *Strongly Connected Components*](@ref scctree)
 * analyze the [distribution of SCC metrics](@ref netstats) at each cutting threshold for
   read data-based SCC tree and the randomized ones
 * [cut the SCC tree at the optimal edge threshold](@ref HierarchicalHotNet.cut) and [export](@ref HierarchicalHotNet.export_flowgraph) the result as a collection of data frames
