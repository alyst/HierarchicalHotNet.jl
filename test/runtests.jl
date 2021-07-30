using Test, Compat
using HierarchicalHotNet, SimpleWeightedGraphs, LightGraphs
using LinearAlgebra, DataFrames

const HHN = HierarchicalHotNet

include("object_pool.jl")
include("partition.jl")
@testset "utils.jl" begin
    include("utils.jl")
end
include("network_permutation.jl")
include("edges_iterator.jl")
include("connected_components.jl")
include("scctree.jl")
include("flows.jl")
include("stats.jl")
include("graph_export.jl")
