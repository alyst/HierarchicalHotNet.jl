using Test
using HierarchicalHotNet, SimpleWeightedGraphs, LightGraphs

const HHN = HierarchicalHotNet

include("array_pool.jl")
include("partition.jl")
@testset "utils.jl" begin
    include("utils.jl")
end
include("network_permutation.jl")
include("edges_iterator.jl")
include("tunnels_matrix.jl")
include("connected_components.jl")
include("scctree.jl")
