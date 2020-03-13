using Test
using HierarchicalHotNet, SimpleWeightedGraphs, LightGraphs

const HHN = HierarchicalHotNet

include("partition.jl")
@testset "utils.jl" begin
    include("utils.jl")
end
include("connected_components.jl")
include("scctree.jl")
