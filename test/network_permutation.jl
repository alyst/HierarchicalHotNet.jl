@testset "vertexbins()" begin
    g = HHN.hhotnet_example_graph()

    bins = HHN.vertexbins(g, nbins=10)
    @test bins isa HHN.IndicesPartition
    @test HHN.nelems(bins) == nv(g)
    @test length(bins) <= 10
end
