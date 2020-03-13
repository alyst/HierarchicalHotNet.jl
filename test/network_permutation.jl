@testset "vertexbins()" begin
    g = HHN.hhotnet_example_graph()
    @test_throws ArgumentError HHN.vertexbins(g, nbins=10, by=:bad)

    @testset "vertexbins(g, by=:by)" for by in [:in, :out]
        bins = HHN.vertexbins(g, nbins=10, by=by)
        @test bins isa HHN.IndicesPartition
        @test HHN.nelems(bins) == nv(g)
        @test length(bins) <= 10
    end
end
