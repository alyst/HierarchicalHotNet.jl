@testset "vertexbins()" begin
    g = HHN.hhotnet_example_graph()
    @test_throws ArgumentError HHN.vertexbins(g, nbins=10, by=:bad)
    @test_throws ArgumentError HHN.vertexbins(g, nbins=10, method=:unknown)

    @testset "vertexbins(g, by=:by, method=:method)" for by in [:in, :out, :outXin], method in [:sort, :tree]
        bins = HHN.vertexbins(g, nbins=10, by=by, method=method)
        @test bins isa HHN.IndicesPartition
        @test HHN.nelems(bins) == nv(g)
        @test length(bins) <= 10
    end
end
