@testset "export_graph" begin
    g = HHN.hhotnet_example_graph()
    gmtx = HHN.stepmatrix(g)
    gwalkmtx = HHN.random_walk_matrix(gmtx, 0.2)
    gtree = HHN.scctree(gwalkmtx, verbose=false)
    flowgraph = HHN.export_flowgraph(gtree, 0.05, gwalkmtx, [1, 2], [3, 4], verbose=true, stepmatrix=gmtx)
    @test !any(flowgraph.edges.has_trace)
    @test !hasproperty(flowgraph.edges, :flowpaths)
    flowgraph2 = HHN.export_flowgraph(gtree, 0.05, gwalkmtx, [1, 2], [3, 4], verbose=true, stepmatrix=gmtx, flowpaths=:stepedges)
    @test_skip any(flowgraph2.edges.has_trace) # doesn't have traces
    @test !hasproperty(flowgraph2.edges, :flowpaths)
    @test !hasproperty(flowgraph2.edges, :flowpaths_rev)
    flowgraph3 = HHN.export_flowgraph(gtree, 0.05, gwalkmtx, [1, 2], [3, 4], verbose=true, stepmatrix=gmtx, flowpaths=:flowattr)
    @test !any(flowgraph3.edges.has_trace)
    @test hasproperty(flowgraph3.edges, :flowpaths)
    @test hasproperty(flowgraph3.edges, :flowpaths_rev)
end