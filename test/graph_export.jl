@testset "export_graph" begin
    g = HHN.hhotnet_example_graph()
    gmtx = HHN.stepmatrix(g)
    gwalkmtx = HHN.random_walk_matrix(gmtx, 0.2)
    gtree = HHN.scctree(gwalkmtx, verbose=false)
    flowgraph = HHN.export_flowgraph(gtree, 0.05, gwalkmtx, [1, 2], [3, 4], verbose=true, stepmatrix=gmtx)
end