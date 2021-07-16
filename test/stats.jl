@testset "treecut_stats()" begin
    g = HHN.tarjan1983_example_graph()
    adjmtx = g.weights
    adjmtx ./= maximum(vec(adjmtx))
    vweights = fill(0.0, nv(g))
    vweights[[2, 4, 7]] .= 1.0
    walkmtx = HHN.random_walk_matrix(adjmtx, 0.5) * Diagonal(vweights)
    tree = HHN.scctree(walkmtx, verbose=false)
    treestats_df = HHN.treecut_stats(tree, walkmatrix=walkmtx)
    @test treestats_df isa DataFrame
    @test !("nflows" in names(treestats_df))

    @testset "with flow stats" begin
        treestats_ex_df = HHN.treecut_stats(tree, walkmatrix=walkmtx,
                                            sources = [2, 4, 7], sinks = [1, 2, 3])
        @test treestats_ex_df isa DataFrame
        @test "nflows" in names(treestats_ex_df)
        @test size(treestats_df, 1) == size(treestats_ex_df, 1)
    end
end
