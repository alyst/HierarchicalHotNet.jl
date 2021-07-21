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
    @test size(treestats_df, 1) == HHN.nthresholds(tree)
    @test treestats_df.threshold == tree.thresholds
    @test treestats_df.ncomponents == [5, 6]
    @test treestats_df.ncomponents_nontrivial == [1, 1]
    @test treestats_df.maxcomponent_size == [3, 2]
    @test treestats_df.topn_components_sizesum == [7, 6]

    @test !("nflows" in names(treestats_df))

    @testset "with flow stats" begin
        treestats_ex_df = HHN.treecut_stats(tree, walkmatrix=walkmtx,
                                            sources = [2, 4, 7], sinks = [1, 2, 3])
        @test treestats_ex_df isa DataFrame
        @test "nflows" in names(treestats_ex_df)
        @test treestats_ex_df.nflows == [9, 6]
        @test treestats_ex_df.ncompflows == [3, 4]
        @test treestats_ex_df.flow_avglen == [2/3, 5/6]
        @test treestats_ex_df.flow_avginvlen ≈ [2/3, 23/54]
        @test size(treestats_ex_df, 1) == size(treestats_df, 1)
        @test treestats_ex_df.threshold == treestats_df.threshold
        @test all(isnan, treestats_ex_df.flow_avgweight)
        @test all(isnan, treestats_ex_df.flow_avghopweight)
        @test treestats_ex_df.ncompsources == [1, 2]
        @test treestats_ex_df.topn_nsources == [3, 3]
        @test treestats_ex_df.ncompsinks == [3, 3]
        @test treestats_ex_df.topn_nsinks == [3, 3]

        @testset "with source-sink weights" begin
            treestats_ex_df = HHN.treecut_stats(tree, walkmatrix=walkmtx, sourcesinkweights=walkmtx,
                                                sources = [2, 4, 7], sinks = [1, 2, 3])
            @test treestats_ex_df isa DataFrame
            @test size(treestats_ex_df, 1) == size(treestats_df, 1)
            @test all(!isnan, treestats_ex_df.flow_avgweight)
            @test all(!isnan, treestats_ex_df.flow_avghopweight)
            @test treestats_ex_df.flow_avgweight ≈ [0.14953506817016599, 0.13783585159418063]
            @test treestats_ex_df.flow_avghopweight ≈ [0.11877205756069342, 0.10962421163253283]
            @test treestats_ex_df.compflow_avgweight ≈ [0.29045570649885494, 0.16933645546777068]
            @test treestats_ex_df.flow_distance == [2/3, 26/9]
        end
    end
end
