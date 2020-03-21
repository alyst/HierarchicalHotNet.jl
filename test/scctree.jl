@testset "scctree()" begin
    g0 = SimpleWeightedDiGraph(0)
    g1 = SimpleWeightedDiGraph(1)

    g = HHN.tarjan1983_example_graph()
    adjmtx = g.weights

    # 4-vertex graph with 2 connected components
    g2 = SimpleWeightedDiGraph(4)
    for (u, v, w) in [(1, 2, 10.), (2, 1, 12.), (3, 4, 11.), (4, 3, 12.)]
        add_edge!(g2, u, v, w)
    end

    # 4-vertex graph with 2 connected components, all weights are the same
    g3 = SimpleWeightedDiGraph(4)
    for (u, v, w) in [(1, 2, 10.), (2, 1, 10.), (3, 4, 10.), (4, 3, 10.)]
        add_edge!(g3, u, v, w)
    end

    @test HHN.weighttype(HHN.SCCTree{Float16}) === Float16

    @testset "scctree(method=:$method)" for method in [:bottomup, :bisect]
        g0tree = HHN.scctree(g0, verbose=false, method=method)
        @test HHN.weighttype(g0tree) === eltype(weights(g0))
        @test HHN.nvertices(g0tree) == 0
        @test isempty(g0tree.thresholds)
        @test isempty(HHN.cut(g0tree, 0.0))

        g1tree = HHN.scctree(g1, verbose=false, method=method)
        @test HHN.nvertices(g1tree) == 1
        @test isempty(g1tree.thresholds)
        @test HHN.cut(g1tree, 0.0) == [[1]]

        gtree = HHN.scctree(adjmtx, verbose=false, method=method)
        @test HHN.nvertices(gtree) == 7
        @test HHN.weighttype(gtree) === eltype(adjmtx)
        @test gtree.thresholds == [8.0, 12, 13, 30]
        groot = gtree.nodes[1]
        @test groot.parent == 0
        @test groot.threshold == 8.0
        @test HHN.nvertices(groot) == 7
        @test HHN.cut(gtree, 0.0) == [1:7]
        @test HHN.cut(gtree, 8.0) == [1:7]
        @test HHN.cut(gtree, 10.0) == [[1, 2, 3, 4, 5, 7], [6]]
        @test HHN.cut(gtree, 10.0, minsize=2) == [[1, 2, 3, 4, 5, 7]]
        @test HHN.cut(gtree, 12.0) == [[1, 2, 3, 4, 5, 7], [6]]
        @test HHN.cut(gtree, 13.0) == [[2, 3, 4, 5, 7], [1], [6]]
        @test HHN.cut(gtree, 16.0) == [[2, 3, 7], [4], [5], [1], [6]]
        @test HHN.cut(gtree, 30.0) == [[2, 3, 7], [4], [5], [1], [6]]
        @test HHN.cut(gtree, 31.0) == [[2], [3], [7], [4], [5], [1], [6]]
        @test HHN.cut(gtree, 100.0) == [[2], [3], [7], [4], [5], [1], [6]]
        @test isempty(HHN.cut(gtree, 100.0, minsize=2))

        gtree_rev = HHN.scctree(adjmtx, verbose=false, rev=true, method=method)
        groot_rev = gtree_rev.nodes[1]
        @test gtree_rev.thresholds == [50., 45, 35, 16, 12]
        @test groot_rev.threshold == 50.0
        @test HHN.nvertices(groot_rev) == 7
        @test HHN.cut(gtree_rev, 0.0) ==  [[1], [2], [7], [3], [4], [5], [6]]
        @test HHN.cut(gtree_rev, 12.0) == [[1, 2], [7], [3], [4], [5], [6]]
        @test HHN.cut(gtree_rev, 16.0) == [[1, 2], [7], [3], [4, 5], [6]]
        @test HHN.cut(gtree_rev, 35.0) == [[1, 2, 7], [3], [4, 5], [6]]
        @test HHN.cut(gtree_rev, 44.0) == [[1, 2, 7], [3], [4, 5], [6]]
        @test HHN.cut(gtree_rev, 45.0) == [[1, 2, 3, 7], [4, 5], [6]]
        @test HHN.cut(gtree_rev, 45.0, minsize=2) == [[1, 2, 3, 7], [4, 5]]
        @test HHN.cut(gtree_rev, 45.0, minsize=3) == [[1, 2, 3, 7]]
        @test isempty(HHN.cut(gtree_rev, 45.0, minsize=5))
        @test HHN.cut(gtree_rev, 50.0) == [1:7]
        @test HHN.cut(gtree_rev, 100.0) == [1:7]
        @test isempty(HHN.cut(gtree_rev, 100.0, minsize=8))

        g2tree = HHN.scctree(g2, verbose=false, method=method)
        @test HHN.nvertices(g2tree) == 4
        @test g2tree.thresholds == [10., 11.]
        g2root = g2tree.nodes[1]
        @test g2root.parent == 0
        @test isnan(g2root.threshold) # because multiple connected components
        @test HHN.nvertices(g2root) == 4
        @test HHN.cut(g2tree, 0.0) ==  [1:2, 3:4]
        @test HHN.cut(g2tree, 10.0) == [1:2, 3:4]
        @test HHN.cut(g2tree, 11.0) == [[1], [2], [3,4]]
        @test HHN.cut(g2tree, 12.0) == [[1], [2], [3], [4]]

        g2tree_rev = HHN.scctree(g2, verbose=false, rev=true, method=method)
        g2root_rev = g2tree_rev.nodes[1]
        @test g2tree_rev.thresholds == [12.]
        @test isnan(g2root_rev.threshold)
        @test HHN.nvertices(g2root_rev) == 4
        @test HHN.cut(g2tree_rev, 0.0) ==  [[1], [2], [3], [4]]
        @test HHN.cut(g2tree_rev, 10.0) == [[1], [2], [3], [4]]
        @test HHN.cut(g2tree_rev, 11.0) == [[1], [2], [3], [4]]
        @test HHN.cut(g2tree_rev, 12.0) == [[1, 2], [3, 4]]
        @test HHN.cut(g2tree_rev, 13.0) == [[1, 2], [3, 4]]

        g3tree = HHN.scctree(g3, verbose=false, method=method)
        @test HHN.nvertices(g3tree) == 4
        @test g3tree.thresholds == [10.]
        g3root = g3tree.nodes[1]
        @test g3root.parent == 0
        @test isnan(g3root.threshold) # because multiple connected components
        @test HHN.nvertices(g3root) == 4
        @test HHN.cut(g3tree, 0.0) ==  [[1, 2], [3, 4]]
        @test HHN.cut(g3tree, 10.0) == [[1, 2], [3, 4]]
        @test HHN.cut(g3tree, 11.0) == [[1], [2], [3], [4]]

        g3tree_rev = HHN.scctree(g3, verbose=false, rev=true, method=method)
        g3root_rev = g3tree_rev.nodes[1]
        @test g3tree_rev.thresholds == [10.]
        @test isnan(g3root_rev.threshold)
        @test HHN.nvertices(g3root_rev) == 4
        @test HHN.cut(g3tree_rev, 0.0) ==  [[1], [2], [3], [4]]
        @test HHN.cut(g3tree_rev, 10.0) == [[1, 2], [3, 4]]
        @test HHN.cut(g3tree_rev, 11.0) == [[1, 2], [3, 4]]
    end

    @testset "scctree(TunnelMatrix)" begin
        mtx = weights(HHN.tarjan1983_example_graph())
        tmtx = HHN.TunnelsMatrix(mtx, [2,1,3], [1,4,6,7])
        tmtx_dense = collect(tmtx)

        @testset "scctree(method=:$method)" for method in [:bottomup, :bisect]
            tmtx_tree = HHN.scctree(tmtx, verbose=false, method=method)
            tmtx_dense_tree = HHN.scctree(tmtx_dense, verbose=false, method=method)
            @test tmtx_tree == tmtx_dense_tree
        end
    end
end
