@testset "TunnelsMatrix" begin
    @testset "empty" begin
        tmtx = @inferred HHN.TunnelsMatrix(Matrix{Int}(undef, 0, 0), Int[], Int[])
        @test eltype(tmtx) === Int
        @test HHN.parenttype(tmtx) === Matrix{Int}
        @test HHN.parenttype(typeof(tmtx)) === Matrix{Int}
        @test size(tmtx) == (0, 0)
        @test size(tmtx, 1) == 0
        @test size(tmtx, 2) == 0
        @test isempty(tmtx)
        @test length(tmtx) == 0
        @test axes(tmtx, 1) == Int[]
        @test axes(tmtx, 2) == Int[]
        @test collect(tmtx) == Matrix{Int}(undef, 0, 0)
        @test HHN.TunnelsMatrix(Matrix{Float64}(undef, 0, 0), Int[], Int[],
                                tunnel_weight=5).tunnel_weight === 5.0
    end

    mtx = [1 2.
           3 4
           5 6
           7 8]

    @testset "no tunnels" begin
        tmtx = @inferred HHN.TunnelsMatrix(mtx, Int[], Int[])
        @test tmtx.tunnel_weight == 8.0
        @test size(tmtx) == (2, 2)
        @test size(tmtx, 1) == 2
        @test size(tmtx, 2) == 2
        @test !isempty(tmtx)
        @test length(tmtx) == 4
        @test axes(tmtx, 1) == 1:2
        @test axes(tmtx, 2) == 1:2
        @test collect(tmtx) == view(mtx, 1:2, 1:2)
    end

    @testset "single tunnel from 4 to 1" begin
        tmtx = @inferred HHN.TunnelsMatrix(mtx, [4], [1])
        @test tmtx.tunnel_weight == 8.0
        @test size(tmtx) == (4, 4)
        @test !isempty(tmtx)
        @test length(tmtx) == 16
        @test axes(tmtx, 1) == 1:4
        @test axes(tmtx, 2) == 1:4
        @test collect(tmtx) ==
            [1 2 0 0.
             3 4 0 3
             7 8 0 0
             0 0 8 0]
    end

    @testset "tunnels from [3,2] to 2" begin
        tmtx = @inferred HHN.TunnelsMatrix(mtx, [4,2], [2])
        @test tmtx.tunnel_weight == 8.0
        @test size(tmtx) == (6, 6)
        @test collect(tmtx) ==
            [1 2 0 0 2 2.
             3 4 0 0 0 0
             7 8 0 0 0 0
             3 0 0 0 0 0
             0 0 8 0 0 0
             0 0 0 8 0 0]
    end

    @testset "tunnels from [3,2] to [1,2] with a view" begin
        tmtx = @inferred HHN.TunnelsMatrix(view(mtx, 2:3, 1:2), [2,1], [1,2])
        @test HHN.parenttype(tmtx) === typeof(view(mtx, 2:3, 1:2))
        @test HHN.parenttype(typeof(tmtx)) === typeof(view(mtx, 2:3, 1:2))
        @test tmtx.tunnel_weight == 6.0
        @test size(tmtx) == (8, 8)
        @test collect(tmtx) ==
            [3 4 0 0 0 4 0 4
             5 6 0 0 5 0 5 0
             5 0 0 0 0 0 0 0
             0 4 0 0 0 0 0 0
             0 0 6 0 0 0 0 0
             0 0 6 0 0 0 0 0
             0 0 0 6 0 0 0 0
             0 0 0 6 0 0 0 0]
    end
end

@testset "TunnelsMatrixOutedgesIterator" begin
    @testset "simple matrix" begin
        mtx = [1 2.
               3 4
               5 6]

        tmtx = HHN.TunnelsMatrix(view(mtx, 2:3, 1:2), [2,1], [1,2])
        @inferred HHN.outedges(tmtx, 1)
        tit = HHN.outedges(tmtx, 1)
        @test tit isa HHN.TunnelsMatrixOutedgesIterator
        @test_throws BoundsError HHN.outedges(tmtx, 0)
        @test_throws BoundsError HHN.outedges(tmtx, 9)

        # test tunnels iterator matches dense iterator
        tmtx_dense = collect(tmtx)
        for i in axes(tmtx, 2)
            @test collect(HHN.outedges(tmtx, i)) == collect(HHN.outedges(tmtx_dense, i))
        end
    end

    @testset "tarjan1983" begin
        mtx = weights(HHN.tarjan1983_example_graph())
        tmtx = HHN.TunnelsMatrix(mtx, [2,1,3], [1,4,6,7])

        tmtx_dense = collect(tmtx)
        for i in axes(tmtx, 2)
            @test collect(HHN.outedges(tmtx, i)) == collect(HHN.outedges(tmtx_dense, i))
        end
    end
end

@testset "indexvalues()" begin
        mtx = [1 2.
               3 4
               5 6]
        tmtx = HHN.TunnelsMatrix(view(mtx, 2:3, 1:2), [2,1], [1,2])
        imtx, weights = HHN.indexvalues(Int, tmtx)
        @test imtx isa HHN.TunnelsMatrix{Int}
        @test imtx.entries == tmtx.entries
        @test imtx.tunnels == tmtx.tunnels
        @test imtx.tunnel_weight == 4
        @test imtx.parent == [1 2; 3 4]

        mtx[1, 2] = 4
        tmtx2 = HHN.TunnelsMatrix(view(mtx, 1:2, 1:2), [2,2], [1,2])
        HHN.indexvalues!(imtx, weights, tmtx2)
        @test imtx.entries == tmtx2.entries
        @test imtx.tunnels == tmtx2.tunnels
        @test imtx.tunnel_weight == 3
        @test imtx.parent == [1 3; 2 3]
end

@testset "condense()" begin
    @testset "simple matrix" begin
        mtx = [1 2.
               3 4
               5 6]

        tmtx = HHN.TunnelsMatrix(view(mtx, 2:3, 1:2), [2,1], [1,2])
        @test HHN.condense(tmtx, HHN.IndicesPartition(8, ngroups=8)) == tmtx
        @test HHN.condense(tmtx, HHN.IndicesPartition(8, ngroups=1)) == fill(6., (1, 1))
        @test HHN.condense(tmtx, [1:2, 3:8]) == HHN.condense(collect(tmtx), [1:2, 3:8])
    end

    @testset "tarjan1983" begin
        mtx = weights(HHN.tarjan1983_example_graph())
        tmtx = HHN.TunnelsMatrix(mtx, [2,1,3], [1,4,6,7])

        ptn = [1:2, 3:8, 9:11, 12:15, 16:21, 22:22]
        @test HHN.condense(tmtx, ptn) == HHN.condense(collect(tmtx), ptn)
    end
end

@testset "subgraph_adjacencymatrix()" begin
    @testset "simple matrix" begin
        mtx = [1 2.
               3 4
               5 6]

        tmtx = HHN.TunnelsMatrix(view(mtx, 2:3, 1:2), [2,1], [1,2])
        tview = @inferred HHN.subgraph_adjacencymatrix(tmtx, 1:2)
        @test tview isa HHN.TunnelsMatrix
        @test tview == view(mtx, 2:3, 1:2)
        @test HHN.subgraph_adjacencymatrix(tmtx, 1:3) == HHN.subgraph_adjacencymatrix(collect(tmtx), 1:3)
        @test_throws ArgumentError HHN.subgraph_adjacencymatrix(tmtx, [1,2,4,5])
        @test HHN.subgraph_adjacencymatrix(tmtx, [1,2,3,5]) == HHN.subgraph_adjacencymatrix(collect(tmtx), [1,2,3,5])
        @test_throws ArgumentError HHN.subgraph_adjacencymatrix(tmtx, [1,2,4,5,7])
        @test HHN.subgraph_adjacencymatrix(tmtx, [1,2,3,4,5,7]) == HHN.subgraph_adjacencymatrix(collect(tmtx), [1,2,3,4,5,7])

        @testset "using ArrayPool" begin
            pool = HHN.ArrayPool{Int}()
            tview2 = @inferred HHN.subgraph_adjacencymatrix(tmtx, [1, 2, 3, 5], pool=pool)
            @test pool.nborrowed > 0
            HHN.release!(pool, tview2)
            @test pool.nborrowed == 0 # it returned all its borrowed
        end
    end

    @testset "tarjan1983" begin
        mtx = weights(HHN.tarjan1983_example_graph())
        tmtx = HHN.TunnelsMatrix(mtx, [2,1,3], [1,4,6,7])

        @test HHN.subgraph_adjacencymatrix(tmtx, 1:10) == HHN.subgraph_adjacencymatrix(collect(tmtx), 1:10)
        @test HHN.subgraph_adjacencymatrix(tmtx, [1:10; 15:20]) == HHN.subgraph_adjacencymatrix(collect(tmtx), [1:10; 15:20])
    end
end
