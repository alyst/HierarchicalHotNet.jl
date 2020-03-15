@testset "strongly_connected_components()" begin
    @test_throws DimensionMismatch HHN.strongly_connected_components(reshape([1, 2], (2, 1)))
    @test HHN.strongly_connected_components(reshape([1], (1, 1))) == [[1]]
    @test HHN.strongly_connected_components([1 0; 0 1]) == [[1], [2]]
    @test sort!.(HHN.strongly_connected_components([1 0; 0 1], skipval=nothing, rev=true)) == [[1, 2]]
    @test sort!.(HHN.strongly_connected_components([1 0; 0 1], skipval=nothing)) == [[1, 2]]
    @test sort!.(HHN.strongly_connected_components([1 0; 0 1], skipval=nothing, threshold=0)) == [[1, 2]]
    @test HHN.strongly_connected_components([1 0; 0 1], skipval=nothing, threshold=1) == [[1], [2]]
    @test HHN.strongly_connected_components([1 2; 2 1], skipval=nothing, threshold=1, rev=true) == [[1], [2]]
    @test sort!.(HHN.strongly_connected_components([1 2; 2 1], skipval=nothing, threshold=2, rev=true)) == [[1, 2]]
    @test sort!.(HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights)) == [1:7]
    @test HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, threshold=100) == [[1], [2], [3], [4], [5], [6], [7]]
    @test sort!.(HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, threshold=20)) == [[1], [5], [2, 3, 7], [4], [6]]
    @test sort!.(HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, threshold=10)) == [[1, 2, 3, 4, 5, 7], [6]]

    pool = HHN.ArrayPool{Int}()
    @test sort!.(HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, pool, threshold=10)) == [[1, 2, 3, 4, 5, 7], [6]]
    @test pool.nborrowed == 0

    @testset "strongly_connected_components(TunnelMatrix)" begin
        mtx = weights(HHN.tarjan1983_example_graph())
        tmtx = HHN.TunnelsMatrix(mtx, [2,1,3], [1,4,6,7])
        @test HHN.strongly_connected_components(tmtx) == HHN.strongly_connected_components(collect(tmtx))
    end
end
