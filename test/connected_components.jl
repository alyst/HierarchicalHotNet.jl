@testset "strongly_connected_components()" begin
    @test_throws DimensionMismatch HHN.strongly_connected_components(reshape([1, 2], (2, 1)))
    @test HHN.strongly_connected_components(reshape([1], (1, 1))) == [[1]]
    @test HHN.strongly_connected_components([1 0; 0 1]) == [[1], [2]]
    @test sort!.(HHN.strongly_connected_components([1 0; 0 1], HHN.EdgeTest{Int}(skipval=nothing, rev=true))) == [[1, 2]]
    @test sort!.(HHN.strongly_connected_components([1 0; 0 1], HHN.EdgeTest{Int}(skipval=nothing))) == [[1, 2]]
    @test sort!.(HHN.strongly_connected_components([1 0; 0 1], HHN.EdgeTest{Int}(skipval=nothing, threshold=0))) == [[1, 2]]
    @test HHN.strongly_connected_components([1 0; 0 1], HHN.EdgeTest{Int}(skipval=nothing, threshold=1)) == [[1], [2]]
    @test HHN.strongly_connected_components([1 2; 2 1], HHN.EdgeTest{Int}(skipval=nothing, threshold=1, rev=true)) == [[1], [2]]
    @test sort!.(HHN.strongly_connected_components([1 2; 2 1], HHN.EdgeTest{Int}(skipval=nothing, threshold=2, rev=true))) == [[1, 2]]
    @test sort!.(HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights)) == [1:7]
    @test HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, HHN.EdgeTest{Float64}(threshold=100)) == [[1], [2], [3], [4], [5], [6], [7]]
    @test sort!.(HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, HHN.EdgeTest{Float64}(threshold=20))) == [[1], [5], [2, 3, 7], [4], [6]]
    @test sort!.(HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, HHN.EdgeTest{Float64}(threshold=10))) == [[1, 2, 3, 4, 5, 7], [6]]

    pools = HHN.ObjectPools()
    @test sort!.(HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, HHN.EdgeTest{Float64}(threshold=10), pools)) == [[1, 2, 3, 4, 5, 7], [6]]
    @test haskey(pools, Vector{Int})
    @test pools[Vector{Int}].nborrowed == 0
    @test haskey(pools, Vector{HHN.DFSState})
    @test pools[Vector{HHN.DFSState}].nborrowed == 0

    # test in-place version
    ptn = HHN.IndicesPartition(3)
    @test HHN.strongly_connected_components!(ptn, HHN.tarjan1983_example_graph().weights, HHN.EdgeTest{Float64}(threshold=10), pools) === ptn
    @test sort!.(ptn) == [[1, 2, 3, 4, 5, 7], [6]]

    @testset "strongly_connected_components(TunnelMatrix)" begin
        mtx = weights(HHN.tarjan1983_example_graph())
        tmtx = HHN.TunnelsMatrix(mtx, [2,1,3], [1,4,6,7])
        @test HHN.strongly_connected_components(tmtx) == HHN.strongly_connected_components(collect(tmtx))
    end
end
