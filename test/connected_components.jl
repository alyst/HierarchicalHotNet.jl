@testset "strongly_connected_components()" begin
    @test_throws DimensionMismatch HHN.strongly_connected_components(reshape([1, 2], (2, 1)))
    @test HHN.strongly_connected_components(reshape([1], (1, 1))) == [[1]]
    @test HHN.strongly_connected_components([1 0; 0 1]) == [[1], [2]]
    @test HHN.strongly_connected_components([1 0; 0 1], skipval=nothing, rev=true) == [[1, 2]]
    @test HHN.strongly_connected_components([1 0; 0 1], skipval=nothing) == [[1, 2]]
    @test HHN.strongly_connected_components([1 0; 0 1], skipval=nothing, threshold=0) == [[1, 2]]
    @test HHN.strongly_connected_components([1 0; 0 1], skipval=nothing, threshold=1) == [[1], [2]]
    @test HHN.strongly_connected_components([1 2; 2 1], skipval=nothing, threshold=1, rev=true) == [[1], [2]]
    @test HHN.strongly_connected_components([1 2; 2 1], skipval=nothing, threshold=2, rev=true) == [[1, 2]]
    @test HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights) == [1:7]
    @test HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, threshold=100) == [[1], [2], [3], [4], [5], [6], [7]]
    @test HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, threshold=20) == [[1], [5], [2, 3, 7], [4], [6]]
    @test HHN.strongly_connected_components(HHN.tarjan1983_example_graph().weights, threshold=10) == [[1, 2, 3, 4, 5, 7], [6]]
end
