@testset "MatrixOutedgesIterator" begin
    mtx = [0.0 2 0;
           0.1 3 0;
             1 0 0]
    @inferred HHN.outedges(mtx, 1)
    @test_throws BoundsError HHN.outedges(mtx, 0)
    @test_throws BoundsError HHN.outedges(mtx, 4)
    @test collect(HHN.outedges(mtx, 1)) == [2 => 0.1, 3 => 1.0]
    @test collect(HHN.outedges(mtx, 2)) == [1 => 2.0, 2 => 3.0]
    @test collect(HHN.outedges(mtx, 3)) == Pair{Int, Float64}[]
    @test collect(HHN.outedges(mtx, 2, HHN.EdgeTest{Float64}(skipval=nothing))) == [1 => 2.0, 2 => 3.0, 3 => 0.0]
    @test collect(HHN.outedges(mtx, 2, HHN.EdgeTest{Float64}(skipval=3))) == [1 => 2.0, 3 => 0.0]
    @test collect(HHN.outedges(mtx, 2, HHN.EdgeTest{Float64}(threshold=3))) == [2 => 3.0]
    @test collect(HHN.outedges(mtx, 2, HHN.EdgeTest{Float64}(threshold=4))) == Pair{Int, Float64}[]
    @test collect(HHN.outedges(mtx, 2, HHN.EdgeTest{Float64}(threshold=2, rev=true))) == [1 => 2.0]
    @test collect(HHN.outedges(mtx, 2, HHN.EdgeTest{Float64}(threshold=0, rev=true))) == Pair{Int, Float64}[]
    @test collect(HHN.outedges(mtx, 2, HHN.EdgeTest{Float64}(threshold=0, skipval=nothing, rev=true))) == [3 => 0.0]
end
