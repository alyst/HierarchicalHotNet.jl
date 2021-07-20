@testset "isweaker/isstronger()" begin
    @test HHN.isweaker(5, 6, rev=false)
    @test !HHN.isweaker(5, 6, rev=true)
    @test HHN.isweaker(5, 6, rev=false)
    @test !HHN.isweaker(5, 5, rev=true)

    @test !HHN.isweaker(6, 5, rev=false)
    @test HHN.isweaker(6, 5, rev=true)
    @test !HHN.isweaker(6, 5, rev=false)
    @test !HHN.isweaker(5, 5, rev=false)

    @test !HHN.isstronger(5, 6, rev=false)
    @test HHN.isstronger(5, 6, rev=true)
    @test !HHN.isstronger(5, 6, rev=false)
    @test !HHN.isstronger(5, 5, rev=true)

    @test HHN.isstronger(6, 5, rev=false)
    @test !HHN.isstronger(6, 5, rev=true)
    @test HHN.isstronger(6, 5, rev=false)
    @test !HHN.isstronger(5, 5, rev=false)
end

@testset "valuegroups" begin
    @test HHN.valuegroups(5, [2, 1, 4, 5, 3]) == [[2], [1], [5], [3], [4]]
    @test HHN.valuegroups(4, [2, 2, 4, 3, 3]) == [Int[], [1, 2], [4, 5], [3]]

    @test HHN.valuegroups([2, 1, 4, 5, 3]) == [[2], [1], [5], [3], [4]]
    @test HHN.valuegroups([2, 2, 4, 3, 3]) == [Int[], [1, 2], [4, 5], [3]]

    @test HHN.valuegroupsizes(5, [2, 1, 4, 5, 3]) == [1, 1, 1, 1, 1]
    @test HHN.valuegroupsizes(4, [2, 2, 4, 3, 3]) == [0, 2, 2, 1]

    @test HHN.valuegroupsizes([2, 1, 4, 5, 3]) == [1, 1, 1, 1, 1]
    @test HHN.valuegroupsizes([2, 2, 4, 3, 3]) == [0, 2, 2, 1]
end

@testset "sortedvalues()" begin
    adjmtx = HHN.tarjan1983_example_graph().weights

    vals_empty = HHN.sortedvalues(Int[])
    @test issorted(vals_empty)

    vals = HHN.sortedvalues(adjmtx)
    @test issorted(vals)

    vals20 = HHN.sortedvalues(adjmtx, HHN.EdgeTest{Float64}(threshold=20))
    @test issorted(vals20)
    @test all(x -> x <= 20, vals20)

    valsm1 = HHN.sortedvalues(adjmtx, HHN.EdgeTest{Float64}(threshold=-1))
    @test isempty(valsm1)

    rvals = HHN.sortedvalues(adjmtx, HHN.EdgeTest{Float64}(rev=true))
    @test issorted(rvals, rev=true)

    rvals20 = HHN.sortedvalues(adjmtx, HHN.EdgeTest{Float64}(rev=true, threshold=20))
    @test issorted(rvals20, rev=true)
    @test all(x -> x >= 20, rvals20)
end

@testset "indexvalues()" begin
    @test HHN.indexvalues(Int, fill(0.0, (0, 0))) == (fill(0, (0, 0)), Float64[])
    @test HHN.indexvalues(Int, fill(2.0, (1, 1))) == (fill(1, (1, 1)), [2.0])
    mtx1 = [0.0 3.0; 2.0 0.0]
    @test HHN.indexvalues(Int, mtx1) == ([0 2; 1 0], [2.0, 3.0])
    @test HHN.indexvalues(Int, mtx1, HHN.EdgeTest{Float64}(rev=true)) == ([0 1; 2 0], [3.0, 2.0])
    @test HHN.indexvalues(Int, mtx1, HHN.EdgeTest{Float64}(skipval=nothing)) == ([1 3; 2 1], [0.0, 2.0, 3.0])
    @test HHN.indexvalues(Int, mtx1, HHN.EdgeTest{Float64}(threshold=2.5)) == ([0 0; 1 0], [2.0])
    @test HHN.indexvalues(Int, mtx1, HHN.EdgeTest{Float64}(threshold=2.0)) == ([0 0; 1 0], [2.0])
    @test HHN.indexvalues(Int, mtx1, HHN.EdgeTest{Float64}(threshold=1.5)) == ([0 0; 0 0], Float64[])
    @test HHN.indexvalues(Int, mtx1, HHN.EdgeTest{Float64}(rev=true, threshold=2.5)) == ([0 1; 0 0], [3.0])
    @test HHN.indexvalues(Int, mtx1, HHN.EdgeTest{Float64}(rev=true, threshold=3.0)) == ([0 1; 0 0], [3.0])
end

@testset "condensed()" begin
    @test_throws DimensionMismatch HHN.condense(reshape([1], (1, 1)), [[1], [2]])
    @test_throws DimensionMismatch HHN.condense(reshape([1, 2], (1, 2)), [[1]])
    @test HHN.condense(reshape([1, 2], (1, 2)), [[1]], [[1], [2]]) == reshape([1, 2], (1, 2))
    @test HHN.condense(reshape([1], (1, 1)), [[1]]) == reshape([1], (1, 1))
    A = [0. 2. 3.;
         1. 0. 2.;
         2. 4. 0.]
    @test HHN.condense(A, [[1], [2, 3]]) == [0. 3; 2 4]
    @test HHN.condense(A, [[1], [2, 3]], HHN.EdgeTest{Float64}(skipval=3)) == [0. 2; 2 4]
    @test HHN.condense(A, [[1], [3, 2]], HHN.EdgeTest{Float64}(rev=true)) == [0. 2; 1 2]
    @test HHN.condense(A, [[1], [2, 3]], HHN.EdgeTest{Float64}(skipval=3, rev=true)) == [0. 2; 1 0]
    @test HHN.condense(A, [[1], [2, 3]], [[1, 2, 3]], HHN.EdgeTest{Float64}(skipval=0, rev=true)) == reshape([2., 1], (2, 1))
    @test HHN.condense(A, [[1], [3, 2]], [[3, 1, 2]]) == reshape([3., 4], (2, 1))
    @test HHN.condense(A, [[1], [2, 3]], [[2, 1], [3]]) == [2. 3; 4 2]
    @test HHN.condense(A, [[1], [2, 3]], [[1, 2], [3]], HHN.EdgeTest{Float64}(skipval=-1, rev=true)) == [0. 3; 0 0]
end

@testset "insertsorted()" begin
    @test HHN.insertsorted!(Vector{Int}(), 0) == [0]
    @test HHN.insertsorted!([1, 2], 0) == [0, 1, 2]
    @test HHN.insertsorted!([1, 2], 3) == [1, 2, 3]
    @test HHN.insertsorted!([1, 3], 2) == [1, 2, 3]
    @test HHN.insertsorted!([1, 2], 2) == [1, 2, 2]
    @test HHN.insertsorted!([1, 2], 2, unique=true) == [1, 2]
    @test HHN.insertsorted!([1, 2, 3], 2, unique=true) == [1, 2, 3]
    @test HHN.insertsorted!([1, 1], 1) == [1, 1, 1]
    @test HHN.insertsorted!([1, 1, 2, 2], 1) == [1, 1, 1, 2, 2]
end

@testset "hilbertorder()" begin
    @test_throws ArgumentError HHN.hilbertorder(1, 1, 0)
    @test_throws ArgumentError HHN.hilbertorder(1, 1, -1)
    @test_throws ArgumentError HHN.hilbertorder(0, 1, 1)
    @test_throws ArgumentError HHN.hilbertorder(1, 0, 1)

    @test HHN.hilbertorder(1, 1, 1) == 1
    @test HHN.hilbertorder.([1, 2, 1, 2], [1, 1, 2, 2], 2) == [1, 2, 4, 3]
    ci4x4 = CartesianIndices((4, 4))
    @test HHN.hilbertorder.(getindex.(ci4x4, 1), getindex.(ci4x4, 2), 4) ==
           [1 2 15 16; 4 3 14 13; 5 8 9 12; 6 7 10 11]
    ci16x16 = CartesianIndices((16, 16))
    @test length(unique!(vec(HHN.hilbertorder.(getindex.(ci16x16, 1), getindex.(ci16x16, 2), 16)))) == 256
end

@testset "check_square()" begin
    @test HHN.check_square(Matrix{Int}(undef, 0, 0))
    @test HHN.check_square(Matrix{Int}(undef, 5, 5))
    @test_throws DimensionMismatch HHN.check_square(Matrix{Int}(undef, 2, 5))
end