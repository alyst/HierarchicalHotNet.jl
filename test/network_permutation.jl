@testset "vertexbins()" begin
    g = HHN.hhotnet_example_graph()
    @test_throws ArgumentError HHN.vertexbins(g, nbins=10, by=:bad)
    @test_throws ArgumentError HHN.vertexbins(g, nbins=10, method=:unknown)

    @testset "vertexbins(g, by=:by, method=:method)" for by in [:in, :out, :outXin], method in [:sort, :tree]
        bins = HHN.vertexbins(g, nbins=10, by=by, method=method)
        @test bins isa HHN.IndicesPartition
        @test HHN.nelems(bins) == nv(g)
        @test length(bins) <= 10
    end

    @testset "vertexbins(g, subset, by=:by, method=:method)" for by in [:in, :out, :outXin], method in [:sort, :tree]
        bins = HHN.vertexbins(g, 2:9, nbins=4, by=by, method=method)
        @test bins isa HHN.IndicesPartition
        @test HHN.nelems(bins) == 8
        @test all(v -> in(v, 2:9), bins.elems)
        @test length(bins) <= 4
    end
end

@testset "randpermgroups()" begin
    @test HHN.randpermgroups(Vector{Int}[]) == Int[]
    x = Int[]
    @test HHN.randpermgroups!(x, Vector{Int}[]) === x
    @test isempty(x)
    @test HHN.randpermgroups!(x, HHN.IndicesPartition(0)) === x

    # some tests disabled since groups could be a subset of the indices
    #@test_throws DimensionMismatch HHN.randpermgroups!([1], Vector{Int}[])
    #@test_throws DimensionMismatch HHN.randpermgroups!([1], HHN.IndicesPartition(0))
    @test_throws DimensionMismatch HHN.randpermgroups!([1], [[1], [2]])
    @test_throws DimensionMismatch HHN.randpermgroups!([1], HHN.IndicesPartition(2))
    @test_throws DimensionMismatch HHN.randpermgroups!(Int[], [[1]])
    @test_throws DimensionMismatch HHN.randpermgroups!(Int[], HHN.IndicesPartition(1))
    @test_throws DimensionMismatch HHN.randpermgroups!(Int[], HHN.IndicesPartition(2))

    @test HHN.randpermgroups([[1], [2], [3], [4], [5]]) == 1:5
    @test HHN.randpermgroups(HHN.IndicesPartition(5)) == 1:5

    for i in 1:10
        @test HHN.randpermgroups([[1, 3], [2, 4]]) ∈ [[1, 2, 3, 4], [3, 2, 1, 4],
                                                      [3, 4, 1, 2], [1, 4, 3, 2]]
    end

    # permutations of a subset
    @test HHN.randpermgroups!([1, 2, 3, 4, 5], Vector{Int}[]) == 1:5
    @test HHN.randpermgroups!([1, 2, 3, 4, 5], [[1], [4]]) == 1:5
    for i in 1:10
        @test HHN.randpermgroups!([1, 2, 3, 4, 5], [[1, 3], [2, 5]]) ∈
            [[1, 2, 3, 4, 5], [3, 2, 1, 4, 5],
             [3, 5, 1, 4, 2], [1, 5, 3, 4, 2]]
    end
end
