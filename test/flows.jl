@testset "componentsflowgraph(Matrix)" begin
    @test_throws DimensionMismatch HHN.componentsflowgraph(fill(0.0, (1, 2)), [[1]], [[1]])
    @test_throws DimensionMismatch HHN.componentsflowgraph(fill(0.0, (2, 2)), [[1]], [[1]])
    @test_throws DimensionMismatch HHN.componentsflowgraph(fill(0.0, (2, 2)), [[1], [2]], [[1]])
    @test_throws DimensionMismatch HHN.componentsflowgraph(fill(0.0, (2, 2)), [[1]], [[1], [2]])
    @test_throws MethodError HHN.componentsflowgraph(fill(0.0, (2, 2)), [[1], [2]], [[1], [2]], HHN.EdgeTest{Int}())

    @testset "empty" begin
        adjmtx = fill(0.0, (0, 0))
        subgraph, flows = @inferred HHN.componentsflowgraph(adjmtx, Vector{Int}[], Vector{Int}[], HHN.EdgeTest{Float64}())
        @test (subgraph, flows) == HHN.componentsflowgraph(adjmtx, HHN.IndicesPartition(), HHN.IndicesPartition())
        @test flows isa Vector{HHN.Flow}
        @test subgraph isa Vector{HHN.Diedge}
    end

    @testset "1x1" begin
        adjmtx = fill(1.0, (1, 1))
        subgraph, flows = HHN.componentsflowgraph(adjmtx, [[]], [[]])
        @test isempty(flows)
        @test isempty(subgraph)

        subgraph, flows = HHN.componentsflowgraph(adjmtx, [[1]], [[]])
        @test isempty(flows)
        @test isempty(subgraph)

        subgraph, flows = HHN.componentsflowgraph(adjmtx, [[]], [[1]])
        @test isempty(flows)
        @test isempty(subgraph)

        subgraph, flows = HHN.componentsflowgraph(adjmtx, [[1]], [[2]])
        @test flows == [(1 => 1, 1)]
        @test subgraph == [1 => 1]

        # self-loops when the matrix doesn't have them
        subgraph, flows = HHN.componentsflowgraph(fill(0.0, (1, 1)), [[1]], [[2]])
        @test flows == [(1 => 1, 1)]
        @test subgraph == [1 => 1]
    end

    @testset "2x2" begin
        subgraph, flows = HHN.componentsflowgraph([1 0; 0 1], [[1], [1]], [[1], []])
        @test flows == [(1 => 1, 1)]
        @test subgraph == [1 => 1]

        subgraph, flows = HHN.componentsflowgraph([1 0; 0 1], [[1], [1]], [[1], [1]])
        @test flows == [(1 => 1, 1), (2 => 2, 1)]
        @test subgraph == [1 => 1, 2 => 2]

        subgraph, flows = HHN.componentsflowgraph([1 1; 0 1], [[1], [1]], [[1], [1]])
        @test sort(flows) == [(1 => 1, 1), (2 => 1, 2), (2 => 2, 1)]
        @test sort(subgraph) == [1 => 1, 2 => 1, 2 => 2]

        # self-loops when there are no edges again
        subgraph, flows = HHN.componentsflowgraph([0 1; 0 0], [[1], [1]], [[1], [1]])
        @test sort(flows) == [(1 => 1, 1), (2 => 1, 2), (2 => 2, 1)]
        @test sort(subgraph) == [1 => 1, 2 => 1, 2 => 2]

        subgraph, flows = HHN.componentsflowgraph([1 1; 0 1], [[1], []], [[], [1]])
        @test isempty(flows)
        @test isempty(subgraph)

        subgraph, flows = HHN.componentsflowgraph([1 0; 1 1], [[1], []], [[], [1]])
        @test flows == [(1 => 2, 2)]
        @test subgraph == [1 => 2]
    end

    @testset "3x3" begin
        subgraph, flows = HHN.componentsflowgraph([0 0 0; 1 0 0; 0 1 0], [[1], [], []], [[], [], [1]])
        @test flows == [(1 => 3, 3)]
        @test sort(subgraph) == [1 => 2, 2 => 3]

        subgraph, flows = HHN.componentsflowgraph([0 0 0; 1 0 0; 0 1 0], [[1], [1], []], [[], [], [1]])
        @test sort(flows) == [(1 => 3, 3), (2 => 3, 2)]
        @test sort(subgraph) == [1 => 2, 2 => 3]

        subgraph, flows = HHN.componentsflowgraph([0 0 0; 1 0 0; 0 1 0], [[], [1], []], [[], [], [1]])
        @test flows == [(2 => 3, 2)]
        @test subgraph == [2 => 3]

        subgraph, flows = HHN.componentsflowgraph([0 0 0; 1 0 0; 1 1 0], [[1], [], []], [[], [], [1]])
        @test sort(flows) == [(1 => 3, 2)]
        @test sort(subgraph) == [1 => 2, 1 => 3, 2 => 3]

        subgraph, flows = HHN.componentsflowgraph([0 0 0; 1 0 1; 0 0 0], [[1], [], [1]], [[], [1], []])
        @test sort(flows) == [(1 => 2, 2), (3 => 2, 2)]
        @test sort(subgraph) == [1 => 2, 3 => 2]

        subgraph, flows = HHN.componentsflowgraph([0 0 0; 1 1 1; 0 0 0], [[1], [1], [1]], [[], [1], []])
        @test sort(flows) == [(1 => 2, 2), (2 => 2, 1), (3 => 2, 2)]
        @test sort(subgraph) == [1 => 2, 2 => 2, 3 => 2]
    end

    @testset "4x4" begin
        subgraph, flows = HHN.componentsflowgraph([0 0 0 0; 1 0 0 0; 0 1 0 0; 0 1 0 0],
                                    [[1], [], [], []],
                                    [[], [], [1], []])
        @test flows == [(1 => 3, 3)]
        @test sort(subgraph) == [1 => 2, 2 => 3]
    end

    @testset "tarjan1983" begin
        graph = HHN.tarjan1983_example_graph()
        etest = HHN.EdgeTest{Float64}(threshold=20.0)
        conncomp = HHN.strongly_connected_components(weights(graph), etest)
        adjmtx = HHN.condense(collect(weights(graph)), conncomp, etest)

        subgraph, flows = HHN.componentsflowgraph(adjmtx,
                                    [[], [1], [], [], []],
                                    [[], [], [], [1], [1]],
                                    etest)
        @test isempty(flows)
        @test isempty(subgraph)

        subgraph, flows = HHN.componentsflowgraph(adjmtx,
                                    [[], [], [1], [], []],
                                    [[1], [], [], [], []],
                                    HHN.EdgeTest{Float64}(threshold=15))
        @test flows == [(3 => 1, 2)]
        @test sort(subgraph) == [2 => 1, 3 => 1, 3 => 2]

        @testset "with ObjectPools" begin
            pools = HHN.ObjectPools()
            subgraph2, flows2 = HHN.componentsflowgraph(adjmtx,
                                        [[], [], [1], [], []],
                                        [[1], [], [], [], []],
                                        HHN.EdgeTest{Float64}(threshold=15),
                                        pools)
            @test flows2 == flows
            @test subgraph2 == subgraph
            @test haskey(pools, Vector{Bool})
            @test pools[Vector{Bool}].nborrowed == 0
            @test haskey(pools, Vector{HHN.DFSState})
            @test pools[Vector{HHN.DFSState}].nborrowed == 0
            @test haskey(pools, Dict{Int,Int})
            @test pools[Dict{Int,Int}].nborrowed == 0
            @test haskey(pools, Vector{Union{Dict{Int,Int},Nothing}})
            @test pools[Vector{Union{Dict{Int,Int},Nothing}}].nborrowed == 0
        end
    end
end

@testset "flowgraph(SCCTree)" begin
    adjmtx = weights(HHN.tarjan1983_example_graph())
    tree = HHN.scctree(adjmtx)
    subgraph, flows, conncomps = HHN.flowgraph(tree, adjmtx, [1, 4, 5, 7], [2, 4, 6],
                                           HHN.EdgeTest{Float64}(threshold=20))
    @test sort(flows) == [(1 => 6, 4 => 5, 2), (4 => 4, 2 => 2, 1), (5 => 2, 3 => 1, 2),
                          (5 => 6, 3 => 5, 3), (7 => 2, 1 => 1, 1), (7 => 6, 1 => 5, 2)]
    @test sort(subgraph) == [(1 => 6, 4 => 5), (2 => 7, 1 => 1), (3 => 2, 1 => 1), (3 => 7, 1 => 1),
                             (5 => 7, 3 => 1), (7 => 3, 1 => 1), (7 => 6, 1 => 5)]
    @test HHN.nflows(tree, adjmtx, [1, 4, 5, 7], [2, 4, 6],
                     HHN.EdgeTest{Float64}(threshold=20)) == (nflows=6, ncompflows=6,
                                                              flowlen_sum=11, compflowlen_sum=11, compflowlen_max=3,
                                                              ncompsources=4, ncompsinks=3)

    @testset "Flows only" begin
        flows2 = Vector{HHN.CompFlow}()
        conncomps2 = HHN.IndicesPartition()
        subgraphres, flowsres, compsres = HHN.flowgraph!(nothing, flows2, conncomps2,
                                                         tree, adjmtx, [1, 4, 5, 7], [2, 4, 6],
                                                         HHN.EdgeTest{Float64}(threshold=20))
        @test isnothing(subgraphres)
        @test flowsres === flows2
        @test flowsres == flows
        @test compsres === conncomps2
        @test compsres == conncomps
    end

    @testset "with ObjectPools" begin
        pools = HHN.ObjectPools()
        subgraph2, flows2, conncomps2 = HHN.flowgraph(tree, adjmtx, [1, 4, 5, 7], [2, 4, 6],
                                                      HHN.EdgeTest{Float64}(threshold=20), pools)
        @test flows2 == flows
        @test subgraph2 == subgraph
        @test conncomps2 == conncomps
        @test haskey(pools, Vector{Bool})
        @test pools[Vector{Bool}].nborrowed == 0
        @test haskey(pools, Vector{HHN.DFSState})
        @test pools[Vector{HHN.DFSState}].nborrowed == 0
        @test haskey(pools, Dict{Int,Int})
        @test pools[Set{Int}].nborrowed == 0
        @test haskey(pools, Vector{Union{Dict{Int,Int},Nothing}})
        @test pools[Vector{Union{Dict{Int,Int},Nothing}}].nborrowed == 0
        @test haskey(pools, HHN.IndicesPartition)
        @test pools[HHN.IndicesPartition].nborrowed == 0
    end
end
