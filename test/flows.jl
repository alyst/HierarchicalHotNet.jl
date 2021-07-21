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
        @test flows == [(1 => 1, HHN.FlowInfo(0, 1.0))]
        @test subgraph == [1 => 1]

        # self-loops when the matrix doesn't have them
        subgraph, flows = HHN.componentsflowgraph(fill(0.0, (1, 1)), [[1]], [[2]])
        @test flows == [(1 => 1, HHN.FlowInfo(0, 0.0))]
        @test subgraph == [1 => 1]
    end

    @testset "2x2" begin
        subgraph, flows = HHN.componentsflowgraph([1 0; 0 1], [[1], [1]], [[1], []])
        @test flows == [(1 => 1, HHN.FlowInfo(0, 1.0))]
        @test subgraph == [1 => 1]

        subgraph, flows = HHN.componentsflowgraph([1 0; 0 1], [[1], [1]], [[1], [1]])
        @test flows == [(1 => 1, HHN.FlowInfo(0, 1.0)), (2 => 2, HHN.FlowInfo(0, 1.0))]
        @test subgraph == [1 => 1, 2 => 2]

        subgraph, flows = HHN.componentsflowgraph([1 0.8; 0 1], [[1], [1]], [[1], [1]])
        @test sort(flows) == [(1 => 1, HHN.FlowInfo(0, 1.0)),
                              (2 => 1, HHN.FlowInfo(1, 0.8)),
                              (2 => 2, HHN.FlowInfo(0, 1.0))]
        @test sort(subgraph) == [1 => 1, 2 => 1, 2 => 2]

        # self-loops when there are no edges again
        subgraph, flows = HHN.componentsflowgraph([0 0.6; 0 0], [[1], [1]], [[1], [1]])
        @test sort(flows) == [(1 => 1, HHN.FlowInfo(0, 0.0)),
                              (2 => 1, HHN.FlowInfo(1, 0.6)),
                              (2 => 2, HHN.FlowInfo(0, 0.0))]
        @test sort(subgraph) == [1 => 1, 2 => 1, 2 => 2]

        subgraph, flows = HHN.componentsflowgraph([1 1; 0 1], [[1], []], [[], [1]])
        @test isempty(flows)
        @test isempty(subgraph)

        subgraph, flows = HHN.componentsflowgraph([1 0; 1 1], [[1], []], [[], [1]])
        @test flows == [(1 => 2, HHN.FlowInfo(1, 1.0))]
        @test subgraph == [1 => 2]
    end

    @testset "3x3" begin
        subgraph, flows = HHN.componentsflowgraph([0 0 0; 1 0 0; 0 1 0], [[1], [], []], [[], [], [1]])
        @test flows == [(1 => 3, HHN.FlowInfo(2, 1.0))]
        @test sort(subgraph) == [1 => 2, 2 => 3]

        subgraph, flows = HHN.componentsflowgraph([0 0 0; 0.5 0 0; 0 0.7 0], [[1], [1], []], [[], [], [1]])
        @test sort(flows) == [(1 => 3, HHN.FlowInfo(2, 0.5)), (2 => 3, HHN.FlowInfo(1, 0.7))]
        @test sort(subgraph) == [1 => 2, 2 => 3]

        subgraph, flows = HHN.componentsflowgraph([0 0 0; 0.5 0 0; 0 0.7 0], [[], [1], []], [[], [], [1]])
        @test flows == [(2 => 3, HHN.FlowInfo(1, 0.7))]
        @test subgraph == [2 => 3]

        subgraph, flows = HHN.componentsflowgraph([0 0 0; 1 0 0; 1 1 0], [[1], [], []], [[], [], [1]])
        @test sort(flows) == [(1 => 3, HHN.FlowInfo(1, 1.0))]
        @test sort(subgraph) == [1 => 2, 1 => 3, 2 => 3]

        subgraph, flows = HHN.componentsflowgraph([0 0 0; 0.4 0 0.6; 0 0 0], [[1], [], [1]], [[], [1], []])
        @test sort(flows) == [(1 => 2, HHN.FlowInfo(1, 0.4)),
                              (3 => 2, HHN.FlowInfo(1, 0.6))]
        @test sort(subgraph) == [1 => 2, 3 => 2]

        subgraph, flows = HHN.componentsflowgraph([0 0 0; 1 1 1; 0 0 0], [[1], [1], [1]], [[], [1], []])
        @test sort(flows) == [(1 => 2, HHN.FlowInfo(1, 1.0)),
                              (2 => 2, HHN.FlowInfo(0, 1.0)),
                              (3 => 2, HHN.FlowInfo(1, 1.0))]
        @test sort(subgraph) == [1 => 2, 2 => 2, 3 => 2]
    end

    @testset "4x4" begin
        subgraph, flows = HHN.componentsflowgraph([0 0 0 0; 1 0 0 0; 0 1 0 0; 0 1 0 0],
                                    [[1], [], [], []],
                                    [[], [], [1], []])
        @test flows == [(1 => 3, HHN.FlowInfo(2, 1.0))]
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
        @test flows == [(3 => 1, HHN.FlowInfo(1, adjmtx[1, 3]))]
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
            @test haskey(pools, Vector{Bool}) && (pools[Vector{Bool}].nborrowed == 0)
            @test haskey(pools, Vector{HHN.DFSState}) && (pools[Vector{HHN.DFSState}].nborrowed == 0)
            @test haskey(pools, Dict{Int,HHN.FlowInfoWIP}) && (pools[Dict{Int,HHN.FlowInfoWIP}].nborrowed == 0)
            @test haskey(pools, Vector{Union{Dict{Int,HHN.FlowInfoWIP},Nothing}}) &&
                  (pools[Vector{Union{Dict{Int,HHN.FlowInfoWIP},Nothing}}].nborrowed == 0)
        end
    end
end

@testset "flowgraph(SCCTree)" begin
    adjmtx = weights(HHN.tarjan1983_example_graph())
    tree = HHN.scctree(adjmtx)
    subgraph, flows, conncomps = HHN.flowgraph(tree, adjmtx, [1, 4, 5, 7], [2, 4, 6],
                                           HHN.EdgeTest{Float64}(threshold=20))
    @test sort(flows) == [(1 => 6, 4 => 5, HHN.FlowInfo(1, adjmtx[6, 1])),
                          (4 => 4, 2 => 2, HHN.FlowInfo(0, adjmtx[4, 4])),
                          (5 => 2, 3 => 1, HHN.FlowInfo(1, adjmtx[7, 5])), # 7 is a member of 1st component that provides highest weight
                          (5 => 6, 3 => 5, HHN.FlowInfo(2, adjmtx[6, 7])),
                          (7 => 2, 1 => 1, HHN.FlowInfo(0, adjmtx[3, 7])), # 3 - 7 has highest walk weight
                          (7 => 6, 1 => 5, HHN.FlowInfo(1, adjmtx[6, 7]))]
    @test sort(subgraph) == [(1 => 6, 4 => 5), (2 => 7, 1 => 1), (3 => 2, 1 => 1), (3 => 7, 1 => 1),
                             (5 => 7, 3 => 1), (7 => 3, 1 => 1), (7 => 6, 1 => 5)]
    @test HHN.flowstats(tree, adjmtx, [1, 4, 5, 7], [2, 4, 6],
                        HHN.EdgeTest{Float64}(threshold=20), sourcesinkweights=adjmtx) ==
                       (nflows=6, ncompflows=6,
                        flowlen_sum=5, compflowlen_sum=5,
                        flowinvlen_sum=23/6, compflowinvlen_sum=23/6,
                        compflowlen_max=2,
                        floweight_sum = 46.0, compfloweight_sum = 161.0,
                        flowavghopweight_sum = 23.0,
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
        @test haskey(pools, Vector{Bool}) && (pools[Vector{Bool}].nborrowed == 0)
        @test haskey(pools, Vector{HHN.DFSState}) && (pools[Vector{HHN.DFSState}].nborrowed == 0)
        @test haskey(pools, Dict{Int,HHN.FlowInfoWIP}) && (pools[Set{Int}].nborrowed == 0)
        @test haskey(pools, Vector{Union{Dict{Int,HHN.FlowInfoWIP},Nothing}}) && (pools[Vector{Union{Dict{Int,HHN.FlowInfoWIP},Nothing}}].nborrowed == 0)
        @test haskey(pools, HHN.IndicesPartition) && (pools[HHN.IndicesPartition].nborrowed == 0)
    end
end

@testset "traceflows()" begin
    FlowPathsDict = Dict{HHN.Diedge, HHN.Partition{Int}}

    @testset "trivial" begin
        @test_throws DimensionMismatch HHN.traceflows(fill(0.0, (1, 2)), HHN.EdgeTest{Float64}(threshold=0.5),
                                                      fill(0.0, (1, 2)), HHN.EdgeTest{Float64}(threshold=0.5))
        @test_throws DimensionMismatch HHN.traceflows(fill(0.0, (1, 1)), HHN.EdgeTest{Float64}(threshold=0.5),
                                                      fill(0.0, (2, 2)), HHN.EdgeTest{Float64}(threshold=0.5))

        @test HHN.traceflows(fill(0.25, (1, 1)), HHN.EdgeTest{Float64}(threshold=0.5),
                             fill(0.5, (1, 1)), HHN.EdgeTest{Float64}(threshold=0.5)) == FlowPathsDict()
        @test HHN.traceflows(fill(0.25, (1, 1)), HHN.EdgeTest{Float64}(threshold=0.25),
                             fill(0.5, (1, 1)), HHN.EdgeTest{Float64}(threshold=0.75)) == FlowPathsDict()
        @test HHN.traceflows(fill(0.25, (1, 1)), HHN.EdgeTest{Float64}(threshold=0.25),
                             fill(0.5, (1, 1)), HHN.EdgeTest{Float64}(threshold=0.5)) == FlowPathsDict() # the flow exist, but we don't create loop edges
    end

    @testset "simple" begin
        stepmtx = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
        walkmtx = [0.0 0.0 0.0; 0.0 0.0 0.0; 1.0 0.0 0.0]
        @test HHN.traceflows(stepmtx, HHN.EdgeTest{Float64}(threshold=0.5),
                             walkmtx, HHN.EdgeTest{Float64}(threshold=0.5)) == Dict((1=>3) => [[2]])
        @test HHN.traceflows(stepmtx, HHN.EdgeTest{Float64}(threshold=0.5),
                             walkmtx, HHN.EdgeTest{Float64}(threshold=0.5),
                             sinks = [2]) == FlowPathsDict()
        @test HHN.traceflows(stepmtx, HHN.EdgeTest{Float64}(threshold=0.5),
                             walkmtx, HHN.EdgeTest{Float64}(threshold=0.5),
                             sources = [2]) == FlowPathsDict()
        @test HHN.traceflows(stepmtx, HHN.EdgeTest{Float64}(threshold=0.5),
                             walkmtx, HHN.EdgeTest{Float64}(threshold=0.5),
                             sinks = [1], sources = [2]) == FlowPathsDict()

        stepmtx = [0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0]
        walkmtx = [0.0 0.0 0.0; 0.0 0.0 0.0; 1.0 0.0 0.0]
        @test HHN.traceflows(stepmtx, HHN.EdgeTest{Float64}(threshold=0.5),
                             walkmtx, HHN.EdgeTest{Float64}(threshold=0.5)) == Dict((1=>3) => [[2], Int[]])
    end
end
