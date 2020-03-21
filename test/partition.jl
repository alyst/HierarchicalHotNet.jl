@testset "Partition" begin
    @testset "Empty partition" begin
        ptn0 = HHN.IndicesPartition(0)
        @test length(ptn0) == 0
        @test isempty(ptn0)
        @test eltype(ptn0) === typeof(view([1,4,2], 2:3))
        @test HHN.nparts(ptn0) == 0
        @test HHN.nelems(ptn0) == 0
        @test HHN.elems(ptn0) == Int[]

        HHN.pushelem!(ptn0, 3)
        @test HHN.nelems(ptn0) == 1
        @test HHN.nparts(ptn0) == 0
        HHN.closepart!(ptn0)
        @test HHN.nelems(ptn0) == 1
        @test HHN.nparts(ptn0) == 1
        @test length(ptn0) == 1
        @test ptn0[1] == [3]
        @test ptn0 == [[3]]

        @test !isempty(ptn0)
        empty!(ptn0)
        @test isempty(ptn0)
    end

    @testset "1-element partition" begin
        ptn1 = HHN.IndicesPartition(1)
        @test length(ptn1) == 1
        @test HHN.nelems(ptn1) == 1
        @test ptn1 == [[1]]
        @test HHN.elems(ptn1) == [1]
        i = 0
        for pt in ptn1
            i += 1
            @test (i == 1) && (pt == [1])
        end
        @test i == 1
    end

    @testset "adding parts" begin
        ptn = HHN.Partition{Char}()
        @test HHN.nparts(ptn) == 0
        HHN.pushelem!(ptn, 'A')
        HHN.closepart!(ptn)
        @test HHN.nparts(ptn) == 1
        HHN.pushelem!(ptn, 'C')
        HHN.pushelem!(ptn, 'B')
        HHN.closepart!(ptn)
        @test HHN.nelems(ptn) == 3
        @test HHN.nparts(ptn) == 2
        @test ptn[1] == ['A']
        @test ptn[2] == ['C', 'B']
        @test HHN.elems(ptn) == ['A', 'C', 'B']

        HHN.pushelem!(ptn, 'E')
        HHN.pushelem!(ptn, 'D')
        HHN.closepart!(ptn)
        @test HHN.nelems(ptn) == 5
        @test HHN.nparts(ptn) == 3
        @test length(ptn) == 3
        @test ptn[1] == ['A']
        @test ptn[2] == ['C', 'B']
        @test ptn[3] == ['E', 'D']
        sort!(ptn[3])
        @test ptn[3] == ['D', 'E']

        ptn_copy = HHN.Partition{Char}()
        copy!(ptn_copy, ptn)
        @test ptn_copy == ptn
        @test ptn_copy.starts !== ptn.starts
        @test ptn_copy.elems !== ptn.elems

        i = 0
        for pt in ptn
            i += 1
            @test (i != 1) || (pt == ['A'])
            @test (i != 2) || (pt == ['C', 'B'])
            @test (i != 3) || (pt == ['D', 'E'])
        end
        @test i == 3
    end

    @testset "1-set partition by default for non-empty" begin
        ptn = HHN.Partition(["A", "B", "C"])
        @test length(ptn) == 1
        @test ptn == [["A", "B", "C"]]

        ptn = HHN.Partition{String}(String[])
        @test length(ptn) == 0
        @test ptn == Vector{String}[]
    end

    @testset "initializing IndicesPartition" begin
        ptn = HHN.IndicesPartition(5)
        @test ptn == [[1], [2], [3], [4], [5]]
        HHN.reset!(ptn, ngroups=1)
        @test ptn == [1:5]
        @test_throws ArgumentError HHN.reset!(ptn, ngroups=3)

        @test HHN.IndicesPartition(3, ngroups=1) == [1:3]
        @test HHN.IndicesPartition(3, ngroups=3) == [[1], [2], [3]]
    end

    @testset "repeat!()" begin
        ptn = HHN.Partition(['A', 'B', 'C'])
        @test_throws ArgumentError HHN.repeat!(ptn, -1)
        @test ptn == [['A', 'B', 'C']]

        @test HHN.repeat!(ptn, 1) === ptn
        @test ptn == [['A', 'B', 'C']]
        @test HHN.repeat!(ptn, 3) === ptn
        @test ptn == [['A', 'B', 'C'], ['A', 'B', 'C'], ['A', 'B', 'C']]
        @test HHN.repeat!(ptn, 0) === ptn
        @test ptn == Vector{Char}[]

        ptn = HHN.Partition(['A', 'B', 'C'], [1, 2, 4])
        HHN.repeat!(ptn, 2)
        @test ptn == [['A'], ['B', 'C'], ['A'], ['B', 'C']]

        ptn = HHN.Partition(['A', 'B', 'C'], [1, 2, 4, 4])
        HHN.repeat!(ptn, 2)
        @test ptn == [['A'], ['B', 'C'], Char[], ['A'], ['B', 'C'], Char[]]

        ptn = HHN.Partition(Char[], [1, 1, 1])
        HHN.repeat!(ptn, 3)
        @test ptn == [Char[], Char[], Char[], Char[], Char[], Char[]]
    end
end
