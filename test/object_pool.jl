@testset "ArrayPool" begin
    begin
        pool = HHN.ArrayPool{Int8}()
        @test eltype(pool) === Int8
        @test pool.nborrowed == 0
        @inferred Vector{Int8} HHN.borrow!(pool, 2)
        @inferred Matrix{Int8} HHN.borrow!(pool, (2, 4))
        @inferred Array{Int8, 3} HHN.borrow!(pool, (2, 4, 2))
        @inferred Matrix{Int8} HHN.borrow!(Int8, pool, (2, 4))
        @test_throws MethodError HHN.borrow!(Int, pool, (2, 4))
    end

    pool1 = HHN.ArrayPool{Int}()
    @test eltype(pool1) === Int
    @test pool1.nborrowed == 0

    begin
        v = HHN.borrow!(pool1, 2)
        @test v isa Vector{Int}
        @test length(v) == 2
        @test pool1.nborrowed == 1
        vv = HHN.borrow!(pool1)
        @test isempty(vv)
        @test pool1.nborrowed == 2
        HHN.release!(pool1, vv)
        vvv = HHN.borrow!(Int, pool1, 1)
        @test vvv isa Vector{Int}
        HHN.release!(pool1, vvv)
        @test_throws MethodError HHN.borrow!(Float64, pool1)
        @test pool1.nborrowed == 1
        @test_throws Exception HHN.release!(pool1, Vector{UInt}())
        HHN.release!(pool1, v)
        @test pool1.nborrowed == 0
        # over-releasing
        @test_throws Exception HHN.release!(pool1, Vector{Int}())
    end

    begin
        @test pool1.nborrowed == 0
        m = HHN.borrow!(pool1, (2, 4))
        @test m isa Matrix{Int}
        @test size(m) == (2, 4)
        HHN.release!(pool1, m)
        @test pool1.nborrowed == 0
    end

    pool2 = HHN.ArrayPool{String}(2)
    @test eltype(HHN.ArrayPool{Char}) === Char
    @test eltype(pool2) === String
    begin
        v1 = HHN.borrow!(pool2, 5)
        @test v1 isa Vector{String}
        v2 = HHN.borrow!(pool2, 15)
        @test_throws Exception HHN.borrow!(pool2, 1)
        HHN.release!(pool2, v2)
        v3 = HHN.borrow!(pool2, 10)
        HHN.release!(pool2, v3)
        HHN.release!(pool2, v1)
        @test pool2.nborrowed == 0
    end
end

@testset "NoopArrayPool" begin
    strnopool = HHN.NoopArrayPool{String}(2)

    @test eltype(HHN.NoopArrayPool{Char}) === Char
    @test eltype(strnopool) === String
    begin
        v1 = @inferred HHN.borrow!(strnopool, 5)
        @test v1 isa Vector{String}
        @test length(v1) == 5
        v2 = @inferred HHN.borrow!(strnopool, (15, 4))
        @test size(v2) == (15, 4)
        HHN.release!(strnopool, v2)
        HHN.release!(strnopool, v1)
        @test_throws MethodError HHN.release!(strnopool, Int[])
    end
end

@testset "ObjectPools" begin
    pools = HHN.ObjectPools()
    @test pools.default_borrow_limit == 0
    @test !haskey(pools, Vector{Int})
    intpool = @inferred pools[Vector{Int}]
    @test intpool isa HHN.ArrayPool{Int}
    @test haskey(pools, Vector{Int})
    intpool2 = pools[Vector{Int}]
    @test intpool2 === intpool
    intpool3 = HHN.arraypool(pools, Int)
    @test intpool3 === intpool

    floatpool = @inferred pools[Vector{Float64}]
    @test floatpool isa HHN.ArrayPool{Float64}

    nointpool = @inferred HHN.arraypool(nothing, Int)
    @test nointpool isa HHN.NoopArrayPool{Int}
end
