const OB = OpticsBase

@testset "RegularGrid construction" begin
    g = RegularGrid((4, 5), (1e-6, 2e-6))
    @test g isa RegularGrid{2, Float64}
    @test g isa AbstractGrid{2}
    @test size(g) == (4, 5)
    @test size(g, 2) == 5
    @test OB.spacing(g) == (1e-6, 2e-6)

    # Promotion to a common floating-point type
    @test RegularGrid((3, 3), (1, 2)) isa RegularGrid{2, Float64}
    @test RegularGrid((3, 3), (1.0f0, 2.0f0)) isa RegularGrid{2, Float32}
    @test RegularGrid((3, 3), (1.0f0, 2.0)) isa RegularGrid{2, Float64}
    @test RegularGrid((3, Int32(3)), (1.0, 1.0)).dims isa NTuple{2, Int}

    # 3D grid
    g3 = RegularGrid((2, 3, 4), (1e-6, 1e-6, 5e-6))
    @test g3 isa AbstractGrid{3}
    @test size(g3) == (2, 3, 4)

    # Invariants
    @test_throws ArgumentError RegularGrid((0, 4), (1e-6, 1e-6))
    @test_throws ArgumentError RegularGrid((4, -1), (1e-6, 1e-6))
    @test_throws ArgumentError RegularGrid((4, 4), (0.0, 1e-6))
    @test_throws ArgumentError RegularGrid((4, 4), (1e-6, -1e-6))
    @test_throws ArgumentError RegularGrid((4, 4), (Inf, 1e-6))
    @test_throws ArgumentError RegularGrid((4, 4), (NaN, 1e-6))
end

@testset "RegularGrid coordinates" begin
    for n in (1, 2, 5, 6, 7, 512, 513), Δ in (1e-6, 20e-6, 0.3)
        g = RegularGrid((n, 3), (Δ, 1.0))
        x = OB.coordinates(g, 1)
        @test length(x) == n
        @test x[n ÷ 2 + 1] == 0
        @test all(i -> x[i] == (i - (n ÷ 2 + 1)) * Δ, 1:n)
        @test step(x) == Δ
    end
    # Even n: one more sample on the negative side (fftshift convention)
    x = OB.coordinates(RegularGrid((4, 4), (1.0, 1.0)), 1)
    @test collect(x) == [-2.0, -1.0, 0.0, 1.0]
    x = OB.coordinates(RegularGrid((5, 5), (1.0, 1.0)), 2)
    @test collect(x) == [-2.0, -1.0, 0.0, 1.0, 2.0]

    # Per-dimension spacing and third dimension
    g3 = RegularGrid((4, 5, 6), (1e-6, 2e-6, 3e-6))
    @test OB.coordinates(g3, 2)[1] == -2 * 2e-6
    @test OB.coordinates(g3, 3)[6 ÷ 2 + 1] == 0
    @test step(OB.coordinates(g3, 3)) == 3e-6

    # Element type and bounds
    g32 = RegularGrid((4, 4), (1.0f-6, 1.0f-6))
    @test eltype(OB.coordinates(g32, 1)) == Float32
    @test_throws BoundsError OB.coordinates(g32, 3)
    @test_throws BoundsError OB.coordinates(g32, 0)
    @inferred OB.coordinates(g32, 1)
end
