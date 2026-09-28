using LinearAlgebra

const OB = OpticsBase

const PORT = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
const λ = 1.064e-6

@testset "SampledField scalar" begin
    g = RegularGrid((8, 6), (1e-6, 2e-6))
    E = rand(ComplexF64, 8, 6)
    f = SampledField(E, g, PORT, λ)
    @test f isa SampledField{1, 2, Float64}
    @test f isa AbstractOpticalField
    @test size(OB.field_array(f)) == (8, 6, 1)
    # No copy: the stored array shares memory with the input
    @test vec(OB.field_array(f)) == vec(E)
    E[3, 4] = 42 + 1im
    @test OB.field_array(f)[3, 4, 1] == 42 + 1im
    OB.field_array(f)[1, 1, 1] = -7im
    @test E[1, 1] == -7im
    @test Base.mightalias(OB.field_array(f), E)

    @test OB.grid(f) === g
    @test OB.port(f) === PORT
    @test OB.wavelength(f) == λ
    @test !is_vectorial(f)
    @test is_coherent(f)

    # Explicit singleton component dimension is also scalar
    E1 = zeros(ComplexF64, 8, 6, 1)
    f1 = SampledField(E1, g, PORT, λ)
    @test f1 isa SampledField{1, 2}
    @test OB.field_array(f1) === E1

    # Wavelength promoted to the element type of E
    f32 = SampledField(zeros(ComplexF32, 8, 6), g, PORT, 1.064e-6)
    @test OB.wavelength(f32) isa Float32
    @test f32 isa SampledField{1, 2, Float32}
end

@testset "SampledField vectorial" begin
    g = RegularGrid((8, 6), (1e-6, 1e-6))
    E = rand(ComplexF64, 8, 6, 3)
    f = SampledField(E, g, PORT, λ)
    @test f isa SampledField{3, 2, Float64}
    @test OB.field_array(f) === E
    @test is_vectorial(f)
    @test is_coherent(f)

    # Generic array types are kept as passed (no collect)
    Ev = view(zeros(ComplexF64, 10, 6, 3), 1:8, :, :)
    fv = SampledField(Ev, g, PORT, λ)
    @test OB.field_array(fv) === Ev

    # 3D grid (volume)
    g3 = RegularGrid((4, 5, 6), (1e-6, 1e-6, 1e-6))
    f3 = SampledField(zeros(ComplexF64, 4, 5, 6, 3), g3, PORT, λ)
    @test f3 isa SampledField{3, 3}
    @test SampledField(zeros(ComplexF64, 4, 5, 6), g3, PORT, λ) isa SampledField{1, 3}
end

@testset "SampledField invariants" begin
    g = RegularGrid((8, 6), (1e-6, 1e-6))
    # N must be 1 or 3
    @test_throws ArgumentError SampledField(zeros(ComplexF64, 8, 6, 2), g, PORT, λ)
    @test_throws ArgumentError SampledField(zeros(ComplexF64, 8, 6, 4), g, PORT, λ)
    # Grid size mismatch
    @test_throws ArgumentError SampledField(zeros(ComplexF64, 6, 8), g, PORT, λ)
    @test_throws ArgumentError SampledField(zeros(ComplexF64, 8, 5, 3), g, PORT, λ)
    # Wrong number of dimensions
    @test_throws ArgumentError SampledField(zeros(ComplexF64, 8), g, PORT, λ)
    @test_throws ArgumentError SampledField(zeros(ComplexF64, 8, 6, 3, 1), g, PORT, λ)
    # Wavelength
    @test_throws ArgumentError SampledField(zeros(ComplexF64, 8, 6), g, PORT, 0.0)
    @test_throws ArgumentError SampledField(zeros(ComplexF64, 8, 6), g, PORT, -1e-6)
end

@testset "SampledField total_power" begin
    # Sampled Gaussian: |E₀|² = 4 Z₀ P / (n π w²) gives power P
    n = 1.5
    port = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0];
        refractive_index = n)
    w = 1e-3
    P = 1.0
    Δ = 20e-6
    g = RegularGrid((512, 512), (Δ, Δ))
    ξ = OB.coordinates(g, 1)
    η = OB.coordinates(g, 2)
    E0 = sqrt(4 * OB.VACUUM_IMPEDANCE * P / (n * π * w^2))
    # Tilted phase must not change the power
    E = [E0 * exp(-(x^2 + y^2) / w^2) * cis(2π * x / 50e-6) for x in ξ, y in η]
    f = SampledField(E, g, port, λ)
    @test abs(total_power(f) - P) < 1e-9

    # Vectorial: the power is summed over all components
    Ev = zeros(ComplexF64, 512, 512, 3)
    Ev[:, :, 1] .= E ./ sqrt(2)
    Ev[:, :, 2] .= im .* E ./ sqrt(2)
    fv = SampledField(Ev, g, port, λ)
    @test abs(total_power(fv) - P) < 1e-9

    # Anisotropic spacing enters as Δξ·Δη
    ga = RegularGrid((4, 4), (1e-6, 3e-6))
    fa = SampledField(ones(ComplexF64, 4, 4), ga, PORT, λ)
    @test total_power(fa) ≈ OB.power_normalization(PORT) * 16 * 3e-12 rtol = 1e-14

    # Volume: power is undefined
    g3 = RegularGrid((4, 4, 4), (1e-6, 1e-6, 1e-6))
    @test_throws ArgumentError total_power(SampledField(ones(ComplexF64, 4, 4, 4), g3, PORT, λ))

    @inferred total_power(f)
end
