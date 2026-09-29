using LinearAlgebra
using StaticArrays

const OB = OpticsBase

const λ = 1.0e-6
const NIDX = 1.5
const PORT = PlanarPort([1e-3, 0.0, -2e-3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0];
    refractive_index = NIDX)

# Unit direction at polar angle θ and azimuth φ about the port normal (z)
oblique(θ, φ) = SVector(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))

const DIR = [SVector(0.0, 0.0, 1.0), oblique(deg2rad(50), deg2rad(30)),
    oblique(deg2rad(10), deg2rad(-100))]
const AMP1 = [1.0e6 + 0im, 2.0e6 * cis(0.4), -0.5e6im]
const WGT = [1e-4, 2e-4, 5e-5]
# Transverse vectorial amplitudes from Jones vectors in the ray basis
const JONES = [SVector(1.0 + 0im, 0), OB.circular_jones(+1), SVector(0.3 + 0im, -0.8im)]
const AMP3 = [1e6 * OB.jones_to_global(PORT, DIR[j], JONES[j]) for j in 1:3]

scalar_spectrum() = PlaneWaveSpectrum(PORT, λ, DIR, AMP1, WGT)
vector_spectrum() = PlaneWaveSpectrum(PORT, λ, DIR, AMP3, WGT)

# Replaces entry i of v by x (copy)
with(v, i, x) = (w = collect(v); w[i] = x; w)

κλ²(port) = OB.power_normalization(port) * (λ / OB.refractive_index(port))^2

@testset "scalar spectrum" begin
    s = scalar_spectrum()
    @test s isa PlaneWaveSpectrum{1, Float64}
    @test s isa AbstractScalarField
    @test PlaneWaveSpectrum{3} <: AbstractVectorField <: AbstractOpticalData
    @test length(s) == 3
    @test !is_vectorial(s)
    @test is_coherent(s)
    @test OB.port(s) === PORT
    @test OB.wavelength(s) == λ
    @test s.amplitude[2] == SVector(2.0e6 * cis(0.4))
    @test s.direction == DIR && s.weight == WGT
    @test s.weight !== WGT            # copied, not aliased
    P_ref = κλ²(PORT) * sum(WGT .* abs2.(AMP1))
    @test total_power(s) ≈ P_ref rtol = 1e-14
end

@testset "vectorial spectrum" begin
    s = vector_spectrum()
    @test s isa PlaneWaveSpectrum{3, Float64}
    @test is_vectorial(s)
    @test is_coherent(s)
    P_ref = κλ²(PORT) * sum(WGT[j] * norm(AMP3[j])^2 for j in 1:3)
    @test total_power(s) ≈ P_ref rtol = 1e-14
    # Plain Vector{Vector} inputs work as well
    s2 = PlaneWaveSpectrum(PORT, λ, collect.(DIR), collect.(AMP3), WGT)
    @test s2 isa PlaneWaveSpectrum{3, Float64}
    @test s2.amplitude == s.amplitude
end

@testset "total_power of one sample" begin
    d = oblique(deg2rad(35), deg2rad(70))
    e = OB.jones_to_global(PORT, d, SVector(0.6 + 0.2im, -0.3im)) * 4.2e7
    w = 3.7e-5
    s = PlaneWaveSpectrum(PORT, λ, [d], [e], [w])
    κ = OB.power_normalization(PORT)
    λm = λ / NIDX
    @test total_power(s) ≈ κ * λm^2 * w * norm(e)^2 rtol = 1e-14
    a = 3.0e7 + 1.0e7im
    s1 = PlaneWaveSpectrum(PORT, λ, [d], [a], [w])
    @test total_power(s1) ≈ κ * λm^2 * w * abs2(a) rtol = 1e-14
end

@testset "empty spectrum" begin
    s = PlaneWaveSpectrum(PORT, λ, SVector{3, Float64}[], ComplexF64[], Float64[])
    @test length(s) == 0
    @test total_power(s) == 0
    @test !is_vectorial(s)
    s3 = PlaneWaveSpectrum(PORT, λ, SVector{3, Float64}[], SVector{3, ComplexF64}[],
        Float64[])
    @test is_vectorial(s3)
    @test total_power(s3) == 0
    # N cannot be inferred from an untyped empty vector
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, [], [], [])
end

@testset "promotion" begin
    s = PlaneWaveSpectrum(PORT, Float32(λ), SVector{3, Float32}.(DIR),
        ComplexF32.(AMP1), Float32.(WGT))
    @test s isa PlaneWaveSpectrum{1, Float32}
    @test total_power(s) isa Float32
    @test OB.wavelength(s) isa Float32
    s3 = PlaneWaveSpectrum(PORT, Float32(λ), SVector{3, Float32}.(DIR),
        SVector{3, ComplexF32}.(AMP3), Float32.(WGT))
    @test s3 isa PlaneWaveSpectrum{3, Float32}
    # Mixed inputs promote to the widest type; integer weights and a real amplitude work
    sm = PlaneWaveSpectrum(PORT, Float32(λ), DIR, [1, 2, 3], [1, 1, 1])
    @test sm isa PlaneWaveSpectrum{1, Float64}
    @test sm.amplitude[3] == SVector(3.0 + 0im)
end

@testset "invariant 1: components and lengths" begin
    amp2 = [SVector(1.0 + 0im, 0.0), SVector(1.0 + 0im, 0.0), SVector(1.0 + 0im, 0.0)]
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, DIR, amp2, WGT)
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, DIR[1:2], AMP1, WGT)
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, DIR, AMP1[1:2], WGT)
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, DIR, AMP1, WGT[1:2])
    # Inconsistent amplitude lengths and 2-element directions
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, DIR,
        [collect(AMP3[1]), collect(AMP3[2])[1:2], collect(AMP3[3])], WGT)
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, with(collect.(DIR), 2, [0.0, 1.0]),
        AMP1, WGT)
end

@testset "invariant 2: wavelength and weights" begin
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, 0.0, DIR, AMP1, WGT)
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, -λ, DIR, AMP1, WGT)
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, DIR, AMP1, with(WGT, 2, 0.0))
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, DIR, AMP1, with(WGT, 3, -1e-4))
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, DIR, AMP1, with(WGT, 1, Inf))
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, DIR, AMP1, with(WGT, 1, NaN))
    err = try
        PlaneWaveSpectrum(PORT, λ, DIR, AMP1, with(WGT, 2, 0.0))
    catch e
        e
    end
    @test occursin("sample 2", err.msg)
end

@testset "invariant 3: directions" begin
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, with(DIR, 2, 1.01 * DIR[2]), AMP1,
        WGT)
    # Backward and grazing directions
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, with(DIR, 3, -DIR[3]), AMP1, WGT)
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ,
        with(DIR, 1, SVector(1.0, 0.0, 0.0)), AMP1, WGT)
    # Unit norm within √eps is accepted
    s = PlaneWaveSpectrum(PORT, λ, with(DIR, 2, (1 + 1e-10) * DIR[2]), AMP1, WGT)
    @test length(s) == 3
end

@testset "invariant 4: transversality" begin
    # Longitudinal component
    bad = with(AMP3, 2, AMP3[2] + 1e3 * DIR[2])
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, DIR, bad, WGT)
    # Transversality is without conjugation: e·d = 0 must hold, not conj(e)·d
    d = DIR[1]
    e = SVector(1.0 + 0im, im, 0.0)
    @test PlaneWaveSpectrum(PORT, λ, [d], [e], [1e-4]) isa PlaneWaveSpectrum{3}
    d2 = oblique(deg2rad(40), 0.0)
    e2 = SVector(im * d2[3], 0.0, -im * d2[1])           # transverse, purely imaginary
    @test PlaneWaveSpectrum(PORT, λ, [d2], [e2], [1e-4]) isa PlaneWaveSpectrum{3}
    e3 = SVector(1.0 + 0im, 0.0, im)                   # e·d2 ≠ 0
    @test_throws ArgumentError PlaneWaveSpectrum(PORT, λ, [d2], [e3], [1e-4])
    # Zero amplitude is allowed
    @test PlaneWaveSpectrum(PORT, λ, [d2], [zero(SVector{3, ComplexF64})], [1e-4]) isa
          PlaneWaveSpectrum{3}
    # The scalar case has no transversality check
    @test PlaneWaveSpectrum(PORT, λ, [d2], [1.0 + 0im], [1e-4]) isa PlaneWaveSpectrum{1}
end
