using LinearAlgebra
using StaticArrays

const OB = OpticsBase

const λ = 633e-9
const PORT = PlanarPort([0.1, 0.2, 0.3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0];
    refractive_index = 1.5)

# Oblique unit direction at polar angle θ and azimuth φ about the port normal (z)
oblique(θ, φ) = SVector(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))

# Three rays on the port plane: on axis and two oblique ones
const POS = [SVector(0.1, 0.2, 0.3), SVector(0.1 + 1e-3, 0.2, 0.3), SVector(0.1, 0.2 - 2e-3, 0.3)]
const DIR = [SVector(0.0, 0.0, 1.0), oblique(deg2rad(64), deg2rad(45)), oblique(deg2rad(20), deg2rad(-120))]
const OPL = [1.0, 1.0 + 1e-6, 1.0 - 2e-6]
const PWR = [1e-3, 2e-3, 0.5e-3]
const PH1 = cis.([0.0, 0.3, -1.2])
# Vectorial phasors from Jones vectors in the ray basis: x-, circular and oblique linear
const JONES = [SVector(1.0 + 0im, 0), OB.circular_jones(+1), SVector(cos(0.4) + 0im, sin(0.4))]
const PH3 = [OB.jones_to_global(PORT, DIR[i], JONES[i]) for i in 1:3]

scalar_bundle(; kw...) = RayBundle(PORT, λ, POS, DIR, OPL, PWR, PH1; kw...)
vector_bundle(; kw...) = RayBundle(PORT, λ, POS, DIR, OPL, PWR, PH3; kw...)

# Replaces entry i of v by x (copy)
with(v, i, x) = (w = collect(v); w[i] = x; w)

@testset "scalar bundle without beamlets" begin
    b = scalar_bundle()
    @test b isa RayBundle{1, Float64}
    @test b isa AbstractOpticalField
    @test length(b) == 3
    @test !is_vectorial(b)
    @test is_coherent(b)
    @test !OB.has_beamlets(b)
    @test b.beamlet === nothing
    @test total_power(b) ≈ 3.5e-3 rtol = 1e-15
    @test OB.port(b) === PORT
    @test OB.wavelength(b) == λ
    @test b.phasor[2] == SVector(cis(0.3))
    @test b.position == POS && b.direction == DIR && b.opl == OPL && b.power == PWR
    # Inputs are copied, not aliased
    @test b.opl !== OPL && b.power !== PWR
end

@testset "vectorial bundle without beamlets" begin
    b = vector_bundle()
    @test b isa RayBundle{3, Float64}
    @test is_vectorial(b)
    @test is_coherent(b)
    @test !OB.has_beamlets(b)
    @test length(b) == 3
    @test total_power(b) ≈ 3.5e-3 rtol = 1e-15
    # Phasors from jones_to_global on oblique rays satisfy invariant 5
    for i in 1:3
        @test abs(transpose(b.phasor[i]) * b.direction[i]) < 1e-14
        @test norm(b.phasor[i]) ≈ 1 atol = 1e-14
    end
    # Plain Vector{Vector} inputs work as well
    b2 = RayBundle(PORT, λ, collect.(POS), collect.(DIR), OPL, PWR, collect.(PH3))
    @test b2 isa RayBundle{3, Float64}
    @test b2.phasor == b.phasor
end

@testset "beamlets" begin
    w = 1e-3
    λm = λ / 1.5
    invq = 1 / 0.5 + im * λm / (π * w^2)
    # Stigmatic 1/q values are turned into Q = (1/q)·I
    b = scalar_bundle(; beamlet = fill(invq, 3))
    @test OB.has_beamlets(b)
    @test b.beamlet isa Vector{SMatrix{2, 2, ComplexF64, 4}}
    @test b.beamlet[1] == SMatrix{2, 2}(invq, 0, 0, invq)
    # Astigmatic, symmetric Q with off-diagonal terms
    Q = SMatrix{2, 2}(2.0 + 1.0im, 0.3 + 0.1im, 0.3 + 0.1im, -1.0 + 2.0im)
    bv = vector_bundle(; beamlet = fill(Q, 3))
    @test OB.has_beamlets(bv)
    @test is_vectorial(bv)
    @test bv.beamlet[2] == Q
    @test total_power(bv) ≈ 3.5e-3 rtol = 1e-15
    # Plain matrices are accepted
    @test vector_bundle(; beamlet = [Matrix(Q) for _ in 1:3]).beamlet[3] == Q
end

@testset "type promotion" begin
    p32 = PlanarPort(Float32[0, 0, 0], Float32[0, 0, 1], Float32[1, 0, 0])
    b = RayBundle(p32, 633.0f-9, [SVector{3, Float32}(0, 0, 0)],
        [SVector{3, Float32}(0, 0, 1)], Float32[0], Float32[1], ComplexF32[1])
    @test b isa RayBundle{1, Float32}
    @test total_power(b) isa Float32
    @test b.beamlet === nothing
    # Integer and mixed inputs promote to a common float type
    b = RayBundle(PORT, λ, [[0.1, 0.2, 0.3]], [[0, 0, 1]], [0], [1], [1])
    @test b isa RayBundle{1, Float64}
    @test b.phasor[1] isa SVector{1, ComplexF64}
    @test typeof(OB.wavelength(b)) == Float64
    # Empty bundle
    e = RayBundle(PORT, λ, SVector{3, Float64}[], SVector{3, Float64}[], Float64[],
        Float64[], SVector{3, ComplexF64}[])
    @test length(e) == 0
    @test total_power(e) === 0.0
    @test is_vectorial(e)
    @test_throws ArgumentError RayBundle(PORT, λ, [], [], Float64[], Float64[], [])
end

@testset "type stability" begin
    @test (@inferred total_power(vector_bundle())) ≈ 3.5e-3
    @test !(@inferred OB.has_beamlets(scalar_bundle()))
    @test @inferred is_vectorial(vector_bundle())
    # Construction from SVector / number inputs is inferred
    @test (@inferred scalar_bundle()) isa RayBundle{1, Float64}
    @test (@inferred vector_bundle()) isa RayBundle{3, Float64}
end

@testset "invariant 1: N and lengths" begin
    # N = 2 (e.g. a Jones vector passed instead of a global field vector)
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR, OPL, PWR, JONES)
    # Mixed phasor lengths
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR, OPL, PWR,
        Any[PH3[1], PH3[2], SVector(1.0 + 0im)])
    # Unequal lengths of each per-ray vector
    @test_throws ArgumentError RayBundle(PORT, λ, POS[1:2], DIR, OPL, PWR, PH1)
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR[1:2], OPL, PWR, PH1)
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR, OPL[1:2], PWR, PH1)
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR, OPL, PWR[1:2], PH1)
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR, OPL, PWR, PH1[1:2])
    @test_throws ArgumentError scalar_bundle(; beamlet = fill(1 + 1im, 2))
    # Position with the wrong number of elements
    @test_throws ArgumentError RayBundle(PORT, λ, with(collect.(POS), 1, [0.1, 0.2]),
        DIR, OPL, PWR, PH1)
end

@testset "invariant 2: wavelength and power" begin
    @test_throws ArgumentError RayBundle(PORT, 0.0, POS, DIR, OPL, PWR, PH1)
    @test_throws ArgumentError RayBundle(PORT, -λ, POS, DIR, OPL, PWR, PH1)
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR, OPL, with(PWR, 2, -1e-9), PH1)
    # Zero power is allowed
    @test total_power(RayBundle(PORT, λ, POS, DIR, OPL, zeros(3), PH1)) == 0
end

@testset "invariant 3: direction" begin
    # Not normalized (no silent normalization)
    @test_throws ArgumentError RayBundle(PORT, λ, POS, with(DIR, 2, 2 * DIR[2]), OPL,
        PWR, PH1)
    # Propagating against the normal, and grazing
    @test_throws ArgumentError RayBundle(PORT, λ, POS, with(DIR, 1, -DIR[1]), OPL, PWR, PH1)
    @test_throws ArgumentError RayBundle(PORT, λ, POS, with(DIR, 1, SVector(1.0, 0, 0)),
        OPL, PWR, PH1)
end

@testset "invariant 4: position on the port plane" begin
    off = POS[2] + SVector(0.0, 0.0, 1e-6)
    @test_throws ArgumentError RayBundle(PORT, λ, with(POS, 2, off), DIR, OPL, PWR, PH1)
    # Rounding-level deviations are accepted
    tiny = POS[2] + SVector(0.0, 0.0, 1e-12)
    @test length(RayBundle(PORT, λ, with(POS, 2, tiny), DIR, OPL, PWR, PH1)) == 3
end

@testset "invariant 5: phasor norm and transversality" begin
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR, OPL, PWR, with(PH1, 3, 2.0 + 0im))
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR, OPL, PWR, with(PH3, 1, 2 * PH3[1]))
    # Unit but longitudinal field vector on the oblique ray
    long = complex.(DIR[2])
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR, OPL, PWR, with(PH3, 2, long))
    # Transversality is checked without conjugation: e = (x̂ + i d)/√2 has e·d = i/√2
    x̂, _ = OB.ray_basis(PORT, DIR[2])
    e = (x̂ + im * DIR[2]) / sqrt(2)
    @test_throws ArgumentError RayBundle(PORT, λ, POS, DIR, OPL, PWR, with(PH3, 2, e))
end

@testset "invariant 6: beamlet matrix" begin
    invq = 1 / 0.5 + 1e-3im
    # Not symmetric
    Qa = SMatrix{2, 2}(1.0 + 1.0im, 0.5, 0.0, 1.0 + 1.0im)
    @test_throws ArgumentError scalar_bundle(; beamlet = [Qa, Qa, Qa])
    # Im Q not positive definite: negative trace, indefinite (det < 0), zero (plane wave)
    @test_throws ArgumentError scalar_bundle(; beamlet = fill(conj(invq), 3))
    Qi = SMatrix{2, 2}(1.0 + 1.0im, 0, 0, 1.0 - 0.5im)
    @test_throws ArgumentError scalar_bundle(; beamlet = fill(Qi, 3))
    @test_throws ArgumentError scalar_bundle(; beamlet = fill(1.0 + 0im, 3))
    # Wrong size
    @test_throws ArgumentError scalar_bundle(; beamlet = fill(Matrix{ComplexF64}(I, 3, 3), 3))
end
