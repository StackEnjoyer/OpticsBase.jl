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
# Polarizations from Jones vectors in the ray basis: x-, circular and oblique linear
const JONES = [SVector(1.0 + 0im, 0), OB.circular_jones(+1), SVector(cos(0.4) + 0im, sin(0.4))]
const POL = [OB.jones_to_global(PORT, DIR[i], JONES[i]) for i in 1:3]

rays() = RayBundle(PORT, λ, POS, DIR)
polarized() = PolarizedRayBundle(PORT, λ, POS, DIR, OPL, PWR, POL)

# Replaces entry i of v by x (copy)
with(v, i, x) = (w = collect(v); w[i] = x; w)

@testset "hierarchy" begin
    @test RayBundle <: AbstractRayBundle <: AbstractOpticalData
    @test PolarizedRayBundle <: AbstractRayBundle
    @test !(RayBundle <: AbstractOpticalField)
    @test !(PolarizedRayBundle <: AbstractOpticalField)
end

@testset "RayBundle: geometry only" begin
    b = rays()
    @test b isa RayBundle{Float64}
    @test length(b) == 3
    @test OB.port(b) === PORT
    @test OB.wavelength(b) == λ
    @test OB.positions(b) == POS && OB.directions(b) == DIR
    # Inputs are copied, not aliased
    @test OB.positions(b) !== POS
    # Plain Vector{Vector} inputs work as well
    @test OB.directions(RayBundle(PORT, λ, collect.(POS), collect.(DIR))) == DIR
    # No power, phase or polarization
    @test !hasmethod(total_power, Tuple{typeof(b)})
    @test !hasmethod(is_vectorial, Tuple{typeof(b)})
end

@testset "PolarizedRayBundle" begin
    b = polarized()
    @test b isa PolarizedRayBundle{Float64}
    @test length(b) == 3
    @test total_power(b) ≈ 3.5e-3 rtol = 1e-15
    @test OB.positions(b) == POS && OB.directions(b) == DIR
    @test b.opl == OPL && b.power == PWR && b.polarization == POL
    @test b.opl !== OPL && b.power !== PWR
    # Polarizations from jones_to_global on oblique rays are unit and transverse
    for i in 1:3
        @test abs(transpose(b.polarization[i]) * b.direction[i]) < 1e-14
        @test norm(b.polarization[i]) ≈ 1 atol = 1e-14
    end
    @test PolarizedRayBundle(PORT, λ, collect.(POS), collect.(DIR), OPL, PWR,
        collect.(POL)).polarization == POL
end

@testset "type promotion" begin
    p32 = PlanarPort(Float32[0, 0, 0], Float32[0, 0, 1], Float32[1, 0, 0])
    b = RayBundle(p32, 633.0f-9, [SVector{3, Float32}(0, 0, 0)], [SVector{3, Float32}(0, 0, 1)])
    @test b isa RayBundle{Float32}
    pb = PolarizedRayBundle(p32, 633.0f-9, [SVector{3, Float32}(0, 0, 0)],
        [SVector{3, Float32}(0, 0, 1)], Float32[0], Float32[1], [ComplexF32[1, 0, 0]])
    @test pb isa PolarizedRayBundle{Float32}
    @test total_power(pb) isa Float32
    # Integer and mixed inputs promote to a common float type
    b = RayBundle(PORT, λ, [[0.1, 0.2, 0.3]], [[0, 0, 1]])
    @test b isa RayBundle{Float64}
    pb = PolarizedRayBundle(PORT, λ, [[0.1, 0.2, 0.3]], [[0, 0, 1]], [0], [1], [[1, 0, 0]])
    @test pb.polarization[1] isa SVector{3, ComplexF64}
    # Empty bundles
    @test length(RayBundle(PORT, λ, SVector{3, Float64}[], SVector{3, Float64}[])) == 0
    e = PolarizedRayBundle(PORT, λ, SVector{3, Float64}[], SVector{3, Float64}[], Float64[],
        Float64[], SVector{3, ComplexF64}[])
    @test total_power(e) === 0.0
end

@testset "type stability" begin
    @test (@inferred rays()) isa RayBundle{Float64}
    @test (@inferred polarized()) isa PolarizedRayBundle{Float64}
    @test (@inferred total_power(polarized())) ≈ 3.5e-3
    @test (@inferred length(rays())) == 3
end

@testset "invariants of RayBundle" begin
    # Unequal lengths, wrong number of elements
    @test_throws ArgumentError RayBundle(PORT, λ, POS[1:2], DIR)
    @test_throws ArgumentError RayBundle(PORT, λ, with(collect.(POS), 1, [0.1, 0.2]), DIR)
    # Wavelength
    @test_throws ArgumentError RayBundle(PORT, 0.0, POS, DIR)
    @test_throws ArgumentError RayBundle(PORT, -λ, POS, DIR)
    # Direction: not normalized (no silent normalization), against the normal, grazing
    @test_throws ArgumentError RayBundle(PORT, λ, POS, with(DIR, 2, 2 * DIR[2]))
    @test_throws ArgumentError RayBundle(PORT, λ, POS, with(DIR, 1, -DIR[1]))
    @test_throws ArgumentError RayBundle(PORT, λ, POS, with(DIR, 1, SVector(1.0, 0, 0)))
    # Position off the port plane; rounding-level deviations are accepted
    @test_throws ArgumentError RayBundle(PORT, λ,
        with(POS, 2, POS[2] + SVector(0.0, 0.0, 1e-6)), DIR)
    @test length(RayBundle(PORT, λ, with(POS, 2, POS[2] + SVector(0.0, 0.0, 1e-12)), DIR)) == 3
    # The same geometric checks apply to polarized rays
    @test_throws ArgumentError PolarizedRayBundle(PORT, λ, POS, with(DIR, 1, -DIR[1]), OPL,
        PWR, POL)
end

@testset "invariants of PolarizedRayBundle" begin
    # One entry per ray
    @test_throws ArgumentError PolarizedRayBundle(PORT, λ, POS, DIR, OPL[1:2], PWR, POL)
    @test_throws ArgumentError PolarizedRayBundle(PORT, λ, POS, DIR, OPL, PWR[1:2], POL)
    @test_throws ArgumentError PolarizedRayBundle(PORT, λ, POS, DIR, OPL, PWR, POL[1:2])
    # Power: negative rejected, zero allowed
    @test_throws ArgumentError PolarizedRayBundle(PORT, λ, POS, DIR, OPL,
        with(PWR, 2, -1e-9), POL)
    @test total_power(PolarizedRayBundle(PORT, λ, POS, DIR, OPL, zeros(3), POL)) == 0
    # Polarization: 3 components, unit norm
    @test_throws ArgumentError PolarizedRayBundle(PORT, λ, POS, DIR, OPL, PWR, JONES)
    @test_throws ArgumentError PolarizedRayBundle(PORT, λ, POS, DIR, OPL, PWR,
        with(POL, 1, 2 * POL[1]))
    # Unit but longitudinal field vector on the oblique ray
    @test_throws ArgumentError PolarizedRayBundle(PORT, λ, POS, DIR, OPL, PWR,
        with(POL, 2, complex.(DIR[2])))
    # Transversality is checked without conjugation: e = (x̂ + i d)/√2 has e·d = i/√2
    x̂, _ = OB.ray_basis(PORT, DIR[2])
    @test_throws ArgumentError PolarizedRayBundle(PORT, λ, POS, DIR, OPL, PWR,
        with(POL, 2, (x̂ + im * DIR[2]) / sqrt(2)))
end

@testset "rays in the problem interface" begin
    struct RayTrace <: AbstractPropagationAlgorithm end
    OB.input_representation(::RayTrace) = RayBundle
    OB.output_representation(::RayTrace) = RayBundle
    prob = PropagationProblem(rays(), nothing, PORT)
    @test prob.field isa RayBundle
    @test is_compatible(rays(), RayTrace())
    @test !is_compatible(polarized(), RayTrace())
    @test_throws MissingConverterError check_compatibility(polarized(), RayTrace())
end
