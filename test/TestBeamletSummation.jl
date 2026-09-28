using LinearAlgebra
using StaticArrays

const OB = OpticsBase

const λ = 1.0e-6
const W = 1.0e-3                       # beamlet waist radius
const P0 = 1.0e-3
const PORT = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])

# 1/q of a stigmatic Gaussian at distance z behind its waist, q = z − i z_R (medium index n)
inv_q(z; w0 = W, n = 1.0) = 1 / (z - im * π * n * w0^2 / λ)

# Grid ±8W with spacing W/16
const GRID = RegularGrid((256, 256), (W / 16, W / 16))

function single_beamlet(; d = SVector(0.0, 0.0, 1.0), invq = inv_q(0.0), phasor = 1.0 + 0im,
        port = PORT, power = P0, opl = 0.0)
    return RayBundle(port, λ, [SVector(0.0, 0.0, 0.0)], [d], [opl], [power], [phasor];
        beamlet = [invq])
end

# Independent scalar reference: stigmatic Gaussian from 1/q at the ray position, evaluated
# at global point r for a ray at the origin with direction d.
function reference_field(r, d, invq, a; n = 1.0)
    k = 2π * n / λ
    s = dot(r, d)
    ρ² = sum(abs2, r - s * d)
    q = 1 / invq
    return a / (1 + s / q) * exp(im * k * s) * exp(im * k * ρ² / (2 * (q + s)))
end

peak_amplitude(P; n = 1.0, w = W) = sqrt(2P / (π * w^2) / OB.power_normalization(
    PlanarPort([0, 0, 0], [0, 0, 1], [1, 0, 0]; refractive_index = n)))

@testset "single beamlet at normal incidence" begin
    f = convert_field(GaussianBeamletSummation(GRID), single_beamlet())
    @test f isa SampledField{1}
    @test OB.port(f) == PORT
    @test OB.wavelength(f) == λ
    @test abs(total_power(f) - P0) / P0 < 1e-9
    ξ = OB.coordinates(GRID, 1)
    η = OB.coordinates(GRID, 2)
    E = OB.field_array(f)[:, :, 1]
    E_ref = [peak_amplitude(P0) * exp(-(x^2 + y^2) / W^2) for x in ξ, y in η]
    @test maximum(abs, E - E_ref) / maximum(abs, E_ref) < 1e-12
end

@testset "phase: opl, phasor and curvature" begin
    # Off-waist beamlet (curved wavefront) with a phasor and an optical path length
    invq = inv_q(0.4)
    opl = 0.123456789
    ph = cis(0.7)
    f = convert_field(GaussianBeamletSummation(GRID), single_beamlet(; invq, opl, phasor = ph))
    E = OB.field_array(f)[:, :, 1]
    k = 2π / λ
    a = sqrt(P0 * k * imag(invq) / (OB.power_normalization(PORT) * π)) * cis(k * opl) * ph
    ξ = OB.coordinates(GRID, 1)
    η = OB.coordinates(GRID, 2)
    E_ref = [reference_field(SVector(x, y, 0.0), SVector(0.0, 0.0, 1.0), invq, a)
             for x in ξ, y in η]
    @test maximum(abs, E - E_ref) / maximum(abs, E_ref) < 1e-10
    @test abs(total_power(f) - P0) / P0 < 1e-9
end

@testset "oblique beamlet" begin
    θ = deg2rad(20)
    d = SVector(sin(θ), 0.0, cos(θ))
    invq = inv_q(0.2)
    f = convert_field(GaussianBeamletSummation(GRID), single_beamlet(; d, invq))
    # No obliquity weighting: the port-plane footprint carries P/cos θ
    @test total_power(f) ≈ P0 / cos(θ) rtol = 1e-4
    # Compare with the scalar tilted-Gaussian formula
    k = 2π / λ
    a = sqrt(P0 * k * imag(invq) / (OB.power_normalization(PORT) * π))
    E = OB.field_array(f)[:, :, 1]
    ξ = OB.coordinates(GRID, 1)
    η = OB.coordinates(GRID, 2)
    E_ref = [reference_field(SVector(x, y, 0.0), d, invq, a) for x in ξ, y in η]
    @test maximum(abs, E - E_ref) / maximum(abs, E_ref) < 1e-10
end

@testset "astigmatic beamlet" begin
    # Different waists along x̂ and ŷ, rotated by 30° in the transverse plane
    qx = inv_q(0.0; w0 = W)
    qy = inv_q(0.0; w0 = 0.5W)
    α = deg2rad(30)
    R = SMatrix{2, 2}(cos(α), sin(α), -sin(α), cos(α))
    Q = R * SMatrix{2, 2, ComplexF64}(qx, 0, 0, qy) * transpose(R)
    b = RayBundle(PORT, λ, [SVector(0.0, 0.0, 0.0)], [SVector(0.0, 0.0, 1.0)], [0.0], [P0],
        [1.0 + 0im]; beamlet = [Q])
    f = convert_field(GaussianBeamletSummation(GRID), b)
    @test abs(total_power(f) - P0) / P0 < 1e-9
    E = OB.field_array(f)[:, :, 1]
    ξ = OB.coordinates(GRID, 1)
    η = OB.coordinates(GRID, 2)
    a = sqrt(P0 * (2π / λ) * sqrt(det(imag(Q))) / (OB.power_normalization(PORT) * π))
    E_ref = [a * exp(im * π / λ * transpose(SVector(x, y)) * Q * SVector(x, y)) for x in ξ, y in η]
    @test maximum(abs, E - E_ref) / maximum(abs, E_ref) < 1e-12
end

@testset "coherent summation" begin
    # Opposite phasors at the same position cancel
    b = RayBundle(PORT, λ, fill(SVector(0.0, 0.0, 0.0), 2), fill(SVector(0.0, 0.0, 1.0), 2),
        [0.0, 0.0], [P0, P0], [1.0 + 0im, -1.0 + 0im]; beamlet = fill(inv_q(0.0), 2))
    f = convert_field(GaussianBeamletSummation(GRID), b)
    E1 = OB.field_array(convert_field(GaussianBeamletSummation(GRID), single_beamlet()))
    @test maximum(abs, OB.field_array(f)) < 1e-12 * maximum(abs, E1)

    # Two broad beamlets at ±θ: fringes with period λ/(2 sin θ); zeros at ±λ/(4 sin θ)
    θ = asin(5e-3)                            # sin θ = 0.005 exactly
    ds = [SVector(sin(θ), 0.0, cos(θ)), SVector(-sin(θ), 0.0, cos(θ))]
    b2 = RayBundle(PORT, λ, fill(SVector(0.0, 0.0, 0.0), 2), ds, [0.0, 0.0], [P0, P0],
        [1.0 + 0im, 1.0 + 0im]; beamlet = fill(inv_q(0.0; w0 = 5e-3), 2))
    g = RegularGrid((257, 1), (5e-6, 5e-6))
    f2 = convert_field(GaussianBeamletSummation(g), b2)
    ξ = OB.coordinates(g, 1)
    E = OB.field_array(f2)[:, 1, 1]
    ξ0 = λ / (4 * sin(θ))                      # 50 µm, a grid point
    i(x) = findfirst(y -> isapprox(y, x; atol = 1e-12), ξ)
    Emax = abs(E[i(0.0)])
    for m in (-3, -1, 1, 3)
        @test abs(E[i(m * ξ0)]) < 1e-6 * Emax
    end
    for m in (-2, 2)
        @test abs(E[i(m * ξ0)]) ≈ Emax rtol = 1e-3
    end
end

@testset "vectorial bundle" begin
    d = normalize(SVector(0.1, -0.05, 1.0))
    e = OB.jones_to_global(PORT, d, OB.circular_jones(+1))
    b = RayBundle(PORT, λ, [SVector(0.0, 0.0, 0.0)], [d], [0.0], [P0], [e];
        beamlet = [inv_q(0.1)])
    f = convert_field(GaussianBeamletSummation(GRID), b)
    @test f isa SampledField{3}
    @test is_vectorial(f)
    # Each component is the scalar beamlet times the constant polarization vector
    fs = convert_field(GaussianBeamletSummation(GRID), single_beamlet(; d, invq = inv_q(0.1)))
    Es = OB.field_array(fs)[:, :, 1]
    for c in 1:3
        @test OB.field_array(f)[:, :, c] ≈ Es .* e[c] rtol = 1e-12
    end
    @test total_power(f) ≈ total_power(fs) rtol = 1e-12
end

@testset "rejected inputs" begin
    pure = RayBundle(PORT, λ, [SVector(0.0, 0.0, 0.0)], [SVector(0.0, 0.0, 1.0)], [0.0], [P0],
        [1.0 + 0im])
    @test_throws MissingConverterError convert_field(GaussianBeamletSummation(GRID), pure)
    sf = convert_field(GaussianBeamletSummation(GRID), single_beamlet())
    @test_throws MissingConverterError convert_field(GaussianBeamletSummation(GRID), sf)
    @test !is_compatible(pure, GaussianBeamletSummation(GRID))
    @test is_compatible(single_beamlet(), GaussianBeamletSummation(GRID))
end

@testset "medium index" begin
    # In a medium the peak amplitude drops by 1/√n for the same power (κ = n/(2Z₀))
    n = 1.5
    pn = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]; refractive_index = n)
    b = single_beamlet(; port = pn, invq = inv_q(0.0; n))
    f = convert_field(GaussianBeamletSummation(GRID), b)
    @test abs(total_power(f) - P0) / P0 < 1e-9
    E = OB.field_array(f)[:, :, 1]
    @test maximum(abs, E) ≈ peak_amplitude(P0; n) rtol = 1e-12
end
