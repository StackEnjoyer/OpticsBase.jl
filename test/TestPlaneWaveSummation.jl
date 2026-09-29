using LinearAlgebra
using StaticArrays

const OB = OpticsBase

const λ = 1.0e-6
const NIDX = 1.5
const λm = λ / NIDX
const K = 2π / λm
const R0 = SVector(1.0e-6, 2.0e-6, -1.0e-6)
const PORT = PlanarPort(R0, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]; refractive_index = NIDX)

# Unit direction at polar angle θ and azimuth φ about the port normal (z)
oblique(θ, φ) = SVector(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))

# Output port tilted by 20° against the spectrum port normal, shifted from r₀
const TILT = deg2rad(20)
const PORT_OUT = PlanarPort(R0 + SVector(0.5e-6, -0.3e-6, 2.0e-6),
    [sin(TILT), 0.0, cos(TILT)], [0.0, 1.0, 0.0]; refractive_index = NIDX)

const GRID2 = RegularGrid((9, 8), (0.3λ, 0.25λ))
const GRID3 = RegularGrid((7, 6, 5), (0.3λ, 0.25λ, 0.4λ))

# Global position of grid index I at port p
function grid_point(p, g, I)
    D = length(size(g))
    ξ = ntuple(d -> OB.coordinates(g, d)[I[d]], D)
    return D == 2 ? OB.to_global(p, SVector(ξ[1], ξ[2])) : OB.to_global(p, SVector(ξ))
end

# Naive reference: every grid point, every plane wave (independent of the implementation)
function naive_sum(pws, p, g)
    N = is_vectorial(pws) ? 3 : 1
    r0 = OB.origin(OB.port(pws))
    E = zeros(ComplexF64, size(g)..., N)
    for I in CartesianIndices(size(g))
        r = grid_point(p, g, Tuple(I))
        acc = zero(SVector{N, ComplexF64})
        for j in 1:length(pws)
            acc += pws.weight[j] * pws.amplitude[j] *
                   exp(im * K * dot(pws.direction[j], r - r0))
        end
        for c in 1:N
            E[I, c] = acc[c]
        end
    end
    return E
end

# Deterministic quasi-random directions in a cone of half angle θmax and transverse
# amplitudes (no Random dependency)
frac(x) = x - floor(x)
function quasi_random_spectrum(M, N; θmax = deg2rad(40))
    dirs = [oblique(acos(1 - (1 - cos(θmax)) * frac(j * 0.6180339887498949)),
                2π * frac(j * 0.7548776662466927)) for j in 1:M]
    jones = [SVector(cis(1.3j), 0.5 * cis(2.1j)) for j in 1:M]
    amp = N == 3 ? [1e6 * OB.jones_to_global(PORT, dirs[j], jones[j]) for j in 1:M] :
          [1e6 * jones[j][1] for j in 1:M]
    w = [1e-4 * (1 + frac(j * 0.414)) for j in 1:M]
    return PlaneWaveSpectrum(PORT, λ, dirs, amp, w)
end

@testset "traits and construction" begin
    conv = PlaneWaveSummation(GRID2)
    @test conv isa OB.AbstractFieldConverter
    @test conv.port === nothing
    @test input_representation(conv) == PlaneWaveSpectrum
    @test output_representation(conv) == SampledField
    @test PlaneWaveSummation(GRID3, PORT_OUT).port === PORT_OUT
    @test_throws ArgumentError PlaneWaveSummation(RegularGrid((4,), (1e-6,)))
end

@testset "single plane wave on a tilted port" begin
    s = oblique(deg2rad(25), deg2rad(40))
    w = 2.5e-4
    ℰ3 = 3e6 * OB.jones_to_global(PORT, s, SVector(0.8 + 0.1im, -0.4im))
    ℰ1 = 3e6 * cis(0.7)
    for (ℰ, N) in ((ℰ1, 1), (ℰ3, 3)), g in (GRID2, GRID3)
        pws = PlaneWaveSpectrum(PORT, λ, [s], [ℰ], [w])
        f = convert_field(PlaneWaveSummation(g, PORT_OUT), pws)
        @test f isa SampledField{N}
        @test OB.port(f) === PORT_OUT
        @test OB.grid(f) === g
        @test OB.wavelength(f) == λ
        E = OB.field_array(f)
        @test size(E) == (size(g)..., N)
        D = length(size(g))
        maxrel = 0.0
        for I in CartesianIndices(size(g))
            r = grid_point(PORT_OUT, g, Tuple(I))
            ref = SVector{N}(w * ℰ * exp(im * K * dot(s, r - R0)))
            val = SVector{N}(ntuple(c -> E[Tuple(I)..., c], N))
            maxrel = max(maxrel, norm(val - ref) / norm(ref))
        end
        @test maxrel < 1e-12
    end
end

@testset "default port is the spectrum port" begin
    pws = quasi_random_spectrum(20, 3)
    f = convert_field(PlaneWaveSummation(GRID2), pws)
    @test OB.port(f) === PORT
    E_ref = naive_sum(pws, PORT, GRID2)
    @test maximum(abs, OB.field_array(f) - E_ref) < 1e-12 * maximum(abs, E_ref)
    # At the grid center (the port origin r₀) the field is Σ w ℰ
    c = size(GRID2) .÷ 2 .+ 1
    Σ = sum(pws.weight .* pws.amplitude)
    @test SVector{3}(OB.field_array(f)[c..., :]) ≈ Σ rtol = 1e-13
end

@testset "chunk independence" begin
    M = 3 * OB._PLANE_WAVE_CHUNK + 1
    g = RegularGrid((5, 4, 3), (0.35λ, 0.3λ, 0.5λ))
    for N in (1, 3)
        pws = quasi_random_spectrum(M, N)
        f = convert_field(PlaneWaveSummation(g, PORT_OUT), pws)
        E = OB.field_array(f)
        E_ref = naive_sum(pws, PORT_OUT, g)
        @test maximum(abs, E - E_ref) < 1e-12 * maximum(abs, E_ref)
        # Other chunk sizes give the same result
        for chunk in (1, 7, M, 2M)
            E_c = OB.field_array(OB._plane_wave_sum(g, PORT_OUT, pws, chunk))
            @test maximum(abs, E_c - E_ref) < 1e-12 * maximum(abs, E_ref)
        end
    end
end

@testset "analytic Gaussian spectrum" begin
    port = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
    k = 2π / λ
    w0 = 20λ
    E0 = 1e3
    Δk = 0.1 / w0
    kmax = 8 / w0
    m = round(Int, kmax / Δk)
    dirs = SVector{3, Float64}[]
    amp = ComplexF64[]
    wgt = Float64[]
    for iu in (-m):m, iv in (-m):m
        ku, kv = iu * Δk, iv * Δk
        kt2 = ku^2 + kv^2
        kt2 <= kmax^2 || continue
        kn = sqrt(k^2 - kt2)
        Ẽ = E0 * w0^2 / (4π) * exp(-kt2 * w0^2 / 4)
        push!(dirs, SVector(ku, kv, kn) / k)
        push!(amp, k * kn * Ẽ)
        push!(wgt, Δk^2 / (k * kn))
    end
    pws = PlaneWaveSpectrum(port, λ, dirs, amp, wgt)
    g = RegularGrid((21, 21), (0.2w0, 0.2w0))
    f = convert_field(PlaneWaveSummation(g), pws)
    E = OB.field_array(f)[:, :, 1]
    ξ = OB.coordinates(g, 1)
    η = OB.coordinates(g, 2)
    E_ref = [E0 * exp(-(x^2 + y^2) / w0^2) for x in ξ, y in η]
    @test maximum(abs, E - E_ref) < 1e-6 * E0
    κ = OB.power_normalization(port)
    @test total_power(pws) ≈ κ * π * w0^2 * E0^2 / 2 rtol = 1e-3
end

@testset "element type" begin
    s = SVector{3, Float32}(oblique(deg2rad(25), deg2rad(40)))
    pws = PlaneWaveSpectrum(PORT, Float32(λ), [s], [ComplexF32(1e6)], [1.0f-4])
    f = convert_field(PlaneWaveSummation(GRID2), pws)
    @test eltype(OB.field_array(f)) == ComplexF32
    empty = PlaneWaveSpectrum(PORT, λ, SVector{3, Float64}[], SVector{3, ComplexF64}[],
        Float64[])
    fe = convert_field(PlaneWaveSummation(GRID3), empty)
    @test size(OB.field_array(fe)) == (size(GRID3)..., 3)
    @test all(iszero, OB.field_array(fe))
end

struct DummyPort <: AbstractPort end

@testset "errors" begin
    pws = quasi_random_spectrum(5, 3)
    # Other refractive index at the output port
    p_n = PlanarPort(R0, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]; refractive_index = 1.0)
    @test_throws ArgumentError convert_field(PlaneWaveSummation(GRID2, p_n), pws)
    # Directions not propagating into the output port (s·n_out ≤ 0)
    p_back = PlanarPort(R0, [1.0, 0.0, -1.0], [0.0, 1.0, 0.0]; refractive_index = NIDX)
    @test_throws ArgumentError convert_field(PlaneWaveSummation(GRID2, p_back), pws)
    s_side = oblique(deg2rad(60), 0.0)
    p_side = PlanarPort(R0, [-sin(deg2rad(35)), 0.0, cos(deg2rad(35))], [0.0, 1.0, 0.0];
        refractive_index = NIDX)                       # s_side·n < 0
    pw_side = PlaneWaveSpectrum(PORT, λ, [s_side], [1.0 + 0im], [1e-4])
    @test_throws ArgumentError convert_field(PlaneWaveSummation(GRID2, p_side), pw_side)
    # Non-planar output port
    @test_throws ArgumentError convert_field(PlaneWaveSummation(GRID2, DummyPort()), pws)
    # Non-spectrum input
    sf = convert_field(PlaneWaveSummation(GRID2), pws)
    @test_throws MissingConverterError convert_field(PlaneWaveSummation(GRID2), sf)
    @test !is_compatible(sf, PlaneWaveSummation(GRID2))
    @test is_compatible(pws, PlaneWaveSummation(GRID2))
end
