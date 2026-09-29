using LinearAlgebra
using StaticArrays

const OB = OpticsBase

const λ = 1.0e-6
const NIDX = 1.5
const λm = λ / NIDX
const K = 2π / λm

# Tilted port (normal not along a global axis), origin off the global origin
const PORT_T = PlanarPort([1.0e-6, -2.0e-6, 0.5e-6], normalize([0.2, -0.3, 1.0]),
    [1.0, 0.0, 0.0]; refractive_index = NIDX)
# Port with normal z and u rotated by 30° about z (E_x has both E_u and E_v parts)
const PORT_Z = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [cosd(30), sind(30), 0.0];
    refractive_index = NIDX)

axes_of(p) = (SVector{3}(OB.local_axes(p)[:, 1]), SVector{3}(OB.local_axes(p)[:, 2]),
    SVector{3}(OB.local_axes(p)[:, 3]))

# Field array of size (size(g)..., N) with component c = vec[c] · f(ξ, η)
function field_on_grid(f, g, vec)
    ξs = OB.coordinates(g, 1)
    ηs = OB.coordinates(g, 2)
    return [vec[c] * f(ξ, η) for ξ in ξs, η in ηs, c in eachindex(vec)]
end

gaussian(w0) = (ξ, η) -> complex(exp(-(ξ^2 + η^2) / w0^2))

@testset "missing FFTW backend" begin
    # Runs only if no other test module has loaded FFTW before.
    if Base.get_extension(OpticsBase, :OpticsBaseFFTWExt) === nothing
        g = RegularGrid((8, 8), (λm / 4, λm / 4))
        f = SampledField(ones(ComplexF64, 8, 8, 1), g, PORT_Z, λ)
        e = try
            convert_field(PlaneWaveDecomposition(), f)
        catch err
            err
        end
        @test e isa MethodError
        @test e.f === OB._centered_fft
        @test occursin("using FFTW", sprint(showerror, e))
    end
end

using FFTW

@testset "construction, traits and errors" begin
    @test PlaneWaveDecomposition().pad_factor == 1
    @test PlaneWaveDecomposition(pad_factor = 2).pad_factor == 2
    @test PlaneWaveDecomposition() isa OB.AbstractFieldConverter
    @test_throws ArgumentError PlaneWaveDecomposition(pad_factor = 0)
    @test_throws ArgumentError PlaneWaveDecomposition(-1)
    conv = PlaneWaveDecomposition()
    @test input_representation(conv) == SampledField{<:Any, 2}
    @test output_representation(conv) == PlaneWaveSpectrum
    # 3D fields are not accepted
    g3 = RegularGrid((4, 4, 3), (λm, λm, λm))
    f3 = SampledField(ones(ComplexF64, 4, 4, 3, 1), g3, PORT_Z, λ)
    @test_throws MissingConverterError convert_field(conv, f3)
end

@testset "plane wave on an FFT bin" begin
    a, b = 3, -2                      # bin indices on the unpadded grid
    for dims in ((16, 12), (15, 13)), pf in (1, 2), N in (1, 3)
        p = PORT_T
        u, v, n = axes_of(p)
        Δ = (λm / 4, λm / 3)
        g = RegularGrid(dims, Δ)
        kξ = a * 2π / (dims[1] * Δ[1])
        kη = b * 2π / (dims[2] * Δ[2])
        kn = sqrt(K^2 - kξ^2 - kη^2)
        s = (kξ * u + kη * v + kn * n) / K
        A = N == 1 ? SVector(2.0e3 * cis(0.4)) :
            2.0e3 * OB.jones_to_global(p, s, SVector(0.6 + 0.2im, -0.5im))
        E = field_on_grid((ξ, η) -> cis(kξ * ξ + kη * η), g, A)
        pws = convert_field(PlaneWaveDecomposition(pad_factor = pf), SampledField(E, g, p, λ))
        @test pws isa PlaneWaveSpectrum{N}
        @test OB.port(pws) === p
        @test OB.wavelength(pws) == λ
        # Sample indices on the padded wave-vector grid
        Δkξ = 2π / (pf * dims[1] * Δ[1])
        Δkη = 2π / (pf * dims[2] * Δ[2])
        npeak = 0
        maxother = 0.0
        for j in 1:length(pws)
            sj = pws.direction[j]
            mξ = round(Int, K * dot(sj, u) / Δkξ)
            mη = round(Int, K * dot(sj, v) / Δkη)
            wℰ = pws.weight[j] * pws.amplitude[j]
            if (mξ, mη) == (pf * a, pf * b)
                npeak += 1
                # Zero padding by pf spreads the plane wave: the peak carries A/pf²
                @test norm(wℰ - A / pf^2) < 1e-12 * norm(A)
                @test norm(sj - s) < 1e-14
            elseif mξ % pf == 0 && mη % pf == 0
                # All other bins of the unpadded grid are empty (all bins for pf = 1)
                maxother = max(maxother, norm(wℰ))
            end
        end
        @test npeak == 1
        @test maxother < 1e-12 * norm(A)
    end
end

@testset "round trip, scalar Gaussian" begin
    # w₀ = 4λₘ, Δ = λₘ/4; the grid reaches ±6w₀ so that the truncated Gaussian has no
    # spectral content at k_t ≥ k (it would be dropped as evanescent).
    g = RegularGrid((192, 192), (λm / 4, λm / 4))
    E = field_on_grid(gaussian(4λm), g, SVector(1.0e3))
    f = SampledField(E, g, PORT_T, λ)
    for pf in (1, 2)
        pws = convert_field(PlaneWaveDecomposition(pad_factor = pf), f)
        f2 = convert_field(PlaneWaveSummation(g), pws)
        @test OB.port(f2) === PORT_T
        @test maximum(abs, OB.field_array(f2) - E) < 1e-10 * maximum(abs, E)
    end
end

@testset "round trip, vectorial E_x" begin
    g = RegularGrid((192, 192), (λm / 4, λm / 4))
    E = field_on_grid(gaussian(4λm), g, SVector(1.0e3, 0.0, 0.0))
    f = SampledField(E, g, PORT_Z, λ)
    pws = convert_field(PlaneWaveDecomposition(), f)
    # The constructor checked transversality; check it again explicitly
    @test all(abs(transpose(pws.amplitude[j]) * pws.direction[j]) ≤
              1e-12 * norm(pws.amplitude[j]) for j in 1:length(pws))
    f2 = convert_field(PlaneWaveSummation(g), pws)
    E2 = OB.field_array(f2)
    Emax = maximum(abs, E)
    @test maximum(abs, E2[:, :, 1] - E[:, :, 1]) < 1e-10 * Emax
    @test maximum(abs, E2[:, :, 2]) < 1e-10 * Emax
    # E_z is recomputed from transversality (nonzero for a finite beam)
    @test maximum(abs, E2[:, :, 3]) > 1e-3 * Emax
    # Idempotent: decomposing the round-trip field gives the same spectrum
    pws2 = convert_field(PlaneWaveDecomposition(), f2)
    @test pws2.direction == pws.direction
    @test pws2.weight == pws.weight
    wℰ = pws.weight .* pws.amplitude
    wℰ2 = pws2.weight .* pws2.amplitude
    @test maximum(norm, wℰ2 - wℰ) < 1e-10 * maximum(norm, wℰ)
end

@testset "power of a paraxial Gaussian" begin
    # w₀ = 10λₘ: the exact flux differs from the paraxial power by ≈ 1/(k w₀)² ≈ 2.5e-4
    g = RegularGrid((160, 160), (λm / 2, λm / 2))
    for N in (1, 3)
        vec = N == 1 ? SVector(1.0e3) : SVector(1.0e3, 0.5e3im, 0.0)
        f = SampledField(field_on_grid(gaussian(10λm), g, vec), g, PORT_Z, λ)
        pws = convert_field(PlaneWaveDecomposition(), f)
        @test total_power(pws) ≈ total_power(f) rtol = 1e-3
    end
end

@testset "evanescent samples are dropped" begin
    # Δ = λₘ/4: k_ξ/k = 4m/M_ξ; propagating iff (4m_ξ/M_ξ)² + (4m_η/M_η)² < 1
    dims = (64, 48)
    g = RegularGrid(dims, (λm / 4, λm / 4))
    f = SampledField(field_on_grid(gaussian(4λm), g, SVector(1.0 + 0im)), g, PORT_Z, λ)
    for pf in (1, 2)
        M = pf .* dims
        count = 0
        for mξ in (-(M[1] ÷ 2)):(M[1] - M[1] ÷ 2 - 1), mη in (-(M[2] ÷ 2)):(M[2] - M[2] ÷ 2 - 1)
            count += (4mξ)^2 * M[2]^2 + (4mη)^2 * M[1]^2 < M[1]^2 * M[2]^2
        end
        pws = convert_field(PlaneWaveDecomposition(pad_factor = pf), f)
        @test length(pws) == count
        @test all(d -> dot(d, SVector(0.0, 0.0, 1.0)) > 0, pws.direction)
    end
end
