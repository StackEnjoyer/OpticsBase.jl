# DebyeWolf: converging rays → PlaneWaveSpectrum (Debye approximation), checked against
# the trace case of the plan, the Richards–Wolf integrals of an aplanatic lens (x and
# radial polarization) and Voronoi solid angles.

using LinearAlgebra
using StaticArrays
const OB = OpticsBase

const λ = 1.0e-6
const n_med = 1.5
const λₘ = λ / n_med
const k = 2π / λₘ
const k0 = 2π / λ
const L = 1.0e-3            # ray port A lies at z = −L, the focus at the origin
const f = 2.0e-3            # focal length (radius of the reference sphere)
const P0 = 1.0e-3
const sinα = 0.9            # aperture half angle in the medium (NA = n sin α = 1.35)
const opl0 = 5.0e-3         # optical path length from the source to the focus

const portF = PlanarPort(SVector(0.0, 0, 0), SVector(0.0, 0, 1), SVector(1.0, 0, 0);
    refractive_index = n_med)
const portA = PlanarPort(SVector(0.0, 0, -L), SVector(0.0, 0, 1), SVector(1.0, 0, 0);
    refractive_index = n_med)
const κ = OB.power_normalization(portF)

# Gauss–Legendre nodes and weights on [−1, 1] (Golub–Welsch)
function gauss_legendre(m)
    β = [i / sqrt(4i^2 - 1) for i in 1:(m - 1)]
    E = eigen(SymTridiagonal(zeros(m), β))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end

# Rays from port A towards the focus along the directions s, in phase at the focus
function converging_rays(s, power, jones; beamlet = nothing)
    positions = [-(L / d[3]) * d for d in s]
    opl = [opl0 - n_med * norm(p) for p in positions]
    phasor = [OB.jones_to_global(portA, d, J) for (d, J) in zip(s, jones)]
    return RayBundle(portA, λ, positions, s, opl, power, phasor; beamlet)
end

direction(θ, φ) = SVector(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))
x_jones(_) = SVector(1.0 + 0im, 0)
radial_jones(φ) = SVector(cos(φ) + 0im, sin(φ))

# Aplanatic lens with a uniform entrance pupil (radius f sin α, power P0): the power per
# solid angle is S₀ f² cos θ. Rays on a Gauss–Legendre (θ) × uniform (φ) grid.
const S0 = P0 / (π * (f * sinα)^2)
function aplanatic_bundle(jones; nθ = 48, nφ = 96)
    x, wx = gauss_legendre(nθ)
    α = asin(sinα)
    s = SVector{3, Float64}[]
    w = Float64[]
    J = SVector{2, ComplexF64}[]
    for (xi, wi) in zip(x, wx), l in 1:nφ
        θ = α / 2 * (xi + 1)
        φ = 2π * (l - 0.5) / nφ
        push!(s, direction(θ, φ))
        push!(w, α / 2 * wi * sin(θ) * 2π / nφ)
        push!(J, jones(φ))
    end
    power = [S0 * f^2 * d[3] * wj for (d, wj) in zip(s, w)]
    return converging_rays(s, power, J), w
end

# Richards–Wolf reference. Bessel functions by the trapezoid rule on
# J_m(x) = (1/2π) ∫₀^{2π} cos(mτ − x sin τ) dτ (exponentially convergent).
besselj_trap(m, x; K = 96) = sum(cos(m * τ - x * sin(τ)) for τ in (2π / K) .* (0:(K - 1))) / K

const A0 = -im / λₘ * sqrt(S0 * f^2 / κ) * cis(k0 * opl0)

function rw_integrals(ρ, z; m = 200)
    x, wx = gauss_legendre(m)
    α = asin(sinα)
    I0 = I1 = I2 = R1 = R0 = zero(ComplexF64)
    for (xi, wi) in zip(x, wx)
        θ = α / 2 * (xi + 1)
        c, s = cos(θ), sin(θ)
        h = α / 2 * wi * sqrt(c) * cis(k * z * c)
        J0, J1, J2 = besselj_trap(0, k * ρ * s), besselj_trap(1, k * ρ * s),
        besselj_trap(2, k * ρ * s)
        I0 += s * (1 + c) * J0 * h
        I1 += s^2 * J1 * h
        I2 += s * (1 - c) * J2 * h
        R1 += c * s * J1 * h
        R0 += s^2 * J0 * h
    end
    return (; I0, I1, I2, R1, R0)
end

function rw_field(ρ, φ, z, pol)
    I = rw_integrals(ρ, z)
    if pol === :x
        return SVector(π * A0 * (I.I0 + I.I2 * cos(2φ)), π * A0 * I.I2 * sin(2φ),
            -2π * im * A0 * I.I1 * cos(φ))
    else
        return SVector(2π * im * A0 * I.R1 * cos(φ), 2π * im * A0 * I.R1 * sin(φ),
            -2π * A0 * I.R0)
    end
end

const grid = RegularGrid((25, 25), (λₘ / 4, λₘ / 4))    # ±3λₘ

function rw_reference(z, pol)
    ξs = OB.coordinates(grid, 1)
    ηs = OB.coordinates(grid, 2)
    E = zeros(ComplexF64, length(ξs), length(ηs), 3)
    for (i, ξ) in enumerate(ξs), (j, η) in enumerate(ηs)
        E[i, j, :] = rw_field(hypot(ξ, η), atan(η, ξ), z, pol)
    end
    return E
end

focal_port(z) = PlanarPort(SVector(0.0, 0, z), SVector(0.0, 0, 1), SVector(1.0, 0, 0);
    refractive_index = n_med)

relerr(E, Eref) = maximum(abs, E - Eref) / maximum(abs, Eref)

@testset "trace case of the plan" begin
    p = SVector(0.5e-3, 0, -1e-3)
    d = normalize(-p)
    e = OB.jones_to_global(portA, d, SVector(1.0, 0))
    bundle = RayBundle(portA, λ, [p], [d], [2e-3], [1e-3], [e])
    pws = convert_field(DebyeWolf(portF; solid_angle = [1e-4]), bundle)
    @test pws isa PlaneWaveSpectrum{3}
    @test OB.port(pws) === portF
    @test pws.direction[1] == d
    @test pws.weight == [1e-4]
    @test e ≈ SVector(2, 0, 1) / sqrt(5)
    expected = -im / λₘ * sqrt(1e-3 / (κ * 1e-4)) * cis(k0 * (2e-3 + n_med * norm(p))) * e
    @test pws.amplitude[1] ≈ expected rtol = 1e-12
    @test norm(pws.amplitude[1]) ≈ 1.0631e8 rtol = 1e-4
    @test total_power(pws) ≈ 1e-3 rtol = 1e-12
end

@testset "scalar rays and beamlets" begin
    s = [direction(0.3, φ) for φ in (0.0, 2.0, 4.0)]
    rays = converging_rays(s, fill(1e-4, 3), fill(SVector(1.0, 0), 3))
    positions = rays.position
    scalar = RayBundle(portA, λ, positions, s, rays.opl, rays.power, fill(1.0 + 0im, 3))
    pws = convert_field(DebyeWolf(portF; solid_angle = fill(1e-3, 3)), scalar)
    @test pws isa PlaneWaveSpectrum{1}
    @test total_power(pws) ≈ 3e-4 rtol = 1e-12
    # Beamlets are converted like their chief rays; Q is ignored
    beamlets = converging_rays(s, fill(1e-4, 3), fill(SVector(1.0, 0), 3);
        beamlet = fill(1e3im, 3))
    a = convert_field(DebyeWolf(portF; solid_angle = fill(1e-3, 3)), beamlets)
    b = convert_field(DebyeWolf(portF; solid_angle = fill(1e-3, 3)), rays)
    @test a.amplitude == b.amplitude
    @test a.weight == b.weight
end

@testset "Richards–Wolf, sin α = 0.9, $pol polarization" for pol in (:x, :radial)
    bundle, w = aplanatic_bundle(pol === :x ? x_jones : radial_jones)
    pws = convert_field(DebyeWolf(portF; solid_angle = w), bundle)
    @test total_power(pws) ≈ total_power(bundle) rtol = 1e-12
    @test total_power(pws) ≈ P0 rtol = 1e-10
    for z in (0.0, λₘ)
        field = convert_field(PlaneWaveSummation(grid, focal_port(z)), pws)
        @test relerr(OB.field_array(field), rw_reference(z, pol)) < 1e-8
    end
    if pol === :radial
        # Strong longitudinal field on the axis
        E = OB.field_array(convert_field(PlaneWaveSummation(grid), pws))
        c = size(E, 1) ÷ 2 + 1
        @test abs(E[c, c, 3]) > 0.9 * maximum(abs, E)
    end
end

@testset "errors" begin
    s = [direction(0.3, φ) for φ in (0.0, 2.0, 4.0)]
    rays = converging_rays(s, fill(1e-4, 3), fill(SVector(1.0, 0), 3))
    w = fill(1e-3, 3)
    # Other medium at the focus port
    port_air = PlanarPort(SVector(0.0, 0, 0), SVector(0.0, 0, 1), SVector(1.0, 0, 0))
    @test_throws ArgumentError convert_field(DebyeWolf(port_air; solid_angle = w), rays)
    # Focus upstream of the rays (diverging bundle)
    behind = focal_port(-2L)
    @test_throws ArgumentError convert_field(DebyeWolf(behind; solid_angle = w), rays)
    # Rays not propagating into the port
    sideways = PlanarPort(SVector(0.0, 0, 0), SVector(1.0, 0, 0), SVector(0.0, 1, 0);
        refractive_index = n_med)
    @test_throws ArgumentError convert_field(DebyeWolf(sideways; solid_angle = w), rays)
    # Invalid explicit solid angles
    @test_throws ArgumentError DebyeWolf(portF; solid_angle = [1e-3, 0.0, 1e-3])
    @test_throws ArgumentError DebyeWolf(portF; solid_angle = [1e-3, Inf, 1e-3])
    @test_throws ArgumentError convert_field(DebyeWolf(portF; solid_angle = [1e-3, 1e-3]),
        rays)
    # Wrong representation
    field = SampledField(zeros(ComplexF64, 4, 4), RegularGrid((4, 4), (1e-6, 1e-6)),
        portF, λ)
    @test_throws MissingConverterError convert_field(DebyeWolf(portF), field)
end

@testset "missing DelaunayTriangulation" begin
    if Base.get_extension(OpticsBase, :OpticsBaseDelaunayTriangulationExt) === nothing
        s = [direction(0.3, φ) for φ in (0.0, 2.0, 4.0)]
        rays = converging_rays(s, fill(1e-4, 3), fill(SVector(1.0, 0), 3))
        err = try
            convert_field(DebyeWolf(portF), rays)
        catch e
            e
        end
        @test err isa MethodError
        @test occursin("using DelaunayTriangulation", sprint(showerror, err))
    end
end

using DelaunayTriangulation

@testset "Voronoi solid angles" begin
    # Regular grid of direction cosines: interior cells Δs², edges Δs²/2, corners Δs²/4
    Δs = 0.05
    g = -0.3:Δs:0.3
    s = [SVector(a, b, sqrt(1 - a^2 - b^2)) for a in g for b in g]
    rays = converging_rays(s, fill(1e-6, length(s)), fill(SVector(1.0, 0), length(s)))
    pws = convert_field(DebyeWolf(portF), rays)
    area = pws.weight .* [d[3] for d in pws.direction]
    interior = [abs(d[1]) < 0.29 && abs(d[2]) < 0.29 for d in s]
    @test all(a -> isapprox(a, Δs^2; rtol = 1e-12), area[interior])
    @test area[1] ≈ Δs^2 / 4 rtol = 1e-12
    @test area[2] ≈ Δs^2 / 2 rtol = 1e-12
    @test sum(area) ≈ 0.6^2 rtol = 1e-12

    # Degenerate direction sets
    @test_throws ArgumentError convert_field(DebyeWolf(portF), converging_rays(s[1:2],
        fill(1e-6, 2), fill(SVector(1.0, 0), 2)))
    dup = [s[1], s[2], s[2], s[20]]
    @test_throws ArgumentError convert_field(DebyeWolf(portF), converging_rays(dup,
        fill(1e-6, 4), fill(SVector(1.0, 0), 4)))
    line = [SVector(a, 0.0, sqrt(1 - a^2)) for a in (-0.2, 0.0, 0.1, 0.3)]
    @test_throws ArgumentError convert_field(DebyeWolf(portF), converging_rays(line,
        fill(1e-6, 4), fill(SVector(1.0, 0), 4)))

    # Fibonacci directions, equal area in (s_u, s_v): equal power for the aplanatic lens
    M = 20_000
    r = sinα .* sqrt.(((0:(M - 1)) .+ 0.5) ./ M)
    ψ = (0:(M - 1)) .* (π * (3 - sqrt(5)))
    s = [SVector(ri * cos(ψi), ri * sin(ψi), sqrt(1 - ri^2)) for (ri, ψi) in zip(r, ψ)]
    rays = converging_rays(s, fill(P0 / M, M), fill(SVector(1.0, 0), M))
    pws = convert_field(DebyeWolf(portF), rays)
    @test total_power(pws) ≈ P0 rtol = 1e-12
    field = convert_field(PlaneWaveSummation(grid), pws)
    @test relerr(OB.field_array(field), rw_reference(0.0, :x)) < 2e-2
end
