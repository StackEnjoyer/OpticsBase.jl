# Second chain: BeamletOptics polarized rays → 90° off-axis parabolic mirror → Detector →
# RayBundle → DebyeWolf (Voronoi solid angles) → PlaneWaveSpectrum → PlaneWaveSummation
# → vectorial focal field, checked against a reference from analytic mirror rays with exact
# solid angles.

import BeamletOptics as BMO
using DelaunayTriangulation
using LinearAlgebra
using StaticArrays
const OB = OpticsBase

const λ = 1.064e-6
const rfl = 25e-3               # reflected focal length (aperture center → focus)
const D = 20e-3                 # diameter of the collimated input beam
const f_p = rfl / 2             # parent focal length of the 90° OAP
const P0 = 1e-3
const E_in = SVector(0.0, 0, 1) # input polarization (perpendicular to the deflection plane)
const y0 = f_p - 30e-3          # start plane of the input rays
const x_det = 0.4 * rfl         # detector plane between mirror and focus, clear of the input
const F = SVector(0.0, f_p, 0)  # focus
const M = 20_000

# Focus port: normal along the central reflected ray (−x), u along the input polarization
const portF = PlanarPort(F, SVector(-1.0, 0, 0), SVector(0.0, 0, 1))
const grid = RegularGrid((64, 64), (3λ / 32, 3λ / 32))      # ±3λ around the focus

# --- BeamletOptics chain ------------------------------------------------------------------

function traced_detector()
    oap = BMO.OffAxisParabolicMirror(rfl, D; angle = 90)
    BMO.translate_to3d!(oap, [rfl, f_p, 0])          # focus at F, aperture center at (rfl, f_p, 0)
    det = BMO.Detector(25e-3)
    BMO.zrotate3d!(det, deg2rad(90))
    BMO.translate_to3d!(det, [x_det, f_p, 0])
    # Fibonacci (equal-area) sampling of the input disc
    beams = map(0:(M - 1)) do j
        r = D / 2 * sqrt((j + 0.5) / M)
        ψ = j * π * (3 - sqrt(5))
        BMO.Beam([rfl + r * cos(ψ), y0, r * sin(ψ)], [0.0, 1, 0], λ, Vector(E_in))
    end
    source = BMO.CollimatedSource(beams, D, [rfl, y0, 0.0], [0.0, 1, 0])
    BMO.solve_system!(BMO.StaticSystem([oap, det]), source)
    return det
end

bundle = RayBundle(traced_detector(); power = P0)
pws = convert_field(DebyeWolf(portF), bundle)
field = convert_field(PlaneWaveSummation(grid), pws)
E = OB.field_array(field)

# --- Reference: analytic mirror rays, Gauss–Legendre pupil, exact solid angles ------------

function gauss_legendre(m)
    β = [i / sqrt(4i^2 - 1) for i in 1:(m - 1)]
    E = eigen(SymTridiagonal(zeros(m), β))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end

# BMO's ideal mirror: Jones matrix diag(−1, 1) in the (s, p) bases of the in and out rays
function reflect_polarization(e, d_in, d_out)
    s = normalize(cross(d_in, d_out))
    p1 = cross(d_in, s)
    p2 = cross(d_out, s)
    return -dot(e, s) * s + dot(e, p1) * p2
end

function reference_bundle(portA; nr = 64, nφ = 128)
    x, wx = gauss_legendre(nr)
    S0 = P0 / (π * (D / 2)^2)
    C = f_p + rfl                         # paraboloid: ‖q − F‖ = C − q_y
    positions, directions, opl, power, phasor, w = SVector{3, Float64}[],
    SVector{3, Float64}[], Float64[], Float64[], SVector{3, ComplexF64}[], Float64[]
    for (xi, wi) in zip(x, wx), l in 1:nφ
        r = D / 4 * (xi + 1)
        φ = 2π * (l - 0.5) / nφ
        dA = D / 4 * wi * r * 2π / nφ
        xin, zin = rfl + r * cos(φ), r * sin(φ)
        t = (C^2 - f_p^2 - xin^2 - zin^2) / (2 * (C - f_p))
        q = SVector(xin, t, zin)
        ρq = norm(F - q)
        d = (F - q) / ρq
        p = q + ((x_det - q[1]) / d[1]) * d
        push!(positions, p)
        push!(directions, d)
        push!(opl, (t - y0) + norm(p - q))
        push!(power, S0 * dA)
        push!(phasor, reflect_polarization(E_in, SVector(0.0, 1, 0), d))
        push!(w, dA / ρq^2)
    end
    return RayBundle(portA, λ, positions, directions, opl, power, phasor), w
end

ref_bundle, w_ref = reference_bundle(OB.port(bundle))
E_ref = OB.field_array(convert_field(PlaneWaveSummation(grid),
    convert_field(DebyeWolf(portF; solid_angle = w_ref), ref_bundle)))

@testset "BMO rays meet in the focus" begin
    @test length(bundle) == M
    @test OB.normal(OB.port(bundle)) ≈ SVector(-1.0, 0, 0)
    miss = maximum(zip(bundle.position, bundle.direction)) do (p, d)
        Δ = F - p
        norm(Δ - dot(Δ, d) * d)
    end
    @test miss < 1e-6 * rfl
    to_focus = [o + norm(F - p) for (o, p) in zip(bundle.opl, bundle.position)]
    @test maximum(to_focus) - minimum(to_focus) < λ / 1000
    # Polarization from BMO agrees with the ideal-mirror formula of the reference
    e_ref = [reflect_polarization(E_in, SVector(0.0, 1, 0), d) for d in bundle.direction]
    @test maximum(norm.(bundle.phasor .- e_ref)) < 1e-9
end

@testset "focal field" begin
    @test total_power(pws) ≈ P0 rtol = 1e-12
    @test maximum(abs, E - E_ref) / maximum(abs, E_ref) < 2e-2
    # Transverse field along u dominates; a longitudinal component (along −x) exists
    @test maximum(abs, E[:, :, 3]) > maximum(abs, E[:, :, 2])
    @test maximum(abs, E[:, :, 1]) > 1e-2 * maximum(abs, E[:, :, 3])
end
