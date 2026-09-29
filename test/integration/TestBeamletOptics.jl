# Integration tests of OpticsBaseBeamletOpticsExt: BMO detector hits → RayBundle,
# PolarizedRayBundle, SampledField and BMO beams from a RayBundle.
# Evaluated in its own module with `using OpticsBase, Test` (see runtests.jl).

import BeamletOptics as BMO
using LinearAlgebra
using StaticArrays

const λ = 1.064e-6
const w0 = 0.5e-3
const P0 = 1e-3
const zdet = 0.2
const NG = 64
const ΔG = w0 / 8

# Detector at distance `z` along +y, optionally rotated by `tilt` about BMO's z axis.
function detector_at(z; tilt = 0.0)
    det = BMO.Detector(0.05)
    tilt == 0 || BMO.zrotate3d!(det, tilt)
    BMO.translate3d!(det, [0, z, 0])
    return det
end

function trace!(det, beams...)
    system = BMO.System([det])
    for b in beams
        BMO.solve_system!(system, b)
    end
    return det
end

# BMO's field on the grid of `grid` (its z axis is our reversed η axis)
function bmo_field(det, grid)
    ξ = collect(OpticsBase.coordinates(grid, 1))
    η = collect(OpticsBase.coordinates(grid, 2))
    _, _, E = BMO.electric_field(det; n = length(ξ), x_min = first(ξ), x_max = last(ξ),
        z_min = -last(η), z_max = -first(η), progress = false)
    return E[:, end:-1:1]
end

const grid = RegularGrid((NG, NG), (ΔG, ΔG))

agb_at(; kw...) = BMO.AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0,
    support = [1.0, 0, 0], kw...)

@testset "RayBundle from beamlet hits" begin
    agb = agb_at()
    det = trace!(detector_at(zdet), agb)
    b = RayBundle(det)
    @test b isa RayBundle
    @test length(b) == 1
    @test OpticsBase.wavelength(b) == λ
    @test b.position[1] ≈ SVector(0.0, zdet, 0.0) atol = 1e-15
    @test b.direction[1] ≈ SVector(0.0, 1.0, 0.0)
    # stigmatic beamlets give the same geometry
    gb = BMO.GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0, support = [1.0, 0, 0])
    bg = RayBundle(trace!(detector_at(zdet), gb))
    @test bg.position ≈ b.position
    @test bg.direction ≈ b.direction
end

@testset "SampledField: single astigmatic beamlet" begin
    det = trace!(detector_at(zdet), agb_at())
    f = SampledField(det, grid)
    @test f isa SampledField{3}
    @test size(OpticsBase.field_array(f)) == (NG, NG, 3)
    @test total_power(f) ≈ P0 rtol = 1e-6
    E = OpticsBase.field_array(f)
    E_bmo = bmo_field(det, grid)
    normE = sqrt.(sum(abs2, E; dims = 3))[:, :, 1]
    @test maximum(abs.(normE .- abs.(E_bmo))) <= 1e-12 * maximum(abs, E_bmo)
    # Polarization along the support axis (port u = +x), scalar part equals BMO's
    @test maximum(abs, E[:, :, 2:3]) <= 1e-12 * maximum(abs, E_bmo)
    @test maximum(abs, E[:, :, 1] - E_bmo) <= 1e-12 * maximum(abs, E_bmo)
    # Direction and value agree with BMO's vector field at a grid point
    agb = agb_at()
    det2 = trace!(detector_at(zdet), agb)
    i, j = NG ÷ 2 + 3, NG ÷ 2 - 2
    ξ = OpticsBase.coordinates(grid, 1)[i]
    η = OpticsBase.coordinates(grid, 2)[j]
    Ep = BMO.polarized_field(agb, [ξ, 0.0, -η], zdet)
    @test maximum(abs.(collect(Ep) .- E[i, j, :])) <= 1e-6 * maximum(abs, E_bmo)
end

@testset "SampledField: tilted detector, tilted polarization, port shift" begin
    tilt = deg2rad(20)
    agb = agb_at()
    det = trace!(detector_at(zdet; tilt), agb)
    f = SampledField(det, grid)
    @test total_power(f) ≈ BMO.optical_power(agb) rtol = 1e-6
    # A port shifted within the detector plane samples BMO's field at the shifted points
    p0 = OpticsBase.port(f)
    u, v = OpticsBase.local_axes(p0)[:, 1], OpticsBase.local_axes(p0)[:, 2]
    p1 = PlanarPort(OpticsBase.origin(p0) + 3ΔG * u - 2ΔG * v, OpticsBase.normal(p0), u)
    f1 = SampledField(det, grid; port = p1)
    E0 = OpticsBase.field_array(f)
    E1 = OpticsBase.field_array(f1)
    @test maximum(abs, E1[1:(end - 3), 3:end, :] - E0[4:end, 1:(end - 2), :]) <=
          1e-9 * maximum(abs, E0)
    @test OpticsBase.port(f1) === p1
end

@testset "SampledField: stigmatic GaussianBeamlet" begin
    gb = BMO.GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0, support = [1.0, 0, 0])
    det = trace!(detector_at(zdet), gb)
    f = SampledField(det, grid)
    @test f isa SampledField{1}
    @test total_power(f) ≈ P0 rtol = 1e-6
    E_bmo = bmo_field(det, grid)
    @test maximum(abs, OpticsBase.field_array(f)[:, :, 1] - E_bmo) <=
          1e-12 * maximum(abs, E_bmo)
    # same scalar field as the astigmatic beamlet
    fa = SampledField(trace!(detector_at(zdet), agb_at()), grid)
    @test maximum(abs, OpticsBase.field_array(f)[:, :, 1] - OpticsBase.field_array(fa)[:, :, 1]) <=
          1e-5 * maximum(abs, OpticsBase.field_array(f))
end

@testset "Pure rays" begin
    rays_at(x) = BMO.Beam(BMO.Ray([x, 0.0, 0.0], [0.0, 1.0, 0.0], λ))
    det = trace!(detector_at(zdet), rays_at(0.0), rays_at(1e-3))
    b = RayBundle(det)
    @test length(b) == 2
    @test b.position[2] ≈ SVector(1e-3, zdet, 0.0)
    @test_throws ArgumentError PolarizedRayBundle(det; power = 1e-3)
    @test_throws ArgumentError SampledField(det, grid)

    E0 = [0.0, 0.0, 2.0im]
    pol_at(x) = BMO.Beam(BMO.PolarizedRay([x, 0.0, 0.0], [0.0, 1.0, 0.0], λ, E0))
    det = trace!(detector_at(zdet), pol_at(0.0), pol_at(1e-3))
    b = PolarizedRayBundle(det; power = 3e-3)
    @test b isa PolarizedRayBundle
    @test total_power(b) ≈ 3e-3
    @test b.power ≈ [1.5e-3, 1.5e-3]
    @test all(e -> e ≈ SVector(0.0, 0.0, 1.0im), b.polarization)
    @test b.opl ≈ [zdet, zdet]
    @test_throws ArgumentError PolarizedRayBundle(det; power = -1.0)
    @test_throws ArgumentError SampledField(det, grid)
end

@testset "Port from the detector" begin
    # Right-handed port with d·n > 0, for a straight and a tilted detector
    for tilt in (0.0, deg2rad(20), deg2rad(-35))
        det = trace!(detector_at(zdet; tilt), agb_at())
        b = RayBundle(det)
        p = OpticsBase.port(b)
        A = OpticsBase.local_axes(p)
        R = BMO.orientation(det)
        @test LinearAlgebra.det(A) ≈ 1
        @test cross(A[:, 1], A[:, 2]) ≈ A[:, 3]
        @test OpticsBase.normal(p) ≈ -R[:, 2]
        @test A[:, 1] ≈ -R[:, 1]
        @test OpticsBase.origin(p) ≈ BMO.position(det)
        @test OpticsBase.refractive_index(p) == 1
        n = OpticsBase.normal(p)
        @test all(d -> dot(d, n) > 0, b.direction)
        @test dot(b.direction[1], n) ≈ cos(tilt)
    end

    # Pure rays on a tilted detector
    det = trace!(detector_at(zdet; tilt = deg2rad(30)),
        BMO.Beam(BMO.Ray([0.0, 0, 0], [0.0, 1.0, 0.0], λ)),
        BMO.Beam(BMO.Ray([1e-3, 0, 0], [0.0, 1.0, 0.0], λ)))
    b = RayBundle(det)
    n = OpticsBase.normal(OpticsBase.port(b))
    @test all(d -> dot(d, n) > 0, b.direction)
    @test all(p -> abs(dot(p - OpticsBase.origin(OpticsBase.port(b)), n)) < 1e-15,
        b.position)
end

@testset "Port keyword and errors" begin
    det = trace!(detector_at(zdet), agb_at())
    p = PlanarPort([1e-3, zdet, 2e-3], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])
    b = RayBundle(det; port = p)
    @test OpticsBase.port(b) === p
    @test_throws ArgumentError RayBundle(det;
        port = PlanarPort([0.0, zdet, 0], [0.0, -1.0, 0.0], [1.0, 0, 0]))
    @test_throws ArgumentError RayBundle(det;
        port = PlanarPort([0.0, zdet + 1e-3, 0], [0.0, 1.0, 0.0], [1.0, 0, 0]))
    @test_throws ArgumentError RayBundle(det;
        port = PlanarPort([0.0, zdet, 0], [0.0, 1.0, 0.0], [1.0, 0, 0];
            refractive_index = 1.5))
    @test_throws ArgumentError RayBundle(det; port = :notaport)
    # SampledField needs the detector axes
    @test_throws ArgumentError SampledField(det, grid; port = p)
    # non-square grid, unequal spacing
    @test_throws ArgumentError SampledField(det, RegularGrid((NG, NG + 2), (ΔG, ΔG)))
    @test_throws ArgumentError SampledField(det, RegularGrid((NG, NG), (ΔG, 2ΔG)))

    @test_throws ArgumentError RayBundle(detector_at(zdet))  # no hits
    @test_throws ArgumentError SampledField(detector_at(zdet), grid)

    # Hits with different wavelengths
    det = trace!(detector_at(zdet),
        BMO.Beam(BMO.Ray([0.0, 0, 0], [0.0, 1.0, 0.0], λ)),
        BMO.Beam(BMO.Ray([1e-3, 0, 0], [0.0, 1.0, 0.0], 2λ)))
    @test_throws ArgumentError RayBundle(det)
end

@testset "BMO → RayBundle → BMO rays → RayBundle" begin
    xs = range(-2e-3, 2e-3; length = 7)
    ds = [normalize([0.02 * x / 1e-3, 1.0, -0.01 * x / 1e-3]) for x in xs]
    beams = [BMO.Beam([x, 0.0, 0.3x], d, λ) for (x, d) in zip(xs, ds)]
    det = trace!(detector_at(zdet; tilt = deg2rad(10)), beams...)
    b1 = RayBundle(det)
    @test length(b1) == 7

    # hand the rays into BMO again and trace them to a second, identical detector
    beams2 = [BMO.Beam(b1, i) for i in 1:length(b1)]
    @test_throws BoundsError BMO.Beam(b1, 8)
    det2 = BMO.Detector(0.05)
    BMO.zrotate3d!(det2, deg2rad(10))
    BMO.translate3d!(det2, [0, zdet + 0.05, 0])
    trace!(det2, beams2...)
    # BMO orders hits by tracing order, which is the beam order here
    b2 = RayBundle(det2)
    n = OpticsBase.normal(OpticsBase.port(b1))
    @test length(b2) == 7
    # b1 positions propagated along their directions to det2's plane equal b2's positions
    plane = OpticsBase.port(b2)
    for i in 1:7
        d = b1.direction[i]
        t = dot(OpticsBase.origin(plane) - b1.position[i], OpticsBase.normal(plane)) /
            dot(d, OpticsBase.normal(plane))
        @test norm(b1.position[i] + t * d - b2.position[i]) <= 1e-12
        @test norm(d - b2.direction[i]) <= 1e-12
    end
    # the reconstructed BMO rays start exactly at the bundle data
    for i in 1:7
        r = BMO.rays(beams2[i])[1]
        @test norm(BMO.position(r) - b1.position[i]) <= 1e-15
        @test norm(BMO.direction(r) - b1.direction[i]) <= 1e-15
        @test BMO.wavelength(r) == λ
    end
end
