# Integration tests of OpticsBaseBeamletOpticsExt: BMO detector hits → RayBundle.
# Evaluated in its own module with `using OpticsBase, Test` (see runtests.jl).

import BeamletOptics as BMO
using LinearAlgebra
using StaticArrays

const λ = 1.064e-6
const w0 = 0.5e-3
const P0 = 1e-3
const zdet = 0.2

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

# Analytic complex beam parameter q = z − i z_R (exp(−iωt)) of a waist w at distance z
q_analytic(z, w) = z - im * π * w^2 / λ

@testset "Astigmatic beamlet (stigmatic waist)" begin
    agb = BMO.AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0,
        support = [1.0, 0, 0])
    det = trace!(detector_at(zdet), agb)
    b = RayBundle(det)

    @test b isa RayBundle{3}
    @test length(b) == 1
    @test OpticsBase.has_beamlets(b)
    @test OpticsBase.wavelength(b) == λ
    @test total_power(b) ≈ BMO.optical_power(agb) rtol = 1e-9
    @test b.position[1] ≈ SVector(0.0, zdet, 0.0) atol = 1e-15
    @test b.direction[1] ≈ SVector(0.0, 1.0, 0.0)
    @test b.opl[1] ≈ zdet rtol = 1e-12

    invq = 1 / q_analytic(zdet, w0)
    Q = b.beamlet[1]
    @test Q[1, 1] ≈ invq rtol = 1e-6
    @test Q[2, 2] ≈ invq rtol = 1e-6
    @test abs(Q[1, 2]) <= 1e-6 * abs(invq)
    @test imag(Q[1, 1]) > 0

    # Polarization along the support axis (port u = +x); phase = Gouy phase
    zR = π * w0^2 / λ
    e = b.phasor[1]
    @test abs(e[1]) ≈ 1
    @test abs(e[2]) + abs(e[3]) <= 1e-12
    @test angle(e[1]) ≈ -atan(zdet / zR) atol = 1e-6
end

@testset "Astigmatic beamlet (w0x ≠ w0y)" begin
    w0x, w0y = 0.5e-3, 0.3e-3
    agb = BMO.AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0x, w0y; P0,
        support = [1.0, 0, 0])
    det = trace!(detector_at(zdet), agb)
    b = RayBundle(det)
    # Support axis s1 = +x = port u, s2 = d × s1 = −z = port v: Q is diagonal
    Q = b.beamlet[1]
    @test Q[1, 1] ≈ 1 / q_analytic(zdet, w0x) rtol = 1e-6
    @test Q[2, 2] ≈ 1 / q_analytic(zdet, w0y) rtol = 1e-6
    @test abs(Q[1, 2]) <= 1e-6 * abs(Q[1, 1])
    @test total_power(b) ≈ BMO.optical_power(agb) rtol = 1e-9
    # Gouy phase of an astigmatic beam: mean of the two axes
    ψ = -(atan(zdet * λ / (π * w0x^2)) + atan(zdet * λ / (π * w0y^2))) / 2
    @test angle(b.phasor[1][1]) ≈ ψ atol = 1e-6
end

@testset "Stigmatic GaussianBeamlet" begin
    gb = BMO.GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0, support = [1.0, 0, 0])
    det = trace!(detector_at(zdet), gb)
    b = RayBundle(det)
    @test b isa RayBundle{1}
    @test length(b) == 1
    @test total_power(b) ≈ BMO.optical_power(gb) rtol = 1e-12
    @test total_power(b) ≈ P0 rtol = 1e-9
    @test b.opl[1] ≈ zdet rtol = 1e-12
    invq = 1 / q_analytic(zdet, w0)
    Q = b.beamlet[1]
    @test Q[1, 1] ≈ invq rtol = 1e-6
    @test Q[2, 2] == Q[1, 1]
    @test Q[1, 2] == 0
    zR = π * w0^2 / λ
    @test angle(b.phasor[1][1]) ≈ -atan(zdet / zR) atol = 1e-6

    # Same Q as the astigmatic beamlet
    agb = BMO.AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0,
        support = [1.0, 0, 0])
    ba = RayBundle(trace!(detector_at(zdet), agb))
    @test ba.beamlet[1] ≈ Q rtol = 1e-6
end

@testset "Pure rays" begin
    rays_at(x) = BMO.Beam(BMO.Ray([x, 0.0, 0.0], [0.0, 1.0, 0.0], λ))
    det = trace!(detector_at(zdet), rays_at(0.0), rays_at(1e-3))
    b = RayBundle(det; power = 2e-3)
    @test b isa RayBundle{1}
    @test !OpticsBase.has_beamlets(b)
    @test length(b) == 2
    @test total_power(b) ≈ 2e-3
    @test b.power ≈ [1e-3, 1e-3]
    @test all(==(SVector(complex(1.0))), b.phasor)
    @test b.opl ≈ [zdet, zdet]
    @test b.position[2] ≈ SVector(1e-3, zdet, 0.0)
    @test_throws ArgumentError RayBundle(det)
    @test_throws ArgumentError RayBundle(det; power = -1.0)

    E0 = [0.0, 0.0, 2.0im]
    pol_at(x) = BMO.Beam(BMO.PolarizedRay([x, 0.0, 0.0], [0.0, 1.0, 0.0], λ, E0))
    det = trace!(detector_at(zdet), pol_at(0.0), pol_at(1e-3))
    b = RayBundle(det; power = 3e-3)
    @test b isa RayBundle{3}
    @test total_power(b) ≈ 3e-3
    @test all(e -> e ≈ SVector(0.0, 0.0, 1.0im), b.phasor)
    @test b.opl ≈ [zdet, zdet]
    @test_throws ArgumentError RayBundle(det)
end

@testset "Port from the detector" begin
    # Right-handed port with d·n > 0, for a straight and a tilted detector
    for tilt in (0.0, deg2rad(20), deg2rad(-35))
        agb = BMO.AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0,
            support = [1.0, 0, 0])
        det = trace!(detector_at(zdet; tilt), agb)
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
        # The beamlet power does not depend on the detector tilt
        @test total_power(b) ≈ BMO.optical_power(agb) rtol = 1e-9
    end

    # Pure rays on a tilted detector
    det = trace!(detector_at(zdet; tilt = deg2rad(30)),
        BMO.Beam(BMO.Ray([0.0, 0, 0], [0.0, 1.0, 0.0], λ)),
        BMO.Beam(BMO.Ray([1e-3, 0, 0], [0.0, 1.0, 0.0], λ)))
    b = RayBundle(det; power = 1e-3)
    n = OpticsBase.normal(OpticsBase.port(b))
    @test all(d -> dot(d, n) > 0, b.direction)
    @test all(p -> abs(dot(p - OpticsBase.origin(OpticsBase.port(b)), n)) < 1e-15,
        b.position)
end

@testset "Port keyword and errors" begin
    agb = BMO.AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0,
        support = [1.0, 0, 0])
    det = trace!(detector_at(zdet), agb)
    # A port in the detector plane with the same normal but another u and origin
    p = PlanarPort([1e-3, zdet, 2e-3], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])
    b = RayBundle(det; port = p)
    @test OpticsBase.port(b) === p
    @test b.beamlet[1] ≈ RayBundle(det).beamlet[1] rtol = 1e-6   # stigmatic: basis-free
    @test total_power(b) ≈ total_power(RayBundle(det)) rtol = 1e-12

    @test_throws ArgumentError RayBundle(det;
        port = PlanarPort([0.0, zdet, 0], [0.0, -1.0, 0.0], [1.0, 0, 0]))
    @test_throws ArgumentError RayBundle(det;
        port = PlanarPort([0.0, zdet + 1e-3, 0], [0.0, 1.0, 0.0], [1.0, 0, 0]))
    @test_throws ArgumentError RayBundle(det;
        port = PlanarPort([0.0, zdet, 0], [0.0, 1.0, 0.0], [1.0, 0, 0];
            refractive_index = 1.5))
    @test_throws ArgumentError RayBundle(det; port = :notaport)
    @test_throws ArgumentError RayBundle(det; power = 1e-3)

    @test_throws ArgumentError RayBundle(detector_at(zdet))  # no hits

    # Hits with different wavelengths
    det = trace!(detector_at(zdet),
        BMO.Beam(BMO.Ray([0.0, 0, 0], [0.0, 1.0, 0.0], λ)),
        BMO.Beam(BMO.Ray([1e-3, 0, 0], [0.0, 1.0, 0.0], 2λ)))
    @test_throws ArgumentError RayBundle(det; power = 1e-3)
end
