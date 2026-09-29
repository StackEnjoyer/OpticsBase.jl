# Integration tests of OpticsBaseOpticSimExt: OpticSim rays → RayBundle → OpticSim rays,
# with the mm/µm ↔ m unit conversion.
# Evaluated in its own module with `using OpticsBase, Test` (see runtests.jl); the runner
# does not load OpticSim, so this file does.

import OpticSim
using LinearAlgebra
using StaticArrays

@test Base.get_extension(OpticsBase, :OpticsBaseOpticSimExt) !== nothing

const λ_um = 0.6328           # OpticSim wavelength unit: µm
const λ = 0.6328e-6           # OpticsBase: m

# A tilted port (normal not along a coordinate axis), origin in m
const n_port = normalize(SVector(0.2, 1.0, -0.3))
const u_port = normalize(cross(SVector(0.0, 0.0, 1.0), n_port))
const r0 = SVector(0.001, 0.02, -0.002)
const port0 = PlanarPort(r0, n_port, u_port)

# Rays in OpticSim units (mm) with origins on the port plane and directions into it
function rays_on_port(; M = 7)
    v = cross(n_port, u_port)
    map(1:M) do j
        s = (j - (M + 1) / 2)                       # mm, transverse offset
        o_mm = 1000 * r0 + s * u_port + 0.5 * s^2 * v
        d = normalize(n_port + 0.02 * s * u_port - 0.01 * v)
        OpticSim.OpticalRay(o_mm, d, 1.0, λ_um)
    end
end

@testset "OpticSim → RayBundle: units and frame" begin
    rays = rays_on_port()
    b = RayBundle(rays, port0)
    @test b isa RayBundle{Float64}
    @test length(b) == length(rays)
    @test OpticsBase.wavelength(b) ≈ λ rtol = 1e-15
    @test OpticsBase.port(b) === port0
    for (r, p, d) in zip(rays, OpticsBase.positions(b), OpticsBase.directions(b))
        @test p ≈ 1e-3 * OpticSim.origin(r) rtol = 1e-14       # mm → m
        @test d ≈ OpticSim.direction(r) rtol = 1e-15
        @test abs(dot(p - r0, n_port)) <= 1e-15
    end
    # explicit wavelength in m, consistent with the rays
    @test RayBundle(rays, port0, λ).position == b.position
    @test_throws ArgumentError RayBundle(rays, port0, 1.064e-6)
    # geometric rays need the wavelength
    geo = [OpticSim.ray(r) for r in rays]
    bg = RayBundle(geo, port0, λ)
    @test bg.position == b.position
    @test bg.direction == b.direction
end

@testset "OpticSim → RayBundle → OpticSim round trip" begin
    rays = rays_on_port()
    b = RayBundle(rays, port0)
    back = [OpticSim.OpticalRay(b, i) for i in 1:length(b)]
    @test back isa Vector{<:OpticSim.OpticalRay{Float64, 3}}
    for (r, s) in zip(rays, back)
        @test norm(OpticSim.origin(s) - OpticSim.origin(r)) <= 1e-12 * norm(OpticSim.origin(r))
        @test norm(OpticSim.direction(s) - OpticSim.direction(r)) <= 1e-12
        @test OpticSim.wavelength(s) ≈ λ_um rtol = 1e-15            # m → µm
        @test OpticSim.power(s) == 1
        @test OpticSim.pathlength(s) == 0
    end
    @test OpticSim.power(OpticSim.OpticalRay(b, 1; power = 0.5)) == 0.5
    @test_throws BoundsError OpticSim.OpticalRay(b, 0)
    @test_throws BoundsError OpticSim.OpticalRay(b, length(b) + 1)
    # and RayBundle → OpticSim → RayBundle
    b2 = RayBundle(back, port0)
    @test maximum(norm.(b2.position - b.position)) <= 1e-12 * maximum(norm.(b.position))
    @test maximum(norm.(b2.direction - b.direction)) <= 1e-12
end

@testset "Rays off the port plane are carried forward" begin
    # origins 5 mm upstream of the port plane: intersected with the plane
    rays = rays_on_port()
    shifted = [OpticSim.OpticalRay(OpticSim.origin(r) - 5 * OpticSim.direction(r),
                   OpticSim.direction(r), 1.0, λ_um) for r in rays]
    b = RayBundle(shifted, port0)
    for (r, p) in zip(rays, b.position)
        @test norm(p - 1e-3 * OpticSim.origin(r)) <= 1e-15
    end
    # origin beyond the port plane: not propagated backwards
    beyond = [OpticSim.OpticalRay(OpticSim.origin(r) + 5 * OpticSim.direction(r),
                  OpticSim.direction(r), 1.0, λ_um) for r in rays]
    @test_throws ArgumentError RayBundle(beyond, port0)
    # direction against the port normal
    against = [OpticSim.OpticalRay(OpticSim.origin(r), -OpticSim.direction(r), 1.0, λ_um)
               for r in rays]
    @test_throws ArgumentError RayBundle(against, port0)
    # mixed wavelengths
    mixed = [rays[1], OpticSim.OpticalRay(OpticSim.origin(rays[2]),
        OpticSim.direction(rays[2]), 1.0, 1.064)]
    @test_throws ArgumentError RayBundle(mixed, port0)
    @test_throws ArgumentError RayBundle(OpticSim.OpticalRay{Float64, 3}[], port0)
end

@testset "Traced OpticSim rays (LensTrace) at a detector plane" begin
    # Biconvex lens, light along +y (OpticSim builds lenses along −z; rotate −z → +y)
    glass = OpticSim.AGFFileReader.CARGILLE.OG0608
    Rot = SMatrix{3, 3}(1.0, 0, 0, 0, 0, 1, 0, -1, 0)
    lens = OpticSim.SphericalLens(glass, 0.0, 50.0, -50.0, 5.0, 12.7)(
        OpticSim.Transform(Rot, SVector(0.0, 20.0, 0.0)))
    det = OpticSim.Rectangle(30.0, 30.0, SVector(0.0, -1.0, 0.0), SVector(0.0, 60.0, 0.0);
        interface = OpticSim.opaqueinterface(Float64))
    sys = OpticSim.CSGOpticalSystem(OpticSim.LensAssembly(lens), det)
    src = [OpticSim.OpticalRay(SVector(x, 0.0, 0.5x), SVector(0.0, 1.0, 0.0), 1.0, λ_um)
           for x in -3.0:1.0:3.0]
    traces = [OpticSim.trace(sys, r; test = true) for r in src]
    @test all(!isnothing, traces)
    traces = identity.(traces)
    portB = PlanarPort([0.0, 0.06, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0])
    b = RayBundle(traces, portB)
    @test OpticsBase.wavelength(b) ≈ λ rtol = 1e-15
    for (t, p, d) in zip(traces, b.position, b.direction)
        @test norm(p - 1e-3 * OpticSim.point(t)) <= 1e-15      # detector hit, mm → m
        @test d == OpticSim.direction(OpticSim.ray(t))
    end
end
