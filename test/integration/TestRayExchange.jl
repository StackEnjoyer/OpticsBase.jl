# Ray exchange between two ray tracers through OpticsBase: BeamletOptics (BMO, SI units,
# light along +y) and OpticSim (mm, µm). A fan of collimated rays is handed over at plane A
# as a `RayBundle`, traced by the other tracer through a biconvex spherical singlet to plane
# B, and compared with the first tracer doing the whole trace itself. Both directions.
# Evaluated in its own module with `using OpticsBase, Test` (see runtests.jl); the runner
# loads BeamletOptics but not OpticSim, so this file does.
#
# Tolerances: positions at B within 1e-9 m (measured ≈ 9e-13 m). Directions within 1e-10
# (measured ≈ 6e-11): BMO lenses are signed distance functions, found by ray marching to
# ~1e-10 m, and their surface normals limit BMO's refracted directions to ~1e-11–1e-10.
# OpticSim's analytic spheres agree with an exact analytic trace to ~1e-16, so the hand-over
# itself adds no measurable error.

import BeamletOptics as BMO
import OpticSim
using LinearAlgebra
using StaticArrays

const λ = 0.6328e-6                  # vacuum wavelength, m (OpticSim: µm, see below)
const λ_um = λ * 1e6

# Geometry in m (global frame, light along +y)
const yA = 0.010                     # hand-over plane A
const yL = 0.020                     # front vertex of the lens
const R1 = 0.050                     # front radius (centre downstream)
const R2 = -0.050                    # back radius (centre upstream): biconvex
const tL = 0.005                     # centre thickness
const dL = 0.0254                    # diameter
const yB = 0.060                     # plane B, before the focus (f ≈ 54 mm)

# Fan of 7 collimated rays at y = 0, offsets along a skew line in the xz plane (m)
const fan = [SVector(s * cosd(30), 0.0, s * sind(30)) for s in (-3:3) .* 1e-3]
const dir0 = SVector(0.0, 1.0, 0.0)

# Lens material: an OpticSim catalog glass; its index at λ is passed to BMO as a constant
const glass = OpticSim.AGFFileReader.CARGILLE.OG0608

# ---------------------------------------------------------------------------------------
# BMO side

function bmo_detector(y)
    det = BMO.Detector(0.05)
    BMO.translate3d!(det, [0, y, 0])
    return det
end

function bmo_lens(n)
    lens = BMO.SphericalLens(R1, R2, tL, dL, n)
    BMO.translate3d!(lens, [0, yL, 0])
    return lens
end

# Traces `beams` through `elements` and returns the detector (last element)
function bmo_trace(elements, beams)
    system = BMO.System(elements)
    for b in beams
        BMO.solve_system!(system, b)
    end
    return last(elements)
end

# ---------------------------------------------------------------------------------------
# OpticSim side (mm)

# OpticSim builds lenses along its −z axis; rotate local −z onto global +y (exact matrix)
const Rot = SMatrix{3, 3}(1.0, 0, 0, 0, 0, 1, 0, -1, 0)

function opticsim_system(y_det)
    lens = OpticSim.SphericalLens(glass, 0.0, 1e3R1, 1e3R2, 1e3tL, 1e3dL / 2)(
        OpticSim.Transform(Rot, SVector(0.0, 1e3yL, 0.0)))
    det = OpticSim.Rectangle(30.0, 30.0, SVector(0.0, -1.0, 0.0),
        SVector(0.0, 1e3y_det, 0.0); interface = OpticSim.opaqueinterface(Float64))
    return OpticSim.CSGOpticalSystem(OpticSim.LensAssembly(lens), det)
end

# `test = true`: no Monte Carlo Fresnel reflections, every ray is refracted
function opticsim_trace(sys, rays)
    traces = [OpticSim.trace(sys, r; test = true) for r in rays]
    @test all(!isnothing, traces)
    return identity.(traces)
end

const sysB = opticsim_system(yB)
const n_glass = OpticSim.AGFFileReader.index(glass, λ_um;
    temperature = OpticSim.temperature(sysB), pressure = OpticSim.pressure(sysB))

# ---------------------------------------------------------------------------------------

maxdiff(a, b) = maximum(norm.(a .- b))

@testset "Refractive index" begin
    @test 1.4 < n_glass < 1.5
end

@testset "BMO → RayBundle → OpticSim" begin
    # BMO traces the fan to plane A; hand over as RayBundle
    detA = bmo_trace([bmo_detector(yA)], [BMO.Beam(p, dir0, λ) for p in fan])
    bA = RayBundle(detA)
    @test length(bA) == length(fan)
    # OpticSim traces from A through the lens to B
    raysA = [OpticSim.OpticalRay(bA, i) for i in 1:length(bA)]
    # Reference: BMO traces the whole path through the same lens to B
    detB = bmo_trace([bmo_lens(n_glass), bmo_detector(yB)],
        [BMO.Beam(p, dir0, λ) for p in fan])
    bB_ref = RayBundle(detB)
    @test length(bB_ref) == length(fan)
    bB = RayBundle(opticsim_trace(sysB, raysA), OpticsBase.port(bB_ref))
    Δp = maxdiff(bB.position, bB_ref.position)
    Δd = maxdiff(bB.direction, bB_ref.direction)
    @info "BMO → OpticSim at plane B" max_position_error_m = Δp max_direction_error = Δd
    @test Δp <= 1e-9
    @test Δd <= 1e-10              # BMO SDF surface normals, see top of file
end

@testset "OpticSim → RayBundle → BMO" begin
    # OpticSim rays of the fan (mm, µm) handed over at plane A
    portA = PlanarPort([0.0, yA, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0])
    rays0 = [OpticSim.OpticalRay(1e3 * p, dir0, 1.0, λ_um) for p in fan]
    bA = RayBundle(rays0, portA)
    @test length(bA) == length(fan)
    # BMO traces from A through the lens to B
    detB = bmo_trace([bmo_lens(n_glass), bmo_detector(yB)],
        [BMO.Beam(bA, i) for i in 1:length(bA)])
    bB = RayBundle(detB)
    # Reference: OpticSim traces the whole path through the same lens to B
    bB_os = RayBundle(opticsim_trace(sysB, rays0), OpticsBase.port(bB))
    Δp = maxdiff(bB.position, bB_os.position)
    Δd = maxdiff(bB.direction, bB_os.direction)
    @info "OpticSim → BMO at plane B" max_position_error_m = Δp max_direction_error = Δd
    @test Δp <= 1e-9
    @test Δd <= 1e-10              # BMO SDF surface normals, see top of file
end
