# Second chain: BeamletOptics → OpticsBase → vectorial high-NA focus
#
# A collimated, linearly polarized beam is traced with BeamletOptics.jl as polarized rays
# over a 90° off-axis parabolic mirror (OAP) to a detector in front of the focus. The hits
# are handed over as a `PolarizedRayBundle`, converted into plane waves by `DebyeWolf` (Debye
# approximation, solid angles from Voronoi cells) and summed on a grid in the focal plane.
# The result is compared with the same computation from analytic mirror rays with exact
# solid angles.
#
# Run from the repository root with the integration environment:
#   julia --project=test/integration -e 'using Pkg; Pkg.instantiate()'
#   julia --project=test/integration examples/bmo_oap_focus.jl

using OpticsBase
using BeamletOptics          # activates OpticsBaseBeamletOpticsExt
using DelaunayTriangulation  # activates OpticsBaseDelaunayTriangulationExt
using LinearAlgebra
using StaticArrays
using Printf

const OB = OpticsBase
const BMO = BeamletOptics

# --- Setup (SI units) ---------------------------------------------------------------------
λ = 1.064e-6                 # vacuum wavelength [m]
rfl = 25e-3                  # reflected focal length of the OAP [m]
D = 20e-3                    # diameter of the collimated input beam [m]
f_p = rfl / 2                # parent focal length of a 90° OAP [m]
P0 = 1e-3                    # power [W]
E_in = [0.0, 0, 1]           # input polarization: perpendicular to the deflection plane
y0 = f_p - 30e-3             # start plane of the rays [m]
x_det = 0.4 * rfl            # detector plane between mirror and focus [m]
F = SVector(0.0, f_p, 0)     # focus: the OAP deflects +y into −x towards F
M = 20_000                   # number of rays

# --- 1. BeamletOptics: polarized rays over the OAP to a detector ---------------------------
oap = BMO.OffAxisParabolicMirror(rfl, D; angle = 90)
BMO.translate_to3d!(oap, [rfl, f_p, 0])
detector = BMO.Detector(25e-3)
BMO.zrotate3d!(detector, deg2rad(90))
BMO.translate_to3d!(detector, [x_det, f_p, 0])
beams = map(0:(M - 1)) do j                       # equal-area (Fibonacci) input sampling
    r = D / 2 * sqrt((j + 0.5) / M)
    ψ = j * π * (3 - sqrt(5))
    BMO.Beam([rfl + r * cos(ψ), y0, r * sin(ψ)], [0.0, 1, 0], λ, E_in)
end
BMO.solve_system!(BMO.StaticSystem([oap, detector]),
    BMO.CollimatedSource(beams, D, [rfl, y0, 0.0], [0.0, 1, 0]))

# --- 2. Handover: detector hits → PolarizedRayBundle (BMO rays carry no power) ---
bundle = PolarizedRayBundle(detector; power = P0)

# --- 3. Debye: rays → plane waves at a port in the focus ----------------------------------
port_F = PlanarPort(F, [-1.0, 0, 0], [0.0, 0, 1])   # n along the central ray, u along E_in
spectrum = convert_field(DebyeWolf(port_F), bundle)

# --- 4. Focal field: plane-wave summation on ±3λ -------------------------------------------
grid = RegularGrid((64, 64), (3λ / 32, 3λ / 32))
E = OB.field_array(convert_field(PlaneWaveSummation(grid), spectrum))

# --- 5. Reference: analytic mirror rays, Gauss–Legendre pupil, exact solid angles ----------
function gauss_legendre(m)
    β = [i / sqrt(4i^2 - 1) for i in 1:(m - 1)]
    G = eigen(SymTridiagonal(zeros(m), β))
    return G.values, 2 .* G.vectors[1, :] .^ 2
end

function reflect_polarization(e, d_in, d_out)       # BMO's ideal mirror, diag(−1, 1) in (s, p)
    s = normalize(cross(d_in, d_out))
    return -dot(e, s) * s + dot(e, cross(d_in, s)) * cross(d_out, s)
end

function reference_spectrum(port_A; nr = 64, nφ = 128)
    x, wx = gauss_legendre(nr)
    C = f_p + rfl                                   # the paraboloid obeys ‖q − F‖ = C − q_y
    rays = map(vec(collect(Iterators.product(zip(x, wx), 1:nφ)))) do ((xi, wi), l)
        r, φ = D / 4 * (xi + 1), 2π * (l - 0.5) / nφ
        dA = D / 4 * wi * r * 2π / nφ               # input area of the ray [m²]
        xin, zin = rfl + r * cos(φ), r * sin(φ)
        t = (C^2 - f_p^2 - xin^2 - zin^2) / (2 * (C - f_p))
        q = SVector(xin, t, zin)                    # hit point on the mirror
        d = normalize(F - q)
        p = q + ((x_det - q[1]) / d[1]) * d         # position on the detector plane
        (; p, d, opl = (t - y0) + norm(p - q), P = P0 * dA / (π * (D / 2)^2),
            e = reflect_polarization(SVector(E_in...), SVector(0.0, 1, 0), d),
            w = dA / norm(F - q)^2)                 # dΩ = dA/‖q − F‖² for a paraboloid
    end
    ref = PolarizedRayBundle(port_A, λ, [r.p for r in rays], [r.d for r in rays],
        [r.opl for r in rays], [r.P for r in rays], [r.e for r in rays])
    return convert_field(DebyeWolf(port_F; solid_angle = [r.w for r in rays]), ref)
end
E_ref = OB.field_array(convert_field(PlaneWaveSummation(grid),
    reference_spectrum(OB.port(bundle))))

# --- 6. Report -----------------------------------------------------------------------------
u, v, n = eachcol(OB.local_axes(port_F))
sin_u = maximum(d -> abs(dot(d, u)), bundle.direction)      # half-aperture sines per plane
sin_v = maximum(d -> abs(dot(d, v)), bundle.direction)
miss = maximum(zip(bundle.position, bundle.direction)) do (p, d)
    norm((F - p) - dot(F - p, d) * d)
end
to_focus = [o + norm(F - p) for (o, p) in zip(bundle.opl, bundle.position)]
κ = OB.power_normalization(port_F)
E_long = [dot(E[i, j, :], n) for i in axes(E, 1), j in axes(E, 2)]
long_fraction = sum(abs2, E_long) / sum(abs2, E)
deviation = maximum(abs, E - E_ref) / maximum(abs, E_ref)

@printf("NA along u, v                        %8.3f %8.3f\n", sin_u, sin_v)
@printf("rays, max. distance to focus [nm]    %8d %8.2e\n", length(bundle), 1e9 * miss)
@printf("OPL spread to focus [λ]              %17.1e\n", (maximum(to_focus) - minimum(to_focus)) / λ)
@printf("power [mW]                           %17.12f\n", 1e3 * total_power(spectrum))
@printf("peak intensity κ|E|² [mW/µm²]        %17.4f\n", 1e-9 * κ * maximum(sum(abs2, E; dims = 3)))
@printf("longitudinal energy fraction (±3λ)   %17.4f\n", long_fraction)
@printf("max |E − E_ref| / max |E_ref|        %17.1e\n", deviation)
