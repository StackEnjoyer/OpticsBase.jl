# First chain: BeamletOptics → OpticsBase → Fourier solver
#
# A Gaussian beamlet is traced with BeamletOptics.jl to a detector, its field is sampled on
# a grid (`SampledField`, using BeamletOptics' own beamlet sum) and propagated further with
# the angular-spectrum method of WaveOpticsPropagation.jl. The result is compared with the
# analytic Gaussian beam and with BeamletOptics' own field at the output plane.
#
# Run from the repository root with the integration environment:
#   julia --project=test/integration -e 'using Pkg; Pkg.instantiate()'
#   julia --project=test/integration examples/bmo_to_fourier.jl

using OpticsBase
using BeamletOptics          # activates OpticsBaseBeamletOpticsExt
using WaveOpticsPropagation  # activates OpticsBaseWaveOpticsPropagationExt
using Printf

const OB = OpticsBase
const BMO = BeamletOptics

# --- Beam and planes (SI units) ---------------------------------------------------------
λ = 1.064e-6                 # vacuum wavelength [m]
w0 = 0.5e-3                  # waist radius [m], waist at the origin
P0 = 1e-3                    # power [W]
zR = π * w0^2 / λ            # Rayleigh length [m]
yA = 0.1                     # handover plane A, 0.1 m behind the waist (BMO's axis is +y)
zAB = zR / 2                 # angular-spectrum propagation A → B [m]
w_analytic(y) = w0 * sqrt(1 + (y / zR)^2)

# --- 1. BeamletOptics: trace an astigmatic Gaussian beamlet to a detector at plane A ----
function traced_detector(y)
    beamlet = BMO.AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0,
        support = [1.0, 0, 0])
    detector = BMO.Detector(0.05)
    BMO.translate3d!(detector, [0, y, 0])
    BMO.solve_system!(BMO.System([detector]), beamlet)
    return detector
end
detector_A = traced_detector(yA)

# --- 2./3. Handover: the detector's beamlet field on a grid at a port in the detector plane -
grid = RegularGrid((256, 256), (w0 / 16, w0 / 16))
field_A = SampledField(detector_A, grid)  # vectorial (astigmatic beamlet), BMO's own sum
port_A = OB.port(field_A)                 # n = +y (downstream), u = +x

# --- 4. Fourier solver: free-space propagation to a parallel port B ----------------------
port_B = PlanarPort(OB.origin(port_A) + zAB * OB.normal(port_A), OB.normal(port_A),
    OB.local_axes(port_A)[:, 1])
sol = solve(PropagationProblem(field_A, FreeSpace(), port_B), AngularSpectrumMethod())
field_B = sol.field

# --- 5. Checks ----------------------------------------------------------------------------
# Beam radius from the second moment of |E|² (for a Gaussian ⟨ξ²⟩ = w²/4)
function beam_radius(f)
    I = dropdims(sum(abs2, OB.field_array(f); dims = 3); dims = 3)
    ξ = OB.coordinates(OB.grid(f), 1)
    return 2 * sqrt(sum(I .* ξ .^ 2) / sum(I))
end

# BeamletOptics' own field on plane B. Its detector frame (x, z) maps to the port frame as
# ξ = x, η = −z, so the η axis is reversed.
detector_B = traced_detector(yA + zAB)
ξ = OB.coordinates(grid, 1)
η = OB.coordinates(grid, 2)
_, _, E_bmo = BMO.electric_field(detector_B; n = 256, x_min = first(ξ), x_max = last(ξ),
    z_min = -last(η), z_max = -first(η), progress = false)
E_ob = OB.field_array(field_B)[:, end:-1:1, 1]
deviation = maximum(abs, E_ob - E_bmo) / maximum(abs, E_bmo)

@printf("%-28s %12s %12s\n", "", "plane A", "plane B")
@printf("%-28s %12.2f %12.2f\n", "distance from waist [mm]", 1e3 * yA, 1e3 * (yA + zAB))
@printf("%-28s %12.5f %12.5f\n", "beam radius, analytic [mm]", 1e3 * w_analytic(yA),
    1e3 * w_analytic(yA + zAB))
@printf("%-28s %12.5f %12.5f\n", "beam radius, OpticsBase [mm]", 1e3 * beam_radius(field_A),
    1e3 * beam_radius(field_B))
@printf("%-28s %12.9f %12.9f\n", "power [mW]", 1e3 * total_power(field_A),
    1e3 * total_power(field_B))
@printf("\nmax |E_OpticsBase − E_BeamletOptics| / max |E| at plane B: %.1e\n", deviation)
