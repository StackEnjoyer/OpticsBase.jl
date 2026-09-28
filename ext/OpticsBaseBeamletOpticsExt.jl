# Glue between BeamletOptics.jl (BMO) and OpticsBase: reads the hits of a BMO `Detector`
# into a `RayBundle` at a port in the detector plane.
#
# Targets BMO 0.13.10 and uses BMO internals (hit records, parabasal ray data, `hits`,
# `hit_point`); the compat bound and the integration tests detect breaking BMO releases.
# BMO itself is not changed.

module OpticsBaseBeamletOpticsExt

using OpticsBase: OpticsBase, PlanarPort, RayBundle, normal, origin, refractive_index,
                  power_normalization, ray_basis
using LinearAlgebra: cross, dot, norm
using StaticArrays: SMatrix, SVector, @SMatrix
import BeamletOptics as BMO

"""
    RayBundle(detector::BeamletOptics.Detector; port = nothing, power = nothing)

Reads the hits of a BeamletOptics (BMO) `Detector` into a [`RayBundle`](@ref) at a port in
the detector plane: one ray per hit, all mutually coherent. Call it after `solve_system!`
and before moving or emptying the detector. Requires `using BeamletOptics` (package
extension).

# Arguments

  - `detector`: a BMO `Detector` with at least one hit. All hits must share the vacuum
    wavelength and the refractive index of the medium in front of the detector.
  - `port` (keyword): the port of the bundle. Default: a [`PlanarPort`](@ref) at the
    detector position with normal `n = −orientation(detector)[:, 2]` (BMO's detector
    normal points against the beam; `n` points downstream), `u = −orientation(detector)[:, 1]`
    (BMO's local detector x axis), so `v = n × u` is BMO's local `−z`, and the refractive
    index of the hits. A given `port` must be a `PlanarPort` in the detector plane with the
    same normal (within `√eps`) and the refractive index of the hits.
  - `power` (keyword): total power in \\[W\\] of pure-ray bundles, split equally between
    the rays (BMO rays carry no power). Required for `RayHit` and `PolarizedRayHit`
    detectors; must be `nothing` for beamlet detectors, whose power comes from BMO.

# Conversion per hit type

  - `AstigmaticGaussianBeamletHit` → vectorial beamlet (`N = 3`). Position: chief-ray hit
    point; direction: chief-ray direction; `opl`: optical path length of the chief ray from
    the BMO source up to the detector (parent beams included). The beamlet matrix is
    `Q = U·H⁻¹` (symmetrized) from BMO's complex parabasal ray heights `h₁, h₂` and slopes
    `u₁, u₂` at the hit plane, projected without conjugation onto the
    [`ray_basis`](@ref) `(x̂, ŷ)` of the port: `H = [x̂·h₁ x̂·h₂; ŷ·h₁ ŷ·h₂]`, `U` likewise.
    The phasor is the chief-ray polarization `E/‖E‖` in the global frame times the unit
    phase of `g = √(area_ref/area)` (Gouy phase, with BMO's complex beam areas in
    \\[m²\\]). The power is `κ·A²·π/(k·√det(Im Q))` with `A = ‖E‖·|g|` in \\[V/m\\],
    `κ =` [`power_normalization`](@ref)`(port)` and `k = 2πn/λ`, evaluated with the exact
    Lagrange invariant `√det(Im Q) = λ/(nπ·|area|)`, i.e. `κ·‖E‖²·π·|area_ref|/2`. In free
    space this equals BMO's `optical_power(beamlet)`; the `Q` from BMO's parabasal rays
    fulfils the invariant only up to `O(θ²)` (`θ` the divergence angle, ~10⁻⁷ relative
    for mm beams).
  - `GaussianBeamletHit` → scalar beamlet (`N = 1`) with `Q = (1/q)·I`,
    `1/q = R + i·λ/(π n w²)` from BMO's `gauss_parameters` at the hit (`w` beam radius
    in \\[m\\], `R` wavefront curvature in \\[1/m\\]). Phasor `exp(i(arg E₀ + ψ))` with BMO's
    Gouy phase `ψ`; power: BMO's `optical_power(beamlet)`; `opl` as above.
  - `PolarizedRayHit` → vectorial pure rays (`N = 3`), phasor `E₀/‖E₀‖` (BMO's global
    polarization vector, including Fresnel phases); `RayHit` → scalar pure rays (`N = 1`),
    phasor 1. `opl` is BMO's optical path length of the hit.

# Conventions

  - All positions and directions are in BMO's global frame, which is the OpticsBase global
    frame, in \\[m\\]. Hit points are projected onto the port plane (they lie on it up to
    rounding). BMO's `opl` is absolute from its source, as OpticsBase requires.
  - BMO also uses exp(−iωt) and peak amplitudes with `Z₀ = 376.730313668 Ω`. The bundle
    reproduces BMO's beamlet field without BMO's `√|cos θ|` detector projection factor
    (OpticsBase fields carry no obliquity weighting).
  - Beamlets with `M² ≠ 1` are not representable by a single `Q`; the parabasal rays then
    give an effective `Q`.

Throws an `ArgumentError` if the detector has no hits, if the hits do not share wavelength
and refractive index, if a hit propagates against the port normal (`d·n ≤ 0`), if `port`
does not fit the detector, or if `power` is missing for pure rays or given for beamlets.
"""
function OpticsBase.RayBundle(detector::BMO.Detector; port = nothing, power = nothing)
    hs = BMO.hits(detector)
    (hs === nothing || isempty(hs)) &&
        throw(ArgumentError("RayBundle: the detector has no hits; run `solve_system!` first"))
    λ, n_ray = _common_wavelength_index(hs)
    p = _detector_port(detector, port, λ, n_ray)
    return _raybundle(p, λ, hs, power)
end

# ---------------------------------------------------------------------------------------
# Hit accessors

_chief_ray(h::BMO.AbstractRayHit) = h.ray
_chief_ray(h::BMO.GaussianBeamletHit) = BMO.rays(h.gauss.chief)[h.id]
_chief_ray(h::BMO.AstigmaticGaussianBeamletHit) = BMO.rays(h.agb.c)[h.id]

_hit_point(h) = SVector{3}(BMO.hit_point(h))
_hit_direction(h) = SVector{3}(BMO.direction(_chief_ray(h)))

function _common_wavelength_index(hs::AbstractVector)
    ray = _chief_ray(first(hs))
    λ = BMO.wavelength(ray)
    n_ray = BMO.refractive_index(ray)
    tol = sqrt(eps(float(typeof(λ))))
    for (i, h) in enumerate(hs)
        r = _chief_ray(h)
        isapprox(BMO.wavelength(r), λ; rtol = tol) ||
            throw(ArgumentError("RayBundle: all detector hits must share the wavelength, got $(BMO.wavelength(r)) m for hit $i and $λ m for hit 1"))
        isapprox(BMO.refractive_index(r), n_ray; rtol = tol) ||
            throw(ArgumentError("RayBundle: all detector hits must share the refractive index, got $(BMO.refractive_index(r)) for hit $i and $n_ray for hit 1"))
    end
    return λ, n_ray
end

# Optical path length of a BMO beam from its source up to (and including) ray `id`.
function _opl_upto(beam, id::Integer)
    parent = beam.parent
    opl = parent === nothing ? zero(BMO.wavelength(first(BMO.rays(beam)))) :
          BMO.optical_path_length(parent)
    for j in 1:id
        opl += BMO.optical_path_length(BMO.rays(beam)[j])
    end
    return opl
end

# ---------------------------------------------------------------------------------------
# Port

function _detector_port(detector::BMO.Detector, ::Nothing, λ, n_ray)
    R = BMO.orientation(detector)
    n = -SVector{3}(R[:, 2])
    u = -SVector{3}(R[:, 1])
    return PlanarPort(SVector{3}(BMO.position(detector)), n, u; refractive_index = n_ray)
end

function _detector_port(detector::BMO.Detector, port::PlanarPort, λ, n_ray)
    ref = _detector_port(detector, nothing, λ, n_ray)
    tol = sqrt(eps(float(typeof(λ))))
    n = normal(ref)
    norm(normal(port) - n) <= tol ||
        throw(ArgumentError("RayBundle: the given port must have the detector normal $(Vector(n)) (pointing downstream), got $(Vector(normal(port)))"))
    Δ = origin(port) - origin(ref)
    abs(dot(Δ, n)) <= tol * max(norm(Δ), λ) ||
        throw(ArgumentError("RayBundle: the origin of the given port must lie in the detector plane, got distance $(dot(Δ, n)) m"))
    isapprox(refractive_index(port), n_ray; rtol = tol) ||
        throw(ArgumentError("RayBundle: the refractive index of the given port must equal that of the detector hits ($n_ray), got $(refractive_index(port))"))
    return port
end

function _detector_port(::BMO.Detector, port, λ, n_ray)
    throw(ArgumentError("RayBundle: the port of a detector bundle must be a PlanarPort, got $(typeof(port))"))
end

# Hit positions (projected onto the port plane) and directions, with the check d·n > 0.
function _geometry(port::PlanarPort, hs::AbstractVector)
    n = normal(port)
    r0 = origin(port)
    positions = map(hs) do h
        p = _hit_point(h)
        p - dot(p - r0, n) * n
    end
    directions = map(_hit_direction, hs)
    for (i, d) in enumerate(directions)
        dot(d, n) > 0 ||
            throw(ArgumentError("RayBundle: detector hit $i propagates against the port normal (direction·normal = $(dot(d, n)) ≤ 0)"))
    end
    return positions, directions
end

# ---------------------------------------------------------------------------------------
# Pure rays

function _pure_power(::Nothing, M)
    throw(ArgumentError("RayBundle: BeamletOptics rays carry no power; pass the total power of the bundle in W via the `power` keyword"))
end
function _pure_power(power::Real, M)
    power >= 0 ||
        throw(ArgumentError("RayBundle: power must be non-negative, got $power"))
    return fill(float(power) / M, M)
end

function _raybundle(port, λ, hs::AbstractVector{<:BMO.RayHit}, power)
    positions, directions = _geometry(port, hs)
    opl = [BMO.optical_path_length(h) for h in hs]
    phasor = ones(Complex{typeof(λ)}, length(hs))
    return RayBundle(port, λ, positions, directions, opl, _pure_power(power, length(hs)),
        phasor)
end

function _raybundle(port, λ, hs::AbstractVector{<:BMO.PolarizedRayHit}, power)
    positions, directions = _geometry(port, hs)
    opl = [BMO.optical_path_length(h) for h in hs]
    phasor = map(hs) do h
        E = SVector{3}(BMO.polarization(h))
        E / norm(E)
    end
    return RayBundle(port, λ, positions, directions, opl, _pure_power(power, length(hs)),
        phasor)
end

# ---------------------------------------------------------------------------------------
# Beamlets

function _raybundle(port, λ, hs::AbstractVector{<:BMO.AbstractBeamletHit}, power)
    power === nothing ||
        throw(ArgumentError("RayBundle: the power of beamlets is taken from BeamletOptics; do not pass the `power` keyword for beamlet detectors"))
    positions, directions = _geometry(port, hs)
    data = map(h -> _beamlet(port, λ, h), hs)
    return RayBundle(port, λ, positions, directions, [b.opl for b in data],
        [b.power for b in data], [b.phasor for b in data];
        beamlet = [b.Q for b in data])
end

function _beamlet(port, λ, h::BMO.GaussianBeamletHit)
    chief = _chief_ray(h)
    n_ray = BMO.refractive_index(chief)
    z = h.l0 + length(chief)
    w, R, ψ, _ = BMO.gauss_parameters(h.gauss, z; hint = (BMO.hit_point(h), h.id))
    invq = R + im * λ / (π * n_ray * w^2)
    Q = SMatrix{2, 2}(invq, zero(invq), zero(invq), invq)
    phasor = cis(angle(BMO.electric_field(h.gauss)) + ψ)
    return (; Q, phasor, power = BMO.optical_power(h.gauss),
        opl = _opl_upto(h.gauss.chief, h.id))
end

function _beamlet(port, λ, h::BMO.AstigmaticGaussianBeamletHit)
    chief = _chief_ray(h)
    d = _hit_direction(h)
    l = length(chief)
    # Parabasal ray heights at the hit plane (slopes are constant along the segment)
    u1 = SVector{3}(h.u1)
    u2 = SVector{3}(h.u2)
    h1 = SVector{3}(h.h1) + l * u1
    h2 = SVector{3}(h.h2) + l * u2
    # Non-conjugating projections onto the (real) ray basis
    x̂, ŷ = ray_basis(port, d)
    H = @SMatrix [sum(x̂ .* h1) sum(x̂ .* h2); sum(ŷ .* h1) sum(ŷ .* h2)]
    U = @SMatrix [sum(x̂ .* u1) sum(x̂ .* u2); sum(ŷ .* u1) sum(ŷ .* u2)]
    Q = U / H
    Q = (Q + transpose(Q)) / 2
    # Complex beam area (h₁ × h₂)·d and amplitude factor √(area_ref/area), as in BMO
    area = sum(cross(h1, h2) .* d)
    g = sqrt(h.area_ref / area)
    E = SVector{3}(BMO.polarization(chief))
    normE = norm(E)
    # Power κ·A²·π/(k·√det(Im Q)) with A = ‖E‖·|g|, evaluated with the exact Lagrange
    # invariant √det(Im Q) = λ/(n π |area|) instead of BMO's parabasal rays, which fulfil
    # it only to O(θ²): κ·‖E‖²·π·|area_ref|/2, BMO's `optical_power` for free space.
    power = power_normalization(port) * normE^2 * π * abs(h.area_ref) / 2
    phasor = (E / normE) * (g / abs(g))
    return (; Q, phasor, power, opl = _opl_upto(h.agb.c, h.id))
end

end # module
