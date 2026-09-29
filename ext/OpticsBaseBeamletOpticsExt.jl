# BMO-side glue between BeamletOptics.jl (BMO) and OpticsBase. This is BMO's code, living
# here temporarily (it may move into BeamletOptics once OpticsBase is registered). It
# reads the hits of a BMO `Detector` into `RayBundle`, `PolarizedRayBundle` and
# `SampledField`, and builds BMO beams from a `RayBundle`.
#
# Targets BMO 0.13.10 and uses BMO internals (hit records, parabasal ray data, `hits`,
# `hit_point`); the compat bound and the integration tests detect breaking BMO releases.
# BMO itself is not changed. No solver numerics live in OpticsBase `src/`: the coherent sum
# of the beamlets is BMO's own `electric_field`.

module OpticsBaseBeamletOpticsExt

using OpticsBase: OpticsBase, PlanarPort, RayBundle, PolarizedRayBundle, SampledField,
                  AbstractGrid, normal, origin, refractive_index, local_axes, spacing
using LinearAlgebra: cross, dot, norm
using StaticArrays: SVector
import BeamletOptics as BMO

const _PREFIX = "BeamletOptics detector"

# ---------------------------------------------------------------------------------------
# Public methods

"""
    RayBundle(detector::BeamletOptics.Detector; port = nothing)

Reads the hits of a BeamletOptics (BMO) `Detector` into a geometry-only
[`RayBundle`](@ref OpticsBase.RayBundle) at a port in the detector plane: one ray per hit (for beamlet hits the
chief ray). Call it after `solve_system!` and before moving or emptying the detector.
Requires `using BeamletOptics` (package extension). This is BMO-side glue living in
OpticsBase temporarily; it will move into BeamletOptics.

# Arguments

  - `detector`: a BMO `Detector` with at least one hit (any hit type). All hits must share
    the vacuum wavelength and the refractive index of the medium in front of the detector.
  - `port` (keyword): the port of the bundle. Default: a [`PlanarPort`](@ref OpticsBase.PlanarPort) at the
    detector position with normal `n = −orientation(detector)[:, 2]` (BMO's detector
    normal points against the beam; `n` points downstream), `u = −orientation(detector)[:, 1]`
    (BMO's local detector x axis), so `v = n × u` is BMO's local `−z`, and the refractive
    index of the hits. A given `port` must be a `PlanarPort` in the detector plane with the
    same normal (within `√eps`) and the refractive index of the hits.

# Conventions

Positions in \\[m\\] and directions are in BMO's global frame, which is the OpticsBase global
frame. Hit points are projected onto the port plane (they lie on it up to rounding). Power,
optical path length and polarization are not part of a `RayBundle`; see
[`PolarizedRayBundle`](@ref OpticsBase.PolarizedRayBundle) for those.

Throws an `ArgumentError` if the detector has no hits, if the hits do not share wavelength
and refractive index, if a hit propagates against the port normal (`d·n ≤ 0`), or if `port`
does not fit the detector.
"""
function OpticsBase.RayBundle(detector::BMO.Detector; port = nothing)
    hs = _hits(detector)
    λ, n_ray = _common_wavelength_index(hs)
    p = _detector_port(detector, port, λ, n_ray)
    positions, directions = _geometry(p, hs)
    return RayBundle(p, λ, positions, directions)
end

"""
    PolarizedRayBundle(detector::BeamletOptics.Detector; port = nothing, power)

Reads the polarized ray hits (`PolarizedRayHit`) of a BeamletOptics (BMO) `Detector` into a
[`PolarizedRayBundle`](@ref OpticsBase.PolarizedRayBundle) at a port in the detector plane: one mutually coherent ray per
hit. Call it after `solve_system!`. Requires `using BeamletOptics`. This is BMO-side glue
living in OpticsBase temporarily; it will move into BeamletOptics.

# Arguments

  - `detector`: a BMO `Detector` whose hits are all `PolarizedRayHit`. Hits of pure rays
    (`RayHit`, no polarization) and of beamlets (whose field BMO sums itself, see
    `SampledField(detector, grid)`) throw an `ArgumentError`. All hits must share the vacuum
    wavelength and the refractive index of the medium in front of the detector.
  - `port` (keyword): as for `RayBundle(detector; port)`.
  - `power` (keyword, required): total power of the bundle in \\[W\\], split equally between
    the rays, since BMO rays carry no power.

# Conventions

  - `opl` is BMO's optical path length of the hit in \\[m\\], absolute from the BMO source.
  - `polarization` is BMO's global polarization vector of the hit `E₀/‖E₀‖` (complex, global
    frame, includes Fresnel phases); only its direction is used, the amplitude of `E₀` is
    replaced by `power`.
  - Positions and directions as for `RayBundle(detector)`.

Throws an `ArgumentError` under the conditions of `RayBundle(detector)`, for unsupported hit
types and for a negative `power`.
"""
function OpticsBase.PolarizedRayBundle(detector::BMO.Detector; port = nothing, power)
    hs = _hits(detector)
    λ, n_ray = _common_wavelength_index(hs)
    p = _detector_port(detector, port, λ, n_ray)
    return _polarized_bundle(p, λ, hs, power)
end

function _polarized_bundle(port, λ, hs::AbstractVector{<:BMO.PolarizedRayHit}, power)
    power isa Real ||
        throw(ArgumentError("$_PREFIX: power must be a real number in W, got $(typeof(power))"))
    power >= 0 ||
        throw(ArgumentError("$_PREFIX: power must be non-negative, got $power"))
    positions, directions = _geometry(port, hs)
    opl = [BMO.optical_path_length(h) for h in hs]
    pol = map(hs) do h
        E = SVector{3}(BMO.polarization(h))
        E / norm(E)
    end
    return PolarizedRayBundle(port, λ, positions, directions, opl,
        fill(float(power) / length(hs), length(hs)), pol)
end

function _polarized_bundle(port, λ, hs::AbstractVector, power)
    throw(ArgumentError("$_PREFIX: PolarizedRayBundle needs polarized ray hits (PolarizedRayHit), got $(eltype(hs))"))
end

"""
    SampledField(detector::BeamletOptics.Detector, grid::AbstractGrid{2}; port = nothing)

Samples the coherent field of the Gaussian beamlets hitting a BeamletOptics (BMO)
`Detector` on `grid` at a port in the detector plane, using BMO's own beamlet field
functions (`BeamletOptics.electric_field`). Call it after `solve_system!` and before moving
or emptying the detector. Requires `using BeamletOptics`. This is BMO-side glue living in
OpticsBase temporarily; it will move into BeamletOptics.

# Arguments

  - `detector`: a BMO `Detector` whose hits are all `GaussianBeamletHit` or all
    `AstigmaticGaussianBeamletHit`. Pure-ray hits (`RayHit`, `PolarizedRayHit`) have no
    scaled field in BMO and throw an `ArgumentError`; use `RayBundle` /
    `PolarizedRayBundle` for them. All hits must share wavelength and refractive index.
  - `grid`: a square `n × n` [`RegularGrid`](@ref OpticsBase.RegularGrid) with equal spacing `Δ` (BMO samples one
    `n` per axis). Sample `i` lies at `(i − (n÷2 + 1))·Δ` from the port origin along the
    port axes, as always.
  - `port` (keyword): default as for `RayBundle(detector)`. A given port must be a
    `PlanarPort` in the detector plane with the detector normal, refractive index of the
    hits and the same axes `(u, v)` as the default port (`u = ` BMO's local detector x,
    `v = ` BMO's local `−z`); its origin may be shifted within the plane.

# Result

  - Stigmatic `GaussianBeamlet` hits give a scalar `SampledField{1}`: BMO carries no
    polarization for them.
  - `AstigmaticGaussianBeamlet` hits give a vectorial `SampledField{3}`. BMO's detector
    method sums only the scalar part. Per hit, this method takes the scalar detector field
    of that single hit and multiplies it by the unit polarization vector of the hit's chief
    ray (`polarization(chief) / E_ref_amp`, global frame); the vector fields of all hits are
    then summed coherently.

# Conventions

  - The field is the physical peak amplitude in \\[V/m\\], global Cartesian components, with
    the full spatial phase and BMO's absolute optical path length; time dependence
    exp(−iωt), as in BMO. BMO uses `Z₀ = 376.730313668 Ω` and peak amplitudes like
    OpticsBase, so `total_power` reproduces BMO's beamlet power without rescaling.
  - BMO includes its detector projection factor `√|cos θ|` (θ the angle between beam and
    detector normal), so the field is the field on the tilted detector plane and its power
    `κ Σ|E|² ΔξΔη` is the power crossing that plane.
  - Sample mapping: `ξ = x`, `η = −z` in BMO's local detector coordinates; the second axis
    is reversed with respect to BMO's output.
  - Accuracy: the sum is exact for the beamlets as BMO defines them (paraxial
    parabasal-ray beamlets); the grid must resolve the field (Nyquist for the phase
    curvature is the user's responsibility).

Throws an `ArgumentError` under the conditions of `RayBundle(detector)`, for a non-square
grid or unequal spacing, for a `port` with other axes or normal than the detector, and for
pure-ray hits.
"""
function OpticsBase.SampledField(detector::BMO.Detector, grid::AbstractGrid{2};
        port = nothing)
    hs = _hits(detector)
    λ, n_ray = _common_wavelength_index(hs)
    p = _detector_port(detector, port, λ, n_ray)
    n, m = size(grid)
    Δ1, Δ2 = spacing(grid)
    (n == m && Δ1 ≈ Δ2) ||
        throw(ArgumentError("$_PREFIX: the grid must be square with equal spacing (BMO samples one n per axis), got size $(size(grid)) and spacing $(spacing(grid))"))
    ref = _detector_port(detector, nothing, λ, n_ray)
    if port !== nothing
        tol = sqrt(eps(float(typeof(λ))))
        norm(local_axes(p) - local_axes(ref)) <= tol ||
            throw(ArgumentError("$_PREFIX: the given port must have the axes of the detector port (u = BMO local x, v = BMO local −z)"))
    end
    # Offset of the port origin within the detector plane, in BMO's local (x, z)
    R = BMO.orientation(detector)
    Δr = origin(p) - origin(ref)
    x_shift = dot(Δr, -SVector{3}(R[:, 1]))
    z_shift = dot(Δr, SVector{3}(R[:, 3]))
    # Sample i at (i − (n÷2 + 1))Δ; BMO's z axis is the reversed η axis
    x_min = -(n ÷ 2) * Δ1
    x_max = (n - 1 - n ÷ 2) * Δ1
    E = _summed_field(detector, hs; n, x_min, x_max, z_min = -x_max, z_max = -x_min,
        x0_shift = x_shift, z0_shift = z_shift, progress = false)
    E = reverse(E; dims = 2)          # η = −z
    return SampledField(E, grid, p, λ)
end

# Stigmatic: BMO's scalar detector field, N = 1
function _summed_field(detector, hs::AbstractVector{<:BMO.GaussianBeamletHit}; kwargs...)
    return BMO.electric_field(detector, hs; kwargs...)[3]
end

# Astigmatic: BMO's scalar field per hit times the unit polarization of the hit, N = 3
function _summed_field(detector, hs::AbstractVector{<:BMO.AstigmaticGaussianBeamletHit};
        kwargs...)
    n = kwargs[:n]
    T = real(eltype(BMO.electric_field(detector, hs[1:1]; kwargs...)[3]))
    E = zeros(Complex{T}, n, n, 3)
    for h in hs
        s = BMO.electric_field(detector, [h]; kwargs...)[3]
        e = SVector{3}(BMO.polarization(_chief_ray(h))) / h.E_ref_amp
        for c in 1:3
            @views E[:, :, c] .+= s .* e[c]
        end
    end
    return E
end

function _summed_field(detector, hs::AbstractVector; kwargs...)
    throw(ArgumentError("$_PREFIX: SampledField needs Gaussian beamlet hits (GaussianBeamletHit or AstigmaticGaussianBeamletHit); pure-ray hits ($(eltype(hs))) have no scaled field in BMO"))
end

"""
    BeamletOptics.Beam(bundle::RayBundle, i::Integer)

Builds the BMO `Beam` of ray `i` of a [`RayBundle`](@ref OpticsBase.RayBundle), for handing rays from another
ray tracer (e.g. OpticSim.jl) into BMO. Requires `using BeamletOptics`. This is BMO-side
glue living in OpticsBase temporarily; it will move into BeamletOptics.

The beam has one BMO `Ray` starting at `position[i]` (\\[m\\], global frame) with direction
`direction[i]` and vacuum wavelength `wavelength(bundle)` (\\[m\\]) in vacuum (BMO rays start
with refractive index 1; the port index is not transferred). Use
`[BMO.Beam(bundle, i) for i in 1:length(bundle)]` for all rays and trace each beam with
`BMO.solve_system!`. The inverse is `RayBundle(detector)` after tracing to a `Detector`.
"""
function BMO.Beam(bundle::RayBundle, i::Integer)
    1 <= i <= length(bundle) ||
        throw(BoundsError(bundle, i))
    return BMO.Beam(bundle.position[i], bundle.direction[i], bundle.wavelength)
end

# ---------------------------------------------------------------------------------------
# Hit accessors

function _hits(detector::BMO.Detector)
    hs = BMO.hits(detector)
    (hs === nothing || isempty(hs)) &&
        throw(ArgumentError("$_PREFIX: the detector has no hits; run `solve_system!` first"))
    return hs
end

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
            throw(ArgumentError("$_PREFIX: all hits must share the wavelength, got $(BMO.wavelength(r)) m for hit $i and $λ m for hit 1"))
        isapprox(BMO.refractive_index(r), n_ray; rtol = tol) ||
            throw(ArgumentError("$_PREFIX: all hits must share the refractive index, got $(BMO.refractive_index(r)) for hit $i and $n_ray for hit 1"))
    end
    return λ, n_ray
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
        throw(ArgumentError("$_PREFIX: the given port must have the detector normal $(Vector(n)) (pointing downstream), got $(Vector(normal(port)))"))
    Δ = origin(port) - origin(ref)
    abs(dot(Δ, n)) <= tol * max(norm(Δ), λ) ||
        throw(ArgumentError("$_PREFIX: the origin of the given port must lie in the detector plane, got distance $(dot(Δ, n)) m"))
    isapprox(refractive_index(port), n_ray; rtol = tol) ||
        throw(ArgumentError("$_PREFIX: the refractive index of the given port must equal that of the detector hits ($n_ray), got $(refractive_index(port))"))
    return port
end

function _detector_port(::BMO.Detector, port, λ, n_ray)
    throw(ArgumentError("$_PREFIX: the port must be a PlanarPort, got $(typeof(port))"))
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
            throw(ArgumentError("$_PREFIX: hit $i propagates against the port normal (direction·normal = $(dot(d, n)) ≤ 0)"))
    end
    return positions, directions
end

end # module
