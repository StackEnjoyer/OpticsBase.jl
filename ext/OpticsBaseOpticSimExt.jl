# OpticSim-side glue between OpticSim.jl and OpticsBase. This is OpticSim's code, living
# here temporarily (like the BeamletOptics glue): it reads OpticSim rays into a
# `RayBundle` at a port and builds OpticSim rays from a `RayBundle`.
#
# Units: OpticSim works in mm (lengths) and µm (wavelengths), OpticsBase in m. Every
# conversion happens explicitly in this file. Targets OpticSim 0.7.1 (main branch, see
# test/integration/Project.toml for the pinned commit); the integration tests detect
# breaking OpticSim changes. OpticSim itself is not changed.

module OpticsBaseOpticSimExt

using OpticsBase: OpticsBase, AbstractPort, RayBundle, normal, origin
using LinearAlgebra: dot, norm
using StaticArrays: SVector
import OpticSim

const _PREFIX = "OpticSim rays"

# Unit conversion factors: OpticSim length unit (mm) and wavelength unit (µm) in m.
const _M_PER_MM = 1 // 1000
const _M_PER_UM = 1 // 1_000_000

# Ray types accepted by `RayBundle`: geometric rays, and rays that carry a wavelength.
const _OpticalRayLike = Union{OpticSim.OpticalRay{<:Real, 3}, OpticSim.LensTrace{<:Real, 3}}
const _AnyRay = Union{OpticSim.Ray{<:Real, 3}, _OpticalRayLike}

"""
    RayBundle(rays::AbstractVector, port::AbstractPort, wavelength::Real)
    RayBundle(rays::AbstractVector, port::AbstractPort)

Reads OpticSim.jl rays into a geometry-only [`RayBundle`](@ref) at `port`, for handing
rays from OpticSim to another ray tracer (e.g. BeamletOptics). Requires `using OpticSim`
(package extension). This is OpticSim-side glue living in OpticsBase temporarily; it will
move into OpticSim.

# Arguments

  - `rays`: OpticSim rays, all in the medium of the port: `OpticSim.Ray`s,
    `OpticSim.OpticalRay`s, or the `OpticSim.LensTrace`s returned by `OpticSim.trace`
    (their final ray is used). Filter out missed rays (`trace` returns `nothing`) first.
  - `port`: the [`AbstractPort`](@ref) at which the bundle is given, in \\[m\\] in the
    OpticsBase global frame. Its `refractive_index` must be the index of the OpticSim
    medium the rays are in (not checked; OpticSim rays do not carry it).
  - `wavelength`: vacuum wavelength in \\[m\\]. Required for geometric `OpticSim.Ray`s.
    Optional for `OpticalRay`s and `LensTrace`s, which carry their wavelength in \\[µm\\]:
    if omitted, it is read from the rays (all must agree within `√eps`); if given, it must
    agree with theirs within `√eps`.

# Conventions

  - Units: OpticSim positions in \\[mm\\] are converted to \\[m\\] (factor 10⁻³),
    OpticSim wavelengths in \\[µm\\] to \\[m\\] (factor 10⁻⁶). Directions are unit vectors
    and need no conversion.
  - Frame: the OpticSim global frame is the OpticsBase global frame (only the length unit
    differs). OpticSim's own conventions (e.g. light along −z in its lens constructors) do
    not apply to OpticsBase; the port defines the frame of the handover.
  - Position: the straight line of each ray (origin, direction) is intersected with the
    port plane, `p = o + t·d` with `t = (origin(port) − o)·n / (d·n)`. The rays are only
    carried forward (`t ≥ 0` within `√eps` relative); a ray whose origin already lies
    beyond the port plane throws. For a `LensTrace` the origin is the last surface hit, so
    the port may be anywhere downstream of it in the same medium.
  - Wavelength: OpticSim's wavelength (documented by OpticSim as the wavelength in air,
    with air index 1 by definition) is taken as the OpticsBase vacuum wavelength; the
    difference (air index ≈ 1.0003) is not converted.
  - Only geometry is transferred: OpticSim's power and optical path length are dropped.

Throws an `ArgumentError` if `rays` is empty, if a ray propagates against the port normal
(`d·n ≤ 0`), lies beyond the port plane, or if the wavelengths disagree. The
[`RayBundle`](@ref) constructor checks its invariants on top.
"""
function OpticsBase.RayBundle(rays::AbstractVector{<:_AnyRay}, port::AbstractPort,
        wavelength::Real)
    isempty(rays) && throw(ArgumentError("$_PREFIX: no rays given"))
    λ = wavelength
    for (i, r) in enumerate(rays)
        λr = _wavelength(r)
        λr === nothing && continue
        tol = sqrt(eps(float(typeof(λr))))
        isapprox(λr, λ; rtol = tol) ||
            throw(ArgumentError("$_PREFIX: ray $i has the wavelength $λr m, the given wavelength is $λ m"))
    end
    return _bundle(rays, port, λ)
end

function OpticsBase.RayBundle(rays::AbstractVector{<:_OpticalRayLike}, port::AbstractPort)
    isempty(rays) && throw(ArgumentError("$_PREFIX: no rays given"))
    λ = _wavelength(first(rays))
    return OpticsBase.RayBundle(rays, port, λ)
end

"""
    OpticSim.OpticalRay(bundle::RayBundle, i::Integer; power = 1)

Builds the OpticSim `OpticalRay` of ray `i` of a [`RayBundle`](@ref), for handing rays from
another ray tracer (e.g. BeamletOptics) into OpticSim. Requires `using OpticSim`. This is
OpticSim-side glue living in OpticsBase temporarily; it will move into OpticSim.

The ray starts at `position[i]` converted from \\[m\\] to OpticSim's \\[mm\\] (factor 10³),
with direction `direction[i]` (unit vector, global frame, unchanged) and wavelength
`wavelength(bundle)` converted from \\[m\\] to OpticSim's \\[µm\\] (factor 10⁶; taken as
OpticSim's wavelength in air without converting vacuum to air). Its optical path length is
0 and its `power` (OpticSim's relative ray power, dimensionless) defaults to 1: a
`RayBundle` carries neither. The ray starts in the medium of the port; build the OpticSim
system accordingly.

Use `[OpticSim.OpticalRay(bundle, i) for i in 1:length(bundle)]` for all rays and trace
each with `OpticSim.trace`. The inverse is `RayBundle(rays, port)`.
"""
function OpticSim.OpticalRay(bundle::RayBundle, i::Integer; power::Real = 1)
    1 <= i <= length(bundle) || throw(BoundsError(bundle, i))
    T = eltype(eltype(bundle.position))
    o = bundle.position[i] ./ T(_M_PER_MM)
    d = bundle.direction[i]
    λ = bundle.wavelength / T(_M_PER_UM)
    return OpticSim.OpticalRay(o, d, T(power), λ)
end

# ---------------------------------------------------------------------------------------
# Ray accessors (OpticSim units in, SI out)

_geometric_ray(r::OpticSim.Ray) = r
_geometric_ray(r::OpticSim.OpticalRay) = OpticSim.ray(r)
_geometric_ray(r::OpticSim.LensTrace) = OpticSim.ray(OpticSim.ray(r))

_wavelength(::OpticSim.Ray) = nothing
_wavelength(r::OpticSim.OpticalRay) = OpticSim.wavelength(r) * _M_PER_UM
_wavelength(r::OpticSim.LensTrace) = OpticSim.wavelength(r) * _M_PER_UM

_origin_m(r) = SVector{3}(OpticSim.origin(_geometric_ray(r))) * _M_PER_MM
_direction(r) = SVector{3}(OpticSim.direction(_geometric_ray(r)))

# Line-plane intersection of every ray with the port plane, forward only.
function _bundle(rays, port::AbstractPort, λ)
    n = normal(port)
    r0 = origin(port)
    T = float(eltype(r0))
    positions = Vector{SVector{3, T}}(undef, length(rays))
    directions = Vector{SVector{3, T}}(undef, length(rays))
    for (i, r) in enumerate(rays)
        o = _origin_m(r)
        d = _direction(r)
        dn = dot(d, n)
        dn > 0 ||
            throw(ArgumentError("$_PREFIX: ray $i propagates against the port normal (direction·normal = $dn ≤ 0)"))
        t = dot(r0 - o, n) / dn
        tol = sqrt(eps(float(typeof(t))))
        t >= -tol * max(norm(o - r0), λ) ||
            throw(ArgumentError("$_PREFIX: ray $i starts beyond the port plane (distance $(-t) m along the ray); rays are only carried forward"))
        positions[i] = o + t * d
        directions[i] = d
    end
    return RayBundle(port, λ, positions, directions)
end

end # module
