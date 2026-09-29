"""
    RayBundle{T <: Real, P <: AbstractPort} <: AbstractRayBundle

A bundle of rays crossing a port, geometry only: position and direction per ray and one
wavelength. It carries no amplitude, phase or polarization and is the format for handing
rays from one ray tracer to another. Rays with power, optical path length and
polarization are a [`PolarizedRayBundle`](@ref).

See the constructor `RayBundle(port, wavelength, position, direction)` for the meaning,
units and frames of all fields and for the invariants that are checked.

# Fields

  - `port::P`: the port at which the rays are given
  - `wavelength::T`: vacuum wavelength, \\[m\\]
  - `position::Vector{SVector{3,T}}`: ray positions on the port plane, global frame, \\[m\\]
  - `direction::Vector{SVector{3,T}}`: unit propagation directions, global frame

# Interface

Implements the [`AbstractRayBundle`](@ref) interface: [`port`](@ref),
[`wavelength`](@ref), [`positions`](@ref), [`directions`](@ref), `length(bundle)`.
"""
struct RayBundle{T <: Real, P <: AbstractPort} <: AbstractRayBundle
    port::P
    wavelength::T
    position::Vector{SVector{3, T}}
    direction::Vector{SVector{3, T}}

    function RayBundle{T, P}(port::P, wavelength::T, position::Vector{SVector{3, T}},
            direction::Vector{SVector{3, T}}) where {T <: Real, P <: AbstractPort}
        _check_rays("RayBundle", port, wavelength, position, direction)
        return new{T, P}(port, wavelength, position, direction)
    end
end

"""
    RayBundle(port, wavelength, position, direction)

Creates a bundle of `M` rays crossing `port`, geometry only. Monochromatic: rays of
several wavelengths are a collection of bundles.

# Arguments

  - `port`: the [`AbstractPort`](@ref) at which the rays are given, e.g. a
    [`PlanarPort`](@ref). The rays propagate downstream, into `direction·normal(port) > 0`,
    in the medium `refractive_index(port)`.
  - `wavelength`: vacuum wavelength in \\[m\\]; the wavelength in the medium is
    `wavelength / refractive_index(port)`.
  - `position`: `M` points in the global frame in \\[m\\], where the rays cross the port
    plane.
  - `direction`: `M` unit propagation directions in the global frame.

# Invariants

Each violation throws an `ArgumentError`, with `tol = √eps(T)`:

 1. `position` and `direction` have the same length; each element has 3 components.
 2. `wavelength > 0`.
 3. `|‖direction‖ − 1| ≤ tol`; `direction·normal(port) > 0`.
 4. The position lies on the port plane:
    `|(p − origin)·n| ≤ tol·max(‖p − origin‖, wavelength)`.

Inputs are not normalized silently. All real inputs are promoted to a common
floating-point type `T`; the storage is `Vector`s of `SVector`s (copied from the inputs).
"""
function RayBundle(port::AbstractPort, wavelength::Real, position::AbstractVector,
        direction::AbstractVector)
    T = float(promote_type(typeof(wavelength), _eltype_real(position),
        _eltype_real(direction)))
    return RayBundle{T, typeof(port)}(port, T(wavelength),
        _svectors3(T, position, "RayBundle", "position"),
        _svectors3(T, direction, "RayBundle", "direction"))
end

"""
    PolarizedRayBundle{T <: Real, P <: AbstractPort} <: AbstractRayBundle

A bundle of polarized rays crossing a port: one monochromatic, mutually coherent
component. Besides the geometry of a [`RayBundle`](@ref), every ray carries its optical
path length, its power and a unit polarization vector, which is what a converter to the
wave world needs (e.g. [`DebyeWolf`](@ref)).

See the constructor `PolarizedRayBundle(port, wavelength, position, direction, opl, power,
polarization)` for the meaning, units and frames of all fields and for the invariants that
are checked.

# Fields

  - `port::P`: the port at which the rays are given
  - `wavelength::T`: vacuum wavelength, \\[m\\]
  - `position::Vector{SVector{3,T}}`: ray positions on the port plane, global frame, \\[m\\]
  - `direction::Vector{SVector{3,T}}`: unit propagation directions, global frame
  - `opl::Vector{T}`: optical path length from the reference point of the component, \\[m\\]
  - `power::Vector{T}`: power per ray, \\[W\\]
  - `polarization::Vector{SVector{3,Complex{T}}}`: unit complex polarization vector per
    ray, global frame, transverse to the ray

# Interface

Implements the [`AbstractRayBundle`](@ref) interface: [`port`](@ref),
[`wavelength`](@ref), [`positions`](@ref), [`directions`](@ref), `length(bundle)`.
Further: [`total_power`](@ref).
"""
struct PolarizedRayBundle{T <: Real, P <: AbstractPort} <: AbstractRayBundle
    port::P
    wavelength::T
    position::Vector{SVector{3, T}}
    direction::Vector{SVector{3, T}}
    opl::Vector{T}
    power::Vector{T}
    polarization::Vector{SVector{3, Complex{T}}}

    function PolarizedRayBundle{T, P}(port::P, wavelength::T,
            position::Vector{SVector{3, T}}, direction::Vector{SVector{3, T}},
            opl::Vector{T}, power::Vector{T},
            polarization::Vector{SVector{3, Complex{T}}}) where {T <: Real, P <: AbstractPort}
        _check_rays("PolarizedRayBundle", port, wavelength, position, direction)
        _check_polarized_rays(direction, opl, power, polarization)
        return new{T, P}(port, wavelength, position, direction, opl, power, polarization)
    end
end

"""
    PolarizedRayBundle(port, wavelength, position, direction, opl, power, polarization)

Creates a bundle of `M` polarized rays crossing `port`: one monochromatic component whose
rays are mutually coherent. Polychromatic or mutually incoherent light is a collection of
bundles.

# Arguments

  - `port`, `wavelength`, `position`, `direction`: as for [`RayBundle`](@ref).
  - `opl`: `M` optical path lengths `∫ n ds` in \\[m\\], absolute from the reference point
    of the coherent component (set by the first solver of the chain), not from the port
    origin. The phase accumulated along the ray is `2π·opl/wavelength`.
  - `power`: `M` powers in \\[W\\] carried by each ray. The bundle power is the sum, see
    [`total_power`](@ref).
  - `polarization`: `M` complex unit 3-vectors `e` in the global frame, transverse to the
    ray (`Σ eᵢ dᵢ = 0`, no conjugation); build them from Jones vectors with
    [`jones_to_global`](@ref). `e` carries all phase not contained in `opl` (Fresnel and
    coating phases).

# Conventions

Field of ray `j` near its position, with time dependence exp(−iωt) not included:
`E ∝ √power · e · exp(i·2π·opl/wavelength)`.

# Invariants

The invariants of [`RayBundle`](@ref), and, with `tol = √eps(T)` (each violation throws an
`ArgumentError`):

 1. `opl`, `power` and `polarization` have one entry per ray.
 2. Every `power ≥ 0`.
 3. `|‖e‖ − 1| ≤ tol` and `|Σ eᵢ dᵢ| ≤ tol`.

Inputs are not normalized silently. All real inputs are promoted to a common
floating-point type `T`; the storage is `Vector`s (of `SVector`s), copied from the inputs.
"""
function PolarizedRayBundle(port::AbstractPort, wavelength::Real, position::AbstractVector,
        direction::AbstractVector, opl::AbstractVector, power::AbstractVector,
        polarization::AbstractVector)
    T = float(promote_type(typeof(wavelength), _eltype_real(position),
        _eltype_real(direction), _eltype_real(opl), _eltype_real(power),
        _eltype_real(polarization)))
    name = "PolarizedRayBundle"
    pol = SVector{3, Complex{T}}[_svector3(x, name, "polarization") for x in polarization]
    return PolarizedRayBundle{T, typeof(port)}(port, T(wavelength),
        _svectors3(T, position, name, "position"), _svectors3(T, direction, name, "direction"),
        Vector{T}(opl), Vector{T}(power), pol)
end

# Geometric invariants shared by all ray bundles, see the `RayBundle` constructor.
function _check_rays(name, port, wavelength::T, position, direction) where {T}
    tol = sqrt(eps(T))
    M = length(position)
    length(direction) == M ||
        throw(ArgumentError("$name: position and direction must have equal lengths, got $M and $(length(direction))"))
    wavelength > 0 ||
        throw(ArgumentError("$name: wavelength must be positive, got $wavelength"))
    n = normal(port)
    r0 = origin(port)
    for i in 1:M
        d = direction[i]
        abs(norm(d) - 1) <= tol ||
            throw(ArgumentError("$name: direction must be a unit vector, got norm $(norm(d)) for ray $i"))
        dot(d, n) > 0 ||
            throw(ArgumentError("$name: direction must propagate into the port (direction·normal > 0), got $(dot(d, n)) for ray $i"))
        Δ = position[i] - r0
        abs(dot(Δ, n)) <= tol * max(norm(Δ), wavelength) ||
            throw(ArgumentError("$name: position must lie on the port plane, got distance $(dot(Δ, n)) m for ray $i"))
    end
    return nothing
end

# Invariants of the per-ray data of `PolarizedRayBundle`.
function _check_polarized_rays(direction, opl, power::Vector{T}, polarization) where {T}
    tol = sqrt(eps(T))
    M = length(direction)
    lengths = (length(opl), length(power), length(polarization))
    all(==(M), lengths) ||
        throw(ArgumentError("PolarizedRayBundle: opl, power and polarization must have one entry per ray ($M), got $lengths"))
    for i in 1:M
        power[i] >= 0 ||
            throw(ArgumentError("PolarizedRayBundle: power must be non-negative, got $(power[i]) for ray $i"))
        e = polarization[i]
        abs(norm(e) - 1) <= tol ||
            throw(ArgumentError("PolarizedRayBundle: polarization must have unit norm, got norm $(norm(e)) for ray $i"))
        # Transversality without conjugation: Σ eᵢ dᵢ = 0 for real and imaginary part
        s = transpose(e) * direction[i]
        abs(s) <= tol ||
            throw(ArgumentError("PolarizedRayBundle: polarization must be transverse to the direction (polarization·direction = 0), got $s for ray $i"))
    end
    return nothing
end

# Real element type of possibly nested / complex inputs, for type promotion.
_eltype_real(x::Number) = real(typeof(x))
_eltype_real(::AbstractArray{S}) where {S <: Number} = real(S)
_eltype_real(::AbstractArray{<:AbstractArray{S}}) where {S <: Number} = real(S)
_eltype_real(x::AbstractArray) = mapreduce(_eltype_real, promote_type, x; init = Bool)

function _svector3(x, name, field)
    length(x) == 3 ||
        throw(ArgumentError("$name: each $field must have 3 elements, got $(length(x))"))
    return SVector{3}(x)
end

_svectors3(::Type{T}, xs, name, field) where {T} = SVector{3, T}[_svector3(x, name, field)
                                                               for x in xs]

Base.length(b::AbstractRayBundle) = length(positions(b))

port(b::Union{RayBundle, PolarizedRayBundle}) = b.port
wavelength(b::Union{RayBundle, PolarizedRayBundle}) = b.wavelength
positions(b::Union{RayBundle, PolarizedRayBundle}) = b.position
directions(b::Union{RayBundle, PolarizedRayBundle}) = b.direction

total_power(b::PolarizedRayBundle) = sum(b.power; init = zero(eltype(b.power)))
