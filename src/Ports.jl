"""
    PlanarPort{T <: Real} <: AbstractPort

A flat port: the plane through `origin` with unit normal `n`, equipped with the
right-handed local frame `(u, v, n)` in global coordinates.

# Fields

  - `origin::SVector{3,T}`: reference point in global coordinates, \\[m\\]
  - `axes::SMatrix{3,3,T,9}`: unit columns `(u, v, n)` in global coordinates
  - `refractive_index::T`: refractive index of the medium at the port

Use the accessors [`origin`](@ref), [`local_axes`](@ref), [`normal`](@ref) and
[`refractive_index`](@ref) instead of the fields.
"""
struct PlanarPort{T <: Real} <: AbstractPort
    origin::SVector{3, T}
    axes::SMatrix{3, 3, T, 9}
    refractive_index::T

    function PlanarPort{T}(origin::AbstractVector, normal::AbstractVector,
            u::AbstractVector, refractive_index::Real) where {T <: Real}
        _check_length3(origin, "origin")
        _check_length3(normal, "normal")
        _check_length3(u, "u")
        refractive_index > 0 ||
            throw(ArgumentError("PlanarPort: refractive_index must be positive, got $refractive_index"))
        n = SVector{3, T}(normal)
        norm_n = norm(n)
        norm_n > 0 || throw(ArgumentError("PlanarPort: normal must not be the zero vector"))
        n = n / norm_n
        u0 = SVector{3, T}(u)
        # Gram–Schmidt: remove the component of u along n
        u_perp = u0 - dot(u0, n) * n
        norm(u_perp) > sqrt(eps(T)) * norm(u0) ||
            throw(ArgumentError("PlanarPort: u must not be zero or parallel to the normal"))
        u_hat = u_perp / norm(u_perp)
        v_hat = cross(n, u_hat)
        return new{T}(SVector{3, T}(origin), hcat(u_hat, v_hat, n), T(refractive_index))
    end
end

_check_length3(x, name) = length(x) == 3 ||
                          throw(ArgumentError("PlanarPort: $name must have 3 elements, got $(length(x))"))

"""
    PlanarPort(origin, normal, u; refractive_index = 1)

Creates a flat port through `origin` with normal `normal` and first transverse axis `u`.

The local frame `(u, v, n)` is built as follows: `n = normal/‖normal‖`; `u` is made
orthogonal to `n` (Gram–Schmidt) and normalized; `v = n × u`, so the frame is right-handed
(`u × v = n`). The element type `T` is the floating-point promotion of all inputs.

# Arguments

  - `origin`: reference point in global coordinates, \\[m\\]. It is the origin of the local
    coordinates `(ξ, η, ζ)`, the center of every grid at this port and the reference point
    r₀ for position-dependent phase terms.
  - `normal`: normal in global coordinates, need not be normalized. It points downstream:
    fields at the port propagate into the half space `n·k > 0`.
  - `u`: first transverse axis in global coordinates, need not be normalized or orthogonal
    to `normal`. It fixes the polarization basis: Jones vectors at the port are given in
    `(u, v)`. There is deliberately no default.
  - `refractive_index`: refractive index of the medium at the port (default 1). It enters
    the power normalization, see [`power_normalization`](@ref).

Throws an `ArgumentError` if an input does not have three elements, if `normal` is the
zero vector, if `u` is zero or parallel to `normal` (`‖u − (u·n)n‖ ≤ √eps(T)·‖u‖`), or if
`refractive_index ≤ 0`.
"""
function PlanarPort(origin::AbstractVector, normal::AbstractVector, u::AbstractVector;
        refractive_index::Real = 1)
    T = float(promote_type(eltype(origin), eltype(normal), eltype(u), typeof(refractive_index)))
    return PlanarPort{T}(origin, normal, u, refractive_index)
end

"""
    origin(port::PlanarPort) -> SVector{3}

Returns the reference point of `port` in global coordinates, \\[m\\].
"""
origin(p::PlanarPort) = p.origin

"""
    local_axes(port::PlanarPort) -> SMatrix{3,3}

Returns the local frame of `port` as a 3×3 matrix with the unit columns `(u, v, n)` in
global coordinates. The frame is right-handed (`u × v = n`) and `n` points downstream.
"""
local_axes(p::PlanarPort) = p.axes

"""
    normal(port::PlanarPort) -> SVector{3}

Returns the unit normal `n` of `port` in global coordinates. It points downstream: fields at
the port propagate into the half space `n·k > 0`.
"""
normal(p::PlanarPort) = p.axes[:, 3]

"""
    refractive_index(port::PlanarPort) -> Real

Returns the refractive index of the medium at `port`.
"""
refractive_index(p::PlanarPort) = p.refractive_index

"""
    to_local(port::PlanarPort, r) -> SVector{3}

Converts the global point `r` in \\[m\\] to the local coordinates `(ξ, η, ζ)` of `port`
along `(u, v, n)`, relative to [`origin`](@ref)`(port)`. Points on the port plane have
`ζ = 0`.
"""
to_local(p::PlanarPort, r::AbstractVector) = transpose(p.axes) * (SVector{3}(r) - p.origin)

"""
    to_global(port::PlanarPort, ξ) -> SVector{3}

Converts local coordinates of `port` to a global point in \\[m\\]. `ξ` has either three
elements `(ξ, η, ζ)` along `(u, v, n)`, or two elements `(ξ, η)` for a point on the port
plane (`ζ = 0`). Inverse of [`to_local`](@ref).
"""
to_global(p::PlanarPort, ξ::StaticVector{3}) = p.origin + p.axes * ξ
to_global(p::PlanarPort, ξ::StaticVector{2}) = p.origin + p.axes[:, 1] * ξ[1] +
                                               p.axes[:, 2] * ξ[2]
function to_global(p::PlanarPort, ξ::AbstractVector)
    length(ξ) == 3 && return to_global(p, SVector{3}(ξ))
    length(ξ) == 2 && return to_global(p, SVector{2}(ξ))
    throw(ArgumentError("to_global: local coordinates must have 2 or 3 elements, got $(length(ξ))"))
end

"""
    ray_basis(port::PlanarPort, dir) -> (x̂, ŷ)

Returns the local transverse basis of a ray with propagation direction `dir` at `port`, as
two unit vectors in global coordinates. Jones vectors of oblique rays at the port are given
in this basis, see [`jones_to_global`](@ref).

The basis is the port basis `(u, v)` carried along by the minimal rotation R that turns the
port normal `n` into the ray direction `d = dir/‖dir‖`, i.e. the rotation about `n × d`:

```math
\\hat{x} = R\\,u, \\quad \\hat{y} = R\\,v, \\quad
R\\,x = c\\,x + k \\times x + \\frac{k \\cdot x}{1 + c}\\,k, \\quad c = n \\cdot d, \\; k = n \\times d.
```

# Conventions

  - `(x̂, ŷ, d)` is right-handed: `x̂ × ŷ = d`.
  - For `d ∥ n` the basis is exactly `(u, v)`.
  - This is the mapping of an aplanatic lens in the Richards–Wolf model: the meridional
    unit vector `ê_ρ` in the port plane maps to `ê_θ`, the azimuthal `ê_φ` stays unchanged.

`dir` need not be normalized. Throws an `ArgumentError` if `dir` is the zero vector or does
not propagate into the port (`d·n ≤ 0`).
"""
function ray_basis(p::PlanarPort, dir::AbstractVector)
    n = normal(p)
    u = p.axes[:, 1]
    v = p.axes[:, 2]
    d0 = SVector{3}(dir)
    norm_d = norm(d0)
    norm_d > 0 || throw(ArgumentError("ray_basis: dir must not be the zero vector"))
    dot(n, d0) > 0 ||
        throw(ArgumentError("ray_basis: dir must propagate into the port, i.e. dir·normal > 0"))
    # Exactly the port basis for rays along the normal
    iszero(cross(n, d0)) && return (u + zero(d0), v + zero(d0))
    d = d0 / norm_d
    c = dot(n, d)
    k = cross(n, d)
    rotate(x) = c * x + cross(k, x) + (dot(k, x) / (1 + c)) * k
    return (rotate(u), rotate(v))
end

"""
    jones_to_global(port::PlanarPort, dir, jones) -> SVector{3,Complex}

Converts the Jones vector `jones = (J₁, J₂)` of a ray with direction `dir` at `port` into a
complex 3D field vector in global coordinates: `J₁ x̂ + J₂ ŷ` with
`(x̂, ŷ) = ray_basis(port, dir)`. The result is transverse to `dir`. Units are those of
`jones`, e.g. \\[V/m\\] for amplitudes or dimensionless for unit polarization states.

For `dir ∥ normal(port)` the Jones vector is given in the port basis `(u, v)`.
"""
function jones_to_global(p::PlanarPort, dir::AbstractVector, jones::AbstractVector)
    length(jones) == 2 ||
        throw(ArgumentError("jones_to_global: jones must have 2 elements, got $(length(jones))"))
    x, y = ray_basis(p, dir)
    return complex.(jones[1] * x + jones[2] * y)
end

"""
    circular_jones(helicity::Integer) -> SVector{2,ComplexF64}

Returns the normalized Jones vector of circular polarization with the given `helicity`.

  - `helicity = +1`: right-circular = positive helicity, `(1, i)/√2`
  - `helicity = −1`: left-circular = negative helicity, `(1, −i)/√2`

# Conventions

Helicity convention (IEEE): with the time dependence exp(−iωt), the real field
`Re(J·exp(−iωt))` of `(1, i)/√2` rotates from the first to the second basis vector, i.e.
counter-clockwise about the propagation direction for an observer facing the source. The
Jones vector must be given in a right-handed basis `(x̂, ŷ)` with `x̂ × ŷ` along the
propagation direction, such as the port basis `(u, v)` or [`ray_basis`](@ref).

Throws an `ArgumentError` for any helicity other than ±1.
"""
function circular_jones(helicity::Integer)
    helicity == 1 || helicity == -1 ||
        throw(ArgumentError("circular_jones: helicity must be +1 or -1, got $helicity"))
    return SVector{2, ComplexF64}(1, helicity * im) / sqrt(2)
end
