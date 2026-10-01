"""
    OpticsBase

Minimal interface package for coupling optical solvers, i.e. numerical solutions of
Maxwell's equations. There is a single exchange format, [`PlaneField`](@ref): the
tangential electric and magnetic field sampled on a plane in 3D space. By the surface
equivalence theorem it determines the field on both sides of the plane, so every solver
(beamlets, BPM, angular spectrum, FDTD, FEM, RCWA, ...) can produce and consume it.

`OpticsBase` contains no solver code and no solver interface. Solvers couple by function
composition: a package offers methods that return a `PlaneField` and methods that take
one.

All quantities are in SI units. See the "Conventions" page of the documentation for the
binding time, phase, power and polarization conventions.
"""
module OpticsBase

using LinearAlgebra: I, det, norm
using StaticArrays: SMatrix, SVector

export PlaneField, VACUUM_IMPEDANCE
export coordinates, reference_phase, power, forward, backward

"""
    VACUUM_IMPEDANCE

Impedance of free space `Z₀ = μ₀ c = 376.730313668` Ω (CODATA 2018, the same literal as
BeamletOptics' `Z_vacuum`).
"""
const VACUUM_IMPEDANCE = 376.730313668

"""
    PlaneField(E, H, spacing, origin, axes, λ; n = 1, R = Inf)
    PlaneField(E, spacing, origin, axes, λ; n = 1, R = Inf)

Monochromatic, coherent electromagnetic field given by its tangential components `E` and
`H` on a regularly sampled plane in 3D space. This is the only exchange format of
OpticsBase: solvers hand fields to each other as `PlaneField`s.

The first form takes both fields and is exact for any field. The second form takes `E`
only and fills in `H` for a local plane wave along the reference direction of each sample
(`n` if `R = Inf`, otherwise the ray of the reference sphere through the sample, see
[`reference_phase`](@ref)), with `E` completed to be transverse to that direction. It is
exact for a wave travelling along the reference direction; for a plane wave at the angle
`θ` to it, `H` is off by the order of `1 − cos θ`. With an `nx × ny` matrix `E` the field
is scalar and is put into `Eu` (`Ev = 0`).

# Arguments

- `E`: `nx × ny × 2` array of `(Eu, Ev)`, the electric field along the plane axes `u`,
  `v`, in \\[V/m\\] (peak amplitude, not RMS).
- `H`: `nx × ny × 2` array of `(Hu, Hv)` in \\[A/m\\] (peak amplitude).
- `spacing`: sample spacings `(Δu, Δv)` in \\[m\\]. Sample `(i, j)` sits at
  `(coordinates(f, 1)[i], coordinates(f, 2)[j])`, with the origin at the fftshift center
  (see [`coordinates`](@ref)).
- `origin`: center of the plane in global coordinates in \\[m\\].
- `axes`: `3 × 3` matrix with the columns `u`, `v`, `n`. Orthonormal to within `1e-6`
  (axes computed in single precision are accepted) and right-handed (`u × v = n`). `n`
  is the positive normal: [`power`](@ref) counts the flux along `+n` as positive, and
  [`forward`](@ref) is the part travelling towards the `+n` side. `u` fixes the
  polarization basis.
- `λ`: vacuum wavelength in \\[m\\].
- `n`: real refractive index of the medium at the plane.
- `R`: radius of the reference sphere in \\[m\\], see [`reference_phase`](@ref). `Inf`:
  no reference sphere. `R > 0`: diverging from the point `origin − R n`; `R < 0`:
  converging to `origin + |R| n`.

# Medium

The plane must lie in a homogeneous, isotropic, lossless and non-magnetic medium of index
`n`, at least across the sampled region. `E` and `H` themselves are exact in any medium,
but `n` enters the E-only constructor, [`forward`](@ref)/[`backward`](@ref),
[`reference_phase`](@ref) and the normal components below. A plane through an
inhomogeneous structure, e.g. a fiber cross-section, has no single `n`: place it in the
surrounding medium instead.

# Conventions

- Time convention `exp(−iωt)`; a plane wave is `exp(i(k·r − ωt))`.
- The physical fields are `E .* reference_phase(f)` and `H .* reference_phase(f)`. They
  contain the full spatial phase, including the absolute optical path from the reference
  point of the chain (set by the first solver, typically the source).
- The normal components are not stored. They follow from Maxwell's curl equations on the
  plane (non-magnetic medium): `En = i Z₀/(k₀ n²) (∂u Hv − ∂v Hu)` and
  `Hn = −i/(k₀ Z₀) (∂u Ev − ∂v Eu)`, with `k₀ = 2π/λ`, applied to the physical fields.
- A field that does not fit a single plane is a `Vector{PlaneField}`: a closed Huygens
  surface (normals outward), a spectrum (one field per line, common time origin `t = 0`)
  or mutually incoherent components (powers add). The vector has no meaning of its own;
  the function that returns or takes it documents which one it is. See the
  "Conventions" page.

# Fields

- `E`, `H`, `spacing`, `origin`, `axes`, `λ`, `n`, `R`: as above. `origin` is an
  `SVector{3}`, `axes` an `SMatrix{3, 3}`. The fields are the public API; there are no
  accessor functions.

Arrays and geometry have independent precision. The geometry type `T` (of `spacing`,
`origin`, `axes`, `λ`, `n`, `R`) follows from `spacing`, `origin`, `axes` and `λ`, so a
`Float32` field can sit at a `Float64` position. `E` and `H` keep their own (complex,
floating point) element type and are stored in the array type `A` of `similar(E)` (e.g.
`Array` or a GPU array); an argument that already has that type is stored without a copy,
any other (e.g. a view) is copied.
"""
struct PlaneField{T <: Real, A <: AbstractArray{<:Complex, 3}}
    E::A
    H::A
    spacing::NTuple{2, T}
    origin::SVector{3, T}
    axes::SMatrix{3, 3, T, 9}
    λ::T
    n::T
    R::T
    function PlaneField{T, A}(E, H, spacing, origin, axes, λ, n, R) where {T, A}
        size(E) == size(H) ||
            throw(DimensionMismatch("E has size $(size(E)), H has size $(size(H))"))
        size(E, 3) == 2 ||
            throw(DimensionMismatch("E and H have 2 components (along u and v), got $(size(E, 3))"))
        all(Δ -> isfinite(Δ) && Δ > 0, spacing) ||
            throw(ArgumentError("spacing must be positive and finite, got $spacing"))
        all(isfinite, origin) || throw(ArgumentError("origin must be finite, got $origin"))
        # 1e-6 admits axes computed in single precision; a tilt that small is negligible
        norm(axes' * axes - I) <= 1e-6 ||
            throw(ArgumentError("the columns of axes (u, v, n) must be orthonormal to within 1e-6"))
        det(axes) > 0 || throw(ArgumentError("axes (u, v, n) must be right-handed"))
        isfinite(λ) && λ > 0 || throw(ArgumentError("λ must be positive and finite, got $λ"))
        isfinite(n) && n > 0 || throw(ArgumentError("n must be positive and finite, got $n"))
        !isnan(R) && R != 0 ||
            throw(ArgumentError("R must be nonzero; use R = Inf for no reference sphere"))
        return new{T, A}(E, H, spacing, origin, axes, λ, n, R)
    end
end

function PlaneField(E::AbstractArray{<:Number, 3}, H::AbstractArray{<:Number, 3},
        spacing, origin, axes, λ; n = 1, R = Inf)
    T = float(promote_type(eltype(spacing), eltype(origin), eltype(axes), typeof(λ)))
    C = complex(float(promote_type(eltype(E), eltype(H))))
    Es = _storage(C, E, E)
    Hs = _storage(C, H, E)
    return PlaneField{T, typeof(Es)}(Es, Hs, T.(Tuple(spacing)), SVector{3, T}(origin),
        SMatrix{3, 3, T}(axes), T(λ), T(n), T(R))
end

function PlaneField(E::AbstractArray{<:Number, 3}, spacing, origin, axes, λ; n = 1, R = Inf)
    size(E, 3) == 2 ||
        throw(DimensionMismatch("E has 2 components (along u and v), got $(size(E, 3))"))
    H = similar(E, complex(float(eltype(E))))
    Rt = real(eltype(H))
    ξ, η = _coordinates(Rt, size(E), spacing)
    Y, Rr = Rt(n / VACUUM_IMPEDANCE), Rt(R)
    Eu, Ev = view(E, :, :, 1), view(E, :, :, 2)
    H[:, :, 1] .= _admittance_u.(Eu, Ev, ξ, η', Rr, Y)
    H[:, :, 2] .= _admittance_v.(Eu, Ev, ξ, η', Rr, Y)
    return PlaneField(E, H, spacing, origin, axes, λ; n, R)
end

function PlaneField(E::AbstractMatrix{<:Number}, spacing, origin, axes, λ; kwargs...)
    E2 = similar(E, size(E)..., 2)
    E2[:, :, 1] .= E
    E2[:, :, 2] .= 0
    return PlaneField(E2, spacing, origin, axes, λ; kwargs...)
end

# `X` as an array of element type `C` and the array type of `similar(prototype)`; no copy
# if it already is one.
function _storage(::Type{C}, X, prototype) where {C}
    target = typeof(similar(prototype, C, ntuple(_ -> 0, ndims(prototype))))
    return X isa target ? X : copyto!(similar(prototype, C, size(X)), X)
end

# Sample coordinates along u and v in the real type `Rt` (see `coordinates`)
function _coordinates(::Type{Rt}, sz, spacing) where {Rt}
    return ntuple(d -> ((0:(sz[d] - 1)) .- sz[d] ÷ 2) .* Rt(spacing[d]), 2)
end

# Local reference direction (du, dv, dn) at the sample (ξ, η): the plane normal for
# R = Inf, otherwise the ray of the reference sphere through the sample.
@inline function _reference_direction(ξ, η, R)
    isinf(R) && return (zero(ξ), zero(ξ), one(ξ))
    r = sqrt(ξ^2 + η^2 + R^2)
    return (sign(R) * ξ / r, sign(R) * η / r, abs(R) / r)
end

# Tangential H of a plane wave along the reference direction d whose tangential E is
# (Eu, Ev), with E completed to be transverse to d: Ht = M Et, M = Y/dn [−du dv, −(1 − du²);
# 1 − dv², du dv]. A wave with the same transverse direction travelling backward has −M.
@inline function _admittance_u(Eu, Ev, ξ, η, R, Y)
    du, dv, dn = _reference_direction(ξ, η, R)
    return -Y / dn * (du * dv * Eu + (1 - du^2) * Ev)
end
@inline function _admittance_v(Eu, Ev, ξ, η, R, Y)
    du, dv, dn = _reference_direction(ξ, η, R)
    return Y / dn * ((1 - dv^2) * Eu + du * dv * Ev)
end

# The inverse, Et = M⁻¹ Ht (det M = Y²)
@inline function _impedance_u(Hu, Hv, ξ, η, R, Y)
    du, dv, dn = _reference_direction(ξ, η, R)
    return (du * dv * Hu + (1 - du^2) * Hv) / (Y * dn)
end
@inline function _impedance_v(Hu, Hv, ξ, η, R, Y)
    du, dv, dn = _reference_direction(ξ, η, R)
    return -((1 - dv^2) * Hu + du * dv * Hv) / (Y * dn)
end

"""
    coordinates(f::PlaneField, d)

Sample positions along `u` (`d = 1`) or `v` (`d = 2`) in \\[m\\], relative to
`f.origin`: sample `i` sits at `(i − 1 − N÷2) Δ`, with `N` the number of samples and `Δ`
the spacing along that axis. The origin is the fftshift center, so there is always a
sample at `0`. The global position of sample `(i, j)` is
`f.origin + u coordinates(f, 1)[i] + v coordinates(f, 2)[j]`.
"""
coordinates(f::PlaneField{T}, d::Integer) where {T} = _coordinates(T, size(f.E), f.spacing)[d]

"""
    reference_phase(f::PlaneField)

`nx × ny` array of the phase factor of the reference sphere,
`exp(i s k₀ n (√(ρ² + R²) − |R|))` with `s = sign(R)`, `ρ` the distance from
`f.origin` in the plane and `k₀ = 2π/λ`. Multiplying `f.E` and `f.H` by it gives the
physical fields. The factor is `1` at the origin and everywhere if `R = Inf`. It has the
element type and array type of `f.E`.

The reference sphere keeps curved wavefronts cheap to sample: only the deviation from the
sphere is stored. A beam converging to a focus at the distance `F` behind the plane has
`R = −F`, one diverging from a point at the distance `F` before it has `R = F`. The rays
of the sphere are also the local reference directions of the E-only constructor and of
[`forward`](@ref)/[`backward`](@ref).
"""
function reference_phase(f::PlaneField)
    C = eltype(f.E)
    P = similar(f.E, C, size(f.E, 1), size(f.E, 2))
    isinf(f.R) && return fill!(P, one(C))
    Rt = real(C)
    ξ, η = _coordinates(Rt, size(f.E), f.spacing)
    a = Rt(sign(f.R) * 2π / f.λ * f.n)
    R = Rt(f.R)
    # √(ρ² + R²) − |R| written without cancellation for ρ ≪ |R|
    P .= cis.(a .* (ξ .^ 2 .+ η' .^ 2) ./ (sqrt.(ξ .^ 2 .+ η' .^ 2 .+ R^2) .+ abs(R)))
    return P
end

"""
    power(f::PlaneField)

Net power through the plane in \\[W\\]: the time-averaged Poynting flux along `n`,
`P = ½ Re ∫ (E × H*)·n dA = ½ Re ∫ (Eu Hv* − Ev Hu*) dA`, as a Riemann sum over the
samples. Exact, with no paraxial factor; negative if more power travels along `−n` than
along `n`. The reference phase cancels, so `R` does not enter.
"""
function power(f::PlaneField)
    Eu, Ev = view(f.E, :, :, 1), view(f.E, :, :, 2)
    Hu, Hv = view(f.H, :, :, 1), view(f.H, :, :, 2)
    return sum(real.(Eu .* conj.(Hv) .- Ev .* conj.(Hu))) * prod(f.spacing) / 2
end

"""
    forward(f::PlaneField)
    backward(f::PlaneField)

Part of `f` travelling forward (along the reference direction, `dn > 0`) or backward, as a
[`PlaneField`](@ref) on the same plane. Each sample is split as a pair of local plane
waves with the transverse direction of the reference direction (see the E-only
constructor of [`PlaneField`](@ref)): with the tangential admittance `M` of the forward
wave (`Ht = M Et`; the backward wave has `−M`), `E± = ½ (Et ± M⁻¹ Ht)` and `H± = ±M E±`.
For `R = Inf` the reference direction is `n` and this is `E± = ½ (Et ∓ (Z₀/n) n × Ht)`.
`forward(f) + backward(f)` restores `f`, and the powers of the parts add up to `power(f)`.

The split is exact for waves along the reference direction. A plane wave at the angle `θ`
to it leaks a spurious backward amplitude of the order `(1 − cos θ)/2` (exactly that for
s-polarization), so `forward` is meant for fields that may contain both directions. A
consumer that knows its input travels one way should take `f.E` directly. An exact split
per plane wave of the angular spectrum needs an FFT and belongs in a solver.
"""
forward(f::PlaneField) = _split(f, 1)

"""
    backward(f::PlaneField)

Part of `f` travelling backward; see [`forward`](@ref).
"""
backward(f::PlaneField) = _split(f, -1)

function _split(f::PlaneField{T}, s) where {T}
    Rt = real(eltype(f.E))
    ξ, η = _coordinates(Rt, size(f.E), f.spacing)
    Y, R = Rt(f.n / VACUUM_IMPEDANCE), Rt(f.R)
    Eu, Ev = view(f.E, :, :, 1), view(f.E, :, :, 2)
    Hu, Hv = view(f.H, :, :, 1), view(f.H, :, :, 2)
    E = similar(f.E)
    H = similar(f.H)
    E[:, :, 1] .= (Eu .+ s .* _impedance_u.(Hu, Hv, ξ, η', R, Y)) ./ 2
    E[:, :, 2] .= (Ev .+ s .* _impedance_v.(Hu, Hv, ξ, η', R, Y)) ./ 2
    Esu, Esv = view(E, :, :, 1), view(E, :, :, 2)
    H[:, :, 1] .= s .* _admittance_u.(Esu, Esv, ξ, η', R, Y)
    H[:, :, 2] .= s .* _admittance_v.(Esu, Esv, ξ, η', R, Y)
    return PlaneField{T, typeof(E)}(E, H, f.spacing, f.origin, f.axes, f.λ, f.n, f.R)
end

end # module OpticsBase
