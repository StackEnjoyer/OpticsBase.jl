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

The first form takes both fields. The second form takes `E` only and fills in `H` for a
local plane wave travelling along `n` in the medium: `H = (n/Z₀) n × E`, i.e.
`Hu = −(n/Z₀) Ev`, `Hv = (n/Z₀) Eu`. It is exact for waves along the normal and a
paraxial approximation otherwise. With an `nx × ny` matrix `E` the field is scalar and is
put into `Eu` (`Ev = 0`).

# Arguments

- `E`: `nx × ny × 2` array of `(Eu, Ev)`, the electric field along the plane axes `u`,
  `v`, in \\[V/m\\] (peak amplitude, not RMS).
- `H`: `nx × ny × 2` array of `(Hu, Hv)` in \\[A/m\\] (peak amplitude).
- `spacing`: sample spacings `(Δu, Δv)` in \\[m\\]. Sample `(i, j)` sits at
  `(coordinates(f, 1)[i], coordinates(f, 2)[j])`, with the origin at the fftshift center
  (see [`coordinates`](@ref)).
- `origin`: center of the plane in global coordinates in \\[m\\].
- `axes`: `3 × 3` matrix with the columns `u`, `v`, `n`. Orthonormal and right-handed
  (`u × v = n`). `n` is the reference direction: [`power`](@ref) counts the flux along
  `n` positive, [`forward`](@ref) is the part travelling along `n`. `u` fixes the
  polarization basis.
- `λ`: vacuum wavelength in \\[m\\].
- `n`: real refractive index of the homogeneous, isotropic medium at the plane.
- `R`: radius of the reference sphere in \\[m\\], see [`reference_phase`](@ref). `Inf`:
  no reference sphere. `R > 0`: diverging from the point `origin − R n`; `R < 0`:
  converging to `origin + |R| n`.

# Conventions

- Time convention `exp(−iωt)`; a plane wave is `exp(i(k·r − ωt))`.
- The physical fields are `E .* reference_phase(f)` and `H .* reference_phase(f)`. They
  contain the full spatial phase, including the absolute optical path from the reference
  point of the chain (set by the first solver, typically the source).
- The normal components are not stored. They follow from Maxwell's curl equations on the
  plane: `En = i Z₀/(k₀ n²) (∂u Hv − ∂v Hu)` and `Hn = −i/(k₀ Z₀) (∂u Ev − ∂v Eu)`, with
  `k₀ = 2π/λ`, applied to the physical fields.
- A field that does not fit a single plane is a `Vector{PlaneField}` by convention: the
  faces of a closed Huygens box, the wavelengths of a pulse (common time origin `t = 0`),
  or mutually incoherent components (powers add).

# Fields

- `E`, `H`, `spacing`, `origin`, `axes`, `λ`, `n`, `R`: as above. `origin` is an
  `SVector{3}`, `axes` an `SMatrix{3, 3}`. The fields are the public API; there are no
  accessor functions.

The element type `T` follows from `E`, `H`, `spacing`, `origin`, `axes` and `λ`; `n` and
`R` are converted to it. `E` and `H` are stored in the array type `A` of `similar(E)` with
element type `Complex{T}` (e.g. `Array` or a GPU array); an argument that already has
that type is stored without a copy, any other (e.g. a view) is copied.
"""
struct PlaneField{T <: Real, A <: AbstractArray{Complex{T}, 3}}
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
        all(>(0), spacing) || throw(ArgumentError("spacing must be positive, got $spacing"))
        tol = 10 * sqrt(eps(float(T)))
        norm(axes' * axes - I) <= tol ||
            throw(ArgumentError("the columns of axes (u, v, n) must be orthonormal"))
        det(axes) > 0 || throw(ArgumentError("axes (u, v, n) must be right-handed"))
        λ > 0 || throw(ArgumentError("λ must be positive, got $λ"))
        n > 0 || throw(ArgumentError("n must be positive, got $n"))
        R != 0 || throw(ArgumentError("R must be nonzero; use R = Inf for no reference sphere"))
        return new{T, A}(E, H, spacing, origin, axes, λ, n, R)
    end
end

function PlaneField(E::AbstractArray{<:Number, 3}, H::AbstractArray{<:Number, 3},
        spacing, origin, axes, λ; n = 1, R = Inf)
    # n and R do not widen T, so their defaults keep Float32 fields Float32
    T = float(promote_type(real(eltype(E)), real(eltype(H)), eltype(spacing),
        eltype(origin), eltype(axes), typeof(λ)))
    Es = _storage(T, E, E)
    Hs = _storage(T, H, E)
    return PlaneField{T, typeof(Es)}(Es, Hs, T.(Tuple(spacing)), SVector{3, T}(origin),
        SMatrix{3, 3, T}(axes), T(λ), T(n), T(R))
end

function PlaneField(E::AbstractArray{<:Number, 3}, spacing, origin, axes, λ; n = 1, R = Inf)
    size(E, 3) == 2 ||
        throw(DimensionMismatch("E has 2 components (along u and v), got $(size(E, 3))"))
    H = similar(E, complex(float(eltype(E))))
    Y = real(eltype(H))(n / VACUUM_IMPEDANCE)
    H[:, :, 1] .= -Y .* view(E, :, :, 2)
    H[:, :, 2] .= Y .* view(E, :, :, 1)
    return PlaneField(E, H, spacing, origin, axes, λ; n, R)
end

function PlaneField(E::AbstractMatrix{<:Number}, spacing, origin, axes, λ; kwargs...)
    E2 = similar(E, size(E)..., 2)
    E2[:, :, 1] .= E
    E2[:, :, 2] .= 0
    return PlaneField(E2, spacing, origin, axes, λ; kwargs...)
end

# `X` as an array of `Complex{T}` of the array type of `similar(prototype)`; no copy if it
# already is one.
function _storage(::Type{T}, X, prototype) where {T}
    Y = similar(prototype, Complex{T}, size(X))
    return X isa typeof(Y) ? X : copyto!(Y, X)
end

"""
    coordinates(f::PlaneField, d)

Sample positions along `u` (`d = 1`) or `v` (`d = 2`) in \\[m\\], relative to
`f.origin`: sample `i` sits at `(i − 1 − N÷2) Δ`, with `N` the number of samples and `Δ`
the spacing along that axis. The origin is the fftshift center, so there is always a
sample at `0`. The global position of sample `(i, j)` is
`f.origin + u coordinates(f, 1)[i] + v coordinates(f, 2)[j]`.
"""
function coordinates(f::PlaneField, d::Integer)
    N = size(f.E, d)
    return ((0:(N - 1)) .- N ÷ 2) .* f.spacing[d]
end

"""
    reference_phase(f::PlaneField)

`nx × ny` array of the phase factor of the reference sphere,
`exp(i s k₀ n (√(ρ² + R²) − |R|))` with `s = sign(R)`, `ρ` the distance from
`f.origin` in the plane and `k₀ = 2π/λ`. Multiplying `f.E` and `f.H` by it gives the
physical fields. The factor is `1` at the origin and everywhere if `R = Inf`.

The reference sphere keeps curved wavefronts cheap to sample: only the deviation from the
sphere is stored. A beam converging to a focus at the distance `F` behind the plane has
`R = −F`, one diverging from a point at the distance `F` before it has `R = F`.
"""
function reference_phase(f::PlaneField{T}) where {T}
    P = similar(f.E, Complex{T}, size(f.E, 1), size(f.E, 2))
    if isinf(f.R)
        return fill!(P, one(Complex{T}))
    end
    a = sign(f.R) * 2π / f.λ * f.n
    R2 = f.R^2
    u, v = coordinates(f, 1), coordinates(f, 2)
    # √(ρ² + R²) − |R| written without cancellation for ρ ≪ |R|
    P .= cis.(a .* (u .^ 2 .+ v' .^ 2) ./ (sqrt.(u .^ 2 .+ v' .^ 2 .+ R2) .+ abs(f.R)))
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

Part of `f` travelling along `+n` (`forward`) or `−n` (`backward`), as a
[`PlaneField`](@ref) on the same plane: `E± = ½ (E ∓ (Z₀/n) n × H)`, with `H` of each
part from the local plane wave rule `H± = ±(n/Z₀) n × E±`. `forward(f) + backward(f)`
restores `E`, and the powers of the parts add up to `power(f)`.

The split is exact for waves along `n` and a paraxial approximation otherwise: each
sample is treated as a local plane wave along `±n`. An exact split per plane wave of the
angular spectrum needs an FFT and belongs in a solver.
"""
forward(f::PlaneField) = _split(f, 1)

"""
    backward(f::PlaneField)

Part of `f` travelling along `−n`; see [`forward`](@ref).
"""
backward(f::PlaneField) = _split(f, -1)

function _split(f::PlaneField{T}, s) where {T}
    Z = VACUUM_IMPEDANCE / f.n
    Eu, Ev = view(f.E, :, :, 1), view(f.E, :, :, 2)
    Hu, Hv = view(f.H, :, :, 1), view(f.H, :, :, 2)
    E = similar(f.E)
    H = similar(f.H)
    E[:, :, 1] .= (Eu .+ s * Z .* Hv) ./ 2
    E[:, :, 2] .= (Ev .- s * Z .* Hu) ./ 2
    H[:, :, 1] .= -s .* view(E, :, :, 2) ./ Z
    H[:, :, 2] .= s .* view(E, :, :, 1) ./ Z
    return PlaneField{T, typeof(E)}(E, H, f.spacing, f.origin, f.axes, f.λ, f.n, f.R)
end

end # module OpticsBase
