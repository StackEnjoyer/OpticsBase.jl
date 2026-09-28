"""
    RayBundle{N, T <: Real, P <: AbstractPort, B} <: AbstractOpticalField

A bundle of rays (optionally Gaussian beamlets) crossing a port: one monochromatic,
mutually coherent component, stored as a struct of arrays with one entry per ray.
`N = 1` is a scalar bundle, `N = 3` a vectorial one. `B` is `Nothing` for pure rays or
`Vector{SMatrix{2,2,Complex{T},4}}` for Gaussian beamlets.

See the constructor `RayBundle(port, wavelength, position, direction, opl, power, phasor)`
for the meaning, units and frames of all fields and for the invariants that are checked.

# Fields

  - `port::P`: the port at which the rays are given
  - `wavelength::T`: vacuum wavelength, \\[m\\]
  - `position::Vector{SVector{3,T}}`: ray positions on the port plane, global frame, \\[m\\]
  - `direction::Vector{SVector{3,T}}`: unit propagation directions, global frame
  - `opl::Vector{T}`: optical path length from the reference point of the component, \\[m\\]
  - `power::Vector{T}`: power per ray, \\[W\\]
  - `phasor::Vector{SVector{N,Complex{T}}}`: unit complex amplitude per ray
  - `beamlet::B`: complex beamlet matrices `Q` in \\[1/m\\], or `nothing`

# Interface

Implements the [`AbstractOpticalField`](@ref) interface: [`port`](@ref),
[`wavelength`](@ref), [`total_power`](@ref), [`is_vectorial`](@ref) (`N == 3`),
[`is_coherent`](@ref) (always `true`). Further: `length(bundle)` (number of rays) and
[`has_beamlets`](@ref).
"""
struct RayBundle{N, T <: Real, P <: AbstractPort, B} <: AbstractOpticalField
    port::P
    wavelength::T
    position::Vector{SVector{3, T}}
    direction::Vector{SVector{3, T}}
    opl::Vector{T}
    power::Vector{T}
    phasor::Vector{SVector{N, Complex{T}}}
    beamlet::B

    function RayBundle{N, T, P, B}(port::P, wavelength::T,
            position::Vector{SVector{3, T}}, direction::Vector{SVector{3, T}},
            opl::Vector{T}, power::Vector{T}, phasor::Vector{SVector{N, Complex{T}}},
            beamlet::B) where {N, T <: Real, P <: AbstractPort, B}
        B === Nothing || B === Vector{SMatrix{2, 2, Complex{T}, 4}} ||
            throw(ArgumentError("RayBundle: beamlet must be `nothing` or a Vector{SMatrix{2,2,Complex{$T},4}}, got $B"))
        _check_raybundle(port, wavelength, position, direction, opl, power, phasor, beamlet)
        return new{N, T, P, B}(port, wavelength, position, direction, opl, power, phasor,
            beamlet)
    end
end

# Invariants of `RayBundle`, see its constructor docstring.
function _check_raybundle(port, wavelength::T, position, direction, opl, power,
        phasor::Vector{SVector{N, Complex{T}}}, beamlet) where {N, T}
    tol = sqrt(eps(T))
    # 1. Number of components and equal lengths
    N == 1 || N == 3 ||
        throw(ArgumentError("RayBundle: phasor must have N = 1 (scalar) or N = 3 (vectorial) components, got N = $N"))
    M = length(position)
    lengths = (length(direction), length(opl), length(power), length(phasor))
    all(==(M), lengths) && (beamlet === nothing || length(beamlet) == M) ||
        throw(ArgumentError("RayBundle: position, direction, opl, power, phasor (and beamlet) must have equal lengths, got $((M, lengths..., beamlet === nothing ? () : length(beamlet))...)"))
    # 2. Wavelength and power
    wavelength > 0 ||
        throw(ArgumentError("RayBundle: wavelength must be positive, got $wavelength"))
    for i in 1:M
        power[i] >= 0 ||
            throw(ArgumentError("RayBundle: power must be non-negative, got $(power[i]) for ray $i"))
    end
    n = normal(port)
    r0 = origin(port)
    for i in 1:M
        d = direction[i]
        # 3. Unit direction propagating into the port
        abs(norm(d) - 1) <= tol ||
            throw(ArgumentError("RayBundle: direction must be a unit vector, got norm $(norm(d)) for ray $i"))
        dot(d, n) > 0 ||
            throw(ArgumentError("RayBundle: direction must propagate into the port (direction·normal > 0), got $(dot(d, n)) for ray $i"))
        # 4. Position on the port plane
        Δ = position[i] - r0
        abs(dot(Δ, n)) <= tol * max(norm(Δ), wavelength) ||
            throw(ArgumentError("RayBundle: position must lie on the port plane, got distance $(dot(Δ, n)) m for ray $i"))
        # 5. Unit, transverse phasor
        e = phasor[i]
        abs(norm(e) - 1) <= tol ||
            throw(ArgumentError("RayBundle: phasor must have unit norm, got norm $(norm(e)) for ray $i"))
        _check_transverse(e, d, tol, i)
        # 6. Beamlet matrix
        beamlet === nothing || _check_beamlet(beamlet[i], tol, i)
    end
    return nothing
end

_check_transverse(::SVector{1}, _, _, _) = nothing
function _check_transverse(e::SVector{3}, d, tol, i)
    # Transversality without conjugation: Σ eᵢ dᵢ = 0 for real and imaginary part
    s = transpose(e) * d
    abs(s) <= tol ||
        throw(ArgumentError("RayBundle: vectorial phasor must be transverse to the direction (phasor·direction = 0), got $s for ray $i"))
    return nothing
end

function _check_beamlet(Q::SMatrix{2, 2}, tol, i)
    abs(Q[1, 2] - Q[2, 1]) <= tol * norm(Q) ||
        throw(ArgumentError("RayBundle: beamlet matrix Q must be symmetric, got Q[1,2] - Q[2,1] = $(Q[1, 2] - Q[2, 1]) for ray $i"))
    ImQ = imag(Q)
    tr(ImQ) > 0 && det(ImQ) > 0 ||
        throw(ArgumentError("RayBundle: imaginary part of beamlet matrix Q must be positive definite (bounded beamlet), got Im Q = $(Matrix(ImQ)) for ray $i"))
    return nothing
end

"""
    RayBundle(port, wavelength, position, direction, opl, power, phasor; beamlet = nothing)

Creates a bundle of `M` rays crossing `port`: one monochromatic component whose rays are
mutually coherent. Polychromatic or mutually incoherent light is a collection of bundles.
Rays are pure rays (`beamlet = nothing`) or Gaussian beamlets.

# Arguments

  - `port`: the [`AbstractPort`](@ref) at which the rays are given, e.g. a
    [`PlanarPort`](@ref). The rays propagate downstream, into `direction·normal(port) > 0`.
  - `wavelength`: vacuum wavelength in \\[m\\]; the wavelength in the medium is
    `wavelength / refractive_index(port)`.
  - `position`: `M` points in the global frame in \\[m\\], where the rays cross the port
    plane.
  - `direction`: `M` unit propagation directions in the global frame.
  - `opl`: `M` optical path lengths `∫ n ds` in \\[m\\], absolute from the reference point
    of the coherent component (set by the first solver of the chain), not from the port
    origin. The phase accumulated along the ray is `2π·opl/wavelength`.
  - `power`: `M` powers in \\[W\\] carried by each ray (for beamlets: the total power of
    the beamlet). The bundle power is the sum, see [`total_power`](@ref).
  - `phasor`: `M` unit complex amplitudes. `N = 1`: a vector of complex numbers (phase
    factors, `|a| = 1`). `N = 3`: complex unit 3-vectors `e` in the global frame giving
    the polarization, transverse to the ray (`Σ eᵢ dᵢ = 0`, no conjugation); build them
    from Jones vectors with [`jones_to_global`](@ref). The phasor carries all phase not
    contained in `opl` (Fresnel and coating phases, Gouy phase).
  - `beamlet` (keyword): `nothing` for pure rays, or `M` beamlet parameters: complex
    symmetric 2×2 matrices `Q` in \\[1/m\\], or complex scalars `1/q` in \\[1/m\\] for
    stigmatic beamlets (converted to `Q = (1/q)·I`).

# Conventions

  - Field of ray `j` near its position, with transverse coordinates `x` in \\[m\\] in the
    [`ray_basis`](@ref) of its direction and `k = 2π·refractive_index(port)/wavelength`:
    `E(x) ∝ √power · phasor · exp(i·2π·opl/wavelength) · exp(i k xᵀQx/2)` (the last
    factor only for beamlets), with time dependence exp(−iωt) not included.
  - `Im Q` must be positive definite, so the beamlet is bounded. For a stigmatic Gaussian
    beam `1/q = 1/R + i λₘ/(π w²)` with `λₘ` the wavelength in the medium, i.e.
    `q = z − i z_R` under exp(−iωt).
  - For `N = 1` the amplitude is a scalar; there is no polarization information.

# Invariants

Each violation throws an `ArgumentError`, with `tol = √eps(T)`:

 1. `N ∈ (1, 3)`; all per-ray vectors (and `beamlet`) have the same length.
 2. `wavelength > 0`; every `power ≥ 0`.
 3. `|‖direction‖ − 1| ≤ tol`; `direction·normal(port) > 0`.
 4. The position lies on the port plane:
    `|(p − origin)·n| ≤ tol·max(‖p − origin‖, wavelength)`.
 5. `|‖phasor‖ − 1| ≤ tol`; for `N = 3` also `|Σ eᵢ dᵢ| ≤ tol`.
 6. `Q` symmetric within `tol·‖Q‖`; `Im Q` positive definite (trace and determinant > 0).

Inputs are not normalized silently. All real inputs are promoted to a common
floating-point type `T`; the storage is `Vector`s of `SVector`s (copied from the inputs).
"""
function RayBundle(port::AbstractPort, wavelength::Real, position::AbstractVector,
        direction::AbstractVector, opl::AbstractVector, power::AbstractVector,
        phasor::AbstractVector; beamlet = nothing)
    T = float(promote_type(typeof(wavelength), _eltype_real(position),
        _eltype_real(direction), _eltype_real(opl), _eltype_real(power),
        _eltype_real(phasor), _eltype_real(beamlet)))
    N = _phasor_components(phasor)
    return _raybundle(Val(N), T, port, wavelength, position, direction, opl, power, phasor,
        beamlet)
end

# Function barrier: N is known statically from here on.
function _raybundle(::Val{N}, ::Type{T}, port::P, wavelength, position, direction, opl,
        power, phasor, beamlet) where {N, T, P}
    pos = SVector{3, T}[_svector3(x, "position") for x in position]
    dir = SVector{3, T}[_svector3(x, "direction") for x in direction]
    ph = SVector{N, Complex{T}}[_phasor(Val(N), x) for x in phasor]
    Q = _beamlets(T, beamlet)
    return RayBundle{N, T, P, typeof(Q)}(port, T(wavelength), pos, dir,
        Vector{T}(opl), Vector{T}(power), ph, Q)
end

# Real element type of possibly nested / complex inputs, for type promotion.
_eltype_real(::Nothing) = Bool
_eltype_real(x::Number) = real(typeof(x))
_eltype_real(::AbstractArray{S}) where {S <: Number} = real(S)
_eltype_real(::AbstractArray{<:AbstractArray{S}}) where {S <: Number} = real(S)
_eltype_real(x::AbstractArray) = mapreduce(_eltype_real, promote_type, x; init = Bool)

# Number of components per phasor: 1 for a vector of numbers, else the (common) length.
_phasor_components(phasor::AbstractVector{<:Number}) = 1
_phasor_components(::AbstractVector{<:StaticVector{N}}) where {N} = N
function _phasor_components(phasor::AbstractVector)
    isempty(phasor) &&
        throw(ArgumentError("RayBundle: cannot infer the number of components N from an empty phasor vector; pass a vector of SVector{N}"))
    # Consistency of the remaining phasors is checked in `_phasor`
    return length(first(phasor))
end

function _svector3(x, name)
    length(x) == 3 ||
        throw(ArgumentError("RayBundle: each $name must have 3 elements, got $(length(x))"))
    return SVector{3}(x)
end

_phasor(::Val{1}, a::Number) = SVector(a)
function _phasor(::Val{N}, e) where {N}
    length(e) == N ||
        throw(ArgumentError("RayBundle: all phasors must have $N components, got $(length(e))"))
    return SVector{N}(e)
end

_beamlets(::Type{T}, ::Nothing) where {T} = nothing
function _beamlets(::Type{T}, beamlet::AbstractVector) where {T}
    return SMatrix{2, 2, Complex{T}, 4}[_beamlet_matrix(T, b) for b in beamlet]
end
_beamlet_matrix(::Type{T}, invq::Number) where {T} = SMatrix{2, 2, Complex{T}, 4}(invq, 0,
    0, invq)
function _beamlet_matrix(::Type{T}, Q::AbstractMatrix) where {T}
    size(Q) == (2, 2) ||
        throw(ArgumentError("RayBundle: beamlet matrix Q must be 2×2, got size $(size(Q))"))
    return SMatrix{2, 2, Complex{T}, 4}(Q)
end

Base.length(b::RayBundle) = length(b.position)

port(b::RayBundle) = b.port
wavelength(b::RayBundle) = b.wavelength
is_vectorial(::RayBundle{N}) where {N} = N == 3
is_coherent(::RayBundle) = true

"""
    total_power(bundle::RayBundle) -> Real

Returns the sum of the per-ray powers of `bundle` in \\[W\\] (zero for an empty bundle).
"""
total_power(b::RayBundle) = sum(b.power; init = zero(eltype(b.power)))

"""
    has_beamlets(bundle::RayBundle) -> Bool

Returns `true` if the rays of `bundle` carry Gaussian beamlet matrices `Q`, `false` for pure
rays (`beamlet = nothing`).
"""
has_beamlets(b::RayBundle) = b.beamlet !== nothing
