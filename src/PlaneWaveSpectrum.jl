"""
    PlaneWaveSpectrum{N, T <: Real, P <: AbstractPort} <: AbstractOpticalField{N}

A discrete spectrum of homogeneous plane waves at a port: one monochromatic, coherent
component, stored as a struct of arrays with one entry per sample (plane wave). `N = 1` is
a scalar spectrum, `N = 3` a vectorial one with 3D field vectors in the global frame.

The spectrum is the central node between the ray and the wave world: a sample is a
propagation direction with a complex amplitude and a solid angle. See the constructor
`PlaneWaveSpectrum(port, wavelength, direction, amplitude, weight)` for the model, the
units and frames of all fields and the invariants that are checked.

# Fields

  - `port::P`: the port at which the spectrum is given; `origin(port)` is the phase
    reference r₀, `refractive_index(port)` the medium
  - `wavelength::T`: vacuum wavelength, \\[m\\]
  - `direction::Vector{SVector{3,T}}`: unit propagation directions `s_j`, global frame
  - `amplitude::Vector{SVector{N,Complex{T}}}`: spectral densities `ℰ_j` per solid angle,
    \\[V/m/sr\\], global frame
  - `weight::Vector{T}`: solid angles `w_j` of the samples, \\[sr\\]

# Interface

Implements the [`AbstractOpticalField`](@ref) interface: [`port`](@ref),
[`wavelength`](@ref), [`total_power`](@ref) (exact flux), [`is_vectorial`](@ref)
(`N == 3`), [`is_coherent`](@ref) (always `true`). Further: `length(spectrum)` (number of
samples). Converters: [`PlaneWaveSummation`](@ref) evaluates the spectrum on a grid.
"""
struct PlaneWaveSpectrum{N, T <: Real, P <: AbstractPort} <: AbstractOpticalField{N}
    port::P
    wavelength::T
    direction::Vector{SVector{3, T}}
    amplitude::Vector{SVector{N, Complex{T}}}
    weight::Vector{T}

    function PlaneWaveSpectrum{N, T, P}(port::P, wavelength::T,
            direction::Vector{SVector{3, T}}, amplitude::Vector{SVector{N, Complex{T}}},
            weight::Vector{T}) where {N, T <: Real, P <: AbstractPort}
        _check_plane_wave_spectrum(port, wavelength, direction, amplitude, weight)
        return new{N, T, P}(port, wavelength, direction, amplitude, weight)
    end
end

# Invariants of `PlaneWaveSpectrum`, see its constructor docstring.
function _check_plane_wave_spectrum(port, wavelength::T, direction,
        amplitude::Vector{SVector{N, Complex{T}}}, weight) where {N, T}
    tol = sqrt(eps(T))
    # 1. Number of components and equal lengths
    N == 1 || N == 3 ||
        throw(ArgumentError("PlaneWaveSpectrum: amplitude must have N = 1 (scalar) or N = 3 (vectorial) components, got N = $N"))
    M = length(direction)
    length(amplitude) == M && length(weight) == M ||
        throw(ArgumentError("PlaneWaveSpectrum: direction, amplitude and weight must have equal lengths, got $((M, length(amplitude), length(weight)))"))
    # 2. Wavelength and weights
    wavelength > 0 ||
        throw(ArgumentError("PlaneWaveSpectrum: wavelength must be positive, got $wavelength"))
    for j in 1:M
        isfinite(weight[j]) && weight[j] > 0 ||
            throw(ArgumentError("PlaneWaveSpectrum: weight must be finite and > 0, got $(weight[j]) for sample $j"))
    end
    n = normal(port)
    for j in 1:M
        s = direction[j]
        # 3. Unit direction into the half space of the port normal
        abs(norm(s) - 1) <= tol ||
            throw(ArgumentError("PlaneWaveSpectrum: direction must be a unit vector, got norm $(norm(s)) for sample $j"))
        dot(s, n) > 0 ||
            throw(ArgumentError("PlaneWaveSpectrum: direction must propagate into the port (direction·normal > 0), got $(dot(s, n)) for sample $j"))
        # 4. Transversality (vectorial only)
        _check_spectrum_transverse(amplitude[j], s, tol, j)
    end
    return nothing
end

_check_spectrum_transverse(::SVector{1}, _, _, _) = nothing
function _check_spectrum_transverse(ℰ::SVector{3}, s, tol, j)
    # Without conjugation: Σ ℰᵢ sᵢ = 0 for real and imaginary part
    t = transpose(ℰ) * s
    abs(t) <= tol * norm(ℰ) ||
        throw(ArgumentError("PlaneWaveSpectrum: vectorial amplitude must be transverse to the direction (amplitude·direction = 0), got $t for sample $j"))
    return nothing
end

"""
    PlaneWaveSpectrum(port, wavelength, direction, amplitude, weight)

Creates a discrete spectrum of `M` homogeneous plane waves at `port`: one monochromatic
component whose samples are mutually coherent. Polychromatic or mutually incoherent light
is a collection of spectra.

# Model

With `r₀ = origin(port)`, the wavelength in the medium `λₘ = wavelength /
refractive_index(port)` and `k = 2π/λₘ`, the spectrum represents the field

```math
\\mathbf{E}(\\mathbf{r}) = \\sum_j w_j\\, \\boldsymbol{\\mathcal{E}}_j\\,
\\exp\\big(i k\\, \\mathbf{s}_j\\cdot(\\mathbf{r} - \\mathbf{r}_0)\\big)
```

in \\[V/m\\] at every global point `r` of the homogeneous medium of the port, with time
dependence exp(−iωt) not included. The sum is the discretization of the angular integral
`∫ ℰ(s) exp(i k s·(r − r₀)) dΩ` with quadrature weights `w_j`.

# Arguments

  - `port`: the [`AbstractPort`](@ref) at which the spectrum is given, e.g. a
    [`PlanarPort`](@ref). Its origin is the phase reference r₀, its refractive index the
    medium the plane waves propagate in. The normal `n` bounds the hemisphere: every
    direction satisfies `s·n > 0`. The axis `u` fixes the Jones bases (port basis and
    [`ray_basis`](@ref)) used to build or read vectorial amplitudes.
  - `wavelength`: vacuum wavelength in \\[m\\].
  - `direction`: `M` unit propagation directions `s_j` in the global frame
    (3-element vectors).
  - `amplitude`: `M` spectral densities `ℰ_j` per solid angle in \\[V/m/sr\\] (peak
    amplitude, not RMS), including the full phase at r₀. `N = 1`: a vector of complex
    numbers. `N = 3`: complex 3-vectors in the global frame, transverse to their
    direction (`Σᵢ ℰᵢ sᵢ = 0`, no conjugation).
  - `weight`: `M` solid angles `w_j` in \\[sr\\], the quadrature weights of the samples.

# Power

[`total_power`](@ref) returns the flux through any plane normal to `n`,
`P = κ λₘ² Σ_j w_j ‖ℰ_j‖²` with κ = [`power_normalization`](@ref)`(port)`: the
quadrature of the exact flux `κ λₘ² ∫ ‖ℰ‖² dΩ` of the continuous spectrum (cross terms
between distinct plane waves integrate to zero over the plane). It differs
from the paraxial power of a [`SampledField`](@ref) of the same field by O(θ²), because
`SampledField` carries no obliquity weighting.

# Invariants

Each violation throws an `ArgumentError` naming the sample index, with `tol = √eps(T)`:

 1. `N ∈ (1, 3)`; `direction`, `amplitude` and `weight` have the same length (`M = 0` is
    allowed).
 2. `wavelength > 0`; every weight is finite and `> 0`.
 3. `|‖s_j‖ − 1| ≤ tol`; `s_j·normal(port) > 0`.
 4. `N = 3`: `|Σᵢ ℰᵢ sᵢ| ≤ tol·‖ℰ‖` (transversality without conjugation).

Inputs are not normalized silently. All real inputs are promoted to a common
floating-point type `T` (the port keeps its own type); the storage is `Vector`s of
`SVector`s (copied from the inputs).
"""
function PlaneWaveSpectrum(port::AbstractPort, wavelength::Real, direction::AbstractVector,
        amplitude::AbstractVector, weight::AbstractVector)
    T = float(promote_type(typeof(wavelength), _eltype_real(direction),
        _eltype_real(amplitude), _eltype_real(weight)))
    N = _spectrum_components(amplitude)
    return _plane_wave_spectrum(Val(N), T, port, wavelength, direction, amplitude, weight)
end

# Function barrier: N is known statically from here on.
function _plane_wave_spectrum(::Val{N}, ::Type{T}, port::P, wavelength, direction,
        amplitude, weight) where {N, T, P}
    dir = SVector{3, T}[_spectrum_direction(s) for s in direction]
    amp = SVector{N, Complex{T}}[_spectrum_amplitude(Val(N), a) for a in amplitude]
    return PlaneWaveSpectrum{N, T, P}(port, T(wavelength), dir, amp, Vector{T}(weight))
end

# Number of components per amplitude: 1 for a vector of numbers, else the (common) length.
_spectrum_components(::AbstractVector{<:Number}) = 1
_spectrum_components(::AbstractVector{<:StaticVector{N}}) where {N} = N
function _spectrum_components(amplitude::AbstractVector)
    isempty(amplitude) &&
        throw(ArgumentError("PlaneWaveSpectrum: cannot infer the number of components N from an empty amplitude vector; pass a vector of numbers or of SVector{N}"))
    # Consistency of the remaining amplitudes is checked in `_spectrum_amplitude`
    return length(first(amplitude))
end

function _spectrum_direction(s)
    length(s) == 3 ||
        throw(ArgumentError("PlaneWaveSpectrum: each direction must have 3 elements, got $(length(s))"))
    return SVector{3}(s)
end

_spectrum_amplitude(::Val{1}, a::Number) = SVector(a)
function _spectrum_amplitude(::Val{N}, a) where {N}
    length(a) == N ||
        throw(ArgumentError("PlaneWaveSpectrum: all amplitudes must have $N components, got $(length(a))"))
    return SVector{N}(a)
end

Base.length(s::PlaneWaveSpectrum) = length(s.direction)

port(s::PlaneWaveSpectrum) = s.port
wavelength(s::PlaneWaveSpectrum) = s.wavelength
is_coherent(::PlaneWaveSpectrum) = true

"""
    total_power(spectrum::PlaneWaveSpectrum) -> Real

Returns the power of `spectrum` in \\[W\\]: the flux through a plane normal to the port
normal, exact for the continuous spectrum that the samples discretize (quadrature with the
weights `w_j`),

```math
P = \\kappa\\, \\lambda_m^2 \\sum_j w_j \\,\\lVert \\boldsymbol{\\mathcal{E}}_j \\rVert^2,
```

with κ = [`power_normalization`](@ref)`(port)`, `λₘ = wavelength / refractive_index(port)`
the wavelength in the medium, `w_j` in \\[sr\\] and `ℰ_j` in \\[V/m/sr\\]. Zero for an
empty spectrum. The result has the element type `T` of the spectrum.

This differs from the paraxial power `κ ∫|E|² dA` of a [`SampledField`](@ref) of the same
field by O(θ²) (no obliquity weighting there).
"""
function total_power(s::PlaneWaveSpectrum{N, T}) where {N, T}
    p = port(s)
    κ = T(power_normalization(p))
    λm = s.wavelength / T(refractive_index(p))
    acc = zero(T)
    for j in eachindex(s.weight)
        acc += s.weight[j] * sum(abs2, s.amplitude[j])
    end
    return κ * λm^2 * acc
end
