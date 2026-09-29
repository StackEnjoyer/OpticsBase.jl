"""
    SampledField{N,D,T,A,G,P} <: AbstractOpticalField{N}

Complex field sampled on a grid in the local coordinates of a port. `N = 1` is a scalar
field, `N = 3` a vectorial field with 3D field vectors; `D` is the grid dimension (2 for a
plane, 3 for a volume).

    SampledField(E::AbstractArray{<:Complex}, grid::AbstractGrid{D}, port::AbstractPort,
                 wavelength::Real)

Creates a field from the complex amplitude array `E` sampled on `grid` at `port`, with
vacuum wavelength `wavelength` in \\[m\\].

  - `ndims(E) == D`: scalar field (`N = 1`); `E` is reshaped to `(size(grid)..., 1)`
    without copying, so the field shares memory with the input.
  - `ndims(E) == D + 1`: `N = size(E, D + 1)`, which must be 1 or 3.

`size(E)[1:D]` must equal `size(grid)`. `E` is stored as passed (no copy, no `collect`),
so any `AbstractArray` including GPU arrays can be used. The real element type `T` is taken
from `eltype(E) == Complex{T}`, and `wavelength` is converted to `T`.

# Conventions

  - `E` is the physical complex field amplitude in \\[V/m\\] (peak, not RMS); the real field
    is `Re(E·exp(−iωt))`. `E` holds the full spatial phase (no carrier removed); position
    dependent phase terms refer to the port origin.
  - For `N = 3` the last array dimension holds the Cartesian components `(Ex, Ey, Ez)` in
    the **global** frame, not in the port frame.
  - For `N = 1` the field is a scalar amplitude with the same units and normalization.
  - Sample positions follow the grid convention of [`RegularGrid`](@ref): port-local
    coordinates `(ξ, η[, ζ])` along `(u, v[, n])`, port origin at sample `n÷2 + 1`.
  - Power: `P = κ Σ|E|² ΔξΔη` with κ = [`power_normalization`](@ref)`(port)`, see
    [`total_power`](@ref).

# Fields

  - `E`: complex amplitudes, size `(size(grid)..., N)`, \\[V/m\\], global frame
  - `grid`: the sampling grid in port-local coordinates
  - `port`: the port at which the field is given
  - `wavelength`: vacuum wavelength in \\[m\\], `> 0`

# Interface

Implements the [`AbstractOpticalField`](@ref) interface: [`port`](@ref),
[`wavelength`](@ref), [`total_power`](@ref) (planar grids only), [`is_vectorial`](@ref)
(`N == 3`) and [`is_coherent`](@ref) (always `true`: a sampled field is one coherent
component). Further accessors: [`grid`](@ref), [`field_array`](@ref).
"""
struct SampledField{N, D, T <: Real, A <: AbstractArray{Complex{T}}, G <: AbstractGrid{D},
    P <: AbstractPort} <: AbstractOpticalField{N}
    E::A
    grid::G
    port::P
    wavelength::T

    function SampledField{N}(E::AbstractArray{Complex{T}}, grid::G, port::P,
            wavelength::Real) where {N, T <: Real, D, G <: AbstractGrid{D},
            P <: AbstractPort}
        N in (1, 3) ||
            throw(ArgumentError("SampledField: number of components N must be 1 or 3, got $N"))
        ndims(E) == D + 1 ||
            throw(ArgumentError("SampledField: E must have $(D + 1) dimensions (grid..., N), got $(ndims(E))"))
        size(E, D + 1) == N ||
            throw(ArgumentError("SampledField: size(E, $(D + 1)) must be N = $N, got $(size(E, D + 1))"))
        ntuple(d -> size(E, d), Val(D)) == size(grid) ||
            throw(ArgumentError("SampledField: size(E)[1:$D] = $(ntuple(d -> size(E, d), Val(D))) does not match size(grid) = $(size(grid))"))
        wavelength > 0 ||
            throw(ArgumentError("SampledField: wavelength must be > 0, got $wavelength"))
        return new{N, D, T, typeof(E), G, P}(E, grid, port, T(wavelength))
    end
end

function SampledField(E::AbstractArray{<:Complex}, grid::AbstractGrid{D},
        port::AbstractPort, wavelength::Real) where {D}
    if ndims(E) == D
        return SampledField{1}(reshape(E, (size(E)..., 1)), grid, port, wavelength)
    elseif ndims(E) == D + 1
        return _sampled_field(Val(size(E, D + 1)), E, grid, port, wavelength)
    else
        throw(ArgumentError("SampledField: E must have $D (scalar) or $(D + 1) dimensions for a $D-dimensional grid, got $(ndims(E))"))
    end
end

_sampled_field(::Val{N}, E, grid, port, λ) where {N} = SampledField{N}(E, grid, port, λ)

"""
    grid(field::SampledField) -> AbstractGrid

Returns the sampling grid of `field`, in the port-local coordinates `(ξ, η[, ζ])` in
\\[m\\].
"""
function grid end

grid(f::SampledField) = f.grid

"""
    field_array(field::SampledField) -> AbstractArray{<:Complex}

Returns the complex amplitude array of `field` (not a copy) in \\[V/m\\], of size
`(size(grid(field))..., N)`; for `N = 3` the last dimension holds `(Ex, Ey, Ez)` in the
global frame.
"""
function field_array end

field_array(f::SampledField) = f.E

port(f::SampledField) = f.port
wavelength(f::SampledField) = f.wavelength
is_coherent(::SampledField) = true

"""
    total_power(field::SampledField) -> Real

Power of a sampled field through its port in \\[W\\]:
`P = κ Σ|E|² Δξ Δη`, summed over all samples and all `N` components, with
κ = [`power_normalization`](@ref)`(port(field))` = n/(2Z₀).

The sum is the rectangle rule for `κ ∫|E|² dA`; it is exact for fields along the port
normal and paraxial otherwise (no obliquity weighting). Only defined for planar grids
(`D = 2`); for `D = 3` an `ArgumentError` is thrown, because the power through a volume is
undefined.
"""
function total_power(f::SampledField{N, 2}) where {N}
    Δξ, Δη = spacing(f.grid)
    s = sum(abs2, f.E)
    return power_normalization(f.port) * s * Δξ * Δη
end

function total_power(::SampledField{N, D}) where {N, D}
    throw(ArgumentError("total_power is only defined for planar (D = 2) SampledFields, got D = $D: the power through a volume is undefined"))
end
