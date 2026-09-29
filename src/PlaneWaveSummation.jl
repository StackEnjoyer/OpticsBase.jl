"""
    PlaneWaveSummation(grid::AbstractGrid, port = nothing) <: AbstractFieldConverter

Converter from a [`PlaneWaveSpectrum`](@ref) to a [`SampledField`](@ref): the plane waves
are summed directly at every sample of `grid` at the output port. The result has the same
wavelength and number of components `N` as the spectrum.

# Model

With the spectrum's phase reference `r₀ = origin(port(spectrum))`, `k = 2π n/λ` in the
medium of the spectrum and the grid point `r = o + ξu + ηv (+ ζn)` of the output port
(origin `o`, axes `(u, v, n)`), the field is

```math
\\mathbf{E}(\\mathbf{r}) = \\sum_j w_j\\, \\boldsymbol{\\mathcal{E}}_j\\,
\\exp\\big(i k\\, \\mathbf{s}_j\\cdot(\\mathbf{r} - \\mathbf{r}_0)\\big)
```

in \\[V/m\\], for `N = 3` in the global frame. A 2D grid samples the port plane in
`(ξ, η)` along `(u, v)`; a 3D grid samples a volume with `ζ` along `n` (the result is
still given at the output port, whose origin is the grid center).

# Approximations

  - None beyond the discrete spectrum itself: the summation is the exact superposition of
    the plane waves in the homogeneous medium of the port, at any distance and for any
    orientation of the output port.
  - A spectrum that was sampled on a regular wave-vector grid (e.g. from an FFT) represents
    a periodic field; the result repeats with the period of that sampling.
  - The cost is O(M · N_points) for `M` plane waves; the work is done in chunks of plane
    waves as matrix products (separable in `ξ` and `η`).

# Arguments

  - `grid`: 2D grid `(ξ, η)` or 3D grid `(ξ, η, ζ)` in the local coordinates of the output
    port in \\[m\\], see [`RegularGrid`](@ref).
  - `port`: the output [`PlanarPort`](@ref), or `nothing` (default) for the port of the
    spectrum. It must have the same refractive index as the spectrum's port (the medium is
    homogeneous) and every direction must satisfy `s·normal(port) > 0`; otherwise
    `convert_field` throws an `ArgumentError`. Its axis `u` fixes the grid orientation.

Inputs other than a [`PlaneWaveSpectrum`](@ref) are rejected by `convert_field` with a
[`MissingConverterError`](@ref).

# Fields

  - `grid`: the sampling grid (`D = 2` or `3`)
  - `port`: the output port, or `nothing`
"""
struct PlaneWaveSummation{G <: AbstractGrid, P} <: AbstractFieldConverter
    grid::G
    port::P

    function PlaneWaveSummation(grid::G, port::P = nothing) where {D, G <: AbstractGrid{D},
            P <: Union{Nothing, AbstractPort}}
        D == 2 || D == 3 ||
            throw(ArgumentError("PlaneWaveSummation: grid must be 2D or 3D, got D = $D"))
        return new{G, P}(grid, port)
    end
end

input_representation(::PlaneWaveSummation) = PlaneWaveSpectrum
output_representation(::PlaneWaveSummation) = SampledField

# Number of plane waves per chunk (bounds the size of the work matrices).
const _PLANE_WAVE_CHUNK = 1024

function __convert_field(conv::PlaneWaveSummation, pws::PlaneWaveSpectrum)
    p = _require_planar(_summation_port(conv.port, pws))
    return _plane_wave_sum(conv.grid, p, pws, _PLANE_WAVE_CHUNK)
end

_summation_port(::Nothing, pws) = port(pws)
_summation_port(p::AbstractPort, _) = p

_require_planar(p::PlanarPort) = p
_require_planar(p::AbstractPort) = throw(ArgumentError("PlaneWaveSummation: the output port must be a PlanarPort, got a $(typeof(p))"))

# ζ coordinates of the grid slices: one slice at ζ = 0 for planar grids.
_normal_coordinates(::AbstractGrid{2}, ::Type{T}) where {T} = (zero(T),)
_normal_coordinates(g::AbstractGrid{3}, ::Type{T}) where {T} = T.(coordinates(g, 3))

# Sum of the plane waves of `pws` on `grid` at port `p`, in chunks of `chunk` waves.
function _plane_wave_sum(grid::AbstractGrid, p::PlanarPort, pws::PlaneWaveSpectrum{N, T},
        chunk::Integer) where {N, T}
    p_in = port(pws)
    n_in = T(refractive_index(p_in))
    n_out = T(refractive_index(p))
    isapprox(n_out, n_in; rtol = sqrt(eps(T))) ||
        throw(ArgumentError("PlaneWaveSummation: the output port must have the refractive index of the spectrum's port ($n_in), got $n_out"))
    ax = local_axes(p)
    u = SVector{3, T}(ax[:, 1])
    v = SVector{3, T}(ax[:, 2])
    nv = SVector{3, T}(ax[:, 3])
    M = length(pws)
    s = pws.direction
    for j in 1:M
        dot(s[j], nv) > 0 ||
            throw(ArgumentError("PlaneWaveSummation: every direction must propagate into the output port (direction·normal > 0), got $(dot(s[j], nv)) for sample $j"))
    end
    Δo = SVector{3, T}(origin(p)) - SVector{3, T}(origin(p_in))
    λ = pws.wavelength
    k = 2 * T(π) * n_in / λ
    ξs = T.(coordinates(grid, 1))
    ηs = T.(coordinates(grid, 2))
    ζs = _normal_coordinates(grid, T)
    Nξ, Nη, Nζ = length(ξs), length(ηs), length(ζs)

    E = zeros(Complex{T}, size(grid)..., N)
    E4 = reshape(E, Nξ, Nη, Nζ, N)     # shares memory with E; Nζ = 1 for planar grids
    for j0 in 1:chunk:M
        js = j0:min(j0 + chunk - 1, M)
        # With s·(r − r₀) = a + ξb + ηc + ζd the phase factor is separable.
        X = [cis(k * ξ * dot(s[j], u)) for ξ in ξs, j in js]               # Nξ × C
        Y = [cis(k * η * dot(s[j], v)) for η in ηs, j in js]               # Nη × C
        A = [pws.weight[j] * pws.amplitude[j][c] * cis(k * dot(s[j], Δo))
             for j in js, c in 1:N]                                         # C × N
        d = [dot(s[j], nv) for j in js]
        XA = similar(X)
        a = similar(A, length(js))
        for (l, ζ) in enumerate(ζs)
            z = cis.(k * ζ .* d)
            for c in 1:N
                a .= view(A, :, c) .* z
                XA .= X .* transpose(a)                 # X · Diagonal(a)
                mul!(view(E4, :, :, l, c), XA, transpose(Y), true, true)
            end
        end
    end
    return SampledField(E, grid, p, λ)
end
